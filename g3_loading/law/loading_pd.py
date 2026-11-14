"""Config-global registered-site SCFs and actual central probes; no well fits.

Native: run-loading.sh zero --config loading-configs/CONFIG.json --mesh 3 --ranks N
Then probe with --probe clean|near, using each config's registered unit axis.
Every displaced state reads the config-global zero.gpw; no symmetry or density transplant.
Frozen records are portable without GPW/PAW files; restart checks are stricter.
"""
import argparse
import fcntl
from importlib.metadata import version
import math
import os
from pathlib import Path
import time
import traceback
import uuid

import numpy as np
from ase.build import bulk
from ase.parallel import broadcast

from dense_pd import (ROOT, calculator, check_profile, commit_json,
                      cpu_s, digest, manifest, read_json, require, root_action)
from loading_lattice import a as lattice_a

SCHEMA = "loading-pd-local-v2"
BASE_SHA = "466191aede61708287ae5433147748a4bf113c3c50afc6251656a5fa277c9a7c"
CODE = ("loading_pd.py", "loading_readout.py", "run-loading.sh", "dense_pd.py",
        "loading_geometry.py", "loading_lattice.py", "loading-lattice-source.json", "Pd/base.json")


def config(path):
    c = read_json(path)
    require(c["schema"] == "loading-pd-config-v2", "Unsupported config schema")
    require(c["config_id"] and all(ch.isalnum() or ch in "-_." for ch in c["config_id"])
            and c["config_id"] not in (".", ".."), "Unsafe config ID")
    require(c["lattice_source"] == lattice_a(c["loading_x"])
            and c["a_A"] == c["lattice_source"]["a_angstrom"], "Lattice source differs")
    require(manifest(ROOT / "Pd/base.json")["sha256"] == BASE_SHA, "Canonical source bytes differ")
    base = read_json(ROOT / "Pd/base.json")
    canonical = bulk("PdH", "rocksalt", a=base["a_A"], cubic=True).repeat((2, 2, 2))
    require(base["atomic_numbers"] == canonical.numbers.tolist()
            and np.allclose(base["reference_positions_A"], canonical.positions, atol=1e-10, rtol=0)
            and np.allclose(base["cell_A"], canonical.cell, atol=1e-10, rtol=0)
            and base["moving_index"] == 25, "Canonical rocksalt geometry/index differs")
    g = c["geometry"]
    original = g["original_atom_indices"]
    vacancies = c["vacancy_original_indices"]
    require(len(set(vacancies)) == len(vacancies)
            and all(0 <= i < 64 and base["atomic_numbers"][i] == 1 for i in vacancies),
            "Vacancies must be distinct original H sites")
    require(original == [i for i in range(64) if i not in vacancies], "Original index map differs")
    require(g["atomic_numbers"] == [base["atomic_numbers"][i] for i in original]
            and g["pbc"] == [True] * 3, "Species/PBC differs")
    scale = c["a_A"] / base["a_A"]
    require(np.allclose(g["reference_positions_A"], canonical.positions[original] * scale, atol=1e-10, rtol=0)
            and np.allclose(g["cell_A"], canonical.cell.array * scale, atol=1e-10, rtol=0),
            "Config is not the scaled registered-site source")
    numbers = np.array(g["atomic_numbers"])
    require((numbers == 46).sum() == 32 and (numbers == 1).sum() / 32 == c["loading_x"],
            "Loading/count mismatch")
    probes = c["probes"]
    require(len(probes) == 2 and {p["probe_id"] for p in probes} == {"clean", "near"}, "Probe set differs")
    for p in probes:
        i = p["moving_index"]
        require(0 <= i < len(original) and numbers[i] == 1 and original[i] == p["original_atom_index"],
                "Probe/source map mismatch")
        offsets = canonical.positions - canonical.positions[p["original_atom_index"]]
        fractional = offsets @ np.linalg.inv(canonical.cell.array)
        radii = np.linalg.norm(offsets - np.rint(fractional) @ canonical.cell.array, axis=1)
        nn = [j for j in range(64) if base["atomic_numbers"][j] == 1
              and abs(radii[j] - base["a_A"] / np.sqrt(2)) < 1e-8]
        vacant = sorted(set(nn) & set(vacancies))
        require(len(nn) == 12 and sorted(p["NN_original_indices"]) == nn
                and sorted(p["vacant_NN_original_indices"]) == vacant
                and p["n_vac_NN"] == len(vacant), "Probe nearest-neighbour inventory differs")
        if p["probe_id"] == "clean":
            require(p["original_atom_index"] == 25 and not vacant, "Clean probe must remain original H25")
        d = np.asarray(p["direction_vector"], dtype=float)
        require(d.shape == (3,) and np.isfinite(d).all()
                and abs(np.linalg.norm(d) - 1) < 1e-12, "Probe direction is not a unit vector")
        if p["probe_id"] == "clean":
            require(np.allclose(d, [1, 0, 0], atol=1e-12, rtol=0)
                    and p["direction_label"] == "100"
                    and all(p[key] is None for key in ("vacancy_facing_original_index",
                        "vacancy_facing_image_shift", "vacancy_facing_offset_over_a")),
                    "Clean probe direction must be +100")
        else:
            require(bool(vacant), "Near probe needs an actual vacant nearest neighbour")
            target = min(vacant)
            delta = canonical.positions[target] - canonical.positions[p["original_atom_index"]]
            shift = -np.rint(delta @ np.linalg.inv(canonical.cell.array)).astype(int)
            delta = delta + shift @ canonical.cell.array
            require(p["vacancy_facing_original_index"] == target
                    and p["vacancy_facing_image_shift"] == shift.tolist()
                    and np.allclose(p["vacancy_facing_offset_over_a"], delta / base["a_A"], atol=1e-12, rtol=0)
                    and np.allclose(d, delta / np.linalg.norm(delta), atol=1e-12, rtol=0)
                    and p["direction_label"] == "vacancy_facing_110_family",
                    "Near probe direction does not point at its selected vacant neighbour")
    for source in c["source_files"]:
        require(manifest(ROOT / source["path"]) == source, "Config source bytes differ")
    return c


def profile(config_path, mesh, *, frozen_execution=None):
    require(mesh == 3, "Only admitted k3 loading profiles are supported")
    c = config(config_path)
    sources = [manifest(ROOT / name) for name in CODE] + [manifest(config_path)]
    if frozen_execution is None:
        sources += [manifest(ROOT / "setups" / name) for name in ("Pd.10.PBE.gz", "H.PBE.gz")]
        software = {name: version(name) for name in ("gpaw", "ase", "numpy", "scipy")}
    else:
        saved = frozen_execution["sources"]
        require(saved[:len(sources)] == sources, "Frozen source/config bytes differ")
        require(len(saved) == len(sources) + 2, "Missing recorded PAW provenance")
        sources, software = saved, frozen_execution["software"]
    value = dict(schema=SCHEMA, operator="LOCAL_one_D_other_nuclei_fixed_registered_sites",
                 material="PdD", electronic_species="H", vibrational_isotope="D",
                 config_id=c["config_id"], loading_x=c["loading_x"], arrangement=c["arrangement"],
                 a_A=c["a_A"], lattice_source=c["lattice_source"],
                 vacancy_original_indices=c["vacancy_original_indices"], probes=c["probes"],
                 geometry=c["geometry"], kmesh=[mesh] * 3, gamma_centered=True,
                 cutoff_eV=500, xc="PBE", setups={"Pd": "10", "H": "PBE"}, spinpol=False,
                 charge=0, expected_electrons=10 * 32 + c["geometry"]["atomic_numbers"].count(1),
                 smearing_eV=0.1, nbands=-20, mixer=[0.05, 5, 50], maxiter=80,
                 eigensolver={"name": "davidson", "niter": 2},
                 symmetry_settings={"point_group": False},
                 convergence={"energy": 1e-7, "density": 1e-4, "forces": 1e-4},
                 force_stopping_rule="GPAW maximum all-atom force change <0.0001 eV/A across last two checked cycles; not site stationarity",
                 calculator_api="gpaw.new.ase_interface.ASECalculator; legacy_gpaw=False",
                 self_consistent_density=True, energy="force_consistent_free_energy_eV",
                 sources=sources, software=software)
    value["profile_id"] = digest(value)
    return value


def folder_for(p):
    return ROOT / "loading_Pd" / p["config_id"] / "k3"


def probe_for(p, probe):
    return next(item for item in p["probes"] if item["probe_id"] == probe)


def direction_vector(p, probe):
    return np.array(probe_for(p, probe)["direction_vector"])


def cases(stage, probe):
    if stage == "zero":
        return [("zero", None, None, 0.0)]
    require(probe in ("clean", "near"), "Displaced stages require --probe clean|near")
    require(stage == "probe", "Only the registered single-axis probe stage is supported")
    return [(f"{probe}_probe_{sign}", probe, "registered_probe_axis", u)
            for sign, u in (("minus", -0.03), ("plus", 0.03))]


def point_positions(p, probe, direction, u):
    xyz = np.array(p["geometry"]["reference_positions_A"])
    if probe is not None:
        require(direction == "registered_probe_axis", "Wrong registered probe direction")
        xyz[probe_for(p, probe)["moving_index"]] += u * direction_vector(p, probe)
    return xyz


def load_point(folder, p, name, probe, direction, u, *, verify_checkpoint=False):
    path = Path(folder) / f"{name}.json"
    row = read_json(path)
    require(row["schema"] == SCHEMA and row["record_type"] == "measured_scf", "Not a loading SCF")
    check_profile(row["profile"], p)
    require((row["point_id"], row["probe_id"], row["direction"], row["displacement_A"])
            == (name, probe, direction, u), "Point identity differs")
    moving = probe_for(p, probe)["moving_index"] if probe else None
    require(row["moving_index"] == moving, "Point moving index differs")
    require(row["direction_vector"] == (probe_for(p, probe)["direction_vector"] if probe else None)
            and row["direction_label"] == (probe_for(p, probe)["direction_label"] if probe else None),
            "Point direction metadata differs")
    require(np.allclose(row["positions_A"], point_positions(p, probe, direction, u), atol=1e-10, rtol=0),
            "Point positions differ")
    force = np.asarray(row["forces_eV_A"])
    require(force.shape == (len(p["geometry"]["atomic_numbers"]), 3)
            and np.isfinite(force).all() and math.isfinite(row["energy_free_eV"]), "Incomplete SCF values")
    require(row["electron_count"] == p["expected_electrons"], "Electron count differs")
    parallel = row["parallel"]
    require(parallel["world_size"] > 0 and parallel["kpt"] == parallel["world_size"]
            and parallel["domain"] == parallel["band"] == 1
            and row["point_group_operations"] == 1, "Actual execution layout/symmetry differs")
    affinity = row["rank_affinity"]
    require(len(affinity) == parallel["world_size"]
            and all(item["rank"] == i and len(item["actual_cpus"]) == 1
                    and set(item["actual_cpus"]) <= set(item["inherited_cpus"]) for i, item in enumerate(affinity))
            and len({item["actual_cpus"][0] for item in affinity}) == len(affinity), "Rank affinity differs")
    if probe:
        require(np.allclose(row["force_source_eV_A"], force[moving], atol=1e-12, rtol=0)
                and abs(row["force_parallel_eV_A"] - force[moving] @ direction_vector(p, probe)) < 1e-12,
                "Force projection differs")
        require(row["restart_point"] == manifest(Path(folder) / "zero.json")
                and row["initialization"] == "full_self_consistent_zero_checkpoint_restart", "Not common-zero restart")
    else:
        require(row["restart_point"] is None and row["initialization"] == "fresh_zero_atomic_density",
                "Zero must be independently initialized")
    if verify_checkpoint:
        require(row["checkpoint"] == manifest(path.with_suffix(".gpw")), "Checkpoint bytes differ")
    return row


def zero_restart(folder, p):
    row = load_point(folder, p, "zero", None, None, 0.0, verify_checkpoint=True)
    path = folder / "zero.json"
    return path.with_suffix(".gpw"), row, manifest(path)


def measure(folder, p, world, case, restart, run_id, affinity):
    name, probe, direction, u = case
    path = folder / f"{name}.json"
    if root_action(world, path.exists):
        return root_action(world, lambda: load_point(folder, p, *case, verify_checkpoint=True))
    root_action(world, lambda: require(not path.with_suffix(".gpw").exists(),
                                      "Orphan checkpoint requires explicit operator disposition"))
    atoms, calc = calculator(p, world, folder / f"{name}.{run_id}.out", restart)
    atoms.positions[:] = point_positions(p, probe, direction, u)
    atoms.calc = calc
    started, cpu0 = time.perf_counter(), cpu_s()
    if world.rank == 0:
        print(f"LOADING_PD_POINT_START config={p['config_id']} point={name} ranks={world.size}", flush=True)
    energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces()
    electrons = float(calc.get_number_of_electrons())
    require(math.isfinite(energy) and np.isfinite(forces).all(), "Nonfinite SCF")
    require(abs(electrons - p["expected_electrons"]) < 1e-8, "Wrong 10*NPd+NH electron count")
    wfs = calc.dft.ibzwfs
    require(wfs.kpt_comm.size == world.size and wfs.domain_comm.size == wfs.band_comm.size == 1,
            "Actual communicators are not k-point-only")
    require(calc.dft.scf_loop.update_density_and_potential, "Frozen-density calculation forbidden")
    operations = len(wfs.ibz.symmetries)
    require(operations == 1, "Spatial point-group reduction must remain disabled")
    total_cpu = float(world.sum_scalar(cpu_s() - cpu0))
    n_ibz = len(calc.get_ibz_k_points())
    temporary = folder / f".{name}.{run_id}.gpw"
    calc.write(str(temporary), mode="all")
    world.barrier()
    moving = probe_for(p, probe)["moving_index"] if probe else None
    def publish():
        checkpoint = path.with_suffix(".gpw")
        os.link(temporary, checkpoint)
        temporary.unlink()
        row = dict(schema=SCHEMA, record_type="measured_scf", point_id=name, profile=p,
                   run_id=run_id, probe_id=probe, moving_index=moving, direction=direction,
                   direction_vector=probe_for(p, probe)["direction_vector"] if probe else None,
                   direction_label=probe_for(p, probe)["direction_label"] if probe else None,
                   displacement_A=u, positions_A=atoms.positions.tolist(), energy_free_eV=energy,
                   forces_eV_A=forces.tolist(), electron_count=electrons,
                   force_source_eV_A=forces[moving].tolist() if probe else None,
                   force_parallel_eV_A=float(forces[moving] @ direction_vector(p, probe)) if probe else None,
                   net_force_eV_A=forces.sum(axis=0).tolist(), elapsed_s=time.perf_counter() - started,
                   cpu_s_rank_processes_scientific_interval=total_cpu,
                   accounting="Diagnostic energy/force rank CPU only; charge whole Native cgroup CPUUsageNSec including failures",
                   parallel=dict(world_size=world.size, kpt=wfs.kpt_comm.size, domain=1, band=1),
                   rank_affinity=affinity, n_ibz_kpoints=n_ibz, point_group_operations=operations,
                   initialization="full_self_consistent_zero_checkpoint_restart" if restart else "fresh_zero_atomic_density",
                   restart_point=restart[2] if restart else None, checkpoint=manifest(checkpoint))
        commit_json(path, row)
        print(f"LOADING_PD_POINT_DONE point={name} diagnostic_rank_cpu_s={total_cpu}", flush=True)
        return row
    return root_action(world, publish)


def bind_ranks(world, ranks):
    require(world.size == ranks and ranks > 0, "Requested/actual world size differs")
    sizes = [int(os.environ[key]) for key in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_SIZE") if key in os.environ]
    require(all(size == world.size for size in sizes), "Launcher/GPAW world differs")
    inherited = sorted(os.sched_getaffinity(0))
    masks = [broadcast(inherited if world.rank == rank else None, root=rank, comm=world)
             for rank in range(world.size)]
    require(all(mask == inherited for mask in masks) and len(inherited) >= world.size,
            "Ranks require the same admitted inherited CPU set with one core each")
    os.sched_setaffinity(0, {inherited[world.rank]})
    actual = sorted(os.sched_getaffinity(0))
    require(actual == [inherited[world.rank]], "Actual pinning differs")
    local = dict(rank=world.rank, inherited_cpus=inherited, actual_cpus=actual)
    affinity = [broadcast(local if world.rank == rank else None, root=rank, comm=world)
                for rank in range(world.size)]
    print(f"LOADING_PD_RANK_AFFINITY rank={world.rank} cpus={actual}", flush=True)
    return affinity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("zero", "probe"))
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--mesh", type=int, choices=(3,), required=True)
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--probe", choices=("clean", "near"))
    args = parser.parse_args()
    require(args.stage != "zero" or args.probe is None, "Zero is config-global; omit --probe")
    selected = cases(args.stage, args.probe)
    os.environ["GPAW_SETUP_PATH"] = str(ROOT / "setups")
    import gpaw
    gpaw.setup_paths[:] = [str(ROOT / "setups")]
    from gpaw.mpi import world
    affinity = bind_ranks(world, args.ranks)
    p = root_action(world, lambda: profile(args.config, args.mesh))
    folder = folder_for(p)
    run_id = os.environ.get("LOADING_PD_RUN_ID", uuid.uuid4().hex)
    require(run_id and all(ch.isalnum() or ch in "-_." for ch in run_id), "Unsafe run ID")
    lock = None
    def prepare():
        nonlocal lock
        folder.mkdir(parents=True, exist_ok=True)
        lock = (folder / ".writer.lock").open("a")
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        path = folder / "profile.json"
        if path.exists():
            check_profile(read_json(path), p)
        else:
            commit_json(path, p)
    root_action(world, prepare)
    try:
        restart = root_action(world, lambda: zero_restart(folder, p)) if args.stage != "zero" else None
        for case in selected:
            measure(folder, p, world, case, restart, run_id, affinity)
        tag = args.stage + (f"-{args.probe}" if args.probe else "")
        root_action(world, lambda: commit_json(folder / f"stage-{tag}-{run_id}.json",
                    dict(schema=SCHEMA, stage=args.stage, probe_id=args.probe, status="completed",
                         profile_id=p["profile_id"], run_id=run_id,
                         points=[manifest(folder / f"{case[0]}.json") for case in selected])))
        if world.rank == 0:
            print(f"LOADING_PD_STAGE_COMPLETED config={p['config_id']} stage={tag} run_id={run_id}", flush=True)
    except BaseException as error:
        if world.rank == 0:
            commit_json(folder / f"failure-{args.stage}-{args.probe}-{run_id}.json",
                        dict(schema=SCHEMA, stage=args.stage, probe_id=args.probe, run_id=run_id,
                             status="failed", profile_id=p["profile_id"],
                             type=type(error).__name__, message=str(error)))
        raise
    finally:
        if lock is not None:
            lock.close()


if __name__ == "__main__":
    try:
        main()
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
