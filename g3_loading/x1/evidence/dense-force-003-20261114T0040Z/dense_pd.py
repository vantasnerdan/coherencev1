"""PdD LOCAL frozen-host 64-site sampling (never a primitive collective mode).

Run with run-dense-pd.sh STAGE --mesh {3,4} --a 4.084 [--ranks N].
STAGE is zero, probe, well100 or well111. Each mesh starts with an actual
undisplaced SCF. Every displaced point restarts that same zero checkpoint;
cross-mesh restart is deliberately unsupported. Spatial point-group reduction
is disabled throughout, so zero's symmetries cannot constrain 100/111 motion.
The probe contains actual -0.03 and +0.03 A source-x SCFs, not reconstructed signs.

Each point is committed as POINT.json after POINT.gpw; both are immutable.
Re-running skips valid committed points, but refuses partial/orphan checkpoints.
A single writer lock covers a mesh folder. No automatic retry, mesh adoption,
resource admission or launch occurs here. SIGKILL is reported by the caller's
systemd ExecStopPost=unit-finish.sh protocol, not by a Python/shell finally.
"""
import argparse
import fcntl
import hashlib
from importlib.metadata import version
import json
import math
import os
from pathlib import Path
import resource
import time
import traceback
import uuid

import numpy as np
from ase import Atoms
from ase.build import bulk

ROOT = Path(__file__).resolve().parent
SCHEMA = "dense-pd-local-v1"
DIRECTIONS = {"100": np.array([1., 0., 0.]), "111": np.ones(3) / np.sqrt(3)}
CASES = {"zero": [("zero", "100", 0.)],
         "probe": [("probe_minus", "100", -0.03), ("probe", "100", 0.03)],
         "well100": [("well100_012", "100", 0.12), ("well100_024", "100", 0.24)],
         "well111": [("well111_012", "111", 0.12), ("well111_024", "111", 0.24)]}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def manifest(path):
    path = Path(path).resolve()
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return dict(path=path.relative_to(ROOT).as_posix(), bytes=path.stat().st_size, sha256=h.hexdigest())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(Path(path).read_text())


def commit_json(path, value):
    """Atomic no-clobber publication; a killed writer cannot expose half JSON."""
    path = Path(path)
    temporary = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with temporary.open("x") as stream:
            stream.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def source_files():
    paths = [ROOT / name for name in ("dense_pd.py", "dense_pd_readout.py", "run-dense-pd.sh",
                                     "well.py", "analyse.py", "Pd/base.json")]
    paths += [ROOT / "setups" / name for name in ("Pd.10.PBE.gz", "H.PBE.gz")]
    return [manifest(path) for path in paths]


def reference(a):
    base = read_json(ROOT / "Pd/base.json")
    canonical = bulk("PdH", "rocksalt", a=base["a_A"], cubic=True).repeat((2, 2, 2))
    require(base["atomic_numbers"] == canonical.numbers.tolist(), "Baseline species/order mismatch")
    require(np.allclose(base["reference_positions_A"], canonical.positions, atol=1e-10, rtol=0)
            and np.allclose(base["cell_A"], canonical.cell, atol=1e-10, rtol=0),
            "Baseline is not the expected 2x2x2 conventional rocksalt cell")
    factor = a / base["a_A"]
    atoms = Atoms(numbers=base["atomic_numbers"],
                  positions=np.array(base["reference_positions_A"]) * factor,
                  cell=np.array(base["cell_A"]) * factor, pbc=True)
    require(len(atoms) == 64 and (atoms.numbers == 1).sum() == 32
            and (atoms.numbers == 46).sum() == 32, "Only stoichiometric Pd32D32 is supported")
    moving = base["moving_index"]
    require(atoms.numbers[moving] == 1, "Source is not an H/D site")
    return dict(atomic_numbers=atoms.numbers.tolist(), cell_A=atoms.cell.tolist(),
                reference_positions_A=atoms.positions.tolist(), pbc=[True] * 3,
                moving_index=moving, H_indices=np.flatnonzero(atoms.numbers == 1).tolist())


def profile(mesh, a, *, frozen_execution=None):
    require(math.isfinite(a) and a > 0, "--a must be finite and positive")
    value = dict(schema=SCHEMA, operator="LOCAL_one_D_other_nuclei_fixed", material="PdD",
                 electronic_species="H", vibrational_isotope="D", loading_x=1,
                 a_A=a, kmesh=[mesh] * 3, gamma_centered=True, cutoff_eV=500,
                 xc="PBE", setups={"Pd": "10", "H": "PBE"}, spinpol=False,
                 smearing_eV=0.1, nbands=-20, mixer=[0.05, 5, 50], maxiter=80,
                 convergence={"energy": 1e-7, "density": 1e-4, "forces": 0.003},
                 eigensolver={"name": "davidson", "niter": 2},
                 symmetry="point_group_disabled; common zero checkpoint for all displaced directions",
                 symmetry_settings={"point_group": False},
                 force_stopping_rule="GPAW Forces: maximum all-atom force change <0.003 eV/A over last two checked cycles",
                 calculator_api="gpaw.new.ase_interface.ASECalculator; legacy_gpaw=False",
                 self_consistent_density=True, energy="force_consistent_free_energy_eV",
                 geometry=reference(a),
                 sources=source_files() if frozen_execution is None else frozen_execution["sources"],
                 software={name: version(name) for name in ("gpaw", "ase", "numpy", "scipy")}
                 if frozen_execution is None else frozen_execution["software"])
    value["profile_id"] = digest(value)
    return value


def check_profile(value, expected):
    unsigned = {key: item for key, item in value.items() if key != "profile_id"}
    require(value.get("profile_id") == digest(unsigned), "Saved profile identity checksum differs")
    require(value == expected, "Profile/geometry/source identity differs; refusing mixed frozen output")


def point_positions(p, direction, u):
    xyz = np.array(p["geometry"]["reference_positions_A"])
    xyz[p["geometry"]["moving_index"]] += u * DIRECTIONS[direction]
    return xyz


def validate_point(row, p, name, direction, u, path, *, verify_checkpoint=True):
    require(row["schema"] == SCHEMA and row["record_type"] == "measured_scf", "Not a measured local SCF")
    check_profile(row["profile"], p)
    require((row["point_id"], row["direction"], row["displacement_A"]) == (name, direction, u),
            "Point name/direction/displacement mismatch")
    require(np.allclose(row["positions_A"], point_positions(p, direction, u), atol=1e-10, rtol=0),
            "Point geometry mismatch")
    force = np.asarray(row["forces_eV_A"])
    require(force.shape == (64, 3) and np.isfinite(force).all()
            and math.isfinite(row["energy_free_eV"]), "Nonfinite or incomplete SCF output")
    moving = p["geometry"]["moving_index"]
    require(np.allclose(row["force_source_eV_A"], force[moving], atol=1e-12, rtol=0)
            and abs(row["force_parallel_eV_A"] - force[moving] @ DIRECTIONS[direction]) < 1e-12,
            "Inconsistent saved force projections")
    if verify_checkpoint:
        checkpoint = Path(path).with_suffix(".gpw")
        require(row["checkpoint"] == manifest(checkpoint), "Checkpoint byte/hash binding differs")
    return row


def load_point(folder, p, name, direction, u, *, verify_checkpoint=False):
    path = Path(folder) / f"{name}.json"
    return validate_point(read_json(path), p, name, direction, u, path, verify_checkpoint=verify_checkpoint)


def inversion_permutation(geometry):
    positions = np.array(geometry["reference_positions_A"])
    cell = np.array(geometry["cell_A"])
    numbers = np.array(geometry["atomic_numbers"])
    source = geometry["moving_index"]
    inverse_cell = np.linalg.inv(cell)
    permutation = []
    for i, position in enumerate(positions):
        target = 2 * positions[source] - position
        delta = (positions - target) @ inverse_cell
        delta -= np.rint(delta)
        matches = np.flatnonzero((np.linalg.norm(delta @ cell, axis=1) < 1e-8) & (numbers == numbers[i]))
        require(len(matches) == 1, f"Reference inversion match not unique for atom {i}")
        permutation.append(int(matches[0]))
    require(sorted(permutation) == list(range(len(positions)))
            and all(permutation[permutation[i]] == i for i in range(len(positions)))
            and permutation[source] == source, "Inversion is not an involutive site permutation")
    return permutation


def root_action(world, function):
    from ase.parallel import broadcast
    result = None
    if world.rank == 0:
        try:
            result = (True, function())
        except Exception as error:
            result = (False, f"{type(error).__name__}: {error}")
    okay, value = broadcast(result, comm=world)
    if not okay:
        raise RuntimeError(value)
    return value


def cpu_s():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_utime + usage.ru_stime


def restart_record(path, p, expected_name, direction, u):
    path = Path(path).resolve()
    row = read_json(path)
    validate_point(row, p, expected_name, direction, u, path)
    return path.with_suffix(".gpw"), row, manifest(path)


def calculator(p, world, logfile, restart=None):
    from gpaw import GPAW, PW, FermiDirac, Mixer
    parallel = dict(kpt=world.size, domain=1, band=1, sl_auto=False)
    if restart is not None:
        checkpoint, row, _ = restart
        # New GPAW allows only parallel layout changes when reading a checkpoint.
        # Physics/mesh are already bound to this exact same-profile point.
        calc = GPAW(str(checkpoint), legacy_gpaw=False,
                    parallel=parallel, communicator=world, txt=str(logfile))
        atoms = calc.get_atoms()
        require(atoms.numbers.tolist() == p["geometry"]["atomic_numbers"]
                and atoms.pbc.all()
                and np.allclose(atoms.cell, p["geometry"]["cell_A"], atol=1e-10, rtol=0)
                and np.allclose(atoms.positions, row["positions_A"], atol=1e-10, rtol=0),
                "Actual checkpoint geometry differs from its bound record")
        require(calc.dft.scf_loop.update_density_and_potential, "Frozen-density restart is forbidden")
        return atoms, calc
    geometry = p["geometry"]
    atoms = Atoms(numbers=geometry["atomic_numbers"], positions=geometry["reference_positions_A"],
                  cell=geometry["cell_A"], pbc=True)
    calc = GPAW(mode=PW(p["cutoff_eV"]), setups={"Pd": "10"}, xc="PBE",
                kpts={"size": tuple(p["kmesh"]), "gamma": True}, nbands=-20,
                occupations=FermiDirac(0.1), spinpol=False, symmetry=p["symmetry_settings"],
                mixer=Mixer(0.05, 5, 50), maxiter=80,
                eigensolver=p["eigensolver"],
                convergence=p["convergence"], parallel=parallel, communicator=world,
                txt=str(logfile), legacy_gpaw=False)
    return atoms, calc


def measure(folder, p, world, name, direction, u, restart, run_id):
    path = folder / f"{name}.json"
    exists = root_action(world, path.exists)
    if exists:
        return root_action(world, lambda: load_point(folder, p, name, direction, u, verify_checkpoint=True))
    root_action(world, lambda: require(not path.with_suffix(".gpw").exists(),
                                      f"Orphan checkpoint {path.with_suffix('.gpw')}; explicit operator disposition required"))
    atoms, calc = calculator(p, world, folder / f"{name}.{run_id}.out", restart)
    atoms.positions[:] = point_positions(p, direction, u)
    atoms.calc = calc
    started, cpu0 = time.perf_counter(), cpu_s()
    if world.rank == 0:
        print(f"DENSE_PD_POINT_START mesh={p['kmesh'][0]} point={name} ranks={world.size}", flush=True)
    energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces()
    require(math.isfinite(energy) and np.isfinite(forces).all(), "SCF produced nonfinite values")
    require(abs(calc.get_number_of_electrons() - 352) < 1e-8, "Expected 352 electrons (10-valence Pd)")
    ibzwfs = calc.dft.ibzwfs
    require(ibzwfs.kpt_comm.size == world.size and ibzwfs.domain_comm.size == 1
            and ibzwfs.band_comm.size == 1, "Actual communicator layout is not k-point-only")
    require(calc.dft.scf_loop.update_density_and_potential, "Full self-consistent density is required")
    local_cpu = cpu_s() - cpu0
    total_cpu = float(world.sum_scalar(local_cpu))
    moving = p["geometry"]["moving_index"]
    temporary = folder / f".{name}.{run_id}.gpw"
    calc.write(str(temporary), mode="all")
    world.barrier()
    checkpoint = path.with_suffix(".gpw")
    def publish():
        os.link(temporary, checkpoint)
        temporary.unlink()
        row = dict(schema=SCHEMA, record_type="measured_scf", point_id=name, profile=p,
                   run_id=run_id, direction=direction, displacement_A=u,
                   positions_A=atoms.positions.tolist(), energy_free_eV=energy,
                   forces_eV_A=forces.tolist(), force_source_eV_A=forces[moving].tolist(),
                   force_parallel_eV_A=float(forces[moving] @ DIRECTIONS[direction]),
                   net_force_eV_A=forces.sum(axis=0).tolist(),
                   elapsed_s=time.perf_counter() - started,
                   cpu_s_rank_processes_scientific_interval=total_cpu,
                   accounting="Diagnostic rank-process CPU during energy/force only; NOT charged Native cost. "
                              "Use external Native/systemd cgroup CPUUsageNSec, including startup and failed runs.",
                   parallel=dict(world_size=world.size, kpt=ibzwfs.kpt_comm.size,
                                 domain=ibzwfs.domain_comm.size, band=ibzwfs.band_comm.size),
                   n_ibz_kpoints=len(calc.get_ibz_k_points()),
                   point_group_operations=len(ibzwfs.ibz.symmetries),
                   initialization="full_self_consistent_zero_checkpoint_restart" if restart else "fresh_zero_atomic_density",
                   restart_point=restart[2] if restart else None, checkpoint=manifest(checkpoint))
        commit_json(path, row)
        print(f"DENSE_PD_POINT_DONE point={name} diagnostic_rank_cpu_s={total_cpu}", flush=True)
        return row
    return root_action(world, publish)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=tuple(CASES))
    parser.add_argument("--mesh", type=int, choices=(3, 4), required=True)
    parser.add_argument("--a", type=float, required=True, help="Explicit conventional lattice parameter in Angstrom")
    parser.add_argument("--ranks", type=int, default=1, help="Required actual GPAW world size; no silent serial duplication")
    args = parser.parse_args()
    require(args.ranks > 0, "--ranks must be positive")
    os.environ["GPAW_SETUP_PATH"] = str(ROOT / "setups")
    import gpaw
    # The CLI may have imported GPAW before this script: pin its live search path too.
    gpaw.setup_paths[:] = [str(ROOT / "setups")]
    from gpaw.mpi import world
    require(world.size == args.ranks, f"Refusing GPAW world.size={world.size}, requested {args.ranks}")
    # Refuse direct mpiexec of plain Python yielding N unrelated serial worlds.
    launcher_sizes = [int(os.environ[key]) for key in ("OMPI_COMM_WORLD_SIZE", "PMI_SIZE", "PMIX_SIZE")
                      if key in os.environ]
    require(all(size == world.size for size in launcher_sizes), "MPI launcher/GPAW world mismatch")
    allowed_cpus = sorted(os.sched_getaffinity(0))
    require(len(allowed_cpus) >= world.size, "Insufficient admitted CPUs for one pinned rank each")
    os.sched_setaffinity(0, {allowed_cpus[world.rank]})
    print(f"DENSE_PD_RANK_AFFINITY rank={world.rank} cpus={sorted(os.sched_getaffinity(0))}", flush=True)
    p = root_action(world, lambda: profile(args.mesh, args.a))
    root_action(world, lambda: inversion_permutation(p["geometry"]))
    folder = ROOT / "dense_Pd" / f"k{args.mesh}"
    run_id = os.environ.get("DENSE_PD_RUN_ID", uuid.uuid4().hex)
    require(run_id and all(c.isalnum() or c in "-_." for c in run_id), "Unsafe run ID")
    lock = None
    def prepare():
        nonlocal lock
        folder.mkdir(parents=True, exist_ok=True)
        lock = (folder / ".writer.lock").open("a")
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        profile_path = folder / "profile.json"
        if profile_path.exists():
            check_profile(read_json(profile_path), p)
        else:
            commit_json(profile_path, p)
    root_action(world, prepare)
    try:
        restart = None
        if args.stage != "zero":
            restart = root_action(world, lambda: restart_record(folder / "zero.json", p, "zero", "100", 0.0))
        for name, direction, u in CASES[args.stage]:
            measure(folder, p, world, name, direction, u, restart, run_id)
        root_action(world, lambda: commit_json(folder / f"stage-{args.stage}-{run_id}.json",
                    dict(schema=SCHEMA, stage=args.stage, status="completed", profile_id=p["profile_id"],
                         run_id=run_id, points=[manifest(folder / f"{name}.json") for name, _, _ in CASES[args.stage]])))
        if world.rank == 0:
            print(f"DENSE_PD_STAGE_COMPLETED stage={args.stage} mesh={args.mesh} run_id={run_id}", flush=True)
    except BaseException as error:
        if world.rank == 0:
            commit_json(folder / f"failure-{args.stage}-{run_id}.json",
                        dict(schema=SCHEMA, stage=args.stage, run_id=run_id, status="failed",
                             profile_id=p["profile_id"], type=type(error).__name__, message=str(error)))
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
