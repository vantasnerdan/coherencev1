"""Actual Ni32H32 zero or signed source-x pair. Caller admits Native resources."""
import argparse
import fcntl
import os
import time
import traceback

import numpy as np
from ase import Atoms

from dense_pd import ROOT, commit_json, cpu_s, manifest, read_json, require, root_action
from loading_pd import bind_ranks
from nih_profile import CASES, FOLDER, SCHEMA, load_point, positions, profile


def calculator(p, world, log, restart):
    from gpaw import GPAW, PW, FermiDirac, Mixer
    parallel = dict(kpt=world.size, domain=1, band=1, sl_auto=False)
    if restart:
        checkpoint, row, _ = restart
        calc = GPAW(str(checkpoint), legacy_gpaw=False, parallel=parallel, communicator=world, txt=str(log))
        atoms = calc.get_atoms()
        require(atoms.numbers.tolist() == p["geometry"]["atomic_numbers"]
                and atoms.pbc.all() and np.allclose(atoms.cell, p["geometry"]["cell_A"], atol=1e-10, rtol=0)
                and np.allclose(atoms.positions, row["positions_A"], atol=1e-10, rtol=0)
                and calc.params.convergence == p["convergence"], "Actual NiH zero checkpoint differs")
        require(calc.dft.scf_loop.update_density_and_potential, "Frozen-density NiH restart forbidden")
        return atoms, calc
    g = p["geometry"]
    atoms = Atoms(numbers=g["atomic_numbers"], positions=g["reference_positions_A"], cell=g["cell_A"], pbc=True)
    calc = GPAW(mode=PW(p["cutoff_eV"]), setups=p["setups"], xc=p["xc"],
        kpts={"size": tuple(p["kmesh"]), "gamma": p["gamma_centered"]}, nbands=p["nbands"],
        occupations=FermiDirac(p["smearing_eV"]), spinpol=p["spinpol"], charge=p["charge"],
        symmetry=p["symmetry_settings"], mixer=Mixer(*p["mixer"]), maxiter=p["maxiter"],
        eigensolver=p["eigensolver"], convergence=p["convergence"], parallel=parallel,
        communicator=world, txt=str(log), legacy_gpaw=False)
    return atoms, calc


def measure(p, world, name, u, restart, run_id, affinity):
    path = FOLDER / f"{name}.json"
    if root_action(world, path.exists):
        return root_action(world, lambda: load_point(p, name, u, verify_checkpoint=True))
    root_action(world, lambda: require(not path.with_suffix(".gpw").exists(), "Orphan NiH checkpoint requires disposition"))
    atoms, calc = calculator(p, world, FOLDER / f"{name}.{run_id}.out", restart)
    atoms.positions[:] = positions(p, u)
    atoms.calc = calc
    started, cpu0 = time.perf_counter(), cpu_s()
    if world.rank == 0:
        print(f"NIH_POINT_START point={name} ranks={world.size}", flush=True)
    energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces()
    electrons = float(calc.get_number_of_electrons())
    wfs = calc.dft.ibzwfs
    require(np.isfinite(energy) and np.isfinite(forces).all() and electrons == p["expected_electrons"], "Invalid NiH SCF")
    require(wfs.kpt_comm.size == world.size and wfs.domain_comm.size == wfs.band_comm.size == 1
            and len(wfs.ibz.symmetries) == 1 and calc.dft.scf_loop.update_density_and_potential,
            "Actual NiH communicator/symmetry/density differs")
    n_ibz = len(calc.get_ibz_k_points())
    total_cpu = float(world.sum_scalar(cpu_s() - cpu0))
    temporary = FOLDER / f".{name}.{run_id}.gpw"
    calc.write(str(temporary), mode="all")
    world.barrier()
    def publish():
        checkpoint = path.with_suffix(".gpw")
        os.link(temporary, checkpoint)
        temporary.unlink()
        row = dict(schema=SCHEMA, record_type="measured_scf", point_id=name, profile=p,
            run_id=run_id, displacement_A=u, positions_A=atoms.positions.tolist(), energy_free_eV=energy,
            forces_eV_A=forces.tolist(), force_parallel_eV_A=float(forces[25, 0]),
            net_force_eV_A=forces.sum(axis=0).tolist(), electron_count=electrons,
            elapsed_s=time.perf_counter()-started, diagnostic_rank_cpu_s=total_cpu,
            accounting="Charge Native cgroup CPUUsageNSec, not the diagnostic rank interval",
            parallel=dict(world_size=world.size, kpt=wfs.kpt_comm.size, domain=1, band=1),
            rank_affinity=affinity, n_ibz_kpoints=n_ibz, point_group_operations=len(wfs.ibz.symmetries),
            initialization="full_self_consistent_zero_checkpoint_restart" if restart else "fresh_zero_atomic_density",
            restart_point=restart[2] if restart else None, checkpoint=manifest(checkpoint))
        commit_json(path, row)
        print(f"NIH_POINT_DONE point={name} diagnostic_rank_cpu_s={total_cpu}", flush=True)
        return row
    return root_action(world, publish)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=CASES)
    parser.add_argument("--ranks", type=int, required=True)
    args = parser.parse_args()
    os.environ["GPAW_SETUP_PATH"] = str(ROOT / "setups")
    import gpaw
    gpaw.setup_paths[:] = [str(ROOT / "setups")]
    from gpaw.mpi import world
    affinity = bind_ranks(world, args.ranks)
    p = root_action(world, profile)
    run_id = os.environ["NIH_RUN_ID"]
    require(run_id and all(c.isalnum() or c in "-_." for c in run_id), "Unsafe NiH run ID")
    lock = {}
    def prepare():
        FOLDER.mkdir(parents=True, exist_ok=True)
        stream = (FOLDER / ".writer.lock").open("a")
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        lock["stream"] = stream
        path = FOLDER / "profile.json"
        if path.exists():
            require(read_json(path) == p, "Existing NiH profile differs")
        else:
            commit_json(path, p)
    root_action(world, prepare)
    try:
        restart = None
        if args.stage == "probe":
            def seed():
                row = load_point(p, "zero", 0., verify_checkpoint=True)
                return FOLDER / "zero.gpw", row, manifest(FOLDER / "zero.json")
            restart = root_action(world, seed)
        for name, u in CASES[args.stage]:
            measure(p, world, name, u, restart, run_id, affinity)
        root_action(world, lambda: commit_json(FOLDER / f"stage-{args.stage}-{run_id}.json",
            dict(schema=SCHEMA, status="completed", stage=args.stage, profile_id=p["profile_id"],
                 points=[manifest(FOLDER / f"{name}.json") for name, _ in CASES[args.stage]])))
        if world.rank == 0:
            print(f"NIH_STAGE_COMPLETED stage={args.stage} run_id={run_id}", flush=True)
    except BaseException as error:
        if world.rank == 0:
            commit_json(FOLDER / f"failure-{args.stage}-{run_id}.json",
                        dict(schema=SCHEMA, error_type=type(error).__name__, message=str(error)))
            traceback.print_exc()
        raise
    finally:
        if world.rank == 0 and "stream" in lock:
            lock["stream"].close()


if __name__ == "__main__":
    main()
