"""Produce exactly one positive-amplitude NiH quartic SCF per invocation."""
import argparse
import fcntl
import os
import time
import traceback

import numpy as np

from dense_pd import ROOT, commit_json, cpu_s, manifest, read_json, require, root_action
from loading_pd import bind_ranks
from nih_dense import calculator
from nih_profile import FOLDER as BASE_FOLDER, load_point as load_base
from nih_quartic_profile import DIRECTIONS, FOLDER, SCHEMA, point_path, positions, profile


def measure(p, world, direction, run_id, affinity, restart):
    path = point_path(direction)
    root_action(world, lambda: require(not path.exists() and not path.with_suffix(".gpw").exists(),
                                      "Existing NiH quartic point requires disposition; no automatic rerun"))
    atoms, calc = calculator(p["base_profile"], world, FOLDER / f"q{direction}.{run_id}.out", restart)
    atoms.positions[:] = positions(p, direction)
    atoms.calc = calc
    started, cpu0 = time.perf_counter(), cpu_s()
    if world.rank == 0:
        print(f"NIH_QUARTIC_POINT_START direction={direction} ranks={world.size}", flush=True)
    energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces()
    electrons = float(calc.get_number_of_electrons())
    wfs = calc.dft.ibzwfs
    require(np.isfinite(energy) and np.isfinite(forces).all() and electrons == 352, "Invalid NiH quartic SCF")
    require(wfs.kpt_comm.size == world.size and wfs.domain_comm.size == wfs.band_comm.size == 1
            and len(wfs.ibz.symmetries) == 1 and calc.dft.scf_loop.update_density_and_potential,
            "Actual NiH quartic communicator/symmetry/density differs")
    n_ibz = len(calc.get_ibz_k_points())
    total_cpu = float(world.sum_scalar(cpu_s() - cpu0))
    temporary = FOLDER / f".q{direction}.{run_id}.gpw"
    calc.write(str(temporary), mode="all")
    world.barrier()
    def publish():
        checkpoint = path.with_suffix(".gpw")
        os.link(temporary, checkpoint)
        temporary.unlink()
        row = dict(schema=SCHEMA, record_type="measured_scf", profile=p, direction=direction,
            run_id=run_id, displacement_A=.24, direction_vector=p["directions"][direction],
            positions_A=atoms.positions.tolist(), energy_free_eV=energy, forces_eV_A=forces.tolist(),
            force_parallel_eV_A=float(forces[25] @ p["directions"][direction]),
            net_force_eV_A=forces.sum(axis=0).tolist(), electron_count=electrons,
            elapsed_s=time.perf_counter()-started, diagnostic_rank_cpu_s=total_cpu,
            accounting="Charge Native cgroup CPUUsageNSec, not the diagnostic rank interval",
            parallel=dict(world_size=world.size, kpt=wfs.kpt_comm.size, domain=1, band=1),
            rank_affinity=affinity, n_ibz_kpoints=n_ibz, point_group_operations=len(wfs.ibz.symmetries),
            initialization="full_self_consistent_zero_checkpoint_restart", restart_point=restart[2],
            checkpoint=manifest(checkpoint))
        commit_json(path, row)
        print(f"NIH_QUARTIC_POINT_DONE direction={direction} diagnostic_rank_cpu_s={total_cpu}", flush=True)
        return row
    return root_action(world, publish)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("direction", choices=DIRECTIONS)
    parser.add_argument("--ranks", type=int, required=True)
    args = parser.parse_args()
    os.environ["GPAW_SETUP_PATH"] = str(ROOT / "setups")
    import gpaw
    gpaw.setup_paths[:] = [str(ROOT / "setups")]
    from gpaw.mpi import world
    affinity = bind_ranks(world, args.ranks)
    p = root_action(world, profile)
    run_id = os.environ["NIH_RUN_ID"]
    require(run_id and all(c.isalnum() or c in "-_." for c in run_id), "Unsafe NiH quartic run ID")
    lock = {}
    def prepare():
        FOLDER.mkdir(parents=True, exist_ok=True)
        stream = (FOLDER / ".writer.lock").open("a")
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        lock["stream"] = stream
        path = FOLDER / "profile.json"
        if path.exists():
            require(read_json(path) == p, "Existing NiH quartic profile differs")
        else:
            commit_json(path, p)
    root_action(world, prepare)
    try:
        def seed():
            row = load_base(p["base_profile"], "zero", 0., verify_checkpoint=True)
            return BASE_FOLDER / "zero.gpw", row, manifest(BASE_FOLDER / "zero.json")
        restart = root_action(world, seed)
        measure(p, world, args.direction, run_id, affinity, restart)
        root_action(world, lambda: commit_json(FOLDER / f"stage-{args.direction}-{run_id}.json",
            dict(schema=SCHEMA, status="completed", direction=args.direction, profile_id=p["profile_id"],
                 point=manifest(point_path(args.direction)))))
        if world.rank == 0:
            print(f"NIH_QUARTIC_STAGE_COMPLETED direction={args.direction} run_id={run_id}", flush=True)
    except BaseException as error:
        if world.rank == 0:
            commit_json(FOLDER / f"failure-{args.direction}-{run_id}.json",
                        dict(schema=SCHEMA, error_type=type(error).__name__, message=str(error)))
            traceback.print_exc()
        raise
    finally:
        if world.rank == 0 and "stream" in lock:
            lock["stream"].close()


if __name__ == "__main__":
    main()
