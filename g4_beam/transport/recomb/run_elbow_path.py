"""Evaluate explicit saved route frames; no optimizer, synthetic energy or retry."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import hashlib
from pathlib import Path
import signal
import time
import traceback
import numpy as np
from ase.io import read, write
from force_bank_results import verify_force_record
from surface_calculator import calculator
from surface_io import manifest, save


class BoundedStop(BaseException):
    pass


def stop(signum, frame):
    raise BoundedStop(f'Native signal{signum}; unfinished frame is not admitted')


def record_point(atoms, directory, index, started_cpu):
    forces = atoms.get_forces(apply_constraint=False)
    energy = float(atoms.get_potential_energy())
    free = float(atoms.get_potential_energy(force_consistent=True))
    checkpoint = directory/f'point-{index:03d}.gpw'
    coordinates = directory/f'point-{index:03d}.traj'
    atoms.calc.write(str(checkpoint), mode='')
    write(coordinates, atoms)
    bank = dict(forces=forces.tolist(), energy_eV=energy,
        positions_A=atoms.positions.tolist(), gpw=manifest(checkpoint), geometry=manifest(coordinates))
    verify_force_record(bank)
    npd = len(atoms)-2
    return dict(schema='coherence-connected-elbow-point-v1', stage='completed',
        route_index=index, bank=bank, free_energy_eV=free,
        separation_A=float(atoms.get_distance(npd, npd+1)),
        height_A=float(atoms.positions[-2:, 2].mean()-atoms.positions[:-2, 2].max()),
        DFT_and_bank_process_cpu_s=time.process_time()-started_cpu)


def run(args):
    if os.sched_getaffinity(0) != {args.cpu}:
        raise RuntimeError('Declared singleton CPU required')
    raw = args.route.read_bytes()
    if hashlib.sha256(raw).hexdigest() != args.route_sha256:
        raise ValueError('Route byte identity differs')
    frames = read(args.route, index=':')
    if (not args.indices or len(set(args.indices)) != len(args.indices)
            or any(index < 0 or index >= len(frames) for index in args.indices)):
        raise ValueError('Explicit distinct existing route indices required')
    atoms = frames[args.indices[0]].copy()
    npd = len(atoms)-2
    if npd not in (6, 8) or list(atoms.symbols) != ['Pd']*npd+['H']*2:
        raise ValueError('Declared3L/4L two-D Pd111 model required')
    fixed = {int(index) for constraint in atoms.constraints for index in constraint.get_indices()}
    if fixed != {0, 1, 2, 3} or tuple(atoms.pbc) != (True, True, False):
        raise ValueError('Actual bottom-two clamp/xy periodicity required')
    for index in args.indices:
        frame = frames[index]
        if (list(frame.symbols) != list(atoms.symbols)
                or not np.array_equal(frame.cell, atoms.cell)
                or not np.array_equal(frame.pbc, atoms.pbc)
                or not np.array_equal(frame.get_masses(), atoms.get_masses())
                or not np.allclose(frame.get_masses()[-2:], 2.014101778, rtol=0, atol=1e-12)
                or not np.allclose(frame.positions[:4], atoms.positions[:4], rtol=0, atol=1e-12)):
            raise ValueError('Route changed atom/cell/isotope/clamped-coordinate identity')
    args.output.mkdir(parents=True, exist_ok=False)
    save(args.output/'input.json', dict(schema='coherence-connected-elbow-run-v1',
        route=manifest(args.route), selected_indices=args.indices, cpu=args.cpu,
        kmesh=args.kmesh, cutoff_eV=300., xc='PBE', Pd_valence=10, H_valence=1,
        convergence=dict(energy=1e-7, density=1e-6, eigenstates=1e-9),
        symmetry=dict(point_group=False, time_reversal=True),
        meaning='Explicit fixed-coordinate path upper bound, not a relaxed saddle or capture calculation',
        code=[manifest(Path(__file__)), manifest(Path(__file__).with_name('surface_calculator.py')),
              manifest(Path(__file__).with_name('force_bank_results.py')),
              manifest(Path(__file__).with_name('surface_io.py'))]))
    atoms.calc = calculator(atoms, args.output/'scf.txt', kmesh=tuple(args.kmesh))
    atoms.calc = atoms.calc.new(
        symmetry={'point_group': False, 'time_reversal': True},
        convergence={'energy': 1e-7, 'density': 1e-6, 'eigenstates': 1e-9},
        legacy_gpaw=False, txt=str(args.output/'scf.txt'))
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    completed = []
    status, error = 'path_complete', None
    print('CONNECTED_ELBOW_START', args.indices, flush=True)
    try:
        for index in args.indices:
            atoms.positions[:] = frames[index].positions
            save(args.output/'inflight.json', dict(route_index=index, admitted=False))
            started = time.process_time()
            point = record_point(atoms, args.output, index, started)
            path = args.output/f'point-{index:03d}.json'
            save(path, point)
            completed.append(manifest(path))
            print('CONNECTED_ELBOW_POINT', index, point['bank']['energy_eV'],
                  point['separation_A'], point['height_A'], flush=True)
    except BoundedStop as exc:
        status, error = 'native_bounded_stop', str(exc)
    except Exception:
        status, error = 'calculation_failure', traceback.format_exc()
        raise
    finally:
        save(args.output/'terminal.json', dict(status=status, error=error,
            completed_points=completed, process_cpu_s=time.process_time(),
            in_flight_reuse_proven=False))
        print('CONNECTED_ELBOW_END', status, len(completed), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--route', type=Path, required=True)
    parser.add_argument('--route-sha256', required=True)
    parser.add_argument('--indices', type=int, nargs='+', required=True)
    parser.add_argument('--kmesh', type=int, nargs=3, required=True)
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
