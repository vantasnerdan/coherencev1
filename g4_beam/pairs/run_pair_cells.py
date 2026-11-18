"""Relax declared cells once; preserve every completed actual compact force bank."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
from pathlib import Path
import signal
import time
import traceback
import numpy as np
from ase import Atom, Atoms
from ase.io import write
from ase.optimize import BFGS
from g4_beam.pairs.force_bank_results import verify_force_record
from g4_beam.pairs.surface_io import geometry, manifest, save


class BoundedStop(BaseException):
    pass


def stop(signum, frame):
    raise BoundedStop(f'Native signal{signum}; unfinished SCF is not admitted')


def prepare_atoms(seed):
    atoms = Atoms(seed['symbols'], positions=seed['positions_A'],
                  cell=seed['cell_A'], pbc=True)
    atoms.set_masses([2.01410177812 if s == 'H' else m
                     for s, m in zip(atoms.symbols, atoms.get_masses())])
    atoms.set_initial_magnetic_moments([0.1 if s == 'H' else 0.05 for s in atoms.symbols])
    return atoms


def record_point(atoms, directory, optimizer, started):
    step = optimizer.nsteps
    forces = atoms.get_forces(apply_constraint=False)
    energy = float(atoms.get_potential_energy())
    free = float(atoms.get_potential_energy(force_consistent=True))
    gpw = directory / f'point-{step:03d}.gpw'
    coordinates = directory / f'point-{step:03d}.traj'
    atoms.calc.write(str(gpw), mode='')
    write(coordinates, atoms)
    bank = dict(forces=forces.tolist(), energy_eV=energy,
                positions_A=atoms.positions.tolist(), gpw=manifest(gpw),
                geometry=manifest(coordinates))
    verify_force_record(bank)
    result = dict(schema='coherence-pairs-force-point-v1', step=step, bank=bank,
                  electronic_free_energy_eV=free,
                  max_unconstrained_force_eV_A=float(np.linalg.norm(forces, axis=1).max()),
                  process_cpu_s=time.process_time() - started)
    save(directory / f'point-{step:03d}.json', result)
    print('PAIR_FORCE', directory.name, step, energy,
          result['max_unconstrained_force_eV_A'], result['process_cpu_s'], flush=True)


def run_cell(atoms, settings, case, output, started):
    from g4_beam.pairs.pair_dft import calculator
    directory = output / case['name']
    directory.mkdir()
    if atoms.constraints or len(atoms) != case['n_Pd'] + case['n_D']:
        raise ValueError('Cluster must have the declared composition and no frozen coordinates')
    if list(atoms.symbols).count('H') != case['n_D']:
        raise ValueError('Incorrect cluster deuterium count')
    save(directory / 'input.json', dict(case=case, settings=settings, geometry=geometry(atoms)))
    atoms.calc = calculator(atoms, settings, directory / 'scf.txt')
    optimizer = BFGS(atoms, trajectory=str(directory / 'relax.traj'),
                     restart=str(directory / 'optimizer.json'),
                     logfile=str(directory / 'relax.log'), maxstep=0.15)
    optimizer.attach(lambda: record_point(atoms, directory, optimizer, started), interval=1)
    converged = bool(optimizer.run(fmax=settings['fmax_eV_A'], steps=settings['relax_steps']))
    result = json.loads((directory / f'point-{optimizer.nsteps:03d}.json').read_text())
    result.update(case=case, settings=settings, geometry=geometry(atoms),
                  converged_relaxation=converged, all_coordinates_mobile=True,
                  neutral_electrons=int(atoms.calc.get_number_of_electrons()))
    save(directory / 'primitive.json', result)
    if not converged:
        raise RuntimeError('Force/step cap reached: no unchanged automatic retry')
    return atoms


def run(args):
    started = time.process_time()
    if os.sched_getaffinity(0) != {args.cpu}:
        raise ValueError('Unexpected Native CPU affinity')
    config = json.loads(Path(args.input).read_text())
    os.environ['GPAW_SETUP_PATH'] = config['setup_directory']
    for name, expected in config['setup_sha256'].items():
        if manifest(Path(config['setup_directory']) / name)['sha256'] != expected:
            raise ValueError('Wrong PAW bytes: ' + name)
    output = Path(args.output)
    output.mkdir(exist_ok=False)
    sources = [manifest(args.input), manifest(__file__),
               manifest(Path(__file__).with_name('pair_dft.py')),
               manifest(verify_force_record.__code__.co_filename),
               manifest(save.__code__.co_filename)]
    save(output / 'execution.json', dict(sources=sources, affinity=[args.cpu],
         cgroup=Path('/proc/self/cgroup').read_text(),
         setup_sha256=config['setup_sha256'],
         start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())))
    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    status, exit_code = 'completed', 0
    try:
        atoms = None
        for case in config['cases']:
            if 'seed' in case:
                atoms = prepare_atoms(case['seed'])
            elif atoms is None:
                atoms = prepare_atoms(config['seed'])
            else:
                atoms.calc = None
                for position in case['add_D_positions_A']:
                    atoms.append(Atom('H', position=position, mass=2.01410177812, magmom=0.1))
            print('PAIR_CASE_BEGIN', case['name'], flush=True)
            atoms = run_cell(atoms, config['settings'], case, output, started)
            print('PAIR_CASE_END', case['name'], flush=True)
    except BoundedStop as error:
        status, exit_code = 'bounded_stop', 124
        print(str(error), flush=True)
    except BaseException:
        status, exit_code = 'failed', 1
        traceback.print_exc()
    finally:
        save(output / 'terminal.json', dict(status=status, process_cpu_s=time.process_time()-started,
             end_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())))
        print('PAIR_CELLS_END', status, flush=True)
    raise SystemExit(exit_code)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--cpu', required=True, type=int)
    run(parser.parse_args())
