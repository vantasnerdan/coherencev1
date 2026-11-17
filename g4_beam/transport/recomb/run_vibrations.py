"""Serial finite-displacement forces from an explicitly bound completed GPW."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
from pathlib import Path
import time
import numpy as np
from vibration_bank import RecordedVibrations
from force_bank_results import verify_force_record
from gpaw import GPAW
from gpaw.mpi import serial_comm
from hessian_modes import central_hessian, mass_spectrum
from surface_io import geometry, manifest, save


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--checkpoint-sha256', required=True)
    parser.add_argument('--indices', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--delta', type=float, default=0.01)
    args = parser.parse_args()
    if os.sched_getaffinity(0) != {args.cpu}:
        raise ValueError('Declared singleton CPU affinity required')
    checkpoint = manifest(args.checkpoint)
    if checkpoint['sha256'] != args.checkpoint_sha256:
        raise ValueError('Completed checkpoint identity changed')
    indices = [int(i) for i in args.indices.split(',')]
    if len(set(indices)) != len(indices) or not indices or args.delta <= 0:
        raise ValueError('Distinct active atoms and positive displacement required')
    args.output.mkdir(parents=True, exist_ok=True)
    source_calc = GPAW(str(args.checkpoint), communicator=serial_comm,
                       legacy_gpaw=False, txt=None)
    atoms = source_calc.get_atoms()
    # Displacements break the saved point-group symmetry. Keep the source
    # Hamiltonian/geometry, but keep only the valid time-reversal reduction.
    atoms.calc = source_calc.new(
        symmetry={'point_group': False, 'time_reversal': True},
        convergence={'energy': 1e-7, 'density': 1e-6, 'eigenstates': 1e-9},
        legacy_gpaw=False, txt=str(args.output / 'scf.txt'))
    del source_calc
    atoms.set_masses([2.014101778 if s == 'H' else m
                      for s, m in zip(atoms.symbols, atoms.get_masses())])
    fixed = set()
    for constraint in atoms.constraints:
        if constraint.todict()['name'] != 'FixAtoms':
            raise ValueError('Only the declared fixed-Cartesian constraint is supported')
        fixed.update(constraint.get_indices())
    if any(i < 0 or i >= len(atoms) or i in fixed for i in indices):
        raise ValueError('Displaced atoms must be in the mobile Cartesian subspace')
    inputs = {'checkpoint': checkpoint, 'geometry': geometry(atoms),
              'active_atoms': indices, 'masses_u': atoms.get_masses()[indices].tolist(),
              'delta_A': args.delta, 'nfree': 2,
              'symmetry': {'point_group': False, 'time_reversal': True},
              'electronic_state': 'fresh SCF from source geometry/parameters',
              'convergence': {'energy': 1e-7, 'density': 1e-6, 'eigenstates': 1e-9},
              'mode_scope': 'Only supplied atoms displaced; all other positions '
                            'held for this local harmonic correction',
              'code': [manifest(Path(__file__).with_name(name)) for name in
                       ('run_vibrations.py', 'hessian_modes.py', 'surface_io.py',
                        'vibration_bank.py', 'force_bank_results.py')]}
    input_path = args.output / 'input.json'
    if input_path.exists() and json.loads(input_path.read_text()) != inputs:
        raise ValueError('Existing displacement cache has different bound inputs')
    if not input_path.exists():
        save(input_path, inputs)
    vibrations = RecordedVibrations(
        atoms, indices=indices, delta=args.delta, nfree=2,
        name=str(args.output / 'cache'), bank_directory=args.output / 'banks')
    # ASE removes undecodable partial locks; no live writer may share this cache.
    vibrations.clean(empty_files=True)
    started = time.process_time_ns()
    center_energy = float(atoms.get_potential_energy())
    center_forces = atoms.get_forces()
    vibrations.run()
    for record in vibrations.cache.values():
        verify_force_record(record)
    result = {'schema': 'coherence-constrained-hessian-v1', 'inputs': manifest(input_path),
              'center_energy_eV': center_energy,
              'center_mobile_fmax_eV_A': float(np.linalg.norm(center_forces[indices], axis=1).max()),
              'active_atoms': indices, 'masses_u': inputs['masses_u'],
              'delta_A': args.delta, 'spectrum': mass_spectrum(
                  central_hessian(vibrations), inputs['masses_u']),
              'force_cache': [manifest(p) for p in
                              sorted(vibrations.cache.directory.glob('cache.*.json'))],
              'process_cpu_ns': str(time.process_time_ns() - started)}
    save(args.output / 'hessian.json', result)
    print('HESSIAN_COMPLETE', args.output.name, flush=True)


if __name__ == '__main__':
    main()
