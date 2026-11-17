"""Relax the explicitly bound geometry/constraints at one declared k mesh."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import hashlib
from pathlib import Path
from ase.io import read
from surface_calculator import calculator
from surface_optimize import optimize
from surface_io import geometry, manifest, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--input-sha256', required=True)
    parser.add_argument('--kmesh', type=int, nargs=3, required=True)
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('--steps', type=int, default=20)
    parser.add_argument('--alpha', type=float, default=70.0)
    args = parser.parse_args()
    if os.sched_getaffinity(0) != {args.cpu}:
        raise ValueError('Actual singleton CPU differs from assignment')
    if hashlib.sha256(args.input.read_bytes()).hexdigest() != args.input_sha256:
        raise ValueError('Prepared geometry binding changed')
    atoms = read(args.input)
    atoms.set_masses([2.014101778 if symbol == 'H' else mass
                      for symbol, mass in zip(atoms.symbols, atoms.get_masses())])
    args.output.mkdir(parents=True, exist_ok=False)
    save(args.output / 'input.json', {
        'schema': 'coherence-surface-control-v1', 'label': args.label,
        'source_geometry': manifest(args.input), 'geometry': geometry(atoms),
        'masses_u': atoms.get_masses().tolist(), 'kmesh': args.kmesh,
        'cutoff_eV': 300, 'fmax_target_eV_A': .05, 'steps_cap': args.steps,
        'BFGS_initial_alpha_eV_A2': args.alpha,
        'symmetry': {'point_group': False, 'time_reversal': True},
        'affinity': sorted(os.sched_getaffinity(0)),
        'code': [manifest(Path(__file__).with_name(name)) for name in
                 ('run_surface_control.py', 'surface_optimize.py',
                  'surface_calculator.py', 'surface_io.py')],
        'claim': 'Minimum only in the explicitly supplied mobile subspace; '
                 'a constrained chart minimum is not a full saddle or unconstrained minimum.',
    })
    calc = calculator(atoms, args.output / 'scf.txt', kmesh=tuple(args.kmesh))
    calc = calc.new(symmetry={'point_group': False, 'time_reversal': True},
                    legacy_gpaw=False, txt=str(args.output / 'scf.txt'))
    print('SURFACE_CONTROL_START', args.label, flush=True)
    optimize(atoms, args.output / 'relax', steps=args.steps,
             restart_calc=calc, fmax=.05, alpha=args.alpha)
    print('SURFACE_CONTROL_COMPLETE', args.label, flush=True)


if __name__ == '__main__':
    main()
