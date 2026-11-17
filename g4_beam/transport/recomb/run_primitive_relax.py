"""Strict primitive stationary point; full two-cell stability remains a separate test."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
import hashlib
from pathlib import Path
from ase.io import read
from primitive_geometry import primitive_from_two_cell
from surface_calculator import calculator
from surface_optimize import optimize
from surface_io import geometry, manifest, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--source-sha256', required=True)
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--label', required=True)
    args = parser.parse_args()
    if os.sched_getaffinity(0) != {args.cpu}:
        raise ValueError('Actual singleton CPU differs from the assigned core')
    if hashlib.sha256(args.source.read_bytes()).hexdigest() != args.source_sha256:
        raise ValueError('Two-cell stationary-point source identity changed')
    primitive, projection = primitive_from_two_cell(read(args.source))
    args.output.mkdir(parents=True, exist_ok=False)
    convergence = {'energy': 1e-7, 'density': 1e-6, 'eigenstates': 1e-9}
    save(args.output / 'input.json', {
        'schema': 'coherence-strict-primitive-stationary-v1', 'label': args.label,
        'source': manifest(args.source), 'geometry': geometry(primitive),
        'projection': projection, 'masses_u': primitive.get_masses().tolist(),
        'kmesh': [10, 10, 1], 'convergence': convergence,
        'symmetry': {'point_group': False, 'time_reversal': True},
        'fmax_target_eV_A': .01, 'steps_cap': 20,
        'code': [manifest(Path(__file__).with_name(name)) for name in
                 ('run_primitive_relax.py', 'primitive_geometry.py',
                  'surface_calculator.py', 'surface_optimize.py', 'surface_io.py')],
        'claim': 'Primitive stationary point only; no assumption that '
                 'the full two-cell zone-boundary sector is stable.',
    })
    calc = calculator(primitive, args.output / 'scf.txt', kmesh=(10, 10, 1))
    calc = calc.new(legacy_gpaw=False, convergence=convergence,
        symmetry={'point_group': False, 'time_reversal': True},
        txt=str(args.output / 'scf.txt'))
    print('STRICT_PRIMITIVE_START', args.label, flush=True)
    optimize(primitive, args.output / 'relax', steps=20,
             restart_calc=calc, fmax=.01)
    print('STRICT_PRIMITIVE_COMPLETE', args.label, flush=True)


if __name__ == '__main__':
    main()
