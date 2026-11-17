"""Root-only serial Native force producer from a settled SHA-bound primitive.

Usage: taskset -c N PYTHON phonon_run.py --primitive settled.traj
 --primitive-sha256 SHA --kmesh 10 5 1 --cpu N --output NEWDIR --delta .01

Input: ASE-readable single Pd3/Pd3+D primitive, explicit masses, xy PBC and
FixAtoms(0,1). --kmesh is for the REPEATED cell. No relaxation is performed.
Outputs: immutable input.json binding geometry, code and settings; ASE phonon
and withheld-x caches containing raw forces/positions/energies and bank hashes;
banks/*.gpw (compact result banks, no wavefunctions) and *.traj; scf.txt;
hessian.json with raw blocks, full signed mass_spectrum and both q sectors.
A >.05 eV/A2 withheld-row mismatch produces refusal.json, never hessian.json.

Resume only with identical input.json, one writer and no live old process.
ASE cache.strip_empties() removes undecodable/partial JSON locks, not just
zero-byte files. Completed decoded entries and their banks are checked, not
repaired. It cannot recover an in-flight SCF, absent/corrupt GPW, stale live
lock, or prove a primitive is relaxed. Root supplies a settled primitive.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
from pathlib import Path
import time

import ase
import numpy as np
from ase.io import read, write
from phonon_core import (RecordedPhonons, SUPERCELL, consume, layout, manifest,
                         run_withholding, save)


def verify_cache(cache, natoms):
    for key, record in cache.items():
        if not isinstance(record, dict):
            raise ValueError(f'Undecodable cache entry {key}; use ASE cleanup without live writers')
        forces = np.asarray(record['forces'])
        positions = np.asarray(record['positions_A'])
        if (forces.shape != (natoms, 3) or positions.shape != forces.shape
                or not np.all(np.isfinite(forces))
                or not np.all(np.isfinite(positions))
                or not np.isfinite(record['energy_eV'])):
            raise ValueError(f'Invalid completed cache record {key}')
        for name in ('gpw', 'geometry'):
            bound = record[name]
            if manifest(bound['path']) != bound:
                raise ValueError(f'Completed {key} {name} bank changed; no automatic repair')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--primitive', type=Path, required=True)
    parser.add_argument('--primitive-sha256', required=True)
    parser.add_argument('--kmesh', type=int, nargs=3, required=True)
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--delta', type=float, default=0.01)
    args = parser.parse_args()
    if os.sched_getaffinity(0) != {args.cpu}:
        raise ValueError('Declared singleton CPU affinity required')
    if (not np.isfinite(args.delta) or args.delta <= 0
            or any(k <= 0 for k in args.kmesh) or args.kmesh[2] != 1):
        raise ValueError('Positive delta and positive repeated-cell xy mesh with kz=1 required')
    source = manifest(args.primitive)
    if source['sha256'] != args.primitive_sha256:
        raise ValueError('Bound settled primitive identity changed')
    images = read(args.primitive, index=':')
    if len(images) != 1:
        raise ValueError('Bind a single settled primitive, not a relaxation trajectory')
    primitive = images[0]
    active, full = layout(primitive)
    repeated = primitive * SUPERCELL
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    from surface_calculator import calculator
    import gpaw
    from ase import phonons
    from ase.utils import filecache
    from hessian_modes import __file__ as modes_file
    from surface_calculator import __file__ as calc_file
    from surface_io import geometry
    inputs = {'schema': 'coherence-two-cell-phonon-input-v1',
              'primitive': source, 'primitive_geometry': geometry(primitive),
              'primitive_masses_u': primitive.get_masses().tolist(),
              'active_primitive': active, 'active_supercell': full,
              'supercell': list(SUPERCELL), 'delta_A': args.delta,
              'kmesh_repeated_cell': args.kmesh, 'cpu': args.cpu,
              'electrons_repeated_cell': sum(10 if s == 'Pd' else 1 for s in repeated.symbols),
              'cutoff_eV': 300, 'xc': 'PBE', 'Pd_setup': '10', 'FD_width_eV': .05,
              'gamma_centered': True, 'dipolelayer': True,
              'legacy_gpaw': False, 'symmetry': {'point_group': False, 'time_reversal': True},
              'convergence': {'energy': 1e-7, 'density': 1e-6, 'eigenstates': 1e-9},
              'ase_version': ase.__version__, 'gpaw_version': gpaw.__version__,
              'code': [manifest(p) for p in (Path(__file__),
                       Path(__file__).with_name('phonon_core.py'), modes_file,
                       calc_file, phonons.__file__, filecache.__file__)]}
    input_path = args.output / 'input.json'
    if input_path.exists():
        if json.loads(input_path.read_text()) != inputs:
            raise ValueError('Existing output is bound to different inputs/code')
        if (args.output / 'hessian.json').exists() or (args.output / 'refusal.json').exists():
            raise ValueError('Finished output is immutable; use a new output directory')
    elif any(args.output.iterdir()):
        raise ValueError('Unbound nonempty output directory refused')
    else:
        save(input_path, inputs)
    banks = args.output / 'banks'
    banks.mkdir(exist_ok=True)
    # Crucially size bands/electrons for the repeated cell, then clone the
    # established slab calculator with force-qualified symmetry and SCF.
    base = calculator(repeated, args.output / 'scf.txt', kmesh=tuple(args.kmesh))
    calc = base.new(legacy_gpaw=False, symmetry=inputs['symmetry'],
                    convergence=inputs['convergence'])
    del base

    def record(atoms, name):
        gpw_path, geometry_path = banks / f'{name}.gpw', banks / f'{name}.traj'
        atoms.calc.write(str(gpw_path), mode='')
        write(geometry_path, atoms)
        return {'gpw': manifest(gpw_path), 'geometry': manifest(geometry_path)}

    ph = RecordedPhonons(primitive, calc, supercell=SUPERCELL, delta=args.delta,
                         name=str(args.output / 'cache'), recorder=record)
    ph.set_atoms(active)
    started = time.process_time()
    try:
        removed = ph.cache.strip_empties()
        verify_cache(ph.cache, len(repeated))
        ph.run()
        withheld, removed_withheld = run_withholding(ph, str(args.output / 'withheld-cache'))
        verify_cache(ph.cache, len(repeated))
        verify_cache(withheld, len(repeated))
        result = consume(ph, withheld)
        result.update(schema='coherence-two-cell-phonon-result-v1',
                      inputs=manifest(input_path),
                      center_energy_eV=float(ph.cache['eq']['energy_eV']),
                      center_mobile_fmax_eV_A=float(np.linalg.norm(
                          ph.cache['eq']['forces'][full], axis=1).max()),
                      process_cpu_s=time.process_time() - started,
                      ase_removed_undecodable_entries=removed + removed_withheld,
                      force_cache=[manifest(p) for directory in (ph.cache.directory, withheld.directory)
                                   for p in sorted(directory.glob('cache.*.json'))])
        save(args.output / ('hessian.json' if result['accepted'] else 'refusal.json'), result)
        if not result['accepted']:
            raise RuntimeError('Translated-cell x-row discrepancy exceeds .05 eV/A2; spectrum refused')
    except Exception as exc:
        save(args.output / 'failure.json', {'type': type(exc).__name__, 'message': str(exc),
             'process_cpu_s': time.process_time() - started,
             'recovery': 'No in-flight SCF recovery. Resume only identical bound input with no live writer.'})
        raise
    print('TWO_CELL_PHONON_COMPLETE', args.output.name, flush=True)


if __name__ == '__main__':
    main()
