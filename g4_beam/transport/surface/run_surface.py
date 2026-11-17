"""One bounded real feasibility pilot; expansion requires measured decision."""
import argparse
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import time
from pathlib import Path
import numpy as np
from ase.io import write
from ase.optimize import FIRE
import gpaw
from surface_calculator import calculator
from surface_geometry import slab
from surface_io import save, geometry, manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    atoms = slab()
    save(out/'initial.json', geometry(atoms))
    atoms.calc = calculator(atoms, out/'scf.txt')
    opt = FIRE(atoms, logfile=str(out/'relax.log'), trajectory=str(out/'relax.traj'))
    started = time.monotonic()
    def point():
        energy = float(atoms.get_potential_energy())
        free = float(atoms.get_potential_energy(force_consistent=True))
        forces = atoms.get_forces()
        raw = atoms.get_forces(apply_constraint=False)
        record = {'ionic_step': opt.nsteps, 'wall_s': time.monotonic()-started,
                  'energy_extrapolated_eV': energy, 'free_energy_eV': free,
                  'forces_eV_A': forces.tolist(), 'raw_forces_eV_A': raw.tolist(),
                  'fmax_eV_A': float(np.linalg.norm(forces, axis=1).max()),
                  'geometry': geometry(atoms)}
        save(out/f'point-{opt.nsteps:03d}.json', record)
        atoms.calc.write(str(out/f'point-{opt.nsteps:03d}.gpw'), mode='all')
        print('PILOT_FORCE_POINT', opt.nsteps, record['wall_s'], energy,
              record['fmax_eV_A'], flush=True)
    opt.attach(point, interval=1)
    converged = opt.run(fmax=0.05, steps=2)
    save(out/'pilot.json', {'schema': 'p4-surface-pilot-v1',
         'gpaw_version': gpaw.__version__, 'converged': bool(converged),
         'wall_s': time.monotonic()-started, 'points':
         [manifest(p) for p in sorted(out.glob('point-*.json'))],
         'interpretation': 'Cost pilot, not a converged barrier or adsorption energy'})
    print('PILOT_COMPLETED', flush=True)

if __name__ == '__main__':
    main()
