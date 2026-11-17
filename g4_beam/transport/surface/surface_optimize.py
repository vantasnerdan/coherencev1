"""Bounded endpoint optimization with actual force and saved-state records."""
import time
import numpy as np
from ase.io import write
from ase.optimize import BFGS
from surface_calculator import calculator
from surface_io import save, geometry, manifest


def optimize(atoms, out, steps=35, restart_calc=None, fmax=0.05):
    out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    atoms.calc = restart_calc or calculator(atoms, out/'scf.txt')
    opt = BFGS(atoms, logfile=str(out/'relax.log'), trajectory=str(out/'relax.traj'), maxstep=0.15)
    points = []
    def checkpoint():
        e = float(atoms.get_potential_energy())
        f = atoms.get_forces()
        p = out/f'point-{opt.nsteps:03d}.json'
        save(p, {'energy_extrapolated_eV': e,
                 'free_energy_eV': float(atoms.get_potential_energy(force_consistent=True)),
                 'forces_eV_A': f.tolist(),
                 'raw_forces_eV_A': atoms.get_forces(apply_constraint=False).tolist(),
                 'fmax_eV_A': float(np.linalg.norm(f, axis=1).max()),
                 'spin_polarized': atoms.calc.get_number_of_spins() == 2,
                 'geometry': geometry(atoms), 'ionic_step': opt.nsteps,
                 'wall_s': time.monotonic()-started})
        atoms.calc.write(str(out/'checkpoint.gpw'), mode='all')
        points.append(manifest(p))
        print('ENDPOINT_FORCE_POINT', out.name, opt.nsteps, e, flush=True)
    opt.attach(checkpoint)
    converged = opt.run(fmax=fmax, steps=steps)
    write(out/'final.traj', atoms)
    save(out/'result.json', {'converged': bool(converged), 'fmax_target_eV_A': fmax,
                           'wall_s': time.monotonic()-started, 'points': points,
                           'final': points[-1], 'checkpoint': manifest(out/'checkpoint.gpw')})
    if not converged:
        raise RuntimeError('Endpoint exceeded fixed ionic-step budget; no automatic scientific retry')
    return atoms
