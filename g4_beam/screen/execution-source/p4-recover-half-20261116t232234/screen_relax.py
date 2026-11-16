"""Explicit bulk constraint/force/stress convergence, separate from spectroscopy."""
import numpy as np
from ase.filters import UnitCellFilter
from ase.optimize import BFGS

from screen_io import require


def relax_bulk(atoms, settings, output):
    mode = settings['mode']
    if mode == 'fixed':
        return dict(mode=mode, steps=0, relaxed=False)
    masks = {'hydrostatic': [1, 1, 1, 0, 0, 0],
             'hexagonal': [1, 1, 1, 0, 0, 0],
             'cell': [1, 1, 1, 1, 1, 1]}
    require(mode in masks, 'Unsupported commissioned bulk relaxation')
    mask = np.asarray(masks[mode], bool)
    filtered = UnitCellFilter(atoms, mask=mask,
                              hydrostatic_strain=(mode == 'hydrostatic'))
    threshold = min(settings['force_tolerance_eV_A'],
                    .4*settings['stress_tolerance_eV_A3']*atoms.get_volume()/len(atoms))
    optimizer = BFGS(filtered, maxstep=.05,
                     logfile=str(output/'relax.log'),
                     trajectory=str(output/'relax.traj'))
    converged = optimizer.run(fmax=threshold, steps=settings['max_steps'])
    forces = atoms.get_forces()
    stress = atoms.get_stress()
    require(converged, 'Bulk relaxation reached its admitted step cap')
    require(float(np.linalg.norm(forces, axis=1).max()) <= settings['force_tolerance_eV_A'],
            'Relaxed internal forces exceed the declared tolerance')
    require(float(np.max(np.abs(stress[mask]))) <= settings['stress_tolerance_eV_A3'],
            'Relaxed active stress exceeds the declared tolerance')
    return dict(mode=mode, steps=optimizer.nsteps, relaxed=True,
                settings=settings, forces_eV_A=forces.tolist(),
                stress_eV_A3=stress.tolist(),
                interpretation='Stationary ordered phase; no global phase or phonon stability claim')
