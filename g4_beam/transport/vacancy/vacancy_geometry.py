"""Neutral fcc-Pd metal vacancies and nearby interstitial H/D geometry."""
import numpy as np
from ase import Atom
from ase.build import bulk
from ase.constraints import FixAtoms
from ase.geometry import find_mic


def make_geometry(config, case):
    a = config['lattice_A']
    atoms = bulk('Pd', 'fcc', a=a).repeat((config['repeat'],) * 3)
    missing_position = atoms.positions[0].copy()
    is_vacancy = case.startswith('vacancy')
    if is_vacancy:
        del atoms[0]  # A missing Pd atom, never an empty D site.
    n = int(case.split('_D')[1]) if '_D' in case and is_vacancy else 0
    if case == 'bulk_D':
        atoms.append(Atom('H', (a / 2, a / 2, a / 2), mass=2.01410177812))
        centre = atoms.positions[-1].copy()
    else:
        centre = missing_position
        seeds = [(1.0, 0, 0), (-1.0, 0, 0), (0, 1.0, 0),
                 (0, -1.0, 0), (0, 0, 1.0), (0, 0, -1.0)]
        if n > len(seeds):
            raise ValueError('No untested vacancy occupancy construction')
        for displacement in seeds[:n]:
            atoms.append(Atom('H', missing_position + displacement,
                              mass=2.01410177812))
    _, radii = find_mic(atoms.positions - centre, atoms.cell, pbc=True)
    frozen = [i for i, atom in enumerate(atoms)
              if atom.symbol == 'Pd' and
              (case == 'bulk' or radii[i] > config['mobile_radius_A'])]
    atoms.set_constraint(FixAtoms(indices=frozen))
    moments = [0.05 if atom.symbol == 'Pd' else 0.1 for atom in atoms]
    atoms.set_initial_magnetic_moments(moments)
    return atoms, dict(missing_Pd=bool(is_vacancy), missing_position_A=missing_position.tolist(),
                       D_count=sum(atom.symbol == 'H' for atom in atoms),
                       frozen_indices=frozen, initial_moments_muB=moments,
                       fixed_lattice=True, relaxation='D and Pd within mobile_radius; remaining Pd fixed')
