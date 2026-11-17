"""Explicit ordered Pd(111) slab and H/D electronic geometry."""
import numpy as np
from ase import Atoms
from ase.build import fcc111, add_adsorbate
from ase.constraints import FixAtoms


def slab(hydrogens=2, layers=3, vacuum=8.0, size=(1, 2)):
    atoms = fcc111('Pd', size=(*size, layers), a=3.890, vacuum=vacuum)
    atoms.pbc = (True, True, False)
    # Only outermost metal layer relaxes; all lower layers remain crystalline.
    atoms.set_constraint(FixAtoms(mask=atoms.get_tags() != 1))
    if hydrogens:
        add_adsorbate(atoms, 'H', 1.0, position='fcc')
    if hydrogens == 2:
        add_adsorbate(atoms, 'H', 1.0, position='fcc', offset=(0, 1))
    atoms.set_masses([2.014101778 if s == 'H' else m
                      for s, m in zip(atoms.symbols, atoms.get_masses())])
    return atoms


def gas(symbol):
    atoms = Atoms(symbol, cell=(10., 10., 10.), pbc=True)
    if len(atoms) == 2:
        atoms.positions[:] = [[5, 5, 4.625], [5, 5, 5.375]]
    else:
        atoms.positions[:] = [[5, 5, 5]]
        atoms.set_initial_magnetic_moments([1.0])
    return atoms
