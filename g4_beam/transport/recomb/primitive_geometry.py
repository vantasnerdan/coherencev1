"""Project only measured sub-milliangstrom two-cell translation residuals."""
import numpy as np
from ase.constraints import FixAtoms
from ase.geometry import find_mic


def primitive_from_two_cell(atoms, tolerance_A=.001):
    symbols = list(atoms.symbols)
    npd = symbols.count('Pd')
    if npd not in (6, 8) or symbols not in (
            ['Pd'] * npd, ['Pd'] * npd + ['H'] * 2):
        raise ValueError('Declared uniformly occupied clean/1ML two-cell ordering required')
    if tuple(atoms.pbc) != (True, True, False):
        raise ValueError('The declared two-dimensional slab is required')
    translation = np.asarray(atoms.cell[1]) / 2
    positions, residuals = [], []
    masses = atoms.get_masses()
    for first in range(0, len(atoms), 2):
        if masses[first] != masses[first + 1]:
            raise ValueError('Translation-equivalent masses differ')
        residual, distance = find_mic(
            atoms.positions[first + 1] - translation - atoms.positions[first],
            atoms.cell, pbc=atoms.pbc)
        if distance > tolerance_A:
            raise ValueError('Broken-symmetry geometry cannot be silently projected')
        positions.append(atoms.positions[first] + residual / 2)
        residuals.append(float(distance))
    primitive = atoms[::2]
    cell = atoms.cell.copy()
    cell[1] /= 2
    primitive.set_cell(cell, scale_atoms=False)
    primitive.positions[:] = positions
    primitive.set_constraint(FixAtoms(indices=[0, 1]))
    primitive.set_masses([2.014101778 if s == 'H' else m
                          for s, m in zip(primitive.symbols, primitive.get_masses())])
    primitive.calc = None
    primitive.info.clear()
    return primitive, {'translation_residuals_A': residuals,
                       'projection_gate_A': tolerance_A,
                       'full_cell_stability': 'Not assumed; both q sectors must be tested'}
