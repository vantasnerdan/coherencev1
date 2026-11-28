"""Index-two PdD supercell and the exact folded primitive electronic grid."""
import numpy as np
from ase.build import make_supercell

from native import A, reference, require

SUPERCELL = np.array([[-1, 1, 1], [1, 0, 0], [0, 1, -1]])
POLARIZATIONS = {'100': np.array([1., 0., 0.]), '010': np.array([0., 1., 0.])}


def mesh(primitive_n):
    require(primitive_n in (12, 16), 'Unregistered primitive sampling')
    # An even Gamma-centered Monkhorst–Pack mesh is the group Z_n^3/n.
    numerator = np.indices((primitive_n,) * 3).reshape(3, -1).T
    mapped = (numerator @ SUPERCELL.T) % primitive_n
    folded, count = np.unique(mapped, axis=0, return_counts=True)
    require(len(folded) == primitive_n**3 // 2 and np.all(count == 2),
            'Folded k points do not have equal two-to-one preimages')
    return ((folded + primitive_n // 2) % primitive_n - primitive_n // 2) / primitive_n


def geometry(mode, polarization):
    require(mode in ('gamma', 'x') and polarization in POLARIZATIONS, 'Unknown PdD mode or polarization')
    atoms = make_supercell(reference(), SUPERCELL, order='atom-major')
    sites = np.flatnonzero(atoms.numbers == 1)
    require(len(atoms) == 4 and sites.tolist() == [2, 3], 'PdD cell does not have exactly two ordered D sites')
    phase_fraction = (atoms.positions[sites, 0] - atoms.positions[sites[0], 0]) / A
    phase = np.cos(2 * np.pi * phase_fraction)
    require(np.allclose(phase, [1, -1], atol=1e-10, rtol=0)
            and np.allclose(atoms.cell.array[:, 0] / A, (atoms.cell.array[:, 0] / A).round(), atol=1e-10, rtol=0),
            'X does not have two opposite periodic site phases')
    signs = np.array([1., -1.]) if mode == 'x' else np.ones(2)
    return atoms, sites, signs, POLARIZATIONS[polarization]
