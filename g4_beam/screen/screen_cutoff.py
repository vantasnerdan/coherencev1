"""Restrict chi0's PW basis, then invert; never restrict the full inverse."""
import numpy as np
from ase.units import Bohr, Hartree

from screen_io import require


def restricted_inverse(Q_Ainv, epsilon_symmetrized, cutoff_eV):
    vectors = np.asarray(Q_Ainv)
    epsilon = np.asarray(epsilon_symmetrized)
    require(vectors.ndim == 2 and vectors.shape[1] == 3,
            'Physical reciprocal vectors must have shape (nG,3)')
    require(epsilon.shape == (len(vectors), len(vectors)),
            'Dielectric and reciprocal bases differ')
    require(cutoff_eV > 0 and np.isfinite(epsilon).all(),
            'Invalid cutoff or dielectric matrix')
    kinetic_eV = 0.5*Hartree*Bohr**2*np.einsum('gv,gv->g', vectors, vectors)
    indices = np.flatnonzero(kinetic_eV < cutoff_eV)
    require(len(indices) > 0 and indices[0] == 0,
            'Restricted PW basis does not retain the physical head')
    selected = epsilon[np.ix_(indices, indices)]
    inverse = np.linalg.inv(selected)
    require(np.linalg.norm(selected@inverse-np.eye(len(indices)), ord=np.inf) < 1e-8,
            'Restricted dielectric inverse failed')
    return indices, inverse
