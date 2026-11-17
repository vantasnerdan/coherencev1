"""Raw central force derivatives and their actual atomic-mass spectrum."""
import numpy as np
from ase import units


def central_hessian(vibrations):
    """Read ASE's completed decoded cache, retaining reciprocity before averaging."""
    if vibrations.nfree != 2:
        raise ValueError('This producer declares two-point central displacements')
    indices = np.asarray(vibrations.indices, dtype=int)
    pairs = {}
    for disp in vibrations.displacements():
        if disp.sign == 0:
            continue
        pairs.setdefault((disp.a, disp.i), {})[disp.sign] = disp.forces()
    rows = []
    for atom in indices:
        for axis in range(3):
            pair = pairs[(atom, axis)]
            minus, plus = (np.asarray(pair[sign], dtype=float) for sign in (-1, 1))
            if minus.shape != (len(vibrations.atoms), 3) or plus.shape != minus.shape:
                raise ValueError('Force cache and bound geometry differ')
            rows.append(((minus - plus)[indices] / (2 * vibrations.delta)).ravel())
    raw = np.asarray(rows)
    if not np.all(np.isfinite(raw)):
        raise ValueError('Nonfinite force derivative')
    return raw


def mass_spectrum(raw_hessian, atomic_masses_u):
    """Preserve signed quanta and full modes; do not clip soft/unstable directions."""
    raw = np.asarray(raw_hessian, dtype=float)
    masses = np.asarray(atomic_masses_u, dtype=float)
    if (raw.shape != (3 * len(masses),) * 2 or not len(masses)
            or np.any(~np.isfinite(raw)) or np.any(~np.isfinite(masses))
            or np.any(masses <= 0)):
        raise ValueError('Finite full Cartesian Hessian and positive isotope masses required')
    symmetric = (raw + raw.T) / 2
    inverse_sqrt_mass = np.repeat(masses ** -0.5, 3)
    weighted = inverse_sqrt_mass[:, None] * symmetric * inverse_sqrt_mass
    eigenvalues, columns = np.linalg.eigh(weighted)
    scale = units._hbar * 1e10 / np.sqrt(units._e * units._amu)
    signed_quanta = np.sign(eigenvalues) * scale * np.sqrt(np.abs(eigenvalues))
    cartesian_rows = (inverse_sqrt_mass[:, None] * columns).T
    raw_norm = float(np.linalg.norm(raw))
    antisymmetric_norm = float(np.linalg.norm(raw - raw.T))
    return {'raw_H_eV_A2': raw.tolist(), 'symmetric_H_eV_A2': symmetric.tolist(),
            'raw_reciprocity_defect_frobenius_eV_A2': antisymmetric_norm,
            'raw_reciprocity_relative_defect': antisymmetric_norm / raw_norm if raw_norm else None,
            'mass_diagonal_u': np.repeat(masses, 3).tolist(),
            'eigenvalues_eV_A2_u': eigenvalues.tolist(),
            'signed_quanta_eV': signed_quanta.tolist(),
            'massweighted_eigenvectors_columns': columns.tolist(),
            'cartesian_mass_normalized_modes_rows': cartesian_rows.tolist()}
