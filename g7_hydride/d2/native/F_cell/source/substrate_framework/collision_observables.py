"""Periodic electron-only potentials and prescribed-trajectory work (issue 244).

Densities are electrons/A^3; potentials are energy/e for a positive unit probe,
in eV.  The omitted G=0 mode defines a zero-cell-mean electron potential.
"""
from __future__ import annotations

import numpy as np
from scipy.constants import elementary_charge, epsilon_0, physical_constants

COULOMB_EV_A = elementary_charge / (4.0 * np.pi * epsilon_0) * 1e10
DEUTERON_KG = physical_constants['deuteron mass'][0]


class PeriodicElectronPotential:
    """Spectral Coulomb readout on the producer's periodic real-space grid."""

    def __init__(self, shape, cell_A):
        self.shape = tuple(int(n) for n in shape)
        self.cell_A = np.asarray(cell_A, dtype=float)
        if len(self.shape) != 3 or min(self.shape) < 3:
            raise ValueError('A three-dimensional periodic density grid is required')
        if self.cell_A.shape != (3, 3) or np.linalg.det(self.cell_A) <= 0:
            raise ValueError('A right-handed nonsingular cell is required')
        reciprocal = 2.0 * np.pi * np.linalg.inv(self.cell_A).T
        axes = [np.fft.fftfreq(n) * n for n in self.shape]
        indices = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1)
        self.vectors = (indices @ reciprocal).reshape(-1, 3)
        q2 = np.einsum('ij,ij->i', self.vectors, self.vectors)
        self.kernel = np.zeros(q2.shape)
        nonzero = q2 > 0
        self.kernel[nonzero] = -4.0 * np.pi * COULOMB_EV_A / q2[nonzero]
        self.voxel_A3 = np.linalg.det(self.cell_A) / np.prod(self.shape)

    def samples(self, density_e_A3, points_A):
        density = np.asarray(density_e_A3, dtype=float)
        points = np.asarray(points_A, dtype=float).reshape(-1, 3)
        if density.shape != self.shape or not np.isfinite(density).all():
            raise ValueError('Density does not match its finite producer grid')
        coefficients = np.fft.fftn(density).ravel() / density.size
        coefficients *= self.kernel
        # Real trigonometric interpolation also specifies off-grid Nyquist modes.
        return np.array([(coefficients * np.exp(1j * (self.vectors @ p))).sum().real
                         for p in points])

    def electrons(self, density_e_A3):
        return float(np.sum(density_e_A3) * self.voxel_A3)


def deuteron_speed_A_as(energy_keV):
    """Nonrelativistic lab speed from the actual deuteron mass."""
    if not np.isfinite(energy_keV) or energy_keV <= 0:
        raise ValueError('Positive finite laboratory energy is required')
    velocity_m_s = np.sqrt(2.0 * energy_keV * 1e3 * elementary_charge / DEUTERON_KG)
    return float(velocity_m_s * 1e-8)


def external_work_eV(positions_A, forces_eV_A):
    """Cumulative external work: minus the trapezoidal projectile force work."""
    positions = np.asarray(positions_A, dtype=float)
    forces = np.asarray(forces_eV_A, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3 or forces.shape != positions.shape:
        raise ValueError('Matched Cartesian projectile positions and forces are required')
    work = -np.einsum('ij,ij->i', 0.5 * (forces[1:] + forces[:-1]), np.diff(positions, axis=0))
    return np.concatenate(([0.0], np.cumsum(work)))


def excitation_and_closure(td_eV, bo_eV, work_eV, initial_excitation_eV=0.0):
    """Retain signed excitation and residual; never clip or rescale a defect."""
    td = np.asarray(td_eV, dtype=float)
    bo = np.asarray(bo_eV, dtype=float)
    work = np.asarray(work_eV, dtype=float)
    if td.shape != bo.shape or td.shape != work.shape or td.size == 0:
        raise ValueError('Matched nonempty TD, BO and work checkpoint arrays are required')
    excitation = td - bo
    residual = work - ((bo - bo[0]) + excitation - initial_excitation_eV)
    return excitation, residual


def extrapolate_contact(separations_A, screening_eV):
    """Report linear/quadratic contact estimates and fit-window sensitivity."""
    r = np.asarray(separations_A, dtype=float)
    u = np.asarray(screening_eV, dtype=float)
    if r.ndim != 1 or u.shape != r.shape or len(r) < 4:
        raise ValueError('At least four distinct short-range samples are required')
    if not np.isfinite(r).all() or not np.isfinite(u).all() or len(np.unique(r)) != len(r):
        raise ValueError('Finite distinct separations are required')
    order = np.argsort(r)
    values = []
    for count in (3, 4):
        selected = order[:count]
        for degree in (1, 2):
            fit = np.polynomial.polynomial.polyfit(r[selected], u[selected], degree)
            values.append(dict(points=count, degree=degree, contact_eV=float(fit[0])))
    contacts = [v['contact_eV'] for v in values]
    return dict(primary_quadratic_three_eV=values[1]['contact_eV'],
                minimum_eV=min(contacts), maximum_eV=max(contacts),
                spread_eV=max(contacts) - min(contacts), fits=values,
                smallest_resolved_separation_A=float(r.min()))
