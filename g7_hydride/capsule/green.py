"""Qualified zero-G Ewald oracle, with stdlib erfc instead of scipy.special.

Constants are the actual production scipy.constants Coulomb value and ASE Bohr;
no inference of either constant from output curves is used.
"""
import math
import numpy as np

COULOMB_EV_A = 14.399645468667815
BOHR_A = 0.5291772105638411


def triplets(shell):
    axis = np.arange(-shell, shell + 1)
    return np.stack(np.meshgrid(axis, axis, axis, indexing='ij'), axis=-1).reshape(-1, 3)


def erfc(values):
    return np.fromiter((math.erfc(float(x)) for x in values), dtype=float, count=len(values))


class PeriodicNuclearGreen:
    def __init__(self, cell_A):
        self.cell = np.asarray(cell_A, dtype=float)
        if self.cell.shape != (3, 3) or not np.isfinite(self.cell).all():
            raise ValueError('Finite 3x3 cell required')
        self.volume = float(np.linalg.det(self.cell))
        if self.volume <= 0:
            raise ValueError('Positive cell volume required')
        self.eta = 2 / self.volume ** (1 / 3)
        self.shifts = triplets(4) @ self.cell
        reciprocal = 2 * np.pi * np.linalg.inv(self.cell).T
        modes = triplets(6)
        self.G = modes[np.any(modes != 0, axis=1)] @ reciprocal
        g2 = np.einsum('ij,ij->i', self.G, self.G)
        self.weights = 4 * np.pi / self.volume * np.exp(-g2 / (4 * self.eta ** 2)) / g2

    def distances(self, displacement):
        point = np.asarray(displacement, dtype=float)
        if point.shape != (3,) or not np.isfinite(point).all():
            raise ValueError('Finite target-relative vector required')
        distances = np.linalg.norm(self.shifts + point, axis=1)
        if np.min(distances) < 1e-12:
            raise ValueError('Point-nuclear coalescence is singular')
        return point, distances

    def point_eV(self, displacement):
        point, distances = self.distances(displacement)
        real = np.sum(erfc(self.eta * distances) / distances)
        reciprocal = np.dot(self.weights, np.cos(self.G @ point))
        return float(COULOMB_EV_A * (real + reciprocal - np.pi / (self.eta ** 2 * self.volume)))

    def point_minus_gaussian_eV(self, displacement, width):
        _, distances = self.distances(displacement)
        if not math.isfinite(width) or width <= 0:
            raise ValueError('Positive Gaussian pair width required')
        return float(COULOMB_EV_A * (np.sum(erfc(width * distances) / distances)
                                    - np.pi / (width ** 2 * self.volume)))
