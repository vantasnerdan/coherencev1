"""Exact wrapped-midpoint S3/inversion classes for a verified cubic site.

The caller must also establish the nonmagnetic ground state's point-group
symmetry.  This module does not infer a dielectric matrix's covariance.
"""
from __future__ import annotations

from itertools import permutations, product

import numpy as np


PERMUTATIONS = tuple(permutations(range(3)))


def midpoint_classes(n: int) -> list[dict]:
    """Partition every midpoint once using exact integer coordinates.

    Fractional q has denominator2n.  Class multiplicity divided by n**3 is
    its normalized BZ weight.  No optical Gamma point is generated.
    """
    if isinstance(n, bool) or not isinstance(n, int) or n < 2:
        raise ValueError('The registered midpoint grid requires integer n>=2')
    wrap = lambda value: (value + n) % (2 * n) - n
    coordinates = tuple(wrap(2 * j + 1) for j in range(n))
    unassigned = set(product(coordinates, repeat=3))
    classes = []
    while unassigned:
        representative = min(unassigned)
        orbit = {tuple(wrap(sign * representative[i]) for i in permutation)
                 for sign in (-1, 1) for permutation in PERMUTATIONS}
        if not orbit <= unassigned:
            raise ArithmeticError('Midpoint orbits overlap or leave the registered grid')
        unassigned.difference_update(orbit)
        classes.append(dict(numerators=list(representative), denominator=2 * n,
                            q_fractional=[value / (2 * n) for value in representative],
                            multiplicity=len(orbit), weight=len(orbit) / n**3,
                            members_numerators=[list(member) for member in sorted(orbit)]))
    return classes


def validate_site_group(cell_A, positions_A, symbols, source_A, kmesh):
    """Refuse S3/inversion reduction unless the actual cell, ions and site fit.

    Primitive-axis permutations must be the same Cartesian orthogonal
    permutations: cell commutes with them.  Positions and source may move
    by direct-lattice vectors, not to a different inequivalent site.
    """
    cell = np.asarray(cell_A, float)
    positions = np.asarray(positions_A, float)
    source = np.asarray(source_A, float)
    mesh = np.asarray(kmesh)
    if (cell.shape != (3, 3) or positions.shape != (len(symbols), 3) or
            source.shape != (3,) or not np.isfinite(cell).all() or
            not np.isfinite(positions).all() or not np.isfinite(source).all()):
        raise ValueError('Malformed physical site geometry')
    if mesh.shape != (3,) or not np.all(mesh == mesh[0]) or mesh[0] <= 0:
        raise ValueError('The site group requires an equal-axis uniform k mesh')
    inverse_cell = np.linalg.inv(cell)
    operations = []
    for sign, permutation in product((-1, 1), PERMUTATIONS):
        operation = sign * np.eye(3, dtype=int)[list(permutation)]
        if not np.allclose(operation @ cell, cell @ operation, rtol=0, atol=1e-9):
            raise ValueError('Cell does not admit the proposed orthogonal site group')
        delta = (source @ operation - source) @ inverse_cell
        if not np.allclose(delta, np.rint(delta), rtol=0, atol=1e-9):
            raise ValueError('Source is not invariant modulo the direct lattice')
        for symbol, position in zip(symbols, positions):
            candidates = positions[np.asarray(symbols) == symbol]
            differences = (position @ operation - candidates) @ inverse_cell
            if not np.any(np.max(np.abs(differences - np.rint(differences)), axis=1) < 1e-9):
                raise ValueError('Ionic species do not admit the proposed site group')
        operations.append(operation.tolist())
    return operations
