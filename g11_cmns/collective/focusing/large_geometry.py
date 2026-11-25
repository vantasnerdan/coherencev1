"""Dense-index lift of the frozen focusing lattice; dynamics stay in engine."""
from functools import cached_property
from itertools import product

import numpy as np
from numba import njit

import engine as E


@njit(cache=True)
def _edge_arrays(grid, side, offsets, strengths, near_indices, near_bonds, near_axes):
    n = len(grid)
    lookup = np.full((side, side, side), -1, dtype=np.int32)
    for i in range(n):
        gx, gy, gz = grid[i]
        lookup[gx, gy, gz] = i

    # For m >= 4 each offset has a distinct residue; d and -d enumerate
    # opposite directions of the same edge. Keep the old i, then stencil order.
    count = n * len(offsets) // 2
    near_count = n * len(near_bonds) // 2
    edges = np.empty((count, 2), dtype=np.int32)
    weights = np.empty(count)
    near = np.empty((near_count, 2), dtype=np.int32)
    bonds = np.empty((near_count, 3))
    directions = np.empty((near_count, 3))
    b, k = 0, 0
    for i in range(n):
        gx, gy, gz = grid[i]
        for s in range(len(offsets)):
            dx, dy, dz = offsets[s]
            j = lookup[(gx + dx) % side, (gy + dy) % side, (gz + dz) % side]
            if j <= i:
                continue
            edges[b, 0], edges[b, 1], weights[b] = i, j, strengths[s]
            b += 1
            ns = near_indices[s]
            if ns >= 0:
                near[k, 0], near[k, 1] = i, j
                for a in range(3):
                    bonds[k, a] = near_bonds[ns, a]
                    directions[k, a] = near_axes[ns, a]
                k += 1
    if b != count or k != near_count:
        raise ValueError('incomplete periodic stencil')
    return edges, weights, near, bonds, directions


class LargeLattice(E.Lattice):
    """Same register order, Hamiltonian and trajectory inputs through 131072 D."""

    def __init__(self, m):
        if m == 2:
            # Distinct shortest images can alias here; retain exact refolding.
            super().__init__(m)
            return
        if m not in (4, 8, 16, 32):
            raise ValueError('conventional cell count must be 2,4,8,16,32')
        self.m, self.box = m, m * E.CELL / 2
        side = 2 * m
        tiles = np.array(list(product(range(m // 2), repeat=3)), dtype=np.int32)
        base = E.BASE['positions']
        self.positions = (base[None, :, :] + E.CELL * tiles[:, None, :]).reshape(-1, 3)
        base_grid = np.rint(base / E.HALF).astype(np.int32)
        self.grid = (base_grid[None, :, :] + 4 * tiles[:, None, :]).reshape(-1, 3)
        self.n = len(self.positions)
        self.stencil = E.stencil()
        offsets = np.array([d for d, _ in self.stencil], dtype=np.int32)
        strengths = np.array([w for _, w in self.stencil])
        nearest = np.einsum('ij,ij->i', offsets, offsets) == 2
        near_indices = np.full(len(offsets), -1, dtype=np.int32)
        near_bonds = offsets[nearest] * E.HALF
        near_indices[nearest] = np.arange(len(near_bonds), dtype=np.int32)
        near_axes = near_bonds / np.linalg.norm(near_bonds, axis=1)[:, None]
        self.edges, self.weights, self.near, self.bonds, self.directions = _edge_arrays(
            self.grid, side, offsets, strengths, near_indices, near_bonds, near_axes,
        )

        self.target = E.TI, E.TJ
        self.target_bond = base[E.TJ] - base[E.TI]
        self.target_axis = self.target_bond / np.linalg.norm(self.target_bond)
        k = np.fft.fftfreq(side) * 2 * np.pi
        kx, ky, kz = k[:, None, None], k[None, :, None], k[None, None, :]
        self.symbol = np.full((side, side, side), E.K)
        phase = np.empty_like(self.symbol)
        for (dx, dy, dz), w in self.stencil:
            np.add(kx * dx + ky * dy, kz * dz, out=phase)
            np.cos(phase, out=phase)
            phase *= w
            self.symbol -= phase
        self.omega = np.sqrt(self.symbol / E.MASS)

    @cached_property
    def _scalar_seed_fft(self):
        # The two target velocities share one polarization, so one spatial
        # transform suffices; restore that polarization after taking site rows.
        seed = np.zeros_like(self.symbol)
        seed[tuple(self.grid[E.TI])] = 1
        seed[tuple(self.grid[E.TJ])] = -1
        return np.fft.fftn(seed)

    @cached_property
    def seed_fft(self):
        """Full vector seed spectrum, if requested by a diagnostic."""
        return self._scalar_seed_fft[..., None] * self.target_axis

    def predicted_kernel(self, t):
        if self.m == 2:
            return super().predicted_kernel(t)
        phase = self.omega * t
        np.cos(phase, out=phase)
        field = np.fft.ifftn(self._scalar_seed_fft * phase).real
        return field[tuple(self.grid.T)][:, None] * self.target_axis
