"""Fixed-time pullback of the production focusing split (eV, Å, fs)."""
from math import isqrt

import numpy as np
from numba import njit

import engine as E
from exterior import energy_force

# Max coordinate displacement of each local force probe, in Å.
_FORCE_PROBE_A = 1e-6


@njit(cache=True)
def _kick(v, force, kick):
    for i in range(len(v)):
        for a in range(3):
            v[i, a] += kick * force[i, a]


@njit(cache=True)
def _kick_drift(x, v, force, kick, half, decay):
    for i in range(len(x)):
        for a in range(3):
            v[i, a] += kick * force[i, a]
            x[i, a] += half * v[i, a]
            v[i, a] *= decay
            x[i, a] += half * v[i, a]


@njit(cache=True)
def _forward(v0, edges, weights, near, directions, step, block, kick, half, decay, stiff):
    x = np.zeros_like(v0)
    v = v0.copy()
    force = np.empty_like(x)
    checkpoints = (step + block - 1) // block
    starts_x = np.empty((checkpoints, len(x), 3))
    starts_v = np.empty_like(starts_x)
    if step:
        energy_force(x, force, edges, weights, near, directions, stiff)
    for p in range(step):
        if p % block == 0:
            starts_x[p // block] = x
            starts_v[p // block] = v
        _kick_drift(x, v, force, kick, half, decay)
        energy_force(x, force, edges, weights, near, directions, stiff)
        _kick(v, force, kick)
    return x, v, starts_x, starts_v


@njit(cache=True)
def _add_force_pullback(x, av, ax, factor, probe, plus, minus, edges, weights, near, directions, stiff):
    # Where differentiable, F = -∇V so DF.T @ av = DF @ av. Use the
    # actual force for a centered local difference, not a surrogate Hessian.
    size = 0.
    for i in range(len(av)):
        for a in range(3):
            size = max(size, abs(av[i, a]))
    if size == 0.:
        return
    width = _FORCE_PROBE_A / size
    for i in range(len(x)):
        for a in range(3):
            probe[i, a] = x[i, a] + width * av[i, a]
    energy_force(probe, plus, edges, weights, near, directions, stiff)
    for i in range(len(x)):
        for a in range(3):
            probe[i, a] = x[i, a] - width * av[i, a]
    energy_force(probe, minus, edges, weights, near, directions, stiff)
    scale = factor * (.5 / width)
    for i in range(len(x)):
        for a in range(3):
            ax[i, a] += scale * (plus[i, a] - minus[i, a])


@njit(cache=True)
def _reverse(endpoint_x, ax, av, starts_x, starts_v, edges, weights, near, directions,
             step, block, kick, half, decay, stiff):
    x, v = np.empty_like(ax), np.empty_like(av)
    history = np.empty((block - 1, len(x), 3))
    probe, plus, minus = np.empty_like(x), np.empty_like(x), np.empty_like(x)
    _add_force_pullback(endpoint_x, av, ax, kick, probe, plus, minus,
                        edges, weights, near, directions, stiff)
    drift = half * (1. + decay)
    for b in range(len(starts_x) - 1, -1, -1):
        begin = b * block
        length = min(block, step - begin)
        if length > 1:
            x[:] = starts_x[b]
            v[:] = starts_v[b]
            energy_force(x, plus, edges, weights, near, directions, stiff)
            for offset in range(1, length):
                _kick_drift(x, v, plus, kick, half, decay)
                history[offset - 1] = x
                if offset + 1 < length:
                    energy_force(x, plus, edges, weights, near, directions, stiff)
                    _kick(v, plus, kick)
        for offset in range(length - 1, -1, -1):
            for i in range(len(av)):
                for a in range(3):
                    av[i, a] = decay * av[i, a] + drift * ax[i, a]
            if begin + offset > 0:
                # The end kick of the previous step and start kick of the
                # next use the same x and the same carried velocity dual.
                pos = history[offset - 1] if offset else starts_x[b]
                _add_force_pullback(pos, av, ax, 2. * kick, probe, plus, minus,
                                    edges, weights, near, directions, stiff)
            # x0 is fixed: its adjoint is not needed.


def _terminal(target_bond, target, x, v, objective):
    i, j = target
    ax = np.zeros_like(x)
    av = np.zeros_like(v)
    relative = v[i] - v[j]
    if objective == 'K':
        av[i] = .5 * E.MASS * relative
        av[j] = -av[i]
        return .25 * E.MASS * float(relative @ relative), ax, av
    bond = target_bond + x[j] - x[i]
    radius = float(np.linalg.norm(bond))
    axis = bond / radius
    if objective == 'distance':
        ax[i], ax[j] = axis, -axis
        return -radius, ax, av
    closing = float(relative @ axis)
    if closing <= 0.:
        return 0., ax, av
    slope = .5 * E.MASS * closing
    av[i] = slope * axis
    av[j] = -av[i]
    ax[j] = (slope / radius) * (relative - closing * axis)
    ax[i] = -ax[j]
    return .25 * E.MASS * closing**2, ax, av


def endpoint_gradient(lattice, v0, dt, step, objective, tau_fs=0., exterior='soft'):
    """Return ``(value, dvalue_dv0, x_end, v_end)`` at step ``step``, with x0=0.

    Maximize 'K' (full relative kinetic eV), 'inward' (strictly closing
    radial kinetic eV), or 'distance' (negative target distance in Å).
    The gradient is unprojected; launch constraints belong to the caller.
    A positive tau_fs applies one amplitude factor exp(-dt/tau_fs) per
    step, between the two half-drifts; zero is conservative.

    Soft/stiff select the declared force continuations, not support rules.
    The scored path must be checked separately for extent and model failure.
    At force-Hessian knots the centered probe samples both sides rather
    than claiming a unique derivative. Checkpoints and one replay block
    use O(N sqrt(step)) storage without storing a full trajectory.
    """
    if objective not in ('K', 'inward', 'distance'):
        raise ValueError(f'unknown endpoint objective: {objective!r}')
    if exterior not in ('soft', 'stiff'):
        raise ValueError(f'unknown exterior: {exterior!r}')
    if step < 0 or dt <= 0. or tau_fs < 0.:
        raise ValueError('expected nonnegative step and tau_fs, and positive dt')
    v0 = np.asarray(v0, dtype=np.float64, order='C')
    if v0.shape != (lattice.n, 3):
        raise ValueError('v0 must have shape (lattice.n, 3)')
    if step == 0:
        x, v = np.zeros_like(v0), v0.copy()
        value, _, av = _terminal(lattice.target_bond, lattice.target, x, v, objective)
        return value, av, x, v
    block = isqrt(2 * step)
    half = .5 * dt
    kick = half / E.MASS
    decay = np.exp(-dt / tau_fs) if tau_fs > 0. else 1.
    force_data = lattice.edges, lattice.weights, lattice.near, lattice.directions
    scheme = step, block, kick, half, decay, exterior == 'stiff'
    x, v, starts_x, starts_v = _forward(v0, *force_data, *scheme)
    value, ax, av = _terminal(lattice.target_bond, lattice.target, x, v, objective)
    if objective != 'inward' or value != 0.:
        _reverse(x, ax, av, starts_x, starts_v, *force_data, *scheme)
    return value, av, x, v
