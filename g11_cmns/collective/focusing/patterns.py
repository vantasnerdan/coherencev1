"""Velocity preparations on the periodic PdD D register; dynamics live in the engine."""

from itertools import product

import numpy as np


def _box_lengths(box):
    box = np.asarray(box, dtype=float)
    return np.diag(box) if box.ndim == 2 else np.broadcast_to(box, (3,))


def _minimum_image(delta, box):
    return delta - box * np.rint(delta / box)


def _bond_geometry(positions, box, target):
    i, j = target
    delta = _minimum_image(positions[j] - positions[i], box)
    return positions[i] + delta / 2, np.linalg.norm(delta)


def launch_diagnostics(v, target, axis):
    """Launch diagnostics, kinetic-only even if the caller supplies a displaced state."""
    i, j = target
    site_norm2 = np.einsum('ij,ij->i', v, v)
    total = float(site_norm2.sum())
    relative = v[i] - v[j]
    mean = v.mean(axis=0)
    return dict(
        participation_at_launch=total**2 / float(site_norm2 @ site_norm2) if total else None,
        inward_relative_speed_at_launch_A_per_fs=float(relative @ axis),
        uniform_velocity_fraction_at_launch=len(v) * float(mean @ mean) / total if total else None,
        max_site_speed_A_per_fs=float(np.sqrt(site_norm2.max())),
        instantaneous_target_fraction_at_launch=float(relative @ relative) / (2 * total) if total else None,
        target_sites_kinetic_fraction_at_launch=float(site_norm2[i] + site_norm2[j]) / total if total else None,
    )


def _perpendicular(axis):
    direction = np.eye(3)[np.argmin(np.abs(axis))]
    direction -= (direction @ axis) * axis
    return direction / np.linalg.norm(direction)


def capped_direct_velocity(field, target, speed, axis, *, cap_multiplier=2.0):
    """Maximize field:v with both target rows zero, fixed total energy and a site cap."""
    n = len(field)
    remaining = np.ones(n, dtype=bool)
    remaining[list(target)] = False
    h = field[remaining]
    a = np.linalg.norm(h, axis=1)
    nonzero = a > 0
    remaining_rms = speed * np.sqrt(n / (n - 2))
    cap = cap_multiplier * remaining_rms
    budget = n * speed**2
    if cap_multiplier < 1:
        raise ValueError('site cap cannot hold the requested energy outside the target')

    r = np.zeros(n - 2)
    if cap_multiplier == 1:
        r.fill(remaining_rms)
    elif nonzero.sum() * cap**2 <= budget:
        # This branch has at least one null coefficient when cap_multiplier > 1.
        r[nonzero] = cap
        r[~nonzero] = np.sqrt((budget - nonzero.sum() * cap**2) / (~nonzero).sum())
    else:
        active = np.flatnonzero(nonzero)
        # Common scaling does not change the active-set solution.
        scaled = a / a.max()
        while active.size:
            amplitudes = scaled[active]
            trial = (np.sqrt(budget) / np.linalg.norm(amplitudes)) * amplitudes
            saturated = trial > cap
            if not saturated.any():
                r[active] = trial
                break
            r[active[saturated]] = cap
            budget -= saturated.sum() * cap**2
            active = active[~saturated]

    directions = np.broadcast_to(_perpendicular(axis), h.shape).copy()
    directions[nonzero] = h[nonzero] / a[nonzero, None]
    # Null coefficients do not change the objective; alternate their signs rather than add a coherent offset.
    null_rows = np.flatnonzero(~nonzero)
    directions[null_rows[1::2]] *= -1
    v = np.zeros_like(field, dtype=float)
    v[remaining] = r[:, None] * directions
    return v


def direct_refinements(v, target, axis):
    """Four nearby polarizations, preserving the total energy and every site cap."""
    if np.any(v[list(target)]):
        raise ValueError('direct refinement requires both target rows zero')
    tangent = _perpendicular(axis)
    for pole in (tangent, np.cross(axis, tangent)):
        parallel = (v @ pole)[:, None] * pole
        tangent_v = v - parallel
        transverse = np.cross(pole, tangent_v)
        theta = np.deg2rad(12)
        base = parallel + np.cos(theta) * tangent_v
        yield base + np.sin(theta) * transverse
        yield base - np.sin(theta) * transverse


def _wave_basis(positions, box, target, axis):
    """Broad, longitudinal reciprocal waves, not independent site controls."""
    center, bond = _bond_geometry(positions, box, target)
    cell_steps = 2 * np.pi / box
    limits = np.floor(np.pi / (bond * cell_steps) + 1e-12).astype(int)
    waves = np.array(list(product(*(range(-m, m + 1) for m in limits)))) * cell_steps
    wavenumbers = np.linalg.norm(waves, axis=1)
    keep = (wavenumbers * bond <= np.pi + 1e-12) & ((waves @ axis) > 0.45 * wavenumbers)
    waves, wavenumbers = waves[keep], wavenumbers[keep]
    order = np.argsort(wavenumbers, kind='stable')
    if len(order) > 16:
        order = order[np.linspace(0, len(order) - 1, 16, dtype=int)]
    waves = waves[order]
    polarizations = waves / wavenumbers[order, None]
    sine = np.sin((positions - center) @ waves.T)
    modes = (sine[:, None, :] * polarizations.T).reshape(3 * len(positions), len(waves))
    basis, _ = np.linalg.qr(modes, mode='reduced')

    # Opposite points relative to the center have opposite velocities. Enforce zero launch
    # on each target row in the same reciprocal subspace, not by shaping a target site impulse.
    restricted = basis.reshape(len(positions), 3, -1)[list(target)].reshape(6, -1)
    _, singular, vh = np.linalg.svd(restricted, full_matrices=True)
    rank = int(np.count_nonzero(singular > 1e-12 * singular[0]))
    return basis @ vh[rank:].T, len(waves)


def _collision_velocity(positions, box, target, axis, speed):
    center, bond = _bond_geometry(positions, box, target)
    delta = _minimum_image(positions - center, box) / bond
    along = delta @ axis
    across = delta - along[:, None] * axis
    across2 = np.einsum('ij,ij->i', across, across) / 1.05**2
    v = np.zeros_like(positions)
    packet_norm = speed * np.sqrt(len(positions) / 2)
    for side in (-1, 1):
        rho2 = ((along - side * 0.95) / 0.36)**2 + across2
        envelope = np.maximum(1 - rho2, 0)**2
        envelope[list(target)] = 0
        norm = np.linalg.norm(envelope)
        if norm == 0:
            raise ValueError('no D sites inside one of the compact collision packets')
        v += ((-side * packet_norm / norm) * envelope[:, None]) * axis
    return v


def velocity_candidates(kind, positions_A, box_A, target, axis, speeds_A_per_fs,
                        predicted_kernel=None, *, focus_times_fs=()):
    """Yield named launch dicts; speed is a scalar RMS for N equal masses and q=0.

    target=(i,j), axis is the unit bond from i toward j, and predicted_kernel(t)
    returns the signed inward (v_i - v_j)·axis adjoint in velocity space at time t (fs).
    Engine owns propagation, support, finite-horizon scoring and reversal of reverse_seed.
    """
    positions = np.asarray(positions_A, dtype=float)
    box = _box_lengths(box_A)
    axis = np.asarray(axis, dtype=float)
    speed = float(speeds_A_per_fs)
    n = len(positions)

    def record(name, v, role, time=None, **extras):
        return dict(name=name, q=np.zeros_like(v), v=v, role=role, focus_time_fs=time,
                    **launch_diagnostics(v, target, axis), **extras)

    if kind in ('reverse_seed', 'preload'):
        v = np.zeros_like(positions)
        sign = -1 if kind == 'reverse_seed' else 1
        v[list(target)] = np.array([sign, -sign])[:, None] * (speed * np.sqrt(n / 2) * axis)
        role = 'outward_seed_for_reversal' if sign < 0 else 'target_preload_control'
        yield record(kind, v, role)
    elif kind == 'collision':
        v = _collision_velocity(positions, box, target, axis, speed)
        yield record(kind, v, 'two_localized_approaching_packets')
    elif kind in ('wave', 'direct'):
        if predicted_kernel is None or not len(focus_times_fs):
            raise ValueError(f'{kind} requires a harmonic adjoint and positive focus times')
        if kind == 'wave':
            basis, wave_count = _wave_basis(positions, box, target, axis)
        for t in focus_times_fs:
            t = float(t)
            if not np.isfinite(t) or t <= 0:
                raise ValueError('focus times must be finite and positive')
            h = np.asarray(predicted_kernel(t), dtype=float)
            if kind == 'direct':
                v = capped_direct_velocity(h, target, speed, axis)
                role, controls = 'capped_distributed_adjoint', {}
            else:
                v = (basis @ (basis.T @ h.ravel())).reshape(n, 3)
                v[list(target)] = 0  # Enforce the reciprocal null constraint to machine precision.
                norm = np.linalg.norm(v)
                if norm <= 1e-12 * np.linalg.norm(h):
                    continue  # No resolved harmonic focusing for this time in the chosen wave subspace.
                v *= speed * np.sqrt(n) / norm
                role = 'phase_matched_longitudinal_waves'
                controls = dict(wave_mode_count=wave_count, wave_subspace_rank=basis.shape[1])
            yield record(f'{kind}_{t:g}fs', v, role, t,
                         harmonic_prediction_A_per_fs=float(np.einsum('ij,ij->', h, v)), **controls)
    else:
        raise ValueError(f'unknown launch kind {kind!r}')
