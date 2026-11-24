"""Measure the first surface window against BO and expose interpolation/step error."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.constants import elementary_charge, physical_constants
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

from charge_readout import CHARGE_FIELDS, write_charge


OBSERVABLES = ('induced_energy_positive_test_eV', 'forces_eV_A', 'internal_energy_eV',
               'projectile_sphere_e', 'vacuum_e_outside_projectile_1p4_excess', *CHARGE_FIELDS)


def load(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def arrays(rows):
    return {key: np.asarray([r[key] for r in rows]) for key in OBSERVABLES}


def bo_at(anchors, edges):
    z = np.asarray([r['z_A'] for r in anchors])
    ordered = np.argsort(z)
    vals = arrays(anchors)
    edges = np.asarray(edges, dtype=float)
    low, high = z[ordered[[0, -1]]]
    # The final analytic position and static grid may round on opposite sides.
    edges = np.where(np.abs(edges-low) <= 4*abs(np.spacing(low)), low, edges)
    edges = np.where(np.abs(edges-high) <= 4*abs(np.spacing(high)), high, edges)
    return {key: PchipInterpolator(z[ordered], value[ordered], axis=0, extrapolate=False)(edges)
            for key, value in vals.items()}


def peak(time_fs, values, sign=None):
    if sign is None:
        sign = np.sign(values[np.argmax(np.abs(values))])
    signed = values * sign
    i = int(np.argmax(signed))
    if signed[i] <= 0:
        return None
    mid = abs(values[i]) / 2
    left = next((j for j in range(i-1, -1, -1) if signed[j] < mid), None)
    right = next((j for j in range(i+1, len(time_fs)) if signed[j] < mid), None)
    lo = None if left is None else np.interp(mid, signed[left:left+2], time_fs[left:left+2])
    hi = None if right is None else np.interp(mid, signed[right-1:right+1][::-1], time_fs[right-1:right+1][::-1])
    return dict(value=float(values[i]), peak_fs=float(time_fs[i]),
                half_height_bracket_fs=[lo, hi],
                half_height_width_fs=None if lo is None or hi is None else float(hi-lo))


def axial_lobes(time, force):
    return [{'outward': peak(time, z, 1), 'inward': peak(time, z, -1)} for z in force[:, :, 2].T]


def profile(rows, anchors):
    z = np.asarray([r['z_A'] for r in rows])
    t = np.asarray([r['time_as'] for r in rows]) * .001
    td = arrays(rows)
    bo = bo_at(anchors, z)
    for key in OBSERVABLES:
        if not np.isfinite(bo[key]).all():
            raise ValueError('BO anchors do not cover the full trajectory')
    return t, z, td, bo


def window(data, start_z):
    t, z, td, bo = data
    selected = z <= start_z
    return (t[selected], z[selected], {k: v[selected] for k, v in td.items()},
            {k: v[selected] for k, v in bo.items()})


def initial_offset(states, anchors):
    """A ground-state BO comparator may start partially neutral while TD starts D+."""
    first = states[0]
    matches = [row for row in anchors if row['z_A'] == first['z_A']]
    if first['time_as'] != 0 or len(matches) != 1:
        raise ValueError('Initial offset requires one exact BO reference at time zero')
    return {k: np.asarray(first[k])-np.asarray(matches[0][k]) for k in OBSERVABLES}


def holdout_errors(anchors, holdouts):
    measured = arrays(holdouts)
    interpolated = bo_at(anchors, [row['z_A'] for row in holdouts])
    return {key: np.abs(measured[key]-predicted) for key, predicted in interpolated.items()}


def impulse_recoil(t, td_force, bo_force, background, speed_A_as, entry):
    """Freeze positions but estimate diagnostic velocity increments from the forces."""
    f, fb = td_force[:, 1:], bo_force[:, 1:]
    if background.shape != f.shape[1:]:
        raise ValueError('Target/neighbour forces lack the host background')
    j = int(np.searchsorted(t, entry, side='right'))
    if not 0 < j < len(t):
        raise ValueError('Trajectory does not bracket the impulse window start')
    weight = (entry-t[j-1])/(t[j]-t[j-1])
    window_time = np.r_[entry, t[j:]]

    def impulse(force):
        start = force[j-1] + weight*(force[j]-force[j-1])
        return cumulative_trapezoid(np.concatenate((start[None], force[j:]))-background,
                                    window_time, axis=0, initial=0)[1:]

    mom, mom_bo = impulse(f), impulse(fb)
    mass = physical_constants['deuteron mass'][0]
    momentum_factor = elementary_charge * 1e-5
    recoil_factor = momentum_factor**2 / (2*mass*elementary_charge)
    recoil = recoil_factor * (mom**2).sum(axis=-1)
    recoil_excess = recoil - recoil_factor*(mom_bo**2).sum(axis=-1)
    relative_v_m_s = -(momentum_factor / mass) * mom
    relative_v_m_s[..., 2] -= speed_A_as * 1e8
    pair_cm_keV = mass * (relative_v_m_s**2).sum(axis=-1) / (4 * elementary_charge * 1e3)
    return window_time[1:], mom-mom_bo, recoil, recoil_excess, pair_cm_keV


def fastest_rate(t, values):
    rates = np.diff(values) / np.diff(t)
    i = int(np.argmax(np.abs(rates)))
    return {'signed_rate_per_fs': float(rates[i]), 'time_bracket_fs': [float(t[i]), float(t[i+1])]}


def run(args):
    rows = load(args.trajectory)
    inputs, *states, end = rows
    if end['kind'] != 'complete' or inputs['dt_as'] != 4 or inputs['energy_lab_keV'] != .5:
        raise ValueError('Expected the completed 0.5-keV, 4-as entry ray')
    geom = inputs['geometry']
    sites = [-1, *geom['receiving_indices']]
    top = geom['surface_z_A']
    bo_rows = load(args.bo)
    anchors = [r for r in bo_rows if not r['holdout']]
    holdouts = [r for r in bo_rows if r['holdout']]
    errors = holdout_errors(anchors, holdouts)
    full = profile(states, anchors)
    t, z, td, bo = window(full, top+.8)
    phi, fkey = 'induced_energy_positive_test_eV', 'forces_eV_A'
    de, df = td[phi]-bo[phi], td[fkey]-bo[fkey]
    offset = initial_offset(states, anchors)
    speed = inputs['speed_A_as']
    entry = (states[0]['z_A']-(top+.8))/(1000*speed)
    ti, impulse, recoil, excess, pair_energy = impulse_recoil(
        full[0], full[2][fkey], full[3][fkey], np.asarray(args.host_force), speed, entry)
    out = {
        'step_as': 4,
        'edge_start_distance_A': .8, 'edge_end_distance_A': float(z[-1]-top),
        'target_site': sites[1], 'potential_force_sites': sites,
        'bo_holdout_indices': [row['index'] for row in holdouts],
        'bo_holdout_potential_and_vector_force_errors': {key: errors[key].tolist() for key in (phi, fkey)},
        'incoming_TD_minus_BO': {key: value.tolist() for key, value in offset.items()},
        'potential_peaks_eV': [peak(t, curve) for curve in de.T],
        'force_z_lobes_eV_per_A': axial_lobes(t, df),
        'target_steepest_resolved_rate': {'potential_eV': fastest_rate(t, de[:, 1]),
                                         'force_z_eV_per_A': fastest_rate(t, df[:, 1, 2])},
        'receiver_impulse_start_fs': entry,
        'receiver_final_impulse_TD_minus_BO_eV_fs_per_A': impulse[-1].tolist(),
        'receiver_peak_TD_recoil_eV': recoil.max(axis=0).tolist(),
        'receiver_signed_peak_excess_recoil_eV': [peak(ti, curve) for curve in excess.T],
        'projectile_charge_approach': write_charge(args.output.with_suffix('.charge.csv'), full, top),
    }
    with args.output.open('x') as fh:
        json.dump(out, fh, allow_nan=False, indent=2)
    with args.output.with_suffix('.csv').open('x', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['time_fs', 'distance_from_outer_Pd_plane_A', 'site', 'phi_td_minus_bo_eV',
                         'force_td_minus_bo_x_eV_A', 'force_td_minus_bo_y_eV_A', 'force_td_minus_bo_z_eV_A',
                         'phi_td_eV', 'phi_bo_eV', 'f_td_x_eV_A', 'f_td_y_eV_A', 'f_td_z_eV_A',
                         'f_bo_x_eV_A', 'f_bo_y_eV_A', 'f_bo_z_eV_A',
                         'phi_td_minus_bo_change_from_incoming_eV', 'force_z_change_from_incoming_eV_A'])
        time_full, z_full, td_full, bo_full = full
        for i, (time, depth) in enumerate(zip(time_full, z_full-top)):
            for j, receiver in enumerate(sites):
                dphi = td_full[phi][i,j]-bo_full[phi][i,j]
                dforce = td_full[fkey][i,j]-bo_full[fkey][i,j]
                writer.writerow([time, depth, receiver, dphi, *dforce, td_full[phi][i,j], bo_full[phi][i,j],
                                 *td_full[fkey][i,j], *bo_full[fkey][i,j], dphi-offset[phi][j],
                                 dforce[2]-offset[fkey][j,2]])
    with args.output.with_suffix('.receivers.csv').open('x', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['time_fs', 'site', 'impulse_TD_minus_BO_x_eV_fs_per_A',
                         'impulse_TD_minus_BO_y_eV_fs_per_A', 'impulse_TD_minus_BO_z_eV_fs_per_A',
                         'frozen_recoil_eV', 'excess_recoil_eV', 'conditional_pair_cm_keV'])
        for i, time in enumerate(ti):
            for j, receiver in enumerate(sites[1:]):
                writer.writerow([time, receiver, *impulse[i,j], recoil[i,j], excess[i,j], pair_energy[i,j]])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('trajectory', type=Path)
    p.add_argument('pilot', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--bo', type=Path, required=True)
    args = p.parse_args()
    args.host_force = next(x['forces_eV_A'] for x in load(args.pilot) if x['phase'] == 'host_BO')
    run(args)
