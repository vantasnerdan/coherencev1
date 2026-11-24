"""Compare projectile-local electron populations along the entire entry ray."""
import csv

import numpy as np


CHARGE_FIELDS = {
    'projectile_1s_overlap_e': 'overlap_1s_e',
    'projectile_host_excess_sphere_e': 'host_subtracted_sphere_1p4_A_e',
}


def summarize(t, distance, td, bo):
    """Return the first half-BO crossing and largest resolved TD secant slope."""
    def location(i, f=0):
        return {'time_fs': float(t[i]+f*(t[i+1]-t[i])),
                'distance_from_outer_Pd_plane_A': float(distance[i]+f*(distance[i+1]-distance[i]))}

    half = None
    gap = td-.5*bo
    if bo[0] > 0 and gap[0] >= 0:
        half = dict(location(0), status='already_at_or_above_at_start')
    else:
        for i in np.flatnonzero((gap[:-1] < 0) & (gap[1:] >= 0)):
            f = -gap[i]/(gap[i+1]-gap[i])
            if bo[i]+f*(bo[i+1]-bo[i]) > 0:
                half = dict(location(i, f), status='crossed',
                            bracket_fs=[float(t[i]), float(t[i+1])])
                break
    rate = np.diff(td)/np.diff(t)
    i = int(np.argmax(rate))
    fastest = (dict(location(i, .5), rate_e_per_fs=float(rate[i]),
                    bracket_fs=[float(t[i]), float(t[i+1])]) if rate[i] > 0 else None)
    return {
        'half_BO': half,
        'fastest_rise': fastest,
        'endpoint': {'time_fs': float(t[-1]), 'distance_from_outer_Pd_plane_A': float(distance[-1]),
                     'TD_e': float(td[-1]), 'BO_e': float(bo[-1]),
                     'BO_minus_TD_e': float(bo[-1]-td[-1]),
                     'TD_fraction_of_BO': float(td[-1]/bo[-1]) if bo[-1] > 0 else None},
    }


def write_charge(path, profile, top):
    t, z, td, bo = profile
    distance = z-top
    populations = [(td[k], bo[k]) if td[k].ndim == 1 else (td[k][:, 1], bo[k][:, 1])
                   for k in CHARGE_FIELDS]
    with path.open('x', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['time_fs', 'distance_from_outer_Pd_plane_A',
                         *[f'{measure}_{label}' for measure in CHARGE_FIELDS.values() for label in ('TD', 'BO')]])
        writer.writerows(zip(t, distance, *(v for pair in populations for v in pair)))
    return {measure: summarize(t, distance, q_td, q_bo)
            for measure, (q_td, q_bo) in zip(CHARGE_FIELDS.values(), populations)}
