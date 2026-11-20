"""Profile a shared two-class occupancy law against published scalar shift summaries.

This is not a refit of reaction counts: the original energy/bin weights are not
available for these summaries. A penetration response at declared CM pivots
relates the scalar fits to site fractions. Use python -m g7_micro.occupancy.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from g4_beam.forward.model import BG as BG_DD

KB_EV_K = 8.617333262e-5
M_D, M_F = 2.01410177812, 18.99840316273
BG_FD = BG_DD * 9 * np.sqrt(2 * M_F / (M_D + M_F))
PD_E_CM = 1e3 * np.array([3.089, 9.200]) * M_D / (M_D + M_F)
# Raiola et al. / EXFOR O1379001 summary. NaN marks an upper limit without sigma.
TI = np.array([[-10., 2.10, .030, np.nan], [50., 1.10, .050, np.nan],
               [100., .26, .250, .040], [150., .23, .295, .040], [200., .20, .290, .065]])
# Cvetinović et al., F on D. x is initial H(D)/Pd, not an isotope assay by site.
PD = np.array([[24., .47, 18.2, 3.3], [24., .70, 3.2, 1.9]])
PD_BASE = 2.19  # their F+D adiabatic reference, not the D+D bulk anchor
SOURCES = {
    'Ti': 'https://raw.githubusercontent.com/IAEA-NDS/exfor_master/master/exforall/O13/O1379.x4',
    'Ti_paper': 'https://doi.org/10.1088/0954-3899/31/11/002',
    'Pd': 'https://doi.org/10.1016/j.physletb.2023.137684',
    'Pd_table': 'https://doi.org/10.1088/1742-6596/2586/1/012101',
    'partitioning': 'https://doi.org/10.1038/s41586-026-11238-4',
}


def occupied_fraction(x, temperature_k, capacity, n_s, binding_eV):
    """Special-site D fraction for x>0, x<=capacity, 0<=n_s<=capacity.

    x and both capacities count atoms or single-occupancy sites per host atom.
    The special capacity replaces ordinary sites. The smaller quadratic root
    applies to extra binding >=0; it satisfies the common chemical potential.
    """
    k = np.exp(-binding_eV / (KB_EV_K * temperature_k))
    a = capacity - n_s - x
    b = n_s + x + k*a
    d = (n_s-x)**2 + k*(2*((n_s+x)*a+2*n_s*x) + k*a*a)
    return 2*n_s / (b + np.sqrt(d))


def capacity_for_fraction(x, temperature_k, capacity, frac, binding_eV):
    """Invert occupied_fraction for site capacity n_s; no fitted f at f=0/1."""
    s = x*frac
    k = np.exp(-binding_eV/(KB_EV_K*temperature_k))
    return s + (capacity-x) * (k*s) / (x-s+k*s)


def log_enhancement(shift, e_cm, gamow):
    """Penetration-only cross-section enhancement, Cvetinović et al. Eq. (2)."""
    low, high = np.sqrt(e_cm), np.sqrt(e_cm+shift)
    return gamow*shift/(low*high*(low+high))


def inverse_log_enhancement(log_factor, e_cm, gamow):
    t = log_factor*np.sqrt(e_cm)/gamow
    return e_cm*t*(2-t)/((1-t)**2)


def mixed_shift(frac, u_base, u_high, e_cm, gamow):
    """Yield-equivalent shift and yield-weighted arithmetic shift, both keV."""
    b = log_enhancement(u_base, e_cm, gamow)
    r = np.expm1(log_enhancement(u_high, e_cm, gamow) - b)
    rel = 1+frac*r
    return (inverse_log_enhancement(b+np.log1p(frac*r), e_cm, gamow),
            u_base + (u_high-u_base)*frac*(1+r)/rel)


def fraction_for_shift(shift, u_base, u_high, e_cm, gamow):
    b = log_enhancement(u_base, e_cm, gamow)
    return np.expm1(log_enhancement(shift, e_cm, gamow)-b) / np.expm1(log_enhancement(u_high, e_cm, gamow)-b)


def ti_profile(e_cm, capacity, u_base, *, curve=False):
    temp, x, obs, sd = TI.T
    temp = temp+273.15
    limits = np.isnan(sd)
    n_base = log_enhancement(u_base, e_cm, BG_DD)
    detected = ~limits
    if np.any(obs[limits] < u_base):
        return {'feasible': False, 'reason': 'ordinary shift exceeds a reported upper limit'}
    trials = []
    u_values = np.geomspace(.08, 100., 91)
    if curve:
        u_values = np.r_[u_values, np.geomspace(110., 1e4, 42)]
    for i, u_high in enumerate(u_values):
        rel = np.expm1(log_enhancement(u_high, e_cm, BG_DD)-n_base)
        f_lim = np.clip(fraction_for_shift(obs[limits], u_base, u_high, e_cm, BG_DD), 0., 1.)
        for binding in np.linspace(0., .5, 51):
            cap = np.min(capacity_for_fraction(x[limits], temp[limits], capacity, f_lim, binding))

            def loss(n_s):
                frac = occupied_fraction(x[detected], temp[detected], capacity, n_s, binding)
                pred = inverse_log_enhancement(n_base+np.log1p(frac*rel), e_cm, BG_DD)
                return np.sum(((pred-obs[detected])/sd[detected])**2)

            if cap:
                # Scale the scalar search to carry small occupancies accurately.
                r = minimize_scalar(lambda t: loss(cap*t), bounds=(0., 1.), method='bounded')
                chi2, n_s = min((float(r.fun), cap*float(r.x)), (loss(cap), cap), (loss(0.), 0.))
            else:
                chi2, n_s = loss(0.), 0.
            trials.append((float(chi2), float(n_s), float(binding), float(u_high), float(cap), i))
    tr = np.asarray(trials)

    def describe(row):
        chi2, n_s, binding, u_high = map(float, row[:4])
        frac = occupied_fraction(x, temp, capacity, n_s, binding)
        eq, av = mixed_shift(frac, u_base, u_high, e_cm, BG_DD)
        return dict(chi2=chi2, N_s_per_Ti=n_s, dE_eV=binding, U_high_keV=u_high,
                    f_occupied=frac.tolist(), U_equivalent_keV=eq.tolist(), U_yield_weighted_keV=av.tolist())

    imin = np.argmin(tr[:, 0])
    best = tr[imin]
    allowed = tr[tr[:, 0] <= best[0]+1., 1:4]
    examples = []
    for try_u in (.3, .6, 1.5, 5., 30., 100.):
        nearest = tr[np.argmin(abs(np.log(tr[:, 3]/try_u))), 3]
        block = tr[tr[:, 3] == nearest]
        examples.append(describe(block[np.argmin(block[:, 0])]))
    fixed_binding = tr[tr[:, 2] == .2]
    shared = fixed_binding[np.argmin(abs(np.log(fixed_binding[:, 3]/30.)))]
    result = dict(feasible=True, E_cm_pivot_keV=e_cm, capacity_per_Ti=capacity, U_base_keV=u_base,
                  best=describe(best), examples=examples, dE_0p2_and_U_near_30_example=describe(shared),
                  delta_chi2_1_profile=dict(keys=['N_s_per_Ti', 'dE_eV', 'U_high_keV'],
                                             min=allowed.min(axis=0).tolist(), max=allowed.max(axis=0).tolist()),
                  scanned=dict(U_high_keV=[.08, float(u_values[-1]), len(u_values)], dE_eV=[0., .5, 51]))
    if not curve:
        return result
    threshold = best[0]+1.
    rows = []
    for i, u_high in enumerate(u_values):
        rel = np.expm1(log_enhancement(u_high, e_cm, BG_DD)-n_base)
        block = tr[tr[:, 5] == i]
        at_min = block[np.argmin(block[:, 0])]
        band = []
        for _, n_s, de, _, cap, _ in block[block[:, 0] <= threshold]:
            def edge(t):
                frac = occupied_fraction(x[detected], temp[detected], capacity, t*cap, de)
                pred = inverse_log_enhancement(n_base+np.log1p(frac*rel), e_cm, BG_DD)
                return np.sum(((pred-obs[detected])/sd[detected])**2)-threshold
            tmin = n_s/cap if cap else 0.
            left = brentq(edge, 0., tmin) if edge(0.) > 0 else 0.
            right = brentq(edge, tmin, 1.) if edge(1.) > 0 else 1.
            band.append((left*cap, right*cap))
        rows.append(dict(U_high_keV=float(u_high), N_s_per_Ti_best=float(at_min[1]),
                         dE_eV_best=float(at_min[2]), chi2=float(at_min[0]),
                         N_s_per_Ti_band_low=min((x[0] for x in band), default=None),
                         N_s_per_Ti_band_high=max((x[1] for x in band), default=None)))
    result['curve'] = rows
    return result


def pd_constraints(e_cm):
    temp, x, obs, sd = PD.T
    temp = temp+273.15
    hard, soft = obs
    b = log_enhancement(PD_BASE, e_cm, BG_FD)
    def excess(shift):
        return np.expm1(log_enhancement(shift, e_cm, BG_FD)-b)

    # Monotonic filling of the special sites at fixed T puts an exact bound on
    # their occupied fraction ratio under unchanged capacities, for any dE.
    max_ratio = x[1]/x[0]
    def max_hard(soft_pred):
        return inverse_log_enhancement(b + np.log1p(max_ratio*excess(soft_pred)), e_cm, BG_FD)

    def edge_loss(soft_pred):
        hard_pred = min(hard, max_hard(soft_pred))
        return ((hard_pred-hard)/sd[0])**2 + ((soft_pred-soft)/sd[1])**2

    edge = minimize_scalar(edge_loss, bounds=(PD_BASE, hard), method='bounded')
    examples = []
    for u_high in (20., 30., 50., 100.):
        frac = fraction_for_shift(obs, PD_BASE, u_high, e_cm, BG_FD)
        _, av = mixed_shift(frac, PD_BASE, u_high, e_cm, BG_FD)
        for binding in (0., .1, .2, .4):
            examples.append(dict(dE_eV=binding, U_high_keV=u_high, f_occupied_hard_soft=frac.tolist(),
                                 N_s_per_Pd_hard_soft=capacity_for_fraction(x, temp, 1., frac, binding).tolist(),
                                 occupied_special_per_Pd_hard_soft=(x*frac).tolist(),
                                 U_yield_weighted_keV_hard_soft=av.tolist()))
    curve = []
    binding_grid = np.linspace(0., .5, 101)[:, None]
    for u_high in np.unique(np.r_[np.geomspace(12., 1e4, 153), hard, soft+sd[1], hard+sd[0]]):
        if u_high < max(PD_BASE, *(obs-sd)):
            continue
        # The source reports separate marginal errors, not the pair covariance.
        # Intersect a rectangular +/-1 sigma sensitivity with physical shifts.
        low = np.maximum(PD_BASE, obs-sd)
        high = np.minimum(u_high, obs+sd)
        nlo = capacity_for_fraction(x, temp, 1., fraction_for_shift(low, PD_BASE, u_high, e_cm, BG_FD), binding_grid)
        nhi = capacity_for_fraction(x, temp, 1., fraction_for_shift(high, PD_BASE, u_high, e_cm, BG_FD), binding_grid)
        central = None
        if u_high >= hard:
            f = fraction_for_shift(obs, PD_BASE, u_high, e_cm, BG_FD)
            ns = capacity_for_fraction(x, temp, 1., f, .2)
            central = dict(N_s_hard=float(ns[0]), N_s_soft=float(ns[1]), delta_N_s=float(ns[0]-ns[1]))
        curve.append(dict(U_high_keV=float(u_high), central_dE_0p2=central,
                          hard_band=[float(nlo[:, 0].min()), float(nhi[:, 0].max())],
                          soft_band=[float(nlo[:, 1].min()), float(nhi[:, 1].max())],
                          delta_band=[float((nlo[:, 0]-nhi[:, 1]).min()), float((nhi[:, 0]-nlo[:, 1]).max())]))
    return dict(E_cm_pivot_keV=e_cm, central_excess_ratio=float(excess(hard)/excess(soft)),
                loading_only_max_excess_ratio=float(max_ratio),
                loading_only_max_hard_U_keV_at_soft_central=float(max_hard(soft)),
                loading_only_relaxed_min_chi2=float(edge.fun),
                loading_only_edge_U_keV_hard_soft=[float(max_hard(edge.x)), float(edge.x)],
                U_high_keV_central_min=float(hard), conditional_density_examples=examples, curve=curve)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('occupancy_results.json'))
    args = parser.parse_args()
    out = dict(sources=SOURCES,
               Ti=dict(observations=[dict(T_C=t, D_per_Ti=x, U_keV=u,
                                          sigma_keV=None if np.isnan(sd) else sd, upper_limit=bool(np.isnan(sd)))
                                     for t, x, u, sd in TI],
                       conditional_fits=[ti_profile(*spec, curve=i < 3) for i, spec in enumerate(((5., 2.2, 0.), (3.85, 2.2, 0.), (10., 2.2, 0.),
                                                                      (5., 2.1, 0.), (5., 3., 0.), (5., 2.2, .025), (5., 2.2, .04)))]),
               Pd=dict(observations=[dict(preparation=name, T_C=t, initial_H_D_per_Pd=x, U_keV=u, sigma_keV=sd)
                                     for name, (t, x, u, sd) in zip(('hard cold rolled','soft commercial'), PD)],
                       U_base_keV=PD_BASE, assumed_capacity_per_Pd=1.,
                       E_cm_incident_range_keV=PD_E_CM.tolist(), gamow_sqrt_keV=float(BG_FD),
                       response_pivots=[pd_constraints(e) for e in (500., *PD_E_CM)]),
               assumptions=[
                   'FD site occupancy for independent sites; common chemical potential fixed by measured loading, no assumed equilibrium across different specimens.',
                   'Ti capacity >=2.1 is a descriptive effective inventory, not identification of one Ti hydride phase. Thermal and loading effects covary.',
                   'Treat Ti upper limits as hard inequalities; no confidence level quoted. Chi-square uses detected summary points only, not a calibrated joint likelihood.',
                   'Ti uses a penetration-matched 5 keV CM summary proxy with pivot sensitivity. No Ti raw curve in O1379; Pt/Hf energy grids do not supply its original fitting weights.',
                   'Pd uses F+D penetration Eq. (2) with pivot sensitivity, not a D+D kernel or original SRIM/normalized thick-target fit. 24 C is loading temperature.',
                   'Pd paper quotes initial total H(D)/Pd; assume equal isotope partition when using total occupancy as a proxy for the D fraction. Fluorine runs lose D; control normalization is not an occupancy assay.',
                   'Curve bands: Ti is the sampled binding-profile envelope where descriptive chi2 <= global sampled minimum+1 (hard limits retained), not confidence coverage. Pd hard/soft +/-1 sigma rectangle, shared dE in [0,.5] eV, with exact within-dE extrema of the inventory difference; curves are sampled at each CM pivot separately.',
                   'Pd shares the site law, not known Ti parameters. Different fitted shifts cannot isolate cold work, supplier, site density, isotope partition, or binding energy.',
                   'Huang Fig. 4 resolves preferential occupancy of coherent L12 versus D022 in Inconel; it does not measure screening or fix binding energy / site counts in Pd or Ti.',
               ])
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
