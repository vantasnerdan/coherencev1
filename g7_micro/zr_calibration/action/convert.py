"""One-off common-support DD action comparison for draft 453206951 (entry_2800)."""
import argparse
import hashlib
import json
from pathlib import Path
from math import acos, expm1, isfinite, log, sqrt

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import quad
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator, PPoly
from scipy.optimize import brentq

from g4_beam.forward.consumer import COULOMB_EV_M, QE, WKB_NUCLEAR_RADIUS_M
from g4_beam.forward.model import BG, reaction_sigma
from g4_beam.forward.thermal import MU_KG
from g6_pairlaw.script.material_pair import read_interaction

C = COULOMB_EV_M * 1e10  # eV Å
R0 = WKB_NUCLEAR_RADIUS_M * 1e10  # Å
K = 2 * sqrt(2 * MU_KG * QE) / 1.054571817e-34 * 1e-10  # (Å sqrt(eV))^-1
FLOOR = 50.0  # eV CM
MAP_SHA = '2d4f74195bfb19f82b1f07be3de2047b19eeeddaeaef6eb8ade309f6ee87e5be'
BASE = '9c0fb426dfe9498e1049eb05c95daeed9b590205'
SHAPES = ('chord_hermite', 'zero_slope_hermite', 'linear_W')
ANCHORS = (.05, .10, .125, .25, .5, 1., 1.25, 2., 2.5, 3., 3.85, 4., 5., 6.5, 10.)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def bind(data, expected):
    if len(data) != expected['bytes'] or digest(data) != expected['sha256']:
        raise ValueError('compact asset or finite-curve binding does not match')


def contact_bound(contact, raw, rel):
    bind(raw, contact['finite_curve'])
    if (contact['finite_curve']['path'] != rel or contact['R_A'] != 0
            or contact['schema'] != 'hydride-signed-contact-sidecar-v1'
            or contact['observable'] != 'signed_regular_contact_limit_of_target_conditioned_BO_total'
            or contact['charge'] != 1 or contact['electronic_branch'] != 'bare-no-added-electron'
            or contact['derivative_eV_per_A'] is not None or not isfinite(contact['signedW0_eV'])
            or not contact['geometry']['passed']):
        raise ValueError('not an admitted q+1 contact with explicitly unqualified slope')


def load(root):
    raw = (root / 'delivery-map.json').read_bytes()
    if digest(raw) != MAP_SHA:
        raise ValueError('expect the exact draft 453206951 delivery map')
    delivery = json.loads(raw)
    records = {row['capsule_path']: row for row in delivery['small_assets'] if row['capsule_path']}

    def get(rel):
        data = (root / rel).read_bytes()
        bind(data, records[rel])
        return data

    sites = {}
    for role in ('plain', 'oxygen'):
        rel = f'curves/{role}_q1_BO_pair.json'
        raw = get(rel)
        contact = json.loads(get(f'{role}_q1_contact.json'))
        contact_bound(contact, raw, rel)
        r, v, _ = read_interaction(json.loads(raw))
        sites[role] = r, v - C/r, contact, raw, rel
    clean, oxy = [row[2] for row in sites.values()]
    if (not np.array_equal(clean['ray'], oxy['ray']) or
            not np.array_equal(clean['geometry']['cell_A'], oxy['geometry']['cell_A']) or
            clean['absolute_convention']['energy'] != oxy['absolute_convention']['energy']):
        raise ValueError('matched ray, cell and energy convention required')
    pair = json.loads(get('paired-comparisons.json'))
    if pair['paired_cell_qualified'] or delivery['paired_cell_qualified']:
        raise ValueError('endpoint controls are not a cell qualification in this snapshot')
    for role, (r, w, contact, _, _) in sites.items():
        base = pair['pairs']['base']
        if (r[0] != .1 or abs(base[f'{role}_R01_lowering_eV'] + w[0]) > 1e-6 or
                contact['signedW0_eV'] != base[f'{role}_signedW0_eV']):
            raise ValueError('compact comparison does not refer to the base curve/contact')
    return delivery, pair, sites


def potential(r, w, w0, shape):
    after = PchipInterpolator(r, w, extrapolate=False)
    a = r[0]
    chord = (w[0] - w0)/a
    if shape == 'linear_W':
        inside = [0., 0., chord, w0]
    else:
        m0 = chord if shape == 'chord_hermite' else 0.
        inside = CubicHermiteSpline([0., a], [w0, w[0]], [m0, float(after.derivative()(a))]).c[:, 0]
    return PPoly(np.column_stack((inside, after.c)), np.r_[0., r], extrapolate=False)


def maximum(poly, lo, hi):
    roots = poly.deriv().roots()
    points = [lo, hi, *(float(t.real) for t in roots if abs(t.imag) < 1e-12 and lo < t.real < hi)]
    return float(np.max(poly(np.array(points))))


def prepare(w):
    # For E>=50 eV: every interval capable of a forbidden region is strictly
    # decreasing in total V, or the conservative upper bound proves it allowed.
    for a, b, coeffs in zip(w.x[:-1], w.x[1:], w.c.T):
        lo = max(a, R0)
        poly = Polynomial(coeffs[::-1])
        if C/lo + maximum(poly, lo-a, b-a) >= FLOOR:
            if maximum(poly.deriv(), lo-a, b-a) >= C/b**2:
                raise ValueError('cannot establish unique supported turning point')
    return w


def turn(w, e):
    end = w.x[-1]
    v = lambda r: C/r + float(w(r)) - e
    if e < FLOOR or C/end > e or v(end) > 0 or v(R0) <= 0:
        return None
    return brentq(v, R0, end, xtol=1e-15, full_output=False)


def gain(w, e):
    ts = turn(w, e)
    if ts is None:
        return None
    tb = C/e
    if tb <= R0:
        return None
    end = max(ts, tb)

    def integrand(t):
        r = np.exp(t)
        a = C/r - e
        correction = float(w(r))
        b = a + correction
        if a > 0 and b > 0:
            return -r*correction/(sqrt(a) + sqrt(b))
        return r*(sqrt(max(a, 0.)) - sqrt(max(b, 0.)))

    splits = sorted({log(x) for x in [ts, tb, *w.x] if R0 < x < end})
    value, err = quad(integrand, log(R0), log(end), points=splits, epsabs=1e-11, epsrel=2e-10, limit=80)
    return dict(G=K*value, R_turn_site_A=ts, R_turn_bare_A=tb, quad_estimate=K*err)


def eq7(e, g):
    if g is None:
        return {'status': 'unsupported', 'U_eV': None}
    e /= 1000.
    x = BG/sqrt(e)
    lower = -log(x)
    f = lambda ell: ell - x*expm1(ell)
    if g > f(lower):
        return {'status': 'no_increasing_branch_root', 'U_eV': None}
    upper = 1.
    while f(upper) > g:
        upper *= 2
    ell = brentq(lambda ell: f(ell)-g, lower, upper, xtol=1e-15, full_output=False)
    return {'status': 'suppression' if g < 0 else 'enhancement' if g > 0 else 'bare',
            'U_eV': 1000*e*expm1(-2*ell) if g else 0.}


def coulomb_action(e):
    x = R0*e/C
    return K*C/sqrt(e)*(acos(sqrt(x)) - sqrt(x*(1-x)))


def check(r, sites, rows):
    seen = {'bare_max_abs_G': 0., 'constant_shift_max_abs_G_delta': 0., 'eq7_roundtrip_max_abs_log_delta': 0.}
    for e in (50., 250., 1250., 6500., 10000.):
        for u in (0., 40., -20.):
            w = prepare(potential(r, np.full_like(r, -u), -u, 'linear_W'))
            readout = gain(w, e)
            assert readout is not None
            measured = readout['G']
            oracle = coulomb_action(e)-coulomb_action(e+u)
            key = 'constant_shift_max_abs_G_delta' if u else 'bare_max_abs_G'
            seen[key] = max(seen[key], abs(measured-oracle))
            for ch in ('p', 'n'):
                ek = e/1000.
                ue = eq7(e, measured)['U_eV']/1000.
                g2 = log(float(reaction_sigma(ek, ue, ch)/reaction_sigma(ek, 0., ch)))
                seen['eq7_roundtrip_max_abs_log_delta'] = max(seen['eq7_roundtrip_max_abs_log_delta'], abs(measured-g2))
    try:
        contact_bound(sites['oxygen'][2], sites['plain'][3], sites['plain'][4])
    except ValueError:
        seen['swapped_real_contact_refused'] = True
    assert seen.get('swapped_real_contact_refused')
    assert gain(w, FLOOR-1.) is None and eq7(FLOOR, None)['status'] == 'unsupported'
    assert eq7(FLOOR, 1000.)['status'] == 'no_increasing_branch_root' and eq7(FLOOR, -.5)['status'] == 'suppression'
    assert max(seen[key] for key in seen if key.endswith(('G', 'delta'))) < 1e-8
    assert all(isfinite(v['delta_G_O_minus_clean']) for v in rows)
    return seen


def evaluate(root):
    delivery, pair, sites = load(root)
    energies = sorted({*np.geomspace(.05, 6.5, 121).tolist(), *ANCHORS})
    prepared = {name: {s: prepare(potential(r, w, c['signedW0_eV'], s)) for s in SHAPES}
                for name, (r, w, c, _, _) in sites.items()}
    rows = []
    for x in energies:
        e = 1000*x
        for shape in SHAPES:
            g = [gain(prepared[n][shape], e) for n in sites]
            results = [dict(v or {}, **eq7(e, v['G'] if v else None)) for v in g]
            clean, oxy = results
            rows.append(dict(E_cm_keV=x, shape=shape, plain=clean, oxygen=oxy,
                             delta_G_O_minus_clean=oxy['G']-clean['G'] if all(g) else None,
                             delta_U_O_minus_clean_eV=oxy['U_eV']-clean['U_eV']
                             if all(r['U_eV'] is not None for r in results) else None))
    guards = check(sites['plain'][0], sites, rows)
    return dict(schema='zr-DD-site-action-readout-v1', base=BASE, draft=delivery['draft_id'],
                map_sha256=MAP_SHA, archived_source_checkpoint=delivery['source_checkpoint'],
                site_bindings={n: dict(curve_sha256=digest(raw), reference=json.loads(raw)['reference'],
                                       contact_sha256=digest((root/f'{n}_q1_contact.json').read_bytes()),
                                       W0_eV=c['signedW0_eV'], R_A=r.tolist(), W_eV=w.tolist(),
                                       ray=c['ray'], target_A=c['geometry']['positions_A'][c['geometry']['target_index']],
                                       fmax_eV_A=c['geometry']['raw_fmax_eV_A'], cell_A=c['geometry']['cell_A'])
                               for n,(r,w,c,raw,_) in sites.items()},
                endpoint_controls=pair, anchors_cm_keV=list(ANCHORS), finite_Rmin_A=R0,
                isotope='DD', rho_convention='same finite 4 fm action; Eq7 S(E) held fixed',
                sensitivity_scope='continuation shapes only; no action-level h/k/cell bound; 10keV reference separate',
                guards=guards, rows=rows)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs', type=Path)
    p.add_argument('result', type=Path)
    args = p.parse_args()
    args.result.write_text(json.dumps(evaluate(args.inputs), indent=2, allow_nan=False)+'\n')
