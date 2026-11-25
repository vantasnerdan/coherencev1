"""PdD .80 Å radial update; the accepted .48 cubic angular closure is unchanged."""
from pathlib import Path
import importlib.util
import json
import types
import numpy as np
from scipy.interpolate import CubicHermiteSpline

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / 'bulk-refit-048'
INPUT = HERE.parent / '.input-080/g10_quantum_site/bulk-wide-wells'

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

prior = load_module('pdd048_well', PRIOR / 'well.py')
K, A, B, K4 = prior.K, prior.A, prior.B, prior.K4
RAYS, VECTORS, R0, RMAX = prior.RAYS, prior.VECTORS, prior.R0, .80
VARIANTS = ('wide_quadratic', 'wide_quartic')
WIDE = json.loads((INPUT / 'bulk-wells.json').read_text())
ANCHORS = []
SPLINES = []
for j, axis in enumerate(RAYS):
    native = [p for p in WIDE['rows'] if p['material'] == 'PdD' and p['direction'] == axis
              and p['displacement_A'] == .8 and p['sample_kind'] == 'native_scf']
    if len(native) != 1:
        raise ValueError(f'Expected one independently native +.80 {axis} row')
    p = native[0]
    points = prior.ANCHORS[j] + [dict(r_A=.8, energy_eV=p['delta_energy_free_eV'],
              slope_eV_A=-p['native_force_parallel_eV_A'], native_cases=[p['case']], independent_native_count=1)]
    ANCHORS.append(points)
    r = [R0] + [p['r_A'] for p in points]
    e = [.5*K*R0**2+K4[j]*R0**4/4] + [p['energy_eV'] for p in points]
    s = [K*R0+K4[j]*R0**3] + [p['slope_eV_A'] for p in points]
    SPLINES.append(CubicHermiteSpline(r, e, s, extrapolate=False))

def radial(r, variant):
    if variant in prior.VARIANTS:
        return prior.radial(r, variant)
    if variant not in VARIANTS:
        raise ValueError(f'Unknown variant {variant!r}')
    r = np.asarray(r, float)
    E, S = np.empty(r.shape+(3,)), np.empty(r.shape+(3,))
    for j, spline in enumerate(SPLINES):
        rc = np.clip(r, R0, RMAX)
        e, s = spline(rc), spline(rc, 1)
        e = np.where(r<R0, .5*K*r*r+K4[j]*r**4/4, e)
        s = np.where(r<R0, K*r+K4[j]*r**3, s)
        t = np.maximum(r-RMAX, 0.)
        er, sr = ANCHORS[j][-1]['energy_eV'], ANCHORS[j][-1]['slope_eV_A']
        if variant == 'wide_quadratic':
            ee, ss = er+sr*t+.5*K*t*t, sr+K*t
        else:
            c = K+3*K4[j]*RMAX**2
            ee = er+sr*t+.5*c*t*t+K4[j]*RMAX*t**3+K4[j]*t**4/4
            ss = sr+c*t+3*K4[j]*RMAX*t*t+K4[j]*t**3
        E[..., j], S[..., j] = np.where(r>RMAX, ee, e), np.where(r>RMAX, ss, s)
    return E, S

# Rebind only the radial callback in the accepted analytic angular operator.
# No second implementation of its angular gradient or origin polynomial.
energy_force = types.FunctionType(prior.energy_force.__code__,
    dict(prior.__dict__, radial=radial), 'energy_force')
