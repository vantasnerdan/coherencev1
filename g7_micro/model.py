"""Site-population averages on the existing frozen CSDA and nuclear kernels.

U is keV; f is occupied fuel-number fraction, not a yield or host-site fraction.
The explicit sheet convention is shared (same law) or bulk40 (40 eV sheet).
"""
from functools import lru_cache

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre

from g4_beam.forward.consumer import QE, csda_nodes
from g4_beam.forward.profile_thermal import transport_stopping
from g6_plateau.forecast import cross_section, dt_stopping, GAUSS_X, GAUSS_W

BULK = .040
U_GRID = np.unique(np.r_[0., BULK, np.geomspace(1e-5, 1e9, 500)])
LOG_GRID = np.log(U_GRID+.005)


@lru_cache(maxsize=4)
def tail_rule(order=96):
    """Log-x exponential rule; omitted mass <1.000003e-10 on [1e-10,36]."""
    gx, gw = roots_legendre(order)
    lo, hi = np.log(1e-10), np.log(36.)
    x = np.exp((lo+hi)/2+(hi-lo)*gx/2)
    return x, gw*(hi-lo)/2*x*np.exp(-x)


def tail_nodes(scale, order=96):
    """U=BULK+scale*x; finite quadrature, not sigma at the mean shift."""
    x, weight = tail_rule(order)
    return BULK+scale*x, weight


def population(family, a, b, order=96):
    if family == 'two_site':
        return np.array([BULK, b]), np.array([1-a, a])
    if family == 'exponential':
        u, w = tail_nodes(b, order)
        return np.r_[BULK, u], np.r_[1-a, a*w]
    raise ValueError('A physical population is two_site or exponential')


class FrozenResponse:
    """Compile only authoritative stopping nodes; interpolate single-site yields.

    The sheet remains separate. No stopping or source profile is recomputed when
    screening changes. Tables are an accelerator; exact=True bypasses them.
    """
    def __init__(self, nominal_dd_cm, profile, *, isotope='DD', measurement=None):
        self.e = float(nominal_dd_cm)
        self.isotope = isotope
        self.measurement = dict(measurement or {})
        self.profile = profile
        layers, self.nads, self.state = profile
        layers = np.asarray(layers, float)
        if (layers.ndim != 2 or layers.shape[1] != 7 or not len(layers) or
                np.any(~np.isfinite(layers)) or not np.isfinite(self.nads) or self.nads < 0):
            raise ValueError('Seven-field finite SI layers and nonnegative separate sheet required')
        if (np.any(layers[:, 1] <= layers[:, 0]) or np.any(layers[:, 2:4] < 0) or
                abs(layers[0, 0]) > 1e-15 or
                np.any(np.abs(layers[1:, 0]-layers[:-1, 1]) > 1e-15)):
            raise ValueError('Contiguous nonnegative-density layers must start at zero')
        self.layers = tuple(tuple(float(v) for v in row) for row in layers)
        f2 = self.measurement.get('D2_fraction', .15)
        f1 = self.measurement.get('D1_fraction', 1-f2)
        spread = self.measurement.get('spread_FWHM', .05)/2.3548
        if (not np.isfinite([f1, f2, spread]).all() or min(f1, f2) < 0 or
                abs(f1+f2-1) > 1e-12 or spread < 0):
            raise ValueError('Explicit current fractions must sum to one; nonnegative spread')
        stop = dt_stopping if isotope == 'DT' else transport_stopping
        self.tracks = []
        for molecular, share in ((1, f1), (2, 2*f2)):
            for x, w in zip(GAUSS_X, GAUSS_W):
                lab = 2*self.e*(1+spread*x)/molecular
                en, _, weights, _ = csda_nodes(lab, self.layers, stop)
                self.tracks.append((lab, share*w, en, weights))
        self.tables = {}
        self.population_tables = {}

    def _exact_track(self, track, u, channel):
        lab, _, en, weights = track
        u = np.atleast_1d(np.asarray(u, float))
        volume = np.zeros(len(u))
        for start in range(0, len(u), 64):
            part = u[start:start+64]
            volume[start:start+64] = weights@cross_section(
                en[:, None], part[None, :], self.isotope, 'proton', channel)
        sheet = self.nads*1e-31*cross_section(lab, u, self.isotope, 'proton', channel)
        return np.array([volume, sheet])

    def _table(self, channel):
        if channel not in self.tables:
            table = []
            total = np.zeros((2, len(U_GRID)))
            for track in self.tracks:
                values = self._exact_track(track, U_GRID, channel)
                total += track[1]*values
                table.append([None if np.all(y == 0) else
                              PchipInterpolator(LOG_GRID, np.log(y), extrapolate=False)
                              for y in values])
            self.tables[channel] = table
            self.population_tables[channel] = [
                None if np.all(y == 0) else PchipInterpolator(LOG_GRID, np.log(y), extrapolate=False)
                for y in total]
        return self.tables[channel]

    def components(self, u, channel='p', *, exact=False) -> np.ndarray:
        u = np.atleast_1d(np.asarray(u, float))
        if np.any(~np.isfinite(u)) or np.any(u < 0) or np.any(u > U_GRID[-1]):
            raise ValueError('U must be finite in [0, 1e9] keV')
        if exact:
            result = np.zeros((2, len(u)))
            for track in self.tracks:
                result += track[1]*self._exact_track(track, u, channel)
            return result
        self._table(channel)
        return np.array([np.zeros(len(u)) if fn is None else np.exp(fn(np.log(u+.005)))
                         for fn in self.population_tables[channel]])

    def yield_per_C(self, family, a, b, *, channel='p', sheet='shared', exact=False,
                    order=96, moment=False, E0=None):
        if sheet not in ('shared', 'bulk40'):
            raise ValueError('Explicit shared or bulk40 sheet convention required')
        if family in ('two_site', 'exponential'):
            if not 0 <= a <= 1 or (family == 'two_site' and b < BULK) or b <= 0:
                raise ValueError('Invalid population weight or screening scale')
            u, probability = population(family, a, b, order)
            values = self.components(u, channel, exact=exact)
            result = values@(probability*u if moment else probability)
        elif family in ('sigmoid', 'anchored_rise'):
            ub, us = (a, b) if family == 'sigmoid' else (BULK, BULK+a)
            knee = E0 if family == 'sigmoid' else b
            result = np.zeros(2)
            for i, track in enumerate(self.tracks):
                u = ub+(us-ub)/(1+np.exp((track[0]/2-knee)/.5))
                if exact:
                    values = self._exact_track(track, [u], channel)[:, 0]
                else:
                    values = np.array([0. if fn is None else float(np.exp(fn(np.log(u+.005))))
                                       for fn in self._table(channel)[i]])
                result += track[1]*values*(u if moment else 1.)
        else:
            raise ValueError('Unknown screening family')
        if sheet == 'bulk40':
            result[1] = self.components([BULK], channel, exact=exact)[1, 0]*(BULK if moment else 1.)
        return float(np.sum(result)/QE)


def raw_yield(nominal_dd_cm, profiles, family, a, b, *, isotope='DD', channel='p',
              sheet='shared', measurement=None, order=96, E0=None):
    """Raw true reactions/C at an actual exact-key source profile; no normalization.

    nominal_dd_cm is lab/2, including for DT: never a DT-CM profile key.
    No empirical stitch, detector efficiency or cup-current correction is applied.
    """
    profile = profiles(float(nominal_dd_cm)) if callable(profiles) else profiles[float(nominal_dd_cm)]
    response = FrozenResponse(nominal_dd_cm, profile, isotope=isotope, measurement=measurement)
    return response.yield_per_C(family, a, b, channel=channel, sheet=sheet, exact=True,
                                order=order, E0=E0)


def depth_components(response, d_nm, shifts, channel='p', *, moment=False):
    """Baseline and linear surface/bulk increments on exactly split CSDA nodes.

    Return arrays [U, component], component=(baseline, top increment,
    deep increment, sheet increment), in reactions/true C. The sheet is not a
    manufactured layer; shared uses f_s there and bulk40 excludes its increment.
    moment=True inserts local U in keV before averaging, giving keV*reactions/C.
    """
    shifts = np.atleast_1d(np.asarray(shifts, float))
    result = np.zeros((len(shifts), 4))
    boundary = float(d_nm)*1e-9
    stop = dt_stopping if response.isotope == 'DT' else transport_stopping
    for lab, share, _, _ in response.tracks:
        en, z, weights, _ = csda_nodes(lab, response.layers, stop, break_m=boundary)
        base = cross_section(en, BULK, response.isotope, 'proton', channel)
        high = cross_section(en[:, None], shifts[None, :], response.isotope, 'proton', channel)
        if moment:
            base = base*BULK
            high = high*shifts[None, :]
        delta = high-base[:, None]
        top = z < boundary
        sheet_base = response.nads*1e-31*cross_section(lab, BULK, response.isotope, 'proton', channel)
        if moment:
            sheet_base *= BULK
        result[:, 0] += share*(float(weights@base)+sheet_base)
        result[:, 1] += share*(weights[top]@delta[top])
        result[:, 2] += share*(weights[~top]@delta[~top])
        sheet_high = cross_section(lab, shifts, response.isotope, 'proton', channel)
        sheet_low = cross_section(lab, BULK, response.isotope, 'proton', channel)
        if moment:
            sheet_high *= shifts
            sheet_low *= BULK
        result[:, 3] += share*response.nads*1e-31*(sheet_high-sheet_low)
    return result/QE


def raw_depth_yield(nominal_dd_cm, profile, *, f_s, f_b, d_nm, U_high_keV,
                    isotope='DD', channel='p', sheet='shared', measurement=None, moment=False):
    """Return Y/C, or its unnormalized U-weighted moment when explicitly requested."""
    if not 0 <= f_b <= f_s <= 1 or d_nm <= 0 or U_high_keV < BULK:
        raise ValueError('Require 0 <= f_b <= f_s <= 1, positive depth and U_high >= 40 eV')
    if sheet not in ('shared', 'bulk40'):
        raise ValueError('Explicit sheet convention required')
    response = FrozenResponse(nominal_dd_cm, profile, isotope=isotope, measurement=measurement)
    base, top, deep, ads = depth_components(response, d_nm, [U_high_keV], channel, moment=moment)[0]
    return float(base+f_s*(top+(ads if sheet == 'shared' else 0.))+f_b*deep)
