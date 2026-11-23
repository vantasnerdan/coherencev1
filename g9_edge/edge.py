"""Finite causal Pd-proxy entry model and immutable transported-profile consumer.

Units are eV, Angstrom and fs except at the existing SI profile boundary.
This is an analysis replay, not a first-principles surface calculation.
"""
import argparse
import csv
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import j0, j1

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from g8_why.common import (DD, ENERGIES, MATERIALS, GX, GW, QE, U,
                          reference, yield_C, source_hashes)
from g4_beam.forward.consumer import csda_nodes
from g4_beam.forward.profile_thermal import transport_stopping

HBAR = .6582119569509067
C = 14.399645351950548
AK = 3.80998212                 # hbar^2/(2 m_e), eV A^2
BOHR = .529177210903
MD = 208.7331208               # deuteron mass, eV fs^2/A^2
RADII = (0., .1, .5, 1., 2.)
MATCH = 1.                     # kinetic / inner-barrier accounting boundary, A
GAMMA = .5                     # conditional Mermin collision width, eV
REG = .1                       # extra spectral realization width, eV; not measured loss


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def cum(y, t):
    return cumulative_trapezoid(y, t, axis=-1, initial=0.)


def gauss_panels(edges, order):
    x, w = np.polynomial.legendre.leggauss(order)
    half = np.diff(edges)/2
    mid = (np.asarray(edges)[1:]+np.asarray(edges)[:-1])/2
    return (mid[:, None]+half[:, None]*x).ravel(), (half[:, None]*w).ravel()


def dsinc(q, r):
    """d sinc(q r)/dr, with the removable origin evaluated analytically."""
    x = q*r
    out = np.empty(np.broadcast_shapes(np.shape(q), np.shape(r)))
    small = abs(x) < 1e-3
    np.divide(x*np.cos(x)-np.sin(x), x*x, out=out, where=~small)
    out[small] = (-x/3+x**3/30-x**5/840)[small]
    return q*out


def loss_moments(energy, loss):
    """Exact E*L and L/E integrals for the ORIGINAL piecewise-linear spectrum."""
    a, b = energy[:-1], energy[1:]
    left = loss[..., :-1]
    slope = np.diff(loss, axis=-1)/(b-a)
    first = np.sum(left*(b*b-a*a)/2+slope*((b**3-a**3)/3-a*(b*b-a*a)/2), axis=-1)
    intercept = left-slope*a
    ratio = np.ones_like(a)
    np.divide(b, a, out=ratio, where=a > 0)
    inverse = np.sum(slope*(b-a)+intercept*np.log(ratio), axis=-1)
    return first, inverse


def lindhard(q, energy, kf):
    """Retarded T=0 electron gas, complex energy in the upper half-plane."""
    x = q/(2*kf)
    u = energy/(2*AK*kf*q)
    plus, minus = u+x, u-x
    f = (.5+(1-plus**2)/(8*x)*np.log((plus+1)/(plus-1))
         -(1-minus**2)/(8*x)*np.log((minus+1)/(minus-1)))
    # Avoid cancellation of O(u) logarithmic terms at high frequency/small q.
    high = abs(u) > 20+abs(x)
    asymptotic = (-1/(3*u*u)-(1/5+x*x/3)/u**4
                  -(1/7+2*x*x/3+x**4/3)/u**6)
    f = np.where(high, asymptotic, f)
    return 1+4*kf/(np.pi*BOHR*q*q)*f


def static_lindhard(q, kf):
    x = q/(2*kf)
    f = .5+(1-x*x)/(4*x)*np.log(abs((1+x)/(1-x)))
    return 1+4*kf/(np.pi*BOHR*q*q)*f


def mermin(q, energy, kf):
    a = lindhard(q, energy+1j*GAMMA, kf)-1
    a0 = static_lindhard(q, kf)-1
    return 1+(energy+1j*GAMMA)*a/(energy+1j*GAMMA*a/a0)


def reduce_sources(matrix_dir, out):
    """Use the saved FULL G inverse, then take its G=0 head (not 1/eps00)."""
    completed = json.loads((ROOT/'g6_response/data/Pd/completed.json').read_text())
    selected = [r for r in completed['points']
                if r['settings']['group'] in ('primary-N3', 'lambda-G')]
    samples, rows, identities = [], [], []
    for record in selected:
        path = matrix_dir/record['arrays']['path']
        if digest(path) != record['arrays']['sha256']:
            raise ValueError('Matrix source identity changed: '+path.name)
        with np.load(path, allow_pickle=False) as data:
            energy = data['energy_eV'].copy()
            inverse = data['inverse_epsilon_coulomb_symmetrized']
            epsilon = data['epsilon_coulomb_symmetrized']
            q = float(np.linalg.norm(data['Q_Ainv'][0]))
            head = inverse[:, 0, 0].copy()
            # This residual contracts the actual full G row/column before reduction.
            residual = float(np.max(abs(np.einsum('wi,wi->w', epsilon[:, 0, :],
                                                   inverse[:, :, 0])-1)))
            scalar_wrong = 1/epsilon[:, 0, 0]
            positive = energy > 0
            loss = -head.imag
            samples.append((q, energy, loss, head))
            for i, en in enumerate(energy):
                rows.append(dict(file=path.name, group=record['settings']['group'],
                    q_Ainv=q, energy_eV=float(en), inverse_head_real=float(head[i].real),
                    loss_raw=float(loss[i]), scalar_inverse_real=float(scalar_wrong[i].real),
                    scalar_inverse_loss=float(-scalar_wrong[i].imag),
                    full_G_head_identity_residual=residual))
            identities.append(dict(file=path.name, sha256=digest(path), q_Ainv=q,
                energy_grid_eV=energy.tolist(), ecut_eV=200, kmesh=12, nbands=96,
                eta_eV=.1, full_G_head_identity_residual=residual,
                minimum_positive_frequency_loss=float(loss[positive].min()),
                zero_frequency_imaginary_artifact=float(loss[0])))
    samples.sort(key=lambda r: r[0])
    energy = samples[0][1]
    if any(not np.array_equal(r[1], energy) for r in samples):
        raise ValueError('Mixed actual energy grids')
    # Isotropize equal-|q| directional heads only after the G inversion.
    qs, losses, statics = [], [], []
    for q, _, loss, head in samples:
        if qs and abs(q-qs[-1]) < 1e-7:
            losses[-1].append(loss)
            statics[-1].append(float(head[0].real))
        else:
            qs.append(q); losses.append([loss]); statics.append([float(head[0].real)])
    losses = np.array([np.mean(a, axis=0) for a in losses])
    statics = np.array([np.mean(a) for a in statics])
    # Preserve negative losses. A real negative positive-frequency specimen fails visibly.
    if losses[:, 1:].min() < -1e-8:
        raise ValueError('Negative physical-frequency loss retained in reduction.csv')
    physical_loss = losses.copy()
    physical_loss[:, 0] = 0.  # raw static imaginary artifacts stay in reduction.csv
    first, _ = loss_moments(energy, physical_loss)
    moments = 2/np.pi*first
    wp = float(np.sqrt(moments[0]))
    density = wp*wp/(8*np.pi*C*AK)
    kf = float((3*np.pi*np.pi*density)**(1/3))
    write_csv(out/'reduction.csv', rows)
    return dict(q=np.array(qs), energy=energy, loss=physical_loss, static=statics,
                wp=wp, kf=kf, density=density, moments=moments, sources=identities)


class Spectrum:
    def __init__(self, source, dt=.0005):
        self.source = source
        self.t = np.arange(round(.8/dt)+1)*dt
        self.q, self.wq = gauss_panels([0., .25, .5, 1., 2., 4., 8., 12., 16.], 16)
        e1, w1 = gauss_panels(source['energy'], 8)
        e2, w2 = gauss_panels(np.arange(96., 2048.01, 2.), 2)
        self.e = np.r_[e1, e2]
        self.we = np.r_[w1, w2]
        loss = -np.imag(1/mermin(self.q[:, None], self.e[None, :], source['kf']))
        covered_q = (self.q >= source['q'][0]) & (self.q <= source['q'][-1])
        covered_e = self.e <= 96
        raw_on_e = np.array([np.interp(self.e[covered_e], source['energy'], row)
                             for row in source['loss']])
        for i in np.flatnonzero(covered_q):
            loss[i, covered_e] = np.array([np.interp(self.q[i], source['q'], raw_on_e[:, j])
                                          for j in range(raw_on_e.shape[1])])
        self.loss = loss
        if not np.isfinite(loss).all() or loss.min() < -1e-7:
            raise ValueError('Extension loss is nonpassive; no clipping is used')
        self.g0 = -2/np.pi*((loss*self.we)@(self.e/(self.e**2+REG**2)))
        self.h = np.empty((len(self.q), len(self.t)))
        # Bound temporary sin-array RSS; one shared spectral realization for all speeds.
        for a in range(0, len(self.t), 128):
            t = self.t[a:a+128]
            self.h[:, a:a+128] = (-2/(np.pi*HBAR)*(loss*self.we)@
                np.sin(self.e[:, None]*t[None, :]/HBAR))*np.exp(-REG*t/HBAR)
        self.H = cum(self.h, self.t)
        wp = source['wp']
        def oscillator(w):
            om = np.sqrt(w*w-GAMMA*GAMMA/4)/HBAR
            return -w*w/(HBAR*HBAR*om)*np.exp(-GAMMA*self.t/(2*HBAR))*np.sin(om*self.t)
        # Exact local long-wavelength Mermin/Drude image: B=-g_bulk+g_surface.
        self.himage = -oscillator(wp)+oscillator(wp/np.sqrt(2))
        self.Himage = cum(self.himage, self.t)
        self.image_static = 0.  # B(0)=0 analytically for the conducting local limit
        self.frozen_raw = np.interp(source['q'], self.q, 1+self.g0)
        self._doppler = {}

    def doppler(self, v):
        """Cache the unchanged spectral operator, not a fitted stopping curve."""
        if v not in self._doppler:
            mu, wm = np.polynomial.legendre.leggauss(48)
            alpha = HBAR*self.q[:, None]*mu[None, :]*v+1j*REG
            g = np.empty(alpha.shape, complex)
            for i in range(len(self.q)):
                g[i] = 2/np.pi*np.sum((self.we*self.loss[i]*self.e)[None, :]/
                    (alpha[i, :, None]**2-self.e[None, :]**2), axis=1)
            self._doppler[v] = mu, wm, g
        return self._doppler[v]

    def moving_wake(self, v, radii):
        """Infinite-prehistory BULK reference of this same damped realization."""
        mu, wm, g = self.doppler(v)
        return {r: float(2*C/np.pi*np.sum(self.wq[:, None]*wm[None, :]/2*
                     np.real(np.exp(1j*self.q[:, None]*mu[None, :]*r)*g))) for r in radii}


class History:
    def __init__(self, spectrum, ecm, kcap=1., qcap=16.):
        s = spectrum
        self.s, self.t, self.ecm = s, s.t, float(ecm)
        self.v = float(np.sqrt(4*1000*ecm/MD))
        self.wq = s.wq*(s.q <= qcap)
        self.k, self.wk = gauss_panels([0., kcap], 64)
        dt = self.t[1]-self.t[0]
        decay = np.exp(-self.k*self.v*dt)
        self.L = np.zeros((len(self.k), len(self.t)))
        # Stable trapezoidal recurrence for integral h_B(u) exp[-k v(t-u)] du.
        for j in range(1, len(self.t)):
            self.L[:, j] = (decay*self.L[:, j-1]+dt/2*
                            (decay*s.himage[j-1]+s.himage[j]))
        r = self.v*self.t
        force_kernel = s.h*dsinc(s.q[:, None], r[None, :])
        self.Fp_bulk = -2*C/np.pi*(self.wq@cum(force_kernel, self.t))
        self.Fp_image = C*((self.wk*self.k)[:, None]*np.exp(-self.k[:, None]*r)*self.L).sum(axis=0)
        mu, wm, g = s.doppler(self.v)
        self.formed_bulk_drag = float(-2*C/np.pi*np.sum(self.wq[:, None]*wm[None, :]/2*
            np.real(1j*s.q[:, None]*mu[None, :]*g)))
        # CSDA ALREADY includes settled drag. Insert only the arrival deficit /
        # excess relative to the formed moving bulk wake, not that drag again.
        self.pp_bulk = cum(self.Fp_bulk-self.formed_bulk_drag, self.t)
        self.pp_image = cum(self.Fp_image, self.t)
        # Impulse to a target held at Z=MATCH+v*t until the matching encounter.
        a = cum(s.h*np.sinc(s.q[:, None]*(MATCH+r)[None, :]/np.pi), self.t)
        causal = 2*C/(np.pi*self.v)*(self.wq@(a-np.sinc(s.q[:, None]*(MATCH+r)[None, :]/np.pi)*s.H))
        frozen = 2*C/(np.pi*self.v)*((self.wq*s.g0)@
            (np.sinc(s.q[:, None]*MATCH/np.pi)-np.sinc(s.q[:, None]*(MATCH+r)[None, :]/np.pi)))
        self.pt_bulk = causal-frozen
        self.pt_image = C/self.v*((self.wk[:, None]*np.exp(-self.k[:, None]*(MATCH+r))*
                                   (s.Himage[None, :]-self.L))).sum(axis=0)

    def potential(self, radius, surface=False):
        s = self.s
        kernel = s.h*np.sinc(s.q[:, None]*(radius+self.v*self.t)[None, :]/np.pi)
        actual = 2*C/np.pi*(self.wq@cum(kernel, self.t))
        frozen = float(2*C/np.pi*np.sum(self.wq*s.g0*np.sinc(s.q*radius/np.pi)))
        if surface:
            actual = actual+C*(self.wk[:, None]*np.exp(-self.k[:, None]*(radius+self.v*self.t))*self.L).sum(axis=0)
        return actual, frozen

    def kinetic(self, depth_A, surface=True):
        tm = max((depth_A-MATCH)/self.v, 0.)
        pp = np.interp(tm, self.t, self.pp_bulk+(self.pp_image if surface else 0.))
        pt = np.interp(tm, self.t, self.pt_bulk+(self.pt_image if surface else 0.))
        return MD/4*((self.v+(pp-pt)/MD)**2-self.v*self.v), pp, pt

    def fixed_sites(self, end):
        """Differentiate at fixed (rho,Z), integrate in LAB time, never along R=const."""
        s = self.s
        stride = max(1, round(.002/(self.t[1]-self.t[0])))
        indices = np.unique(np.r_[np.arange(0, min(len(self.t), np.searchsorted(self.t, end)+1), stride),
                                   min(len(self.t)-1, np.searchsorted(self.t, end))])
        tt = self.t[indices]
        a = 4.084
        sites = [('target', 0., a/2, 1), ('same_plane_4', a/np.sqrt(2), a/2, 4),
                 ('face_4', a/2, 0., 4), ('next_plane_4', a/2, a, 4)]
        rows = []
        for name, rho, z, multiplicity in sites:
            fz, fr = [], []
            for it in indices:
                lag = self.t[:it+1]
                dz = z-self.v*(self.t[it]-lag)
                rr = np.hypot(rho, dz)
                deriv = dsinc(s.q[:, None], rr[None, :])
                direction_z = np.divide(dz, rr, out=np.zeros_like(rr), where=rr != 0)
                direction_r = np.divide(rho, rr, out=np.zeros_like(rr), where=rr != 0)
                hz = (self.wq[:, None]*s.h[:, :it+1]*deriv).sum(axis=0)
                bulkz = -2*C/np.pi*np.trapezoid(hz*direction_z, lag)
                bulkr = -2*C/np.pi*np.trapezoid(hz*direction_r, lag)
                dz0 = z-self.v*self.t[it]
                r0 = np.hypot(rho, dz0)
                radial0 = float(2*C/np.pi*np.sum(self.wq*s.g0*dsinc(s.q, r0)))
                refz = -radial0*dz0/r0 if r0 else 0.
                refr = -radial0*rho/r0 if r0 else 0.
                imagez = C*np.sum(self.wk*self.k*np.exp(-self.k*z)*j0(self.k*rho)*self.L[:, it])
                imager = C*np.sum(self.wk*self.k*np.exp(-self.k*z)*j1(self.k*rho)*self.L[:, it])
                fz.append(bulkz+imagez-refz)
                fr.append(bulkr+imager-refr)
            pz, pr = cum(np.asarray(fz), tt), cum(np.asarray(fr), tt)
            pp = np.interp(tt, self.t, self.pp_bulk+self.pp_image)
            for j, t in enumerate(tt):
                erel = MD/4*((self.v+(pp[j]-pz[j])/MD)**2+(pr[j]/MD)**2-self.v*self.v)
                rows.append(dict(E_cm_keV=self.ecm, site=name, multiplicity=multiplicity,
                    rho_A=rho, fixed_z_A=z, time_fs=float(t),
                    delta_Fz_eV_A=float(fz[j]), delta_Frho_eV_A=float(fr[j]),
                    delta_pz_eV_fs_A=float(pz[j]), delta_prho_eV_fs_A=float(pr[j]),
                    projectile_self_delta_pz_eV_fs_A=float(pp[j]),
                    kick_energy_eV=float((pz[j]**2+pr[j]**2)/(2*MD)),
                    pair_relative_energy_eV=float(1000*self.ecm+erel),
                    pair_delta_kinetic_eV=float(erel),
                    conditional_proton_line_delta_eV_view_z=float(
                        np.sqrt(2*(MD/2)*3.02e6)*(pp[j]+pz[j])/(2*MD))))
        return rows


def lobes(t, values):
    """Initial deficit plus first opposite-sign lobe; interpolated zero crossings."""
    sign = np.sign(values[0])
    crossing = np.flatnonzero(values[:-1]*values[1:] < 0)+1
    ends = []
    for i in crossing[:2]:
        ends.append(float(t[i-1]-values[i-1]*(t[i]-t[i-1])/(values[i]-values[i-1])))
    first = ends[0] if ends else float(t[-1])
    final = ends[1] if len(ends)>1 else float(t[-1])
    initial_mask = t <= first
    later_mask = (t > first) & (t <= final)
    i0 = np.flatnonzero(initial_mask)[np.argmax(abs(values[initial_mask]))]
    i1 = np.flatnonzero(later_mask)[np.argmax(abs(values[later_mask]))] if later_mask.any() else i0
    return dict(initial_peak_eV=float(values[i0]), initial_peak_time_fs=float(t[i0]),
                initial_lobe_duration_fs=first, later_peak_eV=float(values[i1]),
                later_peak_time_fs=float(t[i1]), later_lobe_duration_fs=final-first,
                edge_end_fs=final, two_zero_crossings_resolved=len(ends)>1,
                initial_sign=int(sign))


def rate_history(history, surface=True, turn_radius=None):
    # A common inner-turn constant shift, referenced at MATCH rather than infinity.
    turn = C/(1000*history.ecm) if turn_radius is None else turn_radius
    va, vf = history.potential(turn, surface)
    ma, mf = history.potential(MATCH, surface)
    barrier = -((va-ma)-(vf-mf))
    window = lobes(history.t, barrier)
    return barrier, window, turn


def rate_consumer(spectrum, kcap=1., qcap=16., formed_sheet=False):
    """Same layer/sheet, species, charge, spread, stopping and BH as g8.

    The first finite entry window is split EXACTLY in the existing depth nodes.
    Collision clock is z/v_entry: a declared first-spacings ballistic approximation;
    local collision energy is always the authoritative slowing-down energy.
    """
    rows = []
    for material in MATERIALS:
        for nominal in ENERGIES:
            layers, sheet, _ = reference(material).profiles(float(nominal))
            accum = {ch: np.zeros(6) for ch in ('p', 'n')}
            depths, times, signs = [], [], []
            for mult, frac in ((1, .85), (2, .15)):
                for g, wg in zip(GX, GW):
                    lab = 2*nominal/mult*(1+.05/2.3548*g)
                    history = History(spectrum, lab/2, kcap, qcap)
                    delta_u, window, turn = rate_history(history)
                    end = window['edge_end_fs']
                    edge_depth = history.v*end*1e-10
                    en, z, weight, _ = csda_nodes(float(lab), layers, transport_stopping, break_m=edge_depth)
                    charge = frac*mult*wg/QE
                    en = np.r_[en, lab]
                    z = np.r_[z, 0.]
                    weight = np.r_[weight, sheet*1e-31]*charge
                    clock = z*1e10/history.v
                    inside = clock <= end
                    # t=0 sheet is the gated atom source, not a preformed screening cloud.
                    du = np.interp(clock[inside], history.t, delta_u)
                    dk = np.array([history.kinetic(float(zz*1e10))[0] for zz in z[inside]])
                    if formed_sheet:
                        # A deliberately different prehistory ONLY for the exact sheet.
                        # This is a formed BULK wake sensitivity, not a solved vacuum approach.
                        wake = spectrum.moving_wake(history.v, (turn, MATCH))
                        _, vf = history.potential(turn)
                        _, mf = history.potential(MATCH)
                        du[z[inside] == 0.] = -(wake[turn]-wake[MATCH]-(vf-mf))
                    ec = np.maximum(en/2, .05)
                    # Boundaries are separate: pre-match force work enters kinetic;
                    # only potential DIFFERENCE inside MATCH enters the barrier.
                    ec_new = ec[inside]+dk/1000
                    u_new = U+du/1000
                    if np.any(ec_new <= 0) or np.any(ec_new+u_new <= 0):
                        raise ValueError('Consumer left positive-energy domain; no clipping')
                    for ch in ('p', 'n'):
                        base = DD.sigma(ec, U, ch)
                        transient = DD.sigma(ec_new, u_new, ch)
                        barrier_only = DD.sigma(ec[inside], u_new, ch)
                        kinetic_only = DD.sigma(ec_new, U, ch)
                        a = accum[ch]
                        a[0] += float(weight@base)
                        a[1] += float(weight[~inside]@base[~inside])
                        a[2] += float(weight[inside]@base[inside])
                        a[3] += float(weight[inside]@transient)
                        a[4] += float(weight[inside]@barrier_only)
                        a[5] += float(weight[inside]@kinetic_only)
                    depths.append(edge_depth*1e10); times.append(end); signs.append(float(du.min()))
            for ch, a in accum.items():
                total = a[1]+a[3]
                rows.append(dict(material=material, response_proxy='pure_Pd_RPA_head_plus_Mermin',
                    E_cm_keV=float(nominal), channel=ch, image_kcap_Ainv=kcap,
                    qcap_Ainv=qcap, smooth40_reactions_per_true_C=float(a[0]),
                    settled_outside_edge_reactions_per_true_C=float(a[1]),
                    baseline_inside_edge_reactions_per_true_C=float(a[2]),
                    transient_inside_edge_reactions_per_true_C=float(a[3]),
                    total_reactions_per_true_C=float(total),
                    edge_delta_reactions_per_true_C=float(a[3]-a[2]),
                    raw_yield_ratio=float(total/a[0]), baseline_edge_reaction_share=float(a[2]/a[0]),
                    actual_edge_reaction_share=float(a[3]/total),
                    barrier_only_total_ratio=float((a[1]+a[4])/a[0]),
                    kinetic_only_total_ratio=float((a[1]+a[5])/a[0]),
                    edge_depth_min_A=float(min(depths)), edge_depth_max_A=float(max(depths)),
                    edge_end_min_fs=float(min(times)), edge_end_max_fs=float(max(times)),
                    surface_sheet_source='formed_bulk_wake' if formed_sheet else 'gated_atom_at_t0',
                    local_clock='z/v_incident', Rturn_nominal_A=C/(1000*nominal), match_A=MATCH))
    return rows


def model_results(spectrum, out):
    curves, peaks, fixed, energy_rows = [], [], [], []
    for e in ENERGIES:
        h = History(spectrum, e)
        wake = spectrum.moving_wake(h.v, RADII)
        period = 2*np.pi*HBAR/spectrum.source['wp']
        du, window, turn = rate_history(h)
        energy_rows.append(dict(E_cm_keV=float(e), v_A_fs=h.v, switch_length_A=1.,
            switch_time_fs=1/h.v, moment_energy_eV=spectrum.source['wp'],
            moment_period_fs=period, switch_over_moment_period=1/(h.v*period),
            response_period_over_switch=h.v*period, Rturn_A=turn,
            **window))
        fixed.extend(h.fixed_sites(window['edge_end_fs']))
        for r in RADII:
            for surface in (False, True):
                actual, frozen = h.potential(r, surface)
                delta = actual-frozen
                info = lobes(h.t, delta)
                peaks.append(dict(E_cm_keV=float(e), geometry='entry' if surface else 'bulk',
                    R_ahead_A=r, frozen_adiabatic_eV=frozen,
                    formed_moving_bulk_eV=wake[r], **info))
                for j in range(0, len(h.t), 4):
                    curves.append(dict(E_cm_keV=float(e), geometry='entry' if surface else 'bulk',
                        time_fs=float(h.t[j]), ion_z_A=float(h.v*h.t[j]), R_ahead_A=r,
                        probe_z_A=float(h.v*h.t[j]+r), induced_qD_phi_eV=float(actual[j]),
                        frozen_adiabatic_qD_phi_eV=frozen,
                        delta_from_frozen_eV=float(delta[j]),
                        formed_moving_BULK_qD_phi_eV=wake[r]))
    write_csv(out/'curves.csv', curves)
    write_csv(out/'peaks.csv', peaks)
    write_csv(out/'fixed_sites.csv', fixed)
    write_csv(out/'edge_windows.csv', energy_rows)
    return energy_rows, peaks


def narrow_check(source, spectrum, primary, alternate_surface, alternate_tail, out):
    """One integrated changed-consumer check; no repo-wide validation."""
    # Actual source order: epsilon row times SAVED full inverse column, and wrong
    # scalar-before-inversion witness from reduction.csv, not a second G solve.
    reduction = list(csv.DictReader((out/'reduction.csv').open()))
    source_order = max(abs(float(r['inverse_head_real'])-float(r['scalar_inverse_real']))
                       for r in reduction if float(r['energy_eV']) > 0)
    baseline_error = max(abs(r['smooth40_reactions_per_true_C']/
        yield_C(r['material'], r['E_cm_keV'], r['channel'])-1) for r in primary)
    # t<0 is exactly zero by support, and an empty kernel makes both branches zero.
    causal_zero = float(np.max(abs(spectrum.h[:, 0])))
    h = History(spectrum, .25)
    # At zero speed, the causal convolution tends to its SAME spectral frozen
    # reference (analytic identity); numerically compare the modal static formula.
    modal_static = -2/np.pi*np.sum(spectrum.loss*spectrum.we[None, :]*
        spectrum.e[None, :]/(spectrum.e[None, :]**2+REG**2), axis=1)
    static_identity = float(np.max(abs(modal_static-spectrum.g0)))
    # Local dielectric limit of the stated specular-reflection formula, including
    # inside/vacuum boundary continuity. These are algebraic complex controls.
    eps = 2.3+.7j
    k, z, zp = .8, .4, .9
    es = lambda x: np.exp(-k*abs(x))/eps
    lhs = es(z-zp)+es(z+zp)-2*es(z)*es(zp)/(1+es(0))
    local = (np.exp(-k*abs(z-zp))+(eps-1)/(eps+1)*np.exp(-k*(z+zp)))/eps
    free = np.exp(-k*abs(z-zp))+np.exp(-k*(z+zp))-np.exp(-k*(z+zp))
    image_limit = float(abs(lhs-local))
    free_limit = float(abs(free-np.exp(-k*abs(z-zp))))
    continuity = float(abs(2*es(z)/(1+es(0))-2*np.exp(-k*z)/(eps+1)))
    # Same source/gate in entry and bulk: their difference is just the declared
    # image convolution. This rejects an entry-vs-bulk prehistory mismatch.
    va, _ = h.potential(.5, True); vb, _ = h.potential(.5, False)
    image = C*(h.wk[:, None]*np.exp(-h.k[:, None]*(.5+h.v*h.t))*h.L).sum(axis=0)
    gate_parity = float(np.max(abs(va-vb-image)))
    # Material step-halving: SAME spectra/q, restrict current fine solution to
    # every other time node; only the causal quadrature changes, not the model.
    coarse = object.__new__(Spectrum)
    coarse.__dict__ = dict(spectrum.__dict__)
    coarse.t = spectrum.t[::2]
    coarse.h = spectrum.h[:, ::2]
    coarse.H = cum(coarse.h, coarse.t)
    coarse.himage = spectrum.himage[::2]
    coarse.Himage = cum(coarse.himage, coarse.t)
    c = History(coarse, .25)
    fine_u, fine_window, _ = rate_history(h)
    coarse_u, coarse_window, _ = rate_history(c)
    step_peak = float(np.max(abs(fine_u[::2]-coarse_u)))
    # Inner-turn shift approximation sensitivity: evaluate at HALF the turn radius.
    half_u, _, _ = rate_history(h, turn_radius=C/(250*2))
    inner_shape = float(np.max(abs(half_u-fine_u)))
    # Same-rate point at sheet and a representative first-spacings collision;
    # this is the changed consumer, not merely an unused field diagnostic.
    def selected_rate(history, uu, depth):
        t = depth/history.v
        dk, _, _ = history.kinetic(depth)
        return float(DD.sigma(.25+dk/1000, U+np.interp(t, history.t, uu)/1000)/DD.sigma(.25,U))
    step_rate = max(abs(selected_rate(h, fine_u, z)-selected_rate(c, coarse_u, z))
                    for z in (0., .1, .5, 1., 2.))
    h2 = History(spectrum, .25, kcap=2.)
    peak_sensitivity = []
    for radius in RADII:
        a, frozen = h.potential(radius, True)
        b, _ = h2.potential(radius, True)
        peak_sensitivity.append(dict(R_A=radius,
            kcap1_peak_departure_eV=float(np.max(abs(a-frozen))),
            kcap2_peak_departure_eV=float(np.max(abs(b-frozen))),
            max_curve_change_eV=float(np.max(abs(a-b)))))
    proof = dict(source_full_G_head_identity_max=max(r['full_G_head_identity_residual'] for r in source['sources']),
        scalar_before_G_inversion_wrong_witness=source_order,
        immutable_profile_consumer_max_relative=baseline_error,
        support_t_negative='exact zero: all convolutions integrate only 0<=lag<=t; prehistory not included',
        kernel_t0_max=causal_zero, same_spectral_static_identity_max=static_identity,
        local_image_limit_absolute=image_limit, free_uniform_limit_absolute=free_limit,
        local_interface_continuity_absolute=continuity, identical_gate_entry_bulk_absolute_eV=gate_parity,
        projectile_kinetic_reference='causal self force minus formed moving bulk drag; CSDA drag retained once',
        formed_bulk_drag_at_lowest_energy_eV_A=h.formed_bulk_drag,
        step_halving_dt_fs=[float(coarse.t[1]), float(spectrum.t[1])],
        step_halving_pair_barrier_max_absolute_eV=step_peak,
        step_halving_selected_rate_max_absolute=step_rate,
        step_halving_edge_end_fs=[coarse_window['edge_end_fs'],fine_window['edge_end_fs']],
        half_Rturn_pair_barrier_max_absolute_eV=inner_shape,
        image_kcap_peak_sensitivity_at_lowest_energy=peak_sensitivity,
        image_kcap_1_to_2_max_yield_ratio_change=max(abs(a['raw_yield_ratio']-b['raw_yield_ratio'])
            for a,b in zip(primary,alternate_surface)),
        qcap_16_to_12_max_yield_ratio_change=max(abs(a['raw_yield_ratio']-b['raw_yield_ratio'])
            for a,b in zip(primary,alternate_tail)),
        input_hashes=source_hashes(), script_sha256=digest(__file__))
    if max(baseline_error, image_limit, free_limit, continuity, gate_parity, static_identity) > 1e-8:
        raise AssertionError(proof)
    (out/'check.json').write_text(json.dumps(proof, indent=2)+'\n')
    return proof


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.process_time(); wall = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=False)
    source = reduce_sources(args.matrix_dir, args.output)
    spectrum = Spectrum(source)
    windows, peaks = model_results(spectrum, args.output)
    primary = rate_consumer(spectrum)
    surface = rate_consumer(spectrum, kcap=2.)
    tail = rate_consumer(spectrum, qcap=12.)
    prehistory = rate_consumer(spectrum, formed_sheet=True)
    write_csv(args.output/'yields.csv', primary)
    write_csv(args.output/'sensitivity.csv', [dict(variant=name, **r) for name, rows in
        [('image_kcap_2',surface),('qcap_12',tail),('sheet_formed_bulk_wake',prehistory)] for r in rows])
    proof = narrow_check(source, spectrum, primary, surface, tail, args.output)
    kk_rows = []
    for i,q in enumerate(source['q']):
        loss = source['loss'][i]
        _, invmoment = loss_moments(source['energy'], loss)
        raw_kk = 1-2/np.pi*invmoment
        kk_rows.append(dict(q_Ainv=float(q), raw_static_inverse=float(source['static'][i]),
            sampled_piecewise_KK_inverse=float(raw_kk),
            realized_damped_hybrid_inverse=float(spectrum.frozen_raw[i]),
            moment_energy_eV=float(np.sqrt(source['moments'][i]))))
    write_csv(args.output/'static_moments.csv', kk_rows)
    result = dict(model='conditional pure-Pd isotropic full-inverse-head bulk plus local dielectric image',
        authoritative_base='59a56b0cd55e601f2e95988106550f3e98bcba5d',
        archive_sha256='1be3de744451f7938d44bd3dced26854849643e80e8c84f1f35e77991a384c10',
        sources=source['sources'], q_RPA_support_Ainv=[float(source['q'][0]),float(source['q'][-1])],
        RPA_frequency_support_eV=[0.,96.], actual_RPA_energy_grid_eV=source['energy'].tolist(),
        extended_energy_support_eV=[0.,2048.], q_support_Ainv=[0.,16.],
        moment_clock_eV=source['wp'], fsum_effective_density_A3=source['density'],
        fsum_kF_Ainv=source['kf'], Mermin_gamma_eV=GAMMA, realization_regulator_eV=REG,
        leading_large_q_self_tail_estimate_eV=2*C/np.pi*(8*np.pi*C*source['density']/AK)/(3*16**3),
        edge_windows=windows, check=proof,
        actual_numeric_cpu_s=time.process_time()-start, wall_s=time.monotonic()-wall,
        peak_rss_kB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        nonyield_observable='speed-dependent signed normal impulse and target-kick distribution in fixed_sites.csv; conditional recoil Doppler contribution, not an emission efficiency',
        omissions=['microscopic nonlocal half-space and nonlinear charge-state/prehistory dynamics',
                   'impact-accessible excited-store population and release/receiver coupling',
                   'self-consistent moving nuclei and a time-dependent tunnelling calculation'])
    (args.output/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(cpu_s=result['actual_numeric_cpu_s'], peak_rss_kB=result['peak_rss_kB'],
        windows=windows, check=proof, yield_ratios=[{k:r[k] for k in
        ('material','E_cm_keV','channel','raw_yield_ratio','actual_edge_reaction_share')} for r in primary]),indent=2))


if __name__ == '__main__':
    main()
