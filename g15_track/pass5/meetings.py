"""Whole-foil reduced optical meetings, with separate event and probe Palm measures.

An electrical parent is a mark of the accepted Cox rain, not a packet image,
Bloch component, or molecular constituent. This finite ray/RWA calculation
conditions on physical support before evaluating the common quartic. It is not
a measurement of atomic contact or a probability model for unretained modes.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.special import erf

QE = 1.602176634e-19
MD = 2.0141 * 1.66053906660e-27
KB = 8.617333262145e-5
L = np.sqrt(np.pi) * .012
H = .00025
ALPHA = 11.77e-6
ENERGIES = (.25, .5, 1., 1.5, 2., 2.5)
BEAMS = ('even_steady', 'hot_steady', 'hot_bursts')
ALLOCATIONS = ('onsite', 'coupling')
# Accepted Pass4 source/optical closures; radius is Gaussian amplitude sigma.
ARMS = {
    'fast_Cu': dict(radius=8e-9, ne=15e-15, opt=.0802e-12, dephase=.08e-12, share=.25),
    'nominal_retained': dict(radius=20e-9, ne=.1e-12, opt=.5e-12, dephase=.5e-12, share=.5),
    'long_retained': dict(radius=60e-9, ne=.5e-12, opt=2.068018471814257e-12, dephase=2e-12, share=.92),
}
MANIFESTS = {
    'core': '47f2c543f484419248afa5989d9cb4edbd4ccce9384c1c44d5f9b8e30a495883',
    'evidence': 'fad2f67b94f706dc14a77bfb0d36e561b556d3a8dd9a77c4e3c383d1c1fbcdf6',
}
GEOMETRY_NODES = 49
CONTACT_NODES = 81
MEETING_NODES = 33
PROBE_SAMPLES = 16
EVENT_SAMPLES = 1536
EVENT_RESPONSE_SAMPLES = 16
MAX_LIVE_SAMPLES = 96
SUPPORT = 3.


def _digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


class Accepted:
    """Bind the exact accepted source tables before using the first row."""

    def __init__(self, root):
        self.root = Path(root)
        self.expected, self.used = {}, {}
        for kind, sha in MANIFESTS.items():
            path = self.root / f'manifest_{kind}.json'
            if _digest(path) != sha:
                raise ValueError(f'Wrong accepted Pass4 manifest: {path}')
            metadata = json.loads(path.read_text())
            self.expected.update({r['path'].removeprefix('g15_track/pass4/'): r['sha256']
                                  for r in metadata['files']})

    def __call__(self, name):
        path = self.root / name
        if name not in self.used:
            digest = _digest(path)
            if digest != self.expected[name]:
                raise ValueError(f'Wrong accepted Pass4 input: {name}')
            self.used[name] = digest
        return path

    def json(self, name):
        return json.loads(self(name).read_text())


def _beam(xy, beam):
    x, y = np.moveaxis(xy, -1, 0)
    mask = (abs(x) <= .005) & (abs(y) <= .005)
    if beam == 'even_steady':
        return mask.astype(float)
    out = np.exp(-(x*x + y*y)/(2*.003**2))
    for a, b, g in ((.0025, -.0025, 5.), (-.0025, .0025, 3.)):
        out += g*np.exp(-((x-a)**2 + (y-b)**2)/(2*.0012**2))
    return mask*out


def _beam_norm(beam):
    q, w = np.polynomial.legendre.leggauss(64)
    xy = np.stack(np.meshgrid(.005*q, .005*q, indexing='ij'), axis=-1)
    val = _beam(xy, beam)
    vol = .005**2 * np.outer(w, w)
    return float(np.sum(vol*val)), float(np.sum(vol*val*val))


def _xy(rng, size, beam):
    result = np.empty((size, 2))
    left = np.arange(size)
    while len(left):
        trial = rng.uniform(-.005, .005, (len(left), 2))
        ok = rng.uniform(0, 1 if beam == 'even_steady' else 9, len(left)) < _beam(trial, beam)
        result[left[ok]] = trial[ok]
        left = left[~ok]
    return result


def _disk(rng, size, radius):
    az = rng.uniform(0., 2*np.pi, size)
    return (np.sqrt(rng.random(size))*radius)[:, None]*np.stack((np.cos(az), np.sin(az)), axis=-1)


def _focus_disk(rng, size, radius):
    """Oversample bright focus without losing the finite support/Jacobian."""
    inner = radius*(1.25/SUPPORT)
    point = _disk(rng, size, np.where(rng.random(size) < .65, inner, radius))
    inner_mask = np.einsum('ni,ni->n', point, point) <= inner**2
    pdf = .35/(np.pi*radius**2)+.65/(np.pi*inner**2)*inner_mask
    return point, pdf


class Clock:
    def __init__(self, path):
        with np.load(path) as data:
            self.edges = data['edges_s']
            self.on = data['burst_on']
            self.phase = float(data['phase'])
            self.norm = float(data['normalization'])
            self.seconds = float(data['seconds'])
        self.omega = 2*np.pi*1e5
        a, b = self.edges[:-1], self.edges[1:]
        ua, ub = self.omega*a+self.phase, self.omega*b+self.phase
        d = b-a
        c1 = (np.sin(ub)-np.sin(ua))/self.omega
        c2 = .5*d+(np.sin(2*ub)-np.sin(2*ua))/(4*self.omega)
        c3 = .75*c1+(np.sin(3*ub)-np.sin(3*ua))/(12*self.omega)
        baseline = 1+9*self.on
        self.integrals = {1: self.seconds,
                          2: float(np.sum(baseline**2*d + 1.4*baseline*c1 + .49*c2))/self.norm**2,
                          3: float(np.sum(baseline**3*d + 2.1*baseline**2*c1 +
                                          1.47*baseline*c2 + .343*c3))/self.norm**3}
        self.envelope = (baseline+.7)/self.norm
        self.proposal = {p: np.cumsum(d*self.envelope**p) for p in (2, 3)}

    def value(self, time, beam):
        inside = (time >= 0) & (time <= self.seconds)
        if beam != 'hot_bursts':
            return inside.astype(float)
        ix = np.clip(np.searchsorted(self.edges, time, side='right')-1, 0, len(self.on)-1)
        return inside*(1+.7*np.cos(self.omega*time+self.phase)+9*self.on[ix])/self.norm

    def draw(self, rng, size, beam, power):
        if beam != 'hot_bursts':
            return rng.uniform(0., self.seconds, size), np.full(size, 1/self.seconds)
        cumulative = self.proposal[power]
        result = np.empty(size)
        left = np.arange(size)
        while len(left):
            ix = np.searchsorted(cumulative, rng.uniform(0., cumulative[-1], len(left)))
            trial = rng.uniform(self.edges[ix], self.edges[ix+1])
            r = self.value(trial, beam)/self.envelope[ix]
            ok = rng.random(len(left)) < r**power
            result[left[ok]] = trial[ok]
            left = left[~ok]
        return result, self.value(result, beam)**power/self.integrals[power]


def _rate(ec):
    return float(np.interp(2*ec, np.arange(.5, 5.01, .5),
                           [38, 68, 94, 110, 118, 123, 130, 134, 136, 138]))*1e-6/QE


class Depth:
    def __init__(self, profile):
        self.layers = np.asarray(profile['layers'])
        self.bottom = self.layers[:, 1]
        self.density = self.layers[:, 3]

    def column_integral(self, z, width):
        """Exact near-face Gaussian integral against accepted piecewise D(z)."""
        z, width = np.broadcast_arrays(z, width)
        out = np.empty(z.size)
        # Block the exact layer sum; a homogeneous continuation is never used.
        lo = self.layers[:, 0]
        hi = self.bottom
        for start in range(0, z.size, 1024):
            sl = slice(start, start+1024)
            zz, ww = z.ravel()[sl, None], width.ravel()[sl, None]
            # Beyond 7 sigma erf has already rounded to the exact same ±1.
            first = np.searchsorted(hi, np.min(zz-7*ww))
            last = np.searchsorted(lo, np.max(zz+7*ww), side='right')
            active = slice(first, last)
            out[sl] = (.5*np.sqrt(np.pi)*ww[:, 0] *
                       ((erf((hi[active]-zz)/ww)-erf((lo[active]-zz)/ww)) @ self.density[active]))
        return out.reshape(z.shape)

    def probe(self, rng, size, marks):
        sp = (rng.random(size) < .30/1.15).astype(np.int8)
        xyz = np.empty(size)
        lag, energy, column = np.empty((3, size))
        for j, track in enumerate(marks):
            ix = np.flatnonzero(sp == j)
            z = track.z_bound
            # CDF on the intersections of the exact loading layers and track segments.
            edge = np.unique(np.r_[z, self.bottom[self.bottom < z[-1]]])
            nd = self.density[np.minimum(np.searchsorted(self.bottom, .5*(edge[:-1]+edge[1:])), len(self.bottom)-1)]
            true_column = np.diff(edge)*nd
            cdf = np.cumsum(true_column)
            u = rng.uniform(0., cdf[-1], len(ix))
            k = np.searchsorted(cdf, u, side='right')
            zz = edge[k]+(u-np.r_[0., cdf][k])/nd[k]
            xyz[ix] = zz
            energy[ix] = np.interp(zz, z, track.energy)/1000
            lag[ix] = np.interp(zz, z, track.time) + rng.integers(0, j+1, len(ix))*track.delay
            column[ix] = cdf[-1]
        return xyz, sp+1, energy, column, lag


class Material:
    def __init__(self, accepted, lattice):
        self.mass = lattice['mass_eVc2']*QE/299792458.**2
        self.M = self.mass*1e-20/QE  # eV s² Å⁻²
        self.A = lattice['K4_100_eV_A4']
        self.C = 2*self.A + 3*lattice['K4_111_eV_A4']
        self.g = lattice['collective_strain_slope_eV_A2']
        with accepted('frozen/g16_response.csv').open() as stream:
            tab = sorted((r for r in csv.DictReader(stream) if int(r['kmesh_axis_count']) == 16),
                         key=lambda r: float(r['electronic_excitation_eV_per_D_site']))
        self.u = np.array([float(r['electronic_excitation_eV_per_D_site']) for r in tab])
        self.dK = np.array([float(r['K_change_eV_A2']) for r in tab])
        with np.load(accepted('frozen/native_optical_basis.npz')) as data:
            self.k, self.D, self.vg = data['k'], data['onsite_D'], data['onsite_vg']
            oldD = data['old_D']
        self.shell = np.asarray(lattice['shells'][0]['periodic_alias_offsets_A'])*1e-10
        self.F = np.cos(self.k @ self.shell.T).mean(axis=1)
        self.gradF = -np.sin(self.k @ self.shell.T) @ self.shell / len(self.shell)
        bonds, spring = [], []
        for shell in lattice['shells']:
            aliases = shell['periodic_alias_offsets_A']
            for row in aliases:
                half_cell = np.rint(2*np.asarray(row)/lattice['a_A']).astype(int)
                images = tuple(itertools.product(*[(-2, 2) if abs(v) == 2 else (v,) for v in half_cell]))
                for image in images:
                    bonds.append(np.asarray(image)*lattice['a_A']*.5e-10)
                    spring.append(shell['J_sum_eV_A2']/(len(aliases)*len(images)))
        bonds, spring = np.asarray(bonds), np.asarray(spring)
        qphase = self.k @ bonds.T
        HessK = np.einsum('kb,b,bi,bj->kij', np.cos(qphase), spring, bonds, bonds)
        gamma = int(np.argmin(np.linalg.norm(self.k, axis=1)))
        qc = min(np.linalg.norm(self.k, axis=1)[np.arange(len(self.k)) != gamma])
        eta = oldD[gamma]/self.D[gamma]-1
        # W and its first two derivatives vanish at the boundary of the Gamma bin.
        HessK[gamma] = HessK[gamma]/(1+eta) + np.eye(3)*oldD[gamma]*6*eta/(qc**2*(1+eta)**2)
        w0 = np.sqrt(self.D/self.M)
        self.hessian = HessK/(2*self.M*w0[:, None, None])-np.einsum('ki,kj->kij', self.vg, self.vg)/w0[:, None, None]
        self.spread = np.linalg.norm(self.hessian, axis=(1, 2))/np.sqrt(3.)

    def electronic(self, u):
        if np.any(u < -1e-12) or np.any(u > self.u[-1]*(1+1e-12)):
            raise ValueError('Shared Ue leaves accepted m16 constitutive support')
        return np.interp(u, self.u, self.dK)

    def electronic_slope(self, u):
        ix = np.clip(np.searchsorted(self.u, u, side='right')-1, 0, len(self.u)-2)
        return (self.dK[ix+1]-self.dK[ix])/(self.u[ix+1]-self.u[ix])

    def electronic_project(self, u, mode, allocation):
        if allocation not in ALLOCATIONS:raise ValueError('Native projection requires its declared material rule')
        delta=self.electronic(u);slope=self.electronic_slope(u)
        if allocation=='onsite':return delta,slope,0.
        shape=(len(mode),)+(1,)*(np.ndim(u)-1)
        F=self.F[mode].reshape(shape);gradF=self.gradF[mode].reshape(shape+(3,))
        return delta*F,slope*F,delta[...,None]*gradF

    def source_project(self, delta, mode, allocation):
        if allocation not in ALLOCATIONS:raise ValueError('Unknown incoming source projection')
        if allocation=='onsite':return delta,0.
        shape=(len(mode),)+(1,)*(delta.ndim-1)
        return delta*self.F[mode].reshape(shape),delta[...,None]*self.gradF[mode].reshape(shape+(3,))


    def variance(self, curvature, temperature):
        if np.any(curvature <= 0) or np.any(temperature <= 0):
            raise ValueError('Sample leaves stable Pass4 thermal pole')
        return 2*KB*temperature/(curvature+np.sqrt(curvature**2+4*self.C*KB*temperature))

    def pole(self, curvature):
        if np.any(curvature <= 0):
            raise ValueError('Nonpositive propagated D optical curvature')
        return np.sqrt(curvature/self.M)


class Mark:
    def __init__(self, track, cache, sp, arm, debit, material):
        rows = track['timeline']
        self.delay = .74e-10/np.sqrt(2*track['initial_eV']*QE/MD) if sp else 0.
        self.z_bound = np.r_[0., [r['depth_A']*1e-10 for r in rows]]
        self.time = np.r_[0., [r['time_fs']*1e-15 for r in rows]]
        self.energy = np.r_[track['initial_eV'], [r['energy_eV'] for r in rows]]
        deposited = np.diff(np.r_[0., [r['D_recoil_eV'] for r in rows]])
        self.z_prob = deposited/deposited.sum()
        self.direct = .5*track['components_eV'][2]
        raw = np.asarray(track['components_eV'])*(sp+1)
        raw[:2] -= (debit['electron'], debit['host'])
        share = arm['share']
        self.budget = float(share*raw[0]+.05*raw[1]+(.5+.25*share)*raw[2]+
                            (.25+.3*share)*track['terminal_unresolved_eV']*(sp+1))
        with np.load(cache) as data:
            if not np.array_equal(data['k_m_inverse'], material.k):
                raise ValueError('Pass4 native modal chart differs from accepted optical basis')
            self.t = np.r_[0., data['time_s']]
            self.site = float(data['sites'])
            stiffness = data['base_K']
            selected = data['selected_mode_indices']
            nearest = np.argmin(abs(stiffness[:, None]-stiffness[selected][None, :]), axis=1)
            quiet = material.variance(stiffness[selected], 300.)
            self.q = np.vstack((np.zeros(len(stiffness)), data['variance_q_selected_A2'][:, nearest]-quiet[nearest]))
            self.phase_excess = np.vstack((np.zeros(len(stiffness)), (data['phase_selected_rad']-
                                            data['reference_phase_selected_rad'])[:, nearest]))
            self.dK = np.r_[0., data['projected_deltaK_eVA2']]
            self.theta = np.vstack((np.zeros(len(stiffness)), data['velocity_theta_excess_eV']))
            self.hot = np.r_[0., data['hot_e_eV']]
            self.thermal = np.r_[0., data['thermal_e_excess_eV']]
            self.optical = np.r_[0., data['optical_excess_eV']]
        self.optical_peak = float(self.optical.max())


def _draw_marks(rng, size, marks, material):
    species = (rng.random(size) < .15).astype(np.int8)
    mode = rng.integers(len(material.k), size=size)
    z, birth, delay, end, site, direct, budget = np.empty((7, size))
    for j, mark in enumerate(marks):
        ix = np.flatnonzero(species == j)
        segment = rng.choice(len(mark.z_prob), len(ix), p=mark.z_prob)
        u = rng.random(len(ix))
        z[ix] = mark.z_bound[segment]+u*np.diff(mark.z_bound)[segment]
        birth[ix] = mark.time[segment]+u*np.diff(mark.time)[segment]
        delay[ix], end[ix], site[ix], direct[ix], budget[ix] = mark.delay, mark.t[-1], mark.site, mark.direct, mark.budget
    return dict(sp=species, mode=mode, z=z, formation=birth, delay=delay, end=end,
                duration=end-birth, sites=site, direct=direct, budget=budget, quiet_normal=rng.standard_normal(size))


def _native(a, age, marks, key):
    age = np.asarray(age)
    flat = age.reshape((len(a['sp']), -1))
    result = np.empty_like(flat)
    for j, m in enumerate(marks):
        ids = np.flatnonzero(a['sp'] == j)
        values = getattr(m, key)
        if values.ndim == 1:
            result[ids] = np.interp(flat[ids], m.t, values)
        else:
            for mode in np.unique(a['mode'][ids]):
                ix = ids[a['mode'][ids] == mode]
                result[ix] = np.interp(flat[ix], m.t, values[:, mode])
    return result.reshape(age.shape)


def _response_fraction(a, flight, marks, arm):
    native_peak = np.array([m.optical_peak for m in marks])[a['sp'], None]
    optical = _native(a, a['formation'][:, None]+flight, marks, 'optical')/native_peak
    return np.maximum(optical, np.exp(-flight/arm['opt']))


def _fold(z):
    half = np.floor(z/H).astype(np.int64)
    sign = 1-2*(half % 2)
    return np.where(sign > 0, z-half*H, (half+1)*H-z), sign


def _cumulative(values, step):
    return np.concatenate((np.zeros_like(values[:, :1]),
                           np.cumsum((values[:, 1:]+values[:, :-1])*.5*step, axis=1)), axis=1)


def _strain(strain_sampler, ec, beam, xyz, t, allocation, control):
    signed = np.asarray(strain_sampler(ec, beam, xyz, mount='nominal', time_s=t,
                                      allocation=allocation, control=control))
    if signed.shape != xyz.shape[:-1]+(3, 3):
        raise ValueError('Unexpected shared strain chart shape')
    return np.trace(signed, axis1=-2, axis2=-1)/3


def _rays(a, marks, material, arm, allocation, ec, beam, xy, launch, medium_sampler, strain_sampler,
          controls=False, mean_strain=None):
    """Retarded eikonal at the source chart, with leading refractive gradients.

    The slow common thermal/electronic chart is transported to leading spatial
    order; signed strain is sampled at the source and held during each ps ray.
    D-profile and paths are evaluated in the foil, not a launch-centred box.
    """
    chart = np.c_[xy, a['z']]
    absolute = np.clip(launch+a['formation'], 0., 2.)
    T, U, gradT, gradU = medium_sampler(ec, beam, chart, mount='nominal', time_s=absolute)
    strain = (np.full_like(T, mean_strain) if mean_strain is not None else
              _strain(strain_sampler, ec, beam, chart, absolute, allocation, 'joined'))
    mode = a['mode']
    base = material.D[mode]
    electronic,_,electronic_q=material.electronic_project(U,mode,allocation)
    q0 = material.variance(base, 300.)
    w_frozen = material.pole(base + material.C*q0)
    slow = material.g*(strain+ALPHA*(T-300.))
    duration = a['duration']
    t = duration[:, None]*np.linspace(0., 1., GEOMETRY_NODES)
    dt = (duration/(GEOMETRY_NODES-1))[:, None]
    ages = t+a['formation'][:, None]
    near_dK = _native(a, ages, marks, 'dK')
    near_e,near_q_gradient=material.source_project(near_dK,mode,allocation)
    near_q = _native(a, ages, marks, 'q')
    warm = material.variance(base+slow+electronic, T)
    omega = material.pole(base[:, None]+slow[:, None]+electronic[:, None]+near_e+material.C*(warm[:, None]+near_q))
    vg0 = material.vg[mode]*material.pole(base)[:, None]
    warm_screen = (base+slow+electronic+material.C*warm)/(base+slow+electronic+2*material.C*warm)
    vg = ((vg0+electronic_q/(2*material.M))[:, None, :]*warm_screen[:, None, None]+
          near_q_gradient/(2*material.M))/omega[..., None]
    drift = _cumulative(vg, dt[..., None])
    dz, parity = _fold(a['z'][:, None]+drift[..., 2])
    delta = drift.copy()
    delta[..., 2] = dz-a['z'][:, None]
    T_ray = T[:, None]+np.einsum('nki,ni->nk', delta, gradT)
    dU = np.einsum('nki,ni->nk', delta, gradU)
    # U is a positive mean; first-order spatial chart may leave support only at zero.
    U_ray = np.maximum(U[:, None]+dU, 0.)
    common_e,common_slope,common_q=material.electronic_project(U_ray,mode,allocation)
    K_el = common_e+near_e
    K_slow = material.g*(strain[:, None]+ALPHA*(T_ray-300.))
    bare = base[:, None]+K_slow+K_el-near_e
    q_quiet = material.variance(bare, T_ray)
    q_ray = q_quiet+near_q
    K_pop = material.C*q_ray
    omega = material.pole(base[:, None]+K_slow+K_el+K_pop)
    # Leading Hamiltonian ray kick from the resolved thermal/electronic gradients.
    denominator = bare+2*material.C*q_quiet
    screen = (bare+material.C*q_quiet)/denominator
    dK = ((material.g*ALPHA*screen+material.C*KB/denominator)[..., None]*gradT[:, None, :]+
          (common_slope*screen)[..., None]*gradU[:, None, :])
    dK[..., 2] *= parity
    kick = _cumulative(-dK/(2*material.M*omega[..., None]), dt[..., None])
    dispersion_numerator = ((vg0[:, None, :]+common_q/(2*material.M))*screen[..., None]+
                            near_q_gradient/(2*material.M))
    vg = dispersion_numerator/omega[..., None]+np.einsum('nij,ntj->nti', material.hessian[mode], kick)
    drift = _cumulative(vg, dt[..., None])
    k = material.k[mode, None, :]+kick
    near_only = material.pole(base[:, None]+near_e+material.C*(q0[:, None]+near_q))
    # Archive resolves the single-parent fs quench; integrate only the common
    # warm innovation in addition, rather than adding the same near phase twice.
    phase = _cumulative(omega-near_only-np.einsum('nki,nki->nk', vg, k), dt)+w_frozen[:, None]*t
    phase_frozen = w_frozen[:, None]*t-np.einsum('nki,ni->nk', drift, material.k[mode])
    if controls:
        reference = _strain(strain_sampler, ec, beam, chart, absolute, allocation, 'strain_off')
        if np.any(reference != 0.):
            raise ValueError('Matched named strain-off chart is not zero')
        slow_off = material.g*(reference[:, None]+ALPHA*(T_ray-300.))
        bare_off = base[:, None]+slow_off+K_el-near_e
        q_off = material.variance(bare_off, T_ray)
        omega_off = material.pole(bare_off+near_e+material.C*(q_off+near_q))
        dslow = material.g*(strain-reference)[:, None]
        dq = -dslow*q_off/(bare+material.C*(q_quiet+q_off))
        deltaw = (dslow+material.C*dq)/(material.M*(omega+omega_off))
        dT = T_ray-300.
        dK_T = material.g*ALPHA*dT
        bare_T0 = bare-dK_T
        q_T0 = material.variance(bare_T0, 300.)
        omega_T0 = material.pole(bare_T0+near_e+material.C*(q_T0+near_q))
        dq_T = (KB*dT-dK_T*q_T0)/(bare+material.C*(q_quiet+q_T0))
        dw_T = (dK_T+material.C*dq_T)/(material.M*(omega+omega_T0))
        a.update(strain_phase=_cumulative(deltaw, dt),
                 max_fractional_strain_pole_change=np.max(abs(deltaw/omega_off), axis=1),
                 max_fractional_warm_pole_change=np.max(abs(dw_T/omega_T0), axis=1))
    radius = arm['radius']
    # Isotropized RMS width from the inherited branch Hessian and Gaussian Fourier width.
    width = radius*np.sqrt(1+(material.spread[mode, None]*t/radius**2)**2)
    a.update(duration=duration, drift=drift, k=k, phase=phase, phase_frozen=phase_frozen,
             omega=omega, width=width, T_birth=T, strain_birth=strain,
             K_background=K_slow+(K_el-near_e)+material.C*q_quiet,
             K_parent_e=near_e, K_parent_q=material.C*near_q)
    return a


def _at(a, key, flight):
    arr = a[key]
    shape = (-1,) + (1,)*(flight.ndim-1)
    position = np.clip(flight/a['duration'].reshape(shape), 0., 1.)*(GEOMETRY_NODES-1)
    lo = np.minimum(position.astype(np.int64), GEOMETRY_NODES-2)
    ix = np.arange(len(lo)).reshape(shape)
    weight = position-lo
    if arr.ndim == 3:
        weight = weight[..., None]
    return (1-weight)*arr[ix, lo]+weight*arr[ix, lo+1]


def _centres(a, xy, flight):
    pos = _at(a, 'drift', flight)
    shape = (-1,) + (1,)*(flight.ndim-1)
    pos[..., :2] += xy.reshape(shape+(2,))
    pos[..., 2], reflection = _fold(pos[..., 2]+a['z'].reshape(shape))
    return pos, reflection


def _age_proposal(rng, a, arm):
    u = rng.random(len(a['sp']))
    norm = -np.expm1(-a['duration']/arm['opt'])
    flight = np.where(rng.random(len(u)) < .6, -arm['opt']*np.log1p(-u*norm), u*a['duration'])
    q = .6*np.exp(-flight/arm['opt'])/(arm['opt']*norm)+.4/a['duration']
    return flight, q


def _overlap(depth, center1, center2, width1, width2, active):
    """D-normalized spatial matrix element ∫ n_D |φ1|² |φ2|²."""
    v1, v2 = width1**2, width2**2
    v = v1*v2/(v1+v2)
    d = center1-center2
    distance2 = np.einsum('...i,...i->...', d, d)
    z1, z2 = center1[..., 2], center2[..., 2]
    n1 = np.pi*v1*depth.column_integral(z1, width1)
    n2 = np.pi*v2*depth.column_integral(z2, width2)
    # Action-normalized spatial modes use the same D measure as the pair potential.
    z = (v2*z1+v1*z2)/(v1+v2)
    # Smooth compact support gate is part of the finite RWA projected law;
    # it leaves the appreciable Gaussian overlap and removes unpriced tails.
    overlap = distance2/(SUPPORT*(width1+width2))**2
    gate = np.maximum(1-overlap**4, 0.)**2*active
    return np.pi*v*depth.column_integral(z, np.sqrt(v))*np.exp(-distance2/(v1+v2))*gate/(n1*n2)


def _launch_action(a):
    """Accepted half-direct impulse; molecular echo and quiet p·v never re-marked."""
    omega = a['omega'][:, 0]
    direct, budget = a['direct'], a['budget']
    quiet_kinetic = .5*KB*a['T_birth']*a['quiet_normal']**2
    quiet_cross = 2*np.sign(a['quiet_normal'])*np.sqrt(direct*quiet_kinetic)
    # Stable positive root; enforce the same finite source debit if necessary.
    s1 = np.minimum(1., 2*budget/(quiet_cross+np.sqrt(quiet_cross**2+4*direct*budget)))
    work1 = direct*s1**2+quiet_cross*s1
    echo = omega*a['delay']
    molecular_cross = (quiet_cross+2*direct*s1)*np.cos(echo)
    remaining = np.maximum(budget-work1, np.finfo(float).tiny)
    denominator = molecular_cross+np.sqrt(molecular_cross**2+4*direct*remaining)
    s2 = np.zeros_like(s1)
    is_molecular = a['sp'] == 1
    s2[is_molecular] = np.minimum(1., 2*remaining[is_molecular]/denominator[is_molecular])
    work2 = direct*s2*s2+molecular_cross*s2
    a.update(action1=s1*np.sqrt(direct/omega), action2=s2*np.sqrt(direct/omega),
             launch_work=work1+work2, launch_quiet_cross=quiet_cross*(s1+s2*np.cos(echo)))


def _action(a, flight, arm):
    shape = (-1,) + (1,)*(flight.ndim-1)
    e = a['delay'].reshape(shape)
    w = a['omega'][:, 0].reshape(shape)
    s1, s2 = a['action1'].reshape(shape), a['action2'].reshape(shape)
    relative = np.exp(1j*w*e)*(flight >= e)*s2*np.exp(e/(2*arm['opt']))
    amplitude = (s1+relative)*np.exp(-np.maximum(flight, 0.)/(2*arm['opt']))
    return np.where(flight >= 0., amplitude, 0.)


def _meeting_mesh(duration, *echoes):
    """Keep both sides of each molecular injection in the quadrature."""
    grid = [duration[:, None]*np.linspace(0., 1., MEETING_NODES)**2]
    for echo in echoes:
        near_echo = np.clip(echo, 0., duration)
        eps = 8*np.spacing(duration)
        grid.extend((np.maximum(near_echo-eps, 0.)[:, None], np.minimum(near_echo+eps, duration)[:, None]))
    return np.sort(np.concatenate(grid, axis=1), axis=1)


def _track_quartic(a, b, f1, f2, xy1, xy2, material, depth, arm, marks):
    """Common spatial quartic in action variables: coherent + Hartree population.

    Mode actions change by the accepted source/damping only. Native Gaussian
    population drives a separate coherent phase through the same quartic and
    physically pays its signed reciprocal work; it is not assigned a phase.
    """
    # Integrate each self phase from its own birth, and cross phase from the
    # later birth. Conditioning on the output probe does not change either clock.
    own_phase = []
    for par, xy, end in ((a, xy1, f1), (b, xy2, f2)):
        ti = _meeting_mesh(end, par['delay'])
        center, _ = _centres(par, xy, ti)
        width = _at(par, 'width', ti)
        wc = _at(par, 'omega', ti)
        v = width**2
        normal = np.pi*v*depth.column_integral(center[..., 2], width)
        s = (3*material.A/(8*material.M**2))*.5*np.pi*v*depth.column_integral(center[..., 2], width/np.sqrt(2.))/(normal*wc)**2
        action = abs(_action(par, ti, arm))**2
        own_phase.append(np.trapezoid(2*s*action, ti, axis=1))
    length = np.minimum(f1, f2)
    start1, start2 = f1-length, f2-length
    meeting_time = _meeting_mesh(length, a['delay']-start1, b['delay']-start2)
    t1, t2 = start1[:, None]+meeting_time, start2[:, None]+meeting_time
    c1, _ = _centres(a, xy1, t1)
    c2, _ = _centres(b, xy2, t2)
    w1, w2 = _at(a, 'width', t1), _at(b, 'width', t2)
    active = (t1 >= 0.) & (t2 >= 0.)
    omega1, omega2 = _at(a, 'omega', t1), _at(b, 'omega', t2)
    overlap = _overlap(depth, c1, c2, w1, w2, active)
    cross = (3*material.A/(2*material.M**2))*overlap/(omega1*omega2)
    i1, i2 = abs(_action(a, t1, arm))**2, abs(_action(b, t2, arm))**2
    dt = np.diff(meeting_time, axis=1)
    # The same accepted Gaussian q variance contains the direct population
    # allocation. Remove its retained coherent first-moment share before the
    # Hartree contraction, so a coherent packet is not a second thermal source.
    incoherent = []
    for par, flight, act, pole in ((a, t1, i1, omega1), (b, t2, i2, omega2)):
        native_q = _native(par, par['formation'][:, None]+flight, marks, 'q')
        replacement_q = act*np.exp(-2*flight/arm['dephase'])/(3*material.M*pole*par['sites'][:, None])
        incoherent.append(native_q-replacement_q)
    wpop1 = material.C*overlap*b['sites'][:, None]/(2*material.M*omega1)
    wpop2 = material.C*overlap*a['sites'][:, None]/(2*material.M*omega2)
    cross_phase = np.stack((np.trapezoid(cross*i2+wpop1*incoherent[1], meeting_time, axis=1),
                            np.trapezoid(cross*i1+wpop2*incoherent[0], meeting_time, axis=1)), axis=1)
    product = i1*i2
    cross_energy = cross*product
    # Exact discrete work split: Δ(cP)=P̄Δc+c̄ΔP; P decays twice as fast as I
    # between launches. Charge the finite overlap term at the later birth to
    # that parent's launch port, not to a new copy of the old optical budget.
    coefficient_work = np.sum(np.diff(cross, axis=1)*.5*(product[:, 1:]+product[:, :-1]), axis=1)
    middle = .5*(cross[:, 1:]+cross[:, :-1])
    decay = -np.expm1(-2*dt/arm['opt'])*product[:, :-1]
    damping = np.sum(middle*decay, axis=1)
    launch_work = cross_energy[:, 0]+np.sum(middle*(np.diff(product, axis=1)+decay), axis=1)
    optical_delta = np.zeros((len(length), 2))
    cpop_work = np.zeros_like(length)
    cpop_energy = np.zeros_like(length)
    survival = np.exp(-np.maximum(length[:, None]-.5*(meeting_time[:, 1:]+meeting_time[:, :-1]), 0.)/arm['opt'])
    for donor, coeff, act, cov in ((1, wpop1, i1, incoherent[1]), (0, wpop2, i2, incoherent[0])):
        prod = act*cov
        term = coeff*prod
        cbar = .5*(coeff[:, 1:]+coeff[:, :-1])
        abar = .5*(act[:, 1:]+act[:, :-1])
        qbar = .5*(cov[:, 1:]+cov[:, :-1])
        pop_transfer = cbar*abar*np.diff(cov, axis=1)
        decrement = -np.expm1(-dt/arm['opt'])*act[:, :-1]
        src = cbar*qbar*(np.diff(act, axis=1)+decrement)
        loss = cbar*qbar*decrement
        shape = np.diff(coeff, axis=1)*.5*(prod[:, 1:]+prod[:, :-1])
        optical_delta[:, donor] = -(pop_transfer*survival).sum(axis=1)
        cpop_work += pop_transfer.sum(axis=1)
        coefficient_work += shape.sum(axis=1)
        launch_work += term[:, 0]+src.sum(axis=1)
        damping += loss.sum(axis=1)
        cpop_energy += term[:, -1]
    final_energy = cross_energy[:, -1]+cpop_energy
    defect = final_energy-coefficient_work-launch_work+damping-cpop_work
    return (np.stack(own_phase, axis=1), cross_phase, final_energy, coefficient_work,
            launch_work, damping, defect, cpop_work, optical_delta)


def _velocity(a, flight, xyz, center, reflect, amp, own_phase, depth, mat, arm, marks):
    width = _at(a, 'width', flight)
    delta = xyz-center
    rr = np.einsum('ni,ni->n', delta, delta)/width**2
    focus = np.exp(-.5*rr)*(rr <= SUPPORT**2)
    sites = np.pi*width**2*depth.column_integral(center[:, 2], width)
    omega = _at(a, 'omega', flight)
    k = _at(a, 'k', flight)
    k[:, 2] *= reflect
    spatial = np.einsum('ni,ni->n', k, delta)
    source_phase = _native(a, a['formation']+flight, marks, 'phase_excess')-_native(
        a, a['formation'], marks, 'phase_excess')
    phase = spatial-_at(a, 'phase', flight)-source_phase-own_phase
    frozen_phase = spatial-_at(a, 'phase_frozen', flight)-own_phase
    eta = np.exp(-flight/arm['dephase'])
    magnitude = focus*1e-10*eta*np.sqrt(2*omega/(mat.M*sites))
    velocity = 1j*amp*magnitude*np.exp(1j*phase)
    frozen = 1j*amp*magnitude*np.exp(1j*frozen_phase)
    strain_off = velocity*np.exp(1j*_at(a, 'strain_phase', flight))
    # Kinetic covariance lost to phase memory stays in the same optical bucket.
    # The parent covariance already includes direct D. Replace, do not append,
    # the retained coherent kinetic second moment (phase-averaged, three axes).
    replacement = omega*abs(amp)**2/(3*sites)*(focus*eta)**2
    return velocity, frozen, strain_off, focus**2*a['sites']/sites, replacement, k


def _parameter_work(a, flight, marks, material, arm):
    """Finite eikonal curvature work and opposite coherent/carrier port amounts."""
    time = _meeting_mesh(flight, a['delay'])
    I = abs(_action(a, time, arm))**2
    wi = _at(a, 'omega', time)
    factor = (I[:, 1:]+I[:, :-1])/(2*material.M*(wi[:, 1:]+wi[:, :-1]))
    electronic_steps = np.diff(_at(a, 'K_parent_e', time), axis=1)*factor
    pop_steps = np.diff(_at(a, 'K_parent_q', time), axis=1)*factor
    # Moving through a held background is spatial/refraction work, not an
    # independent thermal pump or a debit from this parent's hot reservoir.
    shape_steps = np.diff(_at(a, 'K_background', time), axis=1)*factor
    hot = _native(a, a['formation'][:, None]+.5*(time[:, 1:]+time[:, :-1]), marks, 'hot')
    thermal = _native(a, a['formation'][:, None]+.5*(time[:, 1:]+time[:, :-1]), marks, 'thermal')
    optical = _native(a, a['formation'][:, None]+.5*(time[:, 1:]+time[:, :-1]), marks, 'optical')
    fraction = np.divide(hot, hot+thermal, out=np.zeros_like(hot), where=hot+thermal>0)
    hot_steps = electronic_steps*fraction
    hot_work = np.sum(hot_steps, axis=1)
    thermal_work = electronic_steps.sum(axis=1)-hot_work
    optical_work = np.sum(pop_steps, axis=1)
    shape_work = np.sum(shape_steps, axis=1)
    average_pole = .5*(wi[:, 1:]+wi[:, :-1])
    decrement = -np.expm1(-np.diff(time, axis=1)/arm['opt'])*I[:, :-1]
    loss = np.sum(decrement*average_pole, axis=1)
    source = np.sum((np.diff(I, axis=1)+decrement)*average_pole, axis=1)
    change = I[:, -1]*wi[:, -1]-I[:, 0]*wi[:, 0]
    defect = change-hot_work-thermal_work-optical_work-shape_work-source+loss
    fractions = (np.maximum(electronic_steps, 0.)/np.maximum(hot+thermal, 1e-300),
                 np.maximum(pop_steps, 0.)/np.maximum(optical, 1e-300))
    # Apply the opposite same-parent donor perturbations through the exact
    # inherited hot-e → optical and optical survival Green functions. Thermal
    # e/Pd are one-way exits for this observable, not extra optical sources.
    lag = np.maximum(flight[:, None]-.5*(time[:, 1:]+time[:, :-1]), 0.)
    survival = np.exp(-lag/arm['opt'])
    hot_to_optical = arm['share']*arm['opt']/(arm['opt']-arm['ne'])*(
        survival-np.exp(-lag/arm['ne']))
    optical_delta = -np.sum(pop_steps*survival+hot_steps*hot_to_optical, axis=1)
    work = np.stack((hot_work, thermal_work, optical_work, shape_work, source, .2*loss, .8*loss, defect,
                     fractions[0].max(axis=1), fractions[1].max(axis=1)), axis=1)
    return work, optical_delta


def _pair_rays(rng, size, marks, mat, arm, allocation, ec, beam, xy, launches, medium, strain, mean_strain):
    a = _draw_marks(rng, size, marks, mat)
    b = _draw_marks(rng, size, marks, mat)
    controls = mean_strain is None
    a = _rays(a, marks, mat, arm, allocation, ec, beam, xy, launches[0], medium, strain,
              controls=controls, mean_strain=mean_strain)
    b = _rays(b, marks, mat, arm, allocation, ec, beam, xy, launches[1], medium, strain,
              controls=controls, mean_strain=mean_strain)
    return a, b


def _support_residence(gap, duration):
    """Piecewise-linear signed-gap crossings; count first pair, not each echo."""
    lo, hi = gap[:, :-1], gap[:, 1:]
    cross_fraction = np.divide(lo, lo-hi, out=np.zeros_like(lo), where=lo != hi)
    fractions = np.where((lo < 0) & (hi < 0), 1.,
                         np.where(lo < 0, cross_fraction, np.where(hi < 0, 1-cross_fraction, 0.)))
    dt = duration/(gap.shape[1]-1)
    duration_inside = dt*np.sum(fractions, axis=1)
    entered = duration_inside > 0.
    i = np.argmax(fractions > 0., axis=1)
    ix = np.arange(len(gap))
    first = dt*(i+np.where(lo[ix, i] < 0., 0., cross_fraction[ix, i]))
    return duration_inside, entered, np.where(entered, first, np.nan)


def _events(rng, size, marks, mat, depth, arm, allocation, ec, beam, clock, norm, medium, strain, mean_strain=None):
    xy = _xy(rng, size, beam)
    time, q_time = clock.draw(rng, size, beam, 2)
    tau = max(m.t[-1] for m in marks)
    dtime = rng.uniform(0., tau, size)   # chronological labels: no extra 1/2
    a, b = _pair_rays(rng, size, marks, mat, arm, allocation, ec, beam, xy, (time, time+dtime),
                      medium, strain, mean_strain)
    # Superset of lateral offsets capable of contact, including both centre shifts.
    radius = (np.linalg.norm(a['drift'][..., :2], axis=-1).max(axis=1)+
              np.linalg.norm(b['drift'][..., :2], axis=-1).max(axis=1)+
              SUPPORT*(a['width'][:, -1]+b['width'][:, -1]))
    xy2 = xy+_disk(rng, size, radius)
    start = np.maximum(a['formation'], dtime+b['formation'])
    stop = np.minimum(np.minimum(a['end'], dtime+b['end']), clock.seconds-time)
    span = np.maximum(stop-start, 0.)
    local = start[:, None]+span[:, None]*np.linspace(0., 1., CONTACT_NODES)
    f1, f2 = local-a['formation'][:, None], local-dtime[:, None]-b['formation'][:, None]
    p1, _ = _centres(a, xy, f1)
    p2, _ = _centres(b, xy2, f2)
    distance = np.linalg.norm(p1-p2, axis=-1)
    summed_width = _at(a, 'width', f1)+_at(b, 'width', f2)
    dwell, near, first = _support_residence(distance-SUPPORT*summed_width, span)
    _, core, _ = _support_residence(distance-summed_width, span)
    response = np.minimum(_response_fraction(a, f1, marks, arm), _response_fraction(b, f2, marks, arm))
    appreciable_gap = np.maximum((distance-summed_width)/summed_width, np.exp(-3.)-response)
    bright_dwell, bright, _ = _support_residence(appreciable_gap, span)
    r = _rate(ec)
    parent2 = _beam(xy2, beam)/norm*clock.value(time+dtime, beam)
    weight = r*r*clock.value(time, beam)*parent2*np.pi*radius**2*tau/(q_time*size)
    contact1 = start+np.nan_to_num(first)-a['formation']
    contact2 = start+np.nan_to_num(first)-dtime-b['formation']
    first1, _ = _centres(a, xy, contact1)
    first2, _ = _centres(b, xy2, contact2)
    s1, s2 = _at(a, 'width', contact1), _at(b, 'width', contact2)
    first_xyz = (s2[:, None]*first1+s1[:, None]*first2)/(s1+s2)[:, None]
    first_xyz[~near] = np.nan
    valid_count = float(weight @ near)
    bright_count = float(weight @ bright)
    result = dict(event_weight=weight, event_near=near, event_core=core, event_appreciable_core=bright,
                event_overlap_residence_s=dwell, event_appreciable_core_residence_s=bright_dwell,
                event_first_xyz_m=first_xyz, event_first_time_s=np.where(near, time+start+first, np.nan),
                event_first_offset_s=np.where(near, start+first, np.nan), event_launch_time_s=time,
                event_parent_lag_s=dtime, event_launch_xy_m=np.stack((xy, xy2), axis=1),
                event_parent_species=np.stack((a['sp']+1, b['sp']+1), axis=1),
                inherited_opportunity_window_s=12*(arm['opt']+arm['ne']),
                unordered_supported_count_2s=valid_count, ordered_supported_count_2s=2*valid_count,
                unordered_core_count_2s=float(weight @ core),
                unordered_appreciable_core_count_2s=bright_count,
                unordered_supported_rate_s=valid_count/clock.seconds,
                jitter_opportunity_count_2s=float(weight.sum()),
                p_near_given_jitter=float((weight @ near)/weight.sum()),
                p_core_given_jitter=float((weight @ core)/weight.sum()),
                p_appreciable_core_given_jitter=float(bright_count/weight.sum()),
                mean_contact_residence_s=float((weight @ dwell)/valid_count) if valid_count else None,
                mean_appreciable_core_residence_s=float((weight @ bright_dwell)/bright_count) if bright_count else None,
                event_effective_near_samples=float(valid_count**2/np.sum((weight*near)**2)) if valid_count else 0.)
    if mean_strain is not None:
        return result
    # Independently sampled live anchors, never counted as additional events.
    # The geometry pass uses their conditional mean signed hydrostatic strain;
    # the response uses the actual signed values at each anchor.
    source_strain = np.stack((a['strain_birth'], b['strain_birth']), axis=1)
    conditional_mean = float(weight@source_strain.mean(axis=1)/weight.sum())
    width_product = _at(a, 'width', f1)**2+_at(b, 'width', f2)**2
    strength = np.where(distance < SUPPORT*summed_width, response*np.exp(-distance**2/width_product), 0.)
    ix = np.arange(size)
    peak = np.argmax(strength, axis=1)
    valid = near & (strength[ix, peak] > 0.) & (weight > 0.)
    event_age1, event_age2 = f1[ix, peak], f2[ix, peak]
    anchor = {
        'event_response_weight': weight*valid,
        'event_response_present': valid,
        'event_response_xyz_m': np.full((size, 3), np.nan),
        'event_response_time_s': np.where(valid, time+local[ix, peak], np.nan),
        'event_response_parent_flight_s': np.where(valid[:, None], np.stack((event_age1, event_age2), axis=1), np.nan),
        'event_response_parent_strain': source_strain,
        'event_response_cross_phase_rad': np.full((size, 2), np.nan),
    }
    scalar_fields = ('pair_energy_eV', 'ray_work_eV', 'pair_launch_work_eV', 'pair_damping_work_eV',
                     'incoherent_work_eV', 'energy_identity_defect_eV', 'v_joined_z_m_s',
                     'v_additive_z_m_s', 'v_frozen_z_m_s', 'v_strain_off_z_m_s')
    anchor.update({f'event_response_{key}': np.full(size, np.nan) for key in scalar_fields})
    if np.any(valid):
        a, b = ({key: val[valid] for key, val in parent.items()} for parent in (a, b))
        xy1, xy2 = xy[valid], xy2[valid]
        ta, tb = event_age1[valid], event_age2[valid]
        _launch_action(a)
        _launch_action(b)
        (own, cross, energy, ray_work, launch_work, damping, defect, pop_work, _) = _track_quartic(
            a, b, ta, tb, xy1, xy2, mat, depth, arm, marks)
        center1, ref1 = _centres(a, xy1, ta)
        center2, ref2 = _centres(b, xy2, tb)
        width1, width2 = _at(a, 'width', ta)[:, None], _at(b, 'width', tb)[:, None]
        xyz = (width2*center1+width1*center2)/(width1+width2)
        v1, vf1, vo1, *_ = _velocity(a, ta, xyz, center1, ref1, _action(a, ta, arm),
                                      own[:, 0], depth, mat, arm, marks)
        v2, vf2, vo2, *_ = _velocity(b, tb, xyz, center2, ref2, _action(b, tb, arm),
                                      own[:, 1], depth, mat, arm, marks)
        j1, j2 = np.exp(-1j*cross.T)
        fields = (energy, ray_work, launch_work, damping, pop_work, defect, (j1*v1+j2*v2).imag,
                  (v1+v2).imag, (j1*vf1+j2*vf2).imag, (j1*vo1+j2*vo2).imag)
        for key, value in zip(scalar_fields, fields):
            anchor[f'event_response_{key}'][valid] = value
        anchor['event_response_xyz_m'][valid] = xyz
        anchor['event_response_cross_phase_rad'][valid] = cross
    return conditional_mean, anchor


def _probe(rng, size, marks, mat, depth, arm, allocation, ec, beam, clock, norm, medium, strain):
    probe_xy = _xy(rng, size, beam)
    z, species, lab, column, probe_age = depth.probe(rng, size, marks)
    beam_time, time_pdf = clock.draw(rng, size, beam, 3)
    probe_time = beam_time+probe_age
    xyz = np.c_[probe_xy, z]
    a = _draw_marks(rng, size, marks, mat)
    b = _draw_marks(rng, size, marks, mat)
    f1, qt1 = _age_proposal(rng, a, arm)
    f2, qt2 = _age_proposal(rng, b, arm)
    age1, age2 = f1+a['formation'], f2+b['formation']
    t1, t2 = probe_time-age1, probe_time-age2
    a = _rays(a, marks, mat, arm, allocation, ec, beam, probe_xy, t1, medium, strain, controls=True)
    b = _rays(b, marks, mat, arm, allocation, ec, beam, probe_xy, t2, medium, strain, controls=True)
    radius1, radius2 = SUPPORT*_at(a, 'width', f1), SUPPORT*_at(b, 'width', f2)
    # Backward Palm proposal for both electrical launch positions. The age and
    # physical lateral area Jacobians are explicit. At fixed probe, the whole
    # two-parent density has its own third Cox intensity, the incident-probe flux.
    offset1, qp1 = _focus_disk(rng, size, radius1)
    offset2, qp2 = _focus_disk(rng, size, radius2)
    xy1 = probe_xy-_at(a, 'drift', f1)[:, :2]+offset1
    xy2 = probe_xy-_at(b, 'drift', f2)[:, :2]+offset2
    p1, ref1 = _centres(a, xy1, f1)
    p2, ref2 = _centres(b, xy2, f2)
    near1 = (np.linalg.norm(xyz-p1, axis=1) <= radius1) & (probe_time <= clock.seconds)
    near2 = (np.linalg.norm(xyz-p2, axis=1) <= radius2) & (probe_time <= clock.seconds)
    r = _rate(ec)
    w1 = r*_beam(xy1, beam)/norm * clock.value(t1, beam)/(qp1*qt1)*near1
    w2 = r*_beam(xy2, beam)/norm * clock.value(t2, beam)/(qp2*qt2)*near2
    probe_factor = clock.value(beam_time, beam)/(clock.seconds*time_pdf*size)
    weight = .5*probe_factor*w1*w2
    _launch_action(a)
    _launch_action(b)
    (self_phase, cross_phase, cross_e, cross_shape, cross_launch, cross_loss, defect,
     pair_pop_work, pair_pop_delta) = _track_quartic(a, b, f1, f2, xy1, xy2, mat, depth, arm, marks)
    amp1, amp2 = _action(a, f1, arm), _action(b, f2, arm)
    c1, c1f, c1o, n1, repl1, k1 = _velocity(a, f1, xyz, p1, ref1, amp1, self_phase[:, 0], depth, mat, arm, marks)
    c2, c2f, c2o, n2, repl2, k2 = _velocity(b, f2, xyz, p2, ref2, amp2, self_phase[:, 1], depth, mat, arm, marks)
    p1_cross, p2_cross = np.exp(-1j*cross_phase.T)
    j1, j2 = c1*p1_cross, c2*p2_cross
    shuffle1=np.exp(1j*rng.uniform(0.,2*np.pi,size));shuffle2=np.exp(1j*rng.uniform(0.,2*np.pi,size))
    s1=j1*shuffle1;s2=j2*shuffle2
    T, U, _, _ = medium(ec, beam, xyz, mount='nominal', time_s=np.minimum(probe_time, clock.seconds))
    quiet = KB*T
    work1, pop_debit1 = _parameter_work(a, f1, marks, mat, arm)
    work2, pop_debit2 = _parameter_work(b, f2, marks, mat, arm)
    cohort_scale1, cohort_scale2 = n1/(3*a['sites']), n2/(3*b['sites'])
    pop1 = n1*_native(a, age1, marks, 'theta')-repl1+cohort_scale1*pop_debit1
    pop2 = n2*_native(b, age2, marks, 'theta')-repl2+cohort_scale2*pop_debit2
    theta = quiet+pop1+pop2
    joined_theta = theta+cohort_scale1*pair_pop_delta[:, 0]+cohort_scale2*pair_pop_delta[:, 1]
    if min(theta.min(), joined_theta.min(), (quiet+pop1).min(), (quiet+pop2).min()) <= 0.:
        raise ValueError('Coherent replacement exhausts physical probe covariance')
    native_pop = np.stack((_native(a, age1, marks, 'optical'), _native(b, age2, marks, 'optical')), axis=1)
    first_donor = native_pop+np.stack((pop_debit1, pop_debit2), axis=1)
    pair_donor = first_donor+pair_pop_delta
    individual_support = np.stack((w1 > 0., w2 > 0.), axis=1)
    if np.any(first_donor[individual_support] < 0.) or np.any(pair_donor[weight > 0.] < 0.):
        raise ValueError('Reciprocal work exhausts its inherited optical donor')

    def vector(c):
        # Pass4 axial direct-D coherent mark; native stochastic population is isotropic.
        return np.c_[np.zeros((len(c), 2)), c.imag]

    # Independent probe, binary mixture of unordered pair occupation and its
    # first-order population. Pair kernels must use inclusion-exclusion with
    # their own conditional additive comparator, not joined minus quiet.
    return dict(xyz_m=xyz, probe_time_s=probe_time, probe_launch_time_s=beam_time,
                probe_species=species, probe_lab_keV=lab, probe_column_weight_m2=column,
                probe_temperature_K=T, probe_Ue_eV_formula=U, occupation_weight=weight,
                probe_time_proposal_Hz=time_pdf,
                single_occupation_weight=.5*probe_factor[:,None]*np.stack((w1, w2), axis=1),
                parent_age_proposal_Hz=np.stack((qt1, qt2), axis=1),
                parent_lateral_proposal_m_inverse2=np.stack((qp1, qp2), axis=1),
                conditional_source_measure_per_proposal_draw=np.stack((w1, w2), axis=1),
                parent1_id=2*np.arange(size), parent2_id=2*np.arange(size)+1,
                parent_species=np.stack((a['sp']+1, b['sp']+1), axis=1),
                parent_impact_age_s=np.stack((age1, age2), axis=1),
                parent_formation_s=np.stack((a['formation'], b['formation']), axis=1),
                parent_launch_time_s=np.stack((t1, t2), axis=1),
                parent_source_xyz_m=np.stack((np.c_[xy1, a['z']], np.c_[xy2, b['z']]), axis=1),
                parent_centroid_xyz_m=np.stack((p1, p2), axis=1), parent_reflection=np.stack((ref1, ref2), axis=1),
                parent_wavevector_m_inverse=np.stack((k1, k2), axis=1),
                parent_amplitude_width_m=np.stack((radius1, radius2), axis=1)/SUPPORT,
                parent_molecular_echo_s=np.stack((a['delay'], b['delay']), axis=1),
                parent_echo_present=np.stack(((a['sp'] == 1) & (f1 >= a['delay']),
                                              (b['sp'] == 1) & (f2 >= b['delay'])), axis=1),
                parent_covariance_projection_eV=np.stack((pop1, pop2), axis=1),
                parent_launch_signed_work_eV=np.stack((a['launch_work'], b['launch_work']), axis=1),
                parent_launch_quiet_cross_eV=np.stack((a['launch_quiet_cross'], b['launch_quiet_cross']), axis=1),
                parent_first_order_carrier_work_eV=np.stack((work1, work2), axis=1),
                parent_reciprocal_optical_population_eV=np.stack((pop_debit1, pop_debit2), axis=1),
                parent_pair_reciprocal_optical_population_eV=pair_pop_delta,
                theta_quiet_eV=quiet, theta_parent1_eV=quiet+pop1, theta_parent2_eV=quiet+pop2,
                theta_additive_eV=theta, theta_joined_eV=joined_theta,
                theta_frozen_eV=joined_theta, theta_shuffled_eV=joined_theta,
                v_parent1_m_s=vector(c1), v_parent2_m_s=vector(c2), v_additive_m_s=vector(c1+c2),
                v_parent1_frozen_m_s=vector(c1f), v_parent2_frozen_m_s=vector(c2f),
                v_parent1_strain_off_m_s=vector(c1o), v_parent2_strain_off_m_s=vector(c2o),
                v_parent1_shuffled_m_s=vector(c1*shuffle1), v_parent2_shuffled_m_s=vector(c2*shuffle2),
                v_joined_m_s=vector(j1+j2), v_frozen_m_s=vector(c1f*p1_cross+c2f*p2_cross),
                v_strain_off_m_s=vector(c1o*p1_cross+c2o*p2_cross),
                v_shuffled_m_s=vector(s1+s2), cross_phase_rad=cross_phase,
                pair_projected_cross_energy_eV=cross_e, optical_shape_work_demand_eV=cross_shape,
                pair_incoherent_to_coherent_work_eV=pair_pop_work,
                optical_pair_launch_work_eV=cross_launch, optical_pair_loss_to_thermal_e_eV=.2*cross_loss,
                optical_pair_loss_to_Pd_eV=.8*cross_loss, phase_hamiltonian_product_defect_eV=defect,
                unresolved_k_beat=(np.linalg.norm(k1-k2, axis=1)*np.sqrt(radius1**2+radius2**2)/SUPPORT < 6.),
                coherent_joined_quadrature_m_s=np.stack((j1, j2), axis=1),
                coherent_uncoupled_quadrature_m_s=np.stack((c1, c2), axis=1),
                joined_minus_off_strain_birth=np.stack((a['strain_birth'], b['strain_birth']), axis=1),
                max_fractional_strain_pole_change=np.stack((a['max_fractional_strain_pole_change'],
                                                          b['max_fractional_strain_pole_change']), axis=1),
                max_fractional_warm_pole_change=np.stack((a['max_fractional_warm_pole_change'],
                                                        b['max_fractional_warm_pole_change']), axis=1))


def _summary(probe, event, prior, ec, beam, allocation, closure, clock, norm, square, source_allocation=None):
    w = probe['occupation_weight']
    occupation = float(w.sum())
    vj, va = probe['v_joined_m_s'][:, 2], probe['v_additive_m_s'][:, 2]
    sw = probe['single_occupation_weight']
    peak = np.flatnonzero(w > 0.)
    donor = probe['parent_first_order_carrier_work_eV']
    relevant_donor = donor[sw > 0.]
    launch = abs(probe['parent_launch_signed_work_eV']).sum(axis=1)
    cross_scale = np.maximum(launch, 1e-300)
    pc = 'max_fractional_strain_pole_change'
    reported = next((r for r in prior['whole_history_spacetime_opportunities']
                     if r['arm'] == f'{closure}_{source_allocation or allocation}' and r['beam'] == beam), None)
    base_rate = _rate(ec)
    clock2 = clock.integrals[2] if beam == 'hot_bursts' else clock.seconds
    historical_opportunity = (base_rate**2 * square/norm**2 * clock2 *
                              np.pi*(SUPPORT*ARMS[closure]['radius'])**2 * event['inherited_opportunity_window_s'])
    return dict(E_cm_keV=ec, beam=beam, allocation=allocation, closure=closure, electrical_parent_rate_s=base_rate,
                Pass4_directed_window_label_2p5_2s=reported['whole_foil_directed_spacetime_window_opportunities_2s'],
                historical_directed_window_recomputed_same_E_2s=float(historical_opportunity),
                clock_square_mean=float(clock2/clock.seconds),
                clock_cube_mean=float(clock.integrals[3]/clock.seconds if beam == 'hot_bursts' else 1.),
                **{k: v for k, v in event.items() if np.ndim(v) == 0},
                unordered_probe_occupation_per_incident_probe=occupation,
                single_probe_occupation_per_incident_probe=float(sw.sum()),
                pair_effective_probe_samples=float(occupation**2/(w@w)) if occupation else 0.,
                pair_occupation_MC_standard_error=float(np.sqrt(len(w)*np.var(w, ddof=1))),
                unordered_count_MC_standard_error=float(np.sqrt(len(event['event_weight'])*np.var(
                    event['event_weight']*event['event_near'], ddof=1))),
                rms_joined_minus_additive_velocity_m_s=float(np.sqrt(w@((vj-va)**2)/occupation)) if occupation else None,
                rms_joined_minus_additive_theta_eV=float(np.sqrt(w@((probe['theta_joined_eV']-
                                                                   probe['theta_additive_eV'])**2)/occupation)) if occupation else None,
                reflected_parent_occupation_fraction=float(np.sum(sw*(probe['parent_reflection']<0))/sw.sum())
                if sw.sum() else None,
                max_supported_optical_flight_m=float(np.max(np.linalg.norm(probe['parent_centroid_xyz_m']-
                                                                           probe['parent_source_xyz_m'], axis=-1)[sw > 0.], initial=0.)),
                largest_supported_cross_phase_rad=float(np.max(abs(probe['cross_phase_rad'][peak]), initial=0.)),
                near_degenerate_beat_occupation_fraction=float(w@probe['unresolved_k_beat']/occupation) if occupation else None,
                weighted_mean_pair_cross_energy_eV=float(w@probe['pair_projected_cross_energy_eV']/occupation) if occupation else None,
                max_pair_incoherent_to_coherent_work_eV=float(np.max(
                    abs(probe['pair_incoherent_to_coherent_work_eV'][peak]), initial=0.)),
                max_overlap_shape_work_demand_eV=float(np.max(abs(probe['optical_shape_work_demand_eV'][peak]), initial=0.)),
                max_pair_launch_work_eV=float(np.max(abs(probe['optical_pair_launch_work_eV'][peak]), initial=0.)),
                max_overlap_shape_work_fraction_of_signed_launch=float(np.max(
                    abs(probe['optical_shape_work_demand_eV'][peak])/cross_scale[peak], initial=0.)),
                max_overlap_action_damping_eV=float(np.max((probe['optical_pair_loss_to_thermal_e_eV']+
                                                            probe['optical_pair_loss_to_Pd_eV'])[peak], initial=0.)),
                max_pair_H_split_defect_eV=float(np.max(abs(probe['phase_hamiltonian_product_defect_eV'][peak]), initial=0.)),
                max_first_order_H_split_defect_eV=float(np.max(abs(relevant_donor[:, -3]), initial=0.)),
                max_first_order_electronic_donor_fraction=float(np.max(relevant_donor[:, -2], initial=0.)),
                max_first_order_incoherent_donor_fraction=float(np.max(relevant_donor[:, -1], initial=0.)),
                max_strain_fractional_pole_change=float(np.max(probe[pc][sw > 0.], initial=0.)),
                max_warm_fractional_pole_change=float(np.max(probe['max_fractional_warm_pole_change'][sw > 0.], initial=0.)),
                rms_joined_minus_strain_off_velocity_m_s=float(np.sqrt(w@((vj-probe['v_strain_off_m_s'][:, 2])**2)/
                                                                        occupation)) if occupation else None)


def build(out, path4, medium_sampler, strain_sampler, *, probe_samples=PROBE_SAMPLES, event_samples=EVENT_SAMPLES,
          event_response_samples=EVENT_RESPONSE_SAMPLES,case_selection=None,prefix='',allocations=ALLOCATIONS,
          source_allocation=None,material_factory=Material,extra_scientific_sources=()):
    """Write the case grid with separate geometry and signed-response designs."""
    if event_samples < 2 or not 2 <= min(probe_samples, event_response_samples) or max(
            probe_samples, event_response_samples) > MAX_LIVE_SAMPLES:
        raise ValueError('Each estimator requires ≥2 proposals; live response counts must be ≤96 per case')
    out = Path(out)
    accepted = Accepted(path4)
    inputs = accepted.json('frozen/pass1_inputs.json')
    prior = accepted.json('decision.json')
    tracks = accepted.json('frozen/tracks.json')
    debits = {a['name']: a['source_debits'] for a in accepted.json('native_summary.json')['arms']}
    clock = Clock(accepted('calendar.npz'))
    material = material_factory(accepted, inputs['lattice'])
    summaries, cases = [], []
    beam_moments = {beam: _beam_norm(beam) for beam in BEAMS}
    epoch_names=('meetings.py','medium.py','live.py','receiver.py','low.py','common.py','medium.npz','low_modes.npz')
    epoch=hashlib.sha256(''.join(_digest(out/name) for name in (*epoch_names,*extra_scientific_sources)).encode()).hexdigest()
    epoch+=f':{probe_samples}:{event_samples}:{event_response_samples}'
    checkpoints=[]
    for ei, ec in enumerate(ENERGIES):
        depth = Depth(next(r for r in inputs['profiles'] if r['case']['E_cm_keV'] == ec))
        for ai, (closure, arm) in enumerate(ARMS.items()):
            for ci, allocation in enumerate(allocations):
                launch_allocation=source_allocation or allocation
                name = f'{closure}_{launch_allocation}'
                marks = [Mark(tracks[f'{ec:g}_{sp}'], accepted(f'native_E{ec:g}_{name}_sp{sp}_m16.npz'),
                              sp, arm, debits[name][f'{ec:g}_{sp}'], material) for sp in (0, 1)]
                for bi, beam in enumerate(BEAMS):
                    if case_selection is not None and (ec,closure,allocation,beam) not in case_selection:continue
                    norm, square = beam_moments[beam]
                    checkpoint=out/f'.{prefix}meeting_E{ec:g}_{closure}_{allocation}_{beam}.npz'
                    checkpoints.append(checkpoint)
                    if checkpoint.exists():
                        with np.load(checkpoint) as cache:
                            if str(cache['_epoch'])==epoch:
                                summaries.append(json.loads(str(cache['_summary'])))
                                cases.append({k:cache[k] for k in cache.files if not k.startswith('_')})
                                continue
                    def case_strain(*args, **kwargs):
                        return strain_sampler(*args, **kwargs, closure=closure)
                    # Independent response anchors and probe Palm proposals;
                    # keep the larger geometry proposal stream and normalization.
                    seed = [975215, ei, ai, bi]
                    geom_rng = np.random.default_rng(np.random.SeedSequence(seed))
                    anchor_rng = np.random.default_rng(np.random.SeedSequence([*seed, 1]))
                    probe_rng = np.random.default_rng(np.random.SeedSequence([*seed, 2]))
                    strain_mean, anchor = _events(anchor_rng, event_response_samples, marks, material, depth, arm,
                                                  allocation, ec, beam, clock, norm, medium_sampler, case_strain)
                    event = _events(geom_rng, event_samples, marks, material, depth, arm, allocation, ec, beam,
                                    clock, norm, medium_sampler, case_strain, mean_strain=strain_mean)
                    probe = _probe(probe_rng, probe_samples, marks, material, depth, arm, allocation, ec, beam,
                                   clock, norm, medium_sampler, case_strain)
                    summary = _summary(probe, event, prior, ec, beam, allocation, closure, clock, norm, square,launch_allocation)
                    summary['launch_source_allocation']=launch_allocation
                    matched_anchor = anchor['event_response_present']
                    summary.update(geometry_conditional_mean_linear_strain=strain_mean,
                                   event_response_live_supported_draws=int(matched_anchor.sum()),
                                   event_response_supported_measure_2s=float(anchor['event_response_weight'].sum()))
                    summaries.append(summary)
                    cases.append({**probe, **anchor, **{k: v for k, v in event.items() if np.ndim(v)},
                                  'E_cm_keV': ec, 'beam': beam, 'closure': closure, 'allocation': allocation,'source_allocation':launch_allocation})
                    np.savez_compressed(checkpoint,**cases[-1],_epoch=np.array(epoch),
                                        _summary=np.array(json.dumps(summary,allow_nan=False)))
                    print(f'meeting E={ec:g} {closure} {allocation} {beam}',flush=True)
    # P4's numbers are directed windows; match its construction only, not an
    # assumption that 1/2 of its asymmetric source windows are physical events.
    selected = [c for c in summaries if c['E_cm_keV'] == 2.5]
    prior_error = max((abs(c['historical_directed_window_recomputed_same_E_2s']/
                          c['Pass4_directed_window_label_2p5_2s']-1) for c in selected),default=None)
    if prior_error is not None and prior_error > 5e-12:
        raise ValueError('The inherited current/beam/Cox opportunity identity does not close')
    prior_by_arm = {}
    for closure in ARMS:
        counts = [c['Pass4_directed_window_label_2p5_2s'] for c in selected if c['closure'] == closure]
        if counts:prior_by_arm[closure] = [min(counts), max(counts)]
    np.savez_compressed(out/(prefix+'meetings.npz'), **{k: np.stack([r[k] for r in cases]) for k in cases[0]})
    report = dict(schema='whole-foil-optical-meetings-v1', accepted_pass4_inputs=accepted.used,
                  historical_opportunity_count_relative_error=prior_error,scientific_source_epoch=epoch,
                  extra_scientific_sources=list(extra_scientific_sources),
                  frozen_incoming_source_allocation=source_allocation,
                  historical_directed_window_ranges_2s_by_closure=prior_by_arm,
                  Pass4_maximum_wave_modulation=prior['maximum_wave_modulation'],
                  Pass4_maximum_thermal_frequency_modulation=prior['maximum_thermal_frequency_modulation'],
                  Pass4_observed_other_parents=prior['observed_independent_optical_overlap_parents'],
                  Pass4_sampled_parent_windows=prior['physical_optical_windows'],
                  definition={
                      'event': 'Campbell factorial measure on chronological electrical launches: R² dt dlag B(x1) B(x2) dx1 dx2 C(t1) C(t2); 2 times chronological count gives ordered count. Count a pair once on first 3-sigma support intersection after causal formation; later reflected episodes add residence, not events. The launch-lag and offset disk are an explicit jitter denominator, not the historical Pass4 source cutoff. 1-sigma core and appreciable core additionally distinguish a normalized native optical or direct-response level e^-3 on both parents, not a measured physical cut.',
                      'sampling': 'Rate/geometry uses independent proposals at the live-anchor pair-Palm mean linear strain, retaining resolved whole-foil temperature/electronic refraction and identical Campbell/Jacobian factors. The mean closure for geometry does not set strain to zero or estimate instantaneous displacement. Independent live anchor draws evolve signed response through the strongest sampled overlap; their weights are not added to the geometry counts. Independent incident-probe draws sample both parents at their live source times and carry their own sample normalization. Geometry Monte Carlo errors are conditional on the estimated mean strain.',
                      'probe': 'Independent incident track probe with own .85/.30 D weighting, exact D(z) column and flight lag. Unordered pair occupation = 1/2 factorial moment of two simultaneously supporting distinct parents at that probe. Proposal factor is C(t_probe_launch)/(2s q(t)), times two launch factors R B(x)C(t_parent)/[q(x)q(age)]; q(x) is a nested-disk area mixture, q(age) is a uniform/physical-decay mixture. Sample mean normalization is in weights. Use probe_column_weight_m2 for G12 and use single_occupation_weight for matched one-parent replacement. No incident probe coincides with a parent electrical mark.',
                      'connected': 'Use F(joined)-F(parent1)-F(parent2)+F(quiet) at identical probe lab energy and T. F(joined)-F(additive) isolates projected nonlinear dynamics; do not count ordinary one-parent population twice. theta is isotropic moving-target thermal population, v is signed axial canonical velocity, not covariance.',
                      'transport': 'Accepted tracking/debits, actual D profile, native m16 phase-space populations by Bloch row, reflected front/back shell centre with unit free optical image energy, true source flight offsets and same-parent D2 echo. Retained angularly isotropized Hessian spreading and first-order refractive ray gradients about the common whole-foil chart; launches within the jitter disk use the same lateral Taylor chart. The archived integrated native phase resolves the fs quench while the warm innovation uses the same transported foil. Phase and arrival clocks are integrated separately to avoid subtracting nearby whole seconds for femtosecond phase.',
                      'interaction': 'Same native K4_100 axial RWA self and cross terms, plus G3 Hartree contraction with the other parent’s actual incoherent q variance. Phase is the common H action derivative, canonical momentum is measured by Ω Im(q), and no pair amplitude gain is inserted. All cross terms use the same D-weighted Gaussian overlap with continuous compact support; close-k exchange and FWM remain unretained.',
                      'controls': 'Additive is the same two responses with cross quartic disabled. Frozen replaces the linear phase by a cold matched phase at equal source phasor, observed geometry, population and attenuation. Shuffled adds independent constant phases to complete parent phasors, not to D2 or reflected paths. Pure dephasing coherently attenuates and is not lost from the covariance.',
                      'work': 'First-order parent-clock curvature work (eV) enters from that parent’s hot e, thermal e or incoherent optical reservoirs with opposite donor amounts; spatial background refraction is separate. Damping routes .2 to thermal e, .8 to Pd. The pair quartic separately charges parent insertion and echoes, physical damping to the same carriers, and coefficient-change work to the optical ray/shape port. Opposite signed native q-population work is propagated back through the optical population Green to the probe covariance, including the inherited hot-e to optical path. Prescribed rays are the leading eikonal approximation: opposite ray-work demand prices omitted backreaction, not independent heat or pump credit.',
                      'limits': 'Conditional model, not realized hits, microscopic phase measurement, nuclear contact probability, or a physical tail bound. Gaussian quartic is phase-averaged; ray centroids reflect under a named free optical scenario rather than a measured reflection law, and this does not resolve simultaneous near-wall image interference. Live acoustic background conditions the optical Palm insertions, whose additional tagged acoustic source is not replayed. Shapes outside 3 sigma, strain gradients along each ps ray, optical mode exchange and microscopic phase-dependent covariance are unretained. The accepted fixed sheet remains quiet once.',
                      'carrier_ledger_columns': ['hot_e_to_coherent', 'thermal_e_to_coherent', 'incoherent_optical_to_coherent',
                                                 'background_refractive_shape_work', 'same_parent_echo_to_coherent',
                                                 'coherent_loss_to_thermal_e', 'coherent_loss_to_Pd', 'coherent_energy_balance_defect',
                                                 'max_positive_electronic_donor_fraction', 'max_positive_optical_donor_fraction'],
                  }, geometry_event_samples_per_case=event_samples, event_response_samples_per_case=event_response_samples,
                  probe_samples_per_case=probe_samples,
                  strain_chart='Whole-beam conditional Gaussian/direct exterior Green plus globally marked physical near rain and retained low prompt/pressure Cox mean. Prompt innovations have a correlation-independent bound; low thermoelastic continuation and exterior return covariance remain unresolved. Matched strain_off is zero. Not a measured phase chart.',
                  cases=summaries)
    (out/(prefix+'meetings.json')).write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    for checkpoint in checkpoints:
        checkpoint.unlink(missing_ok=True)
    return report


def merge(out,prefix=''):
    """Assemble disjoint saved energy jobs, never regenerate a physical case."""
    out=Path(out);arrays=[];reports=[]
    for ec in ENERGIES:
        source_prefix=f'{prefix}E{ec:g}_'
        with np.load(out/(source_prefix+'meetings.npz')) as data:
            if set(np.asarray(data['E_cm_keV'],float))!={ec}:raise ValueError('Meeting energy shard mismatch')
            arrays.append({k:data[k] for k in data.files})
        reports.append(json.loads((out/(source_prefix+'meetings.json')).read_text()))
    if len({r['scientific_source_epoch'] for r in reports})!=1:raise ValueError('Mixed scientific source epochs')
    if any(set(a)!=set(arrays[0]) for a in arrays):raise ValueError('Meeting field contract changed between shards')
    np.savez_compressed(out/(prefix+'meetings.npz'),**{k:np.concatenate([a[k] for a in arrays]) for k in arrays[0]})
    report=reports[-1].copy();report['cases']=[c for r in reports for c in r['cases']]
    report['accepted_pass4_inputs']={k:v for r in reports for k,v in r['accepted_pass4_inputs'].items()}
    report['assembly']='six disjoint exact-energy jobs; identical scientific epoch and proposal counts'
    (out/(prefix+'meetings.json')).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report
