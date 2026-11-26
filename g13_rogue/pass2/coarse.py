"""Whole-foil finite-band variational Pd model; run from any working directory."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import tarfile

import numpy as np
from numpy.polynomial.legendre import leggauss, legval, legder
from scipy.constants import elementary_charge as QE, Boltzmann as KB, atomic_mass
from scipy.linalg import eigh
from scipy.spatial.transform import Rotation
from scipy.special import exp1

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASE = '9c26f8cb688918ae18b63dbd187fdf99793c1e90'
# Equal-area square surrogate for the 24-mm-diameter disk; central 1-cm² spot.
DIMS = np.array([math.sqrt(math.pi)*.012, math.sqrt(math.pi)*.012, .00025])
RHO = 12007.
C11, C12, C44 = 227.1e9, 176.04e9, 71.73e9
REFERENCE = .0006423
SOURCES = ('g4_beam/merit/profiles.json', 'g4_beam/forward/model.py',
           'g11_cmns/collective/part1/stopping.json', 'g13_rogue/depth.py',
           'g12_residual/inputs_ledger.md')


def table(path, rows):
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def profiles_and_stopping():
    profiles = json.loads((ROOT / SOURCES[0]).read_text())['records']
    profiles = {(r['case']['condition'], float(r['case']['E_cm_keV'])): r for r in profiles}
    stopping = json.loads((ROOT / SOURCES[2]).read_text())
    stopping = {(r['target'], float(r['E_cm_keV'])): r for r in stopping}
    return profiles, stopping


def layer_value(layers, z, column=3):
    for l in layers:
        if l[0] <= z < l[1]:
            return l[column]
    if abs(z - layers[-1][1]) < 1e-18:
        return layers[-1][column]
    raise ValueError('Outside fixed profile')


def deposition_depth(track):
    """Host-recoil deposition first moment; timeline depths are actually Angstrom."""
    cumulative = np.array([0.] + [t['host_recoil_eV'] for t in track['timeline']])
    z = np.array([0.] + [t['depth_A']*1e-10 for t in track['timeline']])
    increments = np.diff(cumulative)
    return float(np.dot(increments, (z[1:]+z[:-1])/2)/cumulative[-1])


def column_density(layers, z0, z1, column=3):
    return sum(max(0., min(z1, l[1]) - max(z0, l[0])) * l[column] for l in layers) / (z1-z0)


def basis(points, degree=(5, 5, 2)):
    """Physical L2 mass-orthonormal scalar Legendre basis, with derivatives."""
    values, derivs = [], []
    for axis, n in enumerate(degree):
        xi = 2 * points[:, axis] / DIMS[axis] - 1
        vv, dd = [], []
        for k in range(n+1):
            coeff = np.zeros(k+1)
            coeff[k] = 1
            norm = math.sqrt((2*k+1) / (RHO * DIMS[axis]))
            # Rho is included once in 3-D below, rather than three times.
            norm *= RHO**(1/3)
            vv.append(legval(xi, coeff)*norm)
            dd.append(legval(xi, legder(coeff))*norm*2/DIMS[axis])
        values.append(np.array(vv).T)
        derivs.append(np.array(dd).T)
    phi = np.einsum('qi,qj,qk->qijk', *values).reshape(len(points), -1)
    grad = []
    for axis in range(3):
        factors = [derivs[a] if a == axis else values[a] for a in range(3)]
        grad.append(np.einsum('qi,qj,qk->qijk', *factors).reshape(len(points), -1))
    n = phi.shape[1]
    b = np.zeros((len(points), 6, 3*n))
    b[:, 0, :n] = grad[0]
    b[:, 1, n:2*n] = grad[1]
    b[:, 2, 2*n:] = grad[2]
    # Engineering shear convention: epsilon vector (xx,yy,zz,2yz,2xz,2xy).
    b[:, 3, n:2*n], b[:, 3, 2*n:] = grad[2], grad[1]
    b[:, 4, :n], b[:, 4, 2*n:] = grad[2], grad[0]
    b[:, 5, :n], b[:, 5, n:2*n] = grad[1], grad[0]
    return phi, b


def quadrature(grain, order=6):
    axes, weights = [], []
    for axis, length in enumerate(DIMS):
        cells = max(1, round(length/grain)) if axis < 2 else 1
        x, w = leggauss(max(order, 8 if cells == 1 and axis < 2 else order))
        edges = np.linspace(0., length, cells+1)
        axes.append(np.concatenate([(a+b)/2+(b-a)/2*x for a, b in zip(edges[:-1], edges[1:])]))
        weights.append(np.concatenate([(b-a)/2*w for a, b in zip(edges[:-1], edges[1:])]))
    points = np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1).reshape(-1, 3)
    ww = np.einsum('i,j,k->ijk', *weights).ravel()
    return points, ww


def cubic(rotation=None, isotropic=False):
    if isotropic:
        bulk = (C11+2*C12)/3
        gv = (C11-C12+3*C44)/5
        gr = 5*(C11-C12)*C44/(4*C44+3*(C11-C12))
        mu = (gv+gr)/2
        lam = bulk-2*mu/3
        anis = 0.
    else:
        lam, mu, anis = C12, C44, C11-C12-2*C44
    eye = np.eye(3)
    c = lam*np.einsum('ij,kl->ijkl', eye, eye)
    c += mu*(np.einsum('ik,jl->ijkl', eye, eye)+np.einsum('il,jk->ijkl', eye, eye))
    rot = np.eye(3) if rotation is None else rotation
    c += anis*np.einsum('ia,ja,ka,la->ijkl', rot, rot, rot, rot)
    pairs = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))
    return np.array([[c[i, j, k, l] for k, l in pairs] for i, j in pairs])


def tensors(points, grain, kind):
    if kind != 'grains':
        return np.broadcast_to(cubic(isotropic=kind == 'isotropic'), (len(points), 6, 6)).copy()
    # Explicit synthetic columnar grains; same seed defines the orientation map.
    cells = round(DIMS[0]/grain)
    labels = np.minimum((points[:, :2]/grain).astype(int), cells-1)
    rotations = Rotation.random(cells*cells, random_state=173).as_matrix()
    cs = np.array([cubic(r) for r in rotations])
    return cs[labels[:, 0]*cells+labels[:, 1]]


def assemble(kind, grain, fmin, fmax):
    points, w = quadrature(grain)
    phi, b = basis(points)
    c = tensors(points, grain, kind)
    cb = np.einsum('qij,qjk->qik', c, b)
    k = np.einsum('qij,qik,q->jk', b, cb, w, optimize=True)
    evals, vec = eigh(k, check_finite=False)
    freq = np.sqrt(np.maximum(evals, 0))/(2*np.pi)
    keep = (freq >= fmin) & (freq <= fmax)
    v = vec[:, keep]
    if v.shape[1] < 4:
        raise ValueError('Retained band has too few modes')
    omega = 2*np.pi*freq[keep]
    mass = phi.T @ (w[:, None]*phi)*RHO
    residual = np.linalg.norm(k@v-v*omega**2)/np.linalg.norm(v*omega**2)
    return dict(kind=kind, grain=grain, v=v, omega=omega, points=points, w=w,
                c=c, b=b, residual=residual, mass_error=float(np.max(abs(mass-np.eye(len(mass))))),
                rigid_count=int(np.count_nonzero(freq < 1.)), ndof=len(k))


def source_matrix(model, depths):
    """Retained-band effective signed traction packets with equal specified energy."""
    x = (DIMS[0]-.01)/2+(.5+np.arange(6))*.01/6
    xy = np.stack(np.meshgrid(x, x, indexing='ij'), axis=-1).reshape(-1, 2)
    columns = []
    nscalar = model['v'].shape[0]//3
    for depth in depths:
        pts = np.column_stack([xy, np.full(len(xy), depth)])
        phi, _ = basis(pts)
        s = model['v'][2*nscalar:].T @ phi.T
        norms = np.linalg.norm(s, axis=0)
        if np.any(norms == 0):
            raise ValueError('Unsupported packet source')
        columns.append(s/norms)
    return np.concatenate(columns, axis=1)


def receiver(model, side, z0, z1):
    # Separate exact volume quadrature; does not identify quadrature points as grains.
    ax, aw = leggauss(6)
    bounds = ((DIMS[0]/2-side/2, DIMS[0]/2+side/2),
              (DIMS[1]/2-side/2, DIMS[1]/2+side/2), (z0, z1))
    xyz, ww = [], []
    for axis, (lo, hi) in enumerate(bounds):
        cuts = [lo, hi]
        if axis < 2 and model['kind'] == 'grains':
            cuts += [v for v in np.arange(model['grain'], DIMS[axis], model['grain']) if lo < v < hi]
        cuts = sorted(cuts)
        xyz.append(np.concatenate([(a+b)/2+(b-a)*ax/2 for a, b in zip(cuts[:-1], cuts[1:])]))
        ww.append(np.concatenate([aw*(b-a)/2 for a, b in zip(cuts[:-1], cuts[1:])]))
    pts = np.stack(np.meshgrid(*xyz, indexing='ij'), axis=-1).reshape(-1, 3)
    weights = np.einsum('i,j,k->ijk', *ww).ravel()
    _, b = basis(pts)
    c = tensors(pts, model['grain'], model['kind'])
    bm = b @ model['v']
    kr = np.einsum('qij,qik,q->jk', bm, np.einsum('qij,qjk->qik', c, bm), weights, optimize=True)
    vol = side*side*(z1-z0)
    # A q-coordinate strain-energy density quadratic form, then x=omega*q/sqrt(Eunit).
    a = kr / np.outer(model['omega'], model['omega']) / vol
    return (a+a.T)/2


def localized_strain(model):
    # Front-side bending strain, not epsilon_zz at the neutral plane.
    _, b = basis(np.array([[DIMS[0]/2, DIMS[1]/2, 30e-6]]))
    raw = b[0, 0] @ model['v'] / model['omega']
    norm = np.linalg.norm(raw)
    return raw/norm, norm


def evolve(model, s, kicks, dt, tau, chi, eunit, initial=None):
    """Exact harmonic rotation, positive rank-one quartic Strang kicks, same drive."""
    omega = model['omega']
    a, strain_norm = localized_strain(model)
    ntime, nrep, _ = kicks.shape
    x = np.zeros((nrep, len(omega))) if initial is None else initial[0].copy()
    y = np.zeros_like(x) if initial is None else initial[1].copy()
    history = np.empty((ntime, nrep, len(omega)))
    co, si = np.cos(omega*dt), np.sin(omega*dt)
    loss = math.exp(-dt/(2*tau)) if math.isfinite(tau) else 1.
    for t in range(ntime):
        y += kicks[t] @ s.T / math.sqrt(eunit)
        if chi:
            y -= .5*dt*omega*chi*a*(x@a)[:, None]**3
        x, y = x*co+y*si, y*co-x*si
        if chi:
            y -= .5*dt*omega*chi*a*(x@a)[:, None]**3
        x *= loss
        y *= loss
        history[t] = x
    energy = eunit*(.5*np.sum(x*x+y*y, axis=1)+chi/4*(x@a)**4)
    return history, x, y, energy, strain_norm


def tail_rows(u, thermal_mean, thermal_sd, side, z0, z1, arm, case, dt):
    # u axes time, realization. Connected excursions, no independent-bin fiction.
    mean = float(u.mean())
    exposure = u.size*dt
    rows = []
    for mult in (1., 2., 4., 8., 16.):
        mask = u >= mult*mean
        runs = []
        censored = 0
        for j in range(u.shape[1]):
            padded = np.r_[False, mask[:, j], False].astype(int)
            starts = np.where(np.diff(padded) == 1)[0]
            ends = np.where(np.diff(padded) == -1)[0]
            for start, end in zip(starts, ends):
                if start == 0 or end == len(mask):
                    censored += 1
                else:
                    runs.append((end-start)*dt)
        count = len(runs)
        rows.append(dict(case=case, arm=arm, side_mm=side*1e3, z0_um=z0*1e6, z1_um=z1*1e6,
            threshold_over_beam_mean=mult, beam_mean_J_m3=mean, peak_J_m3=float(u.max()),
            peak_over_mean=float(u.max()/mean), thermal_mean_J_m3=thermal_mean,
            thermal_sd_J_m3=thermal_sd, peak_over_thermal_mean=float(u.max()/thermal_mean),
            peak_over_thermal_sd=float(u.max()/thermal_sd), time_fraction=float(mask.mean()),
            independent_record_time_fraction_min=float(mask.mean(axis=0).min()),
            independent_record_time_fraction_max=float(mask.mean(axis=0).max()),
            independent_record_fraction_standard_error=float(mask.mean(axis=0).std(ddof=1)/math.sqrt(u.shape[1])) if u.shape[1] > 1 else 0.,
            complete_excursions=count, boundary_censored_excursions=censored,
            sampled_receiver_seconds=exposure, rate_per_receiver_s=count/exposure,
            zero_count_95pct_rate_upper_s=(-math.log(.05)/exposure if count == 0 else ''),
            duration_ge_1us_excursions=sum(v >= 1e-6 for v in runs),
            duration_ge_10us_excursions=sum(v >= 10e-6 for v in runs),
            duration_ge_50us_excursions=sum(v >= 50e-6 for v in runs),
            duration_ge_1us_rate_s=sum(v >= 1e-6 for v in runs)/exposure,
            duration_ge_10us_rate_s=sum(v >= 10e-6 for v in runs)/exposure,
            duration_ge_50us_rate_s=sum(v >= 50e-6 for v in runs)/exposure,
            duration_median_us=(float(np.median(runs))*1e6 if runs else ''),
            duration_max_us=(max(runs)*1e6 if runs else ''), duration_resolution_us=dt*1e6))
    return rows


def union_tail_rows(u, arm, case, dt):
    rows = tail_rows(u, 1., 1., 0., 0., 0., arm, case, dt)
    result = []
    for r in rows:
        item = {k: v for k, v in r.items() if k not in (
            'side_mm', 'z0_um', 'z1_um', 'threshold_over_beam_mean', 'beam_mean_J_m3',
            'peak_J_m3', 'thermal_mean_J_m3', 'thermal_sd_J_m3',
            'peak_over_thermal_mean', 'peak_over_thermal_sd')}
        item.update(threshold_over_bank_statistic_mean=r['threshold_over_beam_mean'],
                    bank_statistic_mean_dimensionless=r['beam_mean_J_m3'],
                    bank_statistic_peak_dimensionless=r['peak_J_m3'])
        result.append(item)
    return result


def backward(profiles, stopping, gain):
    rows, frontier, attenuation = [], [], []
    speed, tau_d = 4000., 4.1e-12
    for material in ('Pd', 'PdD'):
        for energy, current in ((.25, 38e-6), (2.5, 138e-6)):
            rec = profiles[material, energy]
            power = current*stopping[material, energy]['host_recoil_J_per_true_C']
            for depth in (30e-6, 50e-6, 200e-6):
                nd = layer_value(rec['layers'], depth)
                # Same hemisphere direct kernel as depth.py; source gain is an explicit sensitivity.
                ls = np.geomspace(1e-7, 1e3, 1800)
                edge = math.hypot(depth, math.sqrt(1e-4/math.pi))
                kernel = exp1(depth/ls)-exp1(edge/ls)
                u = power/1e-4/speed*kernel
                converted = u*speed/ls*tau_d/nd/QE
                imax = int(np.argmax(converted))
                emax = float(converted[imax])
                length0 = 100e-6
                k0 = float(exp1(depth/length0)-exp1(edge/length0))
                e0 = power/1e-4/length0*k0*tau_d/nd/QE
                # eta_band=1, transfer=1 are transparent maximum fractions, not coherent yield.
                required = REFERENCE/e0
                rows.append(dict(material=material, E_cm_keV=energy, depth_um=depth*1e6,
                    n_D_m3=nd, x_D_per_Pd=nd/layer_value(rec['layers'], depth, 2),
                    reference_eV_D=REFERENCE, baseline_L_um=100., eta_band=1., transfer_to_D=1.,
                    baseline_local_conversion_eV_D=e0, concentration_required_alone=required,
                    holding_energy_lifetime_required_s=tau_d*required,
                    attenuation_optimum_L_um=float(ls[imax]*1e6),
                    attenuation_max_eV_D=emax, attenuation_only_residual_factor=REFERENCE/emax,
                    empirical_quadratic_nonlinear_gain_max= gain,
                    nonlinear_scan_residual_factor=required/gain,
                    nonlinear_strength_alone='no extrapolated solution; finite scan gain only',
                    concentration_after_best_attenuation_and_scanned_nonlinearity=REFERENCE/emax/gain))
                for length in (1e-6, 10e-6, 30e-6, 100e-6, 1e-3, .1, 1000.):
                    kk = float(exp1(depth/length)-exp1(edge/length))
                    attenuation.append(dict(material=material, E_cm_keV=energy, depth_um=depth*1e6,
                        L_energy_m=length, geometric_kernel=kk,
                        local_conversion_eV_D=power/1e-4/length*kk*tau_d/nd/QE))
                # Pareto family: C * (tau/tau_D) = residual after best L and finite nonlinear gain.
                residual = REFERENCE/emax/gain
                for holdfactor in np.geomspace(1., residual, 9):
                    frontier.append(dict(material=material, E_cm_keV=energy, depth_um=depth*1e6,
                        attenuation_L_um=float(ls[imax]*1e6), scanned_nonlinear_gain=gain,
                        holding_factor=float(holdfactor), holding_energy_lifetime_s=float(tau_d*holdfactor),
                        concentration_factor=float(residual/holdfactor),
                        log_symmetric_cost_proxy=float(math.log10(holdfactor)**2+math.log10(residual/holdfactor)**2),
                        relation='C*(tau_hold/4.1ps)=reference/(attenuation_max*scanned_gain); eta=transfer=1'))
    return rows, frontier, attenuation


def event_requirements(out):
    """Consume saved observed durations; never replay the numerical producer."""
    peaks = list(csv.DictReader((out/'region_peaks.csv').open()))
    tails = list(csv.DictReader((out/'tails.csv').open()))
    models = list(csv.DictReader((out/'models.csv').open()))
    key = lambda r: (r['case'], r['arm'], r['side_mm'], r['z0_um'], r['z1_um'])
    peak_map = {key(r): r for r in peaks}
    tau_map = {(r['case'], r['arm']): float(r['energy_lifetime_s']) for r in models}
    rows, frontier = [], []
    for t in tails:
        if float(t['threshold_over_beam_mean']) != 4 or not t['duration_max_us']:
            continue
        p = peak_map[key(t)]
        duration = float(t['duration_max_us'])*1e-6
        conversion_rate = float(p['peak_J_m3'])/tau_map[t['case'], t['arm']]/float(p['n_D_m3'])/QE
        needed_seconds = REFERENCE/conversion_rate
        cap = conversion_rate*duration
        item = {k: p[k] for k in ('case', 'arm', 'side_mm', 'z0_um', 'z1_um')}
        rows.append(dict(**item, event_height_over_mean=4., longest_complete_event_s=duration,
            constant_peak_conversion_rate_eV_D_s=conversion_rate,
            hypothetical_sustained_peak_holding_required_s=needed_seconds,
            infinite_holding_event_energy_cap_eV_D=cap,
            concentration_required_even_at_infinite_holding=REFERENCE/cap,
            holding_alone_for_this_event='no finite holding solution' if cap < REFERENCE else 'conditional solution possible',
            envelope='constant observed peak throughout longest observed event; transfer=1; no background added'))
        for hold in np.geomspace(4.1e-12, 10*duration, 9):
            effective_hold = -hold*math.expm1(-duration/hold)
            frontier.append(dict(**item, event_duration_s=duration, holding_lifetime_s=float(hold),
                finite_event_effective_holding_s=effective_hold,
                concentration_factor=needed_seconds/effective_hold,
                infinite_holding_concentration_asymptote=needed_seconds/duration))
    table(out/'finite_event_requirements.csv', rows)
    table(out/'finite_event_joint_frontier.csv', frontier)
    print(json.dumps(dict(finite_event_requirements=len(rows), finite_event_frontier=len(frontier))))


def plate_bounds(stopping):
    rows = []
    mass = RHO*np.prod(DIMS)
    for energy, current in ((.25, 38e-6), (2.5, 138e-6)):
        rate = current/QE
        e = 2*energy*1e3*QE
        impulses = np.array([math.sqrt(2*2.01410*atomic_mass*e),
                             2*math.sqrt(2*2.01410*atomic_mass*e/2)])
        mean_f = rate*np.dot([.85, .15], impulses)
        psd = 2*rate*np.dot([.85, .15], impulses**2)
        for f in (1000., 10000., 100000.):
            for q in (1., 100., 10000.):
                omega = 2*np.pi*f
                # Unit-overlap translational effective mass = foil mass is a declared plate proxy.
                rows.append(dict(E_cm_keV=energy, foil_mass_kg=mass,
                    absorbed_force_N=mean_f, reflected_force_bound_N=2*mean_f,
                    one_sided_absorbed_force_shot_PSD_N2_Hz=psd,
                    hypothetical_supported_mode_Hz=f, assumed_Q=q,
                    effective_mass_kg=mass, overlap=1.,
                    mean_static_displacement_m=mean_f/(mass*omega**2),
                    shot_noise_mode_energy_J=psd*q/(4*mass*omega),
                    status='conditional unit-overlap supported mode; support/Q not measured'))
    return rows


def focused_checks(model, s, dt):
    omega = model['omega']
    a, _ = localized_strain(model)
    rng = np.random.default_rng(13)
    x = rng.normal(size=(1, len(omega))) / math.sqrt(len(omega))
    y = rng.normal(size=x.shape) / math.sqrt(len(omega))
    kicks = np.zeros((1600, 1, s.shape[1]))
    checks = dict(mass_orthonormal_max_error=model['mass_error'], eigen_relative_residual=model['residual'],
                  rigid_body_modes=model['rigid_count'])
    isolated = []
    for chi in (0., .3, 1.):
        h0 = float(.5*np.sum(x*x+y*y)+chi/4*(x@a)[0]**4)
        hist, xx, yy, energy, _ = evolve(model, s, kicks, dt, math.inf, chi, 1., (x, y))
        initial_modes = .5*(x*x+y*y)
        final_modes = .5*(xx*xx+yy*yy)
        isolated.append(dict(chi=chi, final_total_energy=float(energy[0]), initial_total_energy=h0,
            relative_energy_drift=float(abs(energy[0]-h0)/h0),
            mode_quadratic_energy_share_change_L1=float(np.sum(abs(final_modes/final_modes.sum()-initial_modes/initial_modes.sum()))),
            label='unforced random multimode conservation witness'))
    # Signed cancellation rejects energy-only accumulation at the actual evolution boundary.
    kick = np.zeros((3, 1, s.shape[1]))
    kick[0, 0, 0] = 1.
    pos = evolve(model, s, kick, dt, math.inf, 0., 1.)[0]
    neg = evolve(model, s, -kick, dt, math.inf, 0., 1.)[0]
    checks['opposite_signed_wavefield_max_residual'] = float(np.max(abs(pos+neg)))
    checks['max_unforced_relative_energy_drift'] = max(r['relative_energy_drift'] for r in isolated)
    if checks['mass_orthonormal_max_error'] > 1e-10 or checks['eigen_relative_residual'] > 1e-8:
        raise AssertionError('Variational operator failed its physical metric')
    if checks['max_unforced_relative_energy_drift'] > .01:
        raise AssertionError('Nonlinear energy drift exceeds declared 1 percent')
    return checks, isolated


def run(args):
    out = Path(args.output)
    if not out.is_absolute():
        out = ROOT/out
    out.mkdir(parents=True, exist_ok=True)
    profiles, stopping = profiles_and_stopping()
    dt, tau = args.dt, args.tau
    ntime, burn = round(args.exposure/dt), round(5*tau/dt)
    total = ntime+burn
    all_tails, union_tails, summary, input_rows, isolated_rows = [], [], [], [], []
    history_save = {}
    max_nonlinear_gain = 1.
    specs = [('isotropic', DIMS[0]/4), ('single_cubic', DIMS[0]/4),
             ('grains', DIMS[0]/8), ('grains', DIMS[0]/4)]
    for model_index, (kind, grain) in enumerate(specs):
        model = assemble(kind, grain, args.fmin, args.fmax)
        name = kind + ('_' + str(round(grain*1e3, 3)) + 'mm' if kind == 'grains' else '')
        omega = model['omega']
        # Same physical signed parent count realizations in every arm.
        representative = stopping['PdD', 2.5]
        depths = [deposition_depth(t) for t in representative['tracks']]
        s = source_matrix(model, depths)
        if model_index == 2:
            checks, conservation = focused_checks(model, s, dt)
            dump(out/'behavioral_checks.json', checks)
            table(out/'nonlinear_conservation.csv', conservation)
        nonlinear_levels = (0., .3, 1.) if model_index == 2 else (0.,)
        regions = []
        for side in (.0025, .005, .01, float(DIMS[0])):
            for z0, z1 in ((5e-6, 55e-6), (25e-6, 75e-6), (175e-6, 225e-6), (0., .00025)):
                ar = receiver(model, side, z0, z1)
                if side == DIMS[0] and z0 == 0:
                    defect = float(np.max(abs(ar*np.prod(DIMS)-np.eye(len(omega)))))
                    if defect > 1e-8:
                        raise AssertionError('Whole-foil receiver does not recover modal elastic energy')
                regions.append((side, z0, z1, ar))
        for material, energy, current in (('Pd', .25, 38e-6), ('PdD', .25, 38e-6),
                                          ('Pd', 2.5, 138e-6), ('PdD', 2.5, 138e-6)):
            st = stopping[material, energy]
            eh = np.array([st['tracks'][0]['components_eV'][1], 2*st['tracks'][1]['components_eV'][1]])*QE
            host_sum = np.dot([.85, .15], eh)/QE
            # Only central-energy species tracks were saved; accepted total includes
            # the five-node 5%-FWHM beam average. Preserve that mean exactly and
            # retain central-track relative marks as a declared spectral coarsening.
            mark_correction = st['host_recoil_J_per_true_C']/host_sum
            eh *= mark_correction
            actual_depths = [deposition_depth(t) for t in st['tracks']]
            s = source_matrix(model, actual_depths)
            power = current*st['host_recoil_J_per_true_C']*args.eta_band
            eunit = power*tau
            rng = np.random.default_rng(20261126)
            counts = []
            for fraction, eparent in zip((.85, .15), eh):
                mean = current/QE*fraction*dt/36/2
                signed = rng.poisson(mean, (total, args.realizations, 36))-rng.poisson(mean, (total, args.realizations, 36))
                counts.append(signed*math.sqrt(2*args.eta_band*eparent))
            kicks = np.concatenate(counts, axis=2)
            reference_eunit = 138e-6*stopping['PdD', 2.5]['host_recoil_J_per_true_C']*.01*5e-5
            for chi_reference in nonlinear_levels:
                chi = chi_reference*eunit/reference_eunit
                arm = name + '_chi' + str(chi_reference)
                hist, xx, yy, final_energy, snorm = evolve(model, s, kicks, dt, tau, chi, eunit)
                sampled = hist[burn:]
                case = material + '_' + str(energy) + 'keV'
                region_peaks = []
                union = None
                for side, z0, z1, ar in regions:
                    u = eunit/2*np.einsum('tri,ij,trj->tr', sampled, ar, sampled, optimize=True)
                    rec = profiles[material, energy]
                    temp = rec['case']['T_K']
                    # x thermal covariance = kBT/Eunit I, exact finite supported modal correlations.
                    thmean = KB*temp/2*np.trace(ar)
                    thsd = KB*temp/math.sqrt(2)*np.linalg.norm(ar)
                    all_tails += tail_rows(u, thmean, thsd, side, z0, z1, arm, case, dt)
                    standardized = u/u.mean()
                    union = standardized if union is None else np.maximum(union, standardized)
                    nd = column_density(rec['layers'], z0, z1)
                    converted_peak = float(u.max())/tau*4.1e-12/nd/QE
                    # Positive Hnl cannot beat this fixed-total-energy quadratic bound.
                    fixed_energy_bound_u = eunit*float(eigh(ar, eigvals_only=True)[-1])
                    fixed_energy_bound_d = fixed_energy_bound_u/tau*4.1e-12/nd/QE
                    region_peaks.append(dict(case=case, arm=arm, side_mm=side*1e3,
                        z0_um=z0*1e6, z1_um=z1*1e6, n_D_m3=nd,
                        peak_acoustic_quadratic_eV_per_D=float(u.max()/nd/QE),
                        reference_over_peak_acoustic=REFERENCE/(float(u.max())/nd/QE),
                        conditional_peak_conversion_eV_D=converted_peak,
                        conditional_transfer_fraction=1., conditional_holding_lifetime_s=4.1e-12,
                        concentration_required_from_modal_peak=REFERENCE/converted_peak,
                        holding_lifetime_required_from_modal_peak_s=4.1e-12*REFERENCE/converted_peak,
                        fixed_total_energy_budget_J=eunit,
                        any_positive_chi_quadratic_bound_J_m3=fixed_energy_bound_u,
                        any_positive_chi_conversion_bound_eV_D=fixed_energy_bound_d,
                        any_positive_chi_fixed_energy_residual_factor=REFERENCE/fixed_energy_bound_d,
                        mean_J_m3=float(u.mean()), peak_J_m3=float(u.max()),
                        thermal_mean_J_m3=float(thmean), thermal_sd_J_m3=float(thsd),
                        correlated_effective_quadratic_rank=float(np.trace(ar)**2/np.sum(ar*ar))))
                    if side == DIMS[0] and z0 == 0:
                        history_save[case+'_'+arm] = u.astype(np.float32)
                summary += region_peaks
                union_tails += union_tail_rows(union, arm, case, dt)
                if chi == 0:
                    linear_peaks = np.array([r['peak_J_m3'] for r in region_peaks])
                else:
                    max_nonlinear_gain = max(max_nonlinear_gain, float(np.max(np.array([r['peak_J_m3'] for r in region_peaks])/linear_peaks)))
                # Isolated packet shares source vector/norm, not the continuous source's nonlinear initial state.
                if model_index == 2 and energy == 2.5 and material == 'PdD':
                    packet = np.zeros((ntime, 1, 72))
                    packet_energy = args.eta_band*eh[0]
                    packet[0, 0, 14] = math.sqrt(2*packet_energy)
                    ih, _, _, ie, _ = evolve(model, s, packet, dt, math.inf, chi, eunit)
                    ar = regions[0][3]
                    iu = eunit/2*np.einsum('tri,ij,trj->tr', ih, ar, ih, optimize=True)
                    isolated_rows.append(dict(chi=chi, injected_packet_energy_J=packet_energy,
                        receiver_side_mm=2.5, z0_um=5., z1_um=55.,
                        peak_quadratic_J_m3=float(iu.max()), mean_quadratic_J_m3=float(iu.mean()),
                        final_total_energy_J=float(ie[0]), label='same isolated effective packet, unforced/unattenuated after launch'))
                input_rows.append(dict(case=case, arm=arm, retained_modes=len(omega),
                    min_retained_Hz=float(omega.min()/2/np.pi), max_retained_Hz=float(omega.max()/2/np.pi),
                    chart_true_current_uA=current*1e6, host_recoil_J_true_C=st['host_recoil_J_per_true_C'],
                    host_power_W=power/args.eta_band, eta_band=args.eta_band,
                    band_power_W=power, energy_lifetime_s=tau, assumed_energy_attenuation_length_m=4000*tau,
                    central_species_mark_mean_correction=mark_correction,
                    Dplus_host_deposition_depth_m=actual_depths[0],
                    D2_parent_host_deposition_depth_m=actual_depths[1],
                    band_energy_scale_J=eunit, chi_at_reference_PdD_2p5keV=chi_reference, chi=chi,
                    gamma_positive_quartic=chi/(C11*(.0025**2*50e-6)*eunit*snorm**4) if chi else 0.,
                    synthetic_columnar_grain_mm=grain*1e3 if kind == 'grains' else '',
                    energy_error_relative_eigen=model['residual'], rigid_modes_excluded=model['rigid_count']))
                input_rows[-1]['sampled_nonlinear_projection_variance'] = float(np.mean((sampled@localized_strain(model)[0])**2))
    table(out/'tails.csv', all_tails)
    table(out/'receiver_bank_union_tails.csv', union_tails)
    table(out/'region_peaks.csv', summary)
    table(out/'models.csv', input_rows)
    table(out/'isolated_packet.csv', isolated_rows)
    event_requirements(out)
    back, frontier, atten = backward(profiles, stopping, max_nonlinear_gain)
    table(out/'requirements.csv', back)
    table(out/'joint_frontier.csv', frontier)
    table(out/'attenuation_scan.csv', atten)
    table(out/'plate_bounds.csv', plate_bounds(stopping))
    np.savez_compressed(out/'whole_foil_histories.npz', **history_save)
    profile_rows = []
    for material in ('Pd', 'PdD'):
        for energy in (.25, 2.5):
            rec = profiles[material, energy]
            st = stopping[material, energy]
            for z in (.125e-9, 2e-9, 24e-9, 129e-9, 30e-6, 50e-6, 200e-6, 250e-6):
                profile_rows.append(dict(material=material, E_cm_keV=energy, depth_m=z,
                    n_Pd_m3=layer_value(rec['layers'], z, 2), n_D_m3=layer_value(rec['layers'], z),
                    x_D_per_Pd=layer_value(rec['layers'], z)/layer_value(rec['layers'], z, 2),
                    saved_profile_beam_time_s=rec['case']['beam_time_s'], T_K=rec['case']['T_K'],
                    host_recoil_J_true_C=st['host_recoil_J_per_true_C']))
    table(out/'fixed_profile.csv', profile_rows)
    ledger = dict(schema='g13-pass2-whole-foil-v1', base=BASE,
        sources={p: dict(bytes=(ROOT/p).stat().st_size, sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()) for p in SOURCES},
        elasticity=dict(C11_Pa=C11, C12_Pa=C12, C44_Pa=C44, density_kg_m3=RHO,
            citation='J. A. Rayne, Phys. Rev. 118, 1545 (1960), DOI 10.1103/PhysRev.118.1545',
            semantics='pure Pd near-300 K secondary transcription; S7 does not itself label temperature; not measured foil/PdDx tensor',
            transcription='Fukuda et al., Phys. Rev. Applied 20, 034048 (2023), supplemental Table S7; consistent with Malica/Dal Corso arXiv:2109.00459 Fig.1 near 300 K',
            density_semantics='Arblaster, DOI 10.1595/147106712X646113, assessed pure Pd at 293.15 K',
            alternative='Schwarz et al. DOI 10.1016/j.actamat.2004.10.009 and Bach et al. OSTI 977845; 300 K RUS PdDx measured 0.58<x<0.75; figure envelope at x~.6-.65: B165-180,Cprime25-26.5,C4463-65.5 GPa; no full-profile tensor assumed'),
        model=dict(dimensions_m=DIMS.tolist(), degree=[5,5,2], boundary='all six faces natural traction-free; rigid motions projected out',
            geometry='equal-area free square surrogate for 24-mm disk; central 10-mm square beam; mounting/shape not measured operator inputs',
            requested_band_Hz=[args.fmin,args.fmax], grain_maps='synthetic columnar piecewise rotated cubic; no measured grain map',
            time_step_s=dt, sampled_exposure_s=args.exposure, realizations=args.realizations,
            total_receiver_exposure_s=ntime*dt*args.realizations, burn_in_s=burn*dt,
            source='independent Poisson positive/negative parents in 36 equal-area sites; electrical fractions .85/.15; D2 two recoil energies share one parent sign',
            packet='effective z-traction point at host-deposition weighted mean depth; fixed retained-band energy per parent; terminal_depth_A is stored in METRES despite key',
            source_time='bin-boundary impulse approximation; no conditioning/window reweighting; no time compression',
            localization='positive quartic of center epsilon_xx at depth 30um couples retained foil modes; localized receiver sensitivity, not measured full-foil nonlinearity',
            thermal='classical same-band strain quadratic form; mean kBT Tr(A)/2; sd kBT sqrt(Tr(A²)/2); finite modal correlations retained',
            source_fraction=args.eta_band, damping_energy_lifetime_s=tau),
        result=dict(max_scanned_nonlinear_peak_ratio=max_nonlinear_gain,
            statistical_scope='fixed center receivers plus explicit 16-receiver union; union amplitudes normalized to its own mean; no extrapolated spatial event count',
            first_refinement='band-resolved launch and acoustic loss/transfer; nonlinear coefficient is not measured'),
        units=dict(energy_density='J/m3', reference='eV per D; contextual full-loading reference, not x<1 material threshold'))
    dump(out/'inputs.json', ledger)
    print(json.dumps(dict(outputs=str(out.relative_to(ROOT)), models=len(input_rows), tail_rows=len(all_tails),
                         max_nonlinear_peak_ratio=max_nonlinear_gain,
                         PdD_2p5keV_depth50um=[r for r in back if r['material']=='PdD' and r['E_cm_keV']==2.5 and r['depth_um']==50.][0]), indent=2))


def package(args):
    payloads, rows = {}, []
    for p in sorted(HERE.rglob('*')):
        if p.is_symlink():
            raise ValueError('No archive symlinks')
        if p.is_file() and '__pycache__' not in p.parts and p.name not in ('capsule.tar.gz', 'manifest.json', 'delivery.json'):
            data = p.read_bytes()
            name = str(p.relative_to(ROOT))
            payloads[name] = data
            rows.append(dict(path=name, bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
    manifest = dict(schema='g13-pass2-scoped-unpublished-v1', base=BASE, files=rows,
                    apply='new g13_rogue/pass2/ folder only; no parent/top-level replacements')
    raw = (json.dumps(manifest, indent=2)+'\n').encode()
    (HERE/'manifest.json').write_bytes(raw)
    payloads['g13_rogue/pass2/manifest.json'] = raw
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode='wb', filename='', mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode='w') as tf:
            for name, data in sorted(payloads.items()):
                ti = tarfile.TarInfo(name)
                ti.size, ti.mode, ti.mtime = len(data), 0o644, 0
                tf.addfile(ti, io.BytesIO(data))
    data = buf.getvalue()
    (HERE/'capsule.tar.gz').write_bytes(data)
    print(json.dumps(dict(path='g13_rogue/pass2/capsule.tar.gz', bytes=len(data), sha256=hashlib.sha256(data).hexdigest(), members=len(payloads))))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', default='g13_rogue/pass2/output')
    p.add_argument('--fmin', type=float, default=5000.)
    p.add_argument('--fmax', type=float, default=250000.)
    p.add_argument('--dt', type=float, default=5e-8)
    p.add_argument('--exposure', type=float, default=3e-4)
    p.add_argument('--tau', type=float, default=5e-5)
    p.add_argument('--realizations', type=int, default=6)
    p.add_argument('--eta-band', type=float, default=.01)
    p.add_argument('--package', action='store_true')
    p.add_argument('--live', action='store_true',
                   help='run corrected conditional-source, signed-hole and structured-beam sweep')
    p.add_argument('--live-impact-heights', action='store_true',
                   help='aggregate missing predictable impact heights, reusing frozen convergence/tails')
    p.add_argument('--event-requirements', action='store_true',
                   help='consume saved peak/duration tables without replaying propagation')
    a = p.parse_args()
    if a.package:
        package(a)
    elif a.event_requirements:
        out = Path(a.output)
        event_requirements(out if out.is_absolute() else ROOT/out)
    elif a.live_impact_heights:
        import live
        live.height_readout(a)
    elif a.live:
        import live
        live.run(a)
    else:
        if not (0 < a.eta_band <= 1 and a.dt > 0 and a.tau > 0 and a.exposure > 0 and a.realizations > 0):
            p.error('Positive physical settings and 0<eta-band<=1 required')
        run(a)


if __name__ == '__main__':
    main()
