"""Entry_4103/4104: predictable live launches, directional convergence, rogue holes.

The executable entrypoint is coarse.py --live. No damage or dose model is used.
"""
from pathlib import Path
import json
import math
import numpy as np
from scipy.linalg import solve_continuous_lyapunov
from scipy.constants import elementary_charge as QE, Boltzmann as KB
import coarse as c


SITE_SIDE = 6
DT = 5e-8
EXPOSURE = 1e-4
REPS = 2
ETA = .01
TAU_REFERENCE = 50e-6
SPEED_RANGE = 4000.


def sites():
    axis = (c.DIMS[0]-.01)/2+(.5+np.arange(SITE_SIDE))*.01/SITE_SIDE
    xy = np.stack(np.meshgrid(axis, axis, indexing='ij'), axis=-1).reshape(-1, 2)
    group = (xy[:, 0] >= c.DIMS[0]/2).astype(int)*2+(xy[:, 1] >= c.DIMS[1]/2)
    return xy, group

def beam_profile(level):
    xy, _ = sites()
    p = xy-c.DIMS[:2]/2
    if level == 0:
        return np.full(len(xy), 1/len(xy))
    weight = np.exp(-np.sum(p*p, axis=1)/(2*.003**2))
    if level == 2:
        for center, gain in ((np.array([.0025, -.0025]), 5.), (np.array([-.0025, .0025]), 3.)):
            weight += gain*np.exp(-np.sum((p-center)**2, axis=1)/(2*.0012**2))
    return weight/weight.sum()


def beam_timing(steps, level, rng):
    if level == 0:
        return np.ones((steps, REPS))
    phase = rng.uniform(0., 2*np.pi, REPS)
    result = 1+.7*np.cos(2*np.pi*100000.*np.arange(steps)[:, None]*DT+phase)
    if level == 2:
        # Random Markov intensity: ON 2us, OFF 18us in mean, plus current ripple.
        on = rng.random(REPS) < .1
        for t in range(steps):
            result[t] += 9*on
            switch = rng.random(REPS) < np.where(on, -np.expm1(-DT/2e-6), -np.expm1(-DT/18e-6))
            on ^= switch
    return result


def mode_group_velocities(model):
    bulk = (c.C11+2*c.C12)/3
    gv = (c.C11-c.C12+3*c.C44)/5
    gr = 5*(c.C11-c.C12)*c.C44/(4*c.C44+3*(c.C11-c.C12))
    shear = (gv+gr)/2
    young = 9*bulk*shear/(3*bulk+shear)
    nu = (3*bulk-2*shear)/(2*(3*bulk+shear))
    bending = young*c.DIMS[2]**3/(12*(1-nu*nu))
    alpha = math.sqrt(bending/(c.RHO*c.DIMS[2]))
    phi, b = c.basis(model['points'])
    n = phi.shape[1]
    vz = model['v'][2*n:]
    vertical_mass = np.sum(vz*vz, axis=0)
    zz = phi@vz
    dx = b[:, 4, 2*n:]@vz
    dy = b[:, 3, 2*n:]@vz
    kk = np.sum(model['w'][:, None]*(dx*dx+dy*dy), axis=0)/np.maximum(np.sum(model['w'][:, None]*zz*zz, axis=0), np.finfo(float).tiny)
    k = np.sqrt(kk)
    flexural = (vertical_mass > .5) & (k*c.DIMS[2] < .5)
    divergence = (b[:, 0]+b[:, 1]+b[:, 2])@model['v']
    bulk_fraction = bulk*np.sum(model['w'][:, None]*divergence**2, axis=0)/model['omega']**2
    speed = np.where(bulk_fraction > .5, math.sqrt((bulk+4*shear/3)/c.RHO), math.sqrt(shear/c.RHO))
    speed[flexural] = 2*alpha*k[flexural]
    model['group_velocity'] = speed
    return dict(bulk_modulus_Pa=bulk, shear_modulus_Pa=shear, Young_Pa=young,
                poisson=nu, bending_N_m=bending, alpha_m2_s=alpha,
                bulk_L_m_s=math.sqrt((bulk+4*shear/3)/c.RHO),
                bulk_T_m_s=math.sqrt(shear/c.RHO), flexural_modes=int(flexural.sum()),
                rows=[dict(mode=i, frequency_Hz=float(model['omega'][i]/2/np.pi),
                           flexural=bool(flexural[i]), group_speed_proxy_m_s=float(speed[i]),
                           effective_wavelength_m=float(2*np.pi/k[i]) if k[i] else '',
                           bulk_energy_fraction=float(bulk_fraction[i]), kh=float(k[i]*c.DIMS[2]))
                      for i in range(len(k))])

def timeframes(metadata, out):
    from scipy.integrate import quad
    waves = [('Hill_bulk_L', metadata['bulk_L_m_s'], metadata['bulk_L_m_s']/100000., 100000., False),
             ('Hill_bulk_T', metadata['bulk_T_m_s'], metadata['bulk_T_m_s']/100000., 100000., False)]
    for wavelength in (.005, .01, .015, .02):
        k = 2*np.pi/wavelength
        frequency = metadata['alpha_m2_s']*k*k/(2*np.pi)
        waves.append((f'flexural_lambda{wavelength*1000:g}mm', 2*metadata['alpha_m2_s']*k, wavelength, frequency, True))
    paths = [('foil_thickness', .00025), ('spot_equivalent_disk_diameter', 2*math.sqrt(1e-4/math.pi)),
             ('actual_disk_longest_span', .024), ('actual_center_to_edge', .012),
             ('actual_edge_and_return', .024), ('model_square_side', float(c.DIMS[0])),
             ('model_center_to_edge', float(c.DIMS[0]/2)), ('model_edge_and_return', float(c.DIMS[0]))]
    rows, survival = [], []
    xy, _ = sites()
    radius_from_center = np.linalg.norm(xy-c.DIMS[:2]/2, axis=1)
    for name, speed, wavelength, frequency, guided in waves:
        area = 4*quad(lambda x: min(.005, math.sqrt(max(0., wavelength*wavelength-x*x))), 0., min(.005,wavelength), epsabs=1e-14)[0]
        area = min(1e-4, area)
        distances = np.linalg.norm(xy[:, None]-xy[None, :], axis=2)
        max_structured_fraction = float(np.max((distances <= wavelength)@beam_profile(2)))
        for path, distance in paths:
            valid = not (guided and path == 'foil_thickness')
            seconds = distance/speed if valid else None
            for energy, current in ((.25,38e-6),(2.5,138e-6)):
                rows.append(dict(wave=name, E_cm_keV=energy, frequency_Hz=frequency,
                    in_selected_5_to_250kHz_band=5000 <= frequency <= 250000,
                    group_speed_m_s=speed, wavelength_m=wavelength, kh=2*np.pi*.00025/wavelength,
                    path=path, path_m=distance, time_s=seconds if valid else '',
                    whole_spot_mean_parent_impacts=current/QE*seconds if valid else '',
                    wavelength_radius_neighborhood_area_m2=area,
                    uniform_mean_parents_within_wavelength_radius=current/QE*seconds*area/1e-4 if valid else '',
                    structured_center_discrete_profile_fraction=float(beam_profile(2)[radius_from_center <= wavelength].sum()),
                    structured_center_mean_parent_impacts=current/QE*seconds*float(beam_profile(2)[radius_from_center <= wavelength].sum()) if valid else '',
                    structured_peak_neighborhood_mean_parent_impacts=current/QE*seconds*max_structured_fraction if valid else '',
                    structured_peak_neighborhood_profile_fraction=max_structured_fraction,
                    structured_peak_cell_current_density_over_uniform=float(36*beam_profile(2).max()),
                    convention='guided flexure has no independent thickness flight' if not valid else 'group travel; parent charge arrivals, not coherent amplitudes'))
        for length in (.00025,.002,.02,.2,2.):
            for tau in (2e-6,50e-6,200e-6):
                for reflection in (0.,.1,1.):
                    edge = math.exp(-.012/length-.012/speed/tau)
                    returned = reflection*math.exp(-.024/length-.024/speed/tau)
                    survival.append(dict(wave=name, frequency_Hz=frequency, L_energy_m=length,
                        intrinsic_energy_lifetime_s=tau, assumed_flange_energy_reflectivity=reflection,
                        energy_at_actual_edge=edge, energy_at_actual_return=returned,
                        amplitude_at_actual_return=math.sqrt(returned),
                        reaches_edge_with_at_least_10pct_energy=edge >= .1,
                        returns_with_at_least_10pct_energy=returned >= .1,
                        convention='exp(-path/L_E-time/tau_E), one R_E at edge; unmeasured flange; actual disk geometry proxy'))
    c.table(out/'timeframes.csv', rows)
    c.table(out/'edge_survival.csv', survival)
    c.dump(out/'beam_parameters.json', dict(spatial_sites_m=xy.tolist(),
        spatial_weights={str(k):beam_profile(k).tolist() for k in range(3)},
        spatial_hypotheses='0 uniform; 1 Gaussian sigma3mm; 2 adds sigma1.2mm hot spots at (+2.5,-2.5)mm and (-2.5,+2.5)mm, relative gains5 and3',
        temporal_hypotheses='0 steady; 1 amplitude.7 100kHz ripple with random record phase; 2 adds Markov ON/OFF burst gain9, mean ON2us/OFF18us',
        charge='Poisson parent totals drawn once per record and burn interval before profile allocation; same totals across matched space/time arms; random intensity path normalized within each interval',
        measurements='none of the profile/ripple/burst parameters is a measured beam map; chart mean current alone is inherited',
        geometry='timeframes use actual 24mm disk proxy and equal-area square separately; structured neighborhood counts are discrete source-cell quadrature, uniform neighborhood is clipped circle-square area',
        plate={k:v for k,v in metadata.items() if k != 'rows'}))




def probes(model):
    xy, groups = sites()
    projection = [linear_region(model, .0025, 0., 50e-6, center=p) for p in xy]
    return np.array([p[0] for p in projection]), np.array([p[1] for p in projection]), groups


def source_reference(model, s, weights, eps, uz, power):
    """Harmonic same-band covariance under the named quiet launch distribution."""
    omega = model['omega']
    m = len(omega)
    a = np.zeros((2*m, 2*m))
    a[:m, :m] = -np.eye(m)/(2*TAU_REFERENCE)
    a[m:, m:] = a[:m, :m]
    a[:m, m:] = np.diag(omega)
    a[m:, :m] = -np.diag(omega)
    diffusion = np.zeros_like(a)
    diffusion[m:, m:] = 2/TAU_REFERENCE*((s*weights) @ s.T)
    covariance = solve_continuous_lyapunov(a, -diffusion)
    residual = np.linalg.norm(a@covariance+covariance@a.T+diffusion)/np.linalg.norm(diffusion)
    ecanon = power*TAU_REFERENCE
    xx = covariance[:m, :m]
    es = np.sqrt(np.maximum(np.einsum('si,ij,sj->s', eps, xx, eps), 0)*ecanon)
    us = np.sqrt(np.maximum(np.einsum('si,ij,sj->s', uz, xx, uz), 0)*ecanon)
    return covariance, ecanon, es, us, residual


def linear_region(model, side, z0, z1, center=None):
    x, w = np.polynomial.legendre.leggauss(3)
    center = c.DIMS[:2]/2 if center is None else center
    bounds = ((center[0]-side/2, center[0]+side/2),
              (center[1]-side/2, center[1]+side/2), (z0, z1))
    pts = np.stack(np.meshgrid(*[(lo+hi)/2+(hi-lo)*x/2 for lo, hi in bounds], indexing='ij'), axis=-1).reshape(-1, 3)
    weights = np.einsum('i,j,k->ijk', *[w*(hi-lo)/2 for lo, hi in bounds]).ravel()
    phi, b = c.basis(pts)
    n = phi.shape[1]
    vol = side**2*(z1-z0)
    eps = (weights @ b[:, 0])/vol @ model['v']/model['omega']
    uz = (weights @ phi)/vol @ model['v'][2*n:]/model['omega']
    return eps, uz


def off_beam_receiver(model):
    """Finite-band inward energy flux at six actual receiver faces."""
    center = np.array([c.DIMS[0]/2+.0065, c.DIMS[1]/2, .000125])
    half = np.array([.0005, .0005, 50e-6])
    normals = np.r_[np.eye(3), -np.eye(3)]
    points = center+normals*half
    phi, b = c.basis(points)
    n = phi.shape[1]
    bm = b @ model['v']/model['omega']
    stress = np.einsum('qij,qjk->qik', c.tensors(points, model['grain'], model['kind']), bm)
    velocity = np.zeros((6, 3, len(model['omega'])))
    for d in range(3):
        velocity[:, d] = phi @ model['v'][d*n:(d+1)*n]
    _, bc = c.basis(center[None])
    cc = c.tensors(center[None], model['grain'], model['kind'])[0]
    bb = bc[0] @ model['v']/model['omega']
    energy = bb.T @ cc @ bb
    areas = np.array([4*half[1]*half[2], 4*half[0]*half[2], 4*half[0]*half[1]]*2)
    return stress, velocity, normals, areas, energy, center.tolist(), (2*half).tolist()


def convergence(x, y, receiver, ecanon):
    stress_map, velocity_map, normals, areas, energy, _, _ = receiver
    # r=record, g=disjoint parent-origin sector, f=receiver face.
    stress6 = math.sqrt(ecanon)*np.einsum('fij,rgj->rgfi', stress_map, x)
    velocity = math.sqrt(ecanon)*np.einsum('fij,rgj->rgfi', velocity_map, y)
    sig = np.empty((*stress6.shape[:-1], 3, 3))
    sig[..., 0, 0], sig[..., 1, 1], sig[..., 2, 2] = stress6[..., 0], stress6[..., 1], stress6[..., 2]
    sig[..., 1, 2] = sig[..., 2, 1] = stress6[..., 3]
    sig[..., 0, 2] = sig[..., 2, 0] = stress6[..., 4]
    sig[..., 0, 1] = sig[..., 1, 0] = stress6[..., 5]
    flux = -np.einsum('rgfij,rgfj->rgfi', sig, velocity)
    inward = np.maximum(0., -np.einsum('rgfi,fi->rgf', flux, normals))*areas
    floor = .05*inward.sum(axis=(1, 2))
    active = inward > floor[:, None, None]
    directional = np.zeros(len(x), dtype=bool)
    face_separation = normals@normals.T <= 0
    for g in range(4):
        for h in range(g+1, 4):
            directional |= np.any(active[:, g, :, None] & active[:, h, None, :] & face_separation, axis=(1, 2))
    total = x.sum(axis=1)
    coherent = ecanon/2*np.einsum('ri,ij,rj->r', total, energy, total)
    separate = ecanon/2*np.einsum('rgi,ij,rgj->r', x, energy, x)
    constructive = coherent > 1.1*separate
    return directional & constructive, directional, constructive, coherent, separate


def source_masks(eps, uz, eps_sigma, uz_sigma, canonical_eps, canonical_uz):
    # Every predicate consumes only the wave field that existed before this bin.
    a = (np.abs(eps) > .25*canonical_eps) | (np.abs(uz) > .25*canonical_uz)
    ep = eps >= 3*eps_sigma
    en = eps <= -3*eps_sigma
    up = uz >= 3*uz_sigma
    un = uz <= -3*uz_sigma
    hole = ep | en | up | un
    return dict(all=np.ones_like(a), A_disturbed=a, C_strain_positive=ep,
                C_strain_negative=en, C_displacement_positive=up,
                C_displacement_negative=un, C_any=hole,
                A_and_C=a & hole, A_or_C=a | hole, neither=~(a | hole))


def occurrence(mask, dt):
    runs, censored = [], 0
    for r in range(mask.shape[1]):
        d = np.diff(np.r_[False, mask[:, r], False].astype(int))
        for lo, hi in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
            if lo == 0 or hi == len(mask):
                censored += 1
            else:
                runs.append((hi-lo)*dt)
    exposure = mask.size*dt
    return dict(time_fraction=float(mask.mean()), complete_events=len(runs),
        boundary_censored_events=censored, receiver_seconds=exposure,
        initiation_rate_s=len(runs)/exposure,
        zero_event_95pct_Poisson_rate_model_upper_s=(-math.log(.05)/exposure if not runs else ''),
        max_duration_us=(max(runs)*1e6 if runs else ''),
        median_duration_us=(float(np.median(runs))*1e6 if runs else ''),
        duration_ge_1us_events=sum(t >= 1e-6 for t in runs),
        duration_ge_10us_events=sum(t >= 10e-6 for t in runs),
        duration_ge_50us_events=sum(t >= 50e-6 for t in runs))


def signed_tails(values, standardized, thermal_sigma, observable, side, z0, z1, case, arm):
    result = []
    for sign, label in ((1., 'positive_peak'), (-1., 'negative_hole')):
        for height in (1., 3., 5., 10.):
            mask = sign*standardized >= height
            selected = sign*values
            result.append(dict(case=case, arm=arm, observable=observable, sign=label,
                side_mm=side*1e3, z0_um=z0*1e6, z1_um=z1*1e6,
                threshold_over_causal_prior_rms=height, ensemble_mean_reference=0.,
                physical_unit='strain' if observable == 'epsilon_xx' else 'm',
                physical_extreme=float(selected.max()), empirical_rms=float(np.sqrt(np.mean(values**2))),
                mean_signed_value=float(values.mean()), thermal_sigma=thermal_sigma,
                extreme_over_thermal_sigma=float(selected.max()/thermal_sigma),
                standardized_extreme=float((sign*standardized).max()), **occurrence(mask, DT)))
    return result


def launch_statistics(counts, class_mask, eh, fractions, cap, accumulator, state):
    # Counts axes record, site, species, sign, mark. Count-weighted parent statistics.
    ev, uv, ez, uz = state
    values = np.stack([ev, uv, abs(ev), abs(uv), abs(ez), abs(uz)], axis=-1)
    for name, mask in class_mask.items():
        selected = counts*mask[:, :, None, None, None]
        n = float(selected.sum())
        energies = eh[None, None, :, None, None]*fractions[None, None, None, None, :]
        actual = energies*cap[:, None, None, None, None]**2
        v = accumulator[name]
        v['parents'] += n
        v['launch_self_J'] += float(np.sum(selected*actual))
        v['marks'] += selected.sum(axis=(0, 1, 2, 3))
        weights = selected.sum(axis=(2, 3, 4))
        v['parents_by_record'] += weights.sum(axis=1)
        v['state_sums'] += np.sum(weights[:, :, None]*values, axis=(0, 1))
        if n:
            v['max_launch_self_J'] = max(v['max_launch_self_J'], float(np.max(np.where(selected > 0, actual, 0))))
            valid = weights > 0
            v['state_min'] = np.minimum(v['state_min'], np.array([np.where(valid, ev, np.inf).min(), np.where(valid, uv, np.inf).min()]))
            v['state_max'] = np.maximum(v['state_max'], np.array([np.where(valid, ev, -np.inf).max(), np.where(valid, uv, -np.inf).max()]))
            v['standardized_abs_max'] = np.maximum(v['standardized_abs_max'], np.array([np.where(valid, abs(ez), 0.).max(), np.where(valid, abs(uz), 0.).max()]))


def impact_heights(a, all_by_record):
    result = {}
    for i, name in enumerate(('epsilon_xx','u_z','abs_epsilon_xx','abs_u_z','abs_epsilon_over_prior_rms','abs_u_over_prior_rms')):
        result['preimpact_mean_'+name] = float(a['state_sums'][i]/a['parents']) if a['parents'] else ''
    for i, name in enumerate(('epsilon_xx','u_z')):
        result['preimpact_min_'+name] = float(a['state_min'][i]) if a['parents'] else ''
        result['preimpact_max_'+name] = float(a['state_max'][i]) if a['parents'] else ''
        result['preimpact_max_abs_'+name+'_over_prior_rms'] = float(a['standardized_abs_max'][i]) if a['parents'] else ''
    fractions = a['parents_by_record']/all_by_record
    result['fraction_across_independent_records_min'] = float(fractions.min())
    result['fraction_across_independent_records_max'] = float(fractions.max())
    return result


def simulate(model, source, group, eps, uz, reference, receiver, regions, energy, current, tau, length, feedback, spatial, timing, edge_reflectivity, frozen_convergence=None):
    covariance, ecanon, canonical_eps, canonical_uz, _ = reference
    omega = model['omega']
    loss_rate = 1/tau+model['group_velocity']/length-math.log(edge_reflectivity)*model['group_velocity']/c.DIMS[0]
    modal_lifetimes = 1/loss_rate
    tau_eff = float(np.median(modal_lifetimes))
    tau_memory = float(modal_lifetimes.max())
    burn = round(5*tau_memory/DT)
    sample_steps = round(EXPOSURE/DT)
    steps = burn+sample_steps
    co, si = np.cos(omega*DT), np.sin(omega*DT)
    damp = np.exp(-DT/(2*modal_lifetimes))
    eh = source['eh']
    # Exactly mean ETA under p_high=.01 and p_low=.25.
    fractions = np.minimum(1., ETA/.865*np.array([.1, 1., 10.]))
    count_rng = np.random.default_rng(4104)
    mark_rng = np.random.default_rng(8104)
    beam_rng = np.random.default_rng(14109)
    site_weights = beam_profile(spatial)
    time_weights = beam_timing(steps, timing, beam_rng)
    # Shared Poisson record totals, then a normalized random Cox intensity path.
    time_counts = np.empty((steps, REPS), dtype=np.int64)
    total_rng = np.random.default_rng(24109)
    intervals = (slice(0, burn), slice(burn, steps))
    shared_totals = np.array([[total_rng.poisson(current/QE*(part.stop-part.start)*DT) for rep in range(REPS)] for part in intervals], dtype=np.int64)
    for k, part in enumerate(intervals):
        for rep in range(REPS):
            weights_t = time_weights[part, rep]
            time_counts[part, rep] = count_rng.multinomial(int(shared_totals[k, rep]), weights_t/weights_t.sum())
    x = np.zeros((REPS, 4, len(omega)))
    y = np.zeros_like(x)
    variance_eps = np.broadcast_to(canonical_eps**2, (REPS, len(group))).copy()
    variance_uz = np.broadcast_to(canonical_uz**2, variance_eps.shape).copy()
    region_sigma = []
    for r in regions:
        for proj in r['projectors']:
            region_sigma.append(np.full(REPS, ecanon*float(proj@covariance[:len(omega), :len(omega)]@proj)))
    region_sigma = np.array(region_sigma).reshape((-1, REPS))
    projectors = np.array([p for region in regions for p in region['projectors']]).reshape((-1, len(omega)))
    alpha_rms = -math.expm1(-DT/(5*tau_memory))
    sampled_x = np.empty((sample_steps, REPS, len(omega)))
    sampled_signed = np.empty((sample_steps, REPS, len(regions)*2))
    sampled_z = np.empty_like(sampled_signed)
    bmask = np.empty((sample_steps, REPS), dtype=bool) if frozen_convergence is None else np.array(frozen_convergence, copy=True)
    direction_mask = np.empty_like(bmask) if frozen_convergence is None else None
    construct_mask = np.empty_like(bmask) if frozen_convergence is None else None
    off_u = np.empty((sample_steps, REPS)) if frozen_convergence is None else None
    off_separate = np.empty_like(off_u) if frozen_convergence is None else None
    labels = ('all', 'A_disturbed', 'C_strain_positive', 'C_strain_negative',
              'C_displacement_positive', 'C_displacement_negative', 'C_any', 'A_and_C', 'A_or_C', 'neither',
              'A_and_B_remote_window', 'C_and_B_remote_window')
    accumulator = {k: dict(parents=0., launch_self_J=0., max_launch_self_J=0., marks=np.zeros(3),
                           parents_by_record=np.zeros(REPS), state_sums=np.zeros(6),
                           state_min=np.full(2, np.inf), state_max=np.full(2, -np.inf),
                           standardized_abs_max=np.zeros(2)) for k in labels}
    source_work = incoming_host = marked_self = damping_loss = 0.
    initial_energy = 0.
    caps = 0
    max_budget_violation = 0.
    parent_bins = np.empty((REPS, len(group), 2, 2), dtype=np.int64)
    for t in range(steps):
        pre = x.sum(axis=1)
        ev = math.sqrt(ecanon)*(pre @ eps.T)
        uv = math.sqrt(ecanon)*(pre @ uz.T)
        eps_sd = np.sqrt(np.maximum(variance_eps, np.finfo(float).tiny))
        uz_sd = np.sqrt(np.maximum(variance_uz, np.finfo(float).tiny))
        masks = source_masks(ev, uv, eps_sd, uz_sd, canonical_eps, canonical_uz)
        prior_region_values = math.sqrt(ecanon)*(pre@projectors.T)
        prior_region_z = prior_region_values/np.sqrt(np.maximum(region_sigma.T, np.finfo(float).tiny))
        region_sigma = (1-alpha_rms)*region_sigma+alpha_rms*prior_region_values.T**2
        # B has its own receiving-window denominator, not the parent arrival denominator.
        if t >= burn:
            j = t-burn
            if t == burn:
                initial_energy = ecanon/2*float(np.sum(x.sum(axis=1)**2+y.sum(axis=1)**2))
            sampled_x[j] = pre
            sampled_signed[j] = prior_region_values
            sampled_z[j] = prior_region_z
            if frozen_convergence is None:
                bb, dd, kk, uu, ss = convergence(x, y, receiver, ecanon)
                bmask[j], direction_mask[j], construct_mask[j] = bb, dd, kk
                off_u[j], off_separate[j] = uu, ss
            else:
                bb = bmask[j]
            masks['A_and_B_remote_window'] = masks['A_disturbed'] & bb[:, None]
            masks['C_and_B_remote_window'] = masks['C_any'] & bb[:, None]
        variance_eps = (1-alpha_rms)*variance_eps+alpha_rms*ev*ev
        variance_uz = (1-alpha_rms)*variance_uz+alpha_rms*uv*uv
        for rep in range(REPS):
            per_site = count_rng.multinomial(int(time_counts[t, rep]), site_weights)
            dplus = count_rng.binomial(per_site, .85)
            per_species = np.column_stack([dplus, per_site-dplus])
            positive = count_rng.binomial(per_species, .5)
            parent_bins[rep, :, :, 0] = positive
            parent_bins[rep, :, :, 1] = per_species-positive
        signal = .5*(np.tanh(ev/canonical_eps)**2+np.tanh(uv/canonical_uz)**2)
        high_probability = .01+feedback*signal
        # Conditional multinomial marks; the actual spread, not its mean, drives the field.
        low = mark_rng.binomial(parent_bins, .25)
        rest = parent_bins-low
        high = mark_rng.binomial(rest, high_probability[:, :, None, None]/.75)
        mid = rest-high
        counts = np.stack([low, mid, high], axis=-1)
        signed = counts[:, :, :, 0, :]-counts[:, :, :, 1, :]
        amplitudes = np.sum(signed*np.sqrt(2*eh[None, None, :, None]*fractions), axis=-1)
        delta = np.zeros_like(y)
        for g in range(4):
            members = group == g
            for species in range(2):
                delta[:, g] += amplitudes[:, members, species] @ source['s'][species][:, members].T
        old_p = math.sqrt(ecanon)*y.sum(axis=1)
        dp = delta.sum(axis=1)
        aa = .5*np.sum(dp*dp, axis=1)
        bb = np.sum(old_p*dp, axis=1)
        host = np.sum(parent_bins*eh[None, None, :, None], axis=(1, 2, 3))
        # Global coherent cross work is accounted once; new work cannot exceed all incoming host energy.
        cap = np.ones(REPS)
        over = aa+bb > host
        disc = np.sqrt(bb*bb+4*aa*host)
        cap[over] = 2*host[over]/(bb[over]+disc[over])
        actual_work = aa*cap*cap+bb*cap
        max_budget_violation = max(max_budget_violation, float(np.max(actual_work-host)))
        caps += int(over.sum())
        y += delta*cap[:, None, None]/math.sqrt(ecanon)
        if t >= burn:
            launch_statistics(counts, masks, eh, fractions, cap, accumulator, (ev, uv, ev/eps_sd, uv/uz_sd))
            incoming_host += float(host.sum())
            source_work += float(actual_work.sum())
            marked_self += float(np.sum(counts*eh[None, None, :, None, None]*fractions*cap[:, None, None, None, None]**2))
        # Exact unforced finite-band harmonic propagation and explicit energy damping.
        old_modal_energy = ecanon/2*(x.sum(axis=1)**2+y.sum(axis=1)**2)
        x, y = x*co+y*si, y*co-x*si
        x *= damp
        y *= damp
        if t >= burn:
            damping_loss += float(np.sum(old_modal_energy*(1-damp*damp)))
    return dict(x=sampled_x, signed=sampled_signed, z=sampled_z, B=bmask,
                directional=direction_mask, constructive=construct_mask,
                off_u=off_u, off_separate=off_separate, impact=accumulator,
                source_work_J=source_work, incoming_host_J=incoming_host,
                marked_self_J=marked_self, damping_loss_J=damping_loss,
                capped_record_bins=caps, max_budget_violation_J=max_budget_violation,
                energy_balance_residual_J=(ecanon/2*float(np.sum(x.sum(axis=1)**2+y.sum(axis=1)**2))-initial_energy-source_work+damping_loss),
                mode_lifetime_min_s=float(modal_lifetimes.min()), mode_lifetime_max_s=tau_memory,
                spatial_profile_peak_over_mean=float(site_weights.max()*len(group)),
                spatial_profile_cv=float(np.std(site_weights)/np.mean(site_weights)),
                tau_eff=tau_eff, ecanon=ecanon)


def height_readout(args):
    """Replay only predictable impact states; reuse frozen convergence/field tails."""
    import csv
    global DT, EXPOSURE, REPS, ETA
    output = Path(args.output)
    out = (output if output.is_absolute() else c.ROOT/output)/'live'
    config = json.loads((out/'contract.json').read_text())['sample']
    DT, EXPOSURE, REPS, ETA = config['dt_s'], config['exposure_s'], config['independent_records'], args.eta_band
    read_rows = lambda name: list(csv.DictReader((out/name).open()))
    rows = read_rows('impact_classes.csv')
    models = read_rows('models.csv')
    old = {(r['case'],r['arm'],r['event_class']):r for r in rows}
    profiles, stopping = c.profiles_and_stopping()
    model = c.assemble('grains', c.DIMS[0]/8, args.fmin, args.fmax)
    mode_group_velocities(model)
    eps, uz, groups = probes(model)
    st = stopping['PdD',2.5]
    sr = c.source_matrix(model,[c.deposition_depth(t) for t in st['tracks']])
    energy_weights = np.array([.85*st['tracks'][0]['components_eV'][1],.3*st['tracks'][1]['components_eV'][1]])
    energy_weights /= energy_weights.sum()
    weights = np.repeat(energy_weights/36,36)
    quiet = float(np.dot([.25,.74,.01],np.minimum(1.,ETA/.865*np.array([.1,1.,10.]))))
    reference = source_reference(model,sr,weights,eps,uz,quiet*138e-6*st['host_recoil_J_per_true_C'])
    sources = {}
    for energy in (.25,2.5):
        st = stopping['PdD',energy]
        eh = np.array([st['tracks'][0]['components_eV'][1],2*st['tracks'][1]['components_eV'][1]])*QE
        eh *= st['host_recoil_J_per_true_C']/(np.dot([.85,.15],eh)/QE)
        sm = c.source_matrix(model,[c.deposition_depth(t) for t in st['tracks']])
        sources[energy] = dict(eh=eh,s=(sm[:,:36],sm[:,36:]))
    maximum_work_difference = 0.
    with np.load(out/'selected_histories.npz') as frozen:
        for m in models:
            energy = .25 if m['case']=='PdD_0.25keV' else 2.5
            current = 38e-6 if energy==.25 else 138e-6
            r = simulate(model,sources[energy],groups,eps,uz,reference,None,[],energy,current,
                float(m['intrinsic_lifetime_s']),float(m['range_proxy_m']),float(m['feedback']),
                int(m['spatial_clumping']),int(m['temporal_clumping']),float(m['lateral_edge_energy_reflectivity']),
                frozen_convergence=frozen[m['case']+'_'+m['arm']+'_B'])
            difference = abs(r['source_work_J']-float(m['source_work_J']))
            maximum_work_difference = max(maximum_work_difference,difference)
            if difference > 1e-12*float(m['incoming_host_J']):
                raise AssertionError('Height consumer changed its saved source trajectory')
            for label,a in r['impact'].items():
                target = old[m['case'],m['arm'],label]
                if float(target['counted_parents']) != a['parents']:
                    raise AssertionError('Height consumer selected different parent impacts')
                target.update(impact_heights(a,r['impact']['all']['parents_by_record']))
                target['mean_launch_self_eV_parent'] = a['launch_self_J']/a['parents']/QE if a['parents'] else ''
                target['max_supported_launch_self_eV_parent'] = a['max_launch_self_J']/QE if a['parents'] else ''
    c.table(out/'impact_classes.csv',rows)
    distribution = read_rows('launch_distribution.csv')
    for r in distribution:
        if float(old[r['case'],r['arm'],r['event_class']]['counted_parents']) == 0:
            r['fraction_among_class'] = ''
    c.table(out/'launch_distribution.csv',distribution)
    c.dump(out/'impact_height_consumer.json',dict(cases=len(models),impact_rows=len(rows),
        same_parent_counts=True,max_source_work_difference_J=maximum_work_difference,
        reused='frozen convergence labels and all regional signed/energy tails; no tail or physics refinement'))
    postprocess(out, model=model)
    print(json.dumps(dict(impact_height_rows=len(rows),max_source_work_difference_J=maximum_work_difference)))


def postprocess(out, model=None):
    """Consume frozen wave results and reconstruct only their saved source inputs."""
    import csv
    global DT, EXPOSURE, REPS
    config = json.loads((out/'contract.json').read_text())['sample']
    DT, EXPOSURE, REPS = config['dt_s'], config['exposure_s'], config['independent_records']
    read_rows = lambda name: list(csv.DictReader((out/name).open()))
    models = read_rows('models.csv')
    impacts = read_rows('impact_classes.csv')
    convergence_rows = read_rows('convergence.csv')
    if 'thermal_mean_point_J_m3' not in convergence_rows[0]:
        if model is None:
            atlas = read_rows('mode_attenuation_map.csv')
            frequencies = [float(r['frequency_Hz']) for r in atlas]
            model = c.assemble('grains', c.DIMS[0]/8, min(frequencies)*(1-1e-8), max(frequencies)*(1+1e-8))
        point_energy = off_beam_receiver(model)[4]
        profiles, _ = c.profiles_and_stopping()
        for r in convergence_rows:
            energy = .25 if r['case']=='PdD_0.25keV' else 2.5
            temperature = profiles['PdD',energy]['case']['T_K']
            mean = KB*temperature/2*np.trace(point_energy)
            sd = KB*temperature/math.sqrt(2)*np.linalg.norm(point_energy)
            r.update(thermal_mean_point_J_m3=float(mean),thermal_sd_point_J_m3=float(sd),
                peak_B_over_point_beam_mean=float(r['peak_B_quadratic_J_m3'])/float(r['mean_total_quadratic_J_m3']),
                peak_B_over_point_thermal_mean=float(r['peak_B_quadratic_J_m3'])/mean,
                peak_B_over_point_thermal_sd=float(r['peak_B_quadratic_J_m3'])/sd,
                height_observer='receiver-center quadratic density; one-point box cubature, not full region integral')
        c.table(out/'convergence.csv',convergence_rows)
    signed = read_rows('signed_tails.csv')
    # These are the same seeded Cox profiles, not a replay of the field producer.
    temporal = []
    for r in models:
        burn = round(5*float(r['mode_lifetime_max_s'])/DT)
        n = round(EXPOSURE/DT)
        weight = beam_timing(burn+n, int(r['temporal_clumping']), np.random.default_rng(14109))[burn:]
        weight /= weight.mean(axis=0)
        lag = max(1, round(1e-6/DT))
        correlations = [float(np.corrcoef(weight[:-lag, j], weight[lag:, j])[0, 1])
                        for j in range(REPS) if np.std(weight[:, j]) > 1e-12]
        temporal.append(dict(case=r['case'], arm=r['arm'],
            temporal_clumping=r['temporal_clumping'], source_rate_modulation_cv=float(np.std(weight)),
            source_rate_modulation_peak_over_mean=float(weight.max()),
            source_rate_lag1us_correlation=float(np.mean(correlations)) if correlations else '',
            source='reconstructed exact seeded normalized input path; no wave propagation'))
    c.table(out/'temporal_clumping_inputs.csv', temporal)
    impact_map = {(r['case'],r['arm'],r['event_class']):r for r in impacts}
    bmap = {(r['case'],r['arm']):r for r in convergence_rows}
    signed_map = {(r['case'],r['arm'],r['observable'],r['sign']):r for r in signed
                  if float(r['side_mm']) == 2.5 and float(r['z0_um']) == 0.
                  and float(r['threshold_over_causal_prior_rms']) == 3.}
    summary = []
    for r in models:
        if float(r['intrinsic_lifetime_s']) != 50e-6 or float(r['range_proxy_m']) != .2 or float(r['feedback']) != .35 or float(r['lateral_edge_energy_reflectivity']) != 1.:
            continue
        case, arm = r['case'], r['arm']
        b = bmap[case, arm]
        item = dict(case=case, arm=arm, spatial_clumping=r['spatial_clumping'], temporal_clumping=r['temporal_clumping'],
                    spatial_peak_over_mean=r['spatial_profile_peak_over_mean'],
                    counted_parent_charge_impacts=impact_map[case,arm,'all']['counted_parents'],
                    B_receiver_fraction=b['time_fraction'], B_initiation_rate_s=b['initiation_rate_s'],
                    B_peak_quadratic_J_m3=b['peak_B_quadratic_J_m3'])
        for label in ('A_disturbed','C_any','C_strain_positive','C_strain_negative','C_displacement_positive','C_displacement_negative','A_and_C','A_or_C'):
            p = impact_map[case,arm,label]
            item[label+'_fraction'] = p['fraction_among_impacts']
            item[label+'_impact_rate_s'] = p['impact_rate_s']
            item[label+'_mean_launch_eV'] = p['mean_launch_self_eV_parent']
        for observable in ('epsilon_xx','u_z'):
            for sign in ('positive_peak','negative_hole'):
                s = signed_map[case,arm,observable,sign]
                item[observable+'_'+sign+'_initiation_rate_s'] = s['initiation_rate_s']
                item[observable+'_'+sign+'_physical_extreme'] = s['physical_extreme']
        summary.append(item)
    base = {r['case']:r for r in summary if r['spatial_clumping']=='0' and r['temporal_clumping']=='0'}
    for r in summary:
        for metric in ('A_disturbed_impact_rate_s','C_any_impact_rate_s','B_initiation_rate_s',
                       'epsilon_xx_positive_peak_initiation_rate_s','epsilon_xx_negative_hole_initiation_rate_s'):
            denominator = float(base[r['case']][metric])
            r[metric+'_ratio_to_uniform_steady'] = float(r[metric])/denominator if denominator else ''
    c.table(out/'clumping_summary.csv', summary)
    common = []
    for case in sorted({r['case'] for r in models}):
        settings = sorted({(float(r['feedback']),int(r['spatial_clumping']),int(r['temporal_clumping']),float(r['lateral_edge_energy_reflectivity'])) for r in models if r['case']==case})
        for label in ('A_disturbed','C_strain_positive','C_strain_negative','C_displacement_positive','C_displacement_negative','C_any','B_offimpact_convergence'):
            population = convergence_rows if label=='B_offimpact_convergence' else impacts
            fraction_key = 'time_fraction' if label=='B_offimpact_convergence' else 'fraction_among_impacts'
            for feedback,spatial,timing,edge in settings:
                selected = [r for r in population if r['case']==case
                    and float(r['feedback'])==feedback and int(r['spatial_clumping'])==spatial
                    and int(r['temporal_clumping'])==timing and float(r['lateral_edge_energy_reflectivity'])==edge
                    and (label=='B_offimpact_convergence' or r['event_class']==label)]
                for length in sorted({float(r['range_proxy_m']) for r in selected}):
                    subset = [r for r in selected if float(r['range_proxy_m'])==length]
                    reached = [r for r in subset if float(r[fraction_key])>=.5]
                    common.append(dict(case=case,event_class=label,feedback=feedback,spatial_clumping=spatial,
                        temporal_clumping=timing,lateral_edge_energy_reflectivity=edge,range_proxy_m=length,
                        denominator='off-impact receiver time windows' if label=='B_offimpact_convergence' else 'physical parent-charge impacts',
                        common_fraction=.5,tested_lifetimes_s=';'.join(str(v) for v in sorted({float(r['intrinsic_lifetime_s']) for r in subset})),
                        first_tested_common_lifetime_s=min((float(r['intrinsic_lifetime_s']) for r in reached),default=''),
                        max_observed_fraction=max(float(r[fraction_key]) for r in subset)))
    c.table(out/'common_map.csv',common)
    metadata = json.loads((out/'beam_parameters.json').read_text())['plate']
    timeframes(metadata, out)
    print(json.dumps(dict(clumping_summary_rows=len(summary), source_input_rows=len(temporal), field_replayed=False)))


def run(args):
    global DT, EXPOSURE, REPS, ETA
    DT, EXPOSURE, REPS, ETA = args.dt, args.exposure, args.realizations, args.eta_band
    output = Path(args.output)
    if not output.is_absolute():
        output = c.ROOT/output
    out = output/'live'
    out.mkdir(parents=True, exist_ok=True)
    profiles, stopping = c.profiles_and_stopping()
    model = c.assemble('grains', c.DIMS[0]/8, args.fmin, args.fmax)
    mode_metadata = mode_group_velocities(model)
    c.table(out/'mode_attenuation_map.csv', mode_metadata['rows'])
    eps, uz, group = probes(model)
    timeframes(mode_metadata, out)
    st_ref = stopping['PdD', 2.5]
    sr = c.source_matrix(model, [c.deposition_depth(t) for t in st_ref['tracks']])
    eparent_ref = np.array([st_ref['tracks'][0]['components_eV'][1], 2*st_ref['tracks'][1]['components_eV'][1]])
    energy_weights = np.array([.85, .15])*eparent_ref
    energy_weights /= energy_weights.sum()
    weights = np.repeat(energy_weights/36, 36)
    quiet_fraction = float(np.dot([.25,.74,.01], np.minimum(1., ETA/.865*np.array([.1,1.,10.]))))
    reference = source_reference(model, sr, weights, eps, uz, quiet_fraction*138e-6*st_ref['host_recoil_J_per_true_C'])
    receiver = off_beam_receiver(model)
    regions, omitted = [], []
    for side in (.0025, .01, float(c.DIMS[0])):
        for z0, z1 in ((0., 50e-6), (25e-6, 75e-6), (175e-6, 225e-6)):
            e, u = linear_region(model, side, z0, z1)
            supported = np.linalg.norm(u) > 1e-8*np.max(np.linalg.norm(uz, axis=1))
            if not supported:
                omitted.append(dict(observable='u_z', side_m=side, z0_m=z0, z1_m=z1, reason='near-zero supported spatial-average displacement projector; no normalized tail claimed'))
            regions.append(dict(side=side, z0=z0, z1=z1, projectors=(e, u), uz_supported=supported,
                                energy=c.receiver(model, side, z0, z1)))
    impact_rows, distribution, b_rows, hole_rows, energy_tails, peaks, model_rows = [], [], [], [], [], [], []
    histories = {}
    for energy, current in ((2.5, 138e-6), (.25, 38e-6)):
        st = stopping['PdD', energy]
        eh = np.array([st['tracks'][0]['components_eV'][1], 2*st['tracks'][1]['components_eV'][1]])*QE
        eh *= st['host_recoil_J_per_true_C']/(np.dot([.85, .15], eh)/QE)
        sm = c.source_matrix(model, [c.deposition_depth(t) for t in st['tracks']])
        source = dict(eh=eh, s=(sm[:, :36], sm[:, 36:]))
        sweep = [(t, length, feedback, 0, 0) for t in (2e-6, 50e-6, 200e-6) for length in (.002, .02, .2, 2.) for feedback in (0., .35)]
        sweep += [(t, length, .35, 2, 2) for t in (2e-6, 50e-6, 200e-6) for length in (.002, .02, .2, 2.)]
        sweep += [(50e-6, .2, .35, s, t) for s in range(3) for t in range(3)]
        sweep = list(dict.fromkeys(sweep))
        if energy == .25:
            sweep = [(50e-6, .2, .35, 0, 0), (50e-6, .2, .35, 2, 2), (200e-6, 2., .35, 2, 2)]
        sweep = [(*parameters, 1.) for parameters in sweep]
        sweep += [(200e-6, 2., .35, 2, 2, .1)]
        for tau, length, feedback, spatial, timing, edge_reflectivity in sweep:
            arm = f'live_tau{tau:g}_L{length:g}_feedback{feedback:g}_space{spatial}_time{timing}_edge{edge_reflectivity:g}'
            case = f'PdD_{energy}keV'
            r = simulate(model, source, group, eps, uz, reference, receiver, regions, energy, current, tau, length, feedback, spatial, timing, edge_reflectivity)
            if abs(r['energy_balance_residual_J']) > 1e-9*r['incoming_host_J'] or r['max_budget_violation_J'] > 1e-18:
                raise AssertionError('Live source/propagation energy ledger did not close')
            parents = r['impact']['all']['parents']
            for label, a in r['impact'].items():
                impact_rows.append(dict(case=case, arm=arm, event_class=label, intrinsic_lifetime_s=tau,
                    spatial_clumping=spatial, temporal_clumping=timing, feedback=feedback,
                    lateral_edge_energy_reflectivity=edge_reflectivity,
                    range_proxy_m=length, effective_energy_lifetime_s=r['tau_eff'],
                    counted_parents=a['parents'], all_counted_parents=parents,
                    fraction_among_impacts=a['parents']/parents,
                    physical_source_exposure_s=EXPOSURE*REPS,
                    impact_rate_s=a['parents']/(EXPOSURE*REPS),
                    mean_launch_self_eV_parent=a['launch_self_J']/a['parents']/QE if a['parents'] else '',
                    max_supported_launch_self_eV_parent=a['max_launch_self_J']/QE if a['parents'] else '',
                    common_50pct=(a['parents']/parents >= .5)))
                impact_rows[-1].update(impact_heights(a, r['impact']['all']['parents_by_record']))
                for mark, count in enumerate(a['marks']):
                    distribution.append(dict(case=case, arm=arm, event_class=label, mark=('low','middle','high')[mark],
                        count=float(count), fraction_among_class=float(count/a['parents']) if a['parents'] else ''))
            bstat = occurrence(r['B'], DT)
            b_rows.append(dict(case=case, arm=arm, intrinsic_lifetime_s=tau, range_proxy_m=length,
                spatial_clumping=spatial, temporal_clumping=timing, feedback=feedback,
                lateral_edge_energy_reflectivity=edge_reflectivity,
                effective_energy_lifetime_s=r['tau_eff'], off_impact=True,
                directional_only_time_fraction=float(r['directional'].mean()),
                constructive_only_time_fraction=float(r['constructive'].mean()),
                peak_total_quadratic_J_m3=float(r['off_u'].max()),
                mean_total_quadratic_J_m3=float(r['off_u'].mean()),
                peak_total_over_separate_source_energy=float(np.max(r['off_u']/np.maximum(r['off_separate'], np.finfo(float).tiny))),
                peak_B_quadratic_J_m3=float(np.max(np.where(r['B'], r['off_u'], 0))),
                distinct_origin_groups=4, minimum_distinct_origins=2, min_inward_normal_angle_deg=90.,
                common_50pct=(bstat['time_fraction'] >= .5), **bstat))
            for k, region in enumerate(regions):
                side, z0, z1 = region['side'], region['z0'], region['z1']
                for d, observable in enumerate(('epsilon_xx','u_z')):
                    if d == 1 and not region['uz_supported']:
                        continue
                    proj = region['projectors'][d]
                    thermal = math.sqrt(KB*profiles['PdD', energy]['case']['T_K']*np.sum(proj*proj))
                    hole_rows += signed_tails(r['signed'][:, :, 2*k+d], r['z'][:, :, 2*k+d], thermal, observable, side, z0, z1, case, arm)
                ar = region['energy']
                values = r['ecanon']/2*np.einsum('tri,ij,trj->tr', r['x'], ar, r['x'], optimize=True)
                thermal_mean = KB*profiles['PdD', energy]['case']['T_K']/2*np.trace(ar)
                thermal_sd = KB*profiles['PdD', energy]['case']['T_K']/math.sqrt(2)*np.linalg.norm(ar)
                energy_tails += c.tail_rows(values, thermal_mean, thermal_sd, side, z0, z1, arm, case, DT)
                nd = c.column_density(profiles['PdD', energy]['layers'], z0, z1)
                conversion = float(values.max())/r['tau_eff']*4.1e-12/nd/QE
                peaks.append(dict(case=case, arm=arm, side_mm=side*1e3, z0_um=z0*1e6, z1_um=z1*1e6,
                    n_D_m3=nd, peak_J_m3=float(values.max()), mean_J_m3=float(values.mean()),
                    peak_acoustic_quadratic_eV_per_D=float(values.max())/nd/QE,
                    conditional_peak_conversion_eV_D=conversion,
                    concentration_required_from_modal_peak=c.REFERENCE/conversion,
                    hypothetical_sustained_peak_holding_required_s=4.1e-12*c.REFERENCE/conversion))
            model_rows.append(dict(case=case, arm=arm, energy_lifetime_s=r['tau_eff'],
                intrinsic_lifetime_s=tau, range_proxy_m=length, feedback=feedback,
                spatial_clumping=spatial, temporal_clumping=timing,
                lateral_edge_energy_reflectivity=edge_reflectivity,
                mode_lifetime_min_s=r['mode_lifetime_min_s'], mode_lifetime_max_s=r['mode_lifetime_max_s'],
                source_energy_balance_residual_J=r['energy_balance_residual_J'],
                spatial_profile_peak_over_mean=r['spatial_profile_peak_over_mean'], spatial_profile_cv=r['spatial_profile_cv'],
                source_work_J=r['source_work_J'], incoming_host_J=r['incoming_host_J'],
                marked_self_J=r['marked_self_J'], damping_loss_J=r['damping_loss_J'],
                capped_record_bins=r['capped_record_bins'], max_budget_violation_J=r['max_budget_violation_J']))
            histories[case+'_'+arm+'_signed'] = r['signed'][:, :, :2].astype(np.float32)
            histories[case+'_'+arm+'_B'] = r['B']
    c.table(out/'impact_classes.csv', impact_rows)
    c.table(out/'launch_distribution.csv', distribution)
    c.table(out/'convergence.csv', b_rows)
    c.table(out/'signed_tails.csv', hole_rows)
    c.table(out/'tails.csv', energy_tails)
    c.table(out/'region_peaks.csv', peaks)
    c.table(out/'models.csv', model_rows)
    checks = dict(covariance_relative_residual=reference[-1],
        max_source_budget_violation_J=max(r['max_budget_violation_J'] for r in model_rows),
        max_energy_balance_residual_over_incoming_host=max(abs(r['source_energy_balance_residual_J'])/r['incoming_host_J'] for r in model_rows))
    canonical = [r for r in impact_rows if r['event_class']=='all' and r['case']=='PdD_2.5keV'
                 and r['intrinsic_lifetime_s']==50e-6 and r['range_proxy_m']==.2
                 and r['feedback']==.35 and r['lateral_edge_energy_reflectivity']==1.]
    checks['canonical_factorial_unique_total_parent_counts'] = len({r['counted_parents'] for r in canonical})
    if checks['canonical_factorial_unique_total_parent_counts'] != 1:
        raise AssertionError('Beam clumping comparison changed total charge')
    impact_index = {(r['case'],r['arm'],r['event_class']):r for r in impact_rows}
    overlap_error = 0.
    for r in model_rows:
        get = lambda label: impact_index[r['case'],r['arm'],label]['fraction_among_impacts']
        overlap_error = max(overlap_error, abs(get('A_or_C')-get('A_disturbed')-get('C_any')+get('A_and_C')))
    checks['max_A_C_union_identity_error'] = overlap_error
    c.dump(out/'behavioral_checks.json', checks)
    c.event_requirements(out)
    np.savez_compressed(out/'selected_histories.npz', **histories)
    c.dump(out/'contract.json', dict(schema='g13-live-4109-v1', steering=['entry_4097','entry_4103','entry_4104','entry_4109'],
        damage='omitted by direct corrected steering; no damage/dose interpretation of rogue holes',
        model='same 53-mode rotated-cubic-grain/free-face finite variational operator as fixed-source baseline',
        causal_source='pre-bin wave state conditions every physical parent mark; no new event qualifies itself; within 50ns bin state frozen',
        source_spread='three conditional energy-fraction marks and independent signed counts; exact Poisson parent counts and conditional multinomial marks',
        energy='each raw self-energy fraction <=1; coherent kinetic cross work accounted once, batch new work capped by total incoming host energy; negative work returns to outgoing reservoir',
        source_identity='disjoint spatial parent sectors; identities sector/site/bin/species/sign/ordinal; echoes retain original sector',
        convergence='two different origin sectors have inward flux on receiver faces with >=90-degree inward-normal separation; each exceeds 5% total inward power; total strain energy >1.1 sum separate sector energies',
        receiver_center_m=receiver[-2], receiver_size_m=receiver[-1],
        A='pre-impact |epsilon_xx| or |u_z| > .25 fixed quiet canonical rms in 2.5mm-square, depth0-50um neighborhood centered at source',
        C='pre-impact signed epsilon_xx or u_z beyond +/-3 causal-prior rms; positive and negative identical thresholds; ensemble zero reference',
        rms='prior EWMA squared field, time scale 5*max modal lifetime, initialized from quiet same-band covariance; score then update; no future/current launch enters predicate',
        lifetime_range='mode lifetime^-1=tau_intrinsic^-1+v_group_mode/L_energy; energy survival exp(-path/L_energy); thin-plate flexural group proxy and Hill bulk speed mapping, not optical lifetime',
        beam='uniform/peaked/peaked+hotspot space arms; steady/ripple/random-burst Cox time arms; shared Poisson charge totals and normalized intensity paths make matched charge exact',
        sample=dict(dt_s=DT, exposure_s=EXPOSURE, independent_records=REPS, burn='5*max modal lifetime', conditioned_windows=False, charge_conditioning='shared Poisson record totals, not selected wave events'),
        common='50% among physical impacts for A/C; 50% of off-impact receiving time windows for B; not an additive shared event denominator',
        covariance_relative_residual=reference[-1], omitted_projectors=omitted,
        lateral_edges='free variational lateral reflections retained; separate edge=.1 arm adds modal round-trip leakage -log(R_E)*v_group/model_side; boundary shape remains the equal-area square, not measured flange',
        local_conversion='positive quadratic-loss lifetime proxy using median modal lifetime, not a calibrated spatial acoustic-to-D dissipator; per-mode actual total damping is separately budgeted',
        limits='finite-band coarse influx convergence, not resolved atomistic splash fronts; new source law is a bounded parameterized hypothesis, not calibrated acoustic yield'))
    postprocess(out, model=model)
    print(json.dumps(dict(live_cases=len(model_rows), impact_rows=len(impact_rows), signed_tail_rows=len(hole_rows),
        peak_B_fraction=max(r['time_fraction'] for r in b_rows),
        peak_C_fraction=max(r['fraction_among_impacts'] for r in impact_rows if r['event_class']=='C_any'),
        max_source_budget_violation_J=max(r['max_budget_violation_J'] for r in model_rows)), indent=2))
