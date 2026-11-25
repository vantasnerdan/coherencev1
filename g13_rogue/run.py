"""Build the seven-item, finite-model Pass 1 capsule from repository inputs."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import csv
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode = True
import numpy as np
from scipy.constants import atomic_mass, electron_volt, hbar
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize
from types import SimpleNamespace
from scipy.interpolate import CubicHermiteSpline

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / 'g13_rogue'
OUT = HERE / 'output'
HBAR = hbar / electron_volt * 1e15
KB = 8.617333262145e-5
TAU = 8200.  # engine drag time 4100 fs -> amplitude lifetime 8200 fs
AXIS = np.array([1., 1., 0.]) / np.sqrt(2)
BASE_SHA = '02fe7cc19a121706b70ca201195b36fc779e0fb8'
INPUTS = {}


def bind(path):
    path = Path(path)
    if not path.is_absolute():
        path = ROOT / path
    data = path.read_bytes()
    INPUTS[str(path.relative_to(ROOT))] = dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    return data


def module(name, path):
    bind(path)
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def emit(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def table(name, rows):
    with (OUT / name).open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def frozen_rays(model):
    """Scalar restriction of the accepted Hermite anchors, no external overlays."""
    anchors = json.loads(bind('g10_quantum_site/bulk-refit-080/results/potential.json'))['anchors']
    k = float(model['K'][0,0,0])
    quartics = (model['a'],(model['a']+model['b'])/2,(model['a']+2*model['b'])/3)
    splines = []
    for ray,k4 in zip(('100','110','111'),quartics):
        points = sorted((p for p in anchors if p['ray']==ray),key=lambda p:p['r_A'])
        radii = [.03]+[p['r_A'] for p in points]
        values = [.5*k*.03**2+k4*.03**4/4]+[p['energy_eV'] for p in points]
        slopes = [k*.03+k4*.03**3]+[p['slope_eV_A'] for p in points]
        splines.append(CubicHermiteSpline(radii,values,slopes,extrapolate=False))
        assert np.max(abs(splines[-1](radii)-values)) < 1e-12
        assert np.max(abs(splines[-1](radii,1)-slopes)) < 1e-12
    def radial(radius,variant):
        if variant != 'wide_quadratic':
            raise ValueError('Only the accepted quadratic exterior is selected')
        radius = np.asarray(radius)
        values,slopes = [],[]
        for spline,k4 in zip(splines,quartics):
            r = np.clip(radius,.03,.8)
            energy,slope = spline(r),spline(r,1)
            energy = np.where(radius < .03,.5*k*radius**2+k4*radius**4/4,energy)
            slope = np.where(radius < .03,k*radius+k4*radius**3,slope)
            outside = np.maximum(radius-.8,0.)
            energy = np.where(radius > .8,spline(.8)+spline(.8,1)*outside+.5*k*outside**2,energy)
            slope = np.where(radius > .8,spline(.8,1)+k*outside,slope)
            values.append(energy)
            slopes.append(slope)
        return np.stack(values,axis=-1),np.stack(slopes,axis=-1)
    return SimpleNamespace(radial=radial)


def setup():
    m = module('g13_saved_graph', 'g10_quantum_site/model.py')
    bind('g10_quantum_site/bulk-refit-080/well.py')
    w = frozen_rays(m.full())
    beam = module('g13_beam', 'g4_beam/forward/model.py')
    for path in ('g10_quantum_site/bulk-refit-048/well.py',
                 'g10_quantum_site/bulk-refit-048/inputs/old-readout.json',
                 'g10_quantum_site/bulk-refit-048/inputs/bulk-wells.json',
                 'g10_quantum_site/bulk-wide-wells/bulk-wells.json',
                 'g12_residual/inputs_ledger.md', 'g6_response/README.md',
                 'g6_response/inputs-ledger.json', 'g11_cmns/collective/focusing/engine.py',
                 'g11_cmns/collective/part1/motion.py', 'g11_cmns/collective/part1/stopping.json'):
        bind(path)
    models = [m.full()] + m.partial_cases()
    for model in models:
        bind(model['source'])
    mass = m.MASS['D'] * atomic_mass * 1e10 / electron_volt
    return m, w, beam, models, mass


def transitions(w, mass):
    """Actual 1-D [110] nuclear eigenlevels, not a free-energy curvature."""
    t = HBAR**2 / (2*mass)
    dx = .002
    x = np.arange(-2. + dx, 2., dx)
    v = w.radial(abs(x), 'wide_quadratic')[0][:, 1]
    en, vec = eigh_tridiagonal(2*t/dx**2 + v, np.full(len(x)-1, -t/dx**2),
                               select='i', select_range=(0, 7))
    gap = np.diff(en)
    e0, chi = float(gap[0]), float(gap[1]-gap[0])
    rows = [dict(level=i, energy_eV=float(e), excitation_eV=float(e-en[0]),
                 next_gap_eV=float(gap[i]) if i < len(gap) else None,
                 probability_abs_u_gt_080=float(np.sum(vec[abs(x) > .8, i]**2)))
            for i, e in enumerate(en)]
    # Classical action integral on the SAME ray. End-point substitution removes
    # the square-root singularity without subtracting near-equal potentials.
    z, wt = np.polynomial.legendre.leggauss(96)
    theta = (z+1)*np.pi/4
    weights = wt*np.pi/4
    actions = []
    for amp in np.linspace(.01, .8, 160):
        pos = amp*np.sin(theta)
        energy = float(w.radial(np.array(amp), 'wide_quadratic')[0][1])
        diff = energy-w.radial(pos, 'wide_quadratic')[0][:, 1]
        jac = amp*np.cos(theta)
        action = 2/np.pi*np.sum(weights*jac*np.sqrt(np.maximum(diff, 0)/t))
        quantum = np.pi*np.sqrt(t)/np.sum(weights*jac/np.sqrt(diff))
        actions.append(dict(turning_A=amp, classical_energy_eV=energy,
                            action_quanta=float(action), action_frequency_eV=float(quantum)))
    table('transitions.csv', rows)
    table('action_frequency.csv', actions)
    return e0, chi, t, rows, actions


def image_stencil(model):
    """Average source shell mean times every actually occupied bond, once."""
    p, J = model['positions'], model['J']
    cell = np.asarray(json.loads((ROOT/model['source']).read_text())['profile']['geometry']['cell_A'])
    side = np.diag(cell)
    vectors, weights = [], []
    for i in range(len(p)):
        for j in range(len(p)):
            if i == j:
                continue
            d = p[j]-p[i]
            d -= np.rint(d/side)*side
            choices = [(-side[k]/2, side[k]/2) if np.isclose(abs(d[k]), side[k]/2) else (d[k],)
                       for k in range(3)]
            images = list(itertools.product(*choices))
            for delta in images:
                vectors.append(delta)
                weights.append(J[i, j]/(len(p)*len(images)))
    return np.array(vectors), np.array(weights)


def dispersion(k, vectors, weights, e0, t):
    return e0 - t/e0 * (np.cos(np.asarray(k)[..., None]* (vectors@AXIS))@weights)


def loading_and_growth(m, models, e0, chi, t, actions):
    records, growth, curves = [], [], []
    for model in models:
        n = len(model['J'])
        # Preserve all measured F0 and Kii in the finite graph diagnostic.
        def objective(flat):
            energy, force = m.potential(flat.reshape(n, 3), model)
            return float(energy), -force.ravel()
        fit = minimize(objective, np.zeros(3*n), jac=True, method='L-BFGS-B',
                       options=dict(gtol=1e-10, maxiter=1000, ftol=1e-14))
        u = fit.x.reshape(n, 3)
        xx, r2 = u*u, np.sum(u*u, axis=1)
        local = model['K'].copy()
        for i in range(n):
            local[i] += np.diag(3*(model['a']-model['b'])*xx[i] + model['b']*r2[i])
            local[i] += 2*model['b']*np.outer(u[i], u[i])
        kparallel = np.einsum('a,iab,b->i', AXIS, local, AXIS)
        # Partial loading uses the explicit old cubic equilibrium closure. Its
        # translational average is not advertised as a disorder eigenmode.
        effective_e0 = e0 + t/e0*(float(np.mean(kparallel))-float(models[0]['K'][0, 0, 0]))
        vectors, weights = image_stencil(model)
        a = float(np.diag(np.asarray(json.loads((ROOT/model['source']).read_text())['profile']['geometry']['cell_A']))[0]/2)
        k = np.linspace(0., 2*np.pi*np.sqrt(2)/a, 4097)
        eps = dispersion(k, vectors, weights, effective_e0, t)
        top = int(np.argmax(eps))
        k0 = float(k[top])
        curves.extend(dict(case=model['name'], k_A_inverse=float(q), epsilon_eV=float(e))
                      for q, e in zip(k[::16], eps[::16]))
        support = .8 if model['name']=='Pd' else model['support_A']
        records.append(dict(case=model['name'], loading=n/32, occupied_sites=n,
            local_Kii_mean_eV_A2=float(np.trace(model['K'], axis1=1, axis2=2).mean()/3),
            source_sumJ_eV_A2=model['sumJ_source'], actual_graph_mean_sumJ_eV_A2=float(model['J'].sum()/n),
            bare_min_eigenvalue_eV_A2=float(np.linalg.eigvalsh(np.diag(np.einsum('a,iab,b->i', AXIS, model['K'], AXIS))-model['J'])[0]),
            equilibrium_max_radius_A=float(np.linalg.norm(u, axis=1).max()),
            equilibrium_gradient_max_eV_A=float(abs(fit.jac).max()),
            translated_effective_gap_eV=effective_e0, transferred_chi_eV_per_quantum=chi,
            source_support_A=support, top_k_A_inverse=k0,
            top_carrier_wavelength_nm=2*np.pi/k0/10 if k0 else None,
            minimum_epsilon_eV=float(eps.min()), maximum_epsilon_eV=float(eps.max())))
        for tau in (2050., 4100., 8200.):
            for wavelength in np.geomspace(.6, 200., 96):
                p = 2*np.pi/(10*wavelength)
                for carrier, kc in (('top', k0), ('Gamma', 0.)):
                    d = float((dispersion(kc+p,vectors,weights,effective_e0,t)+
                               dispersion(kc-p,vectors,weights,effective_e0,t))/2-
                              dispersion(kc,vectors,weights,effective_e0,t))
                    threshold = ((-d)**2+(HBAR/tau)**2)/(2*chi*(-d)) if d < 0 and chi > 0 else None
                    for action in (0., .001, .01, .03, .1, .3, 1.):
                        radicand = -d*(d+2*chi*action)
                        gain = np.sqrt(max(0.,radicand))/HBAR
                        amp = 2*np.sqrt(t/effective_e0*action)
                        growth.append(dict(case=model['name'], carrier=carrier, amplitude_lifetime_fs=tau,
                            wavelength_nm=float(wavelength), action_quanta=action,
                            energy_meV_per_D=1000*(effective_e0*action+.5*chi*action*action),
                            displacement_peak_A=float(amp), signed_sideband_d_eV=d,
                            growth_per_ps=float(1000*gain), net_growth_per_ps=float(1000*(gain-1/tau)),
                            threshold_action_quanta=float(threshold) if threshold is not None else None,
                            threshold_meV_per_D=float(1000*(effective_e0*threshold+.5*chi*threshold**2)) if threshold is not None else None,
                            threshold_displacement_A=float(2*np.sqrt(t/effective_e0*threshold)) if threshold is not None else None,
                            within_source_displacement=bool(amp+np.linalg.norm(u,axis=1).max() <= support),
                            within_Kerr_low_action=bool(action <= 1.)))
    table('loading.csv', records)
    table('growth.csv', growth)
    table('dispersion.csv', curves)
    return records




def prepared_seed(model, e0, chi, t, mass):
    vectors, weights = image_stencil(model)
    dx = 4.084/(2*np.sqrt(2))
    n = 256
    q = 2*np.pi*np.fft.fftfreq(n, dx)
    eps = dispersion(q,vectors,weights,e0,t)
    k = int(np.argmax(eps))
    kc = q[k]
    # Direct sideband formula avoids indexing ambiguity at a nonsymmetric top.
    d = (dispersion(kc+q,vectors,weights,e0,t)+dispersion(kc-q,vectors,weights,e0,t))/2-eps[k]
    thresholds = np.full(n, np.inf)
    sel = d < -1e-12
    thresholds[sel] = (d[sel]**2+(HBAR/TAU)**2)/(-2*chi*d[sel])
    p = int(np.argmin(thresholds))
    ncrit = float(thresholds[p])
    n0 = min(.3, max(.03, 8*ncrit))
    rng = np.random.default_rng(13013)
    x = np.arange(n)*dx
    carrier = np.exp(1j*kc*x)
    initial = np.sqrt(n0)*carrier*(1+.05*np.cos(q[p]*x)+.002*rng.normal(size=n))
    # A coherent pump sustains the unmodulated carrier; not a beam preparation.
    frame = eps[k]+chi*n0
    pump = np.sqrt(n0)*carrier/TAU
    dt, horizon = .5, 8200.
    linear = np.exp(-1j*(eps-frame)*dt/HBAR)
    field = initial.copy()
    scale = 2*np.sqrt(t/e0)
    rows, snapshots, trace_bank = [], [], []
    selected = None
    best = -1.
    # Pair spacing is 2*dx, exactly the bulk [110] nearest-neighbour separation.
    all_traces = np.empty((int(horizon/dt)+1,n,4))
    for step in range(int(horizon/dt)+1):
        time = step*dt
        physical = field*np.exp(-1j*frame*time/HBAR)
        derivative = -1j*(np.fft.ifft(eps*np.fft.fft(field))+chi*abs(field)**2*field)/HBAR-field/TAU+pump
        derivative *= np.exp(-1j*frame*time/HBAR)
        u, v = scale*physical.real, scale*derivative.real
        radius = 2*dx+np.roll(u,-2)-u
        closing = v-np.roll(v,-2)
        all_traces[step,:,0] = radius
        all_traces[step,:,1] = closing
        all_traces[step,:,2] = mass/4*closing**2
        all_traces[step,:,3] = abs(u).max()
        density = abs(field)**2
        if step % 20 == 0:
            energy = e0*density+.5*chi*density**2
            rows.append(dict(time_fs=time,max_action=float(density.max()),
                max_local_energy_eV=float(energy.max()), mean_local_energy_eV=float(energy.mean()),
                participation_sites=float(density.sum()**2/np.sum(density**2)),
                peak_width_equivalent_nm=float(dx*density.sum()/density.max()/10),
                max_displacement_A=float(abs(u).max())))
            snapshots.append(energy)
        if np.max(mass/4*np.maximum(closing,0)**2) > best:
            best = float(np.max(mass/4*np.maximum(closing,0)**2))
            selected = int(np.argmax(closing))
        if step == int(horizon/dt):
            break
        field = (field+.5*dt*pump)*np.exp(-dt/(2*TAU))
        field *= np.exp(-.5j*chi*abs(field)**2*dt/HBAR)
        field = np.fft.ifft(linear*np.fft.fft(field))
        field *= np.exp(-.5j*chi*abs(field)**2*dt/HBAR)
        field = field*np.exp(-dt/(2*TAU))+.5*dt*pump
    # A physical radial episode comes from this deterministic solved field,
    # not a stochastic bath increment or an invented cosine at its peak.
    data = all_traces[:,selected]
    trace = np.zeros((len(data),12))
    trace[:,0] = np.arange(len(data))*dt
    trace[:,1:3] = data[:,2,None]
    trace[:,3] = data[:,0]
    trace[:,4] = data[:,3]
    trace[:,5] = data[:,0]
    trace[:,8] = data[:,2]
    trace[:,9] = data[:,1]
    invalid = np.flatnonzero((trace[:,4] > .8)|(trace[:,3] < 1.4))
    first = int(invalid[0]) if len(invalid) else len(trace)
    saved_trace = trace[:first+1] if first < len(trace) else trace
    # Preserve the actual first excluded sample when there is an exit.
    replay = module('g13_episode_reader','g11_cmns/collective/focusing/pulse/replay.py')
    spec = dict(scored_samples=first,
                first_original_exit_fs=float(trace[first,0]) if first < len(trace) else None)
    path, turns, segment, summary = replay.episode(saved_trace,spec)
    np.savez_compressed(OUT/'histories.npz', prepared_seed=saved_trace)
    pair_path = 'g11_cmns/collective/contact/bulk_pair.json'
    pair_hash = hashlib.sha256(bind(pair_path)).hexdigest()
    emit('histories.json', dict(schema='g11-radial-histories-v1',pair_sha256=pair_hash,
        histories=[dict(case='prepared_modulated_carrier',launch='sustained_pump_random_seed',
        archived_trace='output/histories.npz',trace_key='prepared_seed',**spec)]))
    peak = max(rows,key=lambda r:r['max_local_energy_eV'])
    half = peak['max_local_energy_eV']/2
    over = np.array([r['max_local_energy_eV'] >= half for r in rows])
    center = rows.index(peak)
    lo = center
    hi = center
    while lo and over[lo-1]: lo -= 1
    while hi+1 < len(rows) and over[hi+1]: hi += 1
    summary.update(carrier_k_A_inverse=float(kc),modulation_wavelength_nm=float(2*np.pi/abs(q[p])/10),
        threshold_action_quanta=ncrit,initial_action_quanta=n0,
        initial_meV_per_D=float(1000*(e0*n0+.5*chi*n0*n0)),
        maximum=peak,half_peak_duration_fs=rows[hi]['time_fs']-rows[lo]['time_fs'],
        half_peak_duration_censored=bool(lo==0 or hi==len(rows)-1),
        longitudinal_projection='transversely uniform [110] slice; not a 3-D localized volume',
        source='explicit coherent pump and seeded carrier, not inferred beam production',
        first_support_exit_fs=float(trace[first,0]) if first < len(trace) else None,
        maximum_action_seen=float(max(r['max_action'] for r in rows)),
        Kerr_support_action_max=1.)
    table('seed.csv',rows)
    emit('seed.json',summary)
    np.savez_compressed(OUT/'seed_field.npz', energy_eV=np.asarray(snapshots),x_A=x,time_fs=np.arange(len(snapshots))*10.)
    return summary


def drive_sources(beam, model):
    # These exact fixed-track splits are read, not recomputed or refitted.
    saved = json.loads(bind('g11_cmns/collective/part1/stopping.json'))
    splits = [r for r in saved if r['target']=='PdD']
    # Extend the energy display, not the accepted track calculation: at 6.5 keV
    # keep the 2.5-keV partition ratios and geometry as an explicit proxy.
    high = json.loads(json.dumps(splits[-1]))
    factor = (13000.-115.)/high['sum_J_per_true_C']
    for key in ('Se_J_per_true_C','host_recoil_J_per_true_C','direct_D_recoil_J_per_true_C','sum_J_per_true_C'):
        high[key] *= factor
    high['initial_J_per_true_C'] = 13000.
    high['E_cm_keV'] = 6.5
    splits.append(high)
    rows = []
    for split in splits:
        ec = split['E_cm_keV']
        ic = float(beam.I_cup(ec))
        provenance = 'fixed saved track split' if ec <= 2.5 else 'scaled 2.5-keV partition and fixed 2.5-keV geometry; not a computed 6.5-keV track'
        fraction = float(beam.collected_fraction(ec,b_hi=1.,b_lo=0.))
        it = ic/fraction
        for current_scale in (.1,1.,10.):
            # Correct the source's wrongly named depth_A field from actual ODE units.
            depth_nm = float(split['tracks'][0]['terminal_depth_A'])*1e9
            assert depth_nm > .8168 and np.isclose(depth_nm, split['tracks'][0]['timeline'][-1]['depth_A']/10.)
            density = len(model['J'])/(.8168**3)
            assert np.isclose(density*.8168**3, 32.)
            for eta in (.001,.1):
                row = dict(E_cm_keV=ec,source_semantics=provenance,current_scale=current_scale,I_cup_uA=ic*current_scale,
                    cup_over_true=fraction,I_true_uA=it*current_scale,spot_cm2=1.,
                    nucleus_flux_nm2_s=it*current_scale*1e-6/electron_volt/1e14*1.15,
                    initial_eV_per_nucleus=split['initial_J_per_true_C']/1.15,
                    electron_eV_per_nucleus=split['Se_J_per_true_C']/1.15,
                    host_recoil_eV_per_nucleus=split['host_recoil_J_per_true_C']/1.15,
                    direct_D_eV_per_nucleus=split['direct_D_recoil_J_per_true_C']/1.15,
                    unresolved_eV_per_nucleus=split['residual_cutoff_J_per_true_C']/1.15,
                    eta_D=eta,eta_acoustic=eta,receiver_conversion=eta,
                    optical_packet_eV=eta*split['direct_D_recoil_J_per_true_C']/1.15,
                    acoustic_packet_eV=eta*split['host_recoil_J_per_true_C']/1.15,
                    track_depth_nm=depth_nm,D_density_nm3=density,
                    optical_mean_eV_per_D=(it*current_scale*1e-6/electron_volt/1e14*1.15)*
                       eta*split['direct_D_recoil_J_per_true_C']/1.15*(TAU/2*1e-15)/(density*depth_nm))
                rows.append(row)
    table('drive.csv',rows)
    transport=[]
    for row in rows:
        for tau_ac in (4100.,8200.,25000.,1250000.,2500000.):
            speed = .004  # nm/fs, a declared acoustic transport scenario
            la = speed*tau_ac
            le = la/2
            f = row['nucleus_flux_nm2_s']
            transport.append(dict(E_cm_keV=row['E_cm_keV'],current_scale=row['current_scale'],eta=row['eta_D'],
                acoustic_amplitude_lifetime_fs=tau_ac,amplitude_path_nm=la,energy_path_nm=le,
                causal_energy_weighted_active_count=2*np.pi*f*speed**2*(tau_ac/2)**3*1e-15,
                acoustic_mean_eV_per_D=f*row['acoustic_packet_eV']*(tau_ac/2*1e-15)/(row['D_density_nm3']*row['track_depth_nm']),
                converted_D_mean_eV_per_D=f*row['acoustic_packet_eV']*row['receiver_conversion']*(TAU/2*1e-15)/(row['D_density_nm3']*row['track_depth_nm'])))
    table('transport.csv',transport)
    nominal = 150e-6/electron_volt/1e14
    causal = []
    for length in (3.,100.,5000.,10000.):
        causal.append(dict(range_nm=length,nominal_charge_flux_nm2_s=nominal,
            hard_causal_interior_mean=float(np.pi*nominal*length**3/(3*.004)*1e-15),
            disk_times_residence_mean=float(np.pi*nominal*length**3/.004*1e-15)))
    table('causal_overlap.csv',causal)
    return [r for r in rows if r['current_scale']==1.]


def statistics(model, e0, chi, t, mass, sources):
    """Stratified single-arrival Palm response, not accelerated physical arrivals."""
    n = len(model['J'])
    L = e0*np.eye(n)-t/e0*model['J']
    energies, modes = np.linalg.eigh(L)
    thermal_n = 1/np.expm1(energies/(KB*300.))
    scale = 2*np.sqrt(t/e0)
    dt, horizon = 2.,4*TAU
    steps = int(horizon/dt)
    linear = (modes*np.exp(-1j*energies*dt/HBAR))@modes.T
    decay = np.exp(-dt/TAU)
    rng = np.random.default_rng(13991)
    samples = 12
    thermal = (rng.normal(size=(samples,n))+1j*rng.normal(size=(samples,n)))*np.sqrt(thermal_n/2)
    initial = thermal@modes.T
    # The same exact bath draws, phases, packet shapes and linear operator in both models.
    bath_sd = np.sqrt(thermal_n*(1-decay**2)/2)
    a = np.stack((initial.copy(),initial.copy()))
    for _ in range(8200):
        noise = ((rng.normal(size=(samples,n))+1j*rng.normal(size=(samples,n)))*bath_sd)@modes.T
        a[1] *= np.exp(-.5j*chi*abs(a[1])**2*dt/HBAR)
        a = (a@linear.T)*decay+noise
        a[1] *= np.exp(-.5j*chi*abs(a[1])**2*dt/HBAR)
    # Three radial area strata retain the near source and remote arrivals.
    side_nm = .8168
    area = side_nm**2
    source_radius = np.sqrt(area/np.pi)
    energy_length = .004*TAU/2
    bounds = (0.,source_radius,energy_length,5*energy_length)
    rmid = np.sqrt((np.square(bounds[:-1])+np.square(bounds[1:]))/2)
    areas = np.pi*(np.square(bounds[1:])-np.square(bounds[:-1]))
    p = model['positions']
    center = model['site']
    distance = p-p[center]
    distance -= np.rint(distance/8.168)*8.168
    local_shape = np.exp(-np.sum(distance**2,axis=1)/(2*2.887824**2)).astype(complex)
    configs=[]
    kicks=[]
    for source in sources:
        for zone,(radius,area_zone) in enumerate(zip(rmid,areas)):
            optical = source['optical_packet_eV']*side_nm/source['track_depth_nm'] if zone==0 else 0.
            # Energy-normalized 2-D annular shell at arrival, uniform through the
            # fixed track depth. The denominator is never smaller than cell area.
            annulus = max(area,2*np.pi*max(radius,source_radius)*side_nm)
            acoustic = source['receiver_conversion']*source['acoustic_packet_eV']*np.exp(-radius/energy_length)*area/annulus*side_nm/source['track_depth_nm']
            configs.append(dict(E_cm_keV=source['E_cm_keV'],eta=source['eta_D'],zone=zone,
                radius_nm=float(radius),annulus_area_nm2=float(annulus),area_stratum_nm2=float(area_zone),
                arrival_fs=float(radius/.004),packet_harmonic_energy_eV=float(optical+acoustic),
                optical_received_eV=float(optical),acoustic_received_eV=float(acoustic),
                flux_nm2_s=source['nucleus_flux_nm2_s']))
            wave = np.exp(1j*(distance@AXIS)*.25)
            shape = np.sqrt(optical)*local_shape/np.sqrt(np.vdot(local_shape,L@local_shape).real)
            shape += np.sqrt(acoustic)*wave/np.sqrt(np.vdot(wave,L@wave).real)
            # Optical/acoustic relative phase is not a closing projection.
            for_phase = np.exp(1j*rng.uniform(0,2*np.pi,samples))[:,None]*shape
            harmonic = np.einsum('bi,ij,bj->b',for_phase.conj(),L,for_phase).real
            for_phase *= np.sqrt((optical+acoustic)/harmonic)[:,None]
            kicks.append(for_phase)
    configs_count = len(configs)
    # index 0 is ambient; all other blocks use that same bath realization.
    fields = np.repeat(a[:,:,None,:],configs_count+1,axis=2)
    fields[:,:,1:,:] += np.moveaxis(np.array(kicks),0,1)[None,:,:,:]
    harmonic_work = np.einsum('bsci,ij,bscj->bsc',fields.conj(),L,fields).real
    total_work = harmonic_work.copy()
    total_work[1] += .5*chi*np.sum(abs(fields[1])**4,axis=-1)
    for c,config in enumerate(configs,1):
        for branch,label in enumerate(('linear','nonlinear')):
            work = total_work[branch,:,c]-total_work[branch,:,0]
            config['mean_injection_work_'+label+'_eV'] = float(work.mean())
            config['min_injection_work_'+label+'_eV'] = float(work.min())
            config['max_injection_work_'+label+'_eV'] = float(work.max())
    thresholds=np.r_[0.,np.geomspace(1e-6,2.,72)]
    pairs=np.array([(i,j) for i in range(n) for j in range(i+1,n)
                    if np.isclose(np.linalg.norm((p[j]-p[i])-np.rint((p[j]-p[i])/8.168)*8.168),2.887824,atol=2e-5)])
    bi,bj=pairs.T
    dr=p[bj]-p[bi];dr-=np.rint(dr/8.168)*8.168
    projection=(dr@AXIS)/np.linalg.norm(dr,axis=1)
    shape=(2,configs_count+1,len(thresholds))
    occupancy=np.zeros(shape); crossings=np.zeros(shape)
    pair_occ=np.zeros(shape);pair_cross=np.zeros(shape)
    previous=None
    max_local=np.zeros((2,samples,configs_count+1));max_close=max_local.copy()
    support_counts=np.zeros((2,configs_count+1));law_counts=support_counts.copy()
    total_samples=0
    for step in range(steps+1):
        action=abs(fields)**2
        local=e0*action+.5*chi*action**2
        # Linear comparison drops chi in both dynamics and reported Hamiltonian.
        local[0]=e0*action[0]
        derivative=-1j*(fields@L.T)/HBAR-fields/TAU
        derivative[1]-=1j*chi*action[1]*fields[1]/HBAR
        velocity=scale*derivative.real
        closing=(velocity[:,:,:,bi]-velocity[:,:,:,bj])*projection
        close=mass/4*np.maximum(closing,0)**2
        max_local=np.maximum(max_local,local.max(axis=-1))
        max_close=np.maximum(max_close,close.max(axis=-1))
        # Histogram bins preserve the whole empirical survival without allocating
        # samples x thresholds at every time step.
        indices=np.searchsorted(thresholds,local,side='right')
        pi=np.searchsorted(thresholds,close,side='right')
        for branch in range(2):
            for config in range(configs_count+1):
                counts=np.bincount(indices[branch,:,config].ravel(),minlength=len(thresholds)+1)
                occupancy[branch,config]+=np.cumsum(counts[::-1])[::-1][1:]
                counts=np.bincount(pi[branch,:,config].ravel(),minlength=len(thresholds)+1)
                pair_occ[branch,config]+=np.cumsum(counts[::-1])[::-1][1:]
                if previous is not None:
                    old,oldpair=previous
                    # Count each discrete upward threshold crossing, not occupied samples.
                    for new,prior,target in ((indices[branch,:,config],old[branch,:,config],crossings),
                                              (pi[branch,:,config],oldpair[branch,:,config],pair_cross)):
                        lower=prior.ravel();upper=new.ravel();rising=upper>lower
                        diff=np.bincount(lower[rising],minlength=len(thresholds)+1)-np.bincount(upper[rising],minlength=len(thresholds)+1)
                        target[branch,config]+=np.cumsum(diff)[:len(thresholds)]
                support_counts[branch,config]+=np.count_nonzero(scale*np.sqrt(action[branch,:,config])>.8)
                law_counts[branch,config]+=np.count_nonzero(action[branch,:,config]>1.)
        previous=(indices,pi)
        total_samples+=samples*n
        if step==steps: break
        noise=((rng.normal(size=(samples,n))+1j*rng.normal(size=(samples,n)))*bath_sd)@modes.T
        fields[1]*=np.exp(-.5j*chi*abs(fields[1])**2*dt/HBAR)
        fields=(fields@linear.T)*decay+noise[None,:,None,:]
        fields[1]*=np.exp(-.5j*chi*abs(fields[1])**2*dt/HBAR)
    window_s=horizon*1e-15
    occupancy/=total_samples
    pair_occ/=samples*len(pairs)*(steps+1)
    # Conditional response integral and count, averaged over phase realizations.
    crossings/=samples*n
    pair_cross/=samples*n  # pair initiation counts per D, not per bond
    records=[];survivals=[];rates=[]
    for branch,label in enumerate(('linear','nonlinear')):
        for c in range(configs_count+1):
            spec=configs[c-1] if c else {key: 0. for key in configs[0]}
            if c == 0:
                spec['zone'] = -1
            for k,threshold in enumerate(thresholds):
                survivals.append(dict(model=label,**spec,threshold_eV=float(threshold),
                    site_energy_occupancy=float(occupancy[branch,c,k]),closing_energy_pair_occupancy=float(pair_occ[branch,c,k]),
                    max_local_survival=float(np.mean(max_local[branch,:,c]>threshold)),
                    max_closing_survival=float(np.mean(max_close[branch,:,c]>threshold)),
                    local_upcrossings_per_D_window=float(crossings[branch,c,k]),
                    pair_upcrossings_per_D_window=float(pair_cross[branch,c,k]),
                    phases=samples,unresolved_max_tail=bool(np.count_nonzero(max_local[branch,:,c]>threshold)==0),
                    max_tail_one_count_resolution=1/samples))
        for source in sources:
            ids=[i+1 for i,c in enumerate(configs) if c['E_cm_keV']==source['E_cm_keV'] and c['eta']==source['eta_D']]
            for current_scale in (.1,1.,10.):
                lambdas=np.array([configs[i-1]['flux_nm2_s']*current_scale*configs[i-1]['area_stratum_nm2'] for i in ids])
                active=float(lambdas.sum()*window_s)
                for k,threshold in enumerate(thresholds):
                    delta=occupancy[branch,ids,k]-occupancy[branch,0,k]
                    pd=pair_occ[branch,ids,k]-pair_occ[branch,0,k]
                    dc=crossings[branch,ids,k]-crossings[branch,0,k]
                    pc=pair_cross[branch,ids,k]-pair_cross[branch,0,k]
                    rates.append(dict(model=label,E_cm_keV=source['E_cm_keV'],eta=source['eta_D'],current_scale=current_scale,
                        threshold_eV=float(threshold),ambient_site_occupancy=float(occupancy[branch,0,k]),
                        single_packet_site_occupancy_excess=float(lambdas@delta*window_s),
                        ambient_pair_closing_occupancy=float(pair_occ[branch,0,k]),
                        single_packet_pair_occupancy_excess=float(lambdas@pd*window_s),
                        ambient_local_initiations_per_D_s=float(crossings[branch,0,k]/window_s),
                        single_packet_local_initiations_excess_per_D_s=float(lambdas@dc),
                        ambient_pair_initiations_per_D_s=float(pair_cross[branch,0,k]/window_s),
                        single_packet_pair_initiations_excess_per_D_s=float(lambdas@pc),
                        two_or_more_arrivals_probability=float(-np.expm1(-active)-active*np.exp(-active)),
                        overlap_tail_coefficient=None,finite_memory_active_mean=active))
        records.append(dict(model=label,maximum_local_eV=float(max_local[branch].max()),
            maximum_closing_eV=float(max_close[branch].max()),
            coordinate_support_fraction=1-support_counts[branch].sum()/((configs_count+1)*total_samples),
            low_action_support_fraction=1-law_counts[branch].sum()/((configs_count+1)*total_samples)))
    table('survivals.csv',survivals)
    table('rates.csv',rates)
    table('packet_geometry.csv',configs)
    emit('statistics.json',dict(ensemble_phases=samples,window_fs=horizon,timestep_fs=dt,
        bath='discrete exact linear Bose OU bath, common complex innovations; Kerr dynamics has the same bath, not an exact nonlinear quantum canonical ensemble',
        event='discrete upward crossings at 2 fs, per D; no count is an independent lifetime or nuclear event',
        origin='single-arrival Palm response, phase-paired ambient subtraction; no accelerated flux',
        radial_energy_fraction_beyond_sampling_radius=float(np.exp(-5.)),
        thermal_energy_eV_per_polarization=float(np.sum(energies*thermal_n)/n),
        thermal_quanta_per_D=float(np.mean(thermal_n)),thermal_vacuum_used=False,
        total_probability_at_least_two='finite-memory Poisson probability only, not a bound on missing tail amplification or an inferred overlap initiation rate',
        sampling='12 paired phases; radial area midpoint quadrature in three explicit strata, not a spatial convergence claim',
        long_time='response window is four amplitude lifetimes; isolated linear residual energy fraction exp(-8)=0.00033546, nonlinear residual is unresolved rather than zero',
        support=records))
    np.savez_compressed(OUT/'maxima.npz',local_eV=max_local,closing_eV=max_close)
    return records


def figures():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def read(name):
        with (OUT/name).open() as f:return list(csv.DictReader(f))
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    disp=read('dispersion.csv')
    for case in dict.fromkeys(r['case'] for r in disp):
        a=[r for r in disp if r['case']==case]
        axes[0,0].plot([float(r['k_A_inverse']) for r in a],[1000*float(r['epsilon_eV']) for r in a],label=case)
    axes[0,0].set(xlabel='[110] wave number (1/Å)',ylabel='signed dispersion (meV)')
    axes[0,0].legend(fontsize=6)
    seed=read('seed.csv')
    axes[0,1].plot([float(r['time_fs'])/1000 for r in seed],[1000*float(r['max_local_energy_eV']) for r in seed])
    axes[0,1].set(xlabel='prepared seed time (ps)',ylabel='maximum local excitation (meV)')
    survival=read('survivals.csv')
    for label in ('linear','nonlinear'):
        for eta in (.001,.1):
            a=[r for r in survival if r['model']==label and float(r['E_cm_keV'])==2.5 and float(r['eta'])==eta and r['zone']=='0']
            axes[1,0].loglog([float(r['threshold_eV']) for r in a],[float(r['max_local_survival']) or np.nan for r in a],label=f'{label}, eta={eta}')
    axes[1,0].set(xlabel='energy threshold (eV)',ylabel='single-packet maximum survival')
    axes[1,0].legend(fontsize=7)
    rates=read('rates.csv')
    for label in ('linear','nonlinear'):
        a=[r for r in rates if r['model']==label and float(r['E_cm_keV'])==2.5 and float(r['eta'])==.1 and float(r['current_scale'])==1.]
        axes[1,1].semilogx([float(r['threshold_eV']) for r in a],[float(r['single_packet_pair_initiations_excess_per_D_s']) for r in a],label=label)
    axes[1,1].set(xlabel='closing-energy threshold (eV)',ylabel='Palm excess initiations / D / s')
    axes[1,1].legend()
    fig.tight_layout()
    fig.savefig(OUT/'waves.png',dpi=150)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spectrum',type=Path,required=True,help='Fresh contact spectrum bound to current bulk_pair.json')
    parser.add_argument('--reuse-solved',action='store_true',help='Reuse saved transitions and deterministic seed; resume the drive/consumer stage')
    args=parser.parse_args()
    OUT.mkdir(exist_ok=True)
    m,w,beam,models,mass=setup()
    if args.reuse_solved:
        with (OUT/'transitions.csv').open() as stream:
            levels = list(csv.DictReader(stream))
        e0 = float(levels[0]['next_gap_eV'])
        chi = float(levels[1]['next_gap_eV'])-e0
        t = HBAR**2/(2*mass)
        actions = None
        seed = json.loads((OUT/'seed.json').read_text())
        bind('g11_cmns/collective/focusing/pulse/replay.py')
        bind('g11_cmns/collective/contact/bulk_pair.json')
    else:
        e0,chi,t,levels,actions=transitions(w,mass)
        seed=prepared_seed(models[0],e0,chi,t,mass)
    loading=loading_and_growth(m,models,e0,chi,t,actions)
    with (OUT/'seed.csv').open() as stream:
        seed_rows = list(csv.DictReader(stream))
    continued = [float(row['time_fs']) for row in seed_rows if float(row['max_action'])>1.]
    seed['first_sampled_Kerr_continuation_fs'] = min(continued) if continued else None
    seed['continuation_time_resolution_fs'] = 10.
    emit('seed.json',seed)
    sources=drive_sources(beam,models[0])
    stats=statistics(models[0],e0,chi,t,mass,sources)
    figures()
    # The existing CLI is the actual final contact consumer. Its preparation,
    # support, joins, flags and energy/gauge checks are not replaced or relaxed.
    spectrum=args.spectrum if args.spectrum.is_absolute() else ROOT/args.spectrum
    bind(spectrum)
    command=[sys.executable,'g11_cmns/collective/focusing/pulse/replay.py',
        '--pair','g11_cmns/collective/contact/bulk_pair.json','--spectrum',str(spectrum.relative_to(ROOT)),
        '--histories','g13_rogue/output/histories.json','--trace-dir','g13_rogue',
        '--output','g13_rogue/output/contact.json']
    subprocess.run(command,cwd=ROOT,check=True,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
    bind('g11_cmns/collective/contact/contact.py')
    bind('g11_cmns/collective/contact/g12_readout.py')
    emit('inputs.json',dict(base=BASE_SHA,files=INPUTS,
        coefficients=dict(hbar_eV_fs=HBAR,mass_eV_fs2_A2=mass,kinetic_eV_A2=t,
                          discrete_gap_eV=e0,Kerr_eV_per_quantum=chi,amplitude_lifetime_fs=TAU),
        electron_density_term=dict(form='g_n * delta_n_e/n_e * psi',coefficient_eV=None,
            required='deformation potential g_n = d(hbar omega_D)/d ln(n_e), with phonon polarization, frequency, loading and ionic-strain boundary condition; and density ripple per acoustic strain',
            reason='g6_response supplies dielectric/drag response, not the derivative of the D vibrational Hessian or transition frequency with electron density'),
        energy_semantics='normal-ordered excitation; no zero-point energy harvested; source proxy from fixed stopping splits, not measured mode conversion',
        geometry_correction=dict(source='g11_cmns/collective/part1/stopping.json',
            fields=['terminal_depth_A'],actual_unit='m',converted_unit='nm',
            timeline_depth_A_unit='angstrom; retained as stored and divided by 10 for the agreement assertion',
            reason='terminal field alone stores solve_ivp SI metres; timeline was correctly converted; earlier results are not modified or assessed'),
        runtime=dict(command='python g13_rogue/run.py --spectrum g13_rogue/input/contact_spectrum.json'+(' --reuse-solved' if args.reuse_solved else ''),
            pulse_consumer='g11_cmns/collective/focusing/pulse/replay.py',exit='successful natural completion',
            preserved_producers='transitions and deterministic seed reused after the coordinate interpretation repair' if args.reuse_solved else None)))
    result=json.loads((OUT/'contact.json').read_text())
    probabilities=[b['log10_probability_proton_per_pair_pulse'] for p in result['pulses'] for v in p['variants'] for b in (v['lower_band'],v['upper_band'])]
    emit('result.json',dict(discrete_gap_eV=e0,nonlinearity_eV_per_quantum=chi,
        seed=seed,statistics=stats,contact_log10_proton_probability_range=[min(probabilities),max(probabilities)],
        scope='conditional normalized l=0 preparation of the modeled pumped seed episode, not an inferred beam reaction chance or a foil rate'))
    print(json.dumps(dict(output='g13_rogue/output/result.json',pulse_consumer_exit=0,
                         proton_probability_log10_range=[min(probabilities),max(probabilities)])))


if __name__=='__main__':
    main()
