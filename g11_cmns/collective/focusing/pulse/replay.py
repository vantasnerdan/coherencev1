"""Replay pair-bound inward/rebound histories in a declared D–D radial model."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import sys

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal

HERE = Path(__file__).resolve().parent
CONTACT = HERE.parent.parent / 'contact'
sys.path.insert(0, str(CONTACT))
from contact import BULK_PAIR, C, HBAR, T, SitePair, contact_log
from g12_readout import _LOG_D_PER_TRUE_C, _LOG_TRUE_CURRENT_A, _PROTON_A3_PER_S, _readout, event_readout

MASS = HBAR**2 / (2*T)   # relative mass in eV fs²/Å²
BULK_HISTORIES = HERE / 'bulk_histories.json'
ENERGY_REPLAY_ATOL_EV = 5e-9
RADIUS = 1.e-5   # Å, instantaneous inside 1 fm, not first passage
P_A = _PROTON_A3_PER_S
DD_A = P_A * (1 + 53.701/55.576)  # same spin convention, p + n branches only


def episode(trace, readout):
    """Nearest enclosing outer–inner–outer sign changes around largest c>0."""
    n = readout['scored_samples']
    exit_fs = readout['first_original_exit_fs']
    if len(trace) != n + int(exit_fs is not None) or (
            exit_fs is not None and not trace[n-1, 0] < exit_fs == trace[n, 0]):
        raise ValueError("Trace does not match the declared scored prefix and first exit")
    c = trace[:n, 9]
    outer = np.flatnonzero((c[:-1] <= 0) & (c[1:] > 0)) + 1
    inner = np.flatnonzero((c[:-1] > 0) & (c[1:] <= 0)) + 1
    peak = int(np.argmax(c))
    before = outer[outer <= peak]
    after = outer[outer > peak]
    minimum = inner[inner > peak]
    if not (c[peak] > 0 and len(before) and len(minimum) and len(after)
            and before[-1] < peak < minimum[0] < after[0]):
        raise ValueError("Expected a complete supported closing-and-rebound episode")
    start, turn, end = before[-1], minimum[0], after[0]
    seg = trace[start-1:end+1]
    t, radius, closing = seg[:, 0], seg[:, 3], seg[:, 9]
    path = CubicSpline(t, radius, bc_type=((1, -closing[0]), (1, -closing[-1])))
    roots = path.derivative().roots(extrapolate=False)
    turns = np.array([roots[(t[i-1] <= roots) & (roots <= t[i])][0]
                      for i in (1, turn-start+1, len(t)-1)])
    assert np.allclose(path(turns, 1), 0., atol=1e-12)
    sample_min = int(np.argmin(radius))
    assert start-1+sample_min in (turn-1, turn)
    outer_samples = (max(radius[0],radius[1]), max(radius[-2],radius[-1]))
    return path, turns, seg, dict(
        signed_turn_brackets_fs=[[float(trace[j-1,0]), float(trace[j,0])] for j in (start,turn,end)],
        interpolated_turns_fs=turns.tolist(), interpolated_turns_R_A=path(turns).tolist(),
        outer_turn_sampled_R_A=list(outer_samples), sampled_min_R_A=float(radius[sample_min]),
        sampled_min_R_fs=float(t[sample_min]), sampled_min_R_index=int(start-1+sample_min),
        first_receding_sample_index=int(turn), closing_peak_fs=float(trace[peak,0]),
        closing_peak_R_A=float(trace[peak,3]), closing_peak_A_fs=float(c[peak]),
        closing_peak_K_radial_eV=float(trace[peak,2]),
        max_c_spline_error_A_fs=float(np.max(np.abs(path(t,1)+closing))),
        first_original_exit_fs=exit_fs)


def drive_profile(r, flat_below, width):
    z = np.clip((r-flat_below)/width, 0, 1)
    return width*z**4*(z*z-3*z+2.5) + np.maximum(r-(flat_below+width), 0.)


def radial_basis(potential, archive, max_cap):
    """Only regenerate retained eigenvectors, verify energies and origin-slope gauge."""
    model = archive['model']
    archived_e = np.array(archive['spectrum']['energy_eV'])
    h = model['h_A']
    r = np.arange(1, model['n_grid']+1)*h
    count = int(np.searchsorted(archived_e, max_cap, side='right'))
    if count == 0 or max_cap > archived_e[-1]:
        raise ValueError("Pulse energy band is not covered by the supplied spectrum")
    en, vec = eigh_tridiagonal(2*T/h**2 + potential(r), np.full(len(r)-1,-T/h**2),
                               select='i', select_range=(0,count-1), tol=1e-11, lapack_driver='stebz')
    lcs = np.array([contact_log(potential, float(e), r, vec[:,i]/math.sqrt(h),
                               model['match_height_eV'], model['ode_rtol'])[:2] for i,e in enumerate(en)])
    vec *= lcs[:,1]
    log_rho = np.array(archive['spectrum']['log10_contact_A_minus3'][:count])
    gauge = dict(max_energy_difference_eV=float(np.max(np.abs(en-archived_e[:count]))),
                 max_log10_contact_difference=float(np.max(np.abs(lcs[:,0]/math.log(10)-log_rho))))
    # A distinct LAPACK subset can change the tiny resolved matching amplitude.
    # Magnitudes consumed by the replay are the original archived ones.
    if gauge['max_energy_difference_eV'] >= ENERGY_REPLAY_ATOL_EV or gauge['max_log10_contact_difference'] >= 2e-5:
        raise ArithmeticError((model['join_A'], model['include_flagged'], gauge))
    return r, archived_e[:count], vec, log_rho, gauge


def cusp_kernel(potential, energy, amp):
    """Gauss factors for sqrt(4π) amp_n r φ_n(r), regular at r=0."""
    nodes, weights = np.polynomial.legendre.leggauss(12)
    x = (nodes+1)/2
    a = potential.core_r[0]
    s = float(potential.electronic(a, 1))
    u0 = float(potential.electronic(a)) - s*a/2 + potential.shift
    coeff = [np.ones_like(energy)]
    for n in range(1,10):
        term = C*RADIUS/T*coeff[n-1]
        if n >= 2:
            term += (u0-energy)*RADIUS**2/T*coeff[n-2]
        if n >= 4:
            term += s/(2*a)*RADIUS**4/T*coeff[n-4]
        coeff.append(term/(n*(n+1)))
    phi = np.polynomial.polynomial.polyval(x, np.array(coeff)).T
    return np.sqrt(2*math.pi*RADIUS**3*weights)[:,None] * x[:,None]*phi*amp


def reaction_log(log_xi, coefficient):
    log_hazard = log_xi + math.log10(coefficient)
    # The correction to the log is beneath float resolution in this regime.
    return log_hazard if log_hazard < -16 else math.log10(-math.expm1(-10.**log_hazard))


def propagate(initial, eps, coupling, times, force):
    """Fourth-order commutator-free Magnus steps, earlier-weighted factor first."""
    width = np.diff(times)
    middle = times[:-1] + width/2
    offset = math.sqrt(3)/6 * width
    f1, f2 = force(middle-offset), force(middle+offset)
    a, b = (3+2*math.sqrt(3))/12, (3-2*math.sqrt(3))/12
    weights = np.stack((a*f1+b*f2, b*f1+a*f2), axis=1)
    energy, modes = np.linalg.eigh(.5*np.diag(eps)-weights[:,:,None,None]*coupling)
    phases = np.exp((-1j/HBAR*width[:,None,None])*energy)
    states = np.empty((len(times),len(initial)), dtype=complex)
    states[0] = z = initial
    for k,(columns,phase) in enumerate(zip(modes,phases),1):
        for q,p in zip(columns,phase):
            z = q@(p*(q.T@z))
        states[k] = z
    return states


def measures(path, turns, segment, peak, pot, spectrum, caps, sigma, drive, nlo=None, step_fs=.0625):
    r, energies, psi, log_rho, _ = spectrum
    count = int(np.searchsorted(energies, caps, side='right'))
    en, vec, lc = energies[:count], psi[:,:count], log_rho[:count]
    n_prepared = nlo or count
    t0, ti, tf = turns
    center = path(t0)
    launch = np.exp(-((r-center)/(2*sigma))**2) - np.exp(-((r+center)/(2*sigma))**2)
    launch /= np.linalg.norm(launch)
    initial = np.zeros(count, dtype=complex)
    initial[:n_prepared] = vec[:,:n_prepared].T@launch
    kept = float(np.vdot(initial,initial).real)
    initial /= math.sqrt(kept)
    R = vec.T@(r[:,None]*vec)
    D = vec.T@(drive_profile(r, *drive)[:,None]*vec)
    reference = float(np.max(lc))
    slopes = np.exp((lc-reference)*(.5*math.log(10)))
    gap = en-en[0]
    v0 = 1j/HBAR*(gap[:,None]-gap[None,:])*R
    v_drive = -1j/HBAR*(D@R-R@D)

    # The accepted splice is C0, not C1. On join crossings integrate from each
    # side separately, so the exact piecewise inverse force is not blurred.
    switches = path.solve(pot.join, extrapolate=False)
    switches = switches[(t0 < switches) & (switches < tf)]
    bounds = sorted(set([t0,peak,ti,tf,*switches]))
    state = initial
    time_parts, state_parts, force_parts = [], [], []
    closing_exposure = total = top_exposure = 0.
    for a,b in zip(bounds,bounds[1:]):
        derivative = pot.lattice if path((a+b)/2) >= pot.join else pot.inner
        times = np.linspace(a,b, math.ceil((b-a)/step_fs)+1)
        def force(t):
            return MASS * path(t,2) + derivative(path(t),1)
        local = propagate(state,gap,D,times,force)
        state = local[-1]
        area = float(simpson(abs(local@slopes)**2,x=times))
        total += area
        if b <= ti:
            closing_exposure += area
        if nlo is not None:
            top_exposure += float(simpson(abs(local[:,nlo:]@slopes[nlo:])**2,x=times))
        skip = int(bool(time_parts))
        time_parts.append(times[skip:])
        state_parts.append(local[skip:])
        force_parts.append(force(times[skip:]))
    tt, state, F = np.concatenate(time_parts), np.concatenate(state_parts), np.concatenate(force_parts)
    scales = np.einsum('ti,ti->t',state.conj(),state).real
    assert np.max(np.abs(scales-1.)) < 5e-10
    density = abs(state@slopes)**2
    def average(op):
        return np.einsum('ti,ti->t',state.conj(),state@op.T).real
    mean_r = average(R)
    current = average(v0)
    truncation_v = F*average(v_drive)
    mean_v = current + truncation_v
    kernel = cusp_kernel(pot,en,slopes)
    occupancy = np.sum(abs(state@kernel.T)**2,axis=1)
    kc, km, kp = int(np.argmax(density)), int(np.argmin(mean_r)), int(np.searchsorted(tt,peak))
    assert tt[kp] == peak
    ki = int(np.searchsorted(tt,ti))
    assert tt[ki] == ti
    assert 0 <= closing_exposure <= total
    assert 0 < total <= (tf-t0)*(slopes@slopes)*(1.+1e-10)
    log_xi = reference + math.log10(total) - 15.
    log_bound_rho = reference + math.log10(float(slopes@slopes))
    baseline = (tf-t0)*10.**(lc[0]-reference)
    excess = total-baseline
    log_excess = reference + math.log10(abs(excess)) - 15. if excess else None
    log_pp = reaction_log(log_xi,P_A)
    log_per_c = log_pp + _LOG_D_PER_TRUE_C
    def sample(k):
        return dict(time_fs=float(tt[k]), mean_R_A=float(mean_r[k]), imposed_R_A=float(path(tt[k])),
                    mean_velocity_A_fs=float(mean_v[k]), integrated_current_A_fs=float(current[k]),
                    imposed_velocity_A_fs=float(path(tt[k],1)))
    classical_e = segment[:,2] + pot(segment[:,3])
    result = dict(energy_cap_eV=caps, last_level_eV=float(en[-1]), levels=count,
                  classical_radial_energy_range_eV=[float(np.min(classical_e)),float(np.max(classical_e))],
                  gaussian_projection_norm2=kept, initial=sample(0), at_classical_closing_peak=sample(kp),
                  at_classical_inner_turn=sample(ki), at_outer_rebound=sample(len(tt)-1),
                  at_quantum_closest=sample(km), at_contact_peak=sample(kc),
                  rms_R_mismatch_A=float(np.sqrt(simpson((mean_r-path(tt))**2,x=tt)/(tf-t0))),
                  rms_velocity_mismatch_A_fs=float(np.sqrt(simpson((mean_v-path(tt,1))**2,x=tt)/(tf-t0))),
                  max_projected_velocity_correction_A_fs=float(np.max(abs(truncation_v))),
                  max_norm_error=float(np.max(np.abs(scales-1.))),
                  cusp_peak_log10_rho_A_minus3=reference+math.log10(float(density[kc])),
                  log10_max_instantaneous_probability_R_lt_1fm=reference+math.log10(float(occupancy.max())),
                  closing_fraction_of_Xi=closing_exposure/total,
                  log10_total_Xi_A_minus3_s=log_xi,
                  log10_undriven_ground_Xi_A_minus3_s=float(lc[0])+math.log10(tf-t0)-15.,
                  signed_excess_Xi=dict(sign=int(np.sign(excess)),log10_magnitude_A_minus3_s=log_excess),
                  log10_probability_proton_per_pair_pulse=log_pp,
                  log10_probability_either_DD_branch_per_pair_pulse=reaction_log(log_xi,DD_A),
                  one_pulse_per_incident_D_total_G12=(
                      _readout(log_per_c+_LOG_TRUE_CURRENT_A,log_per_c) if pot.pair.g12 else None),
                  one_pulse_per_incident_D_excess_G12=(
                      event_readout(log_excess) if excess > 0 and pot.pair.g12 else None),
                  finite_subspace_log10_rho_cap_A_minus3=log_bound_rho,
                  finite_subspace_log10_Xi_cap_A_minus3_s=log_bound_rho+math.log10(tf-t0)-15.,
                  finite_span_log10_snapshot_cap_R_lt_1fm=reference+2*math.log10(float(np.linalg.norm(kernel,ord=2))),
                  drive_join_crossings_fs=switches.tolist())
    if nlo is not None:
        result['upper_only_contact_control'] = dict(
            max_added_level_norm2=float(np.max(np.sum(abs(state[:,nlo:])**2,axis=1))),
            added_level_norm2_at_contact_peak=float(np.sum(abs(state[kc,nlo:])**2)),
            added_levels_alone_Xi_over_total_Xi=top_exposure/total)
    return result


def main():
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('--pair', type=Path, default=BULK_PAIR, help='Site-pair JSON used to solve the spectrum')
    argp.add_argument('--spectrum', type=Path, required=True, help='Spectrum output of contact.py for this pair')
    argp.add_argument('--histories', type=Path, default=BULK_HISTORIES, help='Manifest of pair-bound pulse traces')
    argp.add_argument('--trace-dir', type=Path, default=HERE.parent, help='Root for manifest-relative trace files')
    argp.add_argument('--output', type=Path, required=True, help='Write a new pulse result; archived results are not defaults')
    args = argp.parse_args()
    pair = SitePair(args.pair)
    quantum = math.sqrt(2*T*pair.outward_curvature)
    sigma = math.sqrt(T/quantum)
    profile = pair.pulse_drive
    if profile is None or not (profile['flat_below_A'] > 0 and profile['transition_width_A'] > 0):
        argp.error("Pulse input needs positive flat_below_A and transition_width_A")
    drive = profile['flat_below_A'], profile['transition_width_A']
    hashes = {'pair':pair.sha256}
    def bind(p, name):
        hashes[name] = sha256(p.read_bytes()).hexdigest()
        return p
    source = json.loads(bind(args.histories,'histories').read_text())
    archive = json.loads(bind(args.spectrum,'spectrum').read_text())
    if source['schema'] != 'g11-radial-histories-v1' or source['pair_sha256'] != pair.sha256:
        argp.error("Pulse histories are not bound to this site-pair description")
    if archive['inputs_sha256'] != pair.sha256:
        argp.error("Spectrum was not computed for this site-pair description")
    for name in ('contact.py','g12_readout.py'):
        bind(CONTACT/name,'../contact/'+name)
    rows = {(row['model']['join_A'],row['model']['include_flagged']):row for row in archive['variants']}
    potentials = pair.potentials()
    if len(rows) != len(potentials) or any((p.join, p.include_flagged) not in rows for p in potentials):
        argp.error("Spectrum variants do not match the pair joins and flag treatments")
    models = [(p, rows[p.join,p.include_flagged]) for p in potentials]
    histories = []
    for spec in source['histories']:
        name, key = spec['archived_trace'], spec['trace_key']
        path_to_trace = args.trace_dir/name
        if name not in hashes:
            bind(path_to_trace, name)
        with np.load(path_to_trace) as arrays:
            path, turns, segment, summary = episode(arrays[key], spec)
        if path(turns[1]) < sum(drive):
            argp.error(f"{spec['case']} reaches inside the force's unit-slope region")
        expected = spec.get('expected_inward_peak_eV')
        if expected is not None and not np.isclose(
                summary['closing_peak_K_radial_eV'], expected, atol=1e-12, rtol=0):
            raise ValueError(f"Wrong selected inward trace for {spec['case']} / {spec['launch']}")
        envelope = max(np.max(segment[:,2]+pot(segment[:,3])) for pot,_ in models)
        identity = dict(case=spec['case'], launch=spec['launch'], archived_trace=name, trace_key=key)
        histories.append((path,turns,segment,summary,envelope,identity,spec.get('half_step_control')))
    cut_max = max(row[4] for row in histories) + 4.5*quantum
    spectra = [radial_basis(pot,row,cut_max) for pot,row in models]
    pulses = []
    for path,turns,segment,episode_summary,envelope,identity,control in histories:
        variants = []
        for (pot,archived), spectrum in zip(models,spectra):
            peak = episode_summary['closing_peak_fs']
            lo = measures(path,turns,segment,peak,pot,spectrum,envelope+1.5*quantum,sigma,drive)
            hi = measures(path,turns,segment,peak,pot,spectrum,envelope+4.5*quantum,sigma,drive,nlo=lo['levels'])
            if control == {'join_A':pot.join, 'include_flagged':pot.include_flagged}:
                fine = measures(path,turns,segment,peak,pot,spectrum,envelope+4.5*quantum,sigma,drive,
                                nlo=lo['levels'],step_fs=.03125)
                hi['half_step_control'] = dict(
                    step_fs=.03125,
                    change_log10_Xi=fine['log10_total_Xi_A_minus3_s']-hi['log10_total_Xi_A_minus3_s'],
                    change_log10_peak_rho=fine['cusp_peak_log10_rho_A_minus3']-hi['cusp_peak_log10_rho_A_minus3'])
            variants.append(dict(join_A=pot.join, include_flagged=pot.include_flagged,
                                 lower_band=lo, upper_band=hi,
                                 upper_minus_lower_log10_Xi=hi['log10_total_Xi_A_minus3_s']-lo['log10_total_Xi_A_minus3_s'],
                                 archived_concentrated_event_G12=archived['preparation'].get('excess_event')))
        pulses.append(dict(**identity, **episode_summary, shared_radial_energy_envelope_eV=float(envelope),variants=variants))
    result = dict(schema='pdd-focusing-conditional-pulses-v1', pair_id=pair.id, input_sha256=hashes,
                  finite_model='normalized pure spatial l=0; no inferred material state or pair population',
                  preparation=dict(outward_local_curvature_quantum_eV=quantum, nominal_sigma_A=sigma,
                                   relative_mass_eV_fs2_A_minus2=MASS, near_origin_probability_radius_A=RADIUS,
                                   additive_energy_support_quanta=[1.5,4.5]),
                  spectrum_gauge=[dict(join_A=pot.join, include_flagged=pot.include_flagged, **spectrum[-1])
                                  for (pot,_),spectrum in zip(models,spectra)], pulses=pulses)
    # One compact row at a time (no lossy rounding of conserved inputs/results).
    header = json.dumps({k:v for k,v in result.items() if k != 'pulses'},indent=2,allow_nan=False)
    body = ',\n'.join('    '+json.dumps(row,allow_nan=False,separators=(',',':')) for row in pulses)
    args.output.write_text(header[:-2]+',\n  "pulses": [\n'+body+'\n  ]\n}\n')


if __name__ == '__main__':
    main()
