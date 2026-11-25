"""Actual-trajectory constrained launch optimization; no rate/contact consumer.

Units are eV, Angstrom and fs. Existing engine, inputs and accepted banks are read-only.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time

import numpy as np

import engine as E
import patterns as P
from large_geometry import LargeLattice

HERE = Path(__file__).resolve().parent
HORIZON = 2400.
TAU = 4100.
CASES = {
    'N2048_4.8meV': (8, .0048, .25, 'outside_soft_4.8meV'),
    'N2048_25meV': (8, .025, .25, 'fine_N2048_25meV'),
    'N16384_0.6meV': (16, .0006, .5, 'size_N16384'),
}
OBJECTIVES = ('K', 'inward', 'distance')


def digest(path):
    content = Path(path).read_bytes()
    return dict(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def seed(lattice, key):
    """Recover the accepted launch without rerunning the accepted candidate bank."""
    data = json.loads((HERE/'continuation_results.json').read_text())
    source = CASES[key][3]
    case = data['cases'][source]
    selected = next(r for r in case['selected'] if r['pattern']=='direct_optimized')
    with np.load(HERE/'continuation_states.npz') as states:
        saved = 'outside_soft_'+('4.8' if '4.8' in key else '25')+'meV__launch_v'
        if lattice.n == 2048:
            v = states[saved].copy()
            lineage = dict(source='continuation_states.npz', key=saved)
        else:
            speed = np.sqrt(2*case['energy_per_D_eV']/E.MASS)
            bases = list(P.velocity_candidates('direct', lattice.positions, lattice.box,
                lattice.target, lattice.target_axis, speed, lattice.predicted_kernel,
                focus_times_fs=case['proposed_focus_times_fs']))
            rows = {r['candidate']: r for r in case['candidates']}
            best = max(bases, key=lambda b: rows[b['name']]['peak_relative_kinetic_eV'])
            other = next(b for b in bases if b is not best)
            v = P.capped_direct_velocity(best['v']-.25*other['v'], lattice.target,
                                         speed, lattice.target_axis)
            lineage = dict(source='exact deterministic stored proposals and stored base rankings',
                best=best['name'], other=other['name'], mixture=-.25,
                proposal_times_fs=case['proposed_focus_times_fs'])
    return v, dict(accepted_case=source, accepted_candidate=selected['candidate'],
                   accepted_row=selected, reconstruction=lineage)


def objective_values(trace, objective):
    if objective == 'K':
        return trace[:, 1]
    if objective == 'inward':
        return .25*E.MASS*np.maximum(trace[:, 9], 0.)**2
    if objective == 'distance':
        return np.linalg.norm(E.BASE['positions'][E.TJ]-E.BASE['positions'][E.TI])-trace[:, 3]
    raise ValueError(objective)


def joint(trace, k, objective):
    values = objective_values(trace, objective)
    lo = hi = int(k)
    threshold = values[k]/2
    while lo and values[lo-1] >= threshold:
        lo -= 1
    while hi+1 < len(trace) and values[hi+1] >= threshold:
        hi += 1
    row = trace[k]
    return dict(step=int(k), time_fs=float(row[0]), Krel_eV=float(row[1]),
        inward_K_eV=float(.25*E.MASS*max(0., row[9])**2), Krad_eV=float(row[2]),
        signed_inward_A_fs=float(row[9]), receding=bool(row[9] < 0.), R_A=float(row[3]),
        max_site_radius_A=float(row[4]), min_NN_R_A=float(row[5]),
        total_energy_eV=float(row[6]), objective_value=float(values[k]),
        half_objective_hold_fs=float(trace[hi, 0]-trace[lo, 0]),
        hold_left_censored=bool(lo == 0), hold_right_censored=bool(hi == len(trace)-1))


def summarize(result, objective, dt):
    trace, exit_time = result[2], result[3]
    supported = trace if exit_time < 0 else trace[:-1]
    events = {o: joint(supported, int(np.argmax(objective_values(supported, o))), o)
              for o in OBJECTIVES}
    peak = events[objective]
    return dict(objective=objective, score=peak['objective_value'], selected_event=peak,
        joint_events=events, first_original_exit_fs=None if exit_time < 0 else float(exit_time),
        exit_kind={0: None, 1: 'site_radius', 2: 'bond_distance'}[result[4]],
        full_horizon_supported=bool(exit_time < 0 and trace[-1, 0] == HORIZON),
        scored_through_fs=float(supported[-1, 0]), scored_samples=len(supported),
        max_site_radius_A=float(supported[:, 4].max()), min_NN_R_A=float(supported[:, 5].min()),
        dt_fs=dt, support_rule='strictly before first original all-step crossing')


def constraints(lattice, v, energy):
    norm2 = np.einsum('ij,ij->i', v, v)
    expected = 2*lattice.n*energy/E.MASS
    cap = 2*np.sqrt(lattice.n/(lattice.n-2))*np.sqrt(2*energy/E.MASS)
    return dict(initial_total_Ne_eV=lattice.n*energy,
        actual_initial_kinetic_eV=float(.5*E.MASS*norm2.sum()),
        squared_norm_relative_error=float(abs(norm2.sum()-expected)/expected),
        max_speed_A_fs=float(np.sqrt(norm2.max())), cap_A_fs=float(cap),
        target_velocity_max_A_fs=float(abs(v[list(lattice.target)]).max()),
        initial_displacement_max_A=0., initial_potential_eV=0.,
        kinetic_IPR=float(norm2.sum()**2/(norm2@norm2)),
        minimum_IPR_from_cap=(lattice.n-2)/4,
        valid=bool(abs(norm2.sum()-expected) <= 2e-12*expected and
                   np.sqrt(norm2.max()) <= cap*(1+2e-12) and
                   not np.any(v[list(lattice.target)])))


def forward(lattice, v, dt, tau=0.):
    return E.trajectory(lattice, np.zeros_like(v), v, dt, HORIZON, damping=tau)


def trial_forward(lattice, v, dt, tau, target_step, objective, incumbent):
    """Prefilter the fixed-time branch, then finish every improving trajectory."""
    prefix = E.trajectory(lattice, np.zeros_like(v), v, dt, target_step*dt, damping=tau)
    scored = summarize(prefix, objective, dt)
    if prefix[3] >= 0 or scored['score'] <= incumbent*(1+1e-6):
        return prefix, scored, False
    remaining = int(round(HORIZON/dt))-target_step
    if not remaining:
        return prefix, scored, True
    tail = E.trajectory(lattice, prefix[0], prefix[1], dt, remaining*dt, damping=tau)
    tail[2][:, 0] += target_step*dt
    trace = np.concatenate((prefix[2][:-1], tail[2]))
    result = (tail[0], tail[1], trace,
              tail[3]+target_step*dt if tail[3] >= 0 else -1., tail[4])
    return result, summarize(result, objective, dt), True


def optimize(lattice, v0, energy, dt, objective, tau, prefix, budget_seconds=6000.,
             initial_result=None, prefilter=True):
    """Projected ascent with actual full-horizon pre-exit acceptance and retargeting.

    A nonzero first-order projected direction is retried down to an RMS launch
    displacement of 1e-5. Stopping is measured local stall/plateau, not the cap.
    """
    from adjoint import endpoint_gradient
    started = time.monotonic()
    v = v0.copy()
    result = forward(lattice, v, dt, tau) if initial_result is None else initial_result
    current = summarize(result, objective, dt)
    initial = current.copy()
    accepted = [dict(iteration=0, source='declared seed', **current)]
    trials, costs = [], []
    step_fraction = .03
    stop = 'unresolved_resource_cap'
    stagnant = 0
    gradient_count = 0
    for iteration in range(1, 61):
        if time.monotonic()-started > budget_seconds:
            break
        target_step = current['selected_event']['step']
        tick = time.monotonic()
        value, g, xend, vend = endpoint_gradient(lattice, v, dt, target_step,
            objective, tau_fs=tau, exterior='soft')
        gradient_count += 1
        costs.append(dict(kind='gradient', seconds=time.monotonic()-tick,
                          step=target_step, endpoint_value=float(value)))
        if not np.all(np.isfinite(g)):
            stop = 'unresolved_nonfinite_gradient'
            break
        g[list(lattice.target)] = 0.
        # Remove the sphere normal; cap projection handles the remaining active set.
        g -= v*(np.sum(v*g)/np.sum(v*v))
        gnorm = float(np.linalg.norm(g))
        if gnorm == 0.:
            stop = 'zero_projected_first_order_direction'
            break
        direction = g*(np.linalg.norm(v)/gnorm)
        accepted_trial = None
        alpha = min(.12, step_fraction*1.5)
        scale = max(abs(current['score']), 1e-12)
        timed_out = False
        while alpha >= 1e-5:
            if time.monotonic()-started > budget_seconds:
                timed_out = True
                break
            proposed = P.capped_direct_velocity(v+alpha*direction, lattice.target,
                np.sqrt(2*energy/E.MASS), lattice.target_axis)
            check = constraints(lattice, proposed, energy)
            if not check['valid']:
                raise RuntimeError('projection violated launch contract')
            displacement = float(np.linalg.norm(proposed-v)/np.linalg.norm(v))
            tick = time.monotonic()
            if prefilter:
                candidate_result, candidate, scanned = trial_forward(
                    lattice, proposed, dt, tau, target_step, objective, current['score'])
            else:
                candidate_result = forward(lattice, proposed, dt, tau)
                candidate, scanned = summarize(candidate_result, objective, dt), True
            relative_gain = (candidate['score']-current['score'])/scale
            trial = dict(iteration=iteration, alpha=alpha, relative_launch_step=displacement,
                relative_gain=relative_gain, seconds=time.monotonic()-tick,
                endpoint_target_step=target_step, full_horizon_scanned=scanned,
                projected_linear_gain=float(np.sum(g*(proposed-v))), **candidate)
            trials.append(trial)
            # Never call a lower supported-prefix peak an improvement to a longer trace.
            if relative_gain > 1e-6 and (scanned or candidate_result[3] >= 0):
                accepted_trial = (proposed, candidate_result, candidate, relative_gain, alpha)
                break
            alpha *= .5
        if accepted_trial is None:
            stop = ('unresolved_resource_cap' if timed_out else
                    'feasible_fixed_time_prefilter_stall_to_1e-5_RMS')
            break
        v, result, current, gain, alpha = accepted_trial
        accepted.append(dict(iteration=iteration, alpha=alpha, relative_gain=gain,
                             source='actual nonlinear gradient and full-horizon acceptance', **current))
        step_fraction = alpha
        stagnant = stagnant+1 if gain < .002 else 0
        np.savez_compressed(HERE/(prefix+'_checkpoint.npz'), v=v, trace=result[2])
        (HERE/(prefix+'_checkpoint.json')).write_text(json.dumps(dict(
            accepted=accepted, trials=trials, costs=costs), indent=2, allow_nan=False)+'\n')
        print(json.dumps(dict(run=prefix, iteration=iteration, score=current['score'],
            gain=gain, time_fs=current['selected_event']['time_fs'],
            exit_fs=current['first_original_exit_fs']), allow_nan=False), flush=True)
        if stagnant >= 3:
            # One explicitly larger feasible trial distinguishes a small-step artifact.
            proposed = P.capped_direct_velocity(v+min(.12, 8*alpha)*direction,
                lattice.target, np.sqrt(2*energy/E.MASS), lattice.target_axis)
            tick = time.monotonic()
            candidate_result = forward(lattice, proposed, dt, tau)
            candidate = summarize(candidate_result, objective, dt)
            gain = (candidate['score']-current['score'])/max(abs(current['score']), 1e-12)
            trials.append(dict(iteration=iteration, alpha=min(.12, 8*alpha),
                plateau_challenge=True, relative_gain=gain, seconds=time.monotonic()-tick, **candidate))
            if gain > .002:
                v, result, current = proposed, candidate_result, candidate
                accepted.append(dict(iteration=iteration, relative_gain=gain,
                    source='larger-step plateau challenge', **current))
                stagnant = 0
                step_fraction = min(.12, 8*alpha)
            else:
                stop = 'relative_gain_plateau_below_0.2percent_with_larger_step_challenge'
                break
    else:
        stop = 'unresolved_iteration_safety_cap'
    return v, result, dict(initial=initial, final=current, accepted_history=accepted,
        feasible_trials=trials, costs=costs, stop_reason=stop,
        resolved_local_stop=not stop.startswith('unresolved'), gradient_count=gradient_count,
        runtime_seconds=time.monotonic()-started, launch_constraints=constraints(lattice,v,energy),
        method='fixed-step discrete nonlinear gradient; retarget actual all-step objective peak after acceptance',
        global_optimality=False)


def run_case(key, budget_seconds, prefilter=True):
    m, energy, dt, source = CASES[key]
    lattice = LargeLattice(m)
    v0, lineage = seed(lattice, key)
    old_row = lineage['accepted_row']
    tick = time.monotonic()
    original = forward(lattice, v0, dt)
    # Deterministic original launch reproduction is also the first required scoring pass.
    old_peak = old_row['peak_relative_kinetic_eV']
    reproduced = summarize(original, 'K', dt)
    if abs(reproduced['score']-old_peak) > 2e-9*max(1., old_peak):
        raise RuntimeError('accepted seed did not reproduce accepted selected K peak')
    report = dict(schema='pdd-nonlinear-search-v1', case=key, n_D=lattice.n,
        energy_per_D_eV=energy, dt_fs=dt, horizon_fs=HORIZON,
        fixed_loss_amplitude_tau_fs=TAU, lineage=lineage,
        optimization_settings=dict(endpoint_prefilter=prefilter,
            initial_RMS_step=.03, max_RMS_step=.12, minimum_RMS_step=1e-5,
            accepted_relative_gain_floor=1e-6, meaningful_relative_gain=.002,
            plateau_successes=3, plateau_challenge_multiplier=8,
            per_objective_budget_seconds=budget_seconds/4, iteration_safety_cap=60,
            numerical_threads=1),
        original_baseline={o: summarize(original,o,dt) for o in OBJECTIVES},
        seed_roundtrip_peak_abs_error_eV=abs(reproduced['score']-old_peak),
        seed_roundtrip_seconds=time.monotonic()-tick,
        initial_constraints=constraints(lattice,v0,energy), searches={},
        inputs={n: digest(HERE/n) for n in ('engine.py','patterns.py','exterior.py',
            'large_geometry.py','continuation_results.json',
            'continuation_states.npz','input_manifest.json')})
    arrays = {'original_v':v0, 'original_trace':original[2]}
    per_search = budget_seconds/4
    for objective in OBJECTIVES:
        prefix = 'nonlinear_'+key+'_'+objective
        v,result,record = optimize(lattice,v0,energy,dt,objective,0.,prefix,per_search,
                                 initial_result=original, prefilter=prefilter)
        record['seed_lineage'] = 'accepted selected eleven-candidate original launch'
        report['searches'][objective] = record
        arrays[objective+'_v'],arrays[objective+'_trace'] = v,result[2]
        arrays[objective+'_final_x'], arrays[objective+'_final_v'] = result[:2]
        save_case(key,report,arrays)
    original_loss = forward(lattice,v0,dt,TAU)
    conservative_loss = forward(lattice,arrays['K_v'],dt,TAU)
    report['loss_comparison'] = dict(
        original_seed_under_loss=summarize(original_loss,'K',dt),
        conservative_optimum_under_actual_loss=summarize(conservative_loss,'K',dt))
    arrays['original_loss_trace'] = original_loss[2]
    arrays['conservative_under_loss_trace'] = conservative_loss[2]
    # Reoptimize from the better feasible dynamic-loss launch, not post-hoc scaling.
    loss_seed = arrays['K_v'] if summarize(conservative_loss,'K',dt)['score'] >= summarize(original_loss,'K',dt)['score'] else v0
    report['loss_comparison']['loss_search_seed'] = ('conservative_optimum' if loss_seed is not v0 else 'original')
    # The rejected post-hoc exponential is a distinct comparison, never the dynamics.
    conservative_trace = arrays['K_trace'].copy()
    conservative_trace[:,1:3] *= np.exp(-2*conservative_trace[:,0,None]/TAU)
    conservative_trace[:,9] *= np.exp(-conservative_trace[:,0]/TAU)
    report['loss_comparison']['posthoc_damped_conservative_trace'] = joint(
        conservative_trace,int(np.argmax(conservative_trace[:,1])),'K')
    v,result,record = optimize(lattice,loss_seed,energy,dt,'K',TAU,
        'nonlinear_'+key+'_loss',per_search,
        initial_result=conservative_loss if loss_seed is not v0 else original_loss,
        prefilter=prefilter)
    record['seed_lineage'] = report['loss_comparison']['loss_search_seed']
    report['searches']['loss_K'] = record
    arrays['loss_K_v'],arrays['loss_K_trace'] = v,result[2]
    arrays['loss_K_final_x'], arrays['loss_K_final_v'] = result[:2]
    report['loss_comparison']['reoptimized_under_loss'] = record['final']
    save_case(key,report,arrays)


def save_case(key, report, arrays):
    (HERE/('nonlinear_'+key+'.json')).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(HERE/('nonlinear_'+key+'.npz'),**arrays)


def probe():
    from adjoint import endpoint_gradient
    key = 'N16384_0.6meV'
    lattice = LargeLattice(16)
    v,_ = seed(lattice,key)
    # Compile on a short endpoint, then measure one representative warmed gradient.
    endpoint_gradient(lattice,v,.5,2,'K')
    tick = time.monotonic()
    value,g,x,vend = endpoint_gradient(lattice,v,.5,1600,'K')
    elapsed = time.monotonic()-tick
    tick = time.monotonic()
    result = forward(lattice,v,.5)
    forward_seconds = time.monotonic()-tick
    record = dict(gradient_step=1600,gradient_physical_time_fs=800.,
        warmed_gradient_seconds=elapsed,full_horizon_forward_seconds=forward_seconds,
        estimated_1799fs_gradient_seconds=elapsed*1799/800,
        value=value,gradient_norm=float(np.linalg.norm(g)),
        baseline_peak=summarize(result,'K',.5)['score'],
        numerical_threads=1, parallel_jobs_limit=2, per_job_memory_limit_bytes=2147483648)
    (HERE/'nonlinear_probe.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2),flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('case',choices=(*CASES,'probe'))
    parser.add_argument('--budget-seconds',type=float,default=6000.)
    parser.add_argument('--no-prefilter',action='store_true')
    args = parser.parse_args()
    if args.case == 'probe':
        probe()
    else:
        run_case(args.case,args.budget_seconds,not args.no_prefilter)


if __name__ == '__main__':
    main()
