"""Exact-seed, actual-motion inward continuation; eV, Angstrom, fs.

The accepted banks and stationary contact/rate tables are read-only.  Each
accepted launch is scored through the full horizon or strictly before exit.
"""
from pathlib import Path
import argparse
import json
import time

import numpy as np

import engine as E
import nonlinear as N
import patterns as P
from adjoint import endpoint_gradient
from large_geometry import LargeLattice
from nonlinear_select import selected_result

HERE = Path(__file__).resolve().parent
ENERGY, DT, HORIZON, TAU = .0006, .5, 2400., 4100.


def cap_diagnostics(lattice, launch):
    record = N.constraints(lattice, launch, ENERGY)
    ratios = np.linalg.norm(launch, axis=1)/record['cap_A_fs']
    record.update(max_speed_over_cap=float(ratios.max()),
                  cap_active_sites=int(np.count_nonzero(abs(ratios-1.) < 1e-9)))
    return record


def restore(result_arrays, prefix, readout):
    exit_fs = readout['first_original_exit_fs']
    return (result_arrays[prefix+'_final_x'].copy(),
            result_arrays[prefix+'_final_v'].copy(),
            result_arrays[prefix+'_trace'].copy(),
            -1. if exit_fs is None else exit_fs,
            {None: 0, 'site_radius': 1, 'bond_distance': 2}[readout['exit_kind']])


def save(size, report, arrays):
    stem = HERE/('inward_N'+str(size))
    np.savez_compressed(str(stem)+'.npz', **arrays)
    Path(str(stem)+'.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


def store(arrays, prefix, launch, result):
    arrays[prefix+'_launch_v'] = launch.copy()
    for suffix, value in zip(('_final_x', '_final_v', '_trace'), result[:3]):
        arrays[prefix+suffix] = value.copy()


def exact_seed(lattice):
    if lattice.n == 16384:
        key = 'N16384_0.6meV'
        selected = json.loads((HERE/'nonlinear_selection.json').read_text())[key]['selected']['inward']
        with np.load(HERE/('nonlinear_'+key+'.npz')) as raw, np.load(HERE/'nonlinear_selection.npz') as bound:
            launch, result = selected_result(key, 'inward', selected, raw, bound)
            launch, result = launch.copy(), tuple(a.copy() if isinstance(a, np.ndarray) else a for a in result)
        return launch, result, dict(source='exact previous selected inward launch',
            selection=selected, original_total_criterion_eV=.12483784143540816)
    data = json.loads((HERE/'continuation_results.json').read_text())['cases']['size_N131072']
    bases = list(P.velocity_candidates('direct', lattice.positions, lattice.box,
        lattice.target, lattice.target_axis, np.sqrt(2*ENERGY/E.MASS),
        lattice.predicted_kernel, focus_times_fs=data['proposed_focus_times_fs']))
    rows = {row['candidate']: row for row in data['candidates']}
    best = max(bases, key=lambda base: rows[base['name']]['peak_relative_kinetic_eV'])
    other = next(base for base in bases if base is not best)
    chosen = next(row for row in data['selected'] if row['pattern'] == 'direct_optimized')
    assert chosen['candidate'] == 'direct_mixture_-0.25'
    assert data['proposed_focus_times_fs'] == [1940., 2100., 2260.]
    launch = P.capped_direct_velocity(best['v']-.25*other['v'], lattice.target,
        np.sqrt(2*ENERGY/E.MASS), lattice.target_axis)
    with np.load(HERE/'continuation_traces.npz') as saved:
        expected = saved['size_N131072__direct_mixture_-0.25'].copy()
    # One actual matched initial trajectory, not the old proposal matrix.
    tick = time.monotonic()
    result = N.forward(lattice, launch, DT)
    defect = float(np.max(abs(result[2]-expected)))
    if defect > 2e-9:
        raise RuntimeError('exact original launch/phase restoration failed')
    return launch, result, dict(source='stored three base rankings; continuation.py:128-142',
        best=best['name'], other=other['name'], mixture=-.25,
        stored_proposal_times_fs=data['proposed_focus_times_fs'], selected_candidate=chosen['candidate'],
        restored_trace_max_abs_error=defect, restoration_seconds=time.monotonic()-tick,
        original_total_criterion_eV=chosen['peak_relative_kinetic_eV'])


def optimize(lattice, initial_v, initial_result, tau, budget, checkpoint, check_gradient=False,
             previous=None, max_updates=60, max_trials=20):
    v, result = initial_v.copy(), initial_result
    current = N.summarize(result, 'inward', DT)
    record = previous if previous is not None else dict(
        initial=current, accepted_history=[dict(iteration=0, **current)],
        feasible_trials=[], costs=[], global_optimality=False,
        wall_clock_started_s=time.time())
    started = time.monotonic()-(time.time()-record['wall_clock_started_s'])
    history = record['accepted_history']
    alpha_last = history[-1].get('alpha', .03)
    gains = [row.get('relative_gain', 1.) for row in history[-3:]]
    stagnant = sum(gain < .002 for gain in gains) if all(gain < .002 for gain in gains) else 0
    gradient_cost = record['costs'][-1]['seconds'] if record['costs'] else None
    trial_cost = record['feasible_trials'][-1]['seconds'] if record['feasible_trials'] else None
    stop = 'out_of_proportion_to_further_local_improvement'
    first_iteration = history[-1]['iteration']+1
    for iteration in range(first_iteration, first_iteration+max_updates):
        elapsed = time.monotonic()-started
        # Do not start an update whose measured gradient plus two native trials
        # would exceed this bounded decision.  The first gradient prices itself.
        estimate = 0. if gradient_cost is None else 1.15*gradient_cost+2*trial_cost
        if elapsed+estimate > budget:
            break
        step = current['selected_event']['step']
        tick = time.monotonic()
        value, gradient, xend, vend = endpoint_gradient(lattice, v, DT, step, 'inward', tau_fs=tau)
        gradient_cost = time.monotonic()-tick
        record['costs'].append(dict(kind='gradient', seconds=gradient_cost, step=step,
                                    endpoint_value=float(value)))
        if not np.isfinite(gradient).all():
            raise RuntimeError('nonfinite actual-motion gradient')
        if abs(value-current['score']) > 2e-8*max(1., abs(value)):
            raise RuntimeError('adjoint endpoint differs from actual strongest closing sample')
        if check_gradient:
            endpoint = E.trajectory(lattice, xend, vend, DT, 0., damping=tau)
            diagnostic = endpoint[2][0].copy()
            diagnostic[0] = step*DT
            defect = float(np.max(abs(diagnostic[:10]-result[2][step, :10])))
            record['endpoint_support_witness'] = dict(step=step, time_fs=step*DT,
                native_trace_max_abs_error=defect, native_endpoint_exit_fs=float(endpoint[3]),
                expected_simultaneous_sample=current['selected_event'],
                passed=bool(defect < 2e-8 and endpoint[3] < 0))
            if not record['endpoint_support_witness']['passed']:
                raise RuntimeError('gradient endpoint/support witness failed')
        record.update(final=current, launch_constraints=cap_diagnostics(lattice, v),
                      stop_reason='search_in_progress', elapsed_seconds=time.monotonic()-started)
        checkpoint(v, result, record, gradient_state=(gradient, xend, vend))
        gradient[list(lattice.target)] = 0.
        gradient -= v*(np.sum(v*gradient)/np.sum(v*v))
        norm = np.linalg.norm(gradient)
        if norm == 0.:
            stop = 'zero_projected_first_order_direction'
            break
        direction = gradient*(np.linalg.norm(v)/norm)
        if check_gradient and iteration == 1:
            h = 1e-6
            perturbed = P.capped_direct_velocity(v+h*direction, lattice.target,
                np.sqrt(2*ENERGY/E.MASS), lattice.target_axis)
            tick = time.monotonic()
            native = E.trajectory(lattice, np.zeros_like(v), perturbed, DT, step*DT, damping=tau)
            actual = float(N.objective_values(native[2], 'inward')[-1]-value)
            predicted = float(np.sum(gradient*(perturbed-v)))
            defect = abs(actual-predicted)/max(abs(actual), abs(predicted), 1e-30)
            witness = dict(relative_launch_step=h, actual_fixed_time_gain_eV=actual,
                predicted_gain_eV=predicted, relative_directional_defect=defect,
                closing_sign=float(native[2][-1, 9]), original_prefix_exit_fs=float(native[3]),
                seconds=time.monotonic()-tick, passed=bool(defect < .02 and actual > 0 and native[3] < 0))
            record['directional_witness'] = witness
            if not witness['passed']:
                raise RuntimeError('131k actual inward gradient witness failed: '+str(witness))
        record.update(final=current, launch_constraints=cap_diagnostics(lattice, v),
                      stop_reason='search_in_progress', elapsed_seconds=time.monotonic()-started)
        checkpoint(v, result, record)
        # Price the *feasible* local slope before a costly native line search;
        # a large fixed RMS first step can leave this narrow phase basin.
        local_trial = P.capped_direct_velocity(v+1e-6*direction, lattice.target,
            np.sqrt(2*ENERGY/E.MASS), lattice.target_axis)
        local_slope = float(np.sum(gradient*(local_trial-v))/1e-6)
        accepted = None
        alpha = min(.12, 1.5*alpha_last,
                    max(1e-5, .05*current['score']/max(local_slope, 1e-30)))
        attempts = 0
        while alpha >= 1e-5 and attempts < max_trials:
            if time.monotonic()-started > budget:
                break
            proposal = P.capped_direct_velocity(v+alpha*direction, lattice.target,
                np.sqrt(2*ENERGY/E.MASS), lattice.target_axis)
            assert cap_diagnostics(lattice, proposal)['valid']
            attempts += 1
            tick = time.monotonic()
            candidate_result, candidate, scanned = N.trial_forward(
                lattice, proposal, DT, tau, step, 'inward', current['score'])
            scanned = scanned or candidate_result[3] >= 0  # complete pre-first-exit path
            trial_cost = time.monotonic()-tick
            gain = candidate['score']/current['score']-1.
            row = dict(iteration=iteration, alpha=alpha, seconds=trial_cost,
                relative_gain=gain, full_horizon_scanned=scanned,
                relative_launch_step=float(np.linalg.norm(proposal-v)/np.linalg.norm(v)), **candidate)
            record['feasible_trials'].append(row)
            if scanned and gain > 1e-6:
                accepted = proposal, candidate_result, candidate, gain, alpha
                break
            alpha *= .5
        if accepted is None:
            if alpha < 1e-5:
                stop = 'feasible_fixed_time_backtracking_stall_to_1e-5_RMS'
            if max_updates == 1:
                stop = 'bounded_proposals_no_accepted_improvement_no_plateau'
            break
        v, result, current, gain, alpha_last = accepted
        record['accepted_history'].append(dict(iteration=iteration, relative_gain=gain,
                                              alpha=alpha_last, **current))
        stagnant = stagnant+1 if gain < .002 else 0
        record.update(final=current, launch_constraints=cap_diagnostics(lattice, v),
                      stop_reason='search_in_progress', elapsed_seconds=time.monotonic()-started)
        checkpoint(v, result, record)
        print(json.dumps(dict(n_D=lattice.n, tau_fs=tau, iteration=iteration,
            inward_eV=current['score'], gain=gain, time_fs=current['selected_event']['time_fs'],
            first_exit_fs=current['first_original_exit_fs']), allow_nan=False), flush=True)
        if max_updates == 1:
            stop = 'bounded_actual_update_completed_no_plateau'
            break
        if stagnant >= 3:
            # Challenge the small-step explanation using the same still-useful
            # direction, completing the full trajectory before any acceptance.
            proposal = P.capped_direct_velocity(v+min(.12, 8*alpha_last)*direction,
                lattice.target, np.sqrt(2*ENERGY/E.MASS), lattice.target_axis)
            tick = time.monotonic()
            challenge_result = N.forward(lattice, proposal, DT, tau)
            challenge = N.summarize(challenge_result, 'inward', DT)
            gain = challenge['score']/current['score']-1.
            record['feasible_trials'].append(dict(iteration=iteration, plateau_challenge=True,
                alpha=min(.12, 8*alpha_last), seconds=time.monotonic()-tick,
                relative_gain=gain, full_horizon_scanned=True, **challenge))
            if gain > 1e-6:
                v, result, current = proposal, challenge_result, challenge
                record['accepted_history'].append(dict(iteration=iteration,
                    plateau_challenge=True, relative_gain=gain, **current))
            if gain <= .002:
                stop = 'local_gain_plateau_below_0.2percent_with_larger_step_challenge'
                break
            stagnant = 0
            alpha_last = min(.12, 8*alpha_last)
    record.update(final=current, launch_constraints=cap_diagnostics(lattice, v),
        stop_reason=stop, resolved_local_stop=not stop.startswith(('out_of_proportion', 'bounded')),
        elapsed_seconds=time.monotonic()-started,
        method='actual discrete inward adjoint; projected line search; prefilter then full pre-first-exit acceptance')
    return v, result, record


def run(size, stage, budget, resume=False, max_updates=60, max_trials=20):
    lattice = LargeLattice({16384: 16, 131072: 32}[size])
    path = HERE/('inward_N'+str(size)+'.json')
    if path.exists():
        report = json.loads(path.read_text())
        with np.load(path.with_suffix('.npz')) as saved:
            arrays = {k: saved[k].copy() for k in saved.files}
    else:
        launch, result, lineage = exact_seed(lattice)
        report = dict(schema='pdd-inward-continuation-v1', n_D=size, energy_per_D_eV=ENERGY,
            dt_fs=DT, horizon_fs=HORIZON, fixed_loss_amplitude_tau_fs=TAU, lineage=lineage,
            original=N.summarize(result, 'inward', DT), original_constraints=cap_diagnostics(lattice, launch),
            settings=dict(minimum_RMS_step=1e-5, meaningful_gain=.002,
                          plateau_successes=3, budget_seconds=budget, numerical_threads=1),
            inputs={name: N.digest(HERE/name) for name in
                ('engine.py', 'patterns.py', 'exterior.py', 'large_geometry.py', 'adjoint.py', 'nonlinear.py')})
        arrays = {}
        store(arrays, 'original', launch, result)
        save(size, report, arrays)
    prefix = stage
    previous = report.get(stage)
    if resume and previous is None:
        report['interruption_recovery'] = (
            'saved exact original reused after an interrupted unserialized attempt; '
            'no previous accepted improvement, rejected-trial result or plateau inferred')
    if previous is not None:
        if not resume or previous['stop_reason'] != 'search_in_progress':
            raise RuntimeError('stage already exists; only an interrupted checkpoint may be resumed explicitly')
        launch = arrays[stage+'_launch_v']
        result = restore(arrays, stage, previous['final'])
        tau = 0. if stage == 'conservative' else TAU
    elif stage == 'conservative':
        launch = arrays['original_launch_v']
        result = restore(arrays, 'original', report['original'])
        tau = 0.
    else:
        if 'conservative' not in report or report['conservative']['stop_reason'] == 'search_in_progress':
            raise RuntimeError('finish conservative inward continuation first')
        launch = arrays['conservative_launch_v']
        result = N.forward(lattice, launch, DT, TAU)
        report['passive_loss'] = N.summarize(result, 'inward', DT)
        store(arrays, 'passive_loss', launch, result)
        report['loss_seed'] = 'new conservative inward winner under actual 4100-fs loss'
        if size == 16384:
            with np.load(HERE/'nonlinear_N16384_0.6meV.npz') as earlier:
                old_v = earlier['loss_K_v'].copy()
                old_trace = earlier['loss_K_trace'].copy()
                old_x, old_final_v = earlier['loss_K_final_x'].copy(), earlier['loss_K_final_v'].copy()
            old_result = old_x, old_final_v, old_trace, -1., 0
            old_readout = N.summarize(old_result, 'inward', DT)
            report['previous_loss_designed_inward'] = old_readout
            store(arrays, 'previous_loss', old_v, old_result)
            if old_readout['score'] > report['passive_loss']['score']:
                launch, result = old_v, old_result
                report['loss_seed'] = 'exact previous loss-designed launch; strongest inward pulse, not total K'
        tau = TAU
    report['settings'].update(budget_seconds=budget, max_updates=max_updates,
                              max_trials_per_update=max_trials, first_step_predicted_fractional_gain=.05)
    report['settings']['applies_to_stage'] = stage
    def checkpoint(v, result, record, gradient_state=None):
        report[prefix] = record
        store(arrays, prefix, v, result)
        if gradient_state is not None:
            for suffix, value in zip(('_gradient', '_gradient_endpoint_x', '_gradient_endpoint_v'), gradient_state):
                arrays[prefix+suffix] = value.copy()
        save(size, report, arrays)
    launch, result, record = optimize(lattice, launch, result, tau, budget, checkpoint,
        check_gradient=size == 131072 and stage == 'conservative', previous=previous,
        max_updates=max_updates, max_trials=max_trials)
    checkpoint(launch, result, record)
    print(json.dumps(dict(size=size, stage=stage, final=record['final'],
                          stop=record['stop_reason']), indent=2), flush=True)


def loss_subspace(size, budget, max_trials=2):
    """Actual lossy feedback along already available conservative directions."""
    started = time.monotonic()
    lattice = LargeLattice({16384: 16, 131072: 32}[size])
    path = HERE/('inward_N'+str(size)+'.json')
    report = json.loads(path.read_text())
    if 'loss' in report:
        raise RuntimeError('loss result already exists')
    if report['conservative']['stop_reason'] == 'search_in_progress':
        raise RuntimeError('finish the conservative update first')
    with np.load(path.with_suffix('.npz')) as saved:
        arrays = {key: saved[key].copy() for key in saved.files}
    base = arrays['conservative_launch_v']
    tick = time.monotonic()
    passive = N.forward(lattice, base, DT, TAU)
    native_seconds = time.monotonic()-tick
    current = N.summarize(passive, 'inward', DT)
    report['passive_loss'] = current
    store(arrays, 'passive_loss', base, passive)
    delta = base-arrays['original_launch_v']
    direction_source = 'exact accepted conservative launch increment'
    if np.linalg.norm(delta) < 1e-14*np.linalg.norm(base):
        gradient = arrays['conservative_gradient'].copy()
        gradient[list(lattice.target)] = 0.
        gradient -= base*(np.sum(base*gradient)/np.sum(base*base))
        gradient *= np.linalg.norm(base)/np.linalg.norm(gradient)
        alpha = report['conservative']['feasible_trials'][0]['alpha']
        projected = P.capped_direct_velocity(base+alpha*gradient, lattice.target,
            np.sqrt(2*ENERGY/E.MASS), lattice.target_axis)
        delta = projected-base
        direction_source = 'saved conservative gradient and its first measured feasible proposal scale'
    best_v, best_result = base.copy(), passive
    record = dict(initial=current, accepted_history=[dict(iteration=0, **current)],
        feasible_trials=[], costs=[dict(kind='matched_passive_loss', seconds=native_seconds)],
        searched_subspace=direction_source+'; one-dimensional projected line about the fixed conservative winner',
        seed_lineage='matched conservative inward winner under actual loss',
        global_optimality=False, resolved_local_stop=False,
        combined_allowance_remaining_seconds=budget)
    for iteration, sign in enumerate((1., -1.)[:max_trials], 1):
        if time.monotonic()-started+1.1*native_seconds > budget:
            break
        proposal = P.capped_direct_velocity(base+sign*delta, lattice.target,
            np.sqrt(2*ENERGY/E.MASS), lattice.target_axis)
        assert cap_diagnostics(lattice, proposal)['valid']
        tick = time.monotonic()
        candidate_result = N.forward(lattice, proposal, DT, TAU)
        elapsed = time.monotonic()-tick
        candidate = N.summarize(candidate_result, 'inward', DT)
        gain = candidate['score']/current['score']-1.
        record['feasible_trials'].append(dict(iteration=iteration, direction_sign=sign,
            seconds=elapsed, relative_gain=gain, full_horizon_scanned=True, **candidate))
        if candidate['score'] > current['score']*(1+1e-6):
            best_v, best_result, current = proposal, candidate_result, candidate
            record['accepted_history'].append(dict(iteration=iteration, direction_sign=sign, **current))
        record.update(final=current, launch_constraints=cap_diagnostics(lattice, best_v),
            stop_reason='search_in_progress', elapsed_seconds=time.monotonic()-started)
        report['loss'] = record
        store(arrays, 'loss', best_v, best_result)
        save(size, report, arrays)
    record.update(final=current, launch_constraints=cap_diagnostics(lattice, best_v),
        stop_reason=('bounded_retained_direction_subspace_complete_no_plateau'
                     if len(record['feasible_trials']) == max_trials else
                     'out_of_proportion_remaining_combined_allowance_no_plateau'),
        optimization_attempted=bool(record['feasible_trials']),
        tested_direction_signs=[trial['direction_sign'] for trial in record['feasible_trials']],
        unsearched_direction_signs=[sign for sign in (1., -1.)
            if sign not in [trial['direction_sign'] for trial in record['feasible_trials']]],
        reoptimized=not np.array_equal(best_v, base),
        elapsed_seconds=time.monotonic()-started,
        method='bounded derivative-free actual-loss feedback; no new adjoint or proposal bank')
    report['loss'] = record
    store(arrays, 'loss', best_v, best_result)
    save(size, report, arrays)
    print(json.dumps(dict(size=size, stage='loss', final=current,
        reoptimized=record['reoptimized'], tested_directions=len(record['feasible_trials']),
        stop=record['stop_reason']), indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('size', type=int, choices=(16384, 131072))
    parser.add_argument('stage', choices=('conservative', 'loss'))
    parser.add_argument('--budget-seconds', type=float, default=5400.)
    parser.add_argument('--max-updates', type=int, default=60)
    parser.add_argument('--max-trials', type=int, default=20)
    parser.add_argument('--derivative-free', action='store_true',
                        help='actual loss feedback along retained conservative launch directions')
    parser.add_argument('--resume', action='store_true',
                        help='resume only an interrupted exact accepted checkpoint, retaining its history')
    args = parser.parse_args()
    if args.derivative_free:
        if args.stage != 'loss':
            parser.error('--derivative-free is only for loss')
        loss_subspace(args.size, args.budget_seconds, min(2, args.max_trials))
    else:
        run(args.size, args.stage, args.budget_seconds, args.resume, args.max_updates, args.max_trials)
