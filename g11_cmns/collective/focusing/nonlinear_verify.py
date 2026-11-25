"""One settled changed-consumer exercise and exact saved-output contract check."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

import engine as E
import nonlinear as N
import exterior as X
from adjoint import endpoint_gradient, _terminal
from large_geometry import LargeLattice
from nonlinear_select import selected_result

HERE = Path(__file__).resolve().parent


def endpoint(lattice, v, dt, step, tau, objective):
    result = E.trajectory(lattice,np.zeros_like(v),v,dt,step*dt,damping=tau)
    if result[3] >= 0:
        raise RuntimeError('derivative witness crossed original support')
    x,vf = result[:2]
    return _terminal(lattice.target_bond,lattice.target,x,vf,objective)[0],x,vf


def verify(replay_only=False):
    records, states, geometries = {}, {}, {}
    native_replays = {}
    selections = json.loads((HERE/'nonlinear_selection.json').read_text())
    selected_arrays = dict(np.load(HERE/'nonlinear_selection.npz'))
    for key,(m,energy,dt,source) in N.CASES.items():
        report = json.loads((HERE/('nonlinear_'+key+'.json')).read_text())
        required = {*N.OBJECTIVES,'loss_K'}
        assert set(report['searches']) == required, (key,'missing required independent objective')
        assert set(selections[key]['selected']) == required, (key,'missing delivered objective')
        lattice = LargeLattice(m)
        omitted_bindings = []
        selected_constraints = {}
        with np.load(HERE/('nonlinear_'+key+'.npz')) as arrays:
            for name,search in report['searches'].items():
                v,trace = arrays[name+'_v'],arrays[name+'_trace']
                assert N.constraints(lattice,v,energy)['valid']
                objective = 'K' if name == 'loss_K' else name
                exit_fs = search['final']['first_original_exit_fs']
                result = (None,None,trace,-1. if exit_fs is None else exit_fs,
                    {None:0,'site_radius':1,'bond_distance':2}[search['final']['exit_kind']])
                rebuilt = N.summarize(result,objective,dt)
                assert rebuilt == search['final'], (key,name,'saved readout drift')
                assert all(b['score'] > a['score'] for a,b in
                    zip(search['accepted_history'],search['accepted_history'][1:]))
                for event in rebuilt['joint_events'].values():
                    assert exit_fs is None or event['time_fs'] < exit_fs
                    assert event['inward_K_eV'] == 0. or event['signed_inward_A_fs'] > 0.
            for name,outcome in selections[key]['selected'].items():
                v,stored_result = selected_result(key,name,outcome,arrays,selected_arrays)
                trace = stored_result[2]
                selected_constraints[name] = N.constraints(lattice,v,energy)
                assert selected_constraints[name]['valid']
                row = outcome['readout']
                exit_fs = row['first_original_exit_fs']
                result = (None,None,trace,-1. if exit_fs is None else exit_fs,
                    {None:0,'site_radius':1,'bond_distance':2}[row['exit_kind']])
                objective = 'K' if name == 'loss_K' else name
                assert N.summarize(result,objective,dt) == row
                if name != 'loss_K':
                    for other,run in report['searches'].items():
                        if other != 'loss_K':
                            assert all(row['score']+2e-9 >= t['joint_events'][name]['objective_value']
                                       for t in run['feasible_trials'])
                else:
                    assert all(row['score']+2e-9 >= t['joint_events']['K']['objective_value']
                               for t in report['searches']['loss_K']['feasible_trials'])
                    assert not selections[key]['loss_comparison']['new_conservative_seed_beats_reoptimized_loss'], (
                        key,'winning conservative seed still needs actual-loss feedback')
                selected = row['selected_event']
                tau = N.TAU if name == 'loss_K' else 0.
                replay_key = (key,dt,selected['step'],tau,
                              hashlib.sha256(np.ascontiguousarray(v)).hexdigest())
                if replay_key not in native_replays:
                    native_replays[replay_key] = E.trajectory(
                        lattice,np.zeros_like(v),v,dt,selected['time_fs'],damping=tau)
                replay = native_replays[replay_key]
                assert replay[3] < 0
                saved = trace[selected['step'], :10]
                assert np.allclose(replay[2][-1,:10],saved,rtol=2e-11,atol=2e-12)
                state_key = key+'__'+name
                states[state_key+'__x'],states[state_key+'__v'] = replay[:2]
                geometries[state_key] = dict(
                    selected=X.state_geometry(lattice,replay[0],selected['time_fs']),
                    source_search=outcome['source_search'],
                    source_trial_index=outcome['source_trial_index'],
                    native_replay_bitwise_first10_diagnostics=bool(np.array_equal(replay[2][-1,:10],saved)),
                    native_replay_absolute_differences_first10_diagnostics=np.abs(replay[2][-1,:10]-saved).tolist(),
                    exterior_used_before_scored_target=False)
                if exit_fs is not None:
                    geometries[state_key]['first_original_exit'] = X.state_geometry(
                        lattice,stored_result[0],exit_fs)
                    states[state_key+'__first_exit_x'] = stored_result[0]
                    states[state_key+'__first_exit_v'] = stored_result[1]
            for name,record in report['inputs'].items():
                if replay_only and not (HERE/name).exists():
                    omitted_bindings.append(name)
                else:
                    assert N.digest(HERE/name) == record, ('immutable source changed',name)
        inferred_prefilter = any('full_horizon_scanned' in t
            for s in report['searches'].values() for t in s['feasible_trials'])
        records[key] = dict(saved_launch_constraints=True,all_step_saved_readouts_exact=True,
            selected_launch_constraints=selected_constraints,
            original_bank_comparison={objective:dict(
                original_score=report['original_baseline'][objective]['score'],
                selected_score=selections[key]['selected'][objective]['readout']['score'],
                gain_fraction=selections[key]['selected'][objective]['readout']['score']/
                    report['original_baseline'][objective]['score']-1.,
                score_units='Angstrom approach Rref-R' if objective=='distance' else 'eV')
                for objective in N.OBJECTIVES},
            actual_minimum_target_R_reduction_A=report['original_baseline']['distance']['selected_event']['R_A']-
                selections[key]['selected']['distance']['readout']['selected_event']['R_A'],
            accepted_history_strictly_improving=True,
            every_selected_objective_dominates_every_recorded_feasible_trial=True,
            selected_local_stops={name:s['objective_matched_local_stop'] for name,s in selections[key]['selected'].items()},
            available_original_bindings_unchanged=True,
            intentionally_omitted_original_bindings=omitted_bindings,
            execution_settings=report.get('optimization_settings',dict(
                endpoint_prefilter=inferred_prefilter, initial_RMS_step=.03,
                max_RMS_step=.12, minimum_RMS_step=1e-5,
                accepted_relative_gain_floor=1e-6,meaningful_relative_gain=.002,
                plateau_successes=3,plateau_challenge_multiplier=8,
                iteration_safety_cap=60,numerical_threads=1,
                per_objective_budget_seconds=1425. if m==16 else 1500.)),
            resolved={name:s['resolved_local_stop'] for name,s in report['searches'].items()})
    frozen = []
    if not replay_only:
        old_manifest = json.loads((HERE/'continuation_manifest.json').read_text())
        frozen = [item for item in old_manifest['files'] if Path(item['path']).name
                  not in ('README.md','method.md')]
        for item in frozen:
            assert N.digest(HERE.parents[2]/item['path']) == {
                'bytes':item['bytes'],'sha256':item['sha256']}, item['path']
    np.savez_compressed(HERE/'nonlinear_states.npz',**states)
    (HERE/'nonlinear_geometry.json').write_text(json.dumps(geometries,indent=2,allow_nan=False)+'\n')
    # The actual high-energy saved launch is used at a supported, nonlinear endpoint.
    # This is not a harmonic toy, and the loss omission is a plausible wrong callback.
    key = 'N2048_25meV'
    lattice = LargeLattice(8)
    with np.load(HERE/('nonlinear_'+key+'.npz')) as arrays:
        v = arrays['original_v'].copy()
    dt, step = .25, 2400  # 600 fs, fixed in advance and well inside the original extent.
    rng = np.random.default_rng(37625)
    direction = rng.normal(size=v.shape)
    cap = N.constraints(lattice,v,.025)['cap_A_fs']
    free = np.linalg.norm(v,axis=1) < cap*(1-1e-8)
    free[list(lattice.target)] = False
    direction[~free] = 0.
    direction[free] -= v[free]*(np.sum(v*direction)/np.sum(v[free]**2))
    direction *= np.linalg.norm(v)/np.linalg.norm(direction)
    epsilon = 1e-7
    witnesses = []
    for tau in (0.,N.TAU):
        value,g,x,vf = endpoint_gradient(lattice,v,dt,step,'K',tau_fs=tau)
        actual,xa,va = endpoint(lattice,v,dt,step,tau,'K')
        plus,_,_ = endpoint(lattice,v+epsilon*direction,dt,step,tau,'K')
        minus,_,_ = endpoint(lattice,v-epsilon*direction,dt,step,tau,'K')
        fd = (plus-minus)/(2*epsilon)
        predicted = float(np.sum(g*direction))
        relative = abs(predicted-fd)/max(abs(fd),abs(predicted),1e-10)
        assert np.max(abs(x-xa)) < 2e-11 and np.max(abs(vf-va)) < 2e-12
        assert abs(value-actual) < 2e-11
        assert relative < .015, ('actual nonlinear gradient mismatch',tau,fd,predicted)
        assert abs(predicted-fd) < abs(-predicted-fd), 'opposite gradient sign was not rejected'
        witnesses.append(dict(tau_fs=tau,objective_value=value,directional_fd=fd,
            gradient_dot_direction=predicted,relative_defect=relative,
            endpoint_x_max_abs_difference_A=float(np.max(abs(x-xa))),
            endpoint_v_max_abs_difference_A_fs=float(np.max(abs(vf-va))),
            opposite_sign_rejected=True))
    assert abs(witnesses[0]['directional_fd']-witnesses[1]['directional_fd']) > .01*abs(witnesses[1]['directional_fd'])
    # A genuinely receding endpoint must not be rewarded by the closing objective.
    report = json.loads((HERE/('nonlinear_'+key+'.json')).read_text())
    total = report['original_baseline']['K']['selected_event']
    value,g,x,vf = endpoint_gradient(lattice,v,dt,total['step'],'inward')
    assert total['receding'] and value == 0. and not np.any(g)
    proof = dict(schema='pdd-nonlinear-consumer-proof-v1',cases=records,
        actual_nonlinear_gradient_witnesses=witnesses,
        omitted_loss_callback_rejected=True,actual_receding_inward_value_and_gradient_zero=True,
        exact_original_files_unchanged=None if replay_only else len(frozen),
        immutable_source_check_scope='not repeated; original old outputs intentionally absent' if replay_only else 'all prior manifest entries except the two edited pages',
        actual_selected_endpoint_replays=len(geometries),
        unique_native_endpoint_propagations=len(native_replays),
        exterior_dependent_supported_winners=0,
        finite_difference_scope='unconstrained launch derivative along a target-zero sphere tangent; not a global optimum proof',
        interpretation='finite-model implementation/consumer evidence, not timestep or electronic validation',
        posthoc_scope='kinetic exponential applied to the conservative trace only; its R/extent/total E fields, if present, remain conservative reference geometry and are not a lossy simultaneous state',
        archive_roundtrip='recorded independently by nonlinear_package.py')
    (HERE/'nonlinear_proof.json').write_text(json.dumps(proof,indent=2,allow_nan=False)+'\n')
    print(json.dumps(proof,indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay-only',action='store_true',
                        help='exercise capsule without intentionally omitted old outputs')
    verify(parser.parse_args().replay_only)
