"""Select every observed feasible trial, recovering only an unsaved actual champion."""
from pathlib import Path
import argparse
import json
import time

import numpy as np

import nonlinear as N
import patterns as P
from adjoint import endpoint_gradient
from large_geometry import LargeLattice

HERE = Path(__file__).resolve().parent


def recover(lattice, arrays, report, name, indices):
    """Replay shared accepted-prefix gradients; never repeat the old candidate bank."""
    run = report['searches'][name]
    loss = name == 'loss_K'
    objective, tau = ('K',N.TAU) if loss else (name,0.)
    start = 'K_v' if loss and report['loss_comparison']['loss_search_seed'] == 'conservative_optimum' else 'original_v'
    v = arrays[start].copy()
    current = run['initial']
    trials = run['feasible_trials']
    requested = set(indices)
    last = max(trials[k]['iteration'] for k in indices)
    recovered = {}
    dt, energy = report['dt_fs'], report['energy_per_D_eV']
    tick = time.monotonic()
    for iteration in range(1,last+1):
        _,g,_,_ = endpoint_gradient(lattice,v,dt,current['selected_event']['step'],objective,tau_fs=tau)
        g[list(lattice.target)] = 0.
        g -= v*(np.sum(v*g)/np.sum(v*v))
        direction = g*(np.linalg.norm(v)/np.linalg.norm(g))
        def project(base,alpha):
            return P.capped_direct_velocity(base+alpha*direction,lattice.target,
                np.sqrt(2*energy/N.E.MASS),lattice.target_axis)
        for k in requested:
            trial = trials[k]
            if trial['iteration'] == iteration and not trial.get('plateau_challenge'):
                recovered[k] = project(v,trial['alpha'])
        accepted = [a for a in run['accepted_history'] if a['iteration'] == iteration]
        for event in accepted:
            if event['source'] == 'larger-step plateau challenge':
                challenge = next(t for t in trials if t['iteration'] == iteration and
                                 t.get('plateau_challenge') and t['relative_gain'] > .002)
                v = project(v,challenge['alpha'])
            else:
                v = project(v,event['alpha'])
            current = event
            # Challenge proposals use the just-accepted base but the same carried direction.
            if event['source'] != 'larger-step plateau challenge':
                for k in requested:
                    trial = trials[k]
                    if trial['iteration'] == iteration and trial.get('plateau_challenge'):
                        recovered[k] = project(v,trial['alpha'])
    assert set(recovered) == requested
    return recovered,time.monotonic()-tick


def selected_result(key, name, outcome, case_arrays, bound_arrays):
    """Resolve a final-launch reference or the exact recovered unsaved trial."""
    if outcome['source_trial_index'] is None:
        source = outcome['source_search']
        launch = case_arrays[source+'_v']
        x,v,trace = (case_arrays[source+suffix] for suffix in ('_final_x','_final_v','_trace'))
    else:
        prefix = key+'__'+name+'__'
        launch = bound_arrays[prefix+'launch_v']
        x,v,trace = (bound_arrays[prefix+suffix] for suffix in ('final_x','final_v','trace'))
    row = outcome['readout']
    exit_fs = row['first_original_exit_fs']
    return launch,(x,v,trace,-1. if exit_fs is None else exit_fs,
                   {None:0,'site_radius':1,'bond_distance':2}[row['exit_kind']])


def matched_loss(lattice, launch, dt, key):
    """Reuse only a byte-identical launch and identical declared loss conditions."""
    path = HERE/('nonlinear_'+key+'_matched_loss.npz')
    fingerprint = json.dumps({name:N.digest(HERE/name) for name in
        ('engine.py','exterior.py','patterns.py','large_geometry.py','input_manifest.json')},sort_keys=True)
    if path.exists():
        with np.load(path) as cache:
            if (np.array_equal(cache['launch_v'],launch) and float(cache['dt']) == dt and
                float(cache['tau']) == N.TAU and float(cache['horizon']) == N.HORIZON and
                str(cache['fingerprint']) == fingerprint):
                return (cache['x'].copy(),cache['v'].copy(),cache['trace'].copy(),
                        float(cache['exit_time']),int(cache['exit_kind']))
    result = N.forward(lattice,launch,dt,N.TAU)
    np.savez_compressed(path,launch_v=launch,dt=dt,tau=N.TAU,horizon=N.HORIZON,
        fingerprint=fingerprint,x=result[0],v=result[1],trace=result[2],
        exit_time=result[3],exit_kind=result[4])
    return result


def prepare_loss(key):
    report = json.loads((HERE/('nonlinear_'+key+'.json')).read_text())
    pool = {}
    for name in N.OBJECTIVES:
        run = report['searches'][name]
        pool[name,-1] = run['final']['joint_events']['K']['objective_value']
        pool.update({(name,k):t['joint_events']['K']['objective_value']
                     for k,t in enumerate(run['feasible_trials'])})
    name,index = max(pool,key=lambda ref:(pool[ref],ref[1] == -1))
    if index != -1:
        raise RuntimeError('the winning conservative trial must first have its exact launch bound')
    with np.load(HERE/('nonlinear_'+key+'.npz')) as arrays:
        launch = arrays[name+'_v'].copy()
    lattice = LargeLattice(N.CASES[key][0])
    result = matched_loss(lattice,launch,report['dt_fs'],key)
    print(json.dumps(dict(case=key,source_search=name,
        selected_conservative_under_actual_loss=N.summarize(result,'K',report['dt_fs'])),indent=2),flush=True)


def select(allow_large_replay_seconds=0.,only_case=None):
    selection_path, arrays_path = HERE/'nonlinear_selection.json', HERE/'nonlinear_selection.npz'
    selection = json.loads(selection_path.read_text()) if only_case and selection_path.exists() else {}
    saved = dict(np.load(arrays_path)) if only_case and arrays_path.exists() else {}
    for key,(m,energy,dt,source) in N.CASES.items():
        if only_case is not None and key != only_case:
            continue
        report = json.loads((HERE/('nonlinear_'+key+'.json')).read_text())
        arrays = dict(np.load(HERE/('nonlinear_'+key+'.npz')))
        lattice = LargeLattice(m)
        # A record is (source search, trial index); -1 denotes its saved final launch.
        pool = {}
        for name in N.OBJECTIVES:
            run = report['searches'][name]
            pool[name,-1] = run['final']
            for k,trial in enumerate(run['feasible_trials']):
                pool[name,k] = trial
        champions = {}
        for objective in N.OBJECTIVES:
            # Prefer an already saved final launch on exact scalar ties.
            champions[objective] = max(pool,key=lambda ref:(
                pool[ref]['joint_events'][objective]['objective_value'],ref[1] == -1))
        loss_run = report['searches']['loss_K']
        loss_pool = {('loss_K',-1):loss_run['final']}
        loss_pool.update({('loss_K',k):t for k,t in enumerate(loss_run['feasible_trials'])})
        champions['loss_K'] = max(loss_pool,key=lambda ref:(
            loss_pool[ref]['joint_events']['K']['objective_value'],ref[1] == -1))
        needed = {}
        for name,k in set(champions.values()):
            if k != -1:
                needed.setdefault(name,[]).append(k)
        replay_estimate = 0.
        for name,indices in needed.items():
            last = max(report['searches'][name]['feasible_trials'][k]['iteration'] for k in indices)
            replay_estimate += sum(c['seconds'] for c in report['searches'][name]['costs'][:last])
        if m == 16 and replay_estimate > allow_large_replay_seconds:
            raise RuntimeError(json.dumps(dict(case=key,required_exact_champion_replay=needed,
                measured_prefix_gradient_seconds_estimate=replay_estimate,
                decision='explicit replay allocation needed; this is not a new launch search')))
        actual, replay_records = {}, {}
        for name,indices in needed.items():
            recovered,elapsed = recover(lattice,arrays,report,name,indices)
            replay_records[name] = dict(trial_indices=indices,gradient_replay_seconds=elapsed)
            for k,v in recovered.items():
                loss = name == 'loss_K'
                result = N.forward(lattice,v,dt,N.TAU if loss else 0.)
                actual[name,k] = (v,result)
                destination = loss_pool if loss else pool
                destination[name,k] = N.summarize(result,'K' if loss else name,dt)
                # The original observed joint event must survive exact launch reconstruction.
                old = report['searches'][name]['feasible_trials'][k]
                for objective in (('loss_K',) if loss else N.OBJECTIVES):
                    event = 'K' if loss else objective
                    if champions[objective] == (name,k):
                        assert destination[name,k]['joint_events'][event]['objective_value'] >= old['joint_events'][event]['objective_value']-2e-9
        outcomes = {}
        for objective in N.OBJECTIVES:
            # Full completion of a former prefix champion may improve another objective too.
            ref = max(pool,key=lambda r:(pool[r]['joint_events'][objective]['objective_value'],r[1] == -1))
            name,k = ref
            if k == -1:
                v = arrays[name+'_v']
                row = report['searches'][name]['final']
                exit_fs = row['first_original_exit_fs']
                result = (arrays[name+'_final_x'],arrays[name+'_final_v'],arrays[name+'_trace'],
                    -1. if exit_fs is None else exit_fs,
                    {None:0,'site_radius':1,'bond_distance':2}[row['exit_kind']])
            else:
                v,result = actual[ref]
            row = N.summarize(result,objective,dt)
            run = report['searches'][name]
            matched = name == objective and k == -1
            event = row['selected_event']
            shared = (objective == 'inward' and name == 'K' and k == -1 and
                event['signed_inward_A_fs'] > 0 and
                abs(event['Krel_eV']-event['inward_K_eV']) < 2e-11 and
                event['step'] == run['final']['selected_event']['step'])
            resolved = bool((matched or shared) and run['resolved_local_stop'])
            outcome = dict(source_search=name,source_trial_index=None if k == -1 else k,
                readout=row,objective_matched_local_stop=resolved,
                stop_evidence=run['stop_reason'] if matched else
                    'identical inward radial endpoint gradient and K-trial plateau' if shared else
                    'observed cross-objective champion; objective-matched plateau not established')
            saved[key+'__'+objective+'__launch_v'] = v
            saved[key+'__'+objective+'__trace'] = result[2]
            saved[key+'__'+objective+'__final_x'],saved[key+'__'+objective+'__final_v'] = result[:2]
            outcomes[objective] = outcome
        loss_ref = max(loss_pool,key=lambda ref:(
            loss_pool[ref]['joint_events']['K']['objective_value'],ref[1] == -1))
        loss_index = loss_ref[1]
        if loss_index == -1:
            loss_result = (arrays['loss_K_final_x'],arrays['loss_K_final_v'],
                arrays['loss_K_trace'],-1. if loss_run['final']['first_original_exit_fs'] is None
                else loss_run['final']['first_original_exit_fs'],
                {None:0,'site_radius':1,'bond_distance':2}[loss_run['final']['exit_kind']])
            loss_v = arrays['loss_K_v']
        else:
            loss_v,loss_result = actual[loss_ref]
        outcomes['loss_K'] = dict(source_search='loss_K',
            source_trial_index=None if loss_index == -1 else loss_index,
            readout=N.summarize(loss_result,'K',dt),
            objective_matched_local_stop=bool(loss_index == -1 and loss_run['resolved_local_stop']),
            stop_evidence=loss_run['stop_reason'] if loss_index == -1 else
                'best observed actual-loss trial; fresh objective plateau not established')
        saved[key+'__loss_K__launch_v'],saved[key+'__loss_K__trace'] = loss_v,loss_result[2]
        saved[key+'__loss_K__final_x'],saved[key+'__loss_K__final_v'] = loss_result[:2]
        conservative_v = saved[key+'__K__launch_v']
        changed = not np.array_equal(conservative_v,arrays['K_v'])
        if changed:
            matched_loss_result = matched_loss(lattice,conservative_v,dt,key)
            conservative_under_loss = N.summarize(matched_loss_result,'K',dt)
            saved[key+'__selected_conservative_under_loss__trace'] = matched_loss_result[2]
        else:
            conservative_under_loss = report['loss_comparison']['conservative_optimum_under_actual_loss']
            saved[key+'__selected_conservative_under_loss__trace'] = arrays['conservative_under_loss_trace']
        posthoc = saved[key+'__K__trace'].copy()
        posthoc[:,1:3] *= np.exp(-2*posthoc[:,0,None]/N.TAU)
        posthoc[:,9] *= np.exp(-posthoc[:,0]/N.TAU)
        loss_comparison = dict(tau_fs=N.TAU,
            original_seed_under_loss=report['loss_comparison']['original_seed_under_loss'],
            selected_conservative_under_actual_loss=conservative_under_loss,
            posthoc_damped_selected_conservative=dict(
                N.joint(posthoc,int(np.argmax(posthoc[:,1])),'K'),
                non_dynamical_diagnostic=True,
                geometry_and_total_energy_source='undamped conservative reference trace'),
            reoptimized_under_loss=outcomes['loss_K']['readout'],
            selected_conservative_differs_from_independent_K_search=changed,
            new_conservative_seed_beats_reoptimized_loss=bool(
                conservative_under_loss['score'] > outcomes['loss_K']['readout']['score']*(1+1e-6)))
        selection[key] = dict(n_D=lattice.n,energy_per_D_eV=energy,dt_fs=dt,horizon_fs=N.HORIZON,
            selected=outcomes,exact_champion_replays=replay_records,loss_comparison=loss_comparison,
            original_baseline=report['original_baseline'],
            observed_pool='every recorded feasible conservative trial plus saved final launches',
            source_case_report='nonlinear_'+key+'.json')
        # A selected final is already owned by its case archive; do not copy its
        # large arrays merely to give the same state another objective label.
        for case_key,case in selection.items():
            case['array_storage'] = 'final references resolve in the case NPZ; recovered trial arrays resolve in nonlinear_selection.npz'
            for objective,outcome in case['selected'].items():
                if outcome['source_trial_index'] is None:
                    for suffix in ('launch_v','trace','final_x','final_v'):
                        saved.pop(case_key+'__'+objective+'__'+suffix,None)
        selection_path.write_text(json.dumps(selection,indent=2,allow_nan=False)+'\n')
        np.savez_compressed(arrays_path,**saved)
    print(json.dumps({key:{o:dict(score=s['readout']['score'],resolved=s['objective_matched_local_stop'],source=s['source_search'])
        for o,s in case['selected'].items()} for key,case in selection.items()},indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--allow-large-replay-seconds',type=float,default=0.)
    parser.add_argument('--case',choices=tuple(N.CASES))
    parser.add_argument('--prepare-loss',action='store_true')
    args = parser.parse_args()
    if args.prepare_loss:
        if args.case is None:
            parser.error('--prepare-loss requires --case')
        prepare_loss(args.case)
    else:
        select(args.allow_large_replay_seconds,args.case)
