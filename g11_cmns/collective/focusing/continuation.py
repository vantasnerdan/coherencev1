"""Ordered PdD continuation; no contact/rate adapter is imported or called."""
from pathlib import Path
import argparse
import json
import resource
import time
from multiprocessing import get_context
import numpy as np
import engine as E
import patterns as P

HERE = Path(__file__).resolve().parent
RESULTS = HERE/'continuation_results.json'
TRACES = HERE/'continuation_traces.npz'
HORIZON = 2400.


def focus_times(lattice, speed):
    """The unchanged run.py proposal grid and three-time spacing rule."""
    scored = []
    for t in np.arange(120., HORIZON-119., 20.):
        h = lattice.predicted_kernel(t)
        v = P.capped_direct_velocity(h,lattice.target,speed,lattice.target_axis)
        scored.append((float(np.sum(h*v)),float(t)))
    chosen = []
    for score,t in sorted(scored,reverse=True):
        if all(abs(t-u)>=160 for u in chosen):
            chosen.append(t)
        if len(chosen)==3:
            break
    return sorted(chosen)


def save(data, traces):
    RESULTS.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    np.savez_compressed(TRACES, **traces)


def load():
    data = json.loads(RESULTS.read_text()) if RESULTS.exists() else dict(
        schema='pdd-ordered-continuation-v1', horizon_fs=HORIZON,
        trace_columns=['time_fs','target_Krel_eV','target_Krad_eV','target_R_A',
                       'max_site_radius_A','min_NN_R_A','total_E_eV','potential_E_eV',
                       'same_trace_max_NN_Krel_eV','target_signed_inward_A_fs',
                       'max_site_local_E_eV','max_site_harmonic_E_eV'], cases={}, resources={})
    traces = dict(np.load(TRACES)) if TRACES.exists() else {}
    return data, traces


def _worker_init(lattice,dt):
    global _worker_lattice, _worker_dt
    _worker_lattice, _worker_dt = lattice,dt


def _worker_run(launch):
    tick = time.monotonic()
    result = E.trajectory(_worker_lattice,launch['q'],launch['v'],_worker_dt,HORIZON)
    return result,time.monotonic()-tick


def search(lattice, energy, dt, key, data, traces, waves=False, evolve=None, workers=1):
    """Exactly run.py's three proposals, four phase and four mixture trials."""
    started = time.monotonic()
    data.pop('archive_trace_policy',None)
    speed = np.sqrt(2*energy/E.MASS)
    times = focus_times(lattice, speed)
    rows, launches = [], {}
    propagate = evolve or (lambda q,v: E.trajectory(lattice,q,v,dt,HORIZON))
    pool = get_context('fork').Pool(workers,_worker_init,(lattice,dt)) if workers>1 else None
    def score(launch, timed):
        name,q,v = launch['name'],launch['q'],launch['v']
        result,elapsed = timed
        pattern = 'converging_wave' if name.startswith('wave_') else 'direct_optimized'
        scored = result
        if evolve is not None:
            # Full extended trace, except an invalid model-failure endpoint.
            trace = result[2] if result[5]==0 else result[2][:-1]
            scored = (result[0],result[1],trace,-1.,0)
        row = E.summarize(lattice,q,v,scored,energy,pattern,name)
        if evolve is not None:
            crossing = None if result[3]<0 else float(result[3])
            outside_peak = crossing is not None and row['peak_time_fs']>=crossing
            window = scored[2][abs(scored[2][:,0]-row['peak_time_fs'])<=10.]
            row.update(first_site_or_bond_exit_fs=crossing,
                exit_kind={0:None,1:'site_radius',2:'bond_distance'}[result[4]],
                trajectory_supported=crossing is None,candidate_supported=not outside_peak,
                peak_before_support_exit=not outside_peak,peak_after_original_exit=outside_peak,
                peak_in_original_extent=bool(scored[2][row['peak_index'],4]<=.8 and
                    scored[2][row['peak_index'],5]>=1.28782409436586),
                peak_neighbourhood_in_original_extent=bool(np.all(window[:,4]<=.8) and
                    np.all(window[:,5]>=1.28782409436586)),
                model_valid_full_horizon=result[5]==0,
                model_failure={0:None,1:'near_collision',2:'reference_lineage_crossing',3:'numerical_energy_failure'}[result[5]],
                model_failure_time_fs=None if result[5]==0 else float(result[2][-1,0]),
                interpretation='artificial matched exterior; not computed electronic support')
            status = np.c_[result[2][:,4]>.8,result[2][:,5]<1.28782409436586]
            changes = np.flatnonzero(np.r_[np.any(status[0]),np.any(status[1:]!=status[:-1],axis=1)])
            row['all_step_original_support_transitions'] = [
                dict(time_fs=float(result[2][k,0]),site_outside=bool(status[k,0]),bond_outside=bool(status[k,1]))
                for k in changes]
            row['outside_original_extent_steps'] = int(np.any(status,axis=1).sum())
        row.update(E.launch_energy_diagnostics(lattice,q,v))
        row.update(dt_fs=dt, headline_eligible=bool(
            row['participation_positive_site_energy']>=lattice.n/4 and
            row['target_positive_site_energy_fraction']<.05), runtime_seconds=elapsed)
        for k in ('wave_mode_count','wave_subspace_rank'):
            if k in launch:
                row[k] = launch[k]
        row['proposed_focus_time_fs'] = launch.get('focus_time_fs')
        rows.append(row)
        traces[key+'__'+name] = result[2]
        launches[name] = (q,v)
        print(key,name,row['peak_relative_kinetic_eV'],row['peak_time_fs'],
              row['first_site_or_bond_exit_fs'],elapsed,flush=True)
        return row
    def evaluate(batch):
        if pool:
            # imap preserves declared ordering and therefore exact first-wins ties.
            timed = pool.imap(_worker_run,batch,chunksize=1)
            return [score(launch,result) for launch,result in zip(batch,timed)]
        answer = []
        for launch in batch:
            tick = time.monotonic()
            result = propagate(launch['q'],launch['v'])
            answer.append(score(launch,(result,time.monotonic()-tick)))
        return answer
    try:
        base = []
        for kind in (('wave','direct') if waves else ('direct',)):
            base.extend(P.velocity_candidates(kind,lattice.positions,lattice.box,lattice.target,
                lattice.target_axis,speed,lattice.predicted_kernel,focus_times_fs=times))
        base_rows = evaluate(base)
        direct = [r for r in base_rows if r['pattern']=='direct_optimized']
        best = max(direct,key=lambda r:r['peak_relative_kinetic_eV'])
        q0,v0 = launches[best['candidate']]
        refine = [dict(name=f'direct_phase_{k}',q=q0,v=v)
                  for k,v in enumerate(P.direct_refinements(v0,lattice.target,lattice.target_axis))]
        other = next(r for r in direct if r is not best)
        vb = launches[other['candidate']][1]
        for c in (-.5,-.25,.25,.5):
            v = P.capped_direct_velocity(v0+c*vb,lattice.target,speed,lattice.target_axis)
            refine.append(dict(name=f'direct_mixture_{c:+.2f}',q=q0,v=v))
        evaluate(refine)
    finally:
        if pool:
            pool.close()
            pool.join()
    selected = []
    for pattern in (('converging_wave','direct_optimized') if waves else ('direct_optimized',)):
        options = [r for r in rows if r['pattern']==pattern]
        viable = [r for r in options if r.get('model_valid_full_horizon',True)]
        eligible = [r for r in viable if r['headline_eligible']]
        row = max(eligible or viable or options,key=lambda r:r['peak_relative_kinetic_eV'])
        row.update(selected=True, selection_candidates=len(options))
        selected.append(row)
    data['cases'][key] = dict(n_D=lattice.n, energy_per_D_eV=energy, dt_fs=dt,
        proposed_focus_times_fs=times, candidates=rows, selected=selected,
        runtime_seconds=time.monotonic()-started, numerical_worker_processes=workers,
        selection_rule=('highest full extended model-valid peak; original crossings logged, not censored'
                        if evolve is not None else 'highest supported-prefix peak')+
                       '; eligible launch preferred; exact 3+4+4 direct search')
    if evolve is not None:
        data['cases'][key]['trace_columns'] = data['trace_columns'][:10]
    save(data,traces)
    return selected, launches


def probe(data, traces):
    from large_geometry import LargeLattice
    # Compile the exact int32 production signature before timing the large register.
    warm = LargeLattice(4)
    q = np.zeros((warm.n,3)); v = q.copy()
    E.trajectory(warm,q,v,.5,1.)
    tick = time.monotonic()
    lattice = LargeLattice(32)
    geometry_seconds = time.monotonic()-tick
    speed = np.sqrt(2*.0006/E.MASS)
    tick = time.monotonic()
    times = focus_times(lattice,speed)
    proposal_seconds = time.monotonic()-tick
    launch = next(P.velocity_candidates('direct',lattice.positions,lattice.box,lattice.target,
        lattice.target_axis,speed,lattice.predicted_kernel,focus_times_fs=times))
    tick = time.monotonic()
    result = E.trajectory(lattice,launch['q'],launch['v'],.5,200.)
    seconds = time.monotonic()-tick
    tested_fs = float(result[2][-1,0])
    estimate = seconds*HORIZON/tested_fs*14+proposal_seconds+geometry_seconds
    row = E.summarize(lattice,launch['q'],launch['v'],result,.0006,'direct_probe',launch['name'])
    record = dict(n_D=lattice.n, dt_fs=.5, tested_horizon_fs=tested_fs,
        production_candidates=14, support_policy='stop at first original all-step support exit',
        geometry_seconds=geometry_seconds,proposal_seconds=proposal_seconds,
        trajectory_seconds=seconds, full_search_linear_seconds=estimate,
        rss_peak_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        probe=row, credible_under_90_minutes=bool(estimate<5400),
        caveat='linear measured estimate, not a runtime guarantee; later nonlinear local evaluations may cost more')
    data['resources']['N131072'] = record
    traces['resource_N131072'] = result[2]
    save(data,traces)
    print(json.dumps(record,indent=2),flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('stage',choices=('probe','size16','size32','fine','outside'))
    args = parser.parse_args()
    data,traces = load()
    if args.stage=='probe':
        probe(data,traces)
    elif args.stage.startswith('size'):
        from large_geometry import LargeLattice
        m = int(args.stage[4:])
        search(LargeLattice(m),.0006,.5,f'size_N{4*m**3}',data,traces,waves=True,workers=2)
    elif args.stage=='fine':
        search(E.Lattice(8),.025,.25,'fine_N2048_25meV',data,traces)
    else:
        from exterior import campaign
        campaign(data,traces)


if __name__=='__main__':
    main()
