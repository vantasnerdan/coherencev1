"""Reproduce the fixed 3×3×4 focusing comparison; run from this directory."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import time
import numpy as np
import engine as E
import patterns as P
import rate_mapping as R

HERE = Path(__file__).resolve().parent
DT, HORIZON = .5, 2400.
ENERGIES = (.0006,.0048,.025)


def save_json(path, value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def focus_times(lattice,speed):
    # Bounded, cheap harmonic proposal only; every selected proposal is then
    # propagated in the nonlinear Hamiltonian. No nonlinear search beyond horizon.
    scored = []
    for t in np.arange(120.,HORIZON-119.,20.):
        h = lattice.predicted_kernel(t)
        v = P.capped_direct_velocity(h,lattice.target,speed,lattice.target_axis)
        scored.append((float(np.sum(h*v)),float(t)))
    chosen = []
    for score,t in sorted(scored,reverse=True):
        if all(abs(t-u)>=160 for u in chosen): chosen.append(t)
        if len(chosen)==3: break
    return sorted(chosen)


def run_case(lattice,energy):
    speed = np.sqrt(2*energy/E.MASS)
    times = focus_times(lattice,speed)
    candidates, histories, launches = [],{},{}
    def score(pattern,name,q,v,meta=None,horizon=HORIZON,harmonic=False,damping=0.):
        result = E.trajectory(lattice,q,v,DT,horizon,damping,harmonic)
        row = E.summarize(lattice,q,v,result,energy,pattern,name)
        row.update(E.launch_energy_diagnostics(lattice,q,v))
        row['headline_eligible'] = bool(row['participation_positive_site_energy'] >= lattice.n/4 and row['target_positive_site_energy_fraction']<.05)
        row['selection_rule'] = 'highest supported pre-exit target peak among declared bounded candidates'
        row['loss_tau_amplitude_fs'] = damping or None
        row['hamiltonian'] = 'harmonic' if harmonic else 'updated_080_connected'
        if meta: row.update(meta)
        candidates.append(row)
        histories[name] = result[2]
        launches[name] = (q.copy(),v.copy())
        return row,result
    for kind,pattern in [('wave','converging_wave'),('collision','collision')]:
        for launch in P.velocity_candidates(kind,lattice.positions,lattice.box,lattice.target,lattice.target_axis,speed,
                                             lattice.predicted_kernel,focus_times_fs=times):
            score(pattern,launch['name'],launch['q'],launch['v'],{'proposed_focus_time_fs':launch.get('focus_time_fs')})
    direct_rows = []
    for launch in P.velocity_candidates('direct',lattice.positions,lattice.box,lattice.target,lattice.target_axis,speed,
                                        lattice.predicted_kernel,focus_times_fs=times):
        row,result = score('direct_optimized',launch['name'],launch['q'],launch['v'],{'proposed_focus_time_fs':launch.get('focus_time_fs')})
        direct_rows.append(row)
    best = max(direct_rows,key=lambda r:r['peak_relative_kinetic_eV'])
    q0,v0 = launches[best['candidate']]
    for k,v in enumerate(P.direct_refinements(v0,lattice.target,lattice.target_axis)):
        score('direct_optimized',f'direct_phase_{k}',q0,v)
    # Nonlinear amplitude AND phase refinement in a distributed, target-zero,
    # capped family. Four mixtures, not an assertion of a global optimum.
    others = [r for r in direct_rows if r is not best]
    if others:
        vb = launches[others[0]['candidate']][1]
        for c in (-.5,-.25,.25,.5):
            field = v0+c*vb
            v = P.capped_direct_velocity(field,lattice.target,speed,lattice.target_axis)
            score('direct_optimized',f'direct_mixture_{c:+.2f}',q0,v)
    # Exact reversal is a reachability control. Preparation is stopped at the
    # last supported point; its actual spread/preload decides headline eligibility.
    seed = next(P.velocity_candidates('reverse_seed',lattice.positions,lattice.box,lattice.target,lattice.target_axis,speed))
    prep = E.trajectory(lattice,seed['q'],seed['v'],DT,HORIZON/2)
    prep_time = prep[2][-1,0] if prep[3]<0 else max(0.,prep[3]-DT)
    if prep[3]>=0:
        prep = E.trajectory(lattice,seed['q'],seed['v'],DT,prep_time)
    qr,vr = prep[0],-prep[1]
    row,result = score('reverse_spreading','reverse_spreading',qr,vr,
        {'reverse_preparation_time_fs':float(prep_time),'reverse_preparation_supported':bool(prep[3]<0),
         'reverse_preparation_original_exit_fs':None if prep_time==HORIZON/2 else float(prep_time+DT)},
        horizon=max(DT,prep_time))
    selected = []
    for pattern in ('converging_wave','reverse_spreading','collision','direct_optimized'):
        options = [r for r in candidates if r['pattern']==pattern]
        eligible = [r for r in options if r['headline_eligible']]
        row = max(eligible or options,key=lambda r:r['peak_relative_kinetic_eV'])
        row['selected'] = True
        row['selection_candidates'] = len(options)
        trace = histories[row['candidate']]
        for scope in ('target','whole_cell_max'):
            if scope=='target':
                selected.append(row)
            else:
                q,v = launches[row['candidate']]
                # Reuse the actual selected trajectory, not a separately optimized cell maximum.
                result = (q,v,trace,-1. if row['first_site_or_bond_exit_fs'] is None else row['first_site_or_bond_exit_fs'],
                          {None:0,'site_radius':1,'bond_distance':2}[row['exit_kind']])
                other = E.summarize(lattice,q,v,result,energy,pattern,row['candidate'],scope)
                other.update(headline_eligible=row['headline_eligible'],selected=True,selection_candidates=len(options))
                selected.append(other)
    return selected,candidates,histories,launches,times


def campaign():
    started = time.monotonic()
    rows,candidates,proposals = [],[],[]
    trace_data, launch_data = {},{}
    for m in (2,4,8):
        lattice = E.Lattice(m)
        for energy in ENERGIES:
            rr,cc,hh,ll,tt = run_case(lattice,energy)
            rows.extend(rr);candidates.extend(cc)
            proposals.append(dict(n_D=lattice.n,energy_per_D_eV=energy,focus_times_fs=tt))
            for r in rr:
                if r['scope']!='target':continue
                key=f"n{lattice.n}_e{energy:g}_{r['pattern']}"
                trace_data[key]=hh[r['candidate']][::4]
                launch_data[key+'_q']=ll[r['candidate']][0]
                launch_data[key+'_v']=ll[r['candidate']][1]
            save_json(HERE/'results.json',dict(rows=rows,candidates=candidates,proposals=proposals,complete=False))
            print(json.dumps({'finished_N':lattice.n,'e':energy,'target_peaks':{r['pattern']:r['peak_relative_kinetic_eV'] for r in rr if r['scope']=='target'},'elapsed_s':time.monotonic()-started}),flush=True)
    # Only two matched harmonic consumers: low/high e at the largest cell;
    # same selected nonlinear launch, all-step identical support rules.
    controls=[]
    lattice=E.Lattice(8)
    for energy in (ENERGIES[0],ENERGIES[-1]):
        key=f'n2048_e{energy:g}_direct_optimized'
        q,v=launch_data[key+'_q'],launch_data[key+'_v']
        result=E.trajectory(lattice,q,v,DT,HORIZON,harmonic=True)
        row=E.summarize(lattice,q,v,result,energy,'matched_harmonic','same_nonlinear_selected_launch')
        row['hamiltonian']='harmonic';controls.append(row)
        trace_data[key+'_harmonic']=result[2][::4]
    # One, and only one, fixed optical amplitude-decay sensitivity.
    energy=ENERGIES[-1];key=f'n2048_e{energy:g}_direct_optimized'
    q,v=launch_data[key+'_q'],launch_data[key+'_v']
    result=E.trajectory(lattice,q,v,DT,HORIZON,damping=4136.03694)
    row=E.summarize(lattice,q,v,result,energy,'loss_sensitivity','same_selected_launch_tau4136fs')
    row['loss_tau_amplitude_fs']=4136.03694;controls.append(row)
    trace_data[key+'_loss']=result[2][::4]
    spectra=R.load_contact_spectra()
    readouts=[dict(n_D=r['n_D'],energy_per_D_eV=r['energy_per_D_eV'],pattern=r['pattern'],readout=R.conditional_readout(r,spectra)) for r in rows if r['scope']=='target']
    scalings=[]
    for energy in ENERGIES:
        for pattern in ('converging_wave','direct_optimized'):
            subset=[r for r in rows if r['energy_per_D_eV']==energy and r['pattern']==pattern and r['scope']=='target' and r['headline_eligible']]
            scalings.append(dict(energy_per_D_eV=energy,pattern=pattern,result=R.conditional_size_scaling(subset)))
    np.savez_compressed(HERE/'trajectories.npz',**trace_data)
    np.savez_compressed(HERE/'launches.npz',**launch_data)
    payload=dict(complete=True,dt_fs=DT,horizon_fs=HORIZON,mass_eV_fs2_A2=E.MASS,rows=rows,candidates=candidates,
        proposals=proposals,controls=controls,conditional_readouts=readouts,scalings=scalings,energy_floors=R.energy_budget_floors(),
        trace_columns=['time_fs','target_Krel_eV','target_Kradial_eV','target_R_A','max_site_radius_A','min_bond_R_A','total_energy_eV','potential_eV','wholecell_max_Krel_eV','target_inward_speed_A_fs','most_displaced_site_local_V_eV','most_displaced_site_harmonic_V_eV'],
        runtime_seconds=time.monotonic()-started)
    save_json(HERE/'results.json',payload)
    fields=list(rows[0])
    for r in rows:
        for k in r:
            if k not in fields:fields.append(k)
    with (HERE/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    print('COMPLETE',payload['runtime_seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true',required=True);parser.parse_args()
    campaign()
