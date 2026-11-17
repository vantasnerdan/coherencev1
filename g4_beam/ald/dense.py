"""Selected, resolved ALD1/2 trajectories; no endpoint-as-run-integral shortcut.

Pilot includes reference transport, all nuclear age evaluations and cache I/O.
Only selected designs acquire a published integral. Coarse/reset estimates are
selection diagnostics, not bounds or a replacement for the resolved result.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import csv,gzip,json,os,re,resource,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from g4_beam.merit.api import Target,FOUR,CAUSES,bestfits,evaluate_fit,write_csv
from g4_beam.forward import model as M
from g4_beam.forward.profile_thermal import transport_stopping
from g4_beam.causes.references import (ZERO_SHEET_SURFACE,zero_sheet_fit,without_exposed_pd_sheet,pd_sheet_identity)
from g4_beam.ald.designs import (HERE,MATERIALS,T,profile_tuple,compact_profile,erosion,add_film,
                                seal_model,observed_run,design_target)
WORKERS=int(os.environ.get('G5_DENSE_WORKERS','4'))
if not 1<=WORKERS<=4:raise ValueError('G5_DENSE_WORKERS must be1–4')
SCENARIOS=[('canonical_full_sheet',c,bestfits()[c]) for c in CAUSES]+[(ZERO_SHEET_SURFACE,'surface',zero_sheet_fit())]
BASE=None
SEALS=None
COARSE_TIMES=None


def load(path,default):
    return json.loads(gzip.decompress(path.read_bytes())) if path.exists() else default


def save(path,value):
    temporary=path.with_suffix('.tmp')
    temporary.write_bytes(gzip.compress(json.dumps(value,separators=(',',':')).encode(),compresslevel=1))
    temporary.replace(path)


def designs():return load(HERE/'early_profiles.json.gz',[])


def parameters(design):
    if design['family']==1:
        match=re.fullmatch(r'permeability_x([0-9.e+-]+)',design['scenario'])
        if not match:raise ValueError(design['scenario'])
        return dict(permeability=float(match[1]))
    if design['family']==2:
        match=re.fullmatch(r'D_retained_fraction([0-9.e+-]+);sputter_x([0-9.e+-]+)',design['scenario'])
        if not match:raise ValueError(design['scenario'])
        return dict(retention=float(match[1]),sputter=float(match[2]))
    raise ValueError('This authorized selected rerun set consists of existing ALD1/2 designs')


def duration(design,e):
    profile=design['profiles'][str(float(e))]
    return float(profile['time_s'] if design['family']==1 else profile['state']['beam_time_s'])


def clearance(design,e):
    rate,life=erosion(MATERIALS[design['material']],e,design['thickness_nm'])
    d=design['thickness_nm']
    if design['family']==2:
        factor=parameters(design)['sputter']
        return np.inf if factor==0 else life/factor
    return ((d-.001)+.02*np.log(d/.001))*1e-9/rate


def coarse_times(design,e):
    end=duration(design,e)
    if design['family']==1:return np.array([p['time_s'] for p in design['erosion_snapshots'][str(float(e))]],float)
    life=clearance(design,e)
    requested=np.unique(np.clip([0.,600.,1800.,life/2,.95*life,end],0,end))
    if COARSE_TIMES is not None:
        available=COARSE_TIMES[str(float(e))]
        requested=np.unique([min(available,key=lambda t:abs(t-wanted)) for wanted in requested])
    return requested


def resolved_times(design,e):
    end=duration(design,e);clear=clearance(design,e)
    if np.isinf(clear):return np.unique(np.r_[np.arange(0.,end,300.),0.,end])
    dense_end=min(end,clear+3600.)
    points=list(np.arange(0.,dense_end,300.))+[t for t in np.arange(0.,end,1800.) if t>=dense_end]+[0.,dense_end,end]
    points+=coarse_times(design,e).tolist()
    rate,_=erosion(MATERIALS[design['material']],e,design['thickness_nm']);d=design['thickness_nm']
    for residual in (.1,.05,.02,.005,.001):
        if residual>=d:continue
        age=(d-residual)*1e-9/(rate*parameters(design)['sputter']) if design['family']==2 else ((d-residual)+.02*np.log(d/residual))*1e-9/rate
        points.append(age)
    points += [clear+dt for dt in (-120.,-60.,-30.,0.,30.,60.,120.,3600.)]
    return np.unique(np.clip(points,0.,end))


def scenario_key(name,cause):return name+':'+cause


def raw_values(target,e):
    values={}
    for name,cause,fit in SCENARIOS:
        projected=without_exposed_pd_sheet(target) if name==ZERO_SHEET_SURFACE else target
        values[scenario_key(name,cause)]=float(evaluate_fit(projected,[e],cause,fit,reference_E_cm_keV=None)['physical_yield_per_C'][0])
    return values


def prepare_reference(all_designs,*,energies=FOUR):
    global BASE
    path=HERE/'dense_reference.json.gz';BASE=load(path,{})
    start=time.monotonic()
    for e in energies:
        key=str(float(e));records=BASE.setdefault(key,{})
        desired=np.unique(np.concatenate([resolved_times(d,float(e)) for d in all_designs]))
        times=np.array([t for t in desired if str(float(t)) not in records])
        if not len(times):continue
        model=T.Model(dict(material='Pd',loaded=True,E_cm=float(e)))
        observed_run(model,times);last_save=time.monotonic()
        for i,age in enumerate(model.t_beam):
            profile=compact_profile(model.profile(i));profile['time_s']=float(age)
            target=Target('same_age_transported_PdD',{float(e):profile_tuple(profile)},stopping_fn=transport_stopping)
            records[str(float(age))]=dict(profile=profile,yields=raw_values(target,float(e)))
            if time.monotonic()-last_save>=60:
                save(path,BASE);last_save=time.monotonic()
                print('DENSE_REFERENCE_CHECKPOINT',e,'saved_age_kernels',len(records),flush=True)
        save(path,BASE)
        print('DENSE_REFERENCE',e,'new_ages',len(times),'saved_ages',len(records),'cumulative_wall_s',time.monotonic()-start,flush=True)
    return time.monotonic()-start


def prepare_coarse_reference():
    """Use actual saved source ages; never interpolate a new transport state."""
    global BASE,COARSE_TIMES
    path=HERE/'coarse_reference.json.gz';BASE=load(path,{})
    old=load(HERE/'matched_reference_history.json.gz',{})
    completed_dense=load(HERE/'dense_reference.json.gz',{})
    for e in FOUR:
        key=str(float(e));records=BASE.setdefault(key,{})
        for age,profile in old[key].items():
            if age in records:continue
            if age in completed_dense.get(key,{}):records[age]=completed_dense[key][age];continue
            target=Target('saved_same_age_PdD',{float(e):profile_tuple(profile)},stopping_fn=transport_stopping)
            records[age]=dict(profile=profile,yields=raw_values(target,float(e)))
        save(path,BASE)
    COARSE_TIMES={e:sorted(map(float,rows)) for e,rows in BASE.items()}
    return BASE


def base_at(e,age):return BASE[str(float(e))][str(float(age))]


def cap_profile(design,e,age):
    p=parameters(design);rate,_=erosion(MATERIALS[design['material']],e,design['thickness_nm'])
    remaining=max(0.,design['thickness_nm']-rate*p['sputter']*age*1e9)
    return add_film(profile_tuple(base_at(e,age)['profile']),MATERIALS[design['material']],remaining,2.,p['retention'])


def age_target(design,e,age,profile):
    covered={float(e):profile_tuple(profile)}
    opened={float(e):profile_tuple(base_at(e,age)['profile'])}
    return design_target(design,covered,opened)


def integral(times,numerator,denominator):
    t=np.asarray(times,float);n=np.asarray(numerator,float);d=np.asarray(denominator,float)
    if t[0]!=0 or np.any(d<=0) or np.any(np.diff(t)<=0):raise ValueError('Incomplete/invalid matched trajectory')
    ratio=float(np.trapezoid(n,t)/np.trapezoid(d,t));h=n-10*d;above=0.
    for i,dt in enumerate(np.diff(t)):
        a,b=h[i],h[i+1]
        if a>=0 and b>=0:above+=dt
        elif a<0 and b<0:continue
        elif a>=0:above+=dt*a/(a-b)
        else:above+=dt*b/(b-a)
    return dict(integrated_merit=ratio,end_merit=float(n[-1]/d[-1]),peak_merit=float(np.max(n/d)),
                time_merit_ge10_s=float(above),integral_design_Y_s_per_C=float(np.trapezoid(n,t)),
                integral_reference_Y_s_per_C=float(np.trapezoid(d,t)))


def coarse_one(design):
    start=time.monotonic();cases=[]
    for e in FOUR:
        e=float(e);times=coarse_times(design,e);values=[]
        for age in times:
            profile=next(p for p in design['erosion_snapshots'][str(e)] if p['time_s']==age) if design['family']==1 else cap_profile(design,e,float(age))
            values.append(raw_values(age_target(design,e,float(age),profile),e))
        for name,cause,fit in SCENARIOS:
            key=scenario_key(name,cause);num=np.array([v[key] for v in values]);den=np.array([base_at(e,float(t))['yields'][key] for t in times])
            coarse=integral(times,num,den);clear=clearance(design,e);rt=list(times);rn=list(num);rd=list(den)
            if clear<times[-1]:
                index=int(np.searchsorted(rt,clear))
                reference=float(np.interp(clear,times,den))
                if index<len(rt) and rt[index]==clear:rn[index]=reference
                else:rt.insert(index,clear);rn.insert(index,reference);rd.insert(index,reference)
                rn=[b if t>=clear else a for t,a,b in zip(rt,rn,rd)]
            reset=integral(rt,rn,rd)
            cases.append(dict(reference_scenario=name,cause=cause,E_cm_keV=e,coarse_integrated=coarse['integrated_merit'],
                              reset_integrated=reset['integrated_merit'],sampled_peak=coarse['peak_merit'],
                              disagree_on_10x=bool((coarse['integrated_merit']>=10)!=(reset['integrated_merit']>=10))))
    return dict(target=design['name'],max_sampled_peak=max(r['sampled_peak'] for r in cases),
                tie_break_integral=max(r['coarse_integrated'] for r in cases),disagreement=any(r['disagree_on_10x'] for r in cases),
                cases=cases,wall_s=time.monotonic()-start)


def collect_results(worker,jobs,path,results,label):
    errors=[]
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures={pool.submit(worker,design):design['name'] for design in jobs}
        for future in as_completed(futures):
            try:row=future.result()
            except Exception as error:
                errors.append((futures[future],repr(error)));print('CASE_FAILED',errors[-1],flush=True)
                continue
            results.append(row);save(path,results)
            print(label,len(results),row['target'],'wall_s',row['wall_s'],flush=True)
    if errors:raise RuntimeError('Failed cases remain missing; valid cases checkpointed: '+repr(errors))


def select(all_designs):
    path=HERE/'dense_selection_cache.json.gz';results=load(path,[]);done={r['target'] for r in results}
    collect_results(coarse_one,[d for d in all_designs if d['name'] not in done],path,results,'COARSE_SELECTED_INPUT')
    top={r['target'] for r in sorted(results,key=lambda r:(-r['max_sampled_peak'],-r['tie_break_integral'],r['target']))[:10]}
    canonical_peak={}
    for row in csv.DictReader((HERE/'early_erosion.csv').open()):
        canonical_peak[row['target']]=max(canonical_peak.get(row['target'],0.),float(row['merit_PdD']))
    canonical_top=set(sorted(canonical_peak,key=lambda name:(-canonical_peak[name],name))[:10])
    selected={r['target'] for r in results if r['target'] in canonical_top or r['target'] in top or r['disagreement']}
    cap_recipes={(d['material'],d['thickness_nm']) for d in all_designs if d['name'] in selected and d['family']==2}
    companions={d['name'] for d in all_designs if d['family']==2 and (d['material'],d['thickness_nm']) in cap_recipes}
    selected|=companions
    write_csv(HERE/'selected_designs.csv',[dict(target=r['target'],sampled_peak_merit=r['max_sampled_peak'],
              canonical_sampled_peak_merit=canonical_peak.get(r['target']),
              selected_by_canonical_top10_peak=r['target'] in canonical_top,selected_by_secondary_joint_peak=r['target'] in top,
              selected_by_10x_disagreement=r['disagreement'],selected_by_cap_endpoint_companion=r['target'] in companions,
              selected=r['target'] in selected) for r in results])
    return [d for d in all_designs if d['name'] in selected],{r['target']:r for r in results if r['target'] in selected}


def seal_key(design,e):return '|'.join(map(str,(design['material'],design['thickness_nm'],parameters(design)['permeability'],float(e))))


def prepare_seals(selected):
    global SEALS
    path=HERE/'dense_seal_profiles.json.gz';SEALS=load(path,{})
    for design in selected:
        if design['family']!=1:continue
        for e in FOUR:
            key=seal_key(design,e)
            if key in SEALS:continue
            model,rate,life=seal_model((design['material'],design['thickness_nm'],parameters(design)['permeability'],float(e),None,1.))
            times=resolved_times(design,float(e));observed_run(model,times)
            SEALS[key]={str(float(t)):dict(compact_profile(model.profile(i)),time_s=float(t)) for i,t in enumerate(model.t_beam)}
            save(path,SEALS);print('DENSE_SEAL',key,'ages',len(times),flush=True)


def resolve_one(design,*,energies=FOUR):
    start=time.monotonic();curves=[];summaries=[]
    for e in energies:
        e=float(e);times=resolved_times(design,e);by_scenario={scenario_key(n,c):[] for n,c,f in SCENARIOS}
        for age in times:
            age=float(age)
            profile=SEALS[seal_key(design,e)][str(age)] if design['family']==1 else cap_profile(design,e,age)
            target=age_target(design,e,age,profile);values=raw_values(target,e)
            for name,cause,fit in SCENARIOS:
                key=scenario_key(name,cause);num=values[key];den=base_at(e,age)['yields'][key]
                projected=without_exposed_pd_sheet(target) if name==ZERO_SHEET_SURFACE else target
                row=dict(target=design['name'],family=design['family'],reference_scenario=name,cause=cause,E_cm_keV=e,time_s=age,
                         physical_design_per_true_C=num,physical_reference_per_true_C=den,instantaneous_merit=num/den,
                         **pd_sheet_identity(projected,e))
                curves.append(row);by_scenario[key].append(row)
        for name,cause,fit in SCENARIOS:
            rr=by_scenario[scenario_key(name,cause)];result=integral(times,[r['physical_design_per_true_C'] for r in rr],[r['physical_reference_per_true_C'] for r in rr])
            current=float(M.I_cup(e)*1e-6/M.collected_fraction(np.array([e]),fit['b_hi'],fit['b_lo'])[0])
            summaries.append(dict(target=design['name'],family=design['family'],reference_scenario=name,cause=cause,E_cm_keV=e,
                    **result,assumed_incident_current_A=current,assumed_run_charge_C=current*times[-1],
                    integrated_design_reactions=current*result['integral_design_Y_s_per_C'],
                    integrated_reference_reactions=current*result['integral_reference_Y_s_per_C'],run_end_s=float(times[-1]),
                    ages=len(times),clearance_s=clearance(design,e),interpolation='linear physical yields; ratio of separate integrals',
                    clock='retained provided-IV-as-incident Model clock; matched same-age beam',resolved=True))
    return dict(target=design['name'],curves=curves,integrals=summaries,wall_s=time.monotonic()-start,
                peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


def export(results,checks):
    write_csv(HERE/'resolved_yields.csv',[r for result in results for r in result['curves']])
    write_csv(HERE/'resolved_integrals.csv',[r for result in results for r in result['integrals']])
    rows=[]
    for result in results:
        for item in checks.get(result['target'],{}).get('cases',[]):rows.append(dict(target=result['target'],**item,role='check only; resolved integral is primary'))
    if rows:write_csv(HERE/'selected_integral_checks.csv',rows)


def main():
    global BASE,SEALS
    start=time.monotonic();all_designs=designs()
    if len(all_designs)!=132:raise ValueError('Expected the132 actually completed early designs')
    if '--select' in sys.argv:
        prepare_coarse_reference();selected,checks=select(all_designs)
        existing=load(HERE/'dense_reference.json.gz',{})
        coverage={}
        for e in FOUR:
            required=set(map(float,np.unique(np.concatenate([resolved_times(d,float(e)) for d in selected]))))
            saved=set(map(float,existing.get(str(float(e)),{})))
            coverage[str(float(e))]=dict(required_ages=len(required),reused_age_kernels=len(required&saved),missing_age_kernels=len(required-saved))
        status=dict(selected_designs=len(selected),coverage=coverage,wall_s=time.monotonic()-start,
                    selection_states='actual saved Model ages; cap query snapped to an available age, not an interpolated profile')
        (HERE/'selection_status.json').write_text(json.dumps(status,indent=2));print(json.dumps(status,indent=2),flush=True);return
    cache=HERE/'resolved_cache.json.gz';results=load(cache,[])
    if '--pilot' in sys.argv or '--pilot-seal' in sys.argv:
        name='Al2O3_seal_2.0nm_perm0.01_coverage1.0' if '--pilot-seal' in sys.argv else 'TiD2_cap_20.0nm_ret1.0_erosion0.1'
        chosen=next(d for d in all_designs if d['name']==name)
        reference_wall=prepare_reference([chosen])
        if chosen['family']==1:prepare_seals([chosen])
        if not any(r['target']==name for r in results):results.append(resolve_one(chosen));save(cache,results)
        export(results,{})
        status=dict(pilot_target=name,total_wall_s=time.monotonic()-start,reference_preparation_wall_s=reference_wall,
                    full_design_forward_IO_wall_s=results[-1]['wall_s'],ages=sum(r['ages'] for r in results[-1]['integrals'])//4,
                    threads=1,workers=1,exit_reached=True)
        status_path=HERE/('dense_pilot_seal.json' if '--pilot-seal' in sys.argv else 'dense_pilot.json')
        status_path.write_text(json.dumps(status,indent=2));print(json.dumps(status,indent=2),flush=True);return
    prepare_coarse_reference();selected,checks=select(all_designs)
    reference_wall=prepare_reference(selected);prepare_seals(selected);done={r['target'] for r in results}
    collect_results(resolve_one,[d for d in selected if d['name'] not in done],cache,results,'RESOLVED')
    # A reused pilot is retained as an observed resolved case even if ties later
    # rank it outside the selected10; no valid physical trajectory is replayed.
    export(results,checks)
    status=dict(selected_designs=len(selected),resolved_designs=len(results),wall_s=time.monotonic()-start,workers=WORKERS,threads=1,
                primary='ratio of physical charge integrals',nonrerun_integral_estimates_published=False)
    (HERE/'dense_status.json').write_text(json.dumps(status,indent=2));print(json.dumps(status,indent=2),flush=True)

if __name__=='__main__':main()
