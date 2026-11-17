"""Six laboratory knobs and sealed-target spread, all through shared compare.

Beam-fraction/angle scans are immediate changes on frozen target profiles.
Current/temperature/time also have explicit transported-history responses.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import os
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from g4_beam.merit.api import Target,reference,bestfits,compare,CAUSES,ENERGIES,FOUR,write_csv
from g4_beam.merit.report import scaled_profile,sealed_target
from g4_beam.forward import model as M
from g4_beam.forward.profile_thermal import transport_stopping
from g4_beam.transport import transport as T
HERE=Path(__file__).resolve().parent
WORKERS=int(os.environ.get('G5_RESPONSE_WORKERS','1'))
if not 1<=WORKERS<=2:raise ValueError('G5_RESPONSE_WORKERS must be1–2')


def transport_condition(job):
    kind,value,e=job;case=dict(material='Pd',loaded=True,E_cm=e)
    base=reference().profiles(float(e))[2]
    if kind=='current':
        case.update(flux_scale=value,T=300+(base['T_K']-300)*value)
    elif kind=='temperature':case['dT']=value
    model=T.Model(case)
    te=np.unique(np.clip([600.,1800.,7200.,model.c['t_run']],0,model.c['t_run'])) if kind=='time' else np.array([model.c['t_run']])
    model.run_te(te);records=[]
    for i,t in enumerate(te):
        layers,nads,state=model.profile(i)
        records.append(dict(E_cm_keV=e,time_s=float(t),T_K=model.T,layers=np.asarray(layers).tolist(),nads_m2=float(nads),
                            state={k:np.asarray(v).tolist() for k,v in state.items()}))
    return dict(kind=kind,value=value,E_cm_keV=e,records=records)


def _profile(record):return record['layers'],record['nads_m2'],dict(record['state'],T_K=record['T_K'],beam_time_s=record['time_s'])


def ray_profile(angle_deg):
    cosine=float(np.cos(np.deg2rad(angle_deg)));boundary=bestfits()['surface']['params']['d']*1e-9/cosine
    def profile(e):
        layers,nads,state=reference().profiles(e);result=[]
        for row in layers:
            a=np.array(row,float);a[:2]/=cosine
            if a[0]<boundary<a[1]:
                left=a.copy();left[1]=boundary;right=a.copy();right[0]=boundary;result.extend([left.tolist(),right.tolist()])
            else:result.append(a.tolist())
        return result,nads/cosine,dict(state,ray_depth_coordinate=True,normal_depth_is_ray_depth_times=cosine,angle_deg=angle_deg)
    return profile


def main():
    start=time.monotonic();path=HERE/'knob_profiles.json'
    if path.exists():records=json.loads(path.read_text())
    else:
        jobs=[(kind,value,float(e)) for kind,values in [('current',(.5,2.)),('temperature',(-30.,30.)),('time',(0.,))] for value in values for e in ENERGIES]
        records=[]
        with ProcessPoolExecutor(max_workers=WORKERS) as pool:
            for i,row in enumerate(pool.map(transport_condition,jobs)):
                records.append(row)
                if (i+1)%14==0:
                    path.write_text(json.dumps(records,separators=(',',':')))
                    print('transported knob',i+1,len(jobs),row['kind'],row['value'],flush=True)
        path.write_text(json.dumps(records,separators=(',',':')))
    rows=[]
    baseline=list(csv.DictReader((HERE.parent/'merit/compare.csv').open()))
    baseline=[r for r in baseline if r['target']=='UC_Davis_transported_PdD']
    base={(r['cause'],float(r['E_cm_keV'])):float(r['physical_reactions_per_true_C']) for r in baseline}
    def append(target,knob,value,regime,*,only=None,current_factor=1.):
        for row in compare(target):
            if only is not None and row['cause'] not in only:continue
            row.update(knob=knob,value=value,regime=regime,
                       yield_per_true_C_ratio=row['physical_reactions_per_true_C']/base[row['cause'],row['E_cm_keV']],
                       physical_rate_ratio=current_factor*row['physical_reactions_per_true_C']/base[row['cause'],row['E_cm_keV']],
                       baseline='sameE nominal run-end transported PdD; neutral extracted flux held fixed')
            rows.append(row)
        print('compared',knob,value,regime,flush=True)
    current=lambda e:M.I_cup(e)*1e-6/M.collected_fraction(e,1,0)
    for factor in (.5,2.):
        target=Target('instant_current_x'+str(factor),reference().profiles,stopping_fn=transport_stopping,
                      measurement={'current_A':lambda e,f=factor:current(e)*f})
        append(target,'current',factor,'immediate frozen profile; delivered charged beam changed at fixed extracted neutral flux',current_factor=factor)
    for fraction in (0.,.5,.85):
        target=Target('D2_current_fraction_'+str(fraction),reference().profiles,stopping_fn=transport_stopping,
                      measurement={'D1_fraction':1-fraction,'D2_fraction':fraction})
        append(target,'D2_fraction',fraction,'immediate frozen target; fraction is electrical current, not D atom fraction')
    for angle in (30.,60.):
        profile=ray_profile(angle);target=Target('angle_'+str(angle),profile,stopping_fn=transport_stopping)
        append(target,'angle_deg_from_normal',angle,'immediate frozen profile; straight ray; normal screening thickness preserved; no reflection/straggling',only=('screeningrise','neutrals'))
        fit=bestfits()['surface']['params'];d=fit['d']*1e-9;extended=d/np.cos(np.deg2rad(angle))
        delta=fit['Ues_PdD']-fit['Ue_PdD']
        # Preserve the fitted NORMAL-depth boundary after converting layers to
        # ray-depth. The exact extended boundary was inserted in ray_profile.
        extra=lambda e,z,dd=d,xx=extended,du=delta:du if dd<=z<xx else 0.
        surface=Target('angle_'+str(angle)+'_surface_geometry',profile,stopping_fn=transport_stopping,Ue_shift_keV=extra)
        append(surface,'angle_deg_from_normal',angle,'immediate frozen profile; straight ray; normal screening thickness preserved; no reflection/straggling',only=('surface',))
    for factor in (.1,3.,50.):
        target=Target('top10nm_D_x'+str(factor),scaled_profile(reference().profiles,factor,depth=10e-9),stopping_fn=transport_stopping)
        append(target,'near_surface_D_factor',factor,'top10nm volumetricD only; adsorbate unchanged; physical D stopping')
    for kind,values in [('current',(.5,2.)),('temperature',(-30.,30.))]:
        for value in values:
            selected=[r for r in records if r['kind']==kind and r['value']==value]
            profiles={float(r['E_cm_keV']):_profile(r['records'][-1]) for r in selected}
            measurement={'current_A':lambda e,f=value:current(e)*f} if kind=='current' else {}
            target=Target(f'transported_{kind}_{value}',profiles,stopping_fn=transport_stopping,measurement=measurement)
            append(target,kind,value,'actual Model nominal run history; retained pass4 surface estimates',current_factor=value if kind=='current' else 1.)
    for time_s in (600.,1800.,7200.):
        selected=[r for r in records if r['kind']=='time']
        profiles={float(r['E_cm_keV']):_profile(next(p for p in r['records'] if p['time_s']==time_s)) for r in selected}
        append(Target('beam_time_'+str(time_s),profiles,stopping_fn=transport_stopping),'time_s',time_s,'actual Model after same3600s preload; reference is nominal run end')
    write_csv(HERE/'predictions.csv',rows)
    # Sealed hydride: report BOTH sheet hypotheses, not an implicit retained ML.
    sealrows=[]
    density=sealed_target()
    for keep in (True,False):
        def profiles(e,kept=keep):
            layers,nads,state=density.profiles(e)
            return layers,nads if kept else 0.,dict(state,adsorbed_sheet_kept=kept)
        target=Target('sealed0.7_'+('retained_sheet' if keep else 'no_sheet'),profiles,stopping_fn=transport_stopping)
        sealrows.extend(compare(target))
    write_csv(HERE/'sealed_compare.csv',sealrows)
    # Rank LAB-measurable dose-normalized sameE ratios, retaining one winner per
    # knob so three rows are not simply the same knob at adjacent energies.
    groups={}
    for row in rows:
        key=(row['knob'],row['value'],row['regime'],row['E_cm_keV'])
        groups.setdefault(key,[]).append(row)
    ranked=[]
    for key,rr in groups.items():
        if len(rr)!=3 or any(r['yield_per_true_C_ratio']<=0 for r in rr):continue
        vals={r['cause']:r['yield_per_true_C_ratio'] for r in rr}
        ranked.append(dict(knob=key[0],value=key[1],regime=key[2],E_cm_keV=key[3],
                           separation_factor=max(vals.values())/min(vals.values()),**{c+'_ratio':vals[c] for c in CAUSES},
                           observable='proton counts per independently measured incident charged C, changed/base at sameE and calibrated detector geometry'))
    chosen=[]
    for row in sorted(ranked,key=lambda r:r['separation_factor'],reverse=True):
        if row['knob'] not in {r['knob'] for r in chosen}:chosen.append(row)
        if len(chosen)==3:break
    write_csv(HERE/'laboratory_ratios.csv',chosen)
    # Reuse existing5nm density response; no extra profile fit or simulation.
    with (HERE.parent/'merit/sensitivities.csv').open() as f:
        reused=[r for r in csv.DictReader(f) if r['density_factor']]
    write_csv(HERE/'top5nm_reused_response.csv',reused)
    status=dict(completed_knobs=sorted({r['knob'] for r in rows}),rows=len(rows),wall_s=time.monotonic()-start,
                transported_profile_cases=len(records),worker_threads=1,workers=WORKERS,
                source='853c377 plus accepted merit archive; no Zenodo access',all_designs_through_compare=True)
    (HERE/'response_status.json').write_text(json.dumps(status,indent=2));print(json.dumps(status,indent=2),flush=True)

if __name__=='__main__':main()
