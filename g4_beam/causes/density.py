"""Authorized shared1D and independent Pd/PdD2D density-profile likelihoods.

One physical consumer and its compiled stopping nodes; no per-material amplitude.
Reference remains m=1. All sheet changes and ratio-error assumptions are explicit.
"""
from pathlib import Path
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import json
import os
import sys
import time
import resource
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from g4_beam.merit.api import (Target,reference,bestfits,CAUSES,ENERGIES,FOUR,evaluate,evaluate_fit,
                               _reported,write_csv)
from g4_beam.merit.report import scaled_profile,sealed_target
from g4_beam.merit.responses import Response
from g4_beam.merit.refit import fit as fit_shape
from g4_beam.forward import model as M
from g4_beam.forward.profile_thermal import transport_stopping
HERE=Path(__file__).resolve().parent
DEPTH=10e-9
DATA=M.load_data(str(HERE.parent/'forward/errbars.csv'))
ONE=np.logspace(-1,2,13)
TWO=np.geomspace(.3,30.,7)
WORKERS=int(os.environ.get('G5_DENSITY_WORKERS','4'))
if not 1<=WORKERS<=4:raise ValueError('G5_DENSITY_WORKERS must be1–4')


@lru_cache(maxsize=40)
def changed_target(material,multiplier=1.,sheet=True):
    base=reference(material)
    if material not in ('Pd','PdD'):return base
    volume=scaled_profile(base.profiles,multiplier,depth=DEPTH)
    def profiles(e):
        layers,nads,state=volume(e)
        return layers,nads if sheet else 0.,dict(state,density_multiplier=multiplier,scaled_depth_nm=10.,sheet_retained=sheet)
    return Target(f'{material}_top10nm_x{multiplier}_sheet{sheet}',profiles,material,transport_stopping)


@lru_cache(maxsize=40)
def compiled(material,multiplier=1.,sheet=True):
    target=changed_target(material,multiplier,sheet)
    return {float(e):Response(float(e),target.profiles(float(e))) for e in ENERGIES}


def shape_residual(predictions):
    return np.concatenate([(np.log(np.interp(DATA[m][0],ENERGIES,predictions[m]))-np.log(DATA[m][1]))/DATA[m][2] for m in M.MATS])


def exact_predictions(targets,cause,fit):
    return {m:_reported(targets[m],ENERGIES,cause,fit) for m in M.MATS}


def density_summary(target):
    result={}
    vals=[];local=[]
    for e in ENERGIES:
        layers,_,_=target.profiles(float(e));a=np.asarray(layers,float)
        width=np.maximum(0,np.minimum(a[:,1],DEPTH)-a[:,0]);den=float(np.sum(width*a[:,6]))
        x=float(np.sum(width*a[:,3])/den);vals.append(x)
        local.append(float(np.max(a[width>0,3]/a[width>0,6])))
        if float(e) in FOUR:result[f'D_per_metal_top10_E{e:g}']=x
    result.update(top10_D_per_metal_min=min(vals),top10_D_per_metal_max=max(vals),max_local_D_per_metal=max(local))
    return result


def profile_job(job):
    cause,multiplier,sheet=job;start=time.monotonic();baseline=bestfits()[cause]
    targets={m:changed_target(m,multiplier if m in ('Pd','PdD') else 1.,sheet if m in ('Pd','PdD') else True) for m in M.MATS}
    frozen=exact_predictions(targets,cause,baseline)
    frozen_res=shape_residual(frozen);frozen_chi=float(frozen_res@frozen_res)
    if multiplier==1. and sheet:
        fitted=dict(baseline);fit_seconds=0.
    else:
        profiles={m:{float(e):targets[m].profiles(float(e)) for e in ENERGIES} for m in M.MATS}
        responses={m:compiled(m,multiplier if m in ('Pd','PdD') else 1.,sheet if m in ('Pd','PdD') else True) for m in M.MATS}
        if cause=='neutrals':
            for values in responses.values():
                for response in values.values():
                    if response.neutral_tables is None:response.build_neutral_table()
        fitted=fit_shape(cause,baseline['current_power'],profiles,baseline,responses,initial_params=baseline['params'])
        fitted['pass4_best_chi2']=baseline['pass4_best_chi2'];fitted['pass4_best_meas']=baseline['pass4_best_meas']
        fit_seconds=fitted['wall_s']
    loss=evaluate_fit(targets['PdD'],FOUR,cause,baseline)['physical_yield_per_C']/evaluate(reference(),FOUR,cause)['physical_yield_per_C']
    return dict(kind='density' if sheet else 'sheet0',cause=cause,multiplier=float(multiplier),sheet_retained=sheet,
                profile_chi2=float(fitted['chi2']),frozen_chi2=frozen_chi,current_power=baseline['current_power'],
                success=fitted['success'],fit=fitted,fit_wall_s=fit_seconds,wall_s=time.monotonic()-start,
                frozen_PdD_four_energy_yield_ratio=loss.tolist(),Pd=density_summary(targets['Pd']),PdD=density_summary(targets['PdD']),
                peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)




def joint_job(job):
    mp,md=map(float,job);start=time.monotonic();baseline=bestfits()['screeningrise']
    targets={m:changed_target(m,mp if m=='Pd' else md if m=='PdD' else 1.) for m in M.MATS}
    responses={m:compiled(m,mp if m=='Pd' else md if m=='PdD' else 1.) for m in M.MATS}
    from g4_beam.causes.likelihood import joint_fit
    result=joint_fit('screeningrise',targets,responses,baseline,DATA)
    result.update(Pd_multiplier=mp,PdD_multiplier=md,Pd=density_summary(targets['Pd']),
                  PdD=density_summary(targets['PdD']),wall_s=time.monotonic()-start,
                  peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return result


def flat_profile(row):
    result={k:v for k,v in row.items() if k not in ('fit','Pd','PdD','frozen_PdD_four_energy_yield_ratio')}
    result.update({m+'_'+k:v for m in ('Pd','PdD') for k,v in row[m].items()})
    result.update({f'frozen_yield_ratio_E{e:g}':v for e,v in zip(FOUR,row['frozen_PdD_four_energy_yield_ratio'])})
    return result


def flat_joint(row):
    result={k:v for k,v in row.items() if k not in ('fit','Pd','PdD','model_ratio_by_energy','SI1_E_cm_keV')}
    result.update({m+'_'+k:v for m in ('Pd','PdD') for k,v in row[m].items()})
    result.update({f'ratio_E{e:g}':v for e,v in zip(row['SI1_E_cm_keV'],row['model_ratio_by_energy'])})
    return result


def summaries(profiles,joint):
    one=[]
    for cause in CAUSES:
        rows=[r for r in profiles if r['cause']==cause and r['kind']=='density']
        minimum=min(rows,key=lambda r:r['profile_chi2'])
        result=dict(cause=cause,best_multiplier=minimum['multiplier'],minimum_chi2=minimum['profile_chi2'],
                    best_Pd=minimum['Pd'],best_PdD=minimum['PdD'])
        for delta in (1,4):
            selected=[r['multiplier'] for r in rows if r['profile_chi2']<=minimum['profile_chi2']+delta]
            result[f'delta{delta}_sampled_range']=[min(selected),max(selected)]
            result[f'delta{delta}_sampled_points']=selected
            result[f'delta{delta}_touches_scan_boundary']=min(selected)==float(ONE[0]) or max(selected)==float(ONE[-1])
        one.append(result)
    (HERE/'density1d_summary.json').write_text(json.dumps(one,indent=2))
    if not joint:return
    minimum=min(joint,key=lambda r:r['joint_chi2']);shape_min=min(r['shape_chi2'] for r in joint)
    region=[]
    for row in joint:
        row['shape_delta_chi2']=row['shape_chi2']-shape_min
        row['requested_region']=bool(row['shape_delta_chi2']<=4 and row['ratio_all_energies_in_measured_band'])
        row['geomean_region']=bool(row['shape_delta_chi2']<=4 and row['ratio_geomean_in_measured_band'])
        if row['requested_region']:region.append(flat_joint(row))
    flat=[flat_joint(r) for r in joint]
    write_csv(HERE/'density2d.csv',flat)
    with (HERE/'density2d_region.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(region)
    summary=dict(joint_minimum=flat_joint(minimum),shape_chi2_minimum=shape_min,requested_region_count=len(region),
                 region_criterion='shape component at joint optimum <= grid shape minimum+4 AND every one of7 modeled same-E ratios inside exact SI1 min/max band',
                 alternative_geomean_region_count=sum(r['geomean_region'] for r in joint),
                 measured_ratio_band=[float(min(M.counts_per_C_table('PdD')/M.counts_per_C_table('Pd'))),float(max(M.counts_per_C_table('PdD')/M.counts_per_C_table('Pd')))],
                 grid=TWO.tolist(),all_success=all(r['success'] for r in joint))
    (HERE/'density2d_summary.json').write_text(json.dumps(summary,indent=2))
    (HERE/'density2d_fits.json').write_text(json.dumps(joint,indent=2))


def main():
    start=time.monotonic();one_path=HERE/'density1d_fits.json';two_path=HERE/'density2d_fits.json'
    profiles=json.loads(one_path.read_text()) if one_path.exists() else []
    joint=json.loads(two_path.read_text()) if two_path.exists() else []
    if '--pilot-joint' in sys.argv:
        row=joint_job((float(TWO[0]),float(TWO[0])));joint.append(row)
        two_path.write_text(json.dumps(joint,indent=2));print(json.dumps(flat_joint(row),indent=2));return
    done={(r['cause'],r['multiplier'],r['sheet_retained']) for r in profiles}
    jobs=[(c,float(m),True) for c in CAUSES for m in ONE if (c,float(m),True) not in done]
    jobs += [(c,1.,False) for c in CAUSES if (c,1.,False) not in done]
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures={pool.submit(profile_job,j):j for j in jobs}
        for future in as_completed(futures):
            row=future.result();profiles.append(row);one_path.write_text(json.dumps(profiles,indent=2))
            write_csv(HERE/'density1d.csv',[flat_profile(r) for r in profiles])
            print('PROFILE',row['cause'],row['multiplier'],row['kind'],'chi2',row['profile_chi2'],'frozen',row['frozen_chi2'],'seconds',row['wall_s'],flush=True)
    done={(r['Pd_multiplier'],r['PdD_multiplier']) for r in joint}
    jobs=[(float(a),float(b)) for a in TWO for b in TWO if (float(a),float(b)) not in done]
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures={pool.submit(joint_job,j):j for j in jobs}
        for future in as_completed(futures):
            row=future.result();joint.append(row);two_path.write_text(json.dumps(joint,indent=2))
            print('JOINT',row['Pd_multiplier'],row['PdD_multiplier'],'shape',row['shape_chi2'],'ratio',row['ratio_chi2'],'seconds',row['wall_s'],flush=True)
    summaries(profiles,joint)
    print('DENSITY COMPLETE',len(profiles),len(joint),'wall_s',time.monotonic()-start,flush=True)

if __name__=='__main__':main()
