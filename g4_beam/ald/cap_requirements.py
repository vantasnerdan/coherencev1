"""TiD-counts calibration and retained20nm-TiD2 sputter requirements.
Reuses the resolved physical-yield consumer; no fitted parameters are changed.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from copy import deepcopy
import csv,gzip,json,os,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from g4_beam.ald import dense as D
from g4_beam.ald.designs import HERE,MATERIALS,erosion,add_film,report
from g4_beam.ald.finalize import TiD_calibration,calibrated
from g4_beam.merit.api import ENERGIES,FOUR,reference,write_csv
WORKERS=int(os.environ.get('G5_SPUTTER_WORKERS','4'))
if not 1<=WORKERS<=4:raise ValueError('G5_SPUTTER_WORKERS must be1–4')
BASE_DESIGN=next(d for d in D.designs() if d['name']=='TiD2_cap_20.0nm_ret1.0_erosion0.1')
CACHE=HERE/'cap_sputter_cache.json.gz'
NOMINAL=HERE/'nominal_resolved.json.gz'


def variant(factor):
    design=dict(BASE_DESIGN,name=f'TiD2_cap_20.0nm_ret1.0_erosion{factor}',
                scenario=f'D_retained_fraction1.0;sputter_x{factor}')
    return design


def nominal_design():
    design=deepcopy(variant(1.0));mat=MATERIALS[design['material']]
    for e in ENERGIES:
        base=reference().profiles(float(e));rate,life=erosion(mat,float(e),20.)
        remaining=max(0.,20.-rate*base[2]['beam_time_s']*1e9)
        design['profiles'][str(float(e))]=add_film(base,mat,remaining,2.,1.)
        design['lifetimes'][str(float(e))]['remaining_film_nm']=remaining
    return design


def key(e,factor):return f'{float(e):.17g}|{float(factor):.17g}'


def evaluate_point(job):
    e,factor=job
    result=D.resolve_one(variant(factor),energies=[e])
    if factor==0:
        for row in result['integrals']:
            row['clearance_s']=None
            row['clearance_status']='noneroding control; no finite clearance'
    return e,factor,result


def batch(jobs,points,nominal):
    jobs=sorted(set(jobs)-{(float(r['E_cm_keV']),float(r['sputter_factor'])) for r in points.values()})
    if not jobs:return
    for e in sorted({e for e,factor in jobs}):
        D.prepare_reference([variant(factor) for ee,factor in jobs if ee==e],energies=[e])
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures={pool.submit(evaluate_point,job):job for job in jobs};errors=[]
        for future in as_completed(futures):
            try:e,factor,result=future.result()
            except Exception as error:
                errors.append((futures[future],repr(error)));print('SPUTTER_FAILED',errors[-1],flush=True);continue
            points[key(e,factor)]=dict(E_cm_keV=e,sputter_factor=factor,integrals=result['integrals'],
                                      wall_s=result['wall_s'],source='resolved physical age kernels')
            if factor==1.:
                nominal[str(e)]=result;D.save(NOMINAL,nominal)
            D.save(CACHE,points)
            print('SPUTTER_POINT',e,factor,'wall_s',result['wall_s'],flush=True)
    if errors:raise RuntimeError('Valid sputter points checkpointed; failed points remain missing: '+repr(errors))


def response(points):
    rows=[]
    for point in points.values():
        for item in point['integrals']:
            rows.append(dict(item,sputter_factor=point['sputter_factor'],source=point['source']))
    return rows


def brackets(points,calibration):
    rows=response(points);out=[]
    for e in FOUR:
        for scenario,cause,fit in D.SCENARIOS:
            observed=calibration.get((scenario,cause,float(e)))
            for reading in ('raw_model','SI1_TiD_counts_calibrated'):
                if reading!='raw_model' and observed is None:
                    out.append(dict(reference_scenario=scenario,cause=cause,E_cm_keV=float(e),criterion=reading,
                                    status='no same-energy SI1 anchor',largest_verified_sputter_factor=None,first_failing_factor=None))
                    continue
                factor=float(observed['model_over_observed']) if reading!='raw_model' else 1.
                matching=[r for r in rows if r['reference_scenario']==scenario and r['cause']==cause and float(r['E_cm_keV'])==float(e)]
                limit=next((float(r['integrated_merit'])/factor for r in matching if r['sputter_factor']==0),None)
                cases=sorted((r for r in matching if r['sputter_factor']>0),key=lambda r:r['sputter_factor'])
                successful=[r for r in cases if float(r['integrated_merit'])/factor>=10.]
                base=dict(reference_scenario=scenario,cause=cause,E_cm_keV=float(e),criterion=reading,
                          threshold_merit=10.,TiD_counts_model_over_observed_factor=factor,
                          search_min_factor=min(r['sputter_factor'] for r in cases),search_max_factor=max(r['sputter_factor'] for r in cases),
                          noneroding_integrated_merit=limit)
                if not successful:
                    out.append(dict(base,status='no sampled success; noneroding limit below10' if limit is not None and limit<10 else 'no >=10 point in positive sampled range',largest_verified_sputter_factor=None,
                                    first_failing_factor=None,max_sampled_merit=max(float(r['integrated_merit'])/factor for r in cases)))
                    continue
                low=successful[-1];failed=[r for r in cases if r['sputter_factor']>low['sputter_factor'] and float(r['integrated_merit'])/factor<10.]
                high=failed[0] if failed else None
                out.append(dict(base,status='rightmost sampled crossing bracket' if high else 'upper scan boundary still >=10',
                                largest_verified_sputter_factor=low['sputter_factor'],merit_at_verified_factor=float(low['integrated_merit'])/factor,
                                first_failing_factor=high['sputter_factor'] if high else None,
                                merit_at_failing_factor=float(high['integrated_merit'])/factor if high else None,
                                bracket_factor_ratio=high['sputter_factor']/low['sputter_factor'] if high else None))
    return out


def main():
    start=time.monotonic();points=D.load(CACHE,{});nominal=D.load(NOMINAL,{})
    existing=D.load(HERE/'resolved_cache.json.gz',[])
    for factor in (.1,10.):
        saved=next(r for r in existing if r['target']==f'TiD2_cap_20.0nm_ret1.0_erosion{factor}')
        for e in FOUR:
            k=key(e,factor)
            if k not in points:points[k]=dict(E_cm_keV=float(e),sputter_factor=factor,
                integrals=[r for r in saved['integrals'] if r['E_cm_keV']==float(e)],wall_s=0.,source='previous selected resolved trajectory, not replayed')
    D.save(CACHE,points)
    if '--pilot' in sys.argv:
        batch([(float(e),1.) for e in FOUR],points,nominal)
        design=nominal_design();D.save(HERE/'nominal_design.json.gz',[design]);report([design],'nominal_cap')
        summary=dict(total_wall_s=time.monotonic()-start,workers=WORKERS,threads=1,
                     nominal_point_walls_s={e:r['wall_s'] for e,r in nominal.items()},exit_reached=True)
        (HERE/'cap_requirements_pilot.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True);return
    calibration=TiD_calibration()
    batch([(float(e),f) for e in FOUR for f in (0.,.03,.1,.3,1.,3.,10.,30.)],points,nominal)
    for iteration in range(12):
        requested=set()
        for row in brackets(points,calibration):
            lo=row['largest_verified_sputter_factor'];hi=row['first_failing_factor']
            if lo is None:
                if row.get('noneroding_integrated_merit',0.)>10:
                    minimum=row['search_min_factor']
                    if minimum<=1e-8:raise RuntimeError('Noneroding limit succeeds but no positive factor located; numerical requirement unresolved')
                    requested.add((row['E_cm_keV'],minimum/10))
                continue
            if hi is None:
                if lo>=1000.:raise RuntimeError('Sputter upper boundary remains successful; no finite threshold delivered')
                requested.add((row['E_cm_keV'],min(1000.,3*lo)))
            elif hi/lo>1.10:requested.add((row['E_cm_keV'],float(f'{np.sqrt(lo*hi):.12g}')))
        if not requested:break
        batch(requested,points,nominal)
    else:raise RuntimeError('Sputter brackets did not reach the declared10% factor resolution')
    threshold=brackets(points,calibration)
    for row in threshold:
        row['interpretation']='rightmost sampled crossing; finite scan, not a global monotonicity proof or confidence enclosure'
        row['headline']=row['criterion']==('SI1_TiD_counts_calibrated' if row['E_cm_keV']<=3 else 'raw_model')
    samples=response(points)
    for row in samples:
        cal=calibration.get((row['reference_scenario'],row['cause'],float(row['E_cm_keV'])))
        row['TiD_counts_calibrated_integrated_merit']=calibrated(row['integrated_merit'],cal)
    write_csv(HERE/'cap_sputter_response.csv',samples);write_csv(HERE/'cap_sputter_thresholds.csv',threshold)
    nominal_rows=[dict(r) for result in nominal.values() for r in result['integrals']]
    for row in nominal_rows:
        cal=calibration.get((row['reference_scenario'],row['cause'],float(row['E_cm_keV'])))
        row['TiD_counts_calibrated_integrated_merit']=calibrated(row['integrated_merit'],cal)
    write_csv(HERE/'nominal_cap_integrals.csv',nominal_rows)
    summary=dict(wall_s=time.monotonic()-start,single_energy_points=len(points),threshold_rows=len(threshold),
                 workers=WORKERS,threads=1,bracket_factor_resolution=1.10,source='same frozen fits and explicit physical age kernels')
    (HERE/'cap_requirements_status.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
