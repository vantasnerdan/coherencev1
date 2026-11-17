"""Rerun changed Pd rows; retain the original estimates and unchanged Ti rows."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import time

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'screen'))
from screen_io import manifest,require,write_json
from knob_cases import build_jobs


def run_job(runner,job):
    start=time.monotonic();cpu=time.process_time_ns()
    result=runner.run_case((job['tag'],job['condition'],job['case']))
    require('error' not in result, 'Transport solver failed: '+str(result.get('error')))
    # Original run_case already converts its raw mb*D/nm² yields by1e-13/e.
    output=dict(job=job,T_K=float(result['T']),incident_D_flux_m2_s=float(result['F']),
        t_run_s=float(result['t_run']),end={k:float(v) for k,v in result['end'].items()},
        m30={k:float(v) for k,v in result['m30'].items()},
        physical_yield_per_true_C={k:float(v) for k,v in result['yield'].items()},
        process_cpu_ns=str(time.process_time_ns()-cpu),elapsed_wall_s=time.monotonic()-start,
        yield_scope='DDp reactions/true charged C, original fitted Ue probes held fixed; no detector or cup correction')
    lifetime=float(result['ol_life_s'])
    if job['case'].get('overlayer') is None:
        output.update(overlayer_status='none',overlayer_removal_h=None)
    elif lifetime==float('inf'):
        output.update(overlayer_status='not_removed_by_end',overlayer_removal_h=None,
                      overlayer_removal_lower_bound_h=float(result['t_run'])/3600)
    else:
        require(lifetime>=0,'Invalid overlayer lifetime')
        output.update(overlayer_status='removed_in_run',overlayer_removal_h=lifetime/3600)
    return output


def updated_table(records,original):
    reference={(r['condition'],float(r['E_cm_keV'])):r for r in original if r['knob']=='base'}
    computed={(r['job']['scenario']['id'],r['job']['condition'],r['job']['case']['E_cm']):r
              for r in records if r['job']['knob']=='base'}
    table=[]
    for r in records:
        j=r['job'];scenario=j['scenario']['id'];condition=j['condition'];energy=j['case']['E_cm']
        e,m=r['end'],r['m30'];old=reference[(condition,energy)];base=computed.get((scenario,condition,energy))
        row=dict(parameter_set=scenario,scientific_role=j['scenario']['kind'],condition=condition,E_cm_keV=energy,
            knob=j['knob'],T_K=r['T_K'],DM5_30min=m['DM5'],DM5_end=e['DM5'],DM20_end=e['DM20'],DM70_end=e['DM70'],
            nD5_end=e['nD5'],v5_end=e['v5'],v70_end=e['v70'],ads_ML_end=e['ads_ML'],
            overlayer_nm_30min=m['d_ol_nm'],overlayer_life_h=r['overlayer_removal_h'],
            overlayer_status=r['overlayer_status'],Y_UeP_per_C_end=r['physical_yield_per_true_C']['UeP'])
        denominator=dict(DM5_30min=base['m30']['DM5'],DM5_end=base['end']['DM5'],DM70_end=base['end']['DM70'],
                         v5_end=base['end']['v5'],v70_end=base['end']['v70'],Y_end=base['physical_yield_per_true_C']['UeP']) if base else {
                         k:float(old[k]) for k in ['DM5_30min','DM5_end','DM70_end','v5_end','v70_end']}
        if not base:denominator['Y_end']=float(old['Y_UeP_per_C_end'])
        numerators=dict(DM5_30min=m['DM5'],DM5_end=e['DM5'],DM70_end=e['DM70'],v5_end=e['v5'],v70_end=e['v70'],
                        Y_end=row['Y_UeP_per_C_end'])
        for key,value in numerators.items():
            row['ratio_'+key]=float(value/denominator[key]) if denominator[key]>0 else None
        table.append(row)
    return table


def main():
    p=argparse.ArgumentParser();p.add_argument('--reference-runner',type=Path,required=True)
    p.add_argument('--original-table',type=Path,required=True);p.add_argument('--parameter-sets',type=Path,required=True)
    p.add_argument('--pdrh',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cpu-budget-s',type=float,required=True);p.add_argument('--only-case')
    a=p.parse_args();start=time.process_time_ns();sys.path.insert(0,str(a.reference_runner.resolve().parent))
    spec=importlib.util.spec_from_file_location('knob_reference_runner',a.reference_runner)
    runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
    cache=Path(runner.T.BCA_FILE)
    require(cache.is_file() and manifest(cache)['sha256']=='4fd7ce875442bc6932646d08c4821e5764dba619fdccf3d32b69f22e3a231045',
            'Original BCA cache required; refuse automatic table rebuild')
    scenarios=json.loads(a.parameter_sets.read_bytes());require(scenarios['schema']=='coherence-transport-parameter-scenarios-v1','Wrong parameter-set role')
    rh=json.loads(a.pdrh.read_bytes());require(rh['schema']=='coherence-Pd80Rh20-conditional-isotherm-v1'
        and set(rh['material_overrides'])=={'a300','b300','dH'},'Wrong PdRh conditional bulk-parameter role')
    jobs=build_jobs(runner,scenarios['scenarios'],rh)
    if a.only_case is not None:jobs=[j for j in jobs if j['id']==a.only_case]
    require(jobs,'No requested job');a.output.mkdir()
    with a.original_table.open() as stream:original=list(csv.DictReader(stream))
    producer=dict(reference_runner=manifest(a.reference_runner),transport=manifest(Path(runner.T.__file__)),
        bca=manifest(Path(runner.T.bca.__file__)),cache=manifest(cache),parameter_sets=manifest(a.parameter_sets),
        PdRh=manifest(a.pdrh),original_table=manifest(a.original_table),
        sources=[manifest(Path(__file__).parent/name) for name in ['knob_cases.py','run_updated_knobs.py']],
        selection=('single-case consumer proof' if a.only_case else 'all changed Pd rows; unchanged Ti rows retained from source table'),
        reuse_scope='Original baseline estimates retained, not relabelled new computations. Ti/TiD physics unchanged by Pd Eb/Eads; original rows remain authoritative.',
        cpu_budget_s=a.cpu_budget_s)
    write_json(a.output/'producer.json',producer);records=[]
    for job in jobs:
        require((time.process_time_ns()-start)/1e9<a.cpu_budget_s,'Whole rerun CPU allocation exhausted')
        result=run_job(runner,job);write_json(a.output/(job['id']+'.json'),result);records.append(result)
        print('KNOB_CASE_COMPLETE '+job['id'],flush=True)
        require(result['elapsed_wall_s']<=90,'Threefold30s case-cost deviation; checkpoint and stop before more jobs')
    result=dict(schema='coherence-updated-knob-table-v1',original_estimate_rows=original,
        recomputed_rows=updated_table(records,original),case_outputs=records,
        source=manifest(a.output/'producer.json'),process_cpu_ns=str(time.process_time_ns()-start),
        scope='Conditional electronic-input sensitivities are not validated finite-T transport replacements. Barrier/prefactor never replaces sticking; mtrap remains6.',
        composition_warning=rh['composition_warning'])
    write_json(a.output/'updated-knobs.json',result)
    print('UPDATED_KNOBS_COMPLETE',flush=True)


if __name__=='__main__':main()
