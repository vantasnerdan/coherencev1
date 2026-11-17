"""Join completed compare tables, explicit reference scenarios and resolved runs.
No transport, fitting or primary numerical replay is performed here.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import csv,gzip,json,os,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from g4_beam.merit.api import ENERGIES,FOUR,CAUSES,evaluate_fit,write_csv
from g4_beam.causes.references import zero_sheet_fit,zero_sheet_reference,without_exposed_pd_sheet,pd_sheet_identity,ZERO_SHEET_SURFACE
from g4_beam.ald.designs import HERE,design_target,_write_profiles
from g4_beam.ald.dense import parameters
WORKERS=int(os.environ.get('G5_FINAL_WORKERS','4'))
if not 1<=WORKERS<=4:raise ValueError('G5_FINAL_WORKERS must be1–4')


def secondary(job):
    design,denominator=job;target=design_target(design);projected=without_exposed_pd_sheet(target)
    values=evaluate_fit(projected,ENERGIES,'surface',zero_sheet_fit(),reference_E_cm_keV=None)['physical_yield_per_C']
    return design['name'],{float(e):dict(**pd_sheet_identity(target,float(e)),
            refitted_zero_sheet_surface_merit=float(values[i]/denominator[i]),
            refitted_zero_sheet_surface_design_physical_per_true_C=float(values[i]),
            refitted_zero_sheet_surface_reference_physical_per_true_C=float(denominator[i]),
            refitted_zero_sheet_surface_fit_chi2=zero_sheet_fit()['chi2'],
            zero_sheet_reference_status='shape-only; fixed original bulk; not joined transport',
            zero_sheet_exposed_Pd_removed_on_both_sides=True,
            zero_sheet_retains_Pd_sheet_area_fraction=pd_sheet_identity(projected,float(e))['retains_beam_side_Pd_sheet_area_fraction']) for i,e in enumerate(ENERGIES)}


def TiD_calibration():
    """Same-energy counts calibration, with the SAME fitted reference on each side."""
    from g4_beam.merit.api import reference
    source=list(csv.DictReader((HERE.parent/'merit/transported_material_ratios.csv').open()))
    rows=[dict(r,reference_scenario='canonical_full_sheet') for r in source if r['numerator']=='TiD' and r['denominator']=='PdD']
    energies=np.array(sorted({float(r['E_cm_keV']) for r in rows}))
    observed={float(r['E_cm_keV']):float(r['measured_SI1_counts_ratio']) for r in rows}
    fit=zero_sheet_fit()
    numerator=evaluate_fit(without_exposed_pd_sheet(reference('TiD')),energies,'surface',fit,reference_E_cm_keV=None)['physical_yield_per_C']
    denominator=evaluate_fit(zero_sheet_reference(),energies,'surface',fit,reference_E_cm_keV=None)['physical_yield_per_C']
    for e,ratio in zip(energies,numerator/denominator):
        rows.append(dict(cause='surface',E_cm_keV=float(e),numerator='TiD',denominator='PdD',
                         model_same_E_ratio=float(ratio),measured_SI1_counts_ratio=observed[float(e)],
                         model_over_observed=float(ratio/observed[float(e)]),reference_scenario=ZERO_SHEET_SURFACE,
                         normalization='same-E shared collection/acceptance; explicit refitted zero-Pd-sheet scenario'))
    write_csv(HERE/'cap_TiD_calibration.csv',rows)
    return {(r['reference_scenario'],r['cause'],float(r['E_cm_keV'])):r for r in rows}


def calibrated(value,calibration):
    return float(value)/float(calibration['model_over_observed']) if calibration and value not in (None,'') else None


def cap_spread(designs,resolved):
    caps={d['name']:d for d in designs if d['family']==2};groups={}
    for result in resolved:
        if result['target'] not in caps:continue
        design=caps[result['target']];p=parameters(design)
        if p['sputter'] not in (.1,10.):continue
        label=('retained' if p['retention']==1 else 'lost')+('_slow' if p['sputter']==.1 else '_fast')
        key=(design['material'],design['thickness_nm'],result['reference_scenario'],result['cause'],float(result['E_cm_keV']))
        row=groups.setdefault(key,dict(material=key[0],thickness_nm=key[1],reference_scenario=key[2],cause=key[3],E_cm_keV=key[4]))
        row[label+'_target']=result['target']
        for field in ('integrated_merit','end_merit','peak_merit','time_merit_ge10_s','clearance_s'):
            row[label+'_'+field]=float(result[field])
        row[label+'_TiD_counts_calibrated_integrated_merit']=result['TiD_counts_calibrated_integrated_merit']
        row['TiD_counts_model_over_observed_factor']=result['TiD_counts_model_over_observed_factor']
    for row in groups.values():
        assert all(label+'_integrated_merit' in row for label in ('retained_slow','retained_fast','lost_slow','lost_fast'))
    write_csv(HERE/'cap_endpoint_spread.csv',list(groups.values()))


def laminate_inventory(designs):
    rows=[]
    from g4_beam.merit.api import reference
    for design in (d for d in designs if d['family']==6):
        for e in FOUR:
            profile=design['initial_profiles'][str(float(e))];top=profile['state']['film_total_m']
            film=[r for r in profile['layers'] if r[0]<top];periods=sum(r[4]!=46 for r in film)
            loaded=profile['state']['interface_count']>0
            extra=sum((r[1]-r[0])*(r[3]-.7*r[6]) for r in film if r[4]==46) if loaded else sum((r[1]-r[0])*r[3] for r in film)
            def inventory(layers,origin):
                return sum(max(0.,min(r[1]-origin,.1e-9)-max(r[0]-origin,0.))*r[3] for r in layers)
            substrate_extra=inventory(profile['layers'],top)-inventory(reference().profiles(float(e))[0],0.)
            expected=2*periods*1e19 if loaded else 0.
            actual=extra+substrate_extra
            error=abs(actual-expected)/max(expected,1e19)
            first=film[0];free_extra=(first[1]-first[0])*(first[3]-.7*first[6]) if loaded else (first[1]-first[0])*first[3]
            assert error<1e-10 and abs(free_extra)<1e9
            rows.append(dict(target=design['name'],E_cm_keV=float(e),oxide_Pd_contacts=2*periods,
                             loaded_inventory_endpoint=loaded,expected_interface_D_m2=expected,
                             actual_excess_interface_D_m2=actual,relative_inventory_error=error,
                             free_surface_extra_D_m2=free_extra,
                             status='frozen ideal inventory; not a measured capacity or kinetic prediction'))
    write_csv(HERE/'laminate_inventory.csv',rows)
    return max(r['relative_inventory_error'] for r in rows)


def main():
    start=time.monotonic();all_designs=[];rows=[]
    for stage in ('early','remaining'):
        all_designs+=json.loads(gzip.decompress((HERE/(stage+'_profiles.json.gz')).read_bytes()))
        rows+=list(csv.DictReader((HERE/(stage+'_compare.csv')).open()))
    all_designs+=json.loads(gzip.decompress((HERE/'nominal_design.json.gz').read_bytes()))
    rows+=list(csv.DictReader((HERE/'nominal_cap_compare.csv').open()))
    cachepath=HERE/'secondary_reference.json.gz'
    cache=json.loads(gzip.decompress(cachepath.read_bytes())) if cachepath.exists() else {}
    denominator=evaluate_fit(zero_sheet_reference(),ENERGIES,'surface',zero_sheet_fit(),reference_E_cm_keV=None)['physical_yield_per_C']
    calibration={(r['reference_scenario'],r['cause'],float(r['E_cm_keV'])):r
                 for r in csv.DictReader((HERE/'cap_TiD_calibration.csv').open())}
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures={pool.submit(secondary,(design,denominator)):design['name'] for design in all_designs if design['name'] not in cache}
        for future in as_completed(futures):
            name,values=future.result();cache[name]={str(k):v for k,v in values.items()}
            cachepath.write_bytes(gzip.compress(json.dumps(cache,separators=(',',':')).encode(),compresslevel=1))
            print('SECOND_REFERENCE',len(cache),name,flush=True)
    resolved=list(csv.DictReader((HERE/'resolved_integrals.csv').open()))
    curves=list(csv.DictReader((HERE/'resolved_yields.csv').open()))
    nominal=json.loads(gzip.decompress((HERE/'nominal_resolved.json.gz').read_bytes()))
    nominal_names={r['target'] for r in nominal.values()}
    resolved=[r for r in resolved if r['target'] not in nominal_names]+[row for r in nominal.values() for row in r['integrals']]
    curves=[r for r in curves if r['target'] not in nominal_names]+[row for r in nominal.values() for row in r['curves']]
    write_csv(HERE/'resolved_yields.csv',curves)
    groups={}
    for row in curves:
        key=(row['target'],row['reference_scenario'],row['cause'],float(row['E_cm_keV']))
        groups.setdefault(key,[]).append(row)
    serialization_error=0.;coarsening_changes=[];mean_differences=[]
    for result in resolved:
        key=(result['target'],result['reference_scenario'],result['cause'],float(result['E_cm_keV']))
        rr=sorted(groups[key],key=lambda row:float(row['time_s']))
        t=np.array([float(r['time_s']) for r in rr]);n=np.array([float(r['physical_design_per_true_C']) for r in rr]);d=np.array([float(r['physical_reference_per_true_C']) for r in rr])
        assert t[0]==0 and t[-1]==float(result['run_end_s']) and np.all(d>0)
        quotient=float(np.trapezoid(n,t)/np.trapezoid(d,t))
        serialization_error=max(serialization_error,abs(quotient/float(result['integrated_merit'])-1))
        indices=np.unique(np.r_[np.arange(0,len(t),2),len(t)-1])
        coarsened=float(np.trapezoid(n[indices],t[indices])/np.trapezoid(d[indices],t[indices]))
        change=abs(coarsened/quotient-1);coarsening_changes.append(change)
        mean=float(np.trapezoid(n/d,t)/t[-1]);mean_differences.append(abs(mean/quotient-1))
        result.update(alternate_snapshot_coarsened_merit=coarsened,coarsening_relative_change=change,
                      numerical_uncertainty_status='sampled coarsening diagnostic; not a rigorous enclosure',
                      coarsening_disagrees_on_10x=(coarsened>=10)!=(quotient>=10))
    assert serialization_error<2e-12
    assert len(rows)==len(all_designs)*len(ENERGIES)*len(CAUSES)
    assert len({(r['target'],r['cause'],r['E_cm_keV']) for r in rows})==len(rows)
    for result in resolved:
        cal=calibration.get((result['reference_scenario'],result['cause'],float(result['E_cm_keV'])))
        is_cap=int(result['family'])==2
        result['TiD_counts_model_over_observed_factor']=float(cal['model_over_observed']) if is_cap and cal else None
        result['TiD_counts_calibrated_integrated_merit']=calibrated(result['integrated_merit'],cal) if is_cap else None
    write_csv(HERE/'resolved_integrals.csv',resolved)
    cap_spread(all_designs,resolved)
    inventory_error=laminate_inventory(all_designs)
    (HERE/'verification.json').write_text(json.dumps(dict(
            complete_zero_to_runend=True,positive_same_age_reference=True,
            integrated_serialization_max_relative_error=serialization_error,
            max_alternate_snapshot_coarsening_relative_change=max(coarsening_changes),
            coarsening_10x_disagreement_cases=[{k:r[k] for k in ('target','reference_scenario','cause','E_cm_keV')} for r in resolved if r['coarsening_disagrees_on_10x']],
            max_difference_from_time_average_merit=max(mean_differences),
            canonical_rows=len(rows),designs=len(all_designs),all_requested_energy_cause_rows_present=True,
            laminate_interface_inventory_max_relative_error=inventory_error,
            unselected_integral_estimates_published=False),indent=2))
    index={(r['target'],r['reference_scenario'],r['cause'],float(r['E_cm_keV'])):r for r in resolved}
    selected={r['target'] for r in resolved}
    design_index={d['name']:d for d in all_designs}
    for row in rows:
        row.setdefault('self_R_defined',row.get('self_R') not in ('',None))
        row.setdefault('reported_R_defined',row.get('reported_R') not in ('',None))
        e=float(row['E_cm_keV']);name=row['target'];row.update(cache[name][str(e)])
        row['lifetime_clock']='provided-IV-as-incident retained Model clock; source proxy uncertainty'
        if int(row['family'])==2:
            p=parameters(design_index[name])
            row['scenario_retained_D_fraction']=p['retention'];row['scenario_sputter_factor']=p['sputter']
            row['scenario_sputter_lifetime_h']=float(row['lifetime_nominal_h'])/p['sputter']
        if int(row['family']) in (1,3,4) or row['material']=='BZY20_front':
            d0=float(row['thickness_nm'])
            row['model_soft_clearance_h']=float(row['lifetime_nominal_h'])*((d0-.001)+.02*np.log(d0/.001))/d0
            row['model_clearance_status']='analytic erosion-only1pm cutoff; retained provided-current clock and source sputter proxy'
        for prefix,scenario,cause in [('canonical','canonical_full_sheet',row['cause']),('zero_sheet',ZERO_SHEET_SURFACE,'surface')]:
            result=index.get((name,scenario,cause,e))
            row[prefix+'_integral_status']='resolved' if result else ('not_selected_for_resolved_rerun' if e in FOUR else 'four-energy_integral_not_requested')
            for field in ('integrated_merit','end_merit','peak_merit','time_merit_ge10_s'):
                row[prefix+'_'+field]=result[field] if result else None
        if int(row['family'])==2:
            cal=calibration.get(('canonical_full_sheet',row['cause'],e))
            row['TiD_counts_model_over_observed_factor']=float(cal['model_over_observed']) if cal else None
            row['TiD_counts_measured_same_E_ratio']=float(cal['measured_SI1_counts_ratio']) if cal else None
            row['TiD_counts_calibrated_merit_PdD']=calibrated(row['merit_PdD'],cal)
            row['TiD_counts_calibrated_integrated_merit']=calibrated(row['canonical_integrated_merit'],cal)
            zero_cal=calibration.get((ZERO_SHEET_SURFACE,'surface',e))
            row['zero_sheet_TiD_counts_calibrated_merit_PdD']=calibrated(row['refitted_zero_sheet_surface_merit'],zero_cal)
            row['zero_sheet_TiD_counts_calibrated_integrated_merit']=calibrated(row['zero_sheet_integrated_merit'],zero_cal)
            row['TiD_counts_calibration_status']='SI1 TiD/PdD counts proxy; non-Ti caps require transfer assumption' if cal else 'no same-energy SI1 anchor'
    write_csv(HERE/'all_designs.csv',rows)
    _write_profiles(HERE/'end_run_depth_profiles.csv',all_designs)
    best=[]
    for family in range(1,7):
        for cause in CAUSES:
            for e in FOUR:
                candidates=[r for r in rows if int(r['family'])==family and r['cause']==cause and float(r['E_cm_keV'])==float(e)]
                chosen=max(candidates,key=lambda r:float(r['merit_PdD']))
                best.append(dict(family=family,cause=cause,E_cm_keV=float(e),best_endpoint_target=chosen['target'],
                                 endpoint_merit=float(chosen['merit_PdD']),zero_sheet_surface_endpoint_merit=float(chosen['refitted_zero_sheet_surface_merit']),
                                 canonical_integral_status=chosen['canonical_integral_status'],canonical_integrated_merit=chosen['canonical_integrated_merit'],
                                 feasibility=chosen['feasibility'],surface_status=chosen['surface_status']))
    write_csv(HERE/'best_per_family.csv',best)
    headline=[]
    for scenario,cause in [('canonical_full_sheet',c) for c in CAUSES]+[(ZERO_SHEET_SURFACE,'surface')]:
        for e in FOUR:
            candidates=[r for r in resolved if r['reference_scenario']==scenario and r['cause']==cause and float(r['E_cm_keV'])==float(e)]
            chosen=max(candidates,key=lambda r:float(r['integrated_merit']))
            headline.append(chosen)
    write_csv(HERE/'best_resolved_integrals.csv',headline)
    summary=dict(designs=len(all_designs),canonical_compare_rows=len(rows),resolved_designs=len(selected),
                 reference_scenarios=['canonical_full_sheet',ZERO_SHEET_SURFACE],wall_s=time.monotonic()-start,
                 nonrerun_integral_estimates_published=False,source_unknowns='retention; sputter proxies; material-specific surface kinetics; continuity; true-charge calibration')
    (HERE/'delivery_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
