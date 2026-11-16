"""Exercise actual supplied baseline profiles, not a transport replay.
Run from repo: python g4_beam/forward/consumer_example.py
The fitted Ue/TF inputs are supplied models, not first-principles calculations.
"""
from pathlib import Path
import csv
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from g4_beam.forward.consumer import forward_yields, WKB_NUCLEAR_RADIUS_M
from g4_beam.forward import model as M
from g4_beam.forward.thermal import QE
from g4_beam.forward.profile_thermal import transport_stopping
HERE=Path(__file__).resolve().parent

def tf_potential(lambda_A, first_m=WKB_NUCLEAR_RADIUS_M, last_m=1e-9):
    """Sample the supplied TF model in SI, including an exact inner endpoint."""
    r_m=np.geomspace(first_m,last_m,500)
    r_A=r_m/1e-10
    return r_m,14.3996/r_A*np.exp(-r_A/lambda_A)





def main():
    profiles=json.load(open(HERE/'inputs/profiles_baseline.json'))
    screen=json.load(open(HERE/'inputs/screening_placeholder.json'))
    records=profiles['records']
    energies=np.array([.25,.5,2.5,6.5])
    outputs=[]
    def append(mat,scenario,result,assumptions):
        for i,e in enumerate(result['E_cm_keV']):
            source=next(r['case'] for r in records if r['case']['condition']==mat and r['case']['E_cm_keV']==e)
            outputs.append(dict(material=mat,scenario=scenario,E_cm_keV=float(e),R=float(result['R'][i]),physical_reactions_per_true_C=float(result['physical_yield_per_C'][i]),detected_counts_per_cup_C=float(result['counts_per_C'][i]),ion_reactions_per_true_C=float(result['ion_yield_per_C'][i]),neutral_reactions_per_true_C=float(result['neutral_yield_per_C'][i]),beam_time_s=source['beam_time_s'],T_K=source['T_K'],input_interpretation=assumptions))
    def get_profiles(mat):
        return {float(r['case']['E_cm_keV']):(r['layers'],r['nads_m2'],r['state']) for r in records if r['case']['condition']==mat}
    for mat in M.MATS:
        u=screen['Ue_bulk_fit_keV']['M4_one_power_current'][mat][0]
        result=forward_yields(energies,get_profiles(mat),Ue_keV=u,stopping_fn=transport_stopping)
        append(mat,'bulk supplied M4 Ue',result,'actual original end-of-run profiles and transport/BCA stopping; fitted Ue; efficiency=collection=1, proton channel')
        print(mat,'physical reactions/C',result['physical_yield_per_C'],'R',result['R'],flush=True)
    pp=get_profiles('Pd');u=screen['Ue_bulk_fit_keV']['M4_one_power_current']['Pd'][0]
    meas=dict(efficiency=.0032,collected_fraction=lambda e:M.collected_fraction(e,1,1),current_A=lambda e:M.I_cup(e)*1e-6)
    thermal_fit=min((r for r in json.load(open(HERE/'fits.json')) if r['cand']=='thermal'),key=lambda r:r['chi2'])
    amplitude=thermal_fit['params']['A']
    for cand,p in [('layer',dict(Ue_s=8.1,d_nm=.0035)),('sigmoid_tt',dict(Ue_s=1.9,E0=1.36)),('fast_neutrals',dict(f_n=.28,E_n_lab_keV=5.,delta=1)),('thermal',dict(Ue_hot=.4,E_A_eV=.150,A=amplitude,T_track_K=1080.,E_R_eV=1.))]:
        mm=dict(meas,neutral_reference_current_A=138e-6)
        result=forward_yields(energies,pp,Ue_keV=u,candidate=cand,parameters=p,measurement=mm,stopping_fn=transport_stopping)
        append('Pd',cand,result,'actual saved Pd profiles; representative supplied/declared physical parameters; efficiency and cup collection explicit')
    # Sample the supplied TF analytic MODEL at both stated Pd lambda endpoints.
    for lam in screen['static_thomas_fermi']['lambda_A_range']['Pd']:
        result=forward_yields(energies,pp,potential=tf_potential(lam),stopping_fn=transport_stopping)
        append('Pd','supplied TF lambda_A='+str(lam),result,'TF estimate, NOT DFT; sampled total V_eV/r_m; WKB action-only enhancement of bare BH')
    # Meaningful unit/consumer invariant: uniform no-adsorbate limit against
    # the independently established target integral, including molecular flux.
    t=M.Target('Pd');nd=t.x*t.n*1e6;nh=t.n*1e6
    uniform=lambda e:([(0.,400e-9,nh,nd,46.,106.42,nh)],0.,{})
    ee=np.array([.5,2.5,6.5])
    uniform_result=forward_yields(ee,uniform,Ue_keV=.7)
    expected=M.beam_average(lambda e:t.Y_const(e,.7),2*ee)/QE
    error=np.max(np.abs(uniform_result['physical_yield_per_C']/expected-1))
    assert error<.003,('uniform integral discrepancy',error)
    # Surface-sheet units: nads m^-2 times sigma mb*1e-31, divided by e.
    sheet=lambda e:([(0.,400e-9,nh,0.,46.,106.42,nh)],1e19,{})
    sheet_result=forward_yields(ee,sheet,Ue_keV=.7)
    expected_sheet=M.beam_average(lambda e:1e19*M.reaction_sigma(e/2,.7)*1e-31,2*ee)/QE
    sheet_error=np.max(np.abs(sheet_result['physical_yield_per_C']/expected_sheet-1))
    assert sheet_error<1e-12,('surface-sheet unit error',sheet_error)
    # Measurement changes do not contaminate the physical output.
    shifted=forward_yields(ee,uniform,Ue_keV=.7,measurement=dict(efficiency=.02,collected_fraction=.5,current_A=lambda e:1e-4*e,extra_current_power=1))
    assert np.allclose(shifted['physical_yield_per_C'],uniform_result['physical_yield_per_C'],rtol=1e-12)
    assert np.allclose(shifted['counts_per_C'],.04*shifted['physical_yield_per_C'])
    empty=lambda e:([(0.,400e-9,nh,0.,46.,106.42,nh)],0.,{})
    try:forward_yields(ee,empty,Ue_keV=.7)
    except ValueError as exc:
        assert 'Reference signal is zero' in str(exc)
    else:raise AssertionError('Undefined zero-reference ratio was accepted')
    with open(HERE/'consumer_outputs.csv','w') as f:
        writer=csv.DictWriter(f,fieldnames=list(outputs[0]));writer.writeheader();writer.writerows(outputs)
    summary=dict(profile_source_commit='0aab89af99e2bad8ace81c4b413c033cf87812e2',profile_source_path='g4_beam/transport/baseline/profiles_full.json',profile_source_sha256='f06a1df5fcc98e204e621080e038569a171a19da7fe0cf900de8b3f5291aa51d',uniform_limit_max_relative_error=float(error),surface_sheet_max_relative_error=float(sheet_error),measurement_physical_separation=True,zero_reference_rejected=True,energies_cm_keV=energies.tolist(),screening_input_status='supplied fitted/TF model estimates, not new first-principles results')
    json.dump(summary,open(HERE/'consumer_verification.json','w'),indent=2)
    print(summary,flush=True)

def potential_layer_example():
    """Combined supplied V(r) and a surface-screening shift; no other replay."""
    records=json.load(open(HERE/'inputs/profiles_baseline.json'))['records']
    profile={float(r['case']['E_cm_keV']):(r['layers'],r['nads_m2'],r['state']) for r in records if r['case']['condition']=='Pd'}
    potential=tf_potential(.57)
    result=forward_yields([.5,6.5],profile,potential=potential,candidate='layer',
                          parameters=dict(d_nm=.0035,delta_Ue_surface_keV=3.),
                          stopping_fn=transport_stopping)
    assert np.all(np.isfinite(result['R'])) and np.all(result['physical_yield_per_C']>0)
    print(json.dumps(dict(scenario='supplied TF potential plus3keV surface shift',
                          R=result['R'].tolist(),physical_reactions_per_C=result['physical_yield_per_C'].tolist())))


def potential_domain_example():
    """Reject reporting-only grids, then recover on the actual saved-profile API."""
    records=json.load(open(HERE/'inputs/profiles_baseline.json'))['records']
    profile={float(r['case']['E_cm_keV']):(r['layers'],r['nads_m2'],r['state']) for r in records if r['case']['condition']=='Pd'}
    rejected={}
    for label,potential,expected in [
        ('reporting_window',tf_potential(.57,first_m=.01e-10,last_m=3e-10),'fixed 4 fm WKB nuclear radius'),
        ('missing_outer_turning_point',tf_potential(.57,last_m=.01e-10),'outer turning point')]:
        try:
            forward_yields([.5,6.5],profile,potential=potential,stopping_fn=transport_stopping)
        except ValueError as exc:
            assert expected in str(exc),str(exc)
            rejected[label]=str(exc)
        else:
            raise AssertionError('Incomplete potential domain accepted: '+label)
    full=forward_yields([.5,6.5],profile,potential=tf_potential(.57),stopping_fn=transport_stopping)
    extended=forward_yields([.5,6.5],profile,potential=tf_potential(.57,first_m=1e-15),stopping_fn=transport_stopping)
    assert np.all(np.isfinite(full['R'])) and np.all(full['physical_yield_per_C']>0)
    difference=float(np.max(np.abs(extended['physical_yield_per_C']/full['physical_yield_per_C']-1)))
    assert difference<.003,('First sample moved the physical WKB boundary',difference)
    domain=dict(fixed_WKB_nuclear_radius_m=WKB_NUCLEAR_RADIUS_M,resonance_channel_radius_is_separate_fm=7.,
                rejected_domains=rejected,full_grid_recovered=True,
                grid_below_radius_max_relative_difference=difference,
                full_grid_R=full['R'].tolist(),full_grid_physical_reactions_per_C=full['physical_yield_per_C'].tolist())
    summary=json.load(open(HERE/'consumer_verification.json'))
    summary['potential_domain']=domain
    json.dump(summary,open(HERE/'consumer_verification.json','w'),indent=2)
    print(json.dumps(domain),flush=True)



if __name__=='__main__':
    if '--potential-domain' in sys.argv:potential_domain_example()
    elif '--potential-layer' in sys.argv:potential_layer_example()
    else:main()
