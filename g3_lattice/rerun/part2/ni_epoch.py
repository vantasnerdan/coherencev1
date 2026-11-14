"""One893 Ni quantum/static/pump refit and the authorized110592-site phase case.
Quartic epoch is explicit and frozen once; no automatic late-input rerun.
"""
import argparse
import copy
import json
from math import pi,sqrt
from pathlib import Path
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'pd_x1_quantum'/'engine'))
from anharmonic import spectrum
from refit import consume,static_branch
from pd_x1_quantum.pd_x1_refit import time_and_bath
from lattice import HC,HBAR,run_case
from coherence import PhaseWindows,characterize,linear_guide


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--quartic-input',type=Path,required=True)
    ap.add_argument('--out',type=Path,default=Path('updated/ni_epoch'))
    args=ap.parse_args();folder=args.out;folder.mkdir(parents=True,exist_ok=True)
    quartic=json.loads(args.quartic_input.read_text())
    raw=json.loads((ROOT/'new_inputs/ni_x1_893/source-readout.json').read_text())
    coupling=raw['coupling']
    historical=json.loads((ROOT/'pd_x1_quantum/historical-inputs.json').read_text())
    old_material=next(row for row in historical['materials'] if row['name']=='NiH')
    before=next(row for row in json.loads((ROOT/'pd_x1_quantum/historical-anharmonic.json').read_text())['rows'] if row['material']=='NiH')
    material=dict(old_material,K2_eV_A2=coupling['on_site_k2_eV_A2'],
        K_Gamma_eV_A2=coupling['H_only_q0_stiffness_eV_A2'],J_sum_eV_A2=coupling['sum_neighbour_pulling_J_xx_eV_A2'],
        K4_100_eV_A4=quartic['K4_100_eV_A4'],K4_111_eV_A4=quartic['K4_111_eV_A4'],
        source_entry='893 plus explicitly frozen quartic epoch',reference_INS_meV=89.)
    A=material['K4_100_eV_A4'];D=material['K4_111_eV_A4'];B=(3*D-A)/2
    if min(A,D)<=0:raise ValueError('Actual supplied quartic is not bounded; no forced positive saturation or hidden legacy substitution')
    actual=spectrum(material['K2_eV_A2'],A,B,material['mass_eVc2'],10,material['J_sum_eV_A2'])
    qrow=dict(material='NiH',harmonic_quantum_meV=1000*sqrt(HC*HC*material['K2_eV_A2']/material['mass_eVc2']),spectrum=actual)
    gain=consume(material,{'rows':[qrow]})
    for name in ('N1_proxy_zero_to_one_meV','conditional_critical_loading','conditional_critical_density_m3','threshold_scope'):
        gain.pop(name)
    gamma_hf=material['K2_eV_A2']-material['J_sum_eV_A2']+(3*A+2*B)*actual['local_mean_square_axis_A2_T0']
    gain['Gamma_Hartree_zero_temperature_curvature_eV_A2']=gamma_hf
    gain['Gamma_Hartree_zero_temperature_quantum_meV']=1000*sqrt(HC*HC/material['mass_eVc2']*gamma_hf)
    gain['onresonance_bare_lattice_h_single_RPA_pole_projection']=gain['effective_curvature_modulation_threshold']*gain['effective_pole_spring_eV_A2']/gamma_hf
    primitive=next(row for row in json.loads((ROOT/'pd_x1_quantum/primitive-strains.json').read_text())['model_inputs'] if row['material']=='NiH')
    matched=copy.deepcopy(primitive);matched['K_zero_eV_A2']=material['K2_eV_A2']-material['J_sum_eV_A2']
    old_matched=copy.deepcopy(primitive);old_matched['K_zero_eV_A2']=old_material['K2_eV_A2']-old_material['J_sum_eV_A2']
    quantum=dict(material_inputs=material,quartic_epoch=quartic,spectrum=actual,before_spectrum=before['spectrum'],
        pumped=gain,units_bath=time_and_bath(material,actual),
        static_matched=static_branch(material,actual,matched),
        static_before_matched=static_branch(old_material,before['spectrum'],old_matched),
        static_scope='Matched H-only Kii-sumJ + primitive Ni hydro CONTRAST; frozen local per-axis variance/quartics, extra elastic deltaSumJ0; beyond+/-1% is extrapolation, not SSCHA',
        INS_scope='89meV is a band peak, not an identified Gamma-mode validation',
        attribution='All harmonic/profile and quartic-epoch changes, NOT isolated coupling causality')
    (folder/'quantum.json').write_text(json.dumps(quantum,indent=2,allow_nan=False)+'\n')
    parameters=json.loads((ROOT/'updated/parameters.json').read_text())
    p=parameters['materials']['NiH'];old_p=copy.deepcopy(p)
    p.update(Kii_eV_A2=material['K2_eV_A2'],K4_100_eV_A4=A,K4_111_eV_A4=D,
        zero_point_axis_variance_A2=actual['local_mean_square_axis_A2_T0'],
        collective_strain_slope_eV_A2=primitive['strain']['hydro']['dK_depsilon_eV_A2'])
    for shell in p['shells']:
        matches=[item for item in coupling['shells'] if abs(item['distance_A']-shell['distance_A'])<1e-6]
        if len(matches)!=1:raise ValueError('New Ni spatial shell did not bind uniquely')
        shell['J_sum_eV_A2']=matches[0]['sum_pulling_J_xx_eV_A2']
    if abs(sum(item['J_sum_eV_A2'] for item in p['shells'])-material['J_sum_eV_A2'])>1e-12:
        raise ValueError('New ALL-shell Ni coupling was not preserved')
    p['status']='893 adopted fixed-k3 harmonic correction; quartic epoch '+quartic['source_label']+'; mixed='+str(quartic['mixed_profile'])+'; scalar minimum-image shell continuation'
    (folder/'parameters.json').write_text(json.dumps(parameters,indent=2,allow_nan=False)+'\n')
    source_cases=json.loads((ROOT/'coherence_cases.json').read_text())['cases']
    case=copy.deepcopy(next(row for row in source_cases if row['label']=='cold_x1_h008_L48_t40'))
    old_unit=HBAR/sqrt(HC*HC*old_p['Kii_eV_A2']/old_p['mass_eVc2'])*1e12
    target_ps=case['steps']*case['dt']*old_unit
    new_unit=quantum['units_bath']['time_unit_ps']
    case.update(label='NiH_893_L48_matched_time',material='NiH',arrangement='random',
        steps=round(target_ps/(case['dt']*new_unit)))
    guide=linear_guide(p,(48,48,48),case['pump_h'],case['pump_ratio'],case['dt'])
    stride=case.get('stride',2);samples=round(16/(guide['pump_half_frequency_THz']*new_unit*case['dt']*stride))
    window=PhaseWindows(guide['pump_half_frequency_THz'],samples,p['a_A'])
    legacy_sum=np.zeros((48,48,48),complex);legacy_count=0;first_time=None;last_time=None
    def observe(q,occupied,t):
        nonlocal legacy_count,first_time,last_time
        window(q,occupied,t)
        legacy_sum[:]+=q[0]*np.exp(-2j*pi*guide['pump_half_frequency_THz']*t)
        legacy_count+=1;last_time=t
        if first_time is None:first_time=t
    started=time.perf_counter();summary,fields=run_case(parameters,case,observe=observe)
    joint,joint_arrays=characterize(fields['driven_amplitude_A'],fields['occupied'],p['a_A'])
    legacy=np.zeros_like(fields['driven_amplitude_A']);legacy[0]=2*legacy_sum/legacy_count
    matched_correlation,matched_arrays=characterize(legacy,fields['occupied'],p['a_A'])
    output=dict(summary=summary,matched_legacy_filter=matched_correlation,joint_harmonic_filter=joint,
        equal_time_16cycle_envelopes=window.finish(),linear_guide=guide,
        comparison_protocol=dict(old_case='cold_x1_h008_L48_t40',old_source='731/757 frozen epoch',
            target_elapsed_ps=target_ps,actual_elapsed_ps=case['steps']*case['dt']*new_unit,
            observed_interval_ps=[first_time,last_time],legacy_samples=legacy_count,
            old_whole_record_phase_integral_nm=2.005,old_frozen_Gamma_RPA_meV=90.772,
            clocks='Same h.08,pump_ratio1,dt.3,seed13,T0; carrier retargeted to UPDATED classical Gamma Hartree reference, NOT a fixed laboratoryTHz',
            meaning='One changed-input110592-site trajectory only; no million-site or size/time campaign. Compare all-input profiles, not isolated coupling.'),
        quartic_epoch=quartic,phase_elapsed_s=time.perf_counter()-started)
    (folder/'phase.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(folder/'phase_arrays.npz',**matched_arrays,joint_amplitude_x_A=fields['driven_amplitude_A'][0])
    print(json.dumps({'local01_meV':actual['zero_to_one_meV'],'Gamma_RPA_meV':gain['mode_quantum_meV'],
        'static_strain0_300':[row['critical_linear_strain'] for row in quantum['static_matched']],
        'phase_integral_nm':matched_correlation['phase_only']['integral_correlation'],
        'phase_elapsed_s':output['phase_elapsed_s'],'quartic_epoch':quartic}),flush=True)


if __name__=='__main__':main()
