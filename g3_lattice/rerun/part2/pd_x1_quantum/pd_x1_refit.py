"""Adopted762 local Pd well + same-column Gamma RPA; no lattice dynamics/DFT.
Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python pd_x1_refit.py > result.json
Unchanged original local-well/RPA implementation is imported from engine/.
"""
import hashlib
import json
from math import copysign, pi, sqrt
from pathlib import Path
import resource
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'engine'))
from anharmonic import spectrum
from refit import consume
from domain import HC, HBAR, KB
from strain_tensor_law import calibration_mask, feature_tensor, stiffness_tensor


def load(name):
    return json.loads((ROOT/name).read_text())


def positive_roots(c0,c1,c2):
    """Stable real quadratic roots; no epsilon floor or forced crossing."""
    if c2==0:
        return [-c0/c1] if c1 and -c0/c1>0 else []
    discriminant=c1*c1-4*c2*c0
    if discriminant<0:return []
    q=-.5*(c1+copysign(sqrt(discriminant),c1))
    candidates=[q/c2,c0/q] if q else [-c1/(2*c2)]
    return sorted(set(float(x) for x in candidates if x>0))


def hydro_polynomial(model):
    # Extract formal polynomial coefficients, not physical predictions at unit strain.
    c=np.asarray(model['coefficients_eV_A2'])
    plus=feature_tensor(np.eye(3))[:len(c),0,0]
    minus=feature_tensor(-np.eye(3))[:len(c),0,0]
    return float(c@((plus-minus)/2)),float(c@((plus+minus)/2))


def static_matched(material,spec,linear,quadratic,model,law_name):
    A=material['K4_100_eV_A4'];D=material['K4_111_eV_A4'];B=(3*D-A)/2
    gamma0=material['K2_eV_A2']-material['J_sum_eV_A2']
    output=[]
    for temperature,key in [(0,'local_mean_square_axis_A2_T0'),(300,'local_mean_square_axis_A2_T300')]:
        variance=spec[key]
        correction=(3*A+2*B)*variance  # per-axis cubic Hessian, NOT total 3D MSD
        c0=gamma0+correction
        roots=positive_roots(c0,linear,quadratic)
        critical=roots[0] if roots else None
        def r(epsilon):
            if law_name=='PartA640 quadratic contrast':
                contrast=stiffness_tensor(np.eye(3)*epsilon,model,extrapolate=True)[0,0]-model['K0_eV_A2']
            else:contrast=linear*epsilon
            return float(c0+contrast)
        samples=[f*critical for f in [0.,.5,.9,.99]] if critical is not None else [0.,.01,.025,.05]
        soft=[]
        for epsilon in samples:
            curvature=r(epsilon);energy2=HC*HC/material['mass_eVc2']*curvature
            soft.append({'linear_strain':epsilon,'curvature_eV_A2':curvature,
                'omega_squared_ps_inverse2':energy2/(HBAR*1e12)**2,
                'undamped_frequency_THz':sqrt(energy2)/(2*pi*HBAR*1e12) if energy2>=0 else None,
                'outside_measured_hydro_range':epsilon>.01,
                'outside_PartA_coordinate_box':not bool(calibration_mask(np.eye(3)*epsilon,model))})
        beyond=None
        if critical is not None and critical<.05:
            epsilon=critical+.001;curvature=r(epsilon)
            if curvature<0:
                amplitude=sqrt(-curvature/D)
                beyond={'linear_strain':epsilon,'direction':'111','displacement_magnitude_A':amplitude,
                    'Cartesian_components_A':[amplitude/sqrt(3)]*3,
                    'energy_gain_per_isotope_eV':-curvature*curvature/(4*D),
                    'outside_measured_hydro_range':epsilon>.01}
        vertex=-linear/(2*quadratic) if quadratic>0 else None
        output.append({'T_K':temperature,'local_axis_variance_A2':variance,
            'matched_bare_Gamma_eV_A2':gamma0,'Hartree_correction_eV_A2':correction,
            'zero_strain_renormalized_curvature_eV_A2':c0,
            'law':law_name,'critical_positive_linear_strain':critical,
            'all_positive_roots_of_fitted_polynomial':roots,
            'crosses_below_5percent':critical is not None and critical<.05,
            'outside_measured_hydro_range_at_onset':critical>.01 if critical is not None else None,
            'polynomial_minimum_positive_strain':vertex if vertex is not None and vertex>0 else None,
            'polynomial_minimum_curvature_eV_A2':r(vertex) if vertex is not None and vertex>0 else None,
            'soft_mode_approach':soft,'just_beyond_within_5percent':beyond})
    return output


def time_and_bath(material,spec):
    mass=material['mass_eVc2'];k=material['K2_eV_A2']
    E0=sqrt(HC*HC*k/mass);unit_ps=HBAR/E0*1e12
    length_A=HC/sqrt(mass*E0)
    tau=material['reference_amplitude_decay_ps'];kappa=2/tau
    pole_mass=1/spec['RPA_Gamma_oscillator_strength']
    return {'bare_harmonic_energy_unit_meV':E0*1000,'time_unit_ps':unit_ps,'length_unit_A':length_A,
        'theta_at300K':KB*300/E0,'physical_kappa_friction_ps_inverse':kappa,
        'physical_amplitude_decay_ps':tau,'gamma_atomic_velocity_code':kappa*unit_ps,
        'Gamma_pole_mass_over_bare':pole_mass,
        'gamma_Gamma_pole_bare_force_coefficient_for_same_kappa':pole_mass*kappa*unit_ps,
        'kappa_over_RPA_angular_frequency':kappa*HBAR*1e12/spec['RPA_Gamma_quantum_eV'],
        'if_literal_atomic_gamma_0p01':{'kappa_ps_inverse':.01/unit_ps,'amplitude_decay_ps':2*unit_ps/.01},
        'if_literal_Gamma_bare_force_gamma_0p01':{'kappa_ps_inverse':.01/(pole_mass*unit_ps),'amplitude_decay_ps':2*pole_mass*unit_ps/.01},
        'physical_dt_ps_at_dimensionless_dt0p3':.3*unit_ps}


def main():
    started=time.perf_counter()
    adopted=load('input.json');source=adopted['summary'];column=adopted['source_column']
    historical_inputs=load('historical-inputs.json')
    historical_spectra=load('historical-anharmonic.json');historical=load('historical-refit.json')
    primitive=next(r for r in load('primitive-strains.json')['model_inputs'] if r['material']=='PdD')
    model=load('partA-law.json')
    assert source['on_site_Kii_force_eV_A2']==column['on_site_k2_eV_A2']
    assert source['sum_J_xx_eV_A2']==column['sum_neighbour_pulling_J_xx_eV_A2']
    assert model['reference_a_A']==4.084
    k=column['on_site_k2_eV_A2'];J=column['sum_neighbour_pulling_J_xx_eV_A2']
    A=source['quartics']['100']['K4_joint_eV_A4'];D=source['quartics']['111']['K4_joint_eV_A4']
    if min(A,D)<=0:raise ValueError('Cubic quartic is not bounded; an actual higher term is required')
    B=(3*D-A)/2;linear,quadratic=hydro_polynomial(model)
    rows=[]
    for name in ['PdH','PdD']:
        old=next(r for r in historical_inputs['materials'] if r['name']==name)
        oldq=next(r for r in historical_spectra['rows'] if r['material']==name)
        oldg=next(r for r in historical['rows'] if r['material']==name)
        material=dict(old,K2_eV_A2=k,K_Gamma_eV_A2=column['H_only_q0_stiffness_eV_A2'],
            J_sum_eV_A2=J,K4_100_eV_A4=A,K4_111_eV_A4=D,source_entry='entry_762 adopted fixed-k3 profile')
        small=spectrum(k,A,B,material['mass_eVc2'],8,J)
        actual=spectrum(k,A,B,material['mass_eVc2'],10,J)
        # This is the ONE changed-Hamiltonian truncation comparison; no native gate replay.
        parallel_small={r['block_degeneracy']:r for r in small['parallel_one_to_two_lines']}
        changes={'local01_meV':actual['zero_to_one_meV']-small['zero_to_one_meV'],
            'RPA_Gamma_meV':1000*(actual['RPA_Gamma_quantum_eV']-small['RPA_Gamma_quantum_eV']),
            'parallel12_meV':{str(r['block_degeneracy']):r['one_to_two_meV']-parallel_small[r['block_degeneracy']]['one_to_two_meV'] for r in actual['parallel_one_to_two_lines']},
            'perpendicular12_meV':actual['perpendicular_one_to_two_lines'][0]['one_to_two_meV']-small['perpendicular_one_to_two_lines'][0]['one_to_two_meV'],
            'zero_point_axis_variance_A2':actual['local_mean_square_axis_A2_T0']-small['local_mean_square_axis_A2_T0'],
            'thermal_axis_variance300_A2':actual['local_mean_square_axis_A2_T300']-small['local_mean_square_axis_A2_T300']}
        maximum=max(abs(changes['local01_meV']),abs(changes['RPA_Gamma_meV']),abs(changes['perpendicular12_meV']),*(abs(x) for x in changes['parallel12_meV'].values()))
        changes['max_energy_change_meV']=float(maximum);changes['resolves_0p01meV_by_basis_comparison']=bool(maximum<.01)
        qrow={'material':name,'harmonic_quantum_meV':1000*sqrt(HC*HC*k/material['mass_eVc2']),'spectrum':actual}
        gain=consume(material,{'rows':[qrow]})
        # Legacy N=1 proxy and toy occupancy threshold are not current model outputs.
        gain={key:value for key,value in gain.items() if key not in ['N1_proxy_zero_to_one_meV','conditional_critical_loading','conditional_critical_density_m3','threshold_scope']}
        gamma_hf=k-J+(3*A+2*B)*actual['local_mean_square_axis_A2_T0']
        gain['Gamma_Hartree_zero_temperature_curvature_eV_A2']=gamma_hf
        gain['Gamma_Hartree_zero_temperature_quantum_meV']=1000*sqrt(HC*HC/material['mass_eVc2']*gamma_hf)
        gain['onresonance_bare_lattice_h_single_RPA_pole_projection']=gain['effective_curvature_modulation_threshold']*gain['effective_pole_spring_eV_A2']/gamma_hf
        gain['bare_h_definition']='lattice deltaK=h_code*Gamma_Hartree0; weak-drive single-RPA-pole on resonance, NOT exact nonlinear lattice onset'
        gain['RPA_clock_over_classical_Hartree_clock']=gain['mode_quantum_meV']/gain['Gamma_Hartree_zero_temperature_quantum_meV']
        rows.append({'material':name,'material_inputs':{key:material[key] for key in ['mass_eVc2','a_A','K2_eV_A2','K_Gamma_eV_A2','J_sum_eV_A2','K4_100_eV_A4','K4_111_eV_A4','reference_amplitude_decay_ps']},
            'before_spectrum':oldq['spectrum'],'spectrum':actual,'basis8_to10':changes,
            'before_pumped':oldg,'pumped':gain,
            'before_units_bath':time_and_bath(old,oldq['spectrum']),'units_bath':time_and_bath(material,actual),
            'static_478_published':oldg['static_Hartree'],
            'static_before_matched_PartA':static_matched(old,oldq['spectrum'],linear,quadratic,model,'PartA640 quadratic contrast'),
            'static_after_matched_PartA':static_matched(material,actual,linear,quadratic,model,'PartA640 quadratic contrast'),
            'static_after_matched_old_linear464':static_matched(material,actual,primitive['strain']['hydro']['dK_depsilon_eV_A2'],0.,model,'464 linear hydro contrast'),
            'zero_point_Hartree_once_eV_A2':(3*A+2*B)*actual['local_mean_square_axis_A2_T0']})
    dependencies=['input.json','historical-inputs.json','historical-anharmonic.json','historical-refit.json','primitive-strains.json','partA-law.json','engine/anharmonic.py','engine/refit.py','engine/domain.py','engine/strain_tensor_law.py','pd_x1_refit.py']
    output={'source':'762 fixed-k3 local force curvature+joint100/111 quartics; no lattice or DFT rerun',
        'potential':'Kii sum(u_i^2)/2+A sum(u_i^4)/4+B sum_i<j(u_i^2 u_j^2)/2;B=(3K4_111-A)/2;u inA,m in eV/c^2',
        'assumption':'Fixed-k3 periodically clamped local well+matched scalar-J RPA; PdH mass substitution at4.084A; frozen local variance/quartics under PartA collective strain contrast at fixed extraelastic sumJ; physical reference friction held.',
        'quantum_normalization':'Full finite-basis dipole strengths enter chi_local(E); RPA=1/(chi_local^-1-J0); NO Hartree shift added again to that exact local spectrum.',
        'classical_counterterm':'Only (3A+2B)*variance_axis_T0 enters future Langevin parameters; thermal Hartree used only in the separate static closure, not added to sampled300K noise.',
        'static_reference':'Matched Kii-J0 + [K_PartA(epsilonI)-K_PartA(0)] + (3A+2B)*variance_axis(T). No absolute Gamma as onsite; fixed4.084A reference, not PBE-recentered.',
        'static_interpretation':'Published478 used primitive Gamma intercept with old linear464; same-reference old/new PartA rows isolate the actual force/quartic update. Hydro data only+/-1%; larger roots are extrapolations even inside the mixed-strain coordinate box.',
        'hydro_contrast_polynomial_eV_A2':{'linear':linear,'quadratic':quadratic,'zero_intercept':0.},
        'bath_policy':'Hold physical kappa=2/reference_amplitude_decay_ps, NOT bare gamma=.01. Atomic lattice gamma=kappa*hbar/E0; reduced pole bare-force gamma=M_pole*kappa*hbar/E0. A literal .01 changes physical bath; values for both conventions are explicit.',
        'spectral_mass_scope':'RPA pole mass is for single-pole action/drive normalization; do NOT replace the bare atomic mass in the full nonlinear coordinate lattice.',
        'rows':rows,'local_PdH_over_PdD':rows[0]['spectrum']['zero_to_one_meV']/rows[1]['spectrum']['zero_to_one_meV'],
        'RPA_Gamma_PdH_over_PdD':rows[0]['pumped']['mode_quantum_meV']/rows[1]['pumped']['mode_quantum_meV'],
        'historical_local_PdH_over_PdD':historical_spectra['Pd_H_over_D'],
        'INS_descriptor_ratio_55p8_over37':55.8/37.,
        'source_caveats':adopted['caveats'],
        'next_input':'Actual completed loading stiffness/per-existing-bond J_dd law; joint lattice/map batch remains held. Dynamic pump vertex/loss still not measured.',
        'provenance':{'archive_sha256':adopted['source_archive_sha256'],'profile_id':adopted['profile_id'],'files':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in dependencies}},
        'execution':{'elapsed_s':time.perf_counter()-started,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'numpy':np.__version__,'threads':1,'changed_well_probe':'basis8 versus10 only, same original implementation; no source SCF checks repeated'}}
    print(json.dumps(output,indent=2,allow_nan=False))


if __name__=='__main__':main()
