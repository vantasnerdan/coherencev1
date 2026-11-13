"""Qualified passive-loss scenarios; never updates the accepted lattice inputs."""
import json
import math
from pathlib import Path

HBAR_EV_PS=6.582119569e-4
C_CM_S=2.99792458e10
HC_EV_A=1973.269804


def table(parameters,inputs):
    rows=[]
    for raw in inputs['measured_rows']:
        name=raw['material'];p=parameters['materials'][name]
        A=p['K4_100_eV_A4'];B=(3*p['K4_111_eV_A4']-A)/2
        K=p['Kii_eV_A2']+(3*A+2*B)*p['zero_point_axis_variance_A2']-sum(s['J_sum_eV_A2'] for s in p['shells'])
        energy=math.sqrt(HC_EV_A**2*K/p['mass_eVc2'])
        omega=energy/HBAR_EV_PS
        kappa=2/p['reference_amplitude_decay_ps']
        row=dict(material=name,kind='measured loss UNKNOWN',mode='accepted lattice cold small-amplitude Gamma reference',
            model_energy_meV=1000*energy,model_frequency_THz=omega/(2*math.pi),stiffness_clock_THz=omega/math.pi,
            measured_kappa_ps_inverse=None,measured_Gamma2_ps_inverse=None,measured_threshold=None,
            baseline_assumed_kappa_ps_inverse=kappa,baseline_on_resonance_h_eff=2*kappa/omega,
            baseline_h_code_at_this_mode=2*kappa/omega,source=raw['positive_control']['source'],
            measured_conditions=raw['positive_control']['material_phase_loading'],
            measured_temperature_K=raw['positive_control']['temperature_K'],
            band_observation_not_loss=raw['positive_control']['observed'],
            lifetime_bound=None,reason=raw['why_not_eligible'])
        if name=='PdH':
            width_meV=4.8
            gamma2=width_meV/(2*1000*HBAR_EV_PS)
            row['lifetime_bound']=dict(kind='conditional width-as-ceiling ONLY, not established by the source',
                assumption=inputs['width_ceiling_assumption'],Gamma2_upper_ps_inverse=gamma2,
                T2_lower_ps=1/gamma2,T1_lower_ps=1/(2*gamma2),kappa_friction_upper_ps_inverse=2*gamma2,
                hypothetical_mode_energy_meV=55.8,on_resonance_h_eff_upper=2*width_meV/55.8,
                relation='Gamma2=1/T2=1/(2T1)+gamma_pure, gamma_pure>=0; threshold assumes the same simple Markov oscillator projection.')
        rows.append(row)
    legacy=inputs['legacy_exploratory']
    omega=2*math.pi*legacy['frequency_THz'];kappa=legacy['kappa_per_second']*1e-12
    legacy_h=2*kappa/omega
    rows.append(dict(material='NiH reference',kind='exploratory bare code gamma=0.01, NOT measured',
        mode=legacy['mode'],frequency_THz=legacy['frequency_THz'],stiffness_clock_THz=2*legacy['frequency_THz'],
        kappa_ps_inverse=kappa,Gamma2_ps_inverse=kappa/2,h_eff_threshold=legacy_h,
        h_bare_threshold=legacy['h_bare'],hybrid_metric=legacy['hybrid_metric'],source=legacy['source'],
        convention='lambda is the dimensionless canonical dynamical pole; kappa=-2 Re(lambda)*E0/hbar, Omega=Im(lambda)*E0/hbar. The input gamma is NOT kappa/Omega.'))
    calculated=inputs['calculated_PdH']
    hwhm=calculated['electron_HWHM_cm_inverse']+calculated['phonon_HWHM_cm_inverse']
    frequency=calculated['frequency_cm_inverse_approx']*C_CM_S*1e-12
    kappa=4*math.pi*C_CM_S*hwhm*1e-12
    calc_h=2*kappa/(2*math.pi*frequency)
    rows.append(dict(material='PdH',kind='calculated, NOT measured; conditional on-shell Markov surrogate',
        mode=calculated['mode'],frequency_THz_approx=frequency,stiffness_clock_THz_approx=2*frequency,
        energy_HWHM_cm_inverse=hwhm,energy_FWHM_meV=1000*HBAR_EV_PS*kappa,
        kappa_ps_inverse=kappa,Gamma2_ps_inverse=kappa/2,T1_ps=1/kappa,T2_ps=2/kappa,
        h_eff_threshold_approx=calc_h,h_bare_threshold=None,source=calculated['source'],
        frequency_qualification=calculated['frequency_qualification'],
        transfer='Not a NiH/PdD measurement or a calibrated Gamma loss; multi-peak/frequency-dependent response can invalidate one Markov pole.'))
    return dict(rows=rows,illustrative_scenario_h_eff_range=[legacy_h,calc_h],
        range_meaning='Two named scenarios, NOT confidence bounds or proved limits for actual material loss.',
        oscillator='qddot+kappa*qdot+[Omega_q^2+h_code*Omega_ref^2*P*cos(2*omega_p*t)]q=0',
        threshold='On resonance: h_eff=h_code*Omega_ref^2*P/Omega_q^2; h_eff,th=2*kappa/Omega_q; h_code,th=2*kappa*Omega_q/(Omega_ref^2*P). P is the mass-normalized pump-mode projection. At the lattice Gamma reference P=1 and h_code=h_eff.',
        hybrid_bare='For the earlier hybrid projection h_eff=delta_k/(M_dyn*nu^2), h_bare=delta_k/k_bare. They are not generally equal.',
        linewidth='Isolated homogeneous Lorentzian only: FWHM_E=2*hbar*Gamma2. Pure amplitude-friction case has Gamma2=kappa/2, T2=2/kappa and T1=1/kappa; extra pure dephasing prevents equating T1 and T2.',
        no_model_recalibration=True,incident_field_or_intensity_threshold=None,
        incident_reason='AC stiffness modulation is not an incident field. Transmission, mode-curvature coupling and accompanying DC shift are not supplied.',
        lookup='Bounded primary lookup complete; no qualifying measured homogeneous bulk optic loss established for any of the three materials.')


if __name__=='__main__':
    result=table(json.loads(Path('lattice_inputs.json').read_text()),json.loads(Path('loss_inputs.json').read_text()))
    print(json.dumps(result,indent=2,allow_nan=False))
