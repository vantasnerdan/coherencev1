"""Local Einstein spectrum + same-matrix neighbour-J gain consumer.
Run: python refit.py computed_inputs.json anharmonic.json adopted-strains.json > refit.json
Pumped local+J and primitive-strain Hartree branches are kept distinct.
"""
import json
import sys
from math import sqrt
from domain import HBAR, HC


def consume(row, spectra):
    q=next(r for r in spectra['rows'] if r['material']==row['name'])
    s=q['spectrum']
    E=s['RPA_Gamma_quantum_eV']
    if E is None:
        raise ValueError('No stable RPA Gamma pole for the selected gain consumer')
    Omega=E/HBAR
    Gamma2=1e12/row['reference_amplitude_decay_ps']
    w0=s['raw_zero_to_one_x_strength_A2']
    wG=s['RPA_Gamma_x_strength_A2']
    Z=wG/w0
    singlet=next(x for x in s['parallel_one_to_two_lines'] if x['block_degeneracy']==1)
    doublet=next(x for x in s['parallel_one_to_two_lines'] if x['block_degeneracy']==2)
    UA=singlet['spacing_change_meV']/1000
    UE=doublet['spacing_change_meV']/1000
    Uself=(UA+2*UE)/3
    Upair=(UA-UE)/3
    Ucross=s['cross_manifold_shift_meV']/1000
    # Leading dilute-mode projection, per N cells: U_Gamma=U_local*Z^2/N.
    Uproject=Uself*Z*Z
    rho=sqrt(3)*Gamma2*HBAR/abs(Uproject)  # g=2 Gamma2, delta=0, conservative pulling
    return dict(
        material=row['name'], source_entry=row['source_entry'],
        local_Kii_eV_A2=row['K2_eV_A2'], same_matrix_Jsum_eV_A2=row['J_sum_eV_A2'],
        harmonic_local_quantum_meV=q['harmonic_quantum_meV'],
        local_zero_to_one_meV=s['zero_to_one_meV'],
        parallel_one_to_two=s['parallel_one_to_two_lines'],
        perpendicular_one_to_two=s['perpendicular_one_to_two_lines'],
        N1_proxy_zero_to_one_meV=row['N1_proxy_0_to_1_meV_from448'],
        INS_reference_meV=row['reference_INS_meV'],
        RPA_static_denominator=s['RPA_static_denominator'],
        conditional_critical_loading=s['RPA_critical_J_eV_A2']/row['J_sum_eV_A2'],
        conditional_critical_density_m3=4e30/row['a_A']**3*s['RPA_critical_J_eV_A2']/row['J_sum_eV_A2'],
        threshold_scope='neighbour occupancy rescales J at fixed local well; not a computed loading phase diagram',
        selected_phase_domain_size_nm=None,
        mode_quantum_meV=1000*E, mode_frequency_THz=Omega/(2*3.141592653589793*1e12),
        stiffness_modulation_clock_THz=Omega/(3.141592653589793*1e12),
        incident_carrier_Ey_squared_THz=Omega/(2*3.141592653589793*1e12),
        effective_curvature_modulation_threshold=4*Gamma2/Omega,
        Gamma_mode_strength_A2=wG, Gamma_spectral_mass_over_bare=1/s['RPA_Gamma_oscillator_strength'],
        effective_pole_spring_eV_A2=E/(2*wG),
        action_normalization='n=N*A_site,peak^2/(4*Gamma_mode_strength_A2)',
        local_Uself_meV=1000*Uself, local_Upair_meV=1000*Upair,
        local_Ucross_meV=1000*Ucross, Gamma_Uself_meV_times_N=1000*Uproject,
        Gamma_Upair_meV_times_N=1000*Upair*Z*Z,
        Gamma_Ucross_meV_times_N=1000*Ucross*Z*Z,
        conservative_reference_gain_occupation_per_cell=rho,
        reference_gain_definition='g=2Gamma2, exact pair resonance, scalar self-pulling projection; tensor pair conversion retained separately',
        stationary_stored_energy_per_cell_eV=E*rho+.5*Uproject*rho*rho,
        equilibrium_ordering_energy_gain_eV=0.,
        pump_stop_amplitude_decay_ps=row['reference_amplitude_decay_ps'],
        local_strain_threshold=None, incident_intensity_threshold_W_m2=None, active_linewidth_Hz=None)


def static_branch(row, spectrum, primitive):
    """Frozen local quantum/thermal variance, cubic Hartree contraction."""
    A=row['K4_100_eV_A4']
    D=row['K4_111_eV_A4']
    B=(3*D-A)/2
    K0=primitive['K_zero_eV_A2']
    slope=primitive['strain']['hydro']['dK_depsilon_eV_A2']
    output=[]
    for T,key in ((0,'local_mean_square_axis_A2_T0'),
                  (300,'local_mean_square_axis_A2_T300')):
        variance=spectrum[key]
        correction=(3*A+2*B)*variance
        critical=-(K0+correction)/slope
        soft=[]
        for fraction in (0.,.5,.9,.99):
            eps=fraction*critical
            r=-slope*(critical-eps)
            quantum=sqrt(HC*HC/row['mass_eVc2']*r)
            soft.append(dict(linear_strain=eps, curvature_eV_A2=r,
                             undamped_soft_mode_THz=quantum/(2*3.141592653589793*HBAR*1e12)))
        just=None
        if critical<.05:
            eps=critical+.001  # explicitly 0.1 percentage point beyond onset
            r=slope*.001
            amplitude=sqrt(-r/D)  # weakest cubic direction is111
            just=dict(linear_strain=eps, direction='111',
                      displacement_magnitude_A=amplitude,
                      Cartesian_components_A=[amplitude/sqrt(3)]*3,
                      energy_gain_per_isotope_eV=-r*r/(4*D))
        output.append(dict(T_K=T, local_axis_variance_A2=variance,
            Hartree_correction_eV_A2=correction, effective_B_Hartree_eV_A4=A+2*B/3,
            critical_linear_strain=critical, crosses_below_5percent=critical<.05,
            just_beyond=just, soft_mode_approach=soft))
    return output


if __name__=='__main__':
    data=json.load(open(sys.argv[1]))
    spectra=json.load(open(sys.argv[2]))
    primitive=json.load(open(sys.argv[3]))['model_inputs']
    rows=[consume(r,spectra) for r in data['materials']]
    for row,result,q in zip(data['materials'],rows,spectra['rows']):
        p=next(x for x in primitive if x['material']==('NiH' if row['name']=='NiH' else 'PdD'))
        result['static_Hartree']=static_branch(row,q['spectrum'],p)
    isotope_difference=[rows[1]['static_Hartree'][i]['critical_linear_strain']-
                        rows[2]['static_Hartree'][i]['critical_linear_strain'] for i in (0,1)]
    print(json.dumps(dict(source=data['source'],rows=rows,
        Pd_local_H_over_D=spectra['Pd_H_over_D'],
        Pd_Gamma_H_over_D=rows[1]['mode_quantum_meV']/rows[2]['mode_quantum_meV'],
        Pd_H_minus_D_critical_strain_T0_T300=isotope_difference,
        assumption='401 local quartic wells plus same-matrix J in lattice RPA; leading dilute collective interaction projection, held reference decay; not SSCHA or converged local force constants.',
        strain='Static branch uses primitive K(strain) plus (3A+2B)<ux^2>; local variance frozen versus strain; linear extrapolation and fixed quartics, NOT full self-consistency or local-strain derivatives.',
        incident_formula='Ith=2*c*epsilon0*Keff*heff_th/(|alpha2_eff|*|tE|^2*overlap), for the actually projected finite-frequency Qx^2 Ey^2 vertex.',
        linewidth='Z2 parametric phase locking; needs actual pump noise and phase-switching rates, not prior U(1) linewidth.',
        next='Self-consistent lattice/SSCHA; finite-q dispersion; matched spectral self-energy; dynamic Raman/strain vertex and transmission; geometry/noise.'),
        indent=2,allow_nan=False))
