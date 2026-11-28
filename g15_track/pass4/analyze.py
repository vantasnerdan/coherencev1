"""Decision readout of the saved joined fields; never reruns a producer."""
import csv,json,math
import numpy as np
import run as m
import sources as s
import fields as f
import native


def rows(name):
    with (s.HERE/name).open() as stream:return list(csv.DictReader(stream))


def native_pump():
    """Actual source-supported quench projection, distinct from remote strain."""
    p=s.inputs()['lattice'];C=2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4'];mass=p['mass_eVc2']*s.QE/299792458.**2;out=[]
    for arm in s.ARMS:
        for ec in s.ENERGIES:
            for mesh in (16,12):
                material=native.Material(mesh)
                for sp in (0,1):
                    data=np.load(s.HERE/f'native_E{ec:g}_{arm["name"]}_sp{sp}_m{mesh}.npz');t=data['time_s']
                    selected=data['selected_mode_indices'];j=int(np.argmin(np.linalg.norm(data['k_m_inverse'][selected],axis=1)))
                    index=int(selected[j]);Q=data['variance_q_selected_A2'][:,j];reference=data['base_K'][index]+C*Q;delta=data['projected_deltaK_eVA2']
                    omega=np.sqrt(reference*s.QE/1e-20/mass);phase=data['reference_phase_selected_rad'][:,j];depth=delta/reference
                    duration=float(t[-1]-t[0]);h2=2*abs(np.trapezoid(depth*np.exp(-2j*phase),t))/duration
                    critical=(2/arm['opt']+4/arm['dephase'])/omega;peak=float(data['maximum_U_eV_formula'].max())
                    ks=data['k_m_inverse'][selected];_,v0,D0,_=native.dispersion(p,ks,arm['allocation'])
                    same_covariance=v0[None,:,:]*np.sqrt(D0[None,:]/(D0[None,:]+C*data['variance_q_selected_A2']))[:,:,None]
                    out.append(dict(arm=arm['name'],E_cm_keV=ec,source_species=sp+1,mesh=mesh,local_peak_U_eV_formula=peak,
                        local_peak_U_bare_Gamma_fraction=float(material.value(peak)/native.Material(16).reference),
                        propagated_bare_Gamma_reference_K_eVA2=native.Material(16).reference,source_mesh_native_reference_K_eVA2=material.reference,
                        projected_Gamma_deltaK_min=float(delta.min()),projected_Gamma_deltaK_max=float(delta.max()),
                        instantaneous_electronic_evolved_pole_fraction_max=float(abs(depth).max()),
                        finite_window_chirped_2Omega_stiffness_projection=float(h2),finite_window_exposure_s=duration,
                        periodic_Mathieu_reference_threshold_min=float(critical.min()),instantaneous_depth_over_periodic_threshold_max=float(np.max(abs(depth)/critical)),
                        maximum_actual_source_group_m_s=float(np.linalg.norm(data['group_velocity_selected_m_s'],axis=-1).max()),
                        electronic_group_change_fixed_covariance_max_m_s=float(np.linalg.norm(data['group_velocity_selected_m_s']-same_covariance,axis=-1).max()),
                        native_vs_no_e_velocity_covariance_difference_max_eV=float(abs(data['velocity_theta_excess_eV']-data['velocity_theta_no_e_excess_eV']).max()),
                        finite_pulse_scope='actual native electron quench; Fourier projection is not a periodic instability or independent quartic self-pump',
                        remote_scope='remote acoustic support ends at 6.8THz, below twice this evolved Gamma pole',
                        group_comparator_scope='same instantaneous covariance, electronic constitutive shift removed; not a second evolved off trajectory'))
    return out


def phase_support(arm,beam):
    data=np.load(s.HERE/f'conditioned_{beam}_{arm["name"]}.npz');numerator=0.;denominator=0.;rms=[];peak=(-1.,None)
    for parent,h in enumerate(data['history']):
        r=data['receiver_positions_m'][parent]-data['source_positions_m'][parent]
        focus2=np.exp(-np.sum((r[None,:]-h[:,12:15])**2,axis=1)/arm['radius']**2);weight=h[:,8]*focus2
        n=float(np.trapezoid(weight*h[:,4]**2,h[:,0]));d=float(np.trapezoid(weight,h[:,0]))
        if d<=0:raise ValueError('No supported coherent phase comparison')
        numerator+=n;denominator+=d;rms.append(np.sqrt(n/d))
        index=int(np.argmax(abs(h[:,4])));value=float(abs(h[index,4]))
        if value>peak[0]:peak=(value,dict(phase_peak_parent=parent,phase_peak_age_s=float(h[index,0]),
            phase_peak_receiver_coherent_energy_eV=float(weight[index]),phase_peak_support_relative_to_parent_peak=float(weight[index]/weight.max())))
    return dict(coherent_energy_weighted_phase_rms_rad=float(np.sqrt(numerator/denominator)),
        largest_parent_energy_weighted_phase_rms_rad=float(max(rms)),late_sampled_phase_max_rad=peak[0],**peak[1])


def main():
    p,K,dT,dE=m.mode_parameters();modes=np.load(s.HERE/'frozen/modes.npz');S=modes['strain'][:,:3].sum(axis=1)
    decision=[];source=[];inflight=[];outcomes=[]
    tr=json.loads((s.HERE/'frozen/tracks.json').read_text())
    for ec in s.ENERGIES:
        total=sum(fr*m.source_totals(tr[f'{ec:g}_{sp}'],sp,{}) for sp,fr in enumerate((.85,.15)))
        source.append(dict(E_cm_keV=ec,current_uA=s.current(ec)*1e6,parent_rate_s=s.current(ec)/s.QE,
            nuclei_rate_s=1.15*s.current(ec)/s.QE,incident_power_W=s.current(ec)*2000*ec,
            **{name+'_eV_per_parent':float(total[j]) for j,name in enumerate(f.CHANNELS[:9])},
            native_profile_exposure_s=43200 if ec<=2 else 7200,macro_first_depth_cell_m=2.5e-7))
    for arm in s.ARMS:
        for beam in s.BEAMS:
            fields=np.load(m.field_path(beam,arm['name']));window=np.load(s.HERE/f'phase_{beam}_{arm["name"]}.npz')
            temp=fields['t2_deltaT'][:,:,0];weight=fields['source_fraction'][:,:,0];strain=fields['t2_strain'][:,:,0]
            wave=window['signed_surface_strain'];time=window['map_time_s'];t=time-time.mean()
            design=np.column_stack((np.ones(len(t)),t));detrended=wave.reshape(len(t),-1)-design@np.linalg.lstsq(design,wave.reshape(len(t),-1),rcond=None)[0]
            freq=np.fft.rfftfreq(len(t),float(np.mean(np.diff(time))));spectrum=np.fft.rfft(detrended,axis=0)
            spectrum[(freq<5000)|(freq>250000)]=0.;resolved=np.fft.irfft(spectrum,n=len(t),axis=0).reshape(wave.shape)
            inside=weight>0;outside=~inside;dynamic=float(np.sqrt(np.mean(resolved*resolved)));maximum=float(resolved.max());minimum=float(resolved.min())
            quiet,_,_,_=native.thermal_pole(p,K,300.)
            warm,_,_,_=native.thermal_pole(p,K+s.ALPHA*dE*temp,300+temp)
            strained,_,derivative,_=native.thermal_pole(p,K+s.ALPHA*dE*temp+dE*strain/3,300+temp)
            thermal_fraction=np.sqrt(warm/quiet)-1;static_strain_fraction=np.sqrt(strained/warm)-1
            wave_fraction=dE*resolved*derivative[None,:,:]/(6*strained[None,:,:])
            row=dict(arm=arm['name'],beam=beam,seconds=2.,temperature_max_K=300+float(temp.max()),temperature_beam_weighted_K=300+float(np.sum(temp*weight)),
                thermal_optical_frequency_max_abs=float(abs(thermal_fraction).max()),static_signed_strain_frequency_max_abs=float(abs(static_strain_fraction).max()),
                resolved_8p33_250k_wave_strain_rms=dynamic,detrended_wave_peak=maximum,detrended_wave_hole=minimum,
                wave_frequency_modulation_max_abs=float(abs(wave_fraction).max()),wave_peak_over_rms=max(abs(maximum),abs(minimum))/max(dynamic,1e-300),
                outside_source_wave_rms=float(np.sqrt(np.mean(resolved[:,outside]**2))),inside_source_wave_rms=float(np.sqrt(np.mean(resolved[:,inside]**2))),
                resolved_window_exposure_s=float(window['exposure_s']),saved_phase_dt_s=float(np.mean(np.diff(time))),
                attribution='subtract static warm deformation and affine drift; 120us resolves 8.33–250kHz, not the entire retained 5–250kHz basis; a hot-spot map alone is not interference')
            row.update(phase_support(arm,beam))
            if arm['closure']=='nominal_retained':
                for control in ('temperature_off','strain_off','frozen','native_e_off'):
                    other=np.load(m.field_path(beam,arm['name'],control=control))
                    row[control+'_temperature_difference_max_K']=float(abs(other['t2_deltaT']-fields['t2_deltaT']).max())
                    row[control+'_K_difference_max_eV_A2']=float(abs(other['t2_K']-fields['t2_K']).max())
                    row[control+'_reentrant_strain_difference_max']=float(abs(other['t2_strain']-fields['t2_strain']).max())
            decision.append(row)
    old=rows('frozen/prior_inflight_yields.csv');names={'fast_Cu':'fast_acoustic','nominal_retained':'nominal_mixed','long_retained':'optical_moment_long'}
    for arm in s.ARMS:
        for ec in s.ENERGIES:
            for beam in s.BEAMS:
                selected=[r for r in old if r['scenario']==names[arm['closure']] and float(r['E_cm_keV'])==ec and r['beam']==beam and r['branch']=='p' and r['comparator']=='static40']
                value=sum((.85 if int(r['species'])==1 else .30)/1.15*float(r['struck_D_supported_only_per_D']) for r in selected)
                continued=sum((.85 if int(r['species'])==1 else .30)/1.15*float(r['per_incident_D_struck_D']) for r in selected)
                need=float(selected[0]['G12_required_extra_per_incident_D'])
                inflight.append(dict(E_cm_keV=ec,arm=arm['name'],beam=beam,branch='p',comparator='static40',
                    source_owned_supported_inflight_per_D=value,conditional_BH_continuation_inflight_per_D=continued,
                    supported_inflight_fraction_of_residual=value/need,conditional_inflight_fraction_of_residual=continued/need,
                    source='frozen exact Pass1 four-transfer stopping quadrature, unchanged source; not a new patch/pair consumer',
                    scope='moving-recoil precursor to the same loss-port graph; below0.5keV relative energy is conditional and is not a measured bound'))
    for ec in s.ENERGIES:
        built=[r for r in rows('yield_mixture.csv') if float(r['E_cm_keV'])==ec and r['branch']=='p' and r['comparator']=='static40']
        flight=[r for r in inflight if r['E_cm_keV']==ec]
        outcomes.append(dict(E_cm_keV=ec,built_field_fraction_residual_min=min(float(r['fraction_of_G12_residual']) for r in built),
            built_field_fraction_residual_max=max(float(r['fraction_of_G12_residual']) for r in built),
            supported_inflight_fraction_residual_min=min(r['supported_inflight_fraction_of_residual'] for r in flight),
            supported_inflight_fraction_residual_max=max(r['supported_inflight_fraction_of_residual'] for r in flight),
            interval_kind='named closure/allocation/beam range with explicit source mesh and coupling contrasts; not confidence or a contact bound'))
    shuffle=rows('phase_shuffle.csv');counts=sum(int(r['other_parent_count']) for r in shuffle);expected=sum(float(r['physical_expected_other_parents']) for r in shuffle)
    import yielding
    opportunities=[];calendar=f.Calendar();rate=s.current(2.5)/s.QE
    for arm in s.ARMS:
        for beam in s.BEAMS:
            temporal=yielding.clock_second(calendar) if beam=='hot_bursts' else 1.
            per_parent=rate*s.beam_moments(beam)[1]/1e-4*np.pi*(3*arm['radius'])**2*12*(arm['opt']+arm['ne'])*temporal
            opportunities.append(dict(arm=arm['name'],beam=beam,whole_history_Palm_window_opportunities_per_parent=per_parent,
                whole_foil_directed_spacetime_window_opportunities_2s=2*rate*per_parent,
                scope='same finite 3-radius/12-lifetime window as phase sampler; opportunities, not guaranteed coherent interference or a tail bound'))
    report=dict(decision='Thermal baseline change, static medium strain and nonlocal dynamic modulation are separately measured on the same rain; see discrimination.csv rather than inferring interference from energy maps.',
        source_seconds=2.,physical_optical_windows=len(shuffle),observed_independent_optical_overlap_parents=counts,expected_overlap_parent_count=expected,
        sparse_overlap_scope='finite Palm windows; zero sampled others does not exclude rare encounters elsewhere in the seconds-long foil',
        whole_history_spacetime_opportunities=opportunities,
        maximum_wave_modulation=max(r['wave_frequency_modulation_max_abs'] for r in decision),
        maximum_thermal_frequency_modulation=max(r['thermal_optical_frequency_max_abs'] for r in decision),
        maximum_outside_beam_wave_rms=max(r['outside_source_wave_rms'] for r in decision),
        yield_named_ranges=outcomes,
        high_band_phase_scope='seven bands retain directional energy/moments; microscopic coherent high-band phases are not inferred from source power',
        interpretation='Finite conditional coarse model; no native contact probability, rogue-tail extrapolation or measured 2Omega pump follows from these fields')
    m.table('discrimination.csv',decision);m.table('source_power.csv',source);m.table('inflight_yields.csv',inflight);m.table('yield_outcomes.csv',outcomes);m.emit('decision.json',report)
    m.table('native_pump.csv',native_pump())
    print(json.dumps(report),flush=True)
if __name__=='__main__':main()
