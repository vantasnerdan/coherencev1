"""Decision tables and figures from selected, already computed Pass 5 consumers."""
from pathlib import Path
import csv,json
import numpy as np
from common import QE,MD,KB,RHO,HERE,ENERGIES,BEAMS,inherited,emit,table,digest,calendar
from receiver import Green,receiver_nodes,BULK


def source_ledger(out,path4):
    out=Path(out);s,native,fields,launch,response,acoustic,yielding=inherited(path4);rows=[]
    for path in sorted(out.glob('spectrum_E*.npz')):
        stem=path.stem.removeprefix('spectrum_E');ec_text,closure=stem.split('_',1);ec=float(ec_text);archive=np.load(path)
        for sp in (0,1):
            accepted=np.load(Path(path4)/f'source_E{ec:g}_{closure}_sp{sp}.npz')
            total=float(accepted['energy_J'].sum())
            for pi,ch in enumerate(launch.ports()):
                a=archive[f'a_{sp}_{pi}'];w=archive[f'weight_{sp}_{pi}'];v=a[...,0]+a[...,2]+a[...,3]
                energy=float(np.sum(w[None,:]*(abs(v)**2+abs(a[...,1])**2)))
                mechanical=float(np.sum(w[None,:]*(np.real(np.conj(a[...,0])*v)+abs(a[...,1])**2)))
                heat=float(np.sum(w[None,:]*np.real(np.conj(a[...,2])*v)));electron=float(np.sum(w[None,:]*np.real(np.conj(a[...,3])*v)))
                target=float(accepted['energy_J'][pi].sum());random=archive[f'stochastic_{sp}_{pi}']
                variance=float(np.sum(w[None,:]*abs(a[...,1])**2));sample=float(np.sum(w[None,:]*np.mean(abs(random)**2,axis=-1)))
                rows.append(dict(E_cm_keV=ec,closure=closure,species=sp+1,port=pi,frequency_Hz=ch['frequency'],accepted_energy_J=target,
                    source_energy_fraction=target/total,replayed_energy_J=energy,relative_energy_defect=energy/target-1 if target else 0.,
                    mechanical_work_J=mechanical,late_LTE_work_J=heat,G16_electron_pressure_work_J=electron,
                    signed_work_identity_defect_J=energy-mechanical-heat-electron,stochastic_variance_relative_sample_error=(sample/variance-1) if variance else 0.,
                    source_scope='current exact Pass4 G16 source; old frequency table is shape sensitivity only; 64 common physical recoil realizations'))
    table(out/'source_projection.csv',rows);return rows

def finite_track_contrast(out,path4):
    """One resolved-primary line contrast; native/recoil/late-LTE are untouched.

    Frozen accepted +x radiation projection, exact source depths/flight times,
    path loading and phase. This is NOT an inferred native-pressure interior.
    """
    import copy
    from receiver import wavelet,packet_duration
    from medium import Medium
    from sampling import frequency_readout
    from yield_consumer import difference
    out=Path(out);g=Green(out,path4,2.5);medium=Medium(out)
    tracks=json.loads((Path(path4)/'frozen/tracks.json').read_text())
    normalization=np.load(Path(path4)/'frozen/native_transfer.npz')['normalization']
    nodes,_=np.polynomial.legendre.leggauss(12);rows=[]
    for sp in (0,1):
        tr=tracks[f'2.5_{sp}'];z,t,dep,E,dt,dp=g.launch.track_arrays(tr)
        dz=np.diff(np.r_[0.,[r['depth_A']*1e-10 for r in tr['timeline']]])
        echo=.74e-10/np.sqrt(2*tr['initial_eV']*QE/MD) if sp else 0.
        for pi in (5,6,7):
            ch=g.channels[pi];omega=2*np.pi*ch['frequency'];a=g.bandwidth[pi];duration=packet_duration(a,omega)
            nn=g.data[f'n_{sp}_{pi}'];direction=int(np.argmax(nn[:,0]));polar=g.data[f'polar_{sp}_{pi}']
            scale=np.sqrt(6*normalization[sp,pi]);source_k=g.data[f'k_{sp}_{pi}']
            primary=copy.copy(g);primary.data={k:g.data[k] for k in g.data.files}
            selected=np.zeros_like(g.data[f'a_{sp}_{pi}'])
            for pol in range(3):
                for j,n in enumerate(nn):
                    k=source_k[pol,j];phase=np.exp(1j*omega*t-1j*k*n[2]*z)*np.sinc((omega*dt-k*n[2]*dz)/(2*np.pi))
                    selected[pol,j,0]=scale[pol]*polar[pol,j,2]*np.sum(dp*phase)*np.exp(-.5*(k*2e-9)**2)*(1+np.exp(1j*omega*echo) if sp else 1)
            primary.data[f'a_{sp}_{pi}']=selected
            for radius in (1.3e-6,1e-5,1e-4,1e-3):
                point=np.array([.0025,-.0025,g.source_depth[sp]]);xy=point[:2]-np.array([radius,0.])
                baseline=list(g.components(point,xy,sp,0,False,pi,pi))
                far=list(primary.components(point,xy,sp,None,False,pi,pi))
                entries=[];phase_residual=[];loaded_residual=[];signed_work=[]
                rj=np.sqrt(radius*radius+(point[2]-z)**2)
                Dbar=np.divide(g.density.primitive(point[2])[0]-g.density.primitive(z)[0],point[2]-z,
                    out=g.density.local(z).copy(),where=abs(point[2]-z)>1e-15)
                for pol in range(3):
                    vp,vg,intrinsic,_,_=g.material[pi,pol];phase,group,loss=g.loading[pi,pol];k=omega/vp
                    host=rj/vg;integral=Dbar*host;arrival=host+group*integral
                    attenuation=np.exp(-np.minimum((intrinsic*host+loss*integral)/2,745.))
                    projection=polar[pol,direction,0];form=np.exp(-.5*(source_k[pol,direction]*2e-9)**2)
                    amplitude=scale[pol]*polar[pol,direction,2]*dp*form
                    local_nz=(point[2]-z)/rj
                    smooth=np.sinc((omega*dt-k*local_nz*dz)/(2*np.pi))
                    metric=projection/np.sqrt(4*np.pi*rj*rj*vg*RHO)/vp
                    for delay in ((0.,echo) if sp else (0.,)):
                        C=1j*np.sqrt(2)*metric*attenuation*amplitude*smooth*np.exp(1j*(k*rj-omega*arrival+phase*integral+omega*(t+delay)))
                        entries.extend(zip(arrival,t+delay,C))
                    residual=k*(z-point[2])**2/(rj+radius);phase_residual.append(residual)
                    loaded_residual.append(residual+phase*(integral-float(g.density.local(point[2]))*radius/vg))
                    atom=amplitude*np.exp(1j*omega*t)*np.sinc(omega*dt/(2*np.pi))*(1+np.exp(1j*omega*echo) if sp else 1)
                    total=g.data[f'a_{sp}_{pi}'][pol,direction,[0,2,3]].sum()
                    signed_work.append(np.real(np.conj(atom)*total)/6)
                ages=np.unique(np.concatenate([arrival+emission+(nodes+1)*duration/2 for arrival,emission,C in entries]+
                    [np.asarray(c[3])+(nodes+1)*duration/2 for c in baseline]))
                base=g.parent(point,xy,ages,sp,0,False,pi,pi)[0]
                removed=primary.parent(point,xy,ages,sp,None,False,pi,pi)[0]
                finite=np.zeros_like(base)
                for arrival,emission,C in entries:
                    finite+=C*wavelet(ages-arrival-emission,a,omega)*np.exp(-1j*omega*(ages-arrival))
                corrected=base+finite-removed
                field=medium.sample(2.5,'hot_bursts',point,time_s=2.)
                before=frequency_readout(medium,field,base.real,'onsite')[0]
                after=frequency_readout(medium,field,corrected.real,'onsite')[0]
                coefficient=g.response.covariance(g.arm,ch['frequency'],float(field[0]))['theta_per_linear_strain_eV']
                theta=KB*field[0];v=np.zeros((len(ages),3))
                delta_kernel=difference(np.full(len(ages),5.),theta+np.real(coefficient*base),v,
                    theta+np.real(coefficient*corrected),v,'p',.04)
                W=abs(np.array(signed_work));P=np.array(phase_residual);Q=np.array(loaded_residual)
                rows.append(dict(E_cm_keV=2.5,closure='nominal_retained',allocation='onsite',species=sp+1,port=pi,frequency_Hz=ch['frequency'],
                    radial_separation_m=radius,receiver_depth_m=point[2],source_terminal_depth_m=tr['terminal_depth_A'],
                    exact_geometry_phase_residual_max_rad=float(P.max()),
                    primary_absolute_signed_work_weighted_phase_rms_rad=float(np.sqrt(np.sum(W*P*P)/max(W.sum(),1e-300))),
                    primary_absolute_signed_work_weighted_loaded_phase_rms_rad=float(np.sqrt(np.sum(W*Q*Q)/max(W.sum(),1e-300))),
                    resolved_deposit_energy_weighted_geometry_phase_rms_rad=float(np.sqrt(np.mean(P*P@dep.sum(axis=1)/dep.sum()))),
                    primary_signed_work_at_native_angle_J=float(np.array(signed_work).sum()),
                    baseline_witness_max_abs_strain=float(abs(base.real).max()),corrected_witness_max_abs_strain=float(abs(corrected.real).max()),
                    maximum_known_primary_strain_correction=float(abs((finite-removed).real).max()),
                    baseline_witness_max_abs_Gamma_fraction=float(abs(before).max()),corrected_witness_max_abs_Gamma_fraction=float(abs(after).max()),
                    maximum_Gamma_correction_over_Pass4=float(abs(after-before).max()/2.5054719912489504e-8),
                    maximum_matched_ordinary_G12_kernel_change_m2=float(abs(delta_kernel).max()),retarded_witnesses=len(ages),
                    scope='cached complete source + finite primary stopping line - far primary stopping line, BEFORE squaring; fixed measured +x projection; native pressure, recoil and late-LTE remain uncorrected/unknown; witness contrast is not a rate or total yield'))
    medium.close();table(out/'finite_track_contrast.csv',rows);return rows





def phase_controls(out,path4):
    """One actual radiative impulse sign and guided-gauge operator witness."""
    import copy
    out=Path(out);g=Green(out,path4,2.5,control='zero_native')
    point=np.array([.001,0.,g.source_depth[0]]);source=np.zeros(2)
    unit=copy.copy(g);unit.data={k:g.data[k] for k in g.data.files}
    amplitudes=np.zeros_like(g.data['a_0_5']);amplitudes[2,:,0]=1.
    unit.data['a_0_5']=amplitudes
    c=list(unit.components(point,source,0,None,False,5,5))[2]
    omega=c[1];vp=unit.material[5,2][0]
    impulse_phase=float(np.angle(c[4]*np.exp(-1j*(omega/vp*.001-omega*c[3]))))
    if abs(impulse_phase-np.pi/2)>1e-10:raise ValueError('Wrong radiative force-to-strain impulse phase')
    target=np.array([.009,0.,1e-8]);before=list(g.components(target,source,0,0,False,0,0))[0][4]
    if abs(before)==0 or not g.guide_phase_scope[0,0,0]['known_normal_impulse_angular_nodes']:
        raise ValueError('Guided normal-impulse gauge witness has no supported response')
    profile=g.guides[0,0];rotation=np.exp(.37j)
    rotated=(profile[0],profile[1]*rotation,profile[2]*rotation,*profile[3:6],profile[6]*rotation)
    g.guides[0,0]=rotated
    for sp in (0,1):g.source_phase[0,sp][0]=g.guide_source_phase(0,0,sp,rotated)[0]
    after=list(g.components(target,source,0,0,False,0,0))[0][4]
    gauge_error=float(abs(after/before-1))
    if not np.isfinite(gauge_error) or gauge_error>1e-10:raise ValueError('Guided source/receiver eigenvector gauge is inconsistent')
    from receiver import guided_loading
    from types import SimpleNamespace
    density_value=1.23e28
    uniform_density=SimpleNamespace(edge=np.array([0.,g.s.H]),local=lambda z:np.full_like(z,density_value))
    h,effective=guided_loading(profile[6],uniform_density)
    expected_h=g.s.H*float(np.sum(abs(profile[6])**2/(2*np.arange(len(profile[6]))+1)))
    uniform_error=float(abs(h*effective/(density_value*expected_h)-1))
    if uniform_error>1e-12:raise ValueError('Uniform modal loading does not recover nD h')
    loss=Green(out,path4,2.5,control='zero_real');off=Green(out,path4,2.5,control='zero_native')
    ratios=[]
    for pi in (0,5):
        left=list(loss.components(point,source,0,None,False,pi,pi));right=list(off.components(point,source,0,None,False,pi,pi))
        ratios.extend(a[4]/b[4] for a,b in zip(left,right) if abs(b[4])>0)
    pure_loss_phase=float(max(abs(np.angle(x)) for x in ratios));pure_loss_gain=float(max(abs(x) for x in ratios))
    if pure_loss_phase>1e-10 or pure_loss_gain>1+1e-12:raise ValueError('Pure-loss loading creates phase or gain')
    proof=json.loads((out/'passive_controls.json').read_text())
    proof.update(bulk_radiative_unit_impulse_phase_rad=impulse_phase,guided_actual_operator_gauge_relative_error=gauge_error,
        uniform_guided_loading_relative_error=uniform_error,guided_h_Pa_inverse=h,guided_uniform_density_m_inverse3=density_value,
        guided_loading='15-node Gauss on every actual D(z) layer; h=integral |div(u)/omega|² dz with no /H; effective density applied once',
        zero_real_retained_loss_max_phase_rad=pure_loss_phase,zero_real_retained_loss_max_amplitude_ratio=pure_loss_gain,
        guide_source_scope={f'{pi}_{pol}_{sp}':value for (pi,pol,sp),value in g.guide_phase_scope.items()},
        source_modal_phase_map_unitarity_defect=float(max(np.max(abs(abs(a)-1)) for a in g.source_phase.values())),
        residual_guide_source='normal-impulse phase anchor where supported, gauge-covariant residual source-force scenario otherwise; not a calibrated pressure/recoil/LTE modal projection',
        phase_scope='isolated radiative branch, not an exact source-interior or kR~1 Green; cross-port and reflected interference is retained, not removed')
    emit(out/'passive_controls.json',proof)
    return impulse_phase,gauge_error


def sensitivities(out,path4):
    """Bounded changed-consumer contrasts; reuse the completed nominal PDF."""
    from sampling import marks,distribution,frequency_readout
    from low import Low
    from medium import Medium
    out=Path(out);raw=list(csv.DictReader((out/'E2.5_receiver_statistics.csv').open()))
    selected=next((i,r) for i,r in enumerate(raw) if r['closure']=='nominal_retained' and r['allocation']=='onsite' and r['beam']=='hot_bursts' and r['receiver']=='1')
    index,base=selected
    if base['scalar_green_sha256']!=digest(out/'receiver.py'):raise ValueError('Sensitivity baseline uses a different scalar Green')
    with np.load(out/'E2.5_receiver.npz') as data:history=data['history'][index]
    point=np.array([float(base[k]) for k in ('x_m','y_m','z_m')]);times=history[:,0];n=len(times)
    medium=Medium(out);low=Low(out,path4);field=medium.sample(2.5,'hot_bursts',np.broadcast_to(point,(n,3)),time_s=times)
    low_mean=low.mean(2.5,'hot_bursts',point,times);clock=calendar(path4)
    records=[dict(scenario='nominal_saved',bandwidth_scale=1.,propagation_reference_T_K=300.,quadrature_refinement=1,
        source_nodes=12800,Campbell_high_variance=float(base['Campbell_high_total_second']),
        Campbell_high_fourth=float(base['Campbell_high_total_fourth']),dense_far_third_absolute_ratio=float(base['dense_far_third_absolute_ratio']),
        dense_far_maximum_single_parent_strain=float(base['dense_far_max_parent_strain']),exact_near_mean_active=float(base['exact_near_mean_active']),
        Gamma_frequency_rms=float(base['fractional_Gamma_frequency_rms']),Gamma_frequency_min=float(base['fractional_Gamma_frequency_min']),Gamma_frequency_max=float(base['fractional_Gamma_frequency_max']),
        positive_fraction_Pass4=float(base['positive_time_fraction_pass4']),negative_fraction_Pass4=float(base['negative_time_fraction_pass4']),
        absolute_fraction_warm=float(base['absolute_time_fraction_warm']),physical_seconds=2.,sample_count=n)]
    for name,width,T,refine in (('half_width',.5,300.,1),('double_width',2.,300.,1),
                                ('nominal_warm_reference',1.,312.,1),('weak_mount_warm_reference',1.,335.,1),
                                ('doubled_nested_quadrature',1.,300.,2)):
        g=Green(out,path4,2.5,temperature=T,bandwidth_scale=width)
        mark=marks(g,point,'hot_bursts',refinement=refine);pdf=distribution(mark,clock,'hot_bursts',samples=n)
        fraction=frequency_readout(medium,field,pdf['strain']+low_mean,'onsite')[0]
        records.append(dict(scenario=name,bandwidth_scale=width,propagation_reference_T_K=T,quadrature_refinement=refine,
            source_nodes=len(mark['h']),Campbell_high_variance=float(mark['weight']@(mark['h']**2)),
            Campbell_high_fourth=float(mark['weight']@(mark['h']**4)),dense_far_third_absolute_ratio=pdf['far_standardized_third_absolute'],
            dense_far_maximum_single_parent_strain=pdf['far_maximum_single_parent'],exact_near_mean_active=pdf['exact_near_mean_active'],
            Gamma_frequency_rms=float(np.sqrt(np.mean(fraction*fraction))),Gamma_frequency_min=float(fraction.min()),Gamma_frequency_max=float(fraction.max()),
            positive_fraction_Pass4=float(np.mean(fraction>2.5054719912489504e-8)),negative_fraction_Pass4=float(np.mean(fraction< -2.5054719912489504e-8)),
            absolute_fraction_warm=float(np.mean(abs(fraction)>.00556844834632475)),physical_seconds=2.,sample_count=n))
    # Concrete passive port witness: zero real susceptibility, retained loss.
    source=point[:2]-np.array([1e-3,0.]);loss=Green(out,path4,2.5,control='zero_real');off=Green(out,path4,2.5,control='zero_native')
    left=list(loss.components(point,source,0,None,False,5,5));right=list(off.components(point,source,0,None,False,5,5))
    ratios=[a[4]/b[4] for a,b in zip(left,right) if abs(b[4])>0]
    phase=max(abs(np.angle(x)) for x in ratios);gain=max(abs(x) for x in ratios)
    if phase>1e-10 or gain>1+1e-12:raise ValueError('Pure-loss control creates real phase or gain')
    proof=json.loads((out/'passive_controls.json').read_text())
    proof.update(zero_real_retained_loss_max_phase_rad=float(phase),zero_real_retained_loss_max_amplitude_ratio=float(gain),
        scope='same isolated actual port/path; no claim that cross-port or reflected interference vanishes')
    emit(out/'passive_controls.json',proof)
    table(out/'parameter_sensitivity.csv',records)
    pilot=json.loads((out/'first_integrated.json').read_text())
    pilot_yields=list(csv.DictReader((out/'first_yield_outcomes.csv').open()))
    emit(out/'sampling_budget.json',dict(narrow_seconds=pilot['total_seconds'],live_probe_samples_per_case=16,live_event_response_samples_per_case=16,
        geometry_samples_per_case=1536,receiver_PDF_samples_per_case=n,physical_cores=6,threads_per_job=1,
        corrected_pilot_fraction_of_residual_range=[min(float(r['fraction_of_G12_residual']) for r in pilot_yields),max(float(r['fraction_of_G12_residual']) for r in pilot_yields)],
        corrected_pilot_maximum_fraction_MC_standard_error=max(float(r['fraction_of_residual_MC_standard_error']) for r in pilot_yields),
        decision='The corrected narrow 0–2s incident-probe consumer reaches G12; its conditional scale is far below the residual despite visible thermal sampling error. Preserve physical weights and report that error, rather than enlarging samples to imply source/phase closure.',
        retired_endpoint_pilot='2.83e-7 used heat at 2s, not the incident-probe exposure mean; it is not delivered as the exposure result.',
        refinement_variance_relative_change=records[-1]['Campbell_high_variance']/records[0]['Campbell_high_variance']-1,
        scope='conditional operator sensitivity, not source/phase completeness or a physical tail enclosure'))
    medium.close();return records


def merge_receivers(out):
    from low import Low
    out=Path(out);low=Low(out,out/'inherited');scalar_sha=digest(out/'receiver.py')
    for name in ('receiver_statistics.csv','band_cumulants.csv','signed_history.csv','matched_controls.csv'):
        paths=sorted(out.glob('E*_'+name),key=lambda p:float(p.name[1:].split('_',1)[0]));rows=[]
        for path in paths:
            current=list(csv.DictReader(path.open()))
            if name=='receiver_statistics.csv':
                if any(r['scalar_green_sha256']!=scalar_sha for r in current):raise ValueError('Receiver shard predates the selected signed Green')
                prefix=path.name.removesuffix(name)
                with np.load(out/(prefix+'receiver.npz')) as data:history=data['history']
                for i,row in enumerate(current):
                    ec=float(row['E_cm_keV']);xyz=np.array([float(row[k]) for k in ('x_m','y_m','z_m')])
                    case=low.case(ec,row['beam'],row['closure'],row['allocation'])
                    steady=float(low.projection(xyz)@(case['mean_force']/case['omega']**2))
                    dynamic=history[i,:,8]-steady
                    row.update(low_steady_elastic_beam_load_trace=steady,low_dynamic_Cox_mean_min=float(dynamic.min()),low_dynamic_Cox_mean_max=float(dynamic.max()))
            rows.extend(current)
        if rows:table(out/name,rows)
    return list(csv.DictReader((out/'receiver_statistics.csv').open()))


def limits(out,path4):
    from receiver import packet_duration
    out=Path(out);s,native,fields,launch,response,acoustic,yielding=inherited(path4);tracks=json.loads((Path(path4)/'frozen/tracks.json').read_text());rows=[]
    clock=calendar(path4)
    for arm in s.PHYSICAL_ARMS:
        for ec in ENERGIES:
            shape=Green(out,path4,ec,arm['name'])
            for sp in (0,1):
                tr=tracks[f'{ec:g}_{sp}'];source=np.load(Path(path4)/f'source_E{ec:g}_{arm["name"]}_sp{sp}.npz')
                a=max(2e-9,np.sqrt(2e6*arm['lambda_e']*arm['ne']),tr['terminal_depth_A'])
                for pi,ch in enumerate(launch.ports()):
                    omega=2*np.pi*ch['frequency']
                    represented=np.load(out/f'spectrum_E{ec:g}_{arm["name"]}.npz')
                    kval=represented[f'k_{sp}_{pi}'];energy_nodes=source['energy_J'][pi]
                    k=float(kval[energy_nodes>0].max())
                    radiation_k=min(shape.guides[pi,pol][3] for pol in range(3)) if ch['guided'] else min(omega/shape.material[pi,pol][0] for pol in range(3))
                    radiation_radius=10/radiation_k
                    radius=max(10*a,10*k*a*a,radiation_radius)
                    # Sufficient-radius scenario, not a measured interaction or
                    # a statement that the complete interior is unrepresentable.
                    density_bound=9/s.beam_moments('hot_steady')[0]
                    rate=min(s.current(ec)/QE,s.current(ec)/QE*np.pi*radius*radius*density_bound)
                    window=packet_duration(shape.bandwidth[pi],omega)
                    intensity=rate*window*(15 if ch['guided'] else 21)*11.7/clock.norm
                    volume=4*np.pi*(2e-9)**3/3;energy=(sp+1)*tr['initial_eV']*QE
                    rows.append(dict(E_cm_keV=ec,closure=arm['name'],species=sp+1,port=pi,frequency_Hz=ch['frequency'],
                        source_extent_hypothesis_m=a,sufficient_far_field_radius_scenario_m=radius,
                        radiative_kr10_radius_scenario_m=radiation_radius,
                        geometric_electrical_parent_rate_upper_inside_scenario_s_inverse=rate,
                        declared_packet_support_s=window,conditional_source_copy_active_intensity_upper=intensity,
                        conditional_active_fraction_upper=float(-np.expm1(-intensity)),
                        incident_parent_energy_J=energy,coarse_2nm_volume_mean_strain_incident_energy_cap=float(np.sqrt(2*energy/(BULK*volume))),
                        accepted_mean_port_energy_J=float(source['energy_J'][pi].sum()),
                        effect='No pointwise or physical tail exclusion; radius is sufficient, not an interaction volume. Duty refers only to the declared compact packet scenario.',
                        scope='worst represented positive-energy angular node; source-size/transport hypothesis, not a measured native aperture or unknown-phase bound'))
    table(out/'omitted_scales.csv',rows)
    emit(out/'parameters.json',dict(
        bandwidth_scale_alternatives=[.5,1.,2.],bandwidth_reference='old track integrated-frequency energy variance only, current native pressure phase fixed at each accepted port center',
        source_grain='accepted rotated Christoffel grain 27',recoil_catalog_realizations=64,
        propagation_susceptibility_reference_T_K=[300.,312.,335.],production_reference_T_K=300.,
        bulk_reference='isotropized per-branch sine dispersion using accepted angular speeds; exact anisotropic group-ray caustics are not represented',
        guide_radial_reference='outgoing Hankel radiative asymptotic; kr<10, finite-source interiors and nonradiative terms remain unclosed, not excluded',
        guide_source_phase='normal-impulse modal phase anchor; residual pressure/recoil/LTE phase is an energy-preserving scenario, not a calibrated modal force projection',
        force_to_strain='unconjugated +iwt source coefficients, outgoing +ikR-iwt; bulk +i strain derivative, outgoing guide +pi/4 with receiver derivative already present',
        screening_rate_endpoints_keV=[0.,.04],source_spectral_allocation='same accepted physical-closure source for matched target-native allocation scenarios',
        low_field='prompt and fast-pressure mean/covariance plus separate static elastic beam load; full low thermoelastic continuation unresolved',
        acoustic_intrinsic_laws=['Cu 32 dB/(cm MHz), energy rate vg ln(10)/10*32*f/MHz','retained 1/(50us)+vg/0.2m'],
        mount_edge_W_m2_K=[1e3,1e4,1e5],native_allocations=['onsite','coupling','native_Gamma_X'],physical_closures=['fast_Cu','nominal_retained','long_retained'],
        native_Gamma_X='f3114eb exact longitudinal increments; DeltaK=DeltaGamma+A(1-F1), shell-1 representative, linear low-A closure; frozen onsite incoming source, original pressure, no unique material vg',
        native_source_feedback='Only 2.5-keV spatial U snapshots retained; source regeneration and full U(t,r) donor-work history remain unresolved, not inferred from Umax or projected DeltaK',
        passive_source='same accepted mechanical + late LTE + G16 pressure port; projection is not another acoustic carrier',
        omissions=[
            'Resolved-primary line contrast does not reconstruct native pressure/recoil/LTE source interiors or unknown angular/within-band phase; finite recoil catalog and low thermoelastic continuation are not physical tail bounds.',
            'First-order retained acoustic operator uses reference susceptibility and a finite Gaussian exterior phase closure. Mean/direct variance use actual retarded Cox times; full-parent marginal quadrature retains first returns, but exterior cross-band/return covariance, mode conversion and full warm-path phase are not closed.',
            'Optical meetings are finite Gaussian RWA/eikonal support and leading pair factorial order; strain gradients over ps rays, close-k exchange, microscopic covariance/contact and ray-work backreaction are unretained.'],
        probability_policy='random-time fractions, conditional event Palm rates and per-incident-probe occupations have distinct denominators; no absence/realized nuclear probability inferred'))


def maps(out,path4):
    from live import Strain
    from medium import Medium
    out=Path(out);live=Strain(out,path4);medium=Medium(out);xyz=receiver_nodes();rows=[]
    for allocation in ('onsite','coupling'):
        for beam in BEAMS:
            time=1.999999973
            value=live.sample(2.5,beam,xyz,time_s=time,allocation=allocation);eps=np.trace(value,axis1=-2,axis2=-1)
            T,U,gradT,gradU=medium.sample(2.5,beam,xyz,time_s=time)
            for point,e,t,u in zip(xyz,eps,T,U):rows.append(dict(E_cm_keV=2.5,allocation=allocation,beam=beam,x_m=point[0],y_m=point[1],z_m=point[2],time_s=time,
                volumetric_strain=e,tensile_mean_stress_Pa=BULK*e,temperature_K=t,Ue_eV_per_formula=u,
                scope='same whole-foil Green/conditional exterior, nested actual coordinates to 0.125nm; not a patched propagation region'))
    table(out/'nested_receiver_maps.csv',rows);medium.close();return rows


def report(out=HERE):
    out=Path(out);receivers=list(csv.DictReader((out/'receiver_statistics.csv').open()));yields=list(csv.DictReader((out/'yield_outcomes.csv').open()))
    meetings=json.loads((out/'meetings.json').read_text())
    projected=json.loads((out/'nativeGX_meetings.json').read_text())
    main_cases=meetings['cases']+[r for r in projected['cases'] if r['allocation']=='native_Gamma_X']
    control_cases=[r for r in projected['cases'] if r['allocation']=='coupling_frozen_launch']
    closures=('fast_Cu','nominal_retained','long_retained');allocations=('onsite','coupling','native_Gamma_X')
    expected={(ec,c,a,b) for ec in ENERGIES for c in closures for a in allocations for b in BEAMS}
    def key(row):return float(row['E_cm_keV']),row['closure'],row['allocation'],row['beam']
    if {(key(r),int(r['receiver'])) for r in receivers}!={(k,i) for k in expected for i in range(4)} or len(receivers)!=648:raise ValueError('Incomplete receiver case grid')
    if {(key(r),r['branch'],r['comparator']) for r in yields}!={(k,p,c) for k in expected for p in ('p','n') for c in ('bare','static40')} or len(yields)!=648:raise ValueError('Incomplete G12 case grid')
    if {key(r) for r in main_cases}!=expected or len(main_cases)!=162 or len(control_cases)!=54:raise ValueError('Incomplete optical meeting/control case grid')
    table(out/'meeting_summary.csv',main_cases+control_cases)
    pdf={key(r):r for r in receivers if r['receiver']=='1'};summary=[]
    for r in yields:
        if r['branch']!='p' or r['comparator']!='static40':continue
        denominator=float(r['G12_required_extra_per_D']);field=pdf[key(r)]
        record={k:r[k] for k in ('E_cm_keV','closure','allocation','beam','fraction_of_G12_residual','fraction_of_residual_MC_standard_error','independent_physical_probe_clusters')}
        for name in ('replacement_thermal_delta_per_D','ordinary_signed_acoustic_delta_per_D','replacement_one_parent_delta_per_D','connected_joined_pair_delta_per_D'):
            record[name.removesuffix('_per_D')+'_fraction_of_residual']=float(r[name])/denominator
        for name in ('signed_volumetric_strain_min','signed_volumetric_strain_max','fractional_Gamma_frequency_rms','positive_time_fraction_pass4','negative_time_fraction_pass4','absolute_time_fraction_warm'):
            record['marginal_hotspot_'+name]=field[name]
        record['scope']='p branch, 40-eV comparator; 0–2s incident-probe yield versus separate hotspot marginal PDF, not a common realized crossing census'
        summary.append(record)
    table(out/'decision_summary.csv',summary)
    peak=max(receivers,key=lambda r:max(abs(float(r['fractional_Gamma_frequency_min'])),abs(float(r['fractional_Gamma_frequency_max']))))
    result=dict(base='f6010129f11ba8010815b78a67b3b6bed1db7194',conditional_model=True,
        receiver_rows=len(receivers),yield_rows=len(yields),meeting_cases=len(main_cases),matched_frozen_launch_control_cases=len(control_cases),energies_keV=list(ENERGIES),beam_shapes=list(BEAMS),
        maximum_sampled_abs_Gamma_frequency=max(abs(float(peak['fractional_Gamma_frequency_min'])),abs(float(peak['fractional_Gamma_frequency_max']))),
        maximum_sample_case={k:peak[k] for k in ('E_cm_keV','closure','allocation','beam','receiver','sample_count')},
        rare_tail_status='Unqualified: nominal nested cubature kappa2 ratio .654 and kappa4 ratio 8.33; finite diagnostics only, not P(tall), crossing-rate or material-exclusion bounds',
        maximum_abs_sampled_volumetric_strain=max(max(abs(float(r['signed_volumetric_strain_min'])),abs(float(r['signed_volumetric_strain_max']))) for r in receivers),
        maximum_sampled_over_Pass4=max(max(float(r['plus_over_pass4']),float(r['minus_over_pass4'])) for r in receivers),
        maximum_sampled_over_warm=max(max(float(r['plus_over_warm']),float(r['minus_over_warm'])) for r in receivers),
        G12_fraction_min=min(float(r['fraction_of_G12_residual']) for r in yields),G12_fraction_max=max(float(r['fraction_of_G12_residual']) for r in yields),
        maximum_conditional_fraction_MC_standard_error=max(float(r['fraction_of_residual_MC_standard_error']) for r in yields),
        exposure='incident-probe Cox-weighted time/depth measure over 0–2s; unchanged quiet sheet and own raw 6.5-keV reference',
        uncertainty='within-probe covariance and retained-parent physical clusters; conditional on the selected chart, excluding unclosed physical source/phase uncertainty',
        nuclear_contact='unknown; no assigned probability',acoustic_crossing_rate='not inferred from sparse random-time samples',
        physical_total_energy_holes='not determined from a scalar trace projection; probe_projection_holes.csv reports only cross-band complex-projection cancellation',
        event_rate='chronological distinct-parent first declared 3-sigma support intersection; not nuclear contact or the historical directed-window count',
        source_phase='bulk +i; guide +pi/4 with known normal-impulse modal anchor; residual guide-force phase is a scenario',
        conclusion='The conditional resolved model and the physical unclosed scales are separate. A small conditional yield is not an exclusion of source-interior or unmeasured phase physics.',
        native_projected_allocation='Exact f3114eb Gamma/X increments, named shell-1 continuation and low-A interpolation; same frozen onsite incoming source and original pressure, not a convex endpoint mixture',
        native_source_feedback=json.loads((out/'native_Gamma_X_readout.json').read_text()),
        consumer_identity={p.name:digest(p) for p in sorted(out.glob('*.py'))},
        output_sha256={name:digest(out/name) for name in ('medium.npz','meetings.npz','nativeGX_meetings.npz','acoustic_probe.npz','nativeGX_acoustic_probe.npz','yield_outcomes.csv','frozen_launch_controls.csv','receiver_statistics.csv','meeting_summary.csv','decision_summary.csv','probe_projection_holes.csv','probe_projection_samples.csv','passive_controls.json','finite_track_contrast.csv','native_Gamma_X_controls.json')})
    emit(out/'results.json',result)
    import os
    os.environ.setdefault('MPLCONFIGDIR',str(out/'.mplcache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for allocation in ('onsite','coupling','native_Gamma_X'):
        chosen=[r for r in yields if r['closure']=='nominal_retained' and r['allocation']==allocation and r['beam']=='hot_bursts' and r['branch']=='p' and r['comparator']=='static40']
        axes[0].plot([float(r['E_cm_keV']) for r in chosen],[float(r['fraction_of_G12_residual']) for r in chosen],'o-',label=allocation)
    axes[0].set(xlabel='Ecm (keV)',ylabel='Conditional addition / G12 residual',title='Native third: frozen incoming source');axes[0].legend()
    chosen=[r for r in receivers if float(r['E_cm_keV'])==2.5 and r['beam']=='hot_bursts' and r['receiver']=='1']
    axes[1].bar(np.arange(len(chosen)),[float(r['fractional_Gamma_frequency_rms']) for r in chosen]);axes[1].axhline(2.5054719912489504e-8,color='k',ls='--',label='Pass4 frequency modulation')
    axes[1].set(yscale='log',ylabel='Matched fixed-T Gamma frequency RMS',xlabel='closure / allocation');axes[1].legend();fig.tight_layout();fig.savefig(out/'consumer_comparison.png',dpi=150);plt.close(fig)
    return result


def native_landscape(out,path4):
    """Changed optical observer only; reuse accepted heat and signed maps."""
    from native_gammax import Law,NAME
    from medium import Medium,_trace
    from sampling import frequency_readout
    out=Path(out);medium=Medium(out,native_projection=Law(path4));rows=[]
    original=list(csv.DictReader((out/'medium_refraction.csv').open()))
    for old in original:
        if old['allocation']!='onsite':continue
        beam=old['beam'];at=2. if old['snapshot']=='two_seconds' else np.inf
        axis=2 if old['direction']=='depth' else 0
        r0=np.array([.0025 if axis==2 else .0035,-.0025,50e-9])
        q=np.array([float(old[k]) for k in ('initial_qx_m_inverse','initial_qy_m_inverse','initial_qz_m_inverse')])
        field=medium.sample(2.5,beam,r0,time_s=at);T,U,gT,gU=map(np.asarray,field)
        h,v,g=medium._hamiltonian(q,T,U,gT,gU,0.,0.,NAME,False)
        g_T=medium._hamiltonian(q,T,U,gT,np.zeros(3),0.,0.,NAME,False)[2]
        g_U=medium._hamiltonian(q,T,U,np.zeros(3),gU,0.,0.,NAME,False)[2]
        cold=medium.optical(2.5,beam,r0,q,time_s=at,allocation=NAME,frozen=True)[0]
        scan=np.broadcast_to(r0,(67,3)).copy();scan[:,axis]=np.linspace(0.,250e-9,67) if axis==2 else np.linspace(-.005,.005,67)
        w=medium.optical(2.5,beam,scan,q,time_s=at,allocation=NAME)[0]
        q0=q.copy();q0[axis]=0.
        turn=medium.optical(2.5,beam,scan,q0,time_s=at,allocation=NAME)[0]
        start=float(medium.optical(2.5,beam,r0,q0,time_s=at,allocation=NAME)[0])
        path,qpath,frequencies,reflections=_trace(medium,r0,q,beam,at,NAME,False)
        frozen_path,_,_,frozen_reflections=_trace(medium,r0,q,beam,at,NAME,True)
        vg0=medium.optical(2.5,beam,r0,q,time_s=at,allocation=NAME,frozen=True)[1]
        vg_path=medium.optical(2.5,beam,path,qpath,time_s=at,allocation=NAME)[1];length=np.linalg.norm(q);wave=2*np.pi/length
        row=old.copy();row.update(allocation=NAME,frozen_frequency_Hz=float(cold/(2*np.pi)),evolved_frequency_Hz=float(h/(2*np.pi)),
            scan_frequency_span_Hz=float(np.ptp(w)/(2*np.pi)),local_signed_kdot_normal_m_inv_s=float(g[axis]),kdot_thermal_normal_m_inv_s=float(g_T[axis]),
            kdot_electronic_normal_m_inv_s=float(g_U[axis]),local_evolved_vg_normal_m_s=float(v[axis]),local_frozen_vg_normal_m_s=float(vg0[axis]),
            one_lifetime_travel_m=float(np.linalg.norm(path[16]-r0)),evolved_minus_frozen_one_lifetime_normal_m=float(path[16,axis]-frozen_path[16,axis]),
            five_lifetime_relative_dq=float(np.linalg.norm(qpath[-1]-q)/length),five_lifetime_frequency_fraction_residual=float(np.ptp(frequencies)/h),
            physical_boundary_reflections=reflections,frozen_boundary_reflections=frozen_reflections,local_WKB_lambda_times_log_gradient=float(wave*np.linalg.norm(g)/h),
            local_1D_zero_normal_q_barrier_candidate=bool(np.any((h-start)*(h-turn[scan[:,axis]>=r0[axis]])<0)),
            normal_vg_reversed_in_five_lifetimes=bool(not reflections and np.any(vg_path[:,axis]*v[axis]<0)),
            native_scope='Actual common T/U gradients, exact Gamma/X inputs and named shell-1 representative; not a unique material vg or regenerated source')
        rows.append(row)
    table(out/'native_Gamma_X_refraction.csv',rows)
    maps=list(csv.DictReader((out/'nested_receiver_maps.csv').open()));selected=[(i,r) for i,r in enumerate(maps) if r['allocation']=='onsite']
    T=np.array([float(r['temperature_K']) for _,r in selected]);U=np.array([float(r['Ue_eV_per_formula']) for _,r in selected]);eps=np.array([float(r['volumetric_strain']) for _,r in selected])
    zeros=np.zeros((len(T),3));fraction,off,on,_,_=frequency_readout(medium,(T,U,zeros,zeros),eps,NAME)
    mapped=[]
    for j,(i,r) in enumerate(selected):
        mapped.append(dict(source_nested_map_row=i,allocation=NAME,beam=r['beam'],time_s=r['time_s'],x_m=r['x_m'],y_m=r['y_m'],z_m=r['z_m'],
            Gamma_fractional_frequency=fraction[j],strain_off_Omega_rad_s=off[j],strained_Omega_rad_s=on[j],
            source_scope='Exact alias of the onsite zero-U sound/heat map; recomputed native Gamma pole, not a new acoustic carrier'))
    table(out/'native_Gamma_X_maps.csv',mapped);medium.close();return rows
