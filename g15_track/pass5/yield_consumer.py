"""Last consumer: replacement heat + named one-parent and connected pair Palm laws.

The target drift changes the full relative velocity, including flux normalization
by the ORIGINAL incident speed. The fixed adsorbed sheet is included quietly once.
No optical contact probability is assigned to these conditional distributions.
"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.special import roots_genlaguerre,roots_legendre,gamma
from common import QE,MD,KB,ENERGIES,BEAMS,inherited,emit,table,digest,calendar as frozen_calendar

GX,GW=roots_genlaguerre(8,.5)
MU,MW=roots_legendre(8)
WEIGHT=GW[:,None]/gamma(1.5)*MW[None,:]/2
COEFFICIENTS={'p':(5.5576e4,2.1054e2,-3.2638e-2,1.4987e-6,1.8181e-10),
              'n':(5.3701e4,3.3027e2,-1.2706e-1,2.9327e-5,-2.5151e-9)}


def kinetic(lab,theta,velocity):
    lab=np.asarray(lab,float);theta=np.asarray(theta,float);velocity=np.asarray(velocity,float)
    if np.any(lab<=0) or np.any(theta<=0):raise ValueError('Nonpositive physical probe energy/covariance')
    speed=np.sqrt(2*lab*1000*QE/MD)
    energy=MD*(np.sum(velocity[...,:2]**2,axis=-1)+(speed-velocity[...,2])**2)/(2*QE*1000)
    e=energy[...,None,None]/2;thermal=GX[:,None]*theta[...,None,None]/2000
    return e,thermal,e+thermal-2*np.sqrt(e*thermal)*MU[None,:]


def log_integrand(lab,theta,velocity,branch,U):
    _,_,er=kinetic(lab,theta,velocity)
    S=sum(c*er**i for i,c in enumerate(COEFFICIENTS[branch]))
    if np.any(S<=0) or np.any(er+U<=0):raise ValueError('Relative kernel leaves positive BH chart')
    # sigma(er)*v_rel/v_original: the sqrt(er) factor cancels analytically.
    return np.log(S)-.5*np.log(er+U)-31.397/np.sqrt(er+U)-.5*np.log(np.asarray(lab)[...,None,None]/2)+np.log(1e-31)


def log_ratio(lab,theta0,velocity0,theta1,velocity1,branch,U):
    e0,t0,er0=kinetic(lab,theta0,velocity0);e1,t1,er1=kinetic(lab,theta1,velocity1)
    v0=np.asarray(velocity0);v1=np.asarray(velocity1);dv=v1-v0
    speed=np.sqrt(2*np.asarray(lab)*1000*QE/MD)
    de=MD/(4*QE*1000)*(np.sum(dv*(v1+v0),axis=-1)-2*speed*dv[...,2])
    de=de[...,None,None];dt=GX[:,None]*(np.asarray(theta1)-np.asarray(theta0))[...,None,None]/2000
    denominator=np.sqrt(e1*t1)+np.sqrt(e0*t0)
    root_delta=np.divide(de*t0+e1*dt,denominator,out=np.zeros_like(denominator),where=denominator>0)
    der=de+dt-2*MU[None,:]*root_delta
    S0=sum(c*er0**i for i,c in enumerate(COEFFICIENTS[branch]))
    slope=sum(c*sum(er1**j*er0**(i-1-j) for j in range(i)) for i,c in enumerate(COEFFICIENTS[branch]) if i)
    q0=er0+U;q1=er1+U
    if np.any(q1<=0):raise ValueError('Changed relative kernel leaves positive BH chart')
    root0=np.sqrt(q0);root1=np.sqrt(q1)
    return np.log1p(der*slope/S0)-.5*np.log1p(der/q0)+31.397*der/(root0*root1*(root0+root1))


def rate(lab,theta,velocity,branch,U):
    return np.sum(np.exp(log_integrand(lab,theta,velocity,branch,U))*WEIGHT,axis=(-2,-1))


def difference(lab,theta0,velocity0,theta1,velocity1,branch,U):
    a=log_integrand(lab,theta0,velocity0,branch,U)
    delta=log_ratio(lab,theta0,velocity0,theta1,velocity1,branch,U)
    return np.sum(np.exp(a)*np.expm1(delta)*WEIGHT,axis=(-2,-1))


def connected(lab,quiet,p1,p2,joined,branch,U):
    zero=np.zeros((len(lab),3))
    a=log_ratio(lab,quiet,zero,*p1,branch,U)
    b=log_ratio(lab,quiet,zero,*p2,branch,U)
    c=log_ratio(lab,quiet,zero,*joined,branch,U)
    # Exact identity; do not subtract three first-order rates to recover a
    # second-order pair term.
    innovation=np.exp(a+b)*np.expm1(c-a-b)+np.expm1(a)*np.expm1(b)
    return np.sum(np.exp(log_integrand(lab,quiet,zero,branch,U))*innovation*WEIGHT,axis=(-2,-1))


def standard_error(contribution,groups=None):
    contribution=np.asarray(contribution);n=len(contribution)
    residual=contribution-contribution.sum()/n
    if groups is None:return float(np.sqrt(n/(n-1)*np.sum(residual*residual)))
    sums=np.bincount(groups,weights=residual);count=len(sums)
    if count<2:raise ValueError('Fewer than two independent physical probe clusters')
    return float(np.sqrt(count/(count-1)*np.sum(sums*sums)))


def baseline(model,inputs,tracks,ec,branch,U):
    """The original quiet loading and fixed sheet; no endpoint heat addition."""
    profile=next(r for r in inputs['profiles'] if r['case']['E_cm_keV']==ec)
    layers=np.asarray(profile['layers']);theta0=KB*profile['case']['T_K'];quiet=0.;sheet=0.
    for sp,fraction in ((0,.85/1.15),(1,.30/1.15)):
        track=tracks[f'{ec:g}_{sp}'];z=np.r_[0.,[r['depth_A']*1e-10 for r in track['timeline']]]
        energy=np.r_[track['initial_eV'],[r['energy_eV'] for r in track['timeline']]]/1000
        nodes=np.unique(np.r_[z,layers[:,1][layers[:,1]<z[-1]]]);mid=(nodes[1:]+nodes[:-1])/2
        lab=np.interp(mid,z,energy);nd=layers[np.minimum(np.searchsorted(layers[:,1],mid,side='right'),len(layers)-1),3]
        quiet+=fraction*float(model.moving_rate(lab,theta0,branch,U)@(nd*np.diff(nodes)))
        sheet+=fraction*profile['nads_m2']*float(model.moving_rate(energy[0],theta0,branch,U))
    return quiet+sheet,sheet


def build_acoustic_probe(out,path4,prefix=''):
    """Price the actual signed covariance response once, before cold consumption."""
    from live import Strain
    from receiver import packet_duration
    from common import L,H
    out=Path(out);meetings=np.load(out/(prefix+'meetings.npz'));sampler=Strain(out,path4);values=[];horizons=[];band_fields=[];field_cache={}
    for index,ec in enumerate(meetings['E_cm_keV']):
        ec=float(ec);beam=str(meetings['beam'][index]);allocation=str(meetings['allocation'][index]);closure=str(meetings['closure'][index])
        source_allocation=str(meetings['source_allocation'][index])
        times=np.minimum(meetings['probe_time_s'][index],2.);xyz=meetings['xyz_m'][index];temperature=meetings['probe_temperature_K'][index]
        key=(ec,beam,closure,source_allocation,xyz.tobytes(),times.tobytes(),temperature.tobytes())
        if key not in field_cache:
            field=sampler.bands(ec,beam,xyz,time_s=times,allocation=source_allocation,closure=closure)
            value=sampler.covariance(ec,beam,xyz,temperature,time_s=times,allocation=source_allocation,closure=closure,band_field=field)
            field_cache[key]=(value,field)
        value,field=field_cache[key];values.append(value);band_fields.append(field)
        g=sampler.green(ec,source_allocation,closure);profile=next(r for r in sampler.s.inputs()['profiles'] if r['case']['E_cm_keV']==ec)
        ndmax=float(np.max(np.asarray(profile['layers'])[:,3]))
        inverse=max((1+max(0.,g.loading[key][1])*ndmax)/item[1] for key,item in g.material.items())
        horizons.append((2*np.sqrt(2)*L+2*H)*inverse+H/4460+
            max(packet_duration(g.bandwidth[pi],2*np.pi*ch['frequency']) for pi,ch in enumerate(g.channels)))
        print(f'acoustic_probe E={ec:g} {closure} {allocation} {beam}',flush=True)
    np.savez_compressed(out/(prefix+'acoustic_probe.npz'),theta_delta_eV=np.array(values),bands_trace_complex=np.array(band_fields),physical_cluster_horizon_s=np.array(horizons),meeting_sha256=np.array(digest(out/(prefix+'meetings.npz'))),
        semantics=np.array('signed linear passive G3 covariance response from the same acoustic susceptibility/loss port; no new mean population heating'))
    return np.array(values)


def consume(out,path4,prefix='',sample_contributions=None):
    out=Path(out);s,native,fields,launch,response,acoustic,yielding=inherited(path4)
    model=s.load_file(Path(path4)/'run.py','accepted_g12_velocity_consumer')
    inputs=s.inputs();tracks=json.loads((Path(path4)/'frozen/tracks.json').read_text())
    residual={float(r['E_cm_keV']):r for r in inputs['residual']}
    reference=float(residual[2.5]['model_corrected_raw_reactions_per_true_C'])/float(residual[2.5]['R_model'])
    meetings=np.load(out/(prefix+'meetings.npz'));rows=[];checks=[];quiet_cache={};hole_rows=[];projection_samples=[]
    acoustic=np.load(out/(prefix+'acoustic_probe.npz'))
    if str(acoustic['meeting_sha256'])!=digest(out/(prefix+'meetings.npz')):raise ValueError('Acoustic probe population is not the selected meeting population')
    if 'source_row_index' in acoustic.files:
        if str(acoustic['source_acoustic_sha256'])!=digest(out/'acoustic_probe.npz'):raise ValueError('Frozen source probe arrays changed')
        source_acoustic=np.load(out/'acoustic_probe.npz');source_index=acoustic['source_row_index']
    else:
        source_acoustic=acoustic;source_index=np.arange(len(meetings['E_cm_keV']))
    report=json.loads((out/(prefix+'meetings.json')).read_text())
    import hashlib
    epoch_names=('meetings.py','medium.py','live.py','receiver.py','low.py','common.py','medium.npz','low_modes.npz')
    epoch=hashlib.sha256(''.join(digest(out/name) for name in (*epoch_names,*report.get('extra_scientific_sources',()))).encode()).hexdigest()
    epoch+=f':{report["probe_samples_per_case"]}:{report["geometry_event_samples_per_case"]}:{report["event_response_samples_per_case"]}'
    if report['scientific_source_epoch']!=epoch:raise ValueError('Selected meeting arrays belong to a different scientific operator')
    calendar=frozen_calendar(path4)
    prior=list(csv.DictReader((Path(path4)/'yield_mixture.csv').open()))
    for index,ec in enumerate(meetings['E_cm_keV']):
        ec=float(ec);beam=str(meetings['beam'][index]);closure=str(meetings['closure'][index]);allocation=str(meetings['allocation'][index])
        source_allocation=str(meetings['source_allocation'][index]);source_i=int(source_index[index])
        lab=meetings['probe_lab_keV'][index];q0=meetings['theta_quiet_eV'][index];dtheta=source_acoustic['theta_delta_eV'][source_i];q=q0+dtheta;zero=np.zeros((len(lab),3))
        column=meetings['probe_column_weight_m2'][index];weight=meetings['occupation_weight'][index]*column
        single=meetings['single_occupation_weight'][index]*column[:,None]
        ordinary_weight=column*calendar.value(meetings['probe_launch_time_s'][index],beam)/(2*meetings['probe_time_proposal_Hz'][index]*len(lab))
        required=(float(residual[ec]['R_observed'])-float(residual[ec]['R_model']))*reference*QE/1.15
        profile=next(r for r in inputs['profiles'] if r['case']['E_cm_keV']==ec);quiet_theta=KB*profile['case']['T_K']
        times=meetings['probe_time_s'][index];order=np.argsort(times);groups=np.empty(len(times),int)
        groups[order]=np.cumsum(np.r_[False,np.diff(times[order])>source_acoustic['physical_cluster_horizon_s'][source_i]])
        clusters=len(np.unique(groups))
        band_field=source_acoustic['bands_trace_complex'][source_i]
        coherent=abs(band_field.sum(axis=-1))**2;no_cross=np.sum(abs(band_field)**2,axis=-1)
        ratio=np.divide(coherent,no_cross,out=np.full_like(coherent,np.nan),where=no_cross>0)
        valid=no_cross>0;normalizer=float(ordinary_weight[valid].sum())
        if normalizer<=0:raise ValueError('Empty column-weighted probe projection measure')
        for cutoff in (.1,.5,1.):
            indicator=valid&(ratio<cutoff);fraction=float(ordinary_weight[indicator].sum()/normalizer)
            influence=ordinary_weight*(indicator-fraction*valid)/normalizer
            hole_rows.append(dict(E_cm_keV=ec,closure=closure,allocation=allocation,beam=beam,
                cross_band_projection_ratio_cutoff=cutoff,column_weighted_probe_exposure_fraction=fraction,
                conditional_cluster_MC_standard_error=standard_error(influence,groups),sample_count=len(lab),
                scope='|sum eight complex trace bands|² / sum |band|² on the same chart; column-weighted incident-probe exposure, not a time-only fraction, total-energy hole or nuclear event'))
        for j in range(len(lab)):
            point=meetings['xyz_m'][index,j]
            projection_samples.append(dict(E_cm_keV=ec,closure=closure,allocation=allocation,beam=beam,
                time_s=times[j],x_m=point[0],y_m=point[1],z_m=point[2],post250kHz_signed_trace_strain=float(band_field[j].sum().real),
                complex_trace_projection_intensity=coherent[j],no_cross_band_trace_projection_intensity=no_cross[j],
                ratio=ratio[j],column_weight_m_inverse2=ordinary_weight[j],
                scope='complex trace projection only; a negative real strain is not an intensity hole'))
        for branch in ('p','n'):
            for U,label in ((0.,'bare'),(.04,'static40')):
                key=(ec,branch,U)
                if key not in quiet_cache:quiet_cache[key]=baseline(model,inputs,tracks,ec,branch,U)
                quiet,sheet=quiet_cache[key]
                thermal_contribution=ordinary_weight*difference(lab,np.full(len(lab),quiet_theta),zero,q0,zero,branch,U)
                thermal_add=float(thermal_contribution.sum())
                ac=ordinary_weight*difference(lab,q0,zero,q,zero,branch,U)
                p=[];first=np.zeros(len(lab))
                for parent in (1,2):
                    th=meetings[f'theta_parent{parent}_eV'][index]+dtheta;v=meetings[f'v_parent{parent}_m_s'][index]
                    delta=difference(lab,q,zero,th,v,branch,U);p.append((th,v,delta));first+=single[:,parent-1]*delta
                outputs={}
                for control in ('joined','frozen','shuffled','strain_off','additive'):
                    theta=meetings['theta_'+('joined' if control=='strain_off' else control)+'_eV'][index]+dtheta
                    velocity=meetings['v_'+control+'_m_s'][index]
                    parent_velocity=[meetings[f'v_parent{j}_{control}_m_s'][index] if control in ('frozen','shuffled','strain_off') else p[j-1][1] for j in (1,2)]
                    pair=weight*connected(lab,q,(p[0][0],parent_velocity[0]),(p[1][0],parent_velocity[1]),(theta,velocity),branch,U)
                    if control=='joined':joined_contribution=pair
                    outputs[control]=(float(pair.sum()),standard_error(pair,groups))
                nonlin=weight*difference(lab,meetings['theta_additive_eV'][index]+dtheta,meetings['v_additive_m_s'][index],
                    meetings['theta_joined_eV'][index]+dtheta,meetings['v_joined_m_s'][index],branch,U)
                combined=thermal_contribution+ac+first+joined_contribution
                total=float(combined.sum());total_error=standard_error(combined,groups)
                if sample_contributions is not None:sample_contributions[(ec,closure,allocation,beam,branch,label)]=(combined,groups)
                selected=[r for r in prior if float(r['E_cm_keV'])==ec and r['arm']==closure+'_'+source_allocation and r['beam']==beam and r['mesh']=='16' and r['source_coupling']=='native_on' and r['branch']==branch and r['comparator']==label]
                if len(selected)!=1:raise ValueError('Ambiguous accepted replacement comparator')
                old=float(selected[0]['built_field_addition_per_D'])
                rows.append(dict(E_cm_keV=ec,closure=closure,allocation=allocation,beam=beam,branch=branch,comparator=label,
                    launch_source_allocation=source_allocation,Pass4_reference_allocation=source_allocation,
                    quiet_sheet_and_loading_per_D=quiet,fixed_quiet_sheet_per_D=sheet,replacement_thermal_delta_per_D=thermal_add,
                    thermal_MC_standard_error_per_D=standard_error(thermal_contribution,groups),
                    ordinary_signed_acoustic_delta_per_D=float(ac.sum()),ordinary_acoustic_MC_standard_error_per_D=standard_error(ac,groups),
                    acoustic_covariance_min_eV=float(dtheta.min()),acoustic_covariance_max_eV=float(dtheta.max()),
                    replacement_one_parent_delta_per_D=float(first.sum()),one_parent_MC_standard_error_per_D=standard_error(first,groups),
                    connected_joined_pair_delta_per_D=outputs['joined'][0],pair_MC_standard_error_per_D=outputs['joined'][1],
                    joined_minus_additive_delta_per_D=float(nonlin.sum()),joined_minus_additive_MC_standard_error_per_D=standard_error(nonlin,groups),
                    frozen_pair_delta_per_D=outputs['frozen'][0],shuffled_pair_delta_per_D=outputs['shuffled'][0],phase_only_strain_off_pair_delta_per_D=outputs['strain_off'][0],
                    Pass4_addition_per_D=old,replacement_total_addition_per_D=total,replacement_minus_Pass4_per_D=total-old,
                    G12_required_extra_per_D=required,fraction_of_G12_residual=total/required,G12_raw_6p5keV_reference_per_true_C=reference,
                    replacement_total_MC_standard_error_per_D=total_error,fraction_of_residual_MC_standard_error=total_error/abs(required),
                    independent_physical_probe_clusters=clusters,physical_cluster_horizon_s=float(source_acoustic['physical_cluster_horizon_s'][source_i]),
                    exposure_scope='incident-probe Cox-weighted mean over 0–2s; conditional chart MC error with within-probe covariance and retained-parent-window clustering; physical source/phase uncertainty is separate',
                    per_incident_probe_occupation=float(meetings['occupation_weight'][index].sum()),
                    single_occupation=float(meetings['single_occupation_weight'][index].sum()),
                    interpretation='conditional exterior/optical ray model; no inferred nuclear contact probability',
                    physical_struck_flight='unchanged source-owned accepted struck-D/flight support, separately tabulated'))
                if index==0 and branch=='p' and U==.04:
                    current=rate(lab,q,zero,branch,U);inherited_rate=model.moving_rate(lab,q,branch,U)
                    err=float(np.max(abs(current/inherited_rate-1)))
                    checks.append(dict(check='true_relative_kernel_zero_drift_matches_original',maximum_relative_error=err))
                    if err>5e-12:raise ValueError('Relative velocity/flux normalization changed inherited kernel')
                    # Independent inherited connected covariance implementation.
                    kk=yielding.Kernel(model,np.array([lab[0]]),float(q[0]),0.,float(q[0]),branch,U)
                    delta=float(q[0]*.05);a=float(difference(lab[:1],q[:1],zero[:1],q[:1]+delta,zero[:1],branch,U)[0]);b=float(kk.difference(0,q[0],delta))
                    rel=abs(a-b)/max(abs(a),1e-300);checks.append(dict(check='connected_difference_matches_accepted_Kernel',relative_error=rel))
                    if rel>1e-3:raise ValueError('Connected covariance difference changed accepted consumer')
    table(out/(prefix+'yield_outcomes.csv'),rows)
    table(out/(prefix+'probe_projection_holes.csv'),hole_rows)
    table(out/(prefix+'probe_projection_samples.csv'),projection_samples)
    emit(out/(prefix+'consumer_checks.json'),dict(checks=checks,case_rows=len(rows),meeting_sha256=digest(out/(prefix+'meetings.npz')),medium_sha256=digest(out/'medium.npz'),
        kernel='same 8x8 BH/Maxwell quadrature, stable expm1 connected differences, original incoming speed flux denominator',
        fixed_sheet_once=True,Pass4_native_addition_replaced_not_added=True,nuclear_contact='unknown',
        phase_controls='matched one-parent subtraction for frozen/shuffled/strain_off; phase-only at the same transported population/work, not all-force removal',
        below_0p5keV='accepted conditional BH continuation, not new measured cross section'))
    return rows


def merge_acoustic_probe(out):
    """Join saved disjoint probes and rebind them to the joined meeting arrays."""
    out=Path(out);values=[];horizons=[];band_fields=[]
    for ec in ENERGIES:
        prefix=f'E{ec:g}_'
        with np.load(out/(prefix+'acoustic_probe.npz')) as probe:
            if str(probe['meeting_sha256'])!=digest(out/(prefix+'meetings.npz')):raise ValueError('Acoustic shard changed its meeting population')
            values.append(probe['theta_delta_eV'])
            band_fields.append(probe['bands_trace_complex'])
            horizons.append(probe['physical_cluster_horizon_s'])
    np.savez_compressed(out/'acoustic_probe.npz',theta_delta_eV=np.concatenate(values),bands_trace_complex=np.concatenate(band_fields),physical_cluster_horizon_s=np.concatenate(horizons),
        meeting_sha256=np.array(digest(out/'meetings.npz')),
        semantics=np.array('same signed passive trace-strain susceptibility, disjoint saved energy probes; no mean population reheating'))
