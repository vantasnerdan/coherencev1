"""Actual-path population, depth-conditioned signed motion, and joint pair/flight consumers.
The signed snapshot estimate is NOT promoted to a resolved optical-cycle exposure.
Its separate support envelope is the comparison that can make a replay relevant.
"""
import math
import numpy as np
from scipy.special import roots_legendre
import whole,motion

_PATH_KERNELS={}


def population_kernel(m,ec,sp,s,env,receiver,ch,U,population,N):
    """Reuse the fixed path kernel across force branches and matched controls."""
    key=(ec,sp,s['name'],env,receiver,ch,U,N,population.tobytes())
    if key not in _PATH_KERNELS:
        zz,column,labs,zgrid,egrid=receiver_rule(m,ec,receiver,m.p1.track(ec,sp)['timeline'][-1]['depth_A'])
        theta=m.p1.KB*m.p1.PROFILES[ec]['case']['T_K']
        cap=.36 if env=='bulk' else .24
        kinetic=.5*m.K*cap*cap/1000
        lo=(np.sqrt(labs)-np.sqrt(kinetic))**2
        hi=(np.sqrt(labs)+np.sqrt(kinetic))**2
        quiet=mode_rate(m,labs,0.,theta,ch,U)
        width=np.pi*s['radius']**2/N
        ordinary=np.empty(len(population));envelope=np.empty_like(ordinary)
        for i,p in enumerate(population):
            base=mode_rate(m,labs,p,theta,ch,U)
            ordinary[i]=width*np.sum(column*(base-quiet))
            lower=mode_rate(m,lo,p,theta,ch,U)*np.sqrt(lo/labs)-base
            upper=mode_rate(m,hi,p,theta,ch,U)*np.sqrt(hi/labs)-base
            envelope[i]=width*np.sum(column*np.maximum(abs(lower),abs(upper)))
        _PATH_KERNELS[key]=(zgrid,egrid,theta,ordinary,envelope)
    return _PATH_KERNELS[key]


def mode_rate(m,lab,population,theta,ch,U):
    """Gaussian excess-mode velocity covariance; quiet thermal motion once."""
    lab,population=np.broadcast_arrays(lab,population)
    effective_theta=(theta+population/3)[...,None,None]
    return m.p1.moving_rate(lab,effective_theta,ch,U)


def receiver_rule(m,ec,receiver,source_length_A):
    tr=m.p1.track(ec,receiver)
    z=np.r_[0.,[r['depth_A'] for r in tr['timeline']]]
    energy=np.r_[tr['initial_eV'],[r['energy_eV'] for r in tr['timeline']]]/1000
    limit=min(z[-1],source_length_A)
    edges=np.unique(np.r_[0.,limit,z,[r[1]*1e10 for r in m.p1.PROFILES[ec]['layers']]])
    edges=edges[(edges>=0)&(edges<=limit)];gl,gw=roots_legendre(4)
    zz=((edges[1:]+edges[:-1])[:,None]/2+np.diff(edges)[:,None]*gl/2).ravel()
    dz=(np.diff(edges)[:,None]*gw/2).ravel()*1e-10
    return zz,dz*m.profile_nd(ec,zz),np.interp(zz,z,energy),z,energy


def load_case(m,ec,sp,s,env):
    with np.load(m.HERE/f'peak_{ec:g}_{sp+1}_{s["name"]}_{env}.npz',allow_pickle=False) as archive:
        return {k:archive[k] for k in archive.files}


def velocity_measure(m,d,lane):
    # The two initial-energy moments are used ONLY to recover the physical
    # velocity-energy integral, never as a substitute for a local CSDA yield.
    moment=d['beam_joint_energy_moments_D_fs'][lane]
    a=d['beam_relative_energy_scale_keV']*moment[:,:,1].sum(axis=1)
    vb=np.sqrt(2*d['incoming_energy_keV']*1000/m.M)
    B=(a[1]-a[0])/(vb[0]-vb[1]);C=a[0]+vb[0]*B
    roundoff=64*np.finfo(float).eps*(abs(a[0])+abs(vb[0]*B))
    if C < -max(roundoff,1e-20):raise ValueError('Recovered kinetic-time measure is negative')
    return 1000*max(float(C),0.),float(1000*B/m.M),float(1000*roundoff)


def primary_integrals(m,d,records,ec,sp,s,env,branch,control,receiver,ch,U):
    lane=3*dict((name,i) for i,(name,_) in enumerate(m.BRANCHES))[branch]+control
    N=float(d['occupied_weight'].sum());population=d['population_energy_grid_eV_D']
    moments=d['beam_joint_energy_moments_D_fs'][lane,0,:,0]
    zgrid,egrid,theta,ordinary,envelope=population_kernel(m,ec,sp,s,env,receiver,ch,U,population,N)
    pop=float(moments@ordinary)*1e-15;bound=float(moments@envelope)*1e-15
    energy_time,displacement,roundoff=velocity_measure(m,d,lane)
    if energy_time==0:bound=0.
    depth=np.stack([r['pos'][:,2] for r in records])
    velocity=np.stack([r['v'] for r in records]).astype(float)
    dt=np.array([r['dt'] for r in records])
    weights=np.stack([r['w'] for r in records])*dt[:,None]
    pool=np.broadcast_to(np.array([r['pool'] for r in records])[:,None]/N,depth.shape)
    inside=(depth>=0)&(depth<=zgrid[-1])
    lab=np.interp(depth[inside],zgrid,egrid);v=velocity[inside]
    vb=np.sqrt(2*lab*1000/m.M)
    # Stable relative-energy change retains all observed directional signs.
    delta=m.M*(.5*np.sum(v*v,axis=1)-vb*v[:,2])/1000
    effective=lab+delta
    score=mode_rate(m,effective,pool[inside],theta,ch,U)*np.sqrt(effective/lab)-mode_rate(m,lab,pool[inside],theta,ch,U)
    observed=float(weights[inside]@score)*1e-15
    duration=float(dt[inside.any(axis=1)].sum())
    available=not (observed==0. and energy_time>max(roundoff,1e-25))
    return dict(population=pop,signed_snapshot_estimate=observed if available else None,
        observed_snapshot_value=observed,motion_magnitude_envelope=bound,kinetic_time_eV_fs=energy_time,
        signed_velocity_time_A_D=displacement,kinetic_roundoff_eV_fs=roundoff,
        exposure_fs=float(moments.sum()/N),snapshot_exposure_fs=duration,
        status='depth-conditioned snapshot quadrature, optical-cycle alias uncertainty bounded separately' if available else
            'saved signed samples do not resolve the nonzero per-step exposure; no apparent-zero replacement',
        population_model='phase-symmetric Gaussian excess modes, total energy P, three Cartesian modes per occupied D; quiet thermal once',
        envelope_model='local CSDA energy and accepted cutoff; monotone low-energy BH relative-flux kernel and actual bound-chart speed support')


def pair_rows(m,d,summary,ec,sp,s,env,branch,control,required,denom):
    rows=[];lane=3*dict((name,i) for i,(name,_) in enumerate(m.BRANCHES))[branch]+control
    ref=lane-control+2;weights=d['pair_weights'];Npairs=float(weights.sum())
    crossings=d['pair_crossings_by_block'];measures=d['pair_band_measures']
    source_fraction=.85 if sp==0 else .15
    for k,cut in enumerate(m.CUTS[env]):
        shell=float(cut['D_D_separation_A'] if env=='bulk' else cut['separation_A'])
        barrier=float(cut['total_approach_free_energy_eV'] if env=='bulk' else cut['AB_eV'])
        n=crossings[lane,:,k];nr=crossings[ref,:,k]
        flux=measures[lane,k];unsupported=flux[4]
        kinetic=flux[5]/flux[6] if flux[6]>0 else None
        for channel in ('p','n'):
            for eta in (0.,1e-12,1.):
                if eta==1:prob=(n>0).astype(float);prior=(nr>0).astype(float)
                else:prob=-np.expm1(n*np.log1p(-eta));prior=-np.expm1(nr*np.log1p(-eta))
                count=float(weights@(prob-prior))
                upper=min(Npairs,eta*flux[2]);supported_lower=max(0.,flux[1]-unsupported)
                for beam in m.p1.BEAMS:
                    rate=m.current(ec)/m.EV*source_fraction
                    rows.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,branch=branch,
                        control='built' if control==0 else 'matched_random',beam=beam,nuclear_branch=channel,
                        native_initial_separation_A=m.D0[env],shell_A=shell,native_static_cut_eV=barrier,
                        conditional_channel_hazard_given_available_shell_history=eta,
                        tagged_pairs_per_parent=Npairs,weighted_coherent_inward_crossings=float(weights@n),
                        coherent_history_hazard_reactions_per_parent=count,conditional_coherent_rate_s=rate*count,
                        ordinary_plus_motion_finite_band_flux=flux[0],band_angular_lower=flux[1],band_angular_upper=flux[2],
                        supported_band_flux_lower=supported_lower,unsupported_population_flux_upper=unsupported,
                        conditional_band_channel_union_upper_per_parent=upper,conditional_band_channel_rate_upper_s=rate*upper,
                        finite_band_mean_relative_kinetic_eV=kinetic,band_width_A=float(d['pair_band_width_A']),
                        coherent_per_incident_D=count/(sp+1),coherent_mixture_fraction_of_residual=(source_fraction*count/1.15)/required if required else None,
                        availability='unit tagged-motif prevalence endmember; actual bulk-pair/cage availability and its conditional contact probability remain unmeasured',
                        prevalence_rescaling='capped fields cannot generally be rescaled independently of their common optical resource',
                        event_definition='entry roots and momentum at the same linear-position interpolant time; tangencies excluded; no E/2 instantaneous collision assignment',
                        support='finite native shell is not nuclear contact; p/n channel hazards are alternatives and must not both be assigned unit probability',
                        long_time_validity='unit-hazard numbers are conditional conversion limits, not fixed-profile sustainable rates or measured cross sections'))
    return rows


def run(m,details,summary,samples):
    groups={}
    for r in details:groups.setdefault(r['key'],[]).append(r)
    for records in groups.values():records.sort(key=lambda r:r['t'])
    residual={float(r['E_cm_keV']):r for r in m.p1.DATA['residual']}
    denom=float(residual[2.5]['model_corrected_raw_reactions_per_true_C'])/float(residual[2.5]['R_model'])
    pass1=m.raw_table('yield_consumers.csv');index={(float(r['E_cm_keV']),r['scenario'],int(r['species'])-1,r['beam'],r['branch'],r['comparator']):r for r in pass1}
    work={(r['E_cm_keV'],r['species']-1,r['scenario'],r['environment'],r['branch']):r for r in summary}
    continuation=whole.continued_moving(m,summary,samples)
    last={}
    for r in sorted(samples,key=lambda x:x['time_fs']):
        last[r['E_cm_keV'],r['species']-1,r['scenario'],r['environment'],r['branch'],r['control']]=r
    beam={};pairs=[];rates=[];state=[];yieldrows=[];mixtures={};kernels=[]
    for ec in m.EC:
        required=(float(residual[ec]['R_observed'])-float(residual[ec]['R_model']))*denom*m.EV/1.15
        for s in m.p1.SCENARIOS:
            for env in m.D0:
                for sp in (0,1):
                    d=load_case(m,ec,sp,s,env)
                    for b,(branch,chi) in enumerate(m.BRANCHES):
                        rec=work[ec,sp,s['name'],env,branch]
                        for ctrl in (0,1):
                            records=groups[ec,sp,s['name'],env,branch,ctrl];lane=3*b+ctrl
                            pairs.extend(pair_rows(m,d,rec,ec,sp,s,env,branch,ctrl,required,denom))
                            for ch in ('p','n'):
                                for U,label in ((0.,'bare'),(.04,'static40')):
                                    for receiver in (0,1):
                                        key=(ec,s['name'],env,branch,ctrl,sp,receiver,ch,label)
                                        beam[key]=primary_integrals(m,d,records,ec,sp,s,env,branch,ctrl,receiver,ch,U)
                            for c,name in enumerate(m.p1.CHANNELS):
                                delta=float(d['carrier_integral_eV_fs'][-1,lane,c]-d['original_carrier_integral_eV_fs'][-1,c])
                                state.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,branch=branch,
                                    control='built' if ctrl==0 else 'matched_random',channel=name,
                                    final_original_eV=float(d['original_carriers_eV'][-1,c]),
                                    final_corrected_eV=float(d['corrected_carriers_eV'][-1,lane,c]),
                                    finite_window_inventory_shift_eV_fs=delta,
                                    true_current_mean_inventory_shift_eV=m.current(ec)/m.EV*(.85 if sp==0 else .15)*delta*1e-15,
                                    optical_mean_lifetime_fs=s['tau_opt']*1e15,
                                    thermal='original Pass1 foil thermal/electron/lattice/outflow history unchanged; this is signed phase-routing stock difference',
                                    tail='remaining carrier stocks retained; finite-window residence is not a new DC power source'))
                    for branch,chi in m.BRANCHES:
                        for ctrl in (0,1):
                            for receiver in (0,1):
                                for ch in ('p','n'):
                                    for U,label in ((0.,'bare'),(.04,'static40')):
                                        key=(ec,s['name'],env,branch,ctrl,sp,receiver,ch,label)
                                        kernels.append(dict(E_cm_keV=ec,source_species=sp+1,receiver_species=receiver+1,scenario=s['name'],environment=env,
                                            branch=branch,control='built' if ctrl==0 else 'matched_random',nuclear_branch=ch,comparator=label,**beam[key]))
                for branch,chi in m.BRANCHES:
                    for ctrl in (0,1):
                        for receiver in (0,1):
                            theta=m.p1.KB*m.p1.PROFILES[ec]['case']['T_K'];rec=work[ec,receiver,s['name'],env,branch]
                            observation=last[ec,receiver,s['name'],env,branch,'built' if ctrl==0 else 'matched_random']
                            moving=1-observation['mechanical_signed_work_eV']/max(rec['original_direct_D_eV'],1e-100)
                            if moving < -1e-10:raise ValueError('Direct-D consumer exceeds the existing signed source debit')
                            for ch in ('p','n'):
                                for U,label in ((0.,'bare'),(.04,'static40')):
                                    flight=motion.struck(m,ec,receiver,s,branch,ch,U,theta,env,ctrl)
                                    sibling=motion.sibling(m,ec,s,branch,ch,U,theta,env,ctrl) if receiver else dict(axial=0.,transverse=0.,axial_unsupported_column_fraction=0.,transverse_unsupported_column_fraction=0.)
                                    for beam_name in m.p1.BEAMS:
                                        ref=index[ec,s['name'],receiver,beam_name,ch,label]
                                        population=0.;estimate=0.;envelope=0.;available=True
                                        for source_sp,fraction in ((0,.85),(1,.15)):
                                            item=beam[ec,s['name'],env,branch,ctrl,source_sp,receiver,ch,label]
                                            factor=m.current(ec)/m.EV*fraction/m.p1.AREA*m.beam_factors(beam_name)
                                            population+=factor*item['population'];envelope+=factor*item['motion_magnitude_envelope']
                                            if item['signed_snapshot_estimate'] is None:available=False
                                            else:estimate+=factor*item['signed_snapshot_estimate']
                                        nuclear_weight=(.85 if receiver==0 else .30)/1.15
                                        old_sibling=float(ref['per_incident_D_same_parent_sibling'])
                                        chart=continuation.get((ec,receiver,s['name'],env,branch,ctrl,ch,label,0.),0.)
                                        chart_upper=continuation.get((ec,receiver,s['name'],env,branch,ctrl,ch,label,1.),0.)
                                        secondary=moving*(flight['recomputed_baseline']+flight['delta'])+chart
                                        observed_beam=population+estimate+old_sibling+sibling['axial'] if available else None
                                        row=dict(E_cm_keV=ec,species=receiver+1,scenario=s['name'],environment=env,branch=branch,
                                            control='built' if ctrl==0 else 'matched_random',beam=beam_name,nuclear_branch=ch,comparator=label,
                                            ordinary_population_actual_CSDA_per_D=population,
                                            signed_depth_snapshot_estimate_per_D=estimate if available else None,
                                            bound_chart_velocity_magnitude_envelope_per_D=envelope,
                                            same_parent_supported_axial_increment_per_D=sibling['axial'],
                                            same_parent_supported_transverse_increment_per_D=sibling['transverse'],
                                            molecular_axial_unresolved_column_fraction=sibling['axial_unsupported_column_fraction'],
                                            molecular_transverse_unresolved_column_fraction=sibling['transverse_unsupported_column_fraction'],
                                            pass1_same_parent_random_population_per_D=old_sibling,
                                            beam_supported_estimate_per_D=observed_beam,
                                            beam_supported_magnitude_upper_per_D=population+envelope+abs(old_sibling)+abs(sibling['axial']),
                                            pass1_random_beam_per_D=float(ref['per_incident_D_beam_moving_addition']),
                                            struck_D_supported_target_state_estimate_per_D=secondary,
                                            struck_D_signed_target_increment_per_D=moving*flight['delta'],
                                            struck_D_unresolved_target_baseline_per_D=moving*flight['unsupported_baseline'],
                                            pass1_struck_D_per_D=float(ref['per_incident_D_struck_D']),
                                            conditional_chart_continuation_per_D=chart,conditional_chart_release_bin_max_per_D=chart_upper,
                                            moving_energy_admitted_fraction=moving,struck_flight_max_fs=flight['flight_max_fs'],
                                            largest_target_candidate_speed_A_fs=flight['target_speed_A_fs'],
                                            incoming_nuclear_mixture_weight=nuclear_weight,G12_required_extra_mixture_per_D=required,
                                            beam_envelope_mixture_fraction_of_residual=nuclear_weight*(population+envelope+abs(old_sibling)+abs(sibling['axial']))/required if required else None,
                                            primary_reference='proton; own6.5keV absolute anchor; neutron and bare rows are sensitivities, not a fitted p/n stitch',
                                            normalized_supported_delta_R=nuclear_weight*observed_beam*1.15/m.EV/denom if observed_beam is not None else None,
                                            exposure_hours=12 if ec<=2 else 2,
                                            limitations='unresolved collisional/moving target states are not replaced by zero; supported estimate plus explicit support/envelope diagnostics, not a full joint nuclear contact prediction')
                                        yieldrows.append(row)
                                        key=(ec,s['name'],env,branch,ctrl,beam_name,ch,label)
                                        acc=mixtures.setdefault(key,dict(ordinary=0.,estimate=0.,estimate_available=True,envelope=0.,secondary=0.,chart_upper=0.,required=required))
                                        acc['ordinary']+=nuclear_weight*population;acc['envelope']+=nuclear_weight*row['beam_supported_magnitude_upper_per_D']
                                        acc['secondary']+=nuclear_weight*secondary;acc['chart_upper']+=nuclear_weight*chart_upper
                                        if observed_beam is None:acc['estimate_available']=False
                                        else:acc['estimate']+=nuclear_weight*observed_beam
    mixed=[]
    for key,a in mixtures.items():
        ec,s,env,b,ctrl,beam_name,ch,label=key
        mixed.append(dict(E_cm_keV=ec,scenario=s,environment=env,branch=b,control='built' if ctrl==0 else 'matched_random',beam=beam_name,
            nuclear_branch=ch,comparator=label,population_increment_mixture_per_D=a['ordinary'],
            supported_beam_estimate_mixture_per_D=a['estimate'] if a['estimate_available'] else None,
            supported_beam_magnitude_upper_mixture_per_D=a['envelope'],secondary_estimate_mixture_per_D=a['secondary'],
            conditional_chart_release_upper_mixture_per_D=a['chart_upper'],G12_required_extra_mixture_per_D=a['required'],
            beam_envelope_fraction_of_residual=a['envelope']/a['required'] if a['required'] else None,
            secondary_fraction_of_residual=a['secondary']/a['required'] if a['required'] else None,
            supported_beam_normalized_delta_R=a['estimate']*1.15/m.EV/denom if a['estimate_available'] else None,
            reference_raw_reactions_per_true_C=denom,physical_D_per_electrical_parent=1.15,
            meaning='0.85 D+ and 0.15 D2+ electrical parents; nuclear weights 0.85/1.15 and 0.30/1.15; no second charged-parent tally'))
    m.table('yield_consumers.csv',yieldrows);m.table('yield_mixture.csv',mixed);m.table('pair_consumers.csv',pairs)
    m.table('velocity_support.csv',kernels);m.table('state_difference.csv',state)
    return yieldrows,pairs,rates
