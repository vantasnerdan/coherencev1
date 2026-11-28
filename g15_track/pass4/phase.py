"""Conditioned parent/receiver windows of the same whole-foil marked rain."""
import math
import numpy as np
import sources as s
import native
from fields import VOL,DX
HERE=s.HERE


def conditioned(arm,name,calendar,tracks,parameters):
    p,K,dT,g=parameters;window=np.load(HERE/f'phase_{name}_{arm["name"]}.npz')
    import envelope,launch,run
    network=run.network(arm)
    source_modes=np.load(HERE/'frozen/source_modes.npz');modes=np.load(HERE/'frozen/modes.npz')
    rng=np.random.default_rng(7941);counts=window['source_counts'];prob=counts.sum(axis=(1,2)).astype(float);prob/=prob.sum()
    traces=[];moments=[];shuffle=[];summaries=[];origins=[];receivers=[];arrivals=[];ledger=[];centroid_ages=[]
    for parent in range(64):
        frame=int(rng.choice(len(prob),p=prob));joint=counts[frame].astype(float);joint/=joint.sum()
        sp,node=np.unravel_index(rng.choice(joint.size,p=joint.ravel()),joint.shape)
        lo,hi,yl,yh=window['source_tile_bounds_m'][node]
        while True:
            sx=float(rng.uniform(lo,hi));sy=float(rng.uniform(yl,yh));local=float(s.beam(np.array([sx]),np.array([sy]),name)[0])
            if name=='even_steady' or rng.random()<local/9:break
        a=float(window['time_s'][frame]);b=a+float(window['dt_s']);target=rng.random()*float(calendar.integral(a,b,name))
        left=a;right=b
        for _ in range(32):
            middle=(left+right)/2
            if float(calendar.integral(a,middle,name))<target:left=middle
            else:right=middle
        at=(left+right)/2;tr=tracks[f'2.5_{sp}'];cache=np.load(HERE/f'native_E2.5_{arm["name"]}_sp{sp}_m16.npz')
        sites=float(cache['sites']);mass=p['mass_eVc2']*s.QE/299792458.**2*sites
        nm=len(modes['omega']);state=window['modal_phase_space'][frame];qmode=state[:,:nm];vmode=state[:,nm:]
        # All quantities are read before this bin's source insertion.
        strain_basis=source_modes['strain'][node,:3].sum(axis=0)
        strain=float(strain_basis@qmode[0]);static=float(strain_basis@window['static_reference_q'])
        lattice_velocity=float(source_modes['displacement'][node,2]@vmode[0])
        ix=min(int((sx+s.L/2)/DX),VOL.shape[0]-1);iy=min(int((sy+s.L/2)/DX),VOL.shape[1]-1)
        temperature=float(window['surface_temperature_K'][frame,ix,iy]);T=temperature-300
        selected=parent%len(cache['selected_mode_indices']);index=int(cache['selected_mode_indices'][selected])
        k=cache['k_m_inverse'][index].copy()
        if np.linalg.norm(k)<1:k=np.array([1.,0.,0.])/arm['radius']
        _,vg0,Db,_=native.dispersion(p,k[None,:],arm['allocation'])
        _,_,Dold=s.optical(p,k[None,:]);mass_D=mass/sites;M=mass_D*1e10/s.QE
        slow=g*(strain/3+s.ALPHA*T)
        cold,Qcold,_,_=native.thermal_pole(p,Db+slow,temperature)
        omega0=float(np.sqrt(cold[0]/M)*1e15);direct=.5*float(tr['components_eV'][2])
        raw=np.array(tr['components_eV'])*(sp+1);debit=arm['source_debits'][f'2.5_{sp}']
        raw[0]-=debit['electron'];raw[1]-=debit['host'];cutoff=tr['terminal_unresolved_eV']*(sp+1);share=arm['share']
        available=share*raw[0]+.05*raw[1]+(.5+.25*share)*raw[2]+(.25+.3*share)*cutoff
        impulse=np.sqrt(2*mass*direct*s.QE);delay=.74e-10/np.sqrt(2*tr['initial_eV']*s.QE/s.MD) if sp else 0.
        zs,ts,deposited,_,_,_=launch.track_arrays(tr);direct_weights=deposited[:,2]
        birth=float(ts@direct_weights/direct_weights.sum());depth=float(zs@direct_weights/direct_weights.sum())
        prior_normal=float(rng.standard_normal());vprior=lattice_velocity+np.sqrt(s.KB*temperature*s.QE/mass)*prior_normal
        qprior=np.sqrt(s.KB*temperature/(float(cold[0])*sites))*1e-10*float(rng.standard_normal())
        q=qprior;v=vprior;qr=qprior;vr=vprior;work=0.;control_work=0.
        qcontrol=qprior;vcontrol=np.sqrt(s.KB*temperature*s.QE/mass)*prior_normal
        for constituent in range(sp+1):
            if constituent:
                co=np.cos(omega0*delay);si=np.sin(omega0*delay)
                q,v=q*co+v/omega0*si,v*co-q*omega0*si;qr,vr=qr*co+vr/omega0*si,vr*co-qr*omega0*si
                qcontrol,vcontrol=qcontrol*co+vcontrol/omega0*si,vcontrol*co-qcontrol*omega0*si
            aa=direct;bb=float(v*impulse/s.QE);scale=1.
            if aa+bb>available-work:
                budget=available-work;scale=2*budget/(bb+np.sqrt(bb*bb+4*aa*budget))
            actual=aa*scale*scale+bb*scale;v+=impulse*scale/mass;work+=actual
            control_work+=direct+vcontrol*impulse/s.QE;vcontrol+=impulse/mass
            moments.append(dict(arm=arm['name'],beam=name,parent=parent,species=sp+1,constituent=constituent,
                global_arrival_s=at+birth+constituent*delay,source_x_m=sx,source_y_m=sy,source_depth_m=depth,
                read_before_insertion_velocity_m_s=float(v-impulse*scale/mass),actual_p_dot_v_plus_self_eV=actual,
                frozen_motion_self_eV=direct,impulse_scale=scale,total_parent_optical_budget_eV=available,
                remaining_parent_optical_incoherent_eV=available-work,live_lattice_velocity_m_s=lattice_velocity,
                medium_contribution_to_same_impulse_work_eV=float(impulse*scale*lattice_velocity/s.QE),D2_delay_s=float(delay),
                source_scope='direct-D coherent replacement within the same fixed optical budget; one molecular parent'))
        amplitude=complex(q-qr,(v-vr)/omega0);coherent=mass*omega0**2*abs(amplitude)**2/(2*s.QE)
        moments[-1].update(group_coherent_projection_eV=coherent,group_actual_work_eV=work,
            group_quiet_cross_eV=work-coherent,group_source_balance_eV=coherent+(work-coherent)+(available-work)-available,
            matched_frozen_motion_work_eV=control_work)
        norm=s.beam_moments(name)[0];duration=12*(arm['opt']+arm['ne'])
        expected=s.current(2.5)/s.QE*local/norm*np.pi*(3*arm['radius'])**2*duration*float(calendar.value(at,name))
        other=int(rng.poisson(expected));phases=rng.uniform(0,2*np.pi,other)
        phasor=amplitude+np.sum(np.exp(1j*phases))*np.sqrt(2*direct*s.QE/(mass*omega0**2))
        randomphase=float(rng.uniform(0,2*np.pi));shuffled=abs(amplitude)*np.exp(1j*randomphase)+np.sum(np.exp(1j*phases))*np.sqrt(2*direct*s.QE/(mass*omega0**2))
        incoherent=abs(amplitude)**2+other*2*direct*s.QE/(mass*omega0**2)
        shuffle.append(dict(arm=arm['name'],beam=name,parent=parent,physical_expected_other_parents=expected,other_parent_count=other,
            coherent_cross_impact_energy_eV=float((abs(phasor)**2-incoherent)*mass*omega0**2/(2*s.QE)),
            shuffled_cross_impact_energy_eV=float((abs(shuffled)**2-incoherent)*mass*omega0**2/(2*s.QE)),
            phase_shuffle_scope='one phase per electrical parent; molecular delay never shuffled independently'))
        ages=np.r_[0.,np.sort(rng.uniform(0.,12*arm['opt'],95))];clock=cache['time_s'];birth+=sp*delay
        _,Qcache0,_,_=native.thermal_pole(p,float(cache['base_K'][index]),300.)
        _,Qquiet,_,_=native.thermal_pole(p,Db,300.)
        offset=float(Qcold[0])-float(Qcache0)
        Q=np.interp(birth+ages,clock,cache['variance_q_selected_A2'][:,selected])+offset
        if np.any(Q<=0):raise ValueError('Nonpositive warm conditional covariance')
        delta=np.interp(birth+ages,clock,cache['projected_deltaK_eVA2'])
        shell=np.asarray(p['shells'][0]['periodic_alias_offsets_A'])*1e-10;F=float(np.cos(k@shell.T).mean())
        gradF=-np.einsum('b,bi->i',np.sin(k@shell.T),shell)/len(shell)
        C=2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4'];sensitivity=1. if arm['allocation']=='onsite' else F
        actual,check=envelope.trace(arm,cache,birth,ages,float(Db[0]),sensitivity,offset,slow,q-qr,v-vr,sites,mass_D,p['K4_100_eV_A4'],C,network)
        reference,reference_check=envelope.trace(arm,cache,birth,ages,float(Db[0]),sensitivity,float(Qquiet[0])-float(Qcache0),0.,q-qr,v-vr,sites,mass_D,p['K4_100_eV_A4'],C,network,frozen=True)
        ledger.append(check);coherent=check['initial_coherent_eV']
        moments[-1].update(group_coherent_projection_eV=coherent,group_quiet_cross_eV=work-coherent,
            group_source_balance_eV=coherent+(work-coherent)+(available-work)-available)
        stiffness=float(Db[0])+delta*sensitivity+slow+C*Q
        linear_omega=np.sqrt(stiffness/M)*1e15;omega=actual[:,9]
        grad=2*mass_D*np.sqrt(float(Db[0])/M)*1e15*vg0[0]/(s.QE/1e-20)
        nonlinear=3*p['K4_100_eV_A4']*actual[:,10]/(4*sites*stiffness**2)
        groups=(grad[None,:]+(delta[:,None]*gradF if arm['allocation']=='coupling' else 0.))*s.QE/1e-20/(2*mass_D*linear_omega[:,None])*(1-2*nonlinear[:,None])
        travel=np.vstack((np.zeros(3),np.cumsum((groups[:-1]+groups[1:])/2*np.diff(ages)[:,None],axis=0)))
        origin=np.array([sx,sy,depth]);receiver=origin+np.array([np.interp(arm['opt'],ages,travel[:,axis]) for axis in range(3)])
        r=receiver-origin;focus=np.exp(-np.sum((r[None,:]-travel)**2,axis=1)/(2*arm['radius']**2))
        origins.append(origin);receivers.append(receiver);arrivals.append(at);centroid_ages.append(birth+arm['opt'])
        correlation=np.exp(-ages/arm['dephase'])
        psi=(actual[:,0]+1j*actual[:,1]/linear_omega)*correlation*focus*np.exp(1j*float(k@r))
        canonical_velocity=linear_omega*psi.imag
        phase_difference=actual[:,8]-reference[:,8]
        weight=actual[:,2]*focus**2;weighted_phase=np.sqrt(np.trapezoid(weight*phase_difference**2,ages)/np.trapezoid(weight,ages))
        quiet_cross=(work-coherent)*np.exp(-ages/arm['opt'])
        remainder=actual[:,13]-quiet_cross
        if remainder.min()<-1e-10:raise ValueError('Coherent replacement exhausts its incoherent optical complement')
        history=np.column_stack((birth+ages,psi.real,psi.imag,canonical_velocity,phase_difference,np.linalg.norm(groups,axis=1),stiffness,omega,
            actual[:,2],actual[:,3],quiet_cross,remainder,travel,
            actual[:,4],actual[:,5],actual[:,6],actual[:,7],actual[:,10],actual[:,11]))
        traces.append(history)
        summaries.append(dict(arm=arm['name'],beam=name,parent=parent,global_arrival_s=at,species=sp+1,
            maximum_receiver_displacement_m=float(abs(psi.real).max()),maximum_receiver_velocity_m_s=float(abs(canonical_velocity).max()),
            maximum_medium_phase_rad=float(abs(phase_difference).max()),coherent_energy_weighted_phase_rms_rad=float(weighted_phase),
            maximum_group_speed_m_s=float(np.linalg.norm(groups,axis=1).max()),live_temperature_K=temperature,
            live_signed_strain=strain,static_reference_strain=static,source_birth_age_s=birth,
            maximum_conditional_energy_defect_eV=check['maximum_reciprocal_balance_eV'],minimum_conditional_carrier_eV=check['minimum_source_carrier_eV'],
            electron_to_coherent_work_eV=check['electron_work_eV'],incoherent_optical_to_coherent_work_eV=check['incoherent_optical_work_eV'],
            old_G3_matched_covariance_Gamma_K=float(Dold[0]),source_parent_hash=str(window['source_parent_hash']),
            phase_reference='same source covariance and physical impulse; cold/frozen propagation, native electronic constitutive shift removed',
            velocity_scope='canonical oscillator momentum/mass; coherent first moment, not derivative of an open-system expectation'))
    columns=['age_s','real_m','imag_m','velocity_m_s','medium_phase_rad','group_m_s','instantaneous_K_eVA2','omega_rad_s',
        'coherent_eV','dephased_remnant_eV','quiet_cross_eV','incoherent_remainder_eV','travel_x_m','travel_y_m','travel_z_m',
        'electron_work_eV','incoherent_optical_work_eV','heat_correction_eV','conditional_balance_eV','packet_optical_energy_eV','physical_loss_eV']
    np.savez_compressed(HERE/f'conditioned_{name}_{arm["name"]}.npz',history=np.array(traces),columns=np.array(columns),
        source_positions_m=np.array(origins),receiver_positions_m=np.array(receivers),global_arrival_s=np.array(arrivals),
        phase_ledger_schema=3,group_centroid_arrival_age_s=np.array(centroid_ages),conditional_ledger=np.array([list(x.values()) for x in ledger]),conditional_ledger_columns=np.array(list(ledger[0])),
        correlation_bandwidth_Hz=1/(2*np.pi*arm['dephase']),parent_ids=np.arange(64),physical_window_exposure_s=window['phase_sampling_exposure_s'])
    return summaries,moments,shuffle
