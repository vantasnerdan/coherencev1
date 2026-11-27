"""Consumers of actual saved covariance histories, with matched quiet subtraction."""
import math
import numpy as np
from scipy.special import ndtr,ndtri,roots_genlaguerre,gamma
from scipy.interpolate import CubicHermiteSpline
from scipy.stats import qmc


def integrate(x,t,axis=0):return np.trapezoid(x,t,axis=axis)


def stable_rate_delta(m,lab,theta,delta,ch,U):
    """Actual relative-flux derivative, centered in the small measured covariance change.

    The centered difference step controls the derivative, not the physical pump.
    A direct nonzero extreme comparison is retained in consumer_checks.json.
    """
    lab=np.asarray(lab,float);theta=np.asarray(theta,float);delta=np.asarray(delta,float)
    mid=np.maximum(theta+delta/2,1e-10);h=np.maximum(mid*1e-3,1e-6)
    high=m.moving_rate(lab,mid+h,ch,U);low=m.moving_rate(lab,np.maximum(mid-h,1e-12),ch,U)
    return (high-low)/(2*h)*delta


def field_pairs(m,d,ec,sp,s,env,coupled,receiver):
    whole=receiver.integrate(m,d,ec,sp,s)
    d['whole_receiver']=whole
    t=whole['t'];pair=whole['field_pairs'];w=d['occupied']
    index=np.searchsorted(d['t'],t)
    evidence=m.HERE/f'field_pairs_{ec:g}_{sp+1}_{s["name"]}_{env}_{int(coupled)}.npz'
    np.savez_compressed(evidence,t_fs=t,anomalous_field_projections=pair,
        field_receiver_functions=whole['receiver_functions'],quadrature_position_A=whole['quadrature_position'],
        occupied_weight=whole['quadrature_weight'],source_owner=whole['owner'],retarded_delay_fs=whole['delay'],
        ordered_field_anomalous=whole['signal'],matched_constant_phase=whole['constant_phase'],
        local_mode_omega=d['omega'],local_mode_weight=d['mode_weight'],neighbor_characteristic=d['neighbor_cosine'])
    out=[];reference_omega=math.sqrt(float(d['localK'])/m.M);ordinary=d['optical']@w
    for b,coef in enumerate(m.COEF):
        local=whole['constant_phase'][:,b];global_signal=whole['signal'][:,b]
        phase_fraction=np.divide(abs(global_signal),local,out=np.full_like(local,np.nan),where=local>0.)
        roundoff_scale=32*np.finfo(float).eps*(len(whole['quadrature_weight'])+len(w))*float(local.max(initial=0.))
        pulse=coef*d['electron'];source_time=d['t']
        h2=2*abs(integrate(pulse*np.exp(-2j*reference_omega*source_time[:,None]),source_time))/max(source_time[-1],1.)
        beta=reference_omega/2*abs(integrate(pulse*np.exp(-2j*reference_omega*source_time[:,None]),source_time))
        global_integral=integrate(coef*whole['average_electron']*np.exp(-2j*reference_omega*t),t)
        threshold=2/(s['tau_opt']*1e15*reference_omega)
        last_start=np.max(d['source_clock']+d['source_duration'])+d['source_offsets'][-1]+whole['delay'].max()
        after=t>=last_start;ordered=abs(global_signal[after]);maxafter=ordered.max(initial=0.)
        peak_at=np.flatnonzero(after)[int(np.argmax(ordered))] if np.any(after) else 0
        remaining=np.flatnonzero((t>=t[peak_at])&(abs(global_signal)<=maxafter/math.e)) if maxafter else np.array([],int)
        persistence=t[remaining[0]]-t[peak_at] if len(remaining) else t[-1]-t[peak_at]
        determinant=d['qr']*d['vr']-d['cr']**2
        trace=(d['vr']*d['q'][:,b]+d['qr']*d['v'][:,b]-2*d['cr']*d['c'][:,b])/determinant
        ratio=(d['q'][:,b]*d['v'][:,b]-d['c'][:,b]**2)/determinant
        amplification=(trace+np.sqrt(np.maximum(trace*trace-4*ratio,0.)))/2
        for j in np.unique(np.r_[np.linspace(0,len(t)-1,48).astype(int),peak_at]):
            k=index[j]
            out.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,neighbor_on=coupled,
                coefficient_per_eV=coef,time_fs=t[j],global_anomalous=abs(global_signal[j]),local_anomalous=local[j],
                matched_constant_phase_anomalous=local[j],phase_defined=bool(local[j]>0.),
                phase_aligned_second_moment_fraction=float(phase_fraction[j]) if local[j]>0 else None,
                phase_absolute_roundoff_scale=roundoff_scale,
                phase_support='zero_norm' if local[j]==0 else ('positive_supported' if local[j]>roundoff_scale else 'roundoff_dominated'),
                phase_triangle_excess=max(float(abs(global_signal[j])-local[j]),0.),
                anomalous_field_projection_Frobenius=float(np.linalg.norm(pair[j,b])),
                source_budget_stored_delta_eV=d['energy'][k,b],signed_work_eV=d['work'][k,b],
                ordinary_optical_population_eV=ordinary[k],
                excess_over_ordinary_population=d['energy'][k,b]/ordinary[k] if ordinary[k]>1e-30 else None,
                actual_covariance_amplification=float(amplification[k].max()),reference_omega_per_fs=reference_omega,
                projected_signed_h2_max=float(h2.max()),maintained_Mathieu_h_threshold=threshold,
                projected_h2_over_threshold=float(h2.max()/threshold),finite_pulse_weak_beta_max=float(beta.max()),
                projected_signed_h2_whole=2*abs(global_integral)/max(t[-1],1.),
                finite_pulse_weak_beta_whole=reference_omega/2*abs(global_integral),
                actual_unfiltered_peak_fraction=float(abs(pulse).max()),post_source_1e_persistence_fs=persistence,
                complete_retarded_window_end_fs=whole['complete_window_end_fs'],
                retarded_quadrature_count_error_D=whole['weighted_count_error_D'],
                pure_dephasing_extra_per_fs=0.,independent_extra_dephasing_0p5ps_factor=math.exp(-2*t[j]/500.),
                claim='causal retarded-age whole-footprint ordered second moments; not maintained gain or entanglement'))
    return out


def instantaneous_pairs(m,d,ec,sp,s,env,coupled):
    info=m.PAIR[env];out=[];t=d['t'];weights=d['occupied'];hot=int(np.argmax(d['electron'].max(axis=0)))
    median=int(np.argsort(d['electron'].max(axis=0))[len(weights)//2])
    for site,label in ((hot,'near_track_max'),(median,'median_strength')):
        quietS=d['relq'][:,0,site];quietV=d['relv'][:,0,site];quietC=d['relc'][:,0,site]
        for b,coef in enumerate(m.COEF):
            S=d['relq'][:,b,site];V=d['relv'][:,b,site];C=d['relc'][:,b,site]
            for u,energy in zip(info['u'][1:],info['energy'][1:]):
                closing=2*u
                prob=ndtr(-closing/np.sqrt(S));quiet=ndtr(-closing/np.sqrt(quietS))
                j=int(np.argmax(prob));meanv=C/S*closing
                variance=np.maximum(V-C*C/S,1e-30);sig=np.sqrt(variance)
                inward=sig*np.exp(-meanv*meanv/(2*variance))/math.sqrt(2*np.pi)+meanv*ndtr(meanv/sig)
                density=np.exp(-closing*closing/(2*S))/np.sqrt(2*np.pi*S)
                flux=inward*density
                qmean=quietC/quietS*closing;qvar=np.maximum(quietV-quietC*quietC/quietS,1e-30)
                qsig=np.sqrt(qvar);qflux=(qsig*np.exp(-qmean*qmean/(2*qvar))/math.sqrt(2*np.pi)+qmean*ndtr(qmean/qsig))*np.exp(-closing*closing/(2*quietS))/np.sqrt(2*np.pi*quietS)
                out.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,neighbor_on=coupled,
                    coefficient_per_eV=coef,site=label,radius_A=d['position'][site,0],depth_A=d['position'][site,2],
                    supplied_separation_A=info['d0']-closing,native_approach_energy_eV=energy,
                    maximum_instantaneous_probability=prob[j],quiet_probability_at_same_time=quiet[j],
                    signed_probability_excess=prob[j]-quiet[j],ratio_at_same_time=prob[j]/quiet[j] if quiet[j]>0 else None,
                    time_of_maximum_fs=t[j],conditional_relative_kinetic_eV=m.M*(variance[j]+meanv[j]**2)/4,
                    inward_shell_crossings=integrate(flux,t),matched_quiet_crossings=integrate(qflux,t),
                    signed_shell_crossing_excess=integrate(flux-qflux,t),
                    negative_signed_separation_probability=float(ndtr(-info['d0']/np.sqrt(S)).max()),
                    probability_kind='instantaneous Gaussian relative distribution, not trajectory minimum',
                    native_support=info['minimum_native_A'],continuation_to_contact='unknown; no native contact probability or cage prevalence inferred'))
    return out


def trajectory_minima(m,d,ec,sp,s,env):
    """Paired propagated phase-space draws in the native directional normal form.

    Stop each path at the edge of the supplied cut. Censoring is carried into
    the minimum distribution, not replaced by a guessed contact potential.
    """
    info=m.PAIR[env];t=d['t'];hot=int(np.argmax(d['electron'].max(axis=0)))
    xnodes=2*info['u'];energy=info['energy'];derivative=np.zeros(len(xnodes))
    derivative[1:]=np.gradient(energy,xnodes)[1:]
    if env=='vacancy':
        derivative[1:]=np.array([r['AB_dE_du_eV_A']/2 for r in m.old.VAC['rows']])
    # Bulk total derivatives come from its total curve, never the separately
    # supplied connected-energy derivative.
    curve=CubicHermiteSpline(xnodes,energy,derivative)
    cutoff=xnodes[-1];mu=m.M/2;N=512
    normal=ndtri(np.clip(qmc.Sobol(2,scramble=True,seed=1537).random_base2(9),1e-12,1-1e-12))
    S0=d['relq'][0,0,hot];V0=d['relv'][0,0,hot]
    rng=np.random.default_rng(15371);wr=math.sqrt(info['krel']/mu)
    x=np.broadcast_to(normal[:,0]*math.sqrt(S0),(len(m.COEF),N)).copy()
    vel=np.broadcast_to(normal[:,1]*math.sqrt(V0),(len(m.COEF),N)).copy()
    quiet=x.copy();qvel=vel.copy();minimum=info['d0']-x;quietmin=info['d0']-quiet
    censored=abs(x)>=cutoff;qcensored=censored.copy();kinetic=np.zeros_like(x)
    def force(z):return -np.sign(z)*curve(np.minimum(abs(z),cutoff),1)/mu
    # Explicit Gaussian reference bath and ordinary optical inflow, common to
    # each driven/quiet path. This is a semiclassical native-cut sensitivity,
    # not an invented native anharmonic quantum stationary state.
    onset=float(d['source_clock'][hot])
    first=max(0,int(np.searchsorted(t,onset)-1))
    for j in range(first,len(t)-1):
        left=max(t[j],onset)
        dt=t[j+1]-left;eps=m.COEF*d['electron'][j,hot];dk=.5*float(d['localK'])*eps
        feed=max(float(d['optical'][j+1,hot]-math.exp(-dt/(s['tau_opt']*1e15))*d['optical'][j,hot]),0.)
        for sub in range(max(1,int(math.ceil(dt)))):
            h=dt/max(1,int(math.ceil(dt)))
            active=~censored;qa=~qcensored
            acceleration=force(x)-dk[:,None]*x/mu
            vx=vel+h*acceleration/2;newx=x+h*vx
            qv=qvel+h*force(quiet)/2;newquiet=quiet+h*qv
            crossed=active&(abs(newx)>=cutoff);qcrossed=qa&(abs(newquiet)>=cutoff)
            kinetic[crossed]=mu*vx[crossed]**2/2
            censored|=crossed;qcensored|=qcrossed
            x=np.where(active,newx,x);quiet=np.where(qa,newquiet,quiet)
            vel=np.where(~censored,vx+h*(force(x)-dk[:,None]*x/mu)/2,vel)
            qvel=np.where(~qcensored,qv+h*force(quiet)/2,qvel)
            at=left+(sub+1)*h-onset;co=math.cos(wr*at);si=math.sin(wr*at)
            bs=S0*co*co+V0*si*si/(wr*wr)
            bv=S0*wr*wr*si*si+V0*co*co
            bc=(V0/wr-S0*wr)*co*si
            decay=math.exp(-h/(s['tau_opt']*1e15));fraction=1-decay
            ns=fraction*bs+feed*h/dt/(3*info['krel'])
            nv=fraction*bv+feed*h/dt/(3*mu);nc=fraction*bc
            noise=rng.standard_normal((2,N));nx=math.sqrt(ns)*noise[0]
            nvel=nc/math.sqrt(ns)*noise[0]+math.sqrt(max(nv-nc*nc/ns,0.))*noise[1]
            x=np.where(~censored,math.sqrt(decay)*x+nx,x);vel=np.where(~censored,math.sqrt(decay)*vel+nvel,vel)
            quiet=np.where(~qcensored,math.sqrt(decay)*quiet+nx,quiet);qvel=np.where(~qcensored,math.sqrt(decay)*qvel+nvel,qvel)
            censored|=abs(x)>=cutoff;qcensored|=abs(quiet)>=cutoff
            minimum=np.minimum(minimum,info['d0']-np.minimum(x,cutoff))
            quietmin=np.minimum(quietmin,info['d0']-np.minimum(quiet,cutoff))
    out=[]
    for b,coef in enumerate(m.COEF):
        for level in info['d0']-xnodes[1:]:
            count=np.count_nonzero(minimum[b]<=level);quiet_count=np.count_nonzero(quietmin[b]<=level)
            # Unknown future paths contribute an explicit additional upper mass.
            unresolved=np.count_nonzero(censored[b]&(minimum[b]>level))
            out.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,coefficient_per_eV=coef,
                position_radius_A=d['position'][hot,0],position_depth_A=d['position'][hot,2],supplied_separation_A=level,
                trajectory_minimum_CDF_lower=count/N,trajectory_minimum_CDF_upper=(count+unresolved)/N,
                quiet_trajectory_minimum_CDF=quiet_count/N,paired_probability_excess=(count-quiet_count)/N,
                effective_weight=1/N,draws=N,observed_count=count,zero_count_IID_reference_not_QMC_interval=1-.05**(1/N) if count==0 else None,
                native_support_exit_fraction=float(censored[b].mean()),quiet_support_exit_fraction=float(qcensored[b].mean()),
                minimum_resolved_quantile_A=float(np.quantile(minimum[b],[.01])[0]),
                deterministic_exit_kinetic_contribution_eV_per_censored_draw=float(kinetic[b,censored[b]].mean()) if np.any(censored[b]) else 0.,
                scope='paired Gaussian-Wigner draws, native even directional cut, shared explicit Gaussian reference bath and ordinary optical feed; support exits censored; semiclassical nonlinear sensitivity'))
    return out


def beam_factors(m,beam):
    xy=np.linspace(-.005,.005,101);X,Y=np.meshgrid(xy,xy,indexing='ij');w=m.p1.beam_weight(X,Y,beam)
    spatial=float(np.mean(w*w)/np.mean(w)**2)
    temporal=(1+.7**2/2+18*.1+81*.1)/1.9**2 if beam=='hot_bursts' else 1.
    norm=w.mean()*m.p1.AREA
    return spatial*temporal,norm


def yield_case(m,d,ec,sp,s,env,coupled,checks):
    t=d['t'];pos=d['position'];weights=d['occupied'];out=[]
    residual=next(r for r in m.p1.DATA['residual'] if float(r['E_cm_keV'])==ec)
    prior=[r for r in m.rows(m.ROOT/'g15_track/pass1/yield_consumers.csv') if float(r['E_cm_keV'])==ec and r['scenario']==s['name']]
    needed=float(prior[0]['G12_required_extra_per_incident_D'])
    # Time-conditioned momentum covariance plus exactly ONE ordinary optical population.
    base=m.M*d['vr']
    delta=m.M*d['delta_v'].astype(float)
    for receiver in (0,1):
        tr=m.p1.track(ec,receiver);ts,z,dep,src=m.p1.source(tr,receiver);z=z*1e10
        energy=np.r_[tr['initial_eV'],[x['energy_eV'] for x in tr['timeline']]]/1000
        lab=np.interp(pos[:,2],z,energy);supported=pos[:,2]<=z[-1]
        for ch in ('p','n'):
            for U,label in ((0.,'bare'),(.04,'static40')):
                # The derivative is a local exact-kernel lookup, not a rate inferred from mean energy.
                grid=np.linspace(float(base.min())*.8,float(base.max())*1.05+1e-5,161)
                step=np.maximum(grid*1e-3,1e-6)
                hi=m.moving_rate(lab[:,None],grid[None,:]+step[None,:],ch,U)
                lo=m.moving_rate(lab[:,None],np.maximum(grid[None,:]-step[None,:],1e-12),ch,U)
                derivative=(hi-lo)/(2*step[None,:])
                inc=np.empty((len(m.COEF),len(pos)))
                for b in range(len(m.COEF)):
                    mid=base+delta[:,b]/2
                    response=np.array([np.interp(mid[:,j],grid,derivative[j])*delta[:,b,j] for j in range(len(pos))]).T
                    inc[b]=integrate(response,t)*supported
                rate=m.old.current(ec)/m.EV*(.85 if sp==0 else .15)/m.p1.AREA
                independent=rate*(inc*weights).sum(axis=1)*1e-15
                b=len(m.COEF)-1;j,k=np.unravel_index(np.argmax(abs(delta[:,b])),delta[:,b].shape)
                exact=float(m.moving_rate(lab[k],base[j,k]+delta[j,b,k],ch,U)-m.moving_rate(lab[k],base[j,k],ch,U))
                approx=float(np.interp(base[j,k]+delta[j,b,k]/2,grid,derivative[k])*delta[j,b,k])
                checks.append(dict(E_cm_keV=ec,species=sp+1,receiver=receiver+1,scenario=s['name'],environment=env,neighbor_on=coupled,
                    branch=ch,comparator=label,exact_extreme_delta_m2=exact,consumer_extreme_delta_m2=approx,
                    relative_difference=abs(approx/exact-1) if exact else 0.))
                for beam in m.p1.BEAMS:
                    factor,_=beam_factors(m,beam)
                    oldrow=next(r for r in prior if int(r['species'])==receiver+1 and r['branch']==ch and r['comparator']==label and r['beam']==beam)
                    for b,coef in enumerate(m.COEF):
                        weight=(.85 if receiver==0 else .30)/1.15
                        out.append(dict(E_cm_keV=ec,source_species=sp+1,receiver_species=receiver+1,scenario=s['name'],environment=env,neighbor_on=coupled,
                            coefficient_per_eV=coef,beam=beam,branch=ch,comparator=label,
                            other_parent_breathing_delta_per_incident_receiver_D=independent[b]*factor,
                            mixture_delta_per_incident_D=weight*independent[b]*factor,
                            mixture_delta_per_electrical_parent=1.15*weight*independent[b]*factor,
                            mixture_raw_delta_per_true_C=1.15*weight*independent[b]*factor/m.EV,
                            fraction_of_G12_residual=weight*independent[b]*factor/needed if needed else None,
                            G12_required_extra_per_incident_D=needed,G12_residual_R=float(residual['R_observed'])-float(residual['R_model']),
                            ordinary_pass1_addition_per_receiver_D=float(oldrow['per_incident_D_beam_moving_addition']),
                            ordinary_pass1_struck_per_D=float(oldrow['per_incident_D_struck_D']),
                            profile_exposure_s=43200 if ec<=2 else 7200,native_sheet_once=True,pumped_sheet=False,
                            kernel='actual local CSDA speed; time-conditioned Gaussian velocity covariance; matched same-state quiet plus ordinary optical population once'))
    return out


def struck_case(m,d,ec,sp,s,env,coupled):
    tr=m.p1.track(ec,sp);ts,z,dep,_=m.p1.source(tr,sp);z*=1e10;ts*=1e15
    pos=d['position'];t=d['t'];samples=[]
    gl,gw=np.polynomial.legendre.leggauss(8)
    for j,row in enumerate(dep):
        projectile=tr['initial_eV'] if j==0 else tr['timeline'][j-1]['energy_eV']
        zz=(z[j]+z[j+1])/2;nh,nd=next(layer[2:4] for layer in m.p1.PROFILES[ec]['layers'] if layer[1]*1e10>=zz)
        for frac in (.05,.2,.5,.9):
            energy=projectile*frac;count=row[2]/(4*energy)
            stopping=(nh*(m.p1.bca.se_ls(energy,46)+m.p1.bca.sn_zbl(energy,46,106.42))+nd*(m.p1.bca.se_ls(energy,1)+m.p1.bca.sn_zbl(energy,1,m.p1.bca.M_D)))*1e-18
            length=energy/stopping*1e10;distance=(gl+1)*length/2;remaining=energy*(1-distance/length)
            speed=np.sqrt(2*energy/m.M);age=2*length/speed*(1-np.sqrt(1-distance/length))
            cosine=math.sqrt(frac);rad=distance*math.sqrt(1-frac);depth=zz+distance*cosine
            sites=np.argmin((pos[:,0,None]-rad)**2+(pos[:,2,None]-depth)**2,axis=0)
            at=(ts[j]+ts[j+1])/2+age
            keep=(rad<=s['radius']*1e10)&(depth<=pos[:,2].max()*16/15)
            for k in range(len(gl)):
                if keep[k]:samples.append((sites[k],at[k],remaining[k]/1000,count*nd*length*1e-10*gw[k]/2))
    out=[]
    if not samples:return out
    sites,at,lab,weight=np.array(samples).T;sites=sites.astype(int)
    quiet=np.array([np.interp(tm,t,m.M*d['vr'][:,site]) for site,tm in zip(sites,at)])
    for ch in ('p','n'):
        for U,label in ((0.,'bare'),(.04,'static40')):
            for b,coef in enumerate(m.COEF):
                change=np.array([np.interp(tm,t,m.M*d['delta_v'][:,b,site]) for site,tm in zip(sites,at)])
                delta=stable_rate_delta(m,lab,quiet,change,ch,U)
                value=float(weight@delta)/(sp+1)
                for beam in m.p1.BEAMS:
                    out.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,neighbor_on=coupled,
                        coefficient_per_eV=coef,beam=beam,branch=ch,comparator=label,
                        struck_D_breathing_delta_per_incident_D=value,per_electrical_parent=value*(sp+1),
                        source_count_weight=.85 if sp==0 else .15,
                        scope='same-parent actual four-transfer recoil rays, local slowing and arrival time; supported source footprint only',
                        supported_flight_quadrature_points=len(samples)))
    return out


def acoustic_case(m,d,ec,sp,s,env,coupled):
    t=d['whole_receiver']['t'];dt=np.diff(t);w=d['occupied'];N=w.sum();V=d['formula'].sum()*m.VOLUME*1e-30
    # g is per LINEAR hydrostatic strain. Pressure uses /3 of the volume derivative.
    g=m.P['collective_strain_slope_eV_A2'];B=m.p1.BULK;rho=m.p1.RHO
    cL=math.sqrt(227.1e9/rho);length=max(s['radius'],d['position'][:,2].max()*1e-10/4)
    wa=cL/length*1e-15;tau=length/cL*1e15
    output=[];transfer={}
    generator=np.array([[0.,1.],[-wa*wa,-wa]])
    for b,coef in enumerate(m.COEF):
        deltaq=d['whole_receiver']['average_deltaq'][:,b]
        sigma=g*m.EV*N/(2*V)*deltaq  # 3 axes /6, correct linear-hydro convention
        equilibrium=-sigma/B
        strain=0.;velocity=0.;work=0.;radiated=0.;peak=0.;peak_energy=0.;feedback=0.;energy=0.
        for j,h in enumerate(dt):
            target=(equilibrium[j]+equilibrium[j+1])/2
            if h not in transfer:transfer[h]=m.expm(generator*h)
            state=transfer[h]@np.array([strain-target,velocity])+np.array([target,0.])
            new,newv=state
            stepwork=-V*(sigma[j]+sigma[j+1])/2*(new-strain)/m.EV
            newenergy=B*V*(new*new+(newv/wa)**2)/(2*m.EV)
            # Exact constant-source work minus stored energy is outgoing elastic
            # radiation through the R/c port, NOT a new unnamed heat bin.
            radiated+=stepwork-(newenergy-energy)
            work+=stepwork;strain,velocity=new,newv;energy=newenergy
            peak=max(peak,abs(strain));peak_energy=max(peak_energy,energy)
            feedback=max(feedback,abs(g*strain/3))
        available=float(np.max(np.maximum(d['work'][:,b],0.)))
        demand=peak_energy+max(radiated,0.)
        scale=math.sqrt(available/demand) if demand>available+1e-8 else 1.
        old_source=float(next(row for row in m.p1.DATA['coherent'] if row['scenario']=='track' and row['species']==('D+' if sp==0 else 'D2+'))['acoustic_J'])/m.EV
        direct=abs(coef)*d['electron'].max()*float(d['localK'])
        output.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,neighbor_on=coupled,coefficient_per_eV=coef,
            maximum_signed_pressure_Pa=float((-sigma).max()),minimum_signed_pressure_Pa=float((-sigma).min()),
            quasistatic_volume_strain_peak=float(abs(equilibrium).max()),conditional_propagating_volume_strain_peak=peak*scale,
            elastic_transit_fs=tau,conditional_acoustic_energy_peak_eV=peak_energy*scale*scale,
            acoustic_stored_final_eV=energy*scale*scale,outgoing_elastic_radiation_eV=radiated*scale*scale,
            acoustic_work_integral_eV=work*scale*scale,energy_cap_scale=scale,
            acoustic_energy_closure_eV=(energy+radiated-work)*scale*scale,
            characteristic_frequency_Hz=wa*1e15/(2*np.pi),
            stiffness_feedback_over_direct_peak=feedback*scale/direct if direct else 0.,
            existing_G13_thermal_recoil_launch_eV=old_source,existing_launch_added_again=False,
            clamped_ion_electronic_pressure='not determined by K(Ue); same-state electronic free-energy strain derivative missing',
            conservation='conditional reciprocal elastic work is reallocated from the driven covariance, not a second electronic input; resolved Born backreaction shown separately',
            optical_energy_reallocation_at_end_eV=-(energy+radiated)*scale*scale))
    return output


def whole_rates(m,summaries):
    out=[];marks=[];rng=np.random.default_rng(1533)
    for ec in m.EC:
        for s in m.p1.SCENARIOS:
            for beam in m.p1.BEAMS:
                factor,norm=beam_factors(m,beam);rate=m.old.current(ec)/m.EV
                for receiver,xy in (('center',(0.,0.)),('hotspot',(.0025,-.0025))):
                    density=float(m.p1.beam_weight(np.array(xy[0]),np.array(xy[1]),beam))/norm
                    nearby=rate*density*np.pi*(2*s['radius'])**2
                    span=8*max(s['tau_opt'],s['tau_ne'],s['tau_end'])+1e-12
                    window=512/max(nearby,1e-30);ceiling=10.7/1.9 if beam=='hot_bursts' else 1.
                    count=rng.poisson(nearby*window*ceiling);at=np.sort(rng.uniform(0,window,count))
                    if beam=='hot_bursts':
                        cuts,states=m.p1.cox_times(window,rng);phase=rng.uniform(0,2*np.pi)
                        inten=(1+.7*np.cos(2*np.pi*1e5*at+phase)+9*states[np.searchsorted(cuts,at,side='right')-1])/1.9
                        at=at[rng.random(len(at))<inten/ceiling]
                    species=(rng.random(len(at))>=.85).astype(int)+1
                    prior=np.r_[False,np.diff(at)<2*span]
                    for i,tm in enumerate(at):marks.append(dict(E_cm_keV=ec,scenario=s['name'],beam=beam,receiver=receiver,parent=i,time_s=tm,species=int(species[i]),overlapping_prior=bool(prior[i]),window_s=window))
                    overlap=-math.expm1(-nearby*2*span*ceiling)
                    for env in m.PAIR:
                        for coupled in (False,True):
                            for coef in m.COEF:
                                pair=[next(r for r in summaries if float(r['E_cm_keV'])==ec and int(r['species'])==sp and r['scenario']==s['name'] and r['environment']==env and str(r['neighbor_on'])==str(coupled) and float(r['coefficient_per_eV'])==coef) for sp in (1,2)]
                                peaks=np.array([float(r['maximum_site_excess_eV']) for r in pair])
                                for threshold in (.0001,.001,.01,.1):
                                    probability=float(np.array([.85,.15])@(peaks>=threshold))
                                    out.append(dict(E_cm_keV=ec,scenario=s['name'],environment=env,neighbor_on=coupled,coefficient_per_eV=coef,
                                        beam=beam,receiver=receiver,threshold_eV_D=threshold,one_physical_parent_patch_probability=probability,
                                        whole_foil_patch_rate_s=rate*probability,physical_electrical_parents_s=rate,
                                        source_D_plus_max_eV=peaks[0],source_D2_plus_max_eV=peaks[1],
                                        local_compatible_parent_rate_s=nearby,source_support_window_s=span,
                                        analytic_local_overlap_upper=overlap,actual_window_parent_marks=len(at),actual_window_overlap_count=int(prior.sum()),
                                        marked_window_s=window,foil_exposure_s=43200 if ec<=2 else 7200,
                                        full_exposure_patch_count=rate*probability*(43200 if ec<=2 else 7200),
                                        actual_window_overlap_fraction=float(prior.mean()) if len(prior) else 0.,
                                        zero_overlap_IID_reference_not_Cox_bound=1-.05**(1/len(at)) if len(at) and not prior.any() else None,
                                        Cox_Palm_factor=factor,environment_abundance='bulk=1 or conditional cage footprint; actual cage abundance not supplied',
                                        counting='one maximum per electrical parent, two true siblings; local support/window overlap, not foil current times lifetime'))
    m.table('parent_windows.csv',marks);return out


def quiet_tails(m):
    output=[]
    exact=np.load(m.ROOT/'g10_quantum_site/output/exact/PdD_300K_joint.npz')
    coarse=np.load(m.ROOT/'g10_quantum_site/output/exact/PdD_300K_basis8.npz')
    radial=exact['radial_A'];density=exact['radial_density_A_minus3']
    for radius in (.24,.4,.6,.8):
        j=int(np.argmin(abs(radial-radius)))
        output.append(dict(radius_A=float(radial[j]),native_radial_survival=float(exact['radial_tail'][j]),
            native_density_Aminus3=float(density[j]),basis8_density_Aminus3=float(coarse['radial_density_A_minus3'][j]),
            basis_relative_density_change=float(coarse['radial_density_A_minus3'][j]/density[j]-1),
            Gaussian_radial_survival=float(1-math.erf(radial[j]/math.sqrt(2*m.Q0))+math.sqrt(2/math.pi)*radial[j]/math.sqrt(m.Q0)*math.exp(-radial[j]**2/(2*m.Q0))),
            scope='single-site radial tail, not a pair distribution; Gaussian pair closure does not inherit the native non-Gaussian tail'))
    return output


def load_case(m,name):
    with np.load(m.HERE/name) as opened:d={k:opened[k] for k in opened.files}
    for field,reference in (('q','qr'),('v','vr'),('c','cr')):
        d[field]=d[reference][:,None,:]+d['delta_'+field].astype(float)
    for field in ('relq','relv','relc'):
        d[field]=d['reference_'+field][:,None,:]+d['delta_'+field].astype(float)
    return d

def scientific_modules(m):
    import importlib.util
    spec=importlib.util.spec_from_file_location('pass3_joined',m.HERE/'joined.py')
    join=importlib.util.module_from_spec(spec);spec.loader.exec_module(join)
    spec=importlib.util.spec_from_file_location('pass3_receivers',m.HERE/'receivers.py')
    receiver=importlib.util.module_from_spec(spec);spec.loader.exec_module(receiver)
    return join,receiver

def phase_stats(phases):
    raw=[r['phase_aligned_second_moment_fraction'] for r in phases if r['coefficient_per_eV']!=0. and r['phase_defined']]
    resolved=[r['phase_aligned_second_moment_fraction'] for r in phases if r['coefficient_per_eV']!=0. and r['phase_support']=='positive_supported']
    return dict(maximum_raw_defined_phase_fraction=max(raw,default=None),
        maximum_resolved_matched_phase_fraction=max(resolved,default=None),
        maximum_phase_triangle_excess=max(r['phase_triangle_excess'] for r in phases))

def run(m):
    summaries=m.rows(m.HERE/'phase_summary.csv');phases=[];pairs=[];minima=[];yields=[];struck=[];acoustic=[];checks=[]
    join,receiver=scientific_modules(m)
    siblings=[];neighbors=[];carriers=[];strengths=[]
    for ec in m.EC:
        for sp in (0,1):
            for s in m.p1.SCENARIOS:
                for env in m.PAIR:
                    for coupled in (False,True):
                        name=f'cov_{ec:g}_{sp+1}_{s["name"]}_{env}_{int(coupled)}.npz'
                        d=load_case(m,name)
                        phases+=field_pairs(m,d,ec,sp,s,env,coupled,receiver)
                        pairs+=instantaneous_pairs(m,d,ec,sp,s,env,coupled)
                        if coupled:minima+=trajectory_minima(m,d,ec,sp,s,env)
                        yields+=yield_case(m,d,ec,sp,s,env,coupled,checks)
                        struck+=struck_case(m,d,ec,sp,s,env,coupled)
                        acoustic+=acoustic_case(m,d,ec,sp,s,env,coupled)
                        neighbors+=join.neighbor(m,d,ec,sp,s,env,coupled)
                        carriers+=join.carrier(m,d,ec,sp,s,env,coupled)
                        strengths+=join.strength_and_support(m,d,ec,sp,s,env,coupled)
                        if sp==1:siblings+=join.sibling(m,ec,s,env,coupled)
    m.table('phase_evolution.csv',phases);m.table('pair_instantaneous.csv',pairs);m.table('pair_trajectory_minima.csv',minima)
    m.table('yield_consumers.csv',yields);m.table('struck_consumers.csv',struck);m.table('acoustic_join.csv',acoustic)
    m.table('whole_foil.csv',whole_rates(m,summaries));m.table('native_quiet_tails.csv',quiet_tails(m));m.dump('consumer_checks.json',checks)
    mixture=join.mixture(m,yields,struck,siblings,neighbors)
    m.table('sibling_consumers.csv',siblings);m.table('neighbor_consumers.csv',neighbors)
    m.table('carrier_delta.csv',carriers);m.table('yield_mixture.csv',mixture)
    transmitted=join.acoustic_transmission(m,acoustic)
    m.table('acoustic_transmission.csv',transmitted);m.table('cage_axial_transfer.csv',join.cage_transfer(m))
    m.table('strength_and_support.csv',strengths)
    import importlib.util
    spec=importlib.util.spec_from_file_location('pass3_four_mode',m.HERE/'four_mode.py')
    band=importlib.util.module_from_spec(spec);spec.loader.exec_module(band)
    retained=band.run(m)
    quiet=[r for r in summaries if float(r['coefficient_per_eV'])==0]
    return dict(base=m.BASE,source_cases=len(summaries),coefficients=m.COEF.tolist(),
        maximum_conservation_error_eV=max(float(r['conservation_error_eV']) for r in summaries),
        maximum_occupied_column_relative_error=max(abs(float(r['occupied_column_relative_error'])) for r in summaries),
        maximum_retarded_whole_energy_bound_eV=max(float(r['retarded_whole_conservation_bound_eV']) for r in summaries),
        maximum_retarded_quadrature_count_error_D=max(abs(r['retarded_quadrature_count_error_D']) for r in phases),
        **phase_stats(phases),
        minimum_actual_carrier_eV_formula=min(float(r['min_actual_carrier_eV_formula']) for r in summaries),
        maximum_uncertainty_violation=max(float(r['uncertainty_violation_eV2_fs2']) for r in summaries),
        maximum_zero_drive_variance_error=max(abs(float(r['maximum_variance_ratio'])-1) for r in quiet),
        maximum_zero_drive_work_eV=max(abs(float(r['signed_work_eV'])) for r in quiet),
        maximum_variance_ratio=max(float(r['maximum_variance_ratio']) for r in summaries),
        minimum_variance_ratio=min(float(r['minimum_variance_ratio']) for r in summaries),
        maximum_local_excess_eV=max(float(r['maximum_site_excess_eV']) for r in summaries),
        maximum_abs_other_parent_residual_fraction=max(abs(r['fraction_of_G12_residual'] or 0) for r in yields),
        maximum_abs_joined_residual_fraction=max(abs(r['fraction_of_G12_residual'] or 0) for r in mixture),
        maximum_abs_sibling_delta_per_D=max(abs(r['mixture_delta_per_incident_D']) for r in siblings),
        maximum_acoustic_energy_closure_eV=max(abs(r['acoustic_energy_closure_eV']) for r in acoustic),
        maximum_acoustic_transmission_partition_error_eV=max(abs(r['outgoing_partition_closure_eV']) for r in transmitted),
        maximum_abs_struck_delta_per_D=max(abs(r['struck_D_breathing_delta_per_incident_D']) for r in struck),
        maximum_acoustic_volume_strain=max(r['conditional_propagating_volume_strain_peak'] for r in acoustic),
        maximum_acoustic_feedback_ratio=max(r['stiffness_feedback_over_direct_peak'] for r in acoustic),
        maximum_kernel_derivative_relative_error=max(r['relative_difference'] for r in checks),
        retained_global_modes=retained,
        G16='not supplied at source freeze; six signed comparison branches retained; exact q/occupation/temperature/internal-energy join required',
        no_nuclear_contact_or_cage_abundance_inferred=True)


def replay_whole(m):
    """Rebuild only derived causal phase/elastic/event readouts, never sources or pair paths."""
    import json
    join,receiver=scientific_modules(m);phases=[];acoustic=[]
    for ec in m.EC:
        for sp in (0,1):
            for s in m.p1.SCENARIOS:
                for env in m.PAIR:
                    for coupled in (False,True):
                        d=load_case(m,f'cov_{ec:g}_{sp+1}_{s["name"]}_{env}_{int(coupled)}.npz')
                        phases+=field_pairs(m,d,ec,sp,s,env,coupled,receiver)
                        acoustic+=acoustic_case(m,d,ec,sp,s,env,coupled)
    transmitted=join.acoustic_transmission(m,acoustic)
    m.table('phase_evolution.csv',phases);m.table('acoustic_join.csv',acoustic);m.table('acoustic_transmission.csv',transmitted)
    whole_rates(m,m.rows(m.HERE/'phase_summary.csv'))
    result=json.loads((m.HERE/'results.json').read_text())
    result.pop('maximum_matched_phase_fraction',None);result.update(phase_stats(phases))
    result.update(maximum_acoustic_volume_strain=max(r['conditional_propagating_volume_strain_peak'] for r in acoustic),
        maximum_acoustic_feedback_ratio=max(r['stiffness_feedback_over_direct_peak'] for r in acoustic),
        maximum_acoustic_energy_closure_eV=max(abs(r['acoustic_energy_closure_eV']) for r in acoustic),
        maximum_acoustic_transmission_partition_error_eV=max(abs(r['outgoing_partition_closure_eV']) for r in transmitted))
    m.dump('results.json',result)
    return result
