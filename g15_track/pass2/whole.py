"""Full-footprint G3 modal Green, marked Cox replay and causal flight consumers.
The local nonlinear microblocks are not used as long-distance lattice bonds.
"""
import math
import numpy as np
from scipy.special import roots_legendre


def finite_motion(age,F,impulse,M,omega,kappa,beta,rho,G=None):
    """Exact finite-time damped Green for the signed rise-and-relaxation pulse."""
    age=np.asarray(age,float);s=np.maximum(age,0.)
    alpha=kappa/2;wd=math.sqrt(max(omega*omega-alpha*alpha,1e-16))
    decay=np.exp(-alpha*s);co=np.cos(wd*s);si=np.sin(wd*s)
    def exponential(b):
        den=(b-alpha)**2+wd*wd
        q=(np.exp(-b*s)-decay*co+(b-alpha)/wd*decay*si)/den
        v=(-b*np.exp(-b*s)+decay*(alpha*co+wd*si)+(b-alpha)/wd*decay*(wd*co-alpha*si))/den
        return q,v
    q1,v1=exponential(beta);q2,v2=exponential(beta+rho)
    q=F*(q1-q2)[...,None]/M+impulse*(decay*si/(M*wd))[...,None]
    v=F*(v1-v2)[...,None]/M+impulse*(decay*(co-alpha/wd*si)/M)[...,None]
    if G is not None:
        pulse=np.exp(-beta*s)-np.exp(-(beta+rho)*s)
        q+=G*(v1-v2)[...,None]/M
        v+=G*(pulse-kappa*(v1-v2)-omega*omega*(q1-q2))[...,None]/M
    mask=(age>=0)[...,None]
    return q*mask,v*mask


def dispersion(m,q):
    q=np.asarray(q,float);j=0.;grad=np.zeros(3)
    for shell in m.P['shells']:
        delta=np.array(shell['periodic_alias_offsets_A']);a=delta@q
        value=shell['J_sum_eV_A2']/len(delta)
        j+=value*np.cos(a).sum();grad+=value*(delta*np.sin(a)[:,None]).sum(0)
    omega=math.sqrt(max((m.K-j)/m.M,1e-12))
    return omega,grad/(2*m.M*omega)


def mode_spec(m,ec,sp,closure):
    tr=m.p1.track(ec,sp);v0=math.sqrt(2*tr['initial_eV']/m.M);radius=closure['radius']*1e10
    stripe=2*np.pi/5.3/np.sqrt(2)*np.array([1.,1.,0.])
    names=['Gamma','retarded_clock','opposite_sheet_110','opposite_sheet_1m10','radial']
    waves=np.array([np.zeros(3),[0.,0.,m.OMEGA/v0],stripe,stripe*np.array([1,-1,1]),[np.pi/radius,0.,0.]])
    return names,waves,np.array([dispersion(m,q)[0] for q in waves])

def global_modes(m,summary):
    """Retain motif-sampled probes without promoting their spatial quadrature.
    The resolved full-window electronic projection is written by linear_join.
    """
    rows=[];bounds=[];local=[]
    maxJ=sum(abs(s['J_sum_eV_A2']) for s in m.P['shells'])
    omin=math.sqrt(max((m.K-maxJ)/m.M,1e-12))
    speed_bound=sum(abs(s['J_sum_eV_A2'])*max(np.linalg.norm(x) for x in s['periodic_alias_offsets_A'])
        for s in m.P['shells'])/(2*m.M*omin)
    for ec in m.EC:
        for sp in (0,1):
            for closure in m.p1.SCENARIOS:
                kappa=1/(closure['tau_opt']*1e15)
                bounds.append(dict(E_cm_keV=ec,species=sp+1,scenario=closure['name'],
                    local_period_fs=2*np.pi/m.OMEGA,Gamma_period_fs=2*np.pi/dispersion(m,np.zeros(3))[0],
                    interblock_group_velocity_upper_A_fs=speed_bound,
                    one_energy_lifetime_reach_nm=speed_bound/kappa/10,footprint_radius_nm=closure['radius']*1e9,
                    finite_source_depth_nm=m.p1.track(ec,sp)['timeline'][-1]['depth_A']/10,
                    operator_phase_bound_over_lifetime=maxJ/(2*m.M*omin*kappa),
                    interpretation='linear bulk coupling is included in the resolved full-window Green; this reach alone does not exclude nonlinear interblock feedback',
                    boundary='physical source window, not a million-site nonlinear simulation; cage cross-transport is unmeasured'))
                for env in m.D0:
                    with np.load(m.HERE/f'peak_{ec:g}_{sp+1}_{closure["name"]}_{env}.npz',allow_pickle=False) as archive:
                        keys=['occupied_weight','mode_omega_fs_inverse','mode_q_A_inverse','mode_time_fs','mode_names',
                            'actual_Q_A','actual_V_A_fs','green_Q_A','green_V_A_fs','green_work_eV','green_loss_eV','mode_h',
                            'harmonic_reference_eV','mode_h_trend_leakage','mode_h_window_available','mode_h_trend_robust',
                            'local_modal_amplitude_A','local_modal_fit_valid','local_mode_omega_fs_inverse','local_mode_info',
                            'local_block_clock_phase','corrected_carriers_eV','pair_weights']
                        d={k:archive[k] for k in keys}
                    N=float(d['occupied_weight'].sum());om=d['mode_omega_fs_inverse'];waves=d['mode_q_A_inverse']
                    for b,(label,chi) in enumerate(m.BRANCHES):
                        for control in (0,1):
                            lane=b*3+control
                            for j,at in enumerate(d['mode_time_fs']):
                                amp=d['local_modal_amplitude_A'][j,lane];w=d['pair_weights'];valid=d['local_modal_fit_valid'][j,lane]
                                raw_energy=float(.5*m.M*float(d['local_mode_omega_fs_inverse'])**2*np.sum(w*abs(amp)**2))
                                ww=w*valid;energy=float(.5*m.M*float(d['local_mode_omega_fs_inverse'])**2*np.sum(ww*abs(amp)**2)) if ww.sum() else None
                                denominator=float(ww.sum()*np.sum(ww*abs(amp)**2))
                                fraction=float(abs(np.sum(ww*amp*d['local_block_clock_phase']))**2/denominator) if denominator else None
                                local.append(dict(E_cm_keV=ec,species=sp+1,scenario=closure['name'],environment=env,branch=label,
                                    control='built' if control==0 else 'matched_random',time_fs=at,
                                    frequency_THz=float(d['local_mode_omega_fs_inverse'])*1000/(2*np.pi),
                                    resolved_mode_weight_fraction=float(np.average(valid,weights=w)),
                                    raw_local_fit_reference_eV=raw_energy,resolved_local_closing_reference_eV=energy,clock_matched_fraction=fraction,
                                    projected_normal_mode_number=energy/(m.p1.HB*1e15*float(d['local_mode_omega_fs_inverse'])) if energy is not None else None,
                                    ordinary_population_eV=float(d['corrected_carriers_eV'][j,lane,4]),
                                    ordinary_mean_number_per_polarization_mode=float(d['corrected_carriers_eV'][j,lane,4])/
                                        (3*N*m.p1.HB*1e15*float(d['local_mode_omega_fs_inverse'])),
                                    meaning='initial local Hessian mode, actual fitted trajectory; unresolved/trending intervals are not zero optical response'))
                                for k,name in enumerate(d['mode_names']):
                                    omega=om[k];Q=d['actual_Q_A'][j,lane,k];V=d['actual_V_A_fs'][j,lane,k]
                                    Qg=d['green_Q_A'][j,lane,k];Vg=d['green_V_A_fs'][j,lane,k]
                                    energy=float(.5*m.M*N*(np.sum(abs(V)**2)+omega**2*np.sum(abs(Q)**2)))
                                    ref=float(d['harmonic_reference_eV'][j,lane,k]);h=float(d['mode_h'][j,lane,k])
                                    available=bool(d['mode_h_window_available'][j,lane,k]);robust=bool(d['mode_h_trend_robust'][j,lane,k])
                                    green=float(.5*m.M*N*(np.sum(abs(Vg)**2)+omega**2*np.sum(abs(Qg)**2)))
                                    work=float(d['green_work_eV'][j,lane,k]);loss=float(d['green_loss_eV'][j,lane,k])
                                    rows.append(dict(E_cm_keV=ec,species=sp+1,scenario=closure['name'],environment=env,branch=label,
                                        control='built' if control==0 else 'matched_random',mode=str(name),time_fs=at,
                                        qx_A_inverse=waves[k,0],qy_A_inverse=waves[k,1],qz_A_inverse=waves[k,2],
                                        optical_frequency_THz=omega*1000/(2*np.pi),single_complex_plane_reference_eV=energy,
                                        harmonic_coordinate_reference_eV=ref,canonical_projection_fraction=energy/ref if ref>1e-100 else 0.,
                                        projected_normal_mode_number=energy/(m.p1.HB*1e15*omega),
                                        motif_sampled_linear_receiver_eV=green,receiver_signed_work_eV=work,receiver_loss_eV=loss,
                                        receiver_energy_defect_eV=green+loss-work,h_window_available=available,h_trend_robust=robust,
                                        h_2omega_fit=h if available else None,
                                        h_trend_spread=float(d['mode_h_trend_leakage'][j,lane,k]) if available else None,
                                        h_over_resonant_threshold=h/(2*kappa/omega) if robust else None,
                                        meaning='coarse motif-only quadrature; not a resolved full-window stripe/track overlap, not a sum of independent modes'))
    m.table('motif_mode_probes.csv',rows);m.table('local_modes.csv',local);m.table('phase_transport.csv',bounds)
    return rows,bounds


def peak_map(m, ec, sp, scenario, env):
    with np.load(m.HERE/f'peak_{ec:g}_{sp+1}_{scenario["name"]}_{env}.npz',allow_pickle=False) as opened:
        d={k:opened[k] for k in ('occupied_weight','position_A','site_peak_total_reference_eV',
                                'escape_marks','mechanical_candidate_transfer_eV','mechanical_strong_mark')}
        d['span_s']=float(opened['mode_time_fs'][-1])*1e-15
    peak=d.pop('site_peak_total_reference_eV').reshape(3*len(m.BRANCHES),-1).copy()
    marks=d.pop('escape_marks')
    lane,block,site=marks[:,1:4].astype(int).T
    i=block*int(np.prod(m.SHAPE))+site
    exit_reference=.5*(m.K*np.sum(marks[:,4:7]**2,axis=1)+m.M*np.sum(marks[:,7:10]**2,axis=1))
    np.maximum.at(peak,(lane,i),exit_reference)
    d['peak']=peak
    return d


def patch_level_rates(m):
    """Count one seeded maximum per parent, not one count per occupied target."""
    rows=[];xx=np.linspace(-.005,.005,101);X,Y=np.meshgrid(xx,xx,indexing='ij')
    density_max={beam:(1. if beam=='even_steady' else 9.)/
                 (m.p1.AREA*m.p1.beam_weight(X,Y,beam).mean()) for beam in m.p1.BEAMS}
    for ec in m.EC:
        parents=m.current(ec)/m.EV
        for s in m.p1.SCENARIOS:
            area=np.pi*(2*s['radius'])**2
            for env in m.D0:
                data=[peak_map(m,ec,sp,s,env) for sp in (0,1)]
                maxima=np.array([d['peak'][:,d['occupied_weight']>0].max(axis=1) for d in data])
                span=max(d['span_s'] for d in data)
                for b,(branch,_) in enumerate(m.BRANCHES):
                    for ctrl,control in enumerate(('built','matched_random')):
                        lane=3*b+ctrl
                        for threshold in (.01,.1,1.):
                            fraction=float(np.array([.85,.15])@(maxima[:,lane]>=threshold))
                            for beam in m.p1.BEAMS:
                                ceiling=10.7/1.9 if beam=='hot_bursts' else 1.
                                # Within this represented horizon, starts in
                                # [-span,+span] can meet the original patch.
                                overlap=-math.expm1(-parents*density_max[beam]*ceiling*area*2*span)
                                rows.append(dict(E_cm_keV=ec,scenario=s['name'],environment=env,branch=branch,
                                    control=control,beam=beam,threshold_eV_per_D=threshold,
                                    D_plus_max_eV=float(maxima[0,lane]),D2_plus_max_eV=float(maxima[1,lane]),
                                    physical_electrical_parents_s=parents,isolated_modeled_fraction=fraction,
                                    isolated_modeled_patch_rate_s=parents*fraction,
                                    finite_horizon_s=span,patch_overlap_probability_upper=overlap,
                                    patch_rate_lower_s=parents*fraction*(1-overlap),
                                    patch_rate_upper_s=parents*(fraction+(1-fraction)*overlap),
                                    basis='one maximum over the represented occupied sites/times per parent; two-sided Cox overlap bound; fixed registration and unit environment prevalence'))
    m.table('patch_level_rates.csv',rows)
    return rows


def marked_windows(m):
    """Shared physical-parent Cox marks; per-site peaks include chart-exit marks.
    Overlap is bounded, never synthesized by summing independently capped states.
    """
    cache={};rows=[];marks=[];rng=np.random.default_rng(152711)
    for ec in m.EC:
        for s in m.p1.SCENARIOS:
            disk=np.pi*s['radius']**2;kernel=8*s['tau_opt']+1e-12
            for sp in (0,1):
                for env in m.D0:
                    d=peak_map(m,ec,sp,s,env)
                    d['cdf']=np.cumsum(d['occupied_weight'])/d['occupied_weight'].sum()
                    cache[sp,env]=d
            for beam in m.p1.BEAMS:
                xx=np.linspace(-.005,.005,101);X,Y=np.meshgrid(xx,xx,indexing='ij')
                norm=m.p1.beam_weight(X,Y,beam).mean()*m.p1.AREA
                for receiver,xy in (('center',(0.,0.)),('hotspot',(.0025,-.0025))):
                    density=float(m.p1.beam_weight(np.array(xy[0]),np.array(xy[1]),beam))/norm
                    lam=m.current(ec)/m.EV*disk*density;window=256/(16*lam);events=[]
                    ceiling=10.7/1.9 if beam=='hot_bursts' else 1.
                    for record in range(16):
                        duration=window+kernel
                        cuts,states=m.p1.cox_times(duration,rng) if beam=='hot_bursts' else (None,None)
                        n=rng.poisson(lam*duration*ceiling);t=np.sort(rng.uniform(0,duration,n))
                        phase=rng.uniform(0,2*np.pi)
                        if beam=='hot_bursts':
                            intensity=(1+.7*np.cos(2*np.pi*1e5*t+phase)+9*states[np.searchsorted(cuts,t,side='right')-1])/1.9
                            t=t[rng.random(len(t))<intensity/ceiling]
                        sp=(rng.random(len(t))>=.85).astype(int);u=rng.random(len(t))
                        for j,at in enumerate(t):
                            if at<kernel:continue
                            prior=np.flatnonzero((t<at)&(at-t<kernel))
                            events.append((record,j,at-kernel,int(sp[j]),u[j],sp[prior],u[prior]))
                    overlap_upper=-math.expm1(-lam*kernel*ceiling)
                    for env in m.D0:
                        for b,(branch,chi) in enumerate(m.BRANCHES):
                            for control in (0,1):
                                lane=3*b+control;counts=np.zeros(3);possible=np.zeros(3);shared=np.zeros(3)
                                largest=0.;strong=0;levels=np.array([.01,.1,1.])
                                for record,j,at,sp,u,prior_sp,prior_u in events:
                                    d=cache[sp,env];site=min(np.searchsorted(d['cdf'],u),len(d['cdf'])-1)
                                    own=float(d['peak'][lane,site]);depth=d['position_A'][site,2];amplitudes=0.
                                    for ps,pu in zip(prior_sp,prior_u):
                                        pd=cache[int(ps),env]
                                        if depth>m.p1.track(ec,int(ps))['timeline'][-1]['depth_A']:continue
                                        # A full sampled-footprint peak is conservative for
                                        # unresolved prior ray/depth registration; not a phase state.
                                        amplitudes+=math.sqrt(float(pd['peak'][lane,pd['occupied_weight']>0].max()))
                                    lower=max(math.sqrt(own)-amplitudes,0.)**2
                                    upper=(math.sqrt(own)+amplitudes)**2
                                    counts+=lower>=levels;possible+=upper>=levels;shared+=(upper>=levels)*bool(len(prior_sp))
                                    largest=max(largest,upper);strong+=bool(d['mechanical_strong_mark'][site])
                                    marks.append(dict(E_cm_keV=ec,scenario=s['name'],beam=beam,receiver=receiver,record=record,parent=j,
                                        time_s=at,species=sp+1,environment=env,branch=branch,control=control,occupied_site_quantile=u,
                                        depth_A=depth,isolated_peak_reference_eV=own,overlap_peak_reference_lower_eV=lower,
                                        overlap_peak_reference_upper_eV=upper,prior_independent_parents=len(prior_sp),
                                        unresolved_mechanical_candidate_eV=float(d['mechanical_candidate_transfer_eV'][site]),
                                        energetic_mechanical_mark=bool(d['mechanical_strong_mark'][site])))
                                n=len(events)
                                for k,level in enumerate(levels):
                                    lo=counts[k]/n if n else None;hi=possible[k]/n if n else None
                                    rows.append(dict(E_cm_keV=ec,scenario=s['name'],beam=beam,receiver=receiver,environment=env,branch=branch,
                                        control='built' if control==0 else 'matched_random',level_eV_D=level,physical_parent_marks=n,
                                        independent_receiver_records=16,total_receiver_seconds=16*window,total_receiver_hours=16*window/3600,
                                        foil_exposure_hours=12 if ec<=2 else 2,source_kernel_s=kernel,receiver_parent_rate_s=lam,
                                        sampled_peak_fraction_lower=lo,sampled_peak_fraction_upper=hi,
                                        conditional_whole_beam_level_rate_lower_s=m.current(ec)/m.EV*lo if lo is not None else None,
                                        conditional_whole_beam_level_rate_upper_s=m.current(ec)/m.EV*hi if hi is not None else None,
                                        observed_shared_peak_overlap_fraction=shared[k]/n if n else None,
                                        analytic_Cox_Palm_overlap_upper=overlap_upper,
                                        shared_probability_Frechet_lower=0.,shared_probability_Frechet_upper=min(hi,overlap_upper) if hi is not None else None,
                                        zero_success_geometry_95pct_upper=1-.05**(1/n) if n and possible[k]==0 else None,
                                        largest_sampled_reference_eV=largest,energetic_mechanical_marks=strong,
                                        support='occupied-D/ray-conditioned peak reference, not any-site patch probability; ordinary population is not coherent amplitude; chart-exit q/v retained',
                                        overlap='Cox marks and shared parent IDs; conservative sampled-reference norm bounds, not new physical energy or a superposition of capped nonlinear trajectories',
                                        uncertainty='simulated parent marks, not measured observations; finite quadrature does not infer unsampled rare-event frequency'))
    m.table('shot_windows.csv',rows);m.table('parent_trials.csv',marks);m.table('whole_foil_rates.csv',rows)
    return rows

def continued_moving(m,summary,samples):
    """Conditional per-mark flight, not a demonstrated escape from a fitted well.
    No averaging of the energetic tail before the nuclear kernel. The original
    signed moving-channel debit limits this alternative continuation once.
    """
    rows=[];index={};gl,gw=roots_legendre(12);fraction=(gl+1)/2
    last={}
    for r in sorted(samples,key=lambda x:x['time_fs']):
        last[r['E_cm_keV'],r['species']-1,r['scenario'],r['environment'],r['branch'],r['control']]=r
    for ec in m.EC:
        theta=m.p1.KB*m.p1.PROFILES[ec]['case']['T_K']
        for sp in (0,1):
            for s in m.p1.SCENARIOS:
                for env in m.D0:
                    with np.load(m.HERE/f'peak_{ec:g}_{sp+1}_{s["name"]}_{env}.npz',allow_pickle=False) as data:
                        marks=data['escape_marks'];positions=data['position_A']
                    sites_per=int(np.prod(m.SHAPE))
                    for branch,(name,chi) in enumerate(m.BRANCHES):
                        for control in (0,1,2):
                            mm=marks[marks[:,1]==branch*3+control]
                            if not len(mm):continue
                            q=mm[:,4:7];v=mm[:,7:10];weight=mm[:,-1]
                            kinetic=.5*m.M*np.sum(v*v,axis=1);r2=np.sum(q*q,axis=1)
                            onsite=.5*m.K*r2+.25*m.B*r2*r2+.25*(m.A-m.B)*np.sum(q**4,axis=1)
                            if env=='vacancy':
                                x=q@m.AXIS[env];site=mm[:,3].astype(int);curves,limit=m.CURVES[env]
                                axquart=m.B+(m.A-m.B)*np.sum(m.AXIS[env]**4)
                                for tag,curve,sign in ((0,curves[1],1.),(3,curves[2],-1.)):
                                    take=site==tag;along=sign*x[take]
                                    onsite[take]+=m.curve_value(curve,along,limit)-(.5*m.K*along**2+.25*axquart*along**4)
                            where=mm[:,2].astype(int)*sites_per+mm[:,3].astype(int)
                            z0=positions[where,2];norm=np.linalg.norm(v,axis=1)
                            uz=np.divide(v[:,2],norm,out=np.zeros_like(norm),where=norm>0)
                            nh=np.empty(len(mm));nd=np.empty(len(mm))
                            for i,z in enumerate(z0):
                                layer=next((r for r in m.p1.PROFILES[ec]['layers'] if r[1]*1e10>=z),m.p1.PROFILES[ec]['layers'][-1])
                                nh[i],nd[i]=layer[2:4]
                            ctrl='built' if control==0 else 'matched_random' if control==1 else 'unforced_prior'
                            rec=None if control==2 else last[ec,sp,s['name'],env,name,ctrl]
                            budget=0. if control==2 else rec['mechanical_signed_work_eV']+rec['electronic_work_eV']-rec['D_coordinate_excess_eV']-rec['coherent_loss_eV']+rec['energy_residual_eV']
                            for release in (0.,1.):
                                energy=kinetic+release*np.maximum(onsite,0.);positive=energy>0;stopping=np.ones(len(energy))
                                stopping[positive]=(nh[positive]*(m.p1.bca.se_ls(energy[positive],46)+m.p1.bca.sn_zbl(energy[positive],46,106.42))+
                                    nd[positive]*(m.p1.bca.se_ls(energy[positive],1)+m.p1.bca.sn_zbl(energy[positive],1,m.p1.bca.M_D)))*1e-18
                                R=np.divide(energy,stopping,out=np.zeros_like(energy),where=positive)
                                zz=z0[:,None]+uz[:,None]*R[:,None]*1e10*fraction
                                density=m.profile_nd(ec,zz)*((zz>=0)&(zz<=m.p1.H*1e10))
                                flight_energy=energy[:,None]*(1-fraction)/1000
                                allocated=float(weight@energy);admitted=float(np.clip(budget/max(allocated,1e-100),-1,1))
                                for ch in ('p','n'):
                                    for U,label in ((0.,'bare'),(.04,'static40')):
                                        score=np.zeros_like(flight_energy)
                                        score[positive]=m.p1.moving_rate(flight_energy[positive],theta,ch,U)
                                        per=R/2*np.sum(gw*density*score,axis=1);total=float(weight@per)
                                        key=(ec,sp,s['name'],env,name,control,ch,label,release)
                                        index[key]=admitted*total/(sp+1)
                                        rows.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,branch=name,
                                            control=ctrl,nuclear_branch=ch,comparator=label,raw_chart_exit_marks=len(mm),weighted_D=float(weight.sum()),
                                            maximum_recorded_kinetic_eV=float(kinetic.max()),maximum_onsite_reference_eV=float(onsite.max()),
                                            conditional_onsite_release_fraction=release,conditional_kinetic_allocation_eV=allocated,
                                            conditional_reactions_per_parent=admitted*total,conditional_reactions_per_D=admitted*total/(sp+1),
                                            raw_flight_reactions_per_parent=total,energy_admitted_signed_fraction=admitted,signed_moving_channel_budget_eV=budget,
                                            energy='actual signed q/v and same onsite potential; interaction partition unknown; optional release is not an added energy source',
                                            flight='12 interior nodes per raw mark; original local LS/ZBL constant-stopping range; actual direction and fixed profile along the ray; foil faces truncate reaction path',
                                            state='conditional moving continuation beyond fitted chart, not demonstrated ejection; no physical escape criterion inferred',
                                            thermal='free-flight endmember is not simultaneously promoted to a changed foil heat source; default carrier rates remain unchanged',
                                            vacancy_prevalence='unit cage endmember; abundance and common-resource redistribution remain conditional'))
    for key in list(index):
        ec,sp,s,env,b,ctrl,ch,label,release=key
        if ctrl in (0,1):index[key]-=index.get((ec,sp,s,env,b,2,ch,label,release),0.)
    m.table('continued_moving.csv',rows)
    return index
