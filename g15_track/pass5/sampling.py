"""Campbell quadrature and finite compound-Poisson receiver distributions.

All copies and bands of one parent are summed before cumulants. Near finite
Poisson marks are never replaced by an rms field. This is a finite radiation
projection, not a bound on the unresolved source-interior elastic near field.
"""
from pathlib import Path
import math
import numpy as np
from scipy.special import ndtr
from common import QE, KB, L, H, ENERGIES, BEAMS, emit, table, inherited, digest, calendar as frozen_calendar
from receiver import Green, FREQUENCY_COMPARATORS, BULK


def marks(green, receiver, beam, refinement=1, returns=True):
    """Positive whole-beam radial/age quadrature, with exact physical area units."""
    s=green.s; rng=np.random.default_rng(8801)
    edges=np.geomspace(2e-10,2*L,25*refinement+1)
    nr=len(edges)-1; na=8*refinement; nu=8*refinement
    # Uniform area in each annulus, uniform azimuth and uniform quantile of
    # a causal exponential age proposal. The mixture includes every port.
    rr=np.sqrt((edges[:-1]**2+edges[1:]**2)/2)
    phi=2*np.pi*(np.arange(na)+.5)/na
    radius,angle,u,port=np.meshgrid(rr,phi,(np.arange(nu)+.5)/nu,np.arange(8),indexing='ij')
    radius=radius.ravel();angle=angle.ravel();u=u.ravel();port=port.ravel()
    xy=np.c_[receiver[0]+radius*np.cos(angle),receiver[1]+radius*np.sin(angle)]
    area=np.repeat(np.pi*np.diff(edges**2)/(na*nu*8),na*nu*8)
    # Proposal front uses the same energy-front law as the receiver, not a
    # compressed physical exposure. Sources still occur over [0,2] seconds.
    width=np.asarray(green.bandwidth)
    from scipy.special import gammaincinv
    # Build proposal clocks from the actual same components, including both
    # species and every image/polarization. A central-speed proxy can miss a
    # ps packet after a ns refractive delay even when its energy is large.
    front_lists=[[] for _ in range(8)]
    for sp in (0,1):
        for pi,omega,bandwidth,arrival,coefficient in green.components(receiver,xy,sp,None,returns=returns):
            front_lists[pi].append(arrival)
    fronts=[np.stack(items,axis=1) for items in front_lists]
    chosen=np.empty(len(radius))
    for pi,values in enumerate(fronts):
        ids=np.flatnonzero(port==pi);pick=rng.integers(values.shape[1],size=len(ids))
        chosen[ids]=values[ids,pick]
    age=chosen+gammaincinv(3,u)/(2*width[port])
    proposal=np.zeros(len(age))
    for pi,values in enumerate(fronts):
        delta=age[:,None]-values;x=width[pi]*np.maximum(delta,0.)
        proposal+=np.mean(np.where(delta>=0,4*width[pi]*x*x*np.exp(-np.minimum(2*x,745.)),0.),axis=1)/8
    norm=s.beam_moments(beam)[0]
    density=s.beam(xy[:,0],xy[:,1],beam)/norm
    physical=s.current(green.ec)/QE*density*area/proposal
    species=(rng.random(len(age))>=.85).astype(int);recoil=rng.integers(0,64,len(age))
    total=np.zeros(len(age),complex);intensity=np.zeros(len(age));bands=np.zeros((len(age),8),complex)
    for sp in (0,1):
        take=species==sp
        if not np.any(take):continue
        h,i,b=green.parent(receiver,xy[take],age[take],sp,recoil[take],returns=returns)
        total[take]=h;intensity[take]=i;bands[take]=b
    return dict(age=age,xy=xy,weight=physical,h=total.real,complex=total,intensity=intensity,bands=bands,
        species=species,recoil=recoil,source_area_m2=float(np.pi*np.sum(np.diff(edges**2))),
        green=green,receiver=np.asarray(receiver),returns=returns)


def distribution(mark, calendar, beam, samples=512, seed=9931, control='joined'):
    """Random-time samples over two physical seconds; event rates are separate."""
    rng=np.random.default_rng(seed)
    times=(np.arange(samples)+rng.random(samples))*calendar.seconds/samples
    h=mark['h'];w=mark['weight'];bands=mark['bands'].real
    if control=='parent_resampled':
        # Independent Poisson parent processes for each band, each with the
        # same intensity. No artificial carrier phases are chosen.
        h=bands.ravel();w=np.repeat(w,8);age=np.repeat(mark['age'],8)
    else:age=mark['age']
    k2=float(w@(h*h));sigma=np.sqrt(k2)
    # Explicit absolute STRAIN resolution, not a comparison with deltaOmega/Omega.
    limit=min(max(sigma*.01,1e-300),1e-12)
    near=abs(h)>limit
    values=np.zeros(samples);dense=np.zeros(samples);cumulants=[];rare_intensity=[]
    actual_mean=mark['green'].mean(mark['receiver'],beam,times,calendar,returns=mark['returns'])
    for j,t in enumerate(times):
        lam=w*calendar.value(t-age,beam)
        variance=float(lam[~near]@(h[~near]**2))
        dense[j]=rng.normal()*np.sqrt(variance)
        active=lam[near];jump=h[near]
        # A finite atom with huge occupation is split only by its actual
        # Poisson law; the normal approximation is NOT used for large jumps.
        count=rng.poisson(active)
        values[j]=actual_mean[j]+dense[j]+float((count-active)@jump)
        if j in (0,samples//2,samples-1):
            cumulants.append([t,float(actual_mean[j]),float(lam@(h*h)),float(lam@(h**3)),float(lam@(h**4))])
        rare_intensity.append(float(active.sum()))
    omitted_variance=float(w[~near]@(h[~near]**2));jumpmax=float(abs(h[~near]).max(initial=0.))
    return dict(time_s=times,strain=values,far=dense,cumulants=np.asarray(cumulants),mean=actual_mean,
        far_variance=omitted_variance,far_maximum_single_parent=jumpmax,
        far_standardized_third_absolute=float(w[~near]@(abs(h[~near])**3))/max(omitted_variance**1.5,1e-300),
        exact_near_mean_active=float(np.mean(rare_intensity)),near_nodes=int(near.sum()))


def frequency_readout(medium,fields,eps,allocation):
    """Actual fixed-T Gamma consumer with a cancellation-safe ratio."""
    T,U,gradT,gradU=fields;p=medium.branches;g=p['linear_strain_slope'];C=p['covariance_coefficient']
    electronic=(np.interp(U,p['electronic_abscissa'],p['electronic_deltaK']) if medium.native_projection is None else
                medium.native_projection.project(np.zeros(3),U,allocation)[0])
    D=p['native_reference']+electronic+g*p['alpha']*(T-300.)
    change=g*np.asarray(eps)/3;changed=D+change
    if np.any(changed<=0):raise ValueError('Signed sample leaves the retained positive-curvature optical chart')
    root=np.sqrt(D*D+4*C*KB*T);root1=np.sqrt(changed*changed+4*C*KB*T)
    evolved=(D+root)/2;evolved1=(changed+root1)/2
    delta=.5*change*(1+(changed+D)/(root1+root))
    shifted=delta/(np.sqrt(evolved1*evolved)+evolved)
    zero=medium._hamiltonian(np.zeros(3),T,U,gradT,gradU,0.,np.zeros(3),allocation,False)[0]
    one=medium._hamiltonian(np.zeros(3),T,U,gradT,gradU,np.asarray(eps)/3,np.zeros(3),allocation,False)[0]
    return shifted,zero,one,g/(6*root),g/(6*evolved)


def build(out,path4,energies=ENERGIES,closures=('fast_Cu','nominal_retained','long_retained'),samples=256,refinement=1,prefix=''):
    from medium import Medium
    from low import Low
    out=Path(out);s,native,fields,launch,response,acoustic,yielding=inherited(path4)
    calendar=frozen_calendar(path4);medium=Medium(out);low=Low(out,path4)
    rows=[];bandrows=[];history=[];controls=[];full_history=[];field_cache={}
    points=np.array([[0.,0.,1e-8],[.0025,-.0025,1e-8],[-.0025,.0025,1e-8],[.009,.009,H/2]])
    p4,warm=FREQUENCY_COMPARATORS
    p=medium.branches;G=abs(p['linear_strain_slope'])/(6*np.sqrt(p['covariance_coefficient']*KB*300.))
    for ec in energies:
        for closure in closures:
            for allocation in ('onsite','coupling'):
                green=Green(out,path4,ec,closure,allocation)
                for beam in BEAMS:
                    for site,xyz in enumerate(points):
                        mark=marks(green,xyz,beam,refinement=refinement)
                        dist=distribution(mark,calendar,beam,samples=samples)
                        low_mean=low.mean(ec,beam,xyz,dist['time_s'],closure,allocation)
                        low_sigma=float(low.innovation_rms_bound(ec,beam,xyz,closure,allocation))
                        eps=dist['strain']+low_mean;stress=BULK*eps
                        key=(ec,beam,site)
                        if key not in field_cache:field_cache[key]=medium.sample(ec,beam,np.broadcast_to(xyz,(samples,3)),time_s=dist['time_s'])
                        background=field_cache[key]
                        shifted,zero,one,fixed,frozen=frequency_readout(medium,background,eps,allocation)
                        variance=float(mark['weight']@(mark['h']**2))
                        maximum=sum(float(abs(mark['bands'][:,pi].real).max()) for pi in range(8))
                        vmax=variance*(11.7/calendar.norm if beam=='hot_bursts' else 1.)
                        total_mean=dist['mean']+low_mean
                        def log_bennett(frequency_threshold):
                            strain_threshold=np.log1p(frequency_threshold)/G
                            u=max(0.,strain_threshold-float(abs(total_mean).max()))
                            if not u:return 0.
                            x=maximum*u/max(vmax,1e-300)
                            return -vmax/max(maximum**2,1e-300)*((1+x)*np.log1p(x)-x)
                        margin=.001*p4;low_error=(G*low_sigma/margin)**2
                        fraction_error=float(np.sqrt(np.log(40)/(2*samples)))
                        row=dict(E_cm_keV=ec,closure=closure,allocation=allocation,beam=beam,receiver=site,
                            x_m=xyz[0],y_m=xyz[1],z_m=xyz[2],physical_calendar_seconds=2.,sample_count=samples,
                            signed_volumetric_strain_min=float(eps.min()),signed_volumetric_strain_max=float(eps.max()),strain_rms=float(np.sqrt(np.mean(eps**2))),
                            tensile_mean_stress_min_Pa=float(stress.min()),tensile_mean_stress_max_Pa=float(stress.max()),stress_rms_Pa=float(np.sqrt(np.mean(stress**2))),
                            fractional_Gamma_frequency_min=float(shifted.min()),fractional_Gamma_frequency_max=float(shifted.max()),fractional_Gamma_frequency_rms=float(np.sqrt(np.mean(shifted**2))),
                            plus_over_pass4=float(shifted.max()/p4),minus_over_pass4=float(-shifted.min()/p4),plus_over_warm=float(shifted.max()/warm),minus_over_warm=float(-shifted.min()/warm),
                            positive_time_fraction_pass4=float(np.mean(shifted>p4)),negative_time_fraction_pass4=float(np.mean(shifted< -p4)),absolute_time_fraction_warm=float(np.mean(abs(shifted)>warm)),
                            fraction_sampling_Hoeffding95_halfwidth=fraction_error,
                            positive_fraction_Pass4_conditional_lower=max(0.,float(np.mean(shifted>p4+margin))-low_error-fraction_error),
                            positive_fraction_Pass4_conditional_upper=min(1.,float(np.mean(shifted>p4-margin))+low_error+fraction_error),
                            low_prompt_fast_pressure_mean_abs_max=float(abs(low_mean).max()),low_prompt_fast_pressure_Campbell_rms_bound=low_sigma,
                            low_innovation_frequency_margin=margin,low_innovation_Chebyshev_probability_bound=min(1.,low_error),
                            high_plus_low_centered_rms_upper=float(np.sqrt(vmax)+low_sigma),
                            Gamma_fixed_T_sensitivity_min=float(fixed.min()),Gamma_fixed_T_sensitivity_max=float(fixed.max()),Gamma_frozen_Q_sensitivity_min=float(frozen.min()),Gamma_frozen_Q_sensitivity_max=float(frozen.max()),
                            tail_differential_absolute_error=float(np.max(abs(shifted-fixed*eps))),returned_Omega_ratio_roundoff=float(np.max(abs(one/zero-1-shifted))),
                            dense_far_third_absolute_ratio=dist['far_standardized_third_absolute'],dense_far_max_parent_strain=dist['far_maximum_single_parent'],
                            exact_near_mean_active=dist['exact_near_mean_active'],source_electrical_parent_rate_s_inverse=s.current(ec)/QE,
                            mean_signed_strain_abs_max=float(abs(total_mean).max()),centered_strain_rms=float(np.std(eps-total_mean)),
                            Campbell_high_total_second=variance,Campbell_high_total_fourth=float(mark['weight']@(mark['h']**4)),
                            scalar_green_sha256=digest(Path(__file__).with_name('receiver.py')),
                            stationary_bulk_compression_potential_J_m3=float(BULK*variance/2),
                            stationary_no_cross_cycle_compression_projection_J_m3=float(BULK/4*(mark['weight']@mark['intensity'])),
                            stationary_compression_coherence_ratio=float(2*variance/max(float(mark['weight']@mark['intensity']),1e-300)),
                            finite_recoil64_log_Bennett_warm=log_bennett(warm),finite_recoil64_log_Bennett_Pass4=log_bennett(p4),
                            comparator_scope='evolved fixed-T Gamma at matched T/U including dQ/dK; actual finite-q consumer is separate',
                            low_scope='prompt+fast-pressure Cox mean and same-parent covariance; correlation-independent remainder bound; full thermoelastic continuation unresolved',
                            acoustic_upcrossing_rate='unresolved: random-time samples are not turning-point events',
                            probability_scope='conditional finite marked radiation model; not a physical exclusion, nuclear contact, or crossing rate',
                            source_interior='source/interior/Fresnel and kr~1 radiation remain unclosed; this finite model is not a physical pointwise bound')
                        rows.append(row)
                        full_history.append(np.stack((dist['time_s'],eps,stress,background[0],background[1],shifted,zero,one,low_mean),axis=-1))
                        for pi,ch in enumerate(green.channels):
                            hh=mark['bands'][:,pi].real;ww=mark['weight']
                            bandrows.append(dict(E_cm_keV=ec,closure=closure,allocation=allocation,beam=beam,receiver=site,port=pi,frequency_Hz=ch['frequency'],
                                single_band_Campbell_variance=float(ww@(hh*hh)),single_band_Campbell_fourth_cumulant=float(ww@(hh**4)),maximum_sampled_single_parent_strain=float(abs(hh).max()),
                                variance_scope='same-parent copies summed within band; complete high cumulants use total h; low correlation bounded without independence'))
                        if ec==2.5 and closure=='nominal_retained' and site==1:
                            shuffled=distribution(mark,calendar,beam,samples=samples,control='parent_resampled')
                            direct=marks(green,xyz,beam,refinement=refinement,returns=False)
                            no_return=distribution(direct,calendar,beam,samples=samples)
                            for name,d in (('joined',dist),('parent_resampled',shuffled),('direct_only',no_return)):
                                ee=d['strain']+low_mean;ff=frequency_readout(medium,background,ee,allocation)[0]
                                controls.append(dict(beam=beam,allocation=allocation,control=name,physical_seconds=2.,sample_count=samples,strain_rms=float(np.sqrt(np.mean(ee*ee))),
                                    maximum=float(ee.max()),minimum=float(ee.min()),Gamma_frequency_min=float(ff.min()),Gamma_frequency_max=float(ff.max()),
                                    positive_Gamma_frequency_fraction_pass4=float(np.mean(ff>p4)),negative_Gamma_frequency_fraction_pass4=float(np.mean(ff< -p4)),
                                    positive_Gamma_frequency_fraction_warm=float(np.mean(ff>warm)),negative_Gamma_frequency_fraction_warm=float(np.mean(ff< -warm)),
                                    probability_scope='matched random-time law over 0–2s, not an event/crossing rate'))
                            weight=mark['weight'];quadratic=float(weight@abs(mark['complex'])**2);incoherent=float(weight@mark['intensity'])
                            controls[-3]['complex_projection_to_no_cross_intensity']=quadratic/max(incoherent,1e-300)
                            controls[-3]['intensity_scope']='trace-strain complex quadratic intensity, NOT total acoustic energy and NOT a negative strain trough'
                            for t,e,f in zip(dist['time_s'],eps,shifted):history.append(dict(beam=beam,allocation=allocation,time_s=t,volumetric_strain=e,tensile_mean_stress_Pa=BULK*e,Gamma_fractional_frequency=f))
                    print(f'receiver E={ec:g} {closure} {allocation} {beam}',flush=True)
    table(out/(prefix+'receiver_statistics.csv'),rows);table(out/(prefix+'band_cumulants.csv'),bandrows)
    if history:table(out/(prefix+'signed_history.csv'),history)
    if controls:table(out/(prefix+'matched_controls.csv'),controls)
    np.savez_compressed(out/(prefix+'receiver.npz'),history=np.array(full_history),
        row_closure=np.array([r['closure'] for r in rows]),row_allocation=np.array([r['allocation'] for r in rows]),row_beam=np.array([r['beam'] for r in rows]),row_energy_keV=np.array([r['E_cm_keV'] for r in rows]),row_receiver=np.array([r['receiver'] for r in rows]),
        columns=np.array(['physical_time_s','volumetric_strain','tensile_mean_stress_Pa','T_K','Ue_eV_formula','Gamma_fractional_frequency','strain_off_Gamma_Omega_rad_s','strained_Gamma_Omega_rad_s','low_prompt_fast_pressure_mean']))
    medium.close()
    return rows


