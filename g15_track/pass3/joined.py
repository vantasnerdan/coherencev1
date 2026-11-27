"""Sibling-only, finite-shell and unchanged-carrier joins of the same event."""
import math
import numpy as np
from scipy.linalg import expm
from scipy.special import ndtr


def sibling(m,ec,s,env,coupled):
    """First constituent only, observed BEFORE the second source is inserted."""
    full=m.spatial(ec,1,s);rad=full['position'][:,0];select=np.flatnonzero(rad==rad.min())
    d={k:(value[select] if isinstance(value,np.ndarray) and value.shape and value.shape[0]==len(rad) else value) for k,value in full.items()}
    n=len(select);eig,omega,mw,cosine,local,_,prior=m.modes(env,coupled);nm=len(eig)
    if env=='bulk':q0,v0=m.Q0,m.V0
    else:
        q0=m.HB/(2*m.M*prior)/np.tanh(m.HB*prior/(2*m.p1.KB*300));v0=prior*prior*q0
    q=np.full((len(m.COEF),n,nm),q0);v=np.full_like(q,v0);c=np.zeros_like(q)
    qr=q[0].copy();vr=v[0].copy();cr=c[0].copy();qb=qr.copy();vb=vr.copy();cb=cr.copy()
    y=np.zeros((n,11));A=m.p1.network(s)*1e-15;E=expm(A)
    ext=np.zeros((22,22));ext[:11,:11]=A;ext[:11,11:]=np.eye(11);F=expm(ext)[:11,11:]
    tr=full['tr'];ts,z,_,_=m.p1.source(tr,1);z*=1e10;ts*=1e15
    arrival=np.interp(d['position'][:,2],z,ts)+full['offsets'][1]
    loss=math.exp(-1/(s['tau_opt']*1e15));last=np.zeros((len(m.COEF),n));times=[];quiet=[];delta=[];electrons=[]
    occ=d['occupancy']
    for t in np.arange(0,math.ceil(arrival.max())+2):
        if t:
            k=eig[None,None,:]+last[:,:,None]
            active_dt=np.clip(t-d['start'],0.,1.)[:,None]
            local_loss=np.exp(-active_dt/(s['tau_opt']*1e15))
            q,c,v=m.scalar_transfer(q,c,v,k,active_dt);qr,cr,vr=m.scalar_transfer(qr,cr,vr,eig[None,:],active_dt)
            qb,cb,vb=m.scalar_transfer(qb,cb,vb,eig[None,:],active_dt)
            q=local_loss*q+(1-local_loss)*qb;v=local_loss*v+(1-local_loss)*vb;c=local_loss*c+(1-local_loss)*cb
            qr=local_loss*qr+(1-local_loss)*qb;vr=local_loss*vr+(1-local_loss)*vb;cr=local_loss*cr+(1-local_loss)*cb
            oldopt=y[:,4].copy();y=y@E.T+m.source_rate(d,t-.5,first_only=True)@F.T
            feed=np.maximum(y[:,4]-loss*oldopt,0)/np.maximum(occ,1e-30)
            q+=feed[None,:,None]/(3*k);v+=feed[None,:,None]/(3*m.M)
            qr+=feed[:,None]/(3*eig);vr+=feed[:,None]/(3*m.M)
        last=m.COEF[:,None]*(y[:,3]+y[:,8])[None,:]*local
        times.append(t);quiet.append(vr@mw);delta.append(np.einsum('bnm,m->bn',v-vr[None,:,:],mw));electrons.append(y[:,3]+y[:,8])
    quiet=np.array(quiet);delta=np.array(delta);electrons=np.array(electrons)
    theta=np.array([np.interp(tm,times,m.M*quiet[:,j]) for j,tm in enumerate(arrival)])
    change=np.array([[np.interp(tm,times,m.M*delta[:,b,j]) for j,tm in enumerate(arrival)] for b in range(len(m.COEF))])
    lab=np.interp(d['position'][:,2],z,np.r_[tr['initial_eV'],[x['energy_eV'] for x in tr['timeline']]])/1000
    dz=z[-1]/8*1e-10;column=occ/(m.VOLUME*1e-30)*dz
    out=[]
    for ch in ('p','n'):
        for U,label in ((0.,'bare'),(.04,'static40')):
            for b,coef in enumerate(m.COEF):
                h=np.maximum(theta*1e-3,1e-6);mid=np.maximum(theta+change[b]/2,1e-10)
                derivative=(m.moving_rate(lab,mid+h,ch,U)-m.moving_rate(lab,np.maximum(mid-h,1e-12),ch,U))/(2*h)
                probability=float(column@(derivative*change[b]))
                for beam in m.p1.BEAMS:
                    out.append(dict(E_cm_keV=ec,scenario=s['name'],environment=env,neighbor_on=coupled,coefficient_per_eV=coef,
                        beam=beam,branch=ch,comparator=label,same_parent_delta_per_D2_nucleus=probability/2,
                        mixture_delta_per_incident_D=.15*probability/1.15,delta_per_electrical_D2_parent=probability,
                        electrical_mixture_delta=.15*probability,sibling_delay_fs=full['offsets'][1],axial_separation_A=.74,
                        first_only_max_internal_excitation_eV_formula=max(np.interp(tm,times,electrons[:,j]) for j,tm in enumerate(arrival)),
                        self_response_of_second_included=False,geometry='axial common ray; true half-energy siblings, first constituent only at actual second arrival'))
    return out


def neighbor(m,d,ec,sp,s,env,coupled):
    t=d['t'];w=d['occupied'];info=m.PAIR[env];closing=2*info['u'][1]
    S=d['relq'];V=d['relv'];C=d['relc']
    mean=C/S*closing;variance=np.maximum(V-C*C/S,1e-30);sig=np.sqrt(variance)
    speed=sig*np.exp(-mean*mean/(2*variance))/math.sqrt(2*np.pi)+mean*ndtr(mean/sig)
    density=np.exp(-closing*closing/(2*S))/np.sqrt(2*np.pi*S)
    attempts=np.trapezoid(speed*density,t,axis=0)
    out=[]
    for b,coef in enumerate(m.COEF):
        for eta in (0.,1e-12,1.):
            drive=-np.expm1(-eta*attempts[b]);quiet=-np.expm1(-eta*attempts[0])
            perparent=float((drive-quiet)@w/2)
            for beam in m.p1.BEAMS:
                out.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,neighbor_on=coupled,
                    coefficient_per_eV=coef,beam=beam,native_shell_A=info['native_lowpoint_A'],native_conditional_energy_eV=info['energy'][1],
                    mean_inward_crossings_per_candidate=float(attempts[b]@w/w.sum()),
                    matched_quiet_inward_crossings_per_candidate=float(attempts[0]@w/w.sum()),
                    assumed_contact_transmission=eta,conditional_neighbor_delta_per_electrical_parent=perparent,
                    conditional_mixture_delta_per_incident_D=(.85 if sp==0 else .15)*perparent/1.15,
                    candidates_if_environment_fills_footprint=w.sum()/2,
                    continuation='eta is an unfitted sensitivity, not a native nuclear-contact probability; actual cage abundance and continuation unknown',
                    counting='one disjoint nearest pair per two occupied sites; common translation removed by reciprocal-cell characteristic'))
    return out


def carrier(m,d,ec,sp,s,env,coupled):
    t=d['t']*1e-15;out=[];rate=m.old.current(ec)/m.EV*(.85 if sp==0 else .15)
    for b,coef in enumerate(m.COEF):
        for j,name in enumerate(m.p1.CHANNELS):
            residence=float(np.trapezoid(d['carrier_delta'][:,b,j],t))
            out.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,neighbor_on=coupled,coefficient_per_eV=coef,
                channel=name,changed_residence_eV_s_parent=residence,whole_foil_standing_delta_eV=rate*residence,
                existing_fields='Pass1 whole_foil.csv and foil_*.npz retained; linear signed source correction uses identical onward generator',
                optical_delta_representation='the actual covariance replaces the delta optical reservoir, not an extra copy of it'))
        cov=float(np.trapezoid(d['energy'][:,b],t))
        out.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,neighbor_on=coupled,coefficient_per_eV=coef,
            channel='explicit_covariance_above_matched_quiet',changed_residence_eV_s_parent=cov,whole_foil_standing_delta_eV=rate*cov,
            existing_fields='frozen carrier exits and macroscopic temperature reference unchanged; explicit finite-duration stiffness-work delta',
            optical_delta_representation='contains reversible curvature offset once; not ordinary optical population plus a second initial quiet state'))
    return out


def mixture(m,primary,struck,siblings,neighbor):
    out=[]
    keys=['E_cm_keV','scenario','environment','neighbor_on','coefficient_per_eV','beam','branch','comparator']
    group={}
    for r in primary:
        key=tuple(r[k] for k in keys);group[key]=group.get(key,0.)+r['mixture_delta_per_incident_D']
    for key,beam_delta in group.items():
        k=dict(zip(keys,key))
        matching=lambda r:all(r[x]==k[x] for x in keys if x in r)
        sec=sum((.85 if r['species']==1 else .30)/1.15*r['struck_D_breathing_delta_per_incident_D'] for r in struck if matching(r))
        sib=sum(r['mixture_delta_per_incident_D'] for r in siblings if matching(r))
        required=next(r['G12_required_extra_per_incident_D'] for r in primary if matching(r))
        total=beam_delta+sec+sib
        out.append(dict(**k,other_parent_delta_per_incident_D=beam_delta,struck_delta_per_incident_D=sec,
            sibling_delta_per_incident_D=sib,joined_non_neighbor_delta_per_incident_D=total,
            per_electrical_parent=total*1.15,raw_reactions_per_true_C=total*1.15/m.EV,
            fraction_of_G12_residual=total/required if required else None,
            neighbor_contact_contribution='separate finite-shell transmission and abundance table; not added as a measured nuclear yield'))
    return out


def acoustic_transmission(m,sources):
    """Reuse saved G13 finite-foil Green coefficients, never its producer."""
    saved=m.rows(m.ROOT/'g13_rogue/pass5/readout.csv')
    transport=m.rows(m.ROOT/'g13_rogue/pass5/transport_ledger.csv')
    output=[]
    for source in sources:
        frequency=source['characteristic_frequency_Hz']
        band=int(np.clip(math.floor(math.log10(frequency/2.5e5)),0,6))
        sp=source['species']-1;rate=m.old.current(source['E_cm_keV'])/m.EV*(.85 if sp==0 else .15)
        newpower=rate*source['outgoing_elastic_radiation_eV']*m.EV
        for beam in m.p1.BEAMS:
            for law in ('Cu_extrapolated','Al_transfer'):
                selected=[r for r in transport if r['scenario']=='track' and int(r['band'])==band and r['loss']==law and r['beam']==beam]
                oldpower=sum(float(r['launch_W']) for r in selected)
                own=next(r for r in selected if int(r['species'])==sp)
                escaped=float(own['escape_W'])/float(own['launch_W']) if float(own['launch_W']) else 0.
                converted=float(own['heat_W'])/float(own['launch_W']) if float(own['launch_W']) else 0.
                for r in saved:
                    if r['scenario']!='track' or int(r['band'])!=band or r['loss']!=law or r['beam']!=beam:continue
                    energy_density=float(r['uncensored_mean_J_m3'])*newpower/oldpower if oldpower else 0.
                    output.append(dict(E_cm_keV=source['E_cm_keV'],species=source['species'],scenario=source['scenario'],
                        environment=source['environment'],neighbor_on=source['neighbor_on'],coefficient_per_eV=source['coefficient_per_eV'],
                        beam=beam,loss=law,depth_m=float(r['depth_m']),source_frequency_Hz=frequency,
                        inherited_band_lo_Hz=float(r['lo_Hz']),inherited_band_hi_Hz=float(r['hi_Hz']),
                        new_outgoing_source_W=newpower,conditional_standing_elastic_J_m3=energy_density,
                        conditional_rms_volume_strain=math.sqrt(max(energy_density/m.p1.BULK,0.)),
                        same_event_escaped_eV=source['outgoing_elastic_radiation_eV']*escaped,
                        same_event_onward_Pd_lattice_eV=source['outgoing_elastic_radiation_eV']*converted,
                        outgoing_partition_closure_eV=source['outgoing_elastic_radiation_eV']*(escaped+converted-1),
                        propagation='saved same-band track angular/finite-foil Green transfer, conditional source-shape reuse; no G13 rerun',
                        old_G13_source_added=False))
    return output


def cage_transfer(m):
    out=[]
    for u,energy in zip(m.PAIR['vacancy']['u'][1:],m.PAIR['vacancy']['energy'][1:]):
        k=energy/(2*u*u);omega=math.sqrt(2*k/m.M)
        relative=m.HB/(m.M*omega)/math.tanh(m.HB*omega/(2*m.p1.KB*300))
        out.append(dict(native_displacement_per_D_A=u,native_relative_curvature_eV_A2=k,period_fs=2*np.pi/omega,
            estimated_relative_variance_A2=relative,estimated_single_site_axis_variance_A2=relative/2,
            bulk_moment_transfer_relative_variance_A2=2*m.Q0,
            scope='native finite axial cage harmonic transfer; other four D frozen, transverse Hessian/joint quantum cage state absent; PAW overlap at every point'))
    return out


def strength_and_support(m,d,ec,sp,s,env,coupled):
    out=[];w=d['occupied'];order_weight=w/w.sum()
    for b,coef in enumerate(m.COEF):
        order=np.argsort(d['peak'][b]);cdf=np.cumsum(order_weight[order])
        quantile=np.interp([.5,.9,.99,.999],cdf,d['peak'][b,order])
        maximum=np.unravel_index(np.argmax(d['q'][:,b]),d['q'][:,b].shape)
        ti,site=maximum;variance=float(d['q'][ti,b,site]);quiet=float(d['qr'][ti,site])
        for beam in m.p1.BEAMS:
            out.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,neighbor_on=coupled,
                coefficient_per_eV=coef,beam=beam,one_parent_maximum_site_excess_eV=float(d['peak'][b].max()),
                occupied_site_peak_q50_eV=quantile[0],occupied_site_peak_q90_eV=quantile[1],
                occupied_site_peak_q99_eV=quantile[2],occupied_site_peak_q999_eV=quantile[3],
                maximum_actual_axis_rms_A=math.sqrt(variance),matched_quiet_axis_rms_A=math.sqrt(quiet),
                chart_axis_A=.24,Gaussian_probability_outside_axis_chart=float(2*ndtr(-.24/math.sqrt(variance))),
                matched_quiet_probability_outside_axis_chart=float(2*ndtr(-.24/math.sqrt(quiet))),
                farthest_variance_radius_A=d['position'][site,0],farthest_variance_depth_A=d['position'][site,2],
                farthest_variance_time_fs=d['t'][ti],
                counting='weighted occupied-site strength spread is NOT an event count; the parent maximum is counted once in whole_foil.csv',
                support='positive curvature is not proof that far Gaussian tails stay in the native well; nonlinear cut trajectories are separately censored'))
    return out
