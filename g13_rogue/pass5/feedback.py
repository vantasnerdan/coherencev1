"""Chronological signed source-footprint projection in the two energetic bands.
A finite physical local catalog plus a declared stationary signed projection of
unresolved whole-foil background; no receiver-later state is used for launch.
"""
import hashlib
import json
import numpy as np
import source as s
from run import angle_lookup, make_kernel, source_phase

t,r=s.t,s.r
HERE=s.HERE
ARMS=('quiet','beta035','shuffled035','relative_speed','no_background035','zero_prior035','zero_prior0')


def launch_phase(source,species,band):
    mean=angle_lookup(np.array([[0.,0.,1.]]),source['center_mean'][0,species,band],source['directions'])[0]
    return source_phase(mean)


def source_power_contrast(source,band):
    angular=source['angular_J_sr'][0,:,band]
    variance=source['angular_variance_J_sr'][0,:,band]
    fraction=np.clip(np.divide(variance,angular,out=np.zeros_like(variance),where=angular>0),0,1)
    mean=np.sqrt(1-fraction)*source_phase(source['center_mean'][0,:,band])
    cross=np.einsum('snb,n->sb',2*angular*np.sqrt(fraction)*mean.real,source['angular_weight_sr'])
    return cross/source['energy_J'][0,:,band]


def earlier_field(launches,green,index):
    return np.einsum('jb,jzab->za',launches[:index],green)


def receiver_packets(source,band,law,beam,position,impact,species,pid,launches,signs,normalization,origin):
    """Actual outgoing launch coefficients enter the same whole-foil packets."""
    xy=t.L/2+position
    ker=make_kernel(xy,band,source['energy_J'][0,:,band],law,0,source)
    phase_bank=np.array([launch_phase(source,sp,band) for sp in range(2)])
    nt=512;span=64/ker['f'];nd=r.population()
    time=np.arange(nt)*span/nt
    offsets=origin+20e-9+t.DEPTHS/float(np.mean(ker['vg']))-span/2
    field=np.zeros((len(ARMS),nt,4,3),complex)
    for i,ts in enumerate(impact):
        sp=species[i]
        for d in range(4):
            clock=time+offsets[d]
            for path in range(15):
                b,echo=divmod(path,5);amp=ker['weight'][i,d,path,sp]
                if amp<=0:continue
                arrival=ts+ker['arrival'][i,d,path]+4*ker['sigma']
                arg=(clock-arrival)/ker['sigma']
                env=np.where(abs(arg)<=4,np.exp(-.25*arg**2),0.)
                if not np.any(env):continue
                coefficient=ker['source_mu'][i,d,path,sp]+np.sqrt(ker['source_var'][i,d,path,sp])*signs[i,b]
                coefficient/=phase_bank[sp,b]*normalization[i,b]
                phase=2*np.pi*ker['f']*(clock-ts-4*ker['sigma'])-ker['phase'][i,d,path]+echo*np.pi
                e=source['energy_J'][0,sp,band,b]
                if e<=0:continue
                wave=np.sqrt(amp)*coefficient*env*np.exp(1j*phase)
                field[:,:,d,b]+=launches[:,i,b,None]/np.sqrt(e)*wave[None,:]
    elastic=.5*np.sum(abs(field)**2,axis=-1)
    rows=[]
    for d,z in enumerate(t.DEPTHS):
        threshold=3*float(elastic[0,:,d].mean())
        for a,arm in enumerate(ARMS):
            y=elastic[a,:,d];mask=y>threshold
            initiations=int(np.count_nonzero(mask & ~np.r_[False,mask[:-1]]))
            rows.append(dict(band=band,loss=law,beam=beam,arm=arm,depth_m=z,
                same_physical_parent_count=len(pid),record_start_s=offsets[d],exposure_s=span,
                mean_elastic_J_m3=float(y.mean()),peak_elastic_J_m3=float(y.max()),
                mean_eV_D=float(y.mean()/(t.QE*nd[d])),peak_eV_D=float(y.max()/(t.QE*nd[d])),
                B_threshold_J_m3=threshold,B_fraction=float(mask.mean()),
                B_initiations=initiations,B_rate_s=initiations/span,
                B_boundary_censored=int(mask[0])+int(mask[-1]),
                status='finite catalog retarded packet record' if threshold>0 else 'no quiet packet sampled in this retarded finite record'))
    return rows,dict(time_offset_s=offsets,time_s=time,elastic_J_m3=elastic,source_parent_id=pid)


def run(m):
    source=s.orient_source(s.load());tracks=s.tracks(m)
    xy,prob,_=r.sources(m)
    kernels={};quiet={}
    for band in (5,6):
        for law in ('Cu_extrapolated','Al_transfer'):
            full=make_kernel(xy,band,source['energy_J'][0,:,band],law,0,source)
            kernels[band,law]=full
            mean,_=r.moments(dict(full,weight=np.exp(full['logweight'])),prob[0],None)
            quiet[band,law]=.5*mean[0]
    nearest=np.argmin(np.sum((xy-t.L/2)**2,axis=1))
    hotfactor=prob[1,nearest]/prob[0,nearest]
    output=[];events=[];arrays={};receiver_rows=[]
    p=t.Portable()
    radius=30e-6;exposure=20e-9;burn=20e-9;area=np.pi*radius**2
    for bi,band in enumerate((5,6)):
        phase_bank=np.array([launch_phase(source,sp,band) for sp in range(2)])
        angular_properties=np.stack([source['angular_J_sr'][0,:,band],
            source['center_mean'][0,:,band],source['angular_variance_J_sr'][0,:,band]],axis=-1)
        f=np.sqrt(t.BANDS[band]*t.BANDS[band+1]);omega=2*np.pi*f
        c,ev,_=s.oriented_polarizations(np.array([0.,0.,1.]));k,vg=t.dispersion(f,c)
        sigma=max(.5e-12,30e-9/float(vg.mean()),1/(2*np.pi*(t.BANDS[band+1]-t.BANDS[band])))
        volume=4*np.pi*((20e-9)**2+1/k**2)*vg*np.sqrt(2*np.pi)*sigma
        for beam in ('even_steady','hot_steady','hot_bursts'):
            # Both bands, both losses, and all launch laws share this catalog.
            origin=.5*25e-6
            factor=1. if beam=='even_steady' else hotfactor
            rate=t.RATE*area/1e-4*factor
            edges=np.array([-burn,0.,exposure])
            intensity=np.array([float(r.burst_at(np.array([origin+(a+b)/2]),beam,0)[0]) for a,b in zip(edges[:-1],edges[1:])])
            duration=np.diff(edges);mass=duration*intensity
            n=p.poisson(rate*mass.sum(),0,88001+('even_steady','hot_steady','hot_bursts').index(beam))
            uu=p.uniform(0,'pass5-causal-'+beam,0,n*10).reshape(n,10)
            segment=np.searchsorted(np.cumsum(mass)/mass.sum(),uu[:,0])
            impact=origin+edges[segment]+uu[:,9]*duration[segment]
            order=np.argsort(impact);impact=impact[order];uu=uu[order]
            pos=radius*np.sqrt(uu[:,1,None])*np.column_stack([np.cos(2*np.pi*uu[:,2]),np.sin(2*np.pi*uu[:,2])])
            species=(uu[:,3]>=.85).astype(int)
            pid=np.arange(n,dtype=np.int64)+('even_steady','hot_steady','hot_bursts').index(beam)*1000000
            source_sign=np.where(uu[:,4:7]<.5,-1.,1.)
            normalization=np.sqrt(1+source_power_contrast(source,band)[species]*source_sign)
            history_sha=hashlib.sha256(pid.tobytes()+impact.tobytes()+pos.tobytes()+species.tobytes()).hexdigest()
            for law in ('Cu_extrapolated','Al_transfer'):
                # Twelve signed plane-wave projections reproduce the diagonal
                # quiet covariance exactly under independent phase averaging.
                # They are background state, not new paid launch parents.
                dirs,_=s.angular_grid(2,2)
                cd,ed,_=s.oriented_polarizations(dirs)
                kd,_=t.dispersion(f,cd)
                background_rotation=np.einsum('ia,nib->nab',ev,ed)
                phase0=2*np.pi*p.uniform(0,'pass5-background-'+beam+'-'+law,0,len(dirs)*3).reshape(len(dirs),3)
                full=kernels[band,law]
                outside=dict(full,weight=np.exp(full['logweight']))
                outside['weight'][np.linalg.norm(xy-t.L/2,axis=1)<=radius]=0
                outside_mean,_=r.moments(outside,prob[int(beam!='even_steady')],None)
                amplitude=np.sqrt(outside_mean[0]/(3*len(dirs)))
                vref=np.sqrt(2*quiet[band,law]/(3*t.RHO))
                uref=vref/omega
                launches=np.zeros((len(ARMS),n,3),complex)
                pre_projection=np.zeros_like(launches)
                pre_velocity=np.zeros((len(ARMS),n))
                causal_case=None
                rowdata={a:dict(self_J=0.,cross_J=0.,work_J=0.,heat_J=0.,
                    expected_work_J=0.,expected_self_J=0.,
                    high=0,signal=[],vz=[],relative=[],cross=[],state=[],scale=[]) for a in ARMS}
                for i in range(n):
                    tr=tracks[species[i]];z=tr['z'];rt=impact[i]+tr['t']
                    weight=tr['deposit'][:,1]/tr['deposit'][:,1].sum()
                    # Spatially resolved source-footprint/track-row projection.
                    if i:
                        delta=pos[i]-pos[:i]
                        dist=np.sqrt(np.sum(delta**2,axis=1)[:,None]+z[None,:]**2)
                        direction=np.concatenate([np.broadcast_to(delta[:,None,:],(i,len(z),2)),np.broadcast_to(z[None,:,None],(i,len(z),1))],axis=2)/dist[:,:,None]
                        green=np.empty((i,len(z),3),complex)
                        cp,ep,_=s.oriented_polarizations(direction)
                        kp,vp=t.dispersion(f,cp)
                        gain=np.empty((i,len(z),3))
                        sourcephase=np.empty((i,len(z),3),complex)
                        for sp in range(2):
                            use=species[:i]==sp
                            if not np.any(use):continue
                            nn=direction[use].reshape(-1,3)
                            values=angle_lookup(nn,angular_properties[sp],source['directions']).reshape((-1,len(z),3,3))
                            aj=values[:,:,:,0].real;cm=values[:,:,:,1];vv=values[:,:,:,2].real
                            ej=source['energy_J'][0,sp,band]
                            gain[use]=np.divide(4*np.pi*aj,ej[None,None,:],out=np.zeros_like(aj),where=ej[None,None,:]>0)
                            frac=np.clip(np.divide(vv,aj,out=np.zeros_like(vv),where=aj>0),0,1)
                            sourcephase[use]=np.sqrt(1-frac)*source_phase(cm)+np.sqrt(frac)*source_sign[:i][use,None,:]
                            sourcephase[use]/=phase_bank[sp][None,None,:]*normalization[:i][use,None,:]
                        for b in range(3):
                            travel=dist/vp[:,:,b]
                            arg=(rt[None,:]-impact[:i,None]-travel-4*sigma)/sigma
                            env=np.where((abs(arg)<=4)&(impact[:i,None]<impact[i]),np.exp(-.25*arg**2),0.)
                            optical=(np.log(10)/10*32*f/1e6)*dist if law=='Cu_extrapolated' else t.loss(f,vp[:,:,b],law)*travel
                            geom=4*np.pi*(dist**2+(20e-9)**2+1/kp[:,:,b]**2)*vp[:,:,b]*np.sqrt(2*np.pi)*sigma
                            green[:,:,b]=np.sqrt(np.maximum(gain[:,:,b],0)/geom)*env*np.exp(-.5*optical+1j*(omega*(rt[None,:]-impact[:i,None]-4*sigma)-kp[:,:,b]*dist))*sourcephase[:,:,b]
                        rotation=np.einsum('ia,jzib->jzab',ev,ep)
                        green=green[:,:,None,:]*rotation
                    # Background evaluated at the actual pre-impact row times and
                    # positions; beam-dependent scalar covariance is declared.
                    bg=np.zeros((len(z),3),complex)
                    for b in range(3):
                        spatial=(dirs[:,:2]@pos[i])[:,None]+dirs[:,2,None]*z[None,:]
                        phases=omega*rt[None,:]-kd[:,b,None]*spatial+phase0[:,b,None]
                        bg+=amplitude*np.einsum('nz,na->za',np.exp(1j*phases),background_rotation[:,:,b])
                    E=source['energy_J'][0,species[i],band]
                    phase=phase_bank[species[i]]
                    Ehost=tr['deposit'][:,1].sum()
                    vD=np.sqrt(2*tr['projectile_eV']*t.QE/t.MD)
                    for a,arm in enumerate(ARMS):
                        prior=np.zeros((len(z),3),complex)
                        if i and arm not in ('zero_prior035','zero_prior0'):prior=earlier_field(launches[a],green,i)
                        if arm not in ('no_background035','zero_prior035','zero_prior0'):prior+=bg
                        if arm=='shuffled035':prior*=np.exp(2j*np.pi*uu[i,8])
                        # Sign survives: velocity is real, displacement quadrature
                        # is imaginary. One electrical parent uses one shared mark.
                        vz=np.sqrt(2/t.RHO)*np.sum(prior.real*ev[2][None,:],axis=1)
                        pre_velocity[a,i]=float(weight@vz)
                        uz=np.sqrt(2/t.RHO)*np.sum(prior.imag*ev[2][None,:],axis=1)/omega
                        signal=.5*np.sum(weight*(np.tanh(vz/max(vref,1e-300))**2+np.tanh(uz/max(uref,1e-300))**2))
                        beta=.35 if arm in ('beta035','shuffled035','no_background035','zero_prior035') else 0.
                        phigh=.01+beta*signal
                        mark=.1 if uu[i,7]<.25 else 10. if uu[i,7]>=1-phigh else 1.
                        scale=np.sqrt(mark/.865)
                        relative=float(np.sum(weight*vz/vD))
                        if arm=='relative_speed':scale*=1-relative
                        # Normalized same-footprint projection into the launch mode.
                        old=np.sum(weight[:,None]*prior,axis=0)*np.sqrt(volume)
                        pre_projection[a,i]=old
                        delta=np.sqrt(E)*normalization[i]*scale*phase
                        selfwork=float(np.sum(abs(delta)**2))
                        crosswork=float(2*np.real(np.vdot(old,delta)))
                        # Physical incoming host reservoir, never a paid cross source.
                        capped=selfwork+crosswork>Ehost
                        if capped:
                            aa=selfwork;bb=crosswork
                            cap=2*Ehost/(bb+np.sqrt(bb*bb+4*aa*Ehost))
                            delta*=cap;selfwork*=cap*cap;crosswork*=cap
                        work=selfwork+crosswork
                        launches[a,i]=delta
                        if causal_case is None and i and np.any(green):
                            causal_case=(i,green.copy())
                        if impact[i]<origin:continue
                        d=rowdata[arm];d['self_J']+=selfwork;d['cross_J']+=crosswork;d['work_J']+=work;d['heat_J']+=Ehost-work
                        expected=expected_self=0.
                        for probability,factor in ((.25,.1),(.75-phigh,1.),(phigh,10.)):
                            da=np.sqrt(E*factor/.865)*normalization[i]*phase
                            if arm=='relative_speed':da*=1-relative
                            aa=float(np.vdot(da,da).real);bb=float(2*np.vdot(old,da).real)
                            if aa+bb>Ehost:
                                cap=2*Ehost/(bb+np.sqrt(bb*bb+4*aa*Ehost))
                                aa*=cap*cap;bb*=cap
                            expected+=probability*(aa+bb);expected_self+=probability*aa
                        d['expected_work_J']+=expected;d['expected_self_J']+=expected_self
                        d['high']+=int(mark==10.);d['signal'].append(signal);d['vz'].append(float(weight@vz));d['relative'].append(relative);d['cross'].append(crosswork);d['state'].append(float(weight@vz));d['scale'].append(scale)
                        events.append(dict(band=band,loss=law,beam=beam,arm=arm,parent_id=int(pid[i]),species=int(species[i]),impact_time_s=impact[i],
                            physical_history_sha256=history_sha,earlier_parents=i,pre_vz_m_s=float(weight@vz),relative_speed_fraction=relative,
                            causal_signal=signal,high_probability=phigh,mark=mark,source_self_J=selfwork,cross_work_J=crosswork,
                            actual_work_J=work,residual_host_J=Ehost-work,cap_applied=bool(capped)))
                baseline=rowdata['quiet']
                for a,arm in enumerate(ARMS):
                    d=rowdata[arm];signal=np.array(d['signal']);v=np.array(d['vz']);cw=np.array(d['cross'])
                    high=signal>=np.median(signal) if len(signal) else np.array([],bool)
                    # C: source-parent signed velocity threshold. B: positive
                    # incoming work cross term in the source projection, not a
                    # falsely named remote receiver event.
                    C=abs(v)>2*vref
                    B=cw>0
                    observed=int(np.count_nonzero(impact>=origin))
                    output.append(dict(band=band,loss=law,beam=beam,arm=arm,parents=observed,burn_parents=n-observed,exposure_s=exposure,
                        physical_history_sha256=history_sha,self_J=d['self_J'],cross_J=d['cross_J'],actual_work_J=d['work_J'],
                        work_over_quiet=d['work_J']/baseline['work_J'] if baseline['work_J'] else 0.,
                        conditional_expected_work_J=d['expected_work_J'],
                        conditional_expected_self_J=d['expected_self_J'],
                        expected_work_over_quiet=d['expected_work_J']/baseline['expected_work_J'] if baseline['expected_work_J'] else 0.,
                        residual_host_J=d['heat_J'],high_marks=d['high'],mean_signal=float(signal.mean()) if n else 0.,
                        max_abs_relative_speed=float(np.max(np.abs(d['relative']))) if n else 0.,
                        positive_cross_parents=int(B.sum()),negative_cross_parents=int((cw<0).sum()),
                        positive_cross_rate_s=float(B.sum()/exposure),C_abs_velocity_parents=int(C.sum()),C_rate_s=float(C.sum()/exposure),
                        C_signed_positive=int((v>2*vref).sum()),C_signed_negative=int((v< -2*vref).sum()),
                        high_state_mean_cross_J=float(cw[high].mean()) if np.any(high) else 0.,
                        low_state_mean_cross_J=float(cw[~high].mean()) if np.any(~high) else 0.,
                        projection_volume_min_m3=float(volume.min()),projection_volume_max_m3=float(volume.max()),
                        quiet_velocity_reference_m_s=vref,zero_self_feedback=True))
                key=f'b{band}_{law}_{beam}'
                arrays[key+'_parent_id']=pid;arrays[key+'_time_s']=impact;arrays[key+'_position_m']=pos;arrays[key+'_species']=species
                arrays[key+'_launch_sqrt_J']=launches
                arrays[key+'_source_marks']=source_sign
                arrays[key+'_source_norm']=normalization
                arrays[key+'_signed_background_phase']=phase0
                arrays[key+'_pre_velocity_m_s']=pre_velocity
                arrays[key+'_pre_projection_sqrt_J']=pre_projection
                arrays[key+'_outside_background_total_J_m3']=outside_mean
                arrays[key+'_background_directions']=dirs
                if causal_case is not None:
                    index,gg=causal_case
                    arrays[key+'_causal_green']=gg
                    arrays[key+'_causal_index']=np.array(index)
                received,payload=receiver_packets(source,band,law,beam,pos,impact,species,pid,launches,source_sign,normalization,origin)
                receiver_rows.extend(received)
                for name,value in payload.items():arrays[key+'_receiver_'+name]=value
                print(json.dumps(dict(feedback_band=band,loss=law,beam=beam,physical_parents=n)),flush=True)
    s.table('feedback.csv',output);s.table('feedback_events.csv',events)
    s.table('feedback_receiver.csv',receiver_rows)
    joined={}
    for row in events:
        key=(row['beam'],row['loss'],row['arm'],row['parent_id'])
        joined.setdefault(key,[]).append(row)
    debit=[]
    for (beam,law,arm,pid),rr in joined.items():
        if len(rr)!=2 or {x['band'] for x in rr}!={5,6}:
            raise ValueError('Both energetic bands must share each physical parent')
        sp=rr[0]['species']
        other=float(source['terms'][0,sp].sum()-source['energy_J'][0,sp,5:7].sum())
        acoustic=other+sum(x['actual_work_J'] for x in rr)
        residual=5000*t.QE-acoustic
        if residual<0:raise ValueError('Joint two-band source overspent incident energy')
        debit.append(dict(beam=beam,loss=law,arm=arm,parent_id=pid,incident_J=5000*t.QE,
            fixed_other_acoustic_J=other,two_band_self_J=sum(x['source_self_J'] for x in rr),
            signed_cross_J=sum(x['cross_work_J'] for x in rr),actual_acoustic_work_J=acoustic,
            residual_heat_J=residual,balance_J=5000*t.QE-acoustic-residual))
    s.table('feedback_debits.csv',debit)
    np.savez_compressed(HERE/'feedback.npz',**arrays)
    s.emit('feedback_model.json',dict(bands=[5,6],arms=ARMS,patch_radius_m=radius,exposure_s=exposure,
        parent_rate='physical electrical current times illuminated-area probability; no compressed time or superparents',
        ordering='all track-row states are pre-existing earlier-parent-only at impact+t_row; no same-parent source feedback',
        inherited_driver='Pass2 live.py426-435/live_method9-21: conditional low/mid/high probability sensitivity, not a velocity law; here signed velocity/displacement source projections replace unavailable exact tensor state',
        marks='p_low=.25, p_high=.01+.35 signal, self-energy factors(.1,1,10)/.865; fixed quiet built-wave source covariance reference',
        background='four signed directions times three branch channels with fixed random phases; diagonal covariance from full-foil mean OUTSIDE the explicit 30-um source-history patch; no_background035 sensitivity; no duplicate paid parents',
        projection='same bulk Green local limit, row-weighted host deposition; V=4pi(a_readout^2+k^-2)vg sqrt(2pi)sigma normalizes source footprint mode; a_pre=sqrt(V)A_pre',
        work='DeltaE=2 Re(a_pre* deltaa)+|deltaa|^2; host remainder debited once, negative work returns to reservoir; no extra paid source for cross work',
        relative_speed='separate physical delta p/p=-v_z/v_D row projection, without arbitrary gain',
        B_C='C is source-parent |pre_vz|>2 fixed quiet RMS; feedback_receiver.csv gives actual finite-catalog retarded receiver B above3 times the fixed quiet same-record mean, initiations/exposure and absolute eV/D; positive-source-cross-work is separately named',
        uncertainty='finite source-footprint modal closure and unresolved signed background covariance, not a measured atomistic feedback coefficient; paired phase shuffle and zero-prior physical null included'))
