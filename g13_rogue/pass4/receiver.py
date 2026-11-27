"""Shared-parent receiver Green consumer for the Pass-4 producer.
A finite source quadrature, not independent receiver traces or material patches.
"""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import gammainc, gammaincc
import transport as t


def sources(m):
    # Full illuminated square; logarithmic polar quadrature resolves the track
    # receiver while retaining all distant sources. Each cell has physical area.
    ed=np.r_[0.,np.geomspace(1e-9,.0072,65)]
    xy=[];area=[]
    for r0,r1 in zip(ed[:-1],ed[1:]):
        for q in (.211324865405187,.788675134594813):
            r=np.sqrt(r0*r0+q*(r1*r1-r0*r0))
            for j in range(32):
                a=2*np.pi*(j+.5)/32
                x,y=t.L/2+r*np.array([np.cos(a),np.sin(a)])
                if abs(x-t.L/2)<.005 and abs(y-t.L/2)<.005:
                    xy.append([x,y]);area.append(np.pi*(r1*r1-r0*r0)/64)
    xy=np.array(xy);area=np.array(area)
    f=np.load(t.ROOT/'g13_rogue/pass2/replication/inputs.npz')
    cell=np.minimum(((xy-(t.L-.01)/2)/(.01/6)).astype(int),5)
    label=cell[:,0]*6+cell[:,1]
    pp=xy-t.L/2
    raw=np.exp(-np.sum(pp*pp,axis=1)/(2*.003**2))
    for center,gain in (([.0025,-.0025],5.),([-.0025,.0025],3.)):
        raw+=gain*np.exp(-np.sum((pp-center)**2,axis=1)/(2*.0012**2))
    probs=[]
    for hot in (False,True):
        w=area*(raw if hot else 1.)
        w*=f['site_weights'][2 if hot else 0][label]/np.bincount(label,w,minlength=36)[label]
        probs.append(w)
    return xy,np.array(probs),area


def population():
    prof=t.accepted.profiles_and_stopping()[0]['PdD',2.5]
    layers=prof['layers']
    return np.array([t.accepted.layer_value(layers,z) for z in t.DEPTHS])


def slow_background(ref):
    """Reuse actual accepted source histories; no fresh low-band propagation."""
    arrays={};rows=[];strain=[]
    emin=float(np.linalg.eigvalsh(t.CS).min())
    for beam in ('even_steady','hot_steady','hot_bursts'):
        for law in ('retained','Cu_proxy'):
            path=Path(ref)/f'r2_{beam}_b0.0_{law}_combined_records0-15_histories.npz'
            a=np.load(path)
            u=a['local_J_m3'][:,0,[4,13,22]]
            # A conservative strain norm from the actual local elastic quadratic.
            bound=np.sqrt(2*u/emin)
            arrays[beam+'_'+law]=u
            strain.append(float(bound.max()))
            for j,z in enumerate(t.DEPTHS[1:]):
                rows.append(dict(beam=beam,loss=law,depth_m=z,record=0,mean_elastic_J_m3=float(u[:,j].mean()),max_elastic_J_m3=float(u[:,j].max()),strain_norm_upper=float(bound[:,j].max()),source_sha256=t.sha(path),source_path='g13_rogue/pass3/rung2/'+path.name))
            # Reuse exact sunk low-band source/readout rows, not the 30-arm matrix.
            for suffix in ('energy','same_event','tails'):
                p=Path(ref)/f'r2_{beam}_b0.0_{law}_combined_records0-15_{suffix}.csv'
                if p.exists():
                    (t.HERE/('low_'+beam+'_'+law+'_'+suffix+'.csv')).write_bytes(p.read_bytes())
    t.table('slow_background.csv',rows)
    np.savez_compressed(t.HERE/'slow_background.npz',**arrays)
    return max(strain)


def assumed(m):
    out=np.zeros((2,7,3))
    # One broad recoil-normalized reservoir, not a fraction per decade.
    a=20e-9;tau=.5e-12;c=4460.
    scale=np.hypot(tau,a/c)
    def frac(lo,hi):
        a=(2*np.pi*lo*scale)**2;b=(2*np.pi*hi*scale)**2
        return gammainc(1.5,b)-gammainc(1.5,a) if b<1 else gammaincc(1.5,a)-gammaincc(1.5,b)
    norm=gammaincc(1.5,(2*np.pi*2.5e5*scale)**2)
    for j in range(7):out[:,j,:]=m['heatenergy'][:,1,None]*.99*.1*frac(t.BANDS[j],t.BANDS[j+1])/norm/3
    return out


def plate_character(f,branch,theta,grain,z):
    def mode(k):return t.plate_matrices(k,theta,grain)
    lo=.05;hi=max(20.,2*np.pi*f/500)
    while mode(hi)[0][branch]<2*np.pi*f:hi*=2
    k=brentq(lambda k:mode(k)[0][branch]-2*np.pi*f,lo,hi,xtol=1e-5)
    w,v,_,_,_,_,_=mode(k)
    vg=(mode(k*1.001)[0][branch]-mode(k*.999)[0][branch])/(.002*k)
    phi=np.array([t.legval(2*z/t.H-1,np.eye(15)[j])*np.sqrt((2*j+1)/(t.RHO*t.H)) for j in range(15)])
    u=np.array([phi@v[a*15:(a+1)*15,branch] for a in range(3)])
    return k,float(vg),float(t.RHO*np.vdot(u,u).real)


def characteristics(band,frequency=None):
    f=np.sqrt(t.BANDS[band]*t.BANDS[band+1]) if frequency is None else frequency
    if band==0:
        data=np.array([[plate_character(f,b,0.,27,z) for b in range(3)] for z in t.DEPTHS])
        return f,data[:,:,0],data[:,:,1],data[:,:,2]
    c,e=t.christoffel(np.array([0.,0.,1.]))
    k,vg=t.dispersion(f,c)
    return f,np.tile(k,(4,1)),np.tile(vg,(4,1)),np.ones((4,3))/t.H


def anis_path(xy,z,branch,band,k,vg,frequency=None):
    """Receiver-driven piecewise directional travel, through the fixed grain map.
    Straight characteristics are the weak-interface approximation. No fitted
    focusing or hidden grain absorption is added. The bypassed reflected energy
    is a separately reported conservative redistribution remainder.
    """
    d=np.column_stack([t.L/2-xy,np.full(len(xy),z)])
    r=np.linalg.norm(d,axis=1);n=d/np.maximum(r[:,None],1e-30)
    times=np.zeros(len(xy));phase=np.zeros(len(xy));trans=np.ones(len(xy));vprev=None;interface=np.zeros(len(xy))
    for q in (0.0625,.1875,.3125,.4375,.5625,.6875,.8125,.9375):
        at=xy+q*d[:,:2]
        cell=np.minimum(np.maximum((at/(t.L/8)).astype(int),0),7)
        grain=cell[:,0]*8+cell[:,1]
        g=np.einsum('nijkl,nj,nl->nik',t.CT[grain],n,n)/t.RHO
        val,ev=np.linalg.eigh(g)
        vp=np.sqrt(val[:,branch])
        # Directional derivative of omega gives the physical energy-velocity
        # projection along this stationary-path approximation, not vp=vg at THz.
        e=ev[:,:,branch]
        gv=np.einsum('nijkl,nj,nk,nl->ni',t.CT[grain],e,e,n)/(t.RHO*vp[:,None])
        radial=np.einsum('ni,ni->n',gv,n)
        freq=np.sqrt(t.BANDS[band]*t.BANDS[band+1]) if frequency is None else frequency
        speed=radial*np.sqrt(1-(freq/t.FCAP[branch])**2)
        if band==0:speed=np.full_like(speed,vg)
        times+=r/(8*speed)
        phase+=r/8*(2*np.pi*t.FCAP[branch]/vp)*np.arcsin(freq/t.FCAP[branch])
        if vprev is not None:
            rr=((speed-vprev)/(speed+vprev))**2
            trans*=1-rr;interface+=rr
        vprev=speed
    return r,times,trans,interface,phase


def kernel(xy,band,energies,law,offsets=None,frequency=None,bandwidth=None):
    f,kk,vv,depthnorm=characteristics(band,frequency)
    width=t.BANDS[band+1]-t.BANDS[band] if bandwidth is None else bandwidth
    sigma=max(.5e-12,30e-9/float(np.mean(vv)),1/(2*np.pi*width))
    if offsets is None:offsets=np.zeros((1,2))
    vals=[];arrivals=[];branches=[];reflections=[];phases=[];logvals=[]
    for j,z in enumerate(t.DEPTHS):
        vl=[];tl=[];bl=[];rl=[];pl=[];ll=[]
        for b in range(3):
            for echo in range(5):
                # Direct plus two paired front/back paths. All share parent ID.
                zi=z if echo==0 else (2*((echo+1)//2)*t.H+(-1 if echo%2 else 1)*z)
                r,tt,tr,rr,phase=anis_path(xy,zi,b,band,kk[j,b],vv[j,b],f)
                gamma=t.loss(f,vv[j,b],law)
                if band==0:
                    rad=np.sqrt(np.sum((xy-t.L/2)**2,axis=1)+(20e-9)**2+1/kk[j,b]**2)
                    geom=2*np.pi*rad*vv[j,b]/depthnorm[j,b]
                    tt=rad/vv[j,b];phase=rad*kk[j,b]
                    if echo:tr=np.zeros_like(tr)
                else:geom=4*np.pi*(r*r+(20e-9)**2+1/kk[j,b]**2)*vv[j,b]
                # Cu is a per-length coefficient; Al is a per-time lifetime.
                # Guided radial distance and bulk image-ray length are distinct.
                path_length=rad if band==0 else r
                optical=gamma*path_length/vv[j,b] if law=='Cu_extrapolated' else gamma*tt
                # Integrated branch work is retained; this normal-impulse
                # directivity is a coarse angular transport assumption.
                angle=(zi/np.maximum(r,1e-30))**2
                directivity=3*angle if b==2 else 1.5*(1-angle)
                if band==0:directivity=np.ones_like(angle)
                a=energies[:,b][None,:]*directivity[:,None]*tr[:,None]*np.exp(-optical[:,None])/(geom[:,None]*np.sqrt(2*np.pi)*sigma)
                with np.errstate(divide='ignore'):
                    loga=np.log(energies[:,b])[None,:]+np.log(directivity[:,None])+np.log(tr[:,None])-optical[:,None]-np.log(geom[:,None]*np.sqrt(2*np.pi)*sigma)
                a=np.where((optical<=12)[:,None],a,0.)
                vl.append(a);tl.append(tt);bl.append(b);rl.append(echo);pl.append(phase);ll.append(loga)
        vals.append(np.stack(vl,axis=1));arrivals.append(np.stack(tl,axis=1));branches.append(bl);reflections.append(rl);phases.append(np.stack(pl,axis=1));logvals.append(np.stack(ll,axis=1))
    return dict(f=f,k=kk,vg=vv,sigma=sigma,weight=np.stack(vals,axis=1),logweight=np.stack(logvals,axis=1),arrival=np.stack(arrivals,axis=1),branch=np.array(branches),echo=np.array(reflections),phase=np.stack(phases,axis=1))


def unions(intervals):
    out=[]
    for a,b in sorted(intervals):
        if out and a<=out[-1][1]:out[-1][1]=max(b,out[-1][1])
        else:out.append([a,b])
    return out


def count_masks(ker,prob):
    lam=np.zeros(16);s=ker['sigma']
    # Union before count: a reflected path never becomes a new parent.
    for i,p in enumerate(prob):
        events=[]
        for d in range(4):
            ok=ker['weight'][i,d].sum(axis=-1)>0
            for a,b in unions([(v-4*s,v+4*s) for v in ker['arrival'][i,d,ok]]):
                events.extend([(a,1<<d),(b,-(1<<d))])
        state=0;prev=0.
        for age,delta in sorted(events):
            if state:lam[state]+=t.RATE*p*(age-prev)
            state=state+delta if delta>0 else state+delta
            prev=age
    return lam


def moments(ker,prob,energies):
    a=ker['weight']@t.FRACTION
    tau=ker['sigma'];at=ker['arrival']
    mean=t.RATE*np.einsum('n,ndp->d',prob,a)*np.sqrt(2*np.pi)*tau
    cov=np.zeros((4,4))
    for i in range(4):
        for j in range(i+1):
            dt=at[:,i,:,None]-at[:,j,None,:]
            # Same electrical-parent species and all of its echoes stay paired.
            product=np.einsum('nps,nqs,s->npq',ker['weight'][:,i],ker['weight'][:,j],t.FRACTION)
            cv=t.RATE*np.einsum('n,npq,npq->',prob,product,np.exp(-dt*dt/(4*tau*tau)))*np.sqrt(np.pi)*tau
            cov[i,j]=cov[j,i]=cv
    return mean,cov


def draw_normal(p,record,label,size):
    u=p.uniform(record,label,0,2*((size+1)//2)).reshape(-1,2)
    z=np.sqrt(-2*np.log(u[:,0]))[:,None]*np.stack([np.cos(2*np.pi*u[:,1]),np.sin(2*np.pi*u[:,1])],axis=1)
    return z.ravel()[:size]


def burst_at(time,beam,record):
    if beam!='hot_bursts':return np.ones_like(time)
    f=np.load(t.ROOT/'g13_rogue/pass2/replication/inputs.npz')
    a=f['time_weights'][2,:,record%f['time_weights'].shape[-1]]
    # Reuse the accepted Cox-envelope intensity and its real correlation time.
    return np.interp(np.mod(np.asarray(time)+4704*50e-9, len(a)*50e-9),np.arange(len(a))*50e-9,a)


def sparse_record(p,ker,xy,prob,beam,record,band,scenario,law,nd,slowstrain,path_mode='all',refraction=0.):
    f=ker['f'];tau=ker['sigma'];T=64/f;nt=512
    delay=4*tau;time=np.arange(nt)*T/nt+delay
    weights=ker['weight'];at=ker['arrival']+delay
    active=weights.sum(axis=(1,2,3))>0
    source=np.flatnonzero(active)
    supported=weights.sum(axis=-1)>0
    xmin=np.where(supported[active],at[active],np.inf).min(axis=(1,2))-4*tau
    xmax=np.where(supported[active],at[active],-np.inf).max(axis=(1,2))+4*tau+T
    intens=t.RATE*prob[active]*(xmax-xmin)
    t0=(record+.5)*25e-6
    starts=-xmax+T+delay;ends=-xmin+delay
    if beam=='hot_bursts':
        frozen=np.load(t.ROOT/'g13_rogue/pass2/replication/inputs.npz')
        timing=frozen['time_weights'][2,:,record]
        origin=4704*50e-9+t0;dt=50e-9
        sid=[];slo=[];shi=[];rate=[]
        for i,(lo,hi) in enumerate(zip(starts,ends)):
            for cell in range(int(np.floor((origin+lo)/dt)),int(np.floor((origin+hi)/dt))+1):
                if not 0<=cell<len(timing):raise ValueError('Sparse source emission outside accepted record')
                a=max(lo,cell*dt-origin);b=min(hi,(cell+1)*dt-origin)
                if b>a:
                    sid.append(source[i]);slo.append(a);shi.append(b);rate.append(t.RATE*prob[source[i]]*(b-a)*timing[cell])
        source=np.array(sid,int);starts=np.array(slo);ends=np.array(shi);intens=np.array(rate)
    expected=float(intens.sum())
    n=p.poisson(expected,record,1000+band)
    conditional=n==0
    sample_n=1 if conditional else n
    if sample_n>50000:raise ValueError(f'Exact packet catalog exceeds bounded model route: {sample_n} parents; reconsider support/near-field construction')
    if not len(source):
        return None,dict(beam=beam,band=band,scenario=scenario,loss=law,record=record,finite_exposure_s=T,observed_parents=0,conditional=False,probability_nonempty=0.,omitted_energy_survival_bound=np.exp(-12),status='no path within the explicit 12-energy-length census; attenuated tail is censored, not physical zero')
    u=p.uniform(record,f'pass4-parents-{band}',0,sample_n*4).reshape(sample_n,4)
    segment=np.searchsorted(np.cumsum(intens)/expected,u[:,0])
    ids=source[segment]
    times=starts[segment]+u[:,1]*(ends[segment]-starts[segment])
    species=(u[:,2]>=.85).astype(int)
    shuffle=2*np.pi*u[:,3]
    first_arrival=np.min(np.where(supported[ids],at[ids],np.inf),axis=(1,2))
    order=np.argsort(times+first_arrival)
    parent_id=(np.arange(sample_n,dtype=np.int64)+record*1000000000+band*10000000)[order]
    ids=ids[order];times=times[order];species=species[order];shuffle=shuffle[order]
    # Nine neighbouring receivers share these same parents and timing. Plane
    # stationary-phase continuation carries their common paths; no new traces.
    step=max(20e-9,float(np.min(2*np.pi/ker['k']))/8)
    offsets=np.array([[a*step,b*step] for a in (-1,0,1) for b in (-1,0,1)])
    field=np.zeros((nt,4,9,3),complex);control=np.zeros_like(field)
    counts=np.zeros((nt,4),int);participation=np.zeros((nt,4));sumw=np.zeros_like(participation);sumw2=np.zeros_like(participation)
    prelabels=[];arrived=[];predepth=[];pretime=[];preenergy=[]
    for e,(ii,ts,ss) in enumerate(zip(ids,times,species)):
        own=np.zeros_like(field);en=np.zeros((nt,4));support=np.zeros((nt,4),bool)
        for d,z in enumerate(t.DEPTHS):
            for path in range(15):
                amp=weights[ii,d,path,ss]
                if amp<=0:continue
                b=path//5;echo=path%5;arrival=ts+at[ii,d,path]
                if path_mode=='direct' and echo:continue
                if path_mode=='echo' and not echo:continue
                # Carrier sign: normal momentum has one parent sign; stress
                # reflection phase is pi for each free-surface stress echo.
                zi=z if echo==0 else (2*((echo+1)//2)*t.H+(-1 if echo%2 else 1)*z)
                dist=np.sqrt(np.sum((t.L/2-xy[ii])**2)+zi*zi)
                shifted=np.sqrt(np.sum((t.L/2-xy[ii]+offsets)**2,axis=1)+zi*zi)
                delta=shifted-dist
                travel_shift=delta/ker['vg'][d,b]
                first=max(0,int(np.floor((arrival+travel_shift.min()-4*tau-delay)/(T/nt))))
                last=min(nt,int(np.ceil((arrival+travel_shift.max()+4*tau-delay)/(T/nt)))+1)
                if first>=last:continue
                use=slice(first,last);localtime=time[use]
                arg=(localtime[:,None]-arrival-travel_shift[None,:])/tau
                env=np.where(abs(arg)<=4,np.exp(-.25*arg*arg),0.)
                phase=2*np.pi*f*(localtime-ts-delay)-(1+refraction)*ker['phase'][ii,d,path]+echo*np.pi
                wave=np.sqrt(amp)*env*np.exp(1j*(phase[:,None]-ker['k'][d,b]*delta[None,:]))
                own[use,d,:,b]+=wave
                inside=abs(localtime-arrival)<=4*tau
                en[use,d]+=amp*np.exp(-.5*((localtime-arrival)/tau)**2)*inside
                support[use,d]|=inside
        d0,p0=np.unravel_index(np.argmin(np.where(supported[ii],at[ii],np.inf)),at[ii].shape)
        before=ts+at[ii,d0,p0]-4*tau
        midpoint=int(np.argmin(abs(time-before)))
        old=.5*np.sum(abs(field[midpoint,d0,4])**2)
        prior=float(np.mean(.5*np.sum(abs(field[:,d0,4])**2,axis=-1)))
        tag='censored' if before<time[0] or before>time[-1] else 'empty' if prior==0 else 'peak' if old>2*prior else 'hole' if old<.25*prior else 'ordinary'
        prelabels.append(tag);predepth.append(d0);pretime.append(before);preenergy.append(old if tag!='censored' else np.nan)
        field+=own;control+=own*np.exp(1j*shuffle[e]);counts+=support
        sumw+=en;sumw2+=en*en
        arrived.append(int(np.any(support)))
    np.divide(sumw*sumw,sumw2,out=participation,where=sumw2>0)
    # Analytic signal: time-averaged elastic=.5|A|²; real signed strain proxy
    # has the same norm convention. No claim to a resolved microscopic tensor.
    elastic=.5*np.sum(abs(field)**2,axis=-1)
    shuffled=.5*np.sum(abs(control)**2,axis=-1)
    peak=np.unravel_index(elastic.argmax(),elastic.shape)
    ti,dep,point=peak;u0=elastic[:,dep,point];top=float(u0[ti]);half=u0>=top/2 if top else np.zeros(nt,bool)
    lo=hi=ti
    while lo and half[lo-1]:lo-=1
    while hi<nt-1 and half[hi+1]:hi+=1
    map0=elastic[ti,dep].reshape(3,3)
    supportxy=np.argwhere(map0>=top/2) if top else np.zeros((0,2))
    width=(supportxy.max(axis=0)-supportxy.min(axis=0)+1)*step if len(supportxy) else np.zeros(2)
    payload=dict(time_s=time,depth_m=t.DEPTHS,offset_m=offsets,elastic_J_m3=elastic,shuffled_J_m3=shuffled,signed=field.real,counts=counts,effective_count=participation,parent_id=parent_id,parent_species=species,parent_source_index=ids,parent_time_s=times,conditional=np.array(conditional),expected_parents=np.array(expected))
    payload.update(preexisting_label=np.array(prelabels),preexisting_depth_index=np.array(predepth),preexisting_time_s=np.array(pretime),preexisting_elastic_J_m3=np.array(preenergy),packet_filter_delay_s=np.array(delay))
    row=dict(
        beam=beam,band=band,scenario=scenario,loss=law,record=record,
        finite_exposure_s=T,observed_parents=int(n),
        distinct_arrived_parents=int(sum(arrived)) if not conditional else 0,
        conditional=conditional,probability_nonempty=float(-np.expm1(-expected)),
        conditional_record_weight=float(expected*np.exp(-expected)) if conditional else 1.,
        peak_J_m3=top,peak_depth_m=t.DEPTHS[dep],
        peak_over_own_mean=float(top/u0.mean()) if u0.mean()>0 else 0.,
        max_counts=int(counts.max()) if not conditional else 0,
        max_effective_count=float(participation.max()),peak_elastic_eV_D=top/(t.QE*nd[dep]),
        duration_s=(hi-lo+1)*T/nt,duration_censored=bool(lo==0 or hi==nt-1),
        width_x_m=width[0],width_y_m=width[1],width_censored=bool(np.any(width>=3*step)),
        width_resolution_m=step,shuffled_peak_J_m3=float(shuffled.max()),
        labels_peak=prelabels.count('peak'),labels_hole=prelabels.count('hole'),
        labels_empty=prelabels.count('empty'),
        cause='shared parent direct/echo superposition; no exclusive focusing attribution',
        slow_strain_upper=slowstrain,
        phase_sensitivity_bound_Gruneisen10=float(10*slowstrain*np.max(ker['k'])*min(t.L,12*np.max(ker['vg'])/t.loss(f,np.max(ker['vg']),law))))
    row['record_time_offset_s']=delay
    row['labels_censored']=prelabels.count('censored')
    row['coherence_scope']='mode-diagonal signed complex core within12 energy lengths; not full tensor strain or every weak parent tail'
    return payload,row






def retarded_factor(ker,prob,mean,beam):
    if beam!='hot_bursts':return np.ones((128,4))
    frozen=np.load(t.ROOT/'g13_rogue/pass2/replication/inputs.npz')
    timing=frozen['time_weights'][2,:,0];dt=50e-9
    q,w=np.polynomial.legendre.leggauss(12);q*=4
    weights=w*4*np.exp(-q*q/2)*ker['sigma']
    out=np.ones((128,4))
    for d in range(4):
        age=ker['arrival'][:,d,:,None]+ker['sigma']*q
        value=(ker['weight'][:,d]@t.FRACTION)[:,:,None]*weights
        bins=np.floor(age/dt).astype(int)
        lo=int(bins.min());hi=int(bins.max())
        hist=np.bincount((bins-lo).ravel(),weights=(value*prob[:,None,None]*t.RATE).ravel(),minlength=hi-lo+1)
        emit=4704+np.arange(128)[:,None]*15-np.arange(lo,hi+1)[None,:]
        if emit.min()<0 or emit.max()>=len(timing):raise ValueError('Retarded beam exposure exceeds accepted timing')
        if mean[d]>0:out[:,d]=timing[emit]@hist/mean[d]
    return out


def packet_jump_measure(ker,probs,depth):
    """Positive compound-Poisson marginal; all same-parent echoes combined first."""
    q,w=np.polynomial.legendre.leggauss(12)
    values=[];rates=[];ages=[]
    tau=ker['sigma'];at=ker['arrival'][:,depth];weight=ker['weight'][:,depth]
    for i in range(len(probs[0])):
        use=weight[i].sum(axis=-1)>0
        if not np.any(use):continue
        for a,b in unions([(x-4*tau,x+4*tau) for x in at[i,use]]):
            age=(a+b)/2+(b-a)/2*q
            jump=np.exp(-.5*((age[:,None]-at[i,None,:])/tau)**2)@weight[i]
            rate=t.RATE*(b-a)/2*w[:,None]*t.FRACTION
            values.append(jump.ravel());ages.append(np.repeat(np.floor(age/50e-9).astype(int),2))
            rates.append(np.stack([rate.ravel()*probs[j,i] for j in range(2)]))
    if not values:return None
    value=np.concatenate(values);age=np.concatenate(ages);rate=np.concatenate(rates,axis=1)
    take=value>0;value=value[take];age=age[take];rate=rate[:,take]
    low=float(np.log(value.min()));high=float(np.log(value.max()))
    nbin=256;delta=(high-low)/nbin
    index=np.minimum(((np.log(value)-low)/max(delta,1e-300)).astype(int),nbin-1)
    lam=np.array([np.bincount(index,weights=row,minlength=nbin) for row in rate])
    levels=np.zeros((2,nbin))
    for j in range(2):
        np.divide(np.bincount(index,weights=rate[j]*value,minlength=nbin),lam[j],out=levels[j],where=lam[j]>0)
    alo=int(age.min());ahi=int(age.max());nage=ahi-alo+1
    hot_age=np.bincount(index*nage+age-alo,weights=rate[1],minlength=nbin*nage).reshape(nbin,nage)
    return dict(rate=lam,level=levels,hot_age=hot_age,age=np.arange(alo,ahi+1),max_jump_ratio=float(np.exp(delta)),mean=np.sum(lam*levels,axis=1),second=np.sum(lam*levels*levels,axis=1))


def draw_packet_energy(measure,p,beam,band,depth):
    if measure is None:return np.zeros(128),1.
    profile=int(beam!='even_steady')
    rates=np.tile(measure['rate'][profile],(128,1));level=measure['level'][profile]
    if beam=='hot_bursts':
        timing=np.load(t.ROOT/'g13_rogue/pass2/replication/inputs.npz')['time_weights'][2,:,0]
        emit=4704+np.arange(128)[:,None]*15-measure['age'][None,:]
        if emit.min()<0 or emit.max()>=len(timing):raise ValueError('Packet emission outside accepted source record')
        rates=timing[emit]@measure['hot_age'].T
    result=np.zeros(128)
    for record in range(128):
        for j in np.flatnonzero(np.any(rates>0,axis=0)):
            result[record]+=.5*level[j]*p.poisson(float(rates[record,j]),record,9000000+band*10000+depth*1000+int(j))
    return result,measure['max_jump_ratio']


def run(m,ref,only_beam=None):
    for scenario in ('built','assumed'):
        for path in t.HERE.glob(f'{scenario}_b*_r[0-3].npz'):
            if only_beam is None or f'_{only_beam}_' in path.name:path.unlink()
    xy,prob,area=sources(m);nd=population();slowstrain=slow_background(ref)
    p=t.Portable();built=np.load(t.HERE/'launch.npz')['built'][0]
    scenarios={'built':built,'assumed':assumed(m)}
    ledger=[];stats=[];sparse=[];masks=[];observations={};budgets=[]
    geometry={};jumpbank={}
    joint=np.load(t.HERE/'band_depth_joint.npz')
    beams=(only_beam,) if only_beam else ('even_steady','hot_steady','hot_bursts')
    if only_beam and (t.HERE/'readout.csv').exists():
        def retained(name):
            with (t.HERE/name).open() as f:return [row for row in csv.DictReader(f) if row.get('beam')!=only_beam]
        stats=retained('readout.csv');sparse=retained('sparse_events.csv');masks=retained('overlap_counts.csv')
        with np.load(t.HERE/'joint_distributions.npz') as old:observations={key:old[key] for key in old.files if only_beam not in key}
    t.table('source_quadrature.csv',[dict(x_m=x,y_m=y,area_m2=a,even_probability=prob[0,i],hot_probability=prob[1,i]) for i,((x,y),a) in enumerate(zip(xy,area))])
    for scenario,energy in scenarios.items():
        for species in range(2):
            high=float(energy[species].sum());old=.01*m['heatenergy'][species,1] if scenario=='assumed' else 0.
            if high+old>=5000*t.QE:raise ValueError('Source budget overspent; no clipping')
            budgets.append(dict(scenario=scenario,species=('D+','D2+')[species],incident_J=5000*t.QE,low_assumed_J=old,high_acoustic_J=high,retained_heat_random_motion_and_above_ladder_J=5000*t.QE-high-old,note='supported above-ladder acoustic upper is quantified separately'))
        for band in range(7):
            for law in ('Cu_extrapolated','Al_transfer'):
                if scenario=='built':
                    ker=kernel(xy,band,energy[:,band],law);geometry[band,law]=ker
                else:
                    base=geometry[band,law]
                    ratio=(energy[:,band]/built[:,band]).T[np.repeat(np.arange(3),5)]
                    ker=dict(base,weight=base['weight']*ratio[None,None,:,:],logweight=base['logweight']+np.log(ratio)[None,None,:,:])
                f=ker['f'];vg=float(np.mean(ker['vg']));gamma=float(t.loss(f,vg,law));edge=.1*vg/(t.L/2)
                for s in range(2):
                    for b in range(3):
                        power=t.RATE*t.FRACTION[s]*energy[s,band,b];stored=power/(gamma+edge)
                        ledger.append(dict(scenario=scenario,band=band,species=('D+','D2+')[s],branch=b,loss=law,launch_W=power,stored_J=stored,damping_heat_W=stored*gamma,edge_escape_W=stored*edge,balance_W=power-stored*(gamma+edge),grain_absorption_W=0.,loss_provenance='conditional effective removal; no second grain attenuation'))
                for beam in beams:
                    pp=prob[int(beam!='even_steady')]
                    mean,cov=moments(ker,pp,None);lam=count_masks(ker,pp)
                    countmean=np.array([sum(lam[j] for j in range(1,16) if j&(1<<d)) for d in range(4)])
                    eff=np.divide(mean*mean,np.diag(cov),out=np.zeros(4),where=np.diag(cov)>0)
                    countdraw=joint[law+'_'+beam][:,band]
                    label=f'{scenario}_b{band}_{law}_{beam}'
                    observations[label+'_counts']=countdraw
                    observations[label+'_packet_energy_covariance']=cov
                    for d in range(4):
                        masks.append(dict(scenario=scenario,band=band,loss=law,beam=beam,depth_m=t.DEPTHS[d],mean_distinct_count=countmean[d],poisson_joint_min=int(countdraw[:,d].min()),poisson_joint_median=float(np.median(countdraw[:,d])),poisson_joint_max=int(countdraw[:,d].max()),effective_weighted_count=eff[d],support_sigma_s=ker['sigma'],support_halfwidth_s=4*ker['sigma'],propagation_census_energy_lengths=12,paths_per_parent=15))
                    sparse_depth=countmean<200
                    empirical=[[] for _ in range(4)]
                    if np.any(sparse_depth):
                        selected=dict(ker,weight=ker['weight'].copy());selected['weight'][:,~sparse_depth]=0
                        for record in range(4):
                            payload,row=sparse_record(p,selected,xy,pp,beam,record,band,scenario,law,nd,slowstrain);sparse.append(row)
                            path=t.HERE/(label+f'_r{record}.npz')
                            if payload is None:
                                path.unlink(missing_ok=True)
                                for d in range(4):empirical[d].append(np.zeros(512))
                            else:
                                actual=np.zeros_like(payload['elastic_J_m3']) if bool(payload['conditional']) else payload['elastic_J_m3']
                                payload['unconditional_elastic_J_m3']=actual
                                row['signed_attribution_scope']='all supported parents in selected low-count depths'
                                np.savez_compressed(path,**payload)
                                for d in range(4):empirical[d].append(actual[:,d,4])
                    else:
                        for record in range(4):(t.HERE/(label+f'_r{record}.npz')).unlink(missing_ok=True)
                    factors=retarded_factor(ker,pp,mean,beam)
                    val,vec=np.linalg.eigh(cov);root=vec*np.sqrt(np.maximum(val,0))
                    central=np.array([mean*factors[r]+(root@draw_normal(p,r,label,4))*np.sqrt(factors[r]) for r in range(128)])
                    for d in range(4):
                        jump_ratio=1.
                        if sparse_depth[d]:
                            values=np.concatenate(empirical[d]);model='complete shared-parent signed packet records, plus separately weighted exactly-one-parent evidence'
                        elif eff[d]<100:
                            key=(scenario,band,law,d)
                            if key not in jumpbank:jumpbank[key]=packet_jump_measure(ker,prob,d)
                            values,jump_ratio=draw_packet_energy(jumpbank[key],p,beam,band,d)
                            model='positive compound-Poisson parent-packet marginal; cross-depth packet covariance retained separately'
                        else:
                            values=.5*central[:,d]
                            if np.any(values<0):raise ValueError('Central approximation invalid at declared large participation')
                            model='correlated packet-energy central approximation at measured participation>=100; no signed attribution'
                        observations[label+f'_depth{d}_elastic_samples']=values
                        stats.append(dict(scenario=scenario,band=band,lo_Hz=t.BANDS[band],hi_Hz=t.BANDS[band+1],loss=law,beam=beam,depth_m=t.DEPTHS[d],mean_elastic_J_m3=.5*mean[d],elastic_sd_J_m3=.5*np.sqrt(cov[d,d]),mean_elastic_eV_D=.5*mean[d]/(t.QE*nd[d]),conditional_4p1ps_eV_D=mean[d]*(-np.expm1(-float(t.loss(f,ker['vg'][d].mean(),law))*4.1e-12))/(t.QE*nd[d]),n_D_m3=nd[d],q50=float(np.quantile(values,.5)),q95=float(np.quantile(values,.95)),q99=float(np.quantile(values,.99)),finite_max=float(values.max()),max_over_own_mean=float(values.max()/values.mean()) if values.mean()>0 else 0.,zeros=int(np.count_nonzero(values==0)),samples=len(values),distribution_model=model,max_jump_bin_ratio=jump_ratio,finite_sample_exposure_s=4*64/f if sparse_depth[d] else 128*15*50e-9,carrier_Hz=f,energy_lifetime_s=1/float(t.loss(f,ker['vg'][d].mean(),law)),group_speed_m_s=float(ker['vg'][d].mean()),slow_strain_upper=slowstrain))
                t.table('source_budget.csv',budgets);t.table('transport_ledger.csv',ledger);t.table('readout.csv',stats);t.table('overlap_counts.csv',masks);t.table('sparse_events.csv',sparse)
                print(json.dumps(dict(receiver_band=band,scenario=scenario,law=law)),flush=True)
    np.savez_compressed(t.HERE/'joint_distributions.npz',**observations)
    t.dump('receiver_model.json',dict(source_points=len(xy),source_area_quadrature_m2=float(area.sum()),parent_rate_s=t.RATE,support='four Gaussian sigma; twelve conditional energy lengths',count_model='one marked Poisson measure across all28 band/depth cells in band_depth_joint.npz',sparse_model='complete signed catalog selected by actual parent count, never by participation; shared neighbours and shuffled parent phases',dense_model='positive parent-packet compound-Poisson marginals when weights are rare; correlated first/second packet moments and large-participation central approximation otherwise; no universal Rayleigh',burst='accepted Cox intensity at emission time after burn; energy marginals and full joint count census preserve its delay',interface_model='scalar same-branch impedance R+T=1; grain-reflected redistribution is a distinct receiver remainder, not extra heat',edge_model='conditional ten-percent escape per mean half-foil traverse; finite image receiver versus residence remainder kept distinct',phase_model='wave-vector carrier phase; group-velocity envelope transit; stress-echo pi sign',limitations=['guided/bulk and late lateral cavity receiver remainders are conditional, not resonance reproduction','energy jump marginals do not claim a reconstructed cross-depth signed field','spatial widths are shared-neighbour finite/censored measurements','selected frequency transfer refinement supersedes false band-centre tail exclusions']))
