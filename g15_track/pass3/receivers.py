"""Positive retarded-age quadrature inside unchanged source cells.

This integrates projections of the LOCAL anomalous-covariance field and their
matched constant-phase control with identical physical weights.
"""
import numpy as np


def rule(m,d,ec,sp,s):
    tr=m.p1.track(ec,sp);ts,z,_,_=m.p1.source(tr,sp);ts*=1e15;z*=1e10
    L=z[-1];R=s['radius']*1e10
    radial_edges=np.unique([0.,.8,3.,R/4,R/2,.75*R,R])
    rr=np.sqrt((radial_edges[:-1]**2+radial_edges[1:]**2)/2)
    zg,zw=np.polynomial.legendre.leggauss(3);rg,rw=np.polynomial.legendre.leggauss(2)
    owners=[];positions=[];weights=[];delays=[];starts=[]
    profile_edges=np.array([x[1]*1e10 for x in m.p1.PROFILES[ec]['layers']])
    for n,pos in enumerate(d['position']):
        cell=int(d['source_segment'][n]);za=cell*L/8;zb=(cell+1)*L/8
        edges=np.unique(np.r_[za,zb,z[(z>za)&(z<zb)],profile_edges[(profile_edges>za)&(profile_edges<zb)]])
        zz=((edges[:-1]+edges[1:])[:,None]+np.diff(edges)[:,None]*zg)/2
        dz=np.diff(edges)[:,None]*zw/2
        nd=m.old.profile_nd(ec,zz.ravel()).reshape(zz.shape)
        ir=int(np.argmin(abs(rr-pos[0])))
        lo,hi=radial_edges[ir:ir+2]**2
        radii=np.sqrt((lo+hi+(hi-lo)*rg)/2);area=np.pi*(hi-lo)*rw/2*1e-20
        local_weights=[];local_pos=[];local_delays=[]
        for radius,ar in zip(radii,area):
            for depth,weight in zip(zz.ravel(),(nd*dz*1e-10*ar).ravel()):
                local_weights.append(weight);local_pos.append([radius,0.,depth])
                local_delays.append(float(np.interp(depth,z,ts)-np.interp(pos[2],z,ts)+(radius-pos[0])/10.))
        local_weights=np.array(local_weights)
        if local_weights.sum()>0:local_weights*=d['occupied'][n]/local_weights.sum()
        for position,weight,delay in zip(local_pos,local_weights,local_delays):
            if weight<=0:continue
            owners.append(n);positions.append(position);weights.append(weight);delays.append(delay)
            starts.append(d['source_clock'][n]+delay)
    return np.array(owners),np.array(positions),np.array(weights),np.array(delays),np.array(starts)


def integrate(m,d,ec,sp,s):
    owner,pos,weight,delay,start=rule(m,d,ec,sp,s)
    h=float(d['t'][1]-d['t'][0]);uniform=np.isclose(d['t']/h,np.rint(d['t']/h),atol=1e-10)
    t=d['t'][uniform];anomaly=d['anomalous'][uniform].astype(complex)
    deltaq=d['delta_q'][uniform].astype(float)
    electron=d['electron'][uniform]
    T,branches,nodes=anomaly.shape;N=d['occupied'].sum();R=s['radius']*1e10
    L=d['position'][:,2].max()*16/15
    phi=np.array([np.ones(len(owner)),np.cos(np.pi*pos[:,2]/L),np.sin(np.pi*pos[:,2]/L),
        np.cos(np.pi*pos[:,0]/R),np.exp(2j*m.old.OMEGA*start)]).T
    gram=phi.conj().T@(weight[:,None]*phi)/N
    eig,U=np.linalg.eigh(gram);keep=eig>1e-10;change=U[:,keep]/np.sqrt(eig[keep])
    outer=phi.conj()[:,:,None]*phi.conj()[:,None,:]
    offset=np.floor(delay/h).astype(int);fraction=delay/h-offset
    lo=int(offset.min());hi=int(offset.max()+1);count=hi-lo+1
    kernels=np.zeros((count,nodes,25),complex);scalar=np.zeros((count,nodes))
    for side,factor in ((0,1-fraction),(1,fraction)):
        slot=offset+side-lo
        np.add.at(scalar,(slot,owner),weight*factor)
        np.add.at(kernels,(slot,owner),(weight*factor)[:,None]*outer.reshape(-1,25))
    pair=np.zeros((T,branches,25),complex);signal=np.zeros((T,branches),complex)
    constant=np.zeros((T,branches));qmean=np.zeros((T,branches))
    emean=np.zeros(T)
    for slot,shift in enumerate(range(lo,hi+1)):
        a=max(0,shift);b=min(T,T+shift)
        if b<=a:continue
        source=slice(a-shift,b-shift);target=slice(a,b)
        pair[target]+=np.einsum('tbn,na->tba',anomaly[source],kernels[slot],optimize=True)
        signal[target]+=np.einsum('tbn,n->tb',anomaly[source],scalar[slot],optimize=True)
        constant[target]+=np.einsum('tbn,n->tb',abs(anomaly[source]),scalar[slot],optimize=True)
        qmean[target]+=np.einsum('tbn,n->tb',deltaq[source],scalar[slot],optimize=True)
        emean[target]+=np.einsum('tn,n->t',electron[source],scalar[slot],optimize=True)
    # Fractional-delay interpolation could otherwise use a post-arrival value
    # just before a subcell front. Replace that short startup band by the exact
    # causal positive quadrature, identically for numerator and norm.
    for n in range(nodes):
        chosen=np.flatnonzero(owner==n)
        if not len(chosen):continue
        earliest=start[chosen].min();latest=start[chosen].max()
        ti=np.flatnonzero((t>=earliest-h)&(t<=latest+h))
        if not len(ti):continue
        oldpair=np.zeros((len(ti),branches,25),complex);oldsignal=np.zeros((len(ti),branches),complex)
        oldnorm=np.zeros((len(ti),branches));oldq=np.zeros((len(ti),branches))
        olde=np.zeros(len(ti))
        for slot,shift in enumerate(range(lo,hi+1)):
            ix=ti-shift;valid=(ix>=0)&(ix<T);safe=np.clip(ix,0,T-1)
            values=anomaly[safe,:,n]*valid[:,None]
            oldpair+=values[:,:,None]*kernels[slot,n]
            oldsignal+=values*scalar[slot,n];oldnorm+=abs(values)*scalar[slot,n]
            oldq+=deltaq[safe,:,n]*valid[:,None]*scalar[slot,n]
            olde+=electron[safe,n]*valid*scalar[slot,n]
        newpair=np.zeros_like(oldpair);newsignal=np.zeros_like(oldsignal);newnorm=np.zeros_like(oldnorm);newq=np.zeros_like(oldq)
        newe=np.zeros(len(ti))
        for sub in chosen:
            query=t[ti]-delay[sub];causal=t[ti]>=start[sub]
            values=np.column_stack([np.interp(query,t,anomaly[:,b,n],left=0.,right=0.) for b in range(branches)])*causal[:,None]
            norm=np.column_stack([np.interp(query,t,abs(anomaly[:,b,n]),left=0.,right=0.) for b in range(branches)])*causal[:,None]
            qq=np.column_stack([np.interp(query,t,deltaq[:,b,n],left=0.,right=0.) for b in range(branches)])*causal[:,None]
            newpair+=values[:,:,None]*(weight[sub]*outer[sub].reshape(25))
            newsignal+=values*weight[sub];newnorm+=norm*weight[sub];newq+=qq*weight[sub]
            newe+=np.interp(query,t,electron[:,n],left=0.,right=0.)*causal*weight[sub]
        pair[ti]+=newpair-oldpair;signal[ti]+=newsignal-oldsignal
        constant[ti]+=newnorm-oldnorm;qmean[ti]+=newq-oldq
        emean[ti]+=newe-olde
    pair=pair.reshape(T,branches,5,5)/N
    pair=np.einsum('ia,tbij,jc->tbac',change.conj(),pair,change.conj(),optimize=True)
    # Empty causal support is analytically zero. Enforce that identity after
    # kernel-minus-startup accumulation, rather than dividing roundoff residue.
    # Positive-support amplitudes, however small, are left untouched.
    empty=t<start.min()
    pair[empty]=0.;signal[empty]=0.;constant[empty]=0.;qmean[empty]=0.;emean[empty]=0.
    # The common captured interval contains no extrapolated post-window field.
    complete=t<=t[-1]+delay.min()-h
    return dict(t=t[complete],field_pairs=pair[complete],signal=signal[complete]/N,
        constant_phase=constant[complete]/N,average_deltaq=qmean[complete]/N,
        average_electron=emean[complete]/N,
        quadrature_position=pos,quadrature_weight=weight,owner=owner,delay=delay,
        receiver_functions=phi@change,first_arrival_fs=float(start.min()),
        complete_window_end_fs=float(t[complete][-1]),
        weighted_count_error_D=float(weight.sum()-N),
        source_depth_samples=len(pos),same_positive_weights_for_phase_and_norm=True)
