"""Causal finite-window stiffness phasors of the actual accepted Hamiltonian.
No canonical q+i v/omega square is identified with a 2-omega pump during a
moving well-centre transient. Polynomial-trend ambiguity remains explicit.
"""
from collections import deque
import numpy as np


def fit_window(t,weights,values,omega,now,harmonic=2):
    """Trailing weighted fit; sample weights are elapsed fs, never observations."""
    period=2*np.pi/omega
    use=t>=now-2*period
    ts=t[use];w=weights[use];y=values[use]
    if len(ts)<8 or ts[-1]-ts[0]<1.6*period:
        return None
    x=(ts-now)/(2*period);phase=harmonic*omega*(ts-now)
    root=np.sqrt(w/max(w.sum(),1e-100))
    result=[]
    for order in (1,2):
        design=np.array([x**k for k in range(order+1)]+[np.cos(phase),np.sin(phase)]).T
        X=design*root[:,None]
        coef,_,rank,s=np.linalg.lstsq(X,y.reshape(len(ts),-1)*root[:,None],rcond=1e-10)
        if rank!=order+3 or s[0]/s[-1]>1e7:return None
        residual=np.sqrt(np.sum(abs((design@coef-y.reshape(len(ts),-1))*root[:,None])**2,axis=0))
        result.append((coef[-2]+1j*coef[-1],coef[-2]-1j*coef[-1],residual,coef[0]))
    a,b,res,offset=result[1]
    # Deterministic trend-model spread, NOT a confidence interval or sampling error.
    spread=abs(a-result[0][0])+np.sqrt(2)*res
    return tuple(z.reshape(y.shape[1:]) for z in (a,b,spread,offset))


def local_mode(m,env):
    """Initial mass-normalized local eigenmode with maximal tagged closing overlap."""
    nsite=8 if env=='bulk' else 6;nd=3*nsite
    H=np.zeros((nd,nd));eps=1e-6
    active=np.zeros((1,1)+m.SHAPE,bool);active.reshape(-1)[:nsite]=True
    for j in range(nd):
        site,c=divmod(j,3);q=np.zeros((1,1,3)+m.SHAPE)
        q.reshape(1,1,3,-1)[0,0,c,site]=eps
        plus=m.lattice_force(q,env,active)[0].reshape(3,-1)[:,:nsite].T.ravel()
        minus=m.lattice_force(-q,env,active)[0].reshape(3,-1)[:,:nsite].T.ravel()
        H[:,j]=-(plus-minus)/(2*eps)
    symmetry=float(abs(H-H.T).max());eigen,U=np.linalg.eigh((H+H.T)/2)
    if eigen.min()<=0:raise ValueError('Declared reference local chart is not a stable well')
    ia,ib=(0,4) if env=='bulk' else (0,3)
    target=np.zeros((nsite,3));target[ia]=m.AXIS[env]/np.sqrt(2);target[ib]=-m.AXIS[env]/np.sqrt(2)
    groups=[];start=0
    for j in range(1,nd+1):
        if j==nd or abs(eigen[j]-eigen[start])>1e-8*max(1.,eigen.max()):
            groups.append(np.arange(start,j));start=j
    group=max(groups,key=lambda ids:np.sum((U[:,ids].T@target.ravel())**2))
    projected=U[:,group]@(U[:,group].T@target.ravel());overlap=float(projected@projected)
    vector=np.zeros((8,3));vector[:nsite]=projected.reshape(nsite,3)/np.sqrt(overlap)
    return vector,np.sqrt(eigen[group].mean()/m.M),dict(
        eigenvalue_eV_A2=eigen.tolist(),frequency_THz=(np.sqrt(eigen/m.M)*1000/(2*np.pi)).tolist(),
        tagged_closing_subspace_overlap=overlap,selected_degenerate_dimension=len(group),
        Hessian_symmetry_defect_eV_A2=symmetry,
        meaning='physical mass metric and this declared local Hamiltonian; not a measured vacancy spectrum')


class CurvatureObserver:
    def __init__(self,m,env,position,weights,plane,omega,clock):
        self.m=m;self.env=env;self.weights=weights;self.N=float(weights.sum())
        self.plane2=plane**2*weights[:,None];self.omega=np.asarray(omega)
        self.vector,self.local_omega,self.local_info=local_mode(m,env)
        self.nb=len(clock);self.block_w=weights.reshape(self.nb,8).sum(1)/(8 if env=='bulk' else 6)
        self.local_phase=np.exp(-1j*self.local_omega*np.asarray(clock))
        self.buffer=deque();self.maximum_window=2*2*np.pi/min(self.local_omega,min(omega))
        self.ia,self.ib=(0,4) if env=='bulk' else (0,3)
        self.norm=np.r_[m.M*self.N*self.omega**2,m.M*self.block_w.sum()*self.local_omega**2]
        self.pair_cross=(plane.reshape(self.nb,8,-1)[:,self.ia]*
            plane.reshape(self.nb,8,-1)[:,self.ib]*self.block_w[:,None])

    def add(self,time,dt,q,active):
        m=self.m;env=self.env
        flat=q.transpose(0,1,3,4,5,2).reshape(len(q),self.nb,8,3)
        L=len(flat);r2=np.sum(flat*flat,axis=-1)
        H=2*m.B*flat[..., :,None]*flat[...,None,:]
        diagonal=m.B*r2[...,None]+3*(m.A-m.B)*flat*flat
        H+=np.einsum('lbna,ac->lbnac',diagonal,np.eye(3))
        reference=np.repeat(H[2::3],3,axis=0)
        H-=reference
        host=(H.reshape(L,-1,9).transpose(0,2,1)@self.plane2).transpose(0,2,1).reshape(L,-1,3,3)
        e=self.vector
        local=np.einsum('na,lbnac,nc->lb',e,H,e,optimize=True)
        ia,ib=self.ia,self.ib;n=m.AXIS[env]
        x=flat[:,:,ia]@n;y=flat[:,:,ib]@n;u=(x-y)/2
        g=m.connected_g(u,env);gp=m.connected_g(u,env,1);gpp=m.connected_g(u,env,2)
        hxx=-y*gp-.25*x*y*gpp;hyy=x*gp-.25*x*y*gpp
        hxy=-g-u*gp+.25*x*y*gpp
        alive=active.reshape(L,self.nb,8);pair_alive=alive[:,:,ia]&alive[:,:,ib]
        corrections=[]
        for h in (hxx,hyy,hxy):
            corrections.append((h-np.repeat(h[2::3],3,axis=0))*pair_alive)
        hxx,hyy,hxy=corrections
        ea=e[ia]@n;eb=e[ib]@n
        local+=hxx*ea*ea+hyy*eb*eb+2*hxy*ea*eb
        phase=self.plane2.reshape(self.nb,8,-1)
        wa=phase[:,ia];wb=phase[:,ib]
        # Equal physical weights at the two physical sites in each microblock.
        pair_tensor=(np.einsum('lb,bm->lm',hxx,wa)+np.einsum('lb,bm->lm',hyy,wb)
            +2*np.einsum('lb,bm->lm',hxy,self.pair_cross))
        host+=pair_tensor[...,None,None]*np.outer(n,n)
        if env=='vacancy':
            curves,limit=m.CURVES[env];axquart=m.B+(m.A-m.B)*np.sum(n**4)
            for idx,c,sgn in ((ia,curves[1],1.),(ib,curves[2],-1.)):
                xx=(flat[:,:,idx]@n)*sgn;xr=np.repeat(xx[2::3],3,axis=0)
                dh=(m.curve_value(c,xx,limit,2)-m.curve_value(c,xr,limit,2)-3*axquart*(xx*xx-xr*xr))*alive[:,:,idx]
                local+=dh*(e[idx]@n)**2
                host+=np.einsum('lb,bm->lm',dh,phase[:,idx])[...,None,None]*np.outer(n,n)
        host/=self.norm[:-1][None,:,None,None]
        local=(local*self.block_w[None,:]*self.local_phase[None,:]**2).sum(1)/self.norm[-1]
        S=np.concatenate([host.reshape(L,len(self.omega),9),np.zeros((L,1,9),complex)],axis=1)
        S[:,-1,0]=local
        corr=flat-np.repeat(flat[2::3],3,axis=0)
        Q=np.einsum('na,lbna->lb',e,corr,optimize=True)
        self.buffer.append((time,dt,S,Q))
        while len(self.buffer)>1 and self.buffer[1][0]<time-self.maximum_window:self.buffer.popleft()

    def read(self,now):
        t=np.array([x[0] for x in self.buffer]);w=np.array([x[1] for x in self.buffer])
        S=np.array([x[2] for x in self.buffer]);Q=np.array([x[3] for x in self.buffer])
        L=S.shape[1];nm=len(self.omega)+1
        h=np.zeros((L,nm));leak=np.zeros_like(h);valid=np.zeros_like(h,dtype=bool)
        available=np.zeros_like(h,dtype=bool)
        for k,om in enumerate(np.r_[self.omega,self.local_omega]):
            fit=fit_window(t,w,S[:,:,k,:],om,now)
            if fit is None:continue
            available[:,k]=True
            a,b,spread,offset=fit
            if k<nm-1:
                h[:,k]=np.linalg.svd(a.reshape(L,3,3),compute_uv=False)[:,0]
                leak[:,k]=np.linalg.norm(spread,axis=1)
            else:h[:,k]=abs(a[:,0]);leak[:,k]=spread[:,0].real
            valid[:,k]=leak[:,k]<=np.maximum(.5*h[:,k],1e-10)
        fit=fit_window(t,w,Q,self.local_omega,now,harmonic=1)
        amplitude=np.zeros_like(Q[-1],complex);spread=np.zeros_like(Q[-1]);offset=np.zeros_like(Q[-1])
        modal_valid=np.zeros_like(Q[-1],bool)
        if fit is not None:
            amplitude,_,spread,offset=fit
            modal_valid=spread<=np.maximum(.5*abs(amplitude),1e-10)
        oscillatory=.5*self.m.M*self.local_omega**2*np.sum(self.block_w[None,:]*abs(amplitude)**2,axis=1)
        total_amp=np.sum(self.block_w[None,:]*abs(amplitude)**2,axis=1)
        coherent=abs(np.sum(self.block_w[None,:]*amplitude*self.local_phase[None,:],axis=1))**2
        fraction=np.divide(coherent,self.block_w.sum()*total_amp,out=np.zeros(L),where=total_amp>1e-100)
        return dict(h=h,leak=leak,valid=valid,window_available=available,local_amplitude=amplitude,
            local_fit_spread=spread,local_fit_valid=modal_valid,local_slow_offset=offset,
            local_oscillatory_reference_eV=oscillatory,local_clock_mode_fraction=fraction)


def discriminator():
    omega=.08;t=np.arange(0.,400.,.5);weights=np.full(len(t),.5)
    # Static translated well curvature, including a slow offset: no 2Omega signal.
    static=np.full((len(t),1),.23)
    zero=fit_window(t,weights,static,omega,t[-1])[0][0]
    expected=.04+.03j
    signal=static+.0001*t[:,None]+.04*np.cos(2*omega*(t-t[-1]))[:,None]+.03*np.sin(2*omega*(t-t[-1]))[:,None]
    recovered=fit_window(t,weights,signal,omega,t[-1])[0][0]
    if abs(zero)>1e-12 or abs(recovered-expected)>1e-12:raise ValueError('Curvature observer discriminator failed')
    return dict(static_translated_curvature_false_h=float(abs(zero)),injected_h=float(abs(expected)),
        recovered_h=float(abs(recovered)),complex_phasor_error=float(abs(recovered-expected)),
        window='causal trailing two mode periods; weighted offset+trend+2omega; quadratic/linear trend spread explicit',
        limitation='observer identity only; injected coefficient is not a measured PdD coupling')
