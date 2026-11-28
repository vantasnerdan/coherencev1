"""Native G16 constitutive law and positive G3 soft-branch packet projection.
No electronic-structure producer, predecessor producer or publication is run.
"""
from functools import lru_cache
import csv,json,math
import numpy as np
import sources as s
HERE=s.HERE
V0=17.029316176e-30
ALLOCATIONS=('onsite','coupling')


class Material:
    def __init__(self,mesh=16):
        with (HERE/'frozen/g16_response.csv').open() as stream:
            rows=[r for r in csv.DictReader(stream) if int(r['kmesh_axis_count'])==mesh]
        rows.sort(key=lambda r:float(r['electronic_excitation_eV_per_D_site']))
        self.mesh=mesh;self.u=np.array([float(r['electronic_excitation_eV_per_D_site']) for r in rows])
        self.dk=np.array([float(r['K_change_eV_A2']) for r in rows]);self.pressure=np.array([float(r['pressure_change_GPa'])*1e9 for r in rows])
        self.reference=float(rows[0]['K_Gamma_100_eV_A2']);self.p0=self.pressure[1]/self.u[1]
    def value(self,u,kind='K'):
        u=np.asarray(u,float)
        if np.any(u< -1e-13) or np.any(u>self.u[-1]*(1+1e-12)):raise ValueError('Electronic field outside native G16 energy support')
        values=self.dk if kind=='K' else self.pressure
        return np.interp(u,self.u,values)
    def slope(self,u,kind='K'):
        values=self.dk if kind=='K' else self.pressure
        index=np.minimum(np.maximum(np.searchsorted(self.u,u,side='right')-1,0),len(self.u)-2)
        return np.diff(values)[index]/np.diff(self.u)[index]


def primitive_basis(p):
    cell=.5*np.array([[1,-1,1],[1,1,0],[0,0,1]],float)
    reciprocal=np.linalg.inv(cell).T*2*np.pi/(p['a_A']*1e-10)
    fractions=np.array(list(__import__('itertools').product((-.5,-.25,0.,.25),repeat=3)))
    translations=np.array(list(__import__('itertools').product((-1,0,1),repeat=3)))
    candidates=(fractions[:,None,:]+translations[None,:,:])@reciprocal.T
    pick=np.argmin(np.sum(candidates*candidates,axis=2),axis=1)
    return candidates[np.arange(len(fractions)),pick]


def dispersion(p,k,allocation='onsite',u=0.,mesh=16,frozen=False):
    """Positive one-Gamma-bin bridge; neither this cutoff nor native q-law is measured."""
    k=np.asarray(k);oldw,oldv,oldD=s.optical(p,k);_,_,Gamma=s.optical(p,np.zeros((1,3)))
    basis=primitive_basis(p);radius=float(np.min(np.linalg.norm(basis,axis=1)[np.linalg.norm(basis,axis=1)>1.]))
    r2=np.sum(k*k,axis=-1)/(radius*radius);inside=r2<1
    W=np.where(inside,(1-r2)**3,0.);gradW=np.where(inside[...,None],-6*(1-r2)[...,None]**2*k/(radius*radius),0.)
    table=Material(mesh);anchor=Material(16).reference;eta=float(Gamma[0]/anchor-1);denom=1+eta*W
    mass=p['mass_eVc2']*s.QE/299792458.**2
    gradold=2*mass*oldw[...,None]*oldv/(s.QE/1e-20)
    D=oldD.copy() if frozen else oldD/denom
    grad=gradold.copy() if frozen else gradold/denom[...,None]-oldD[...,None]*eta*gradW/(denom[...,None]**2)
    first=p['shells'][0];r=np.asarray(first['periodic_alias_offsets_A'])*1e-10
    argument=k@r.T;F=np.cos(argument).mean(axis=-1);gradF=-np.einsum('...b,bi->...i',np.sin(argument),r)/len(r)
    delta=table.value(u)
    if not frozen:
        if allocation=='onsite':D=D+delta
        elif allocation=='coupling':D=D+delta*F;grad=grad+delta[...,None]*gradF
        else:raise ValueError(allocation)
    if np.any(D<=0):raise ValueError('Unstable sampled optical branch; no clipping')
    omega=np.sqrt(D*s.QE/1e-20/mass);vg=grad*s.QE/1e-20/(2*mass*omega[...,None])
    return omega,vg,D,dict(cutoff_q_m_inverse=radius,minimum_wavelength_m=2*np.pi/radius,native_Gamma_K=anchor,
        frozen_Gamma_K=float(Gamma[0]),bridge='Kold/(1+eta W), W=(1-|q/qc|²)^3 inside first nonzero 4³ reciprocal bin; zero outside',
        onsite_K_scale=p['Kii_eV_A2']+(3*p['K4_100_eV_A4']+(3*p['K4_111_eV_A4']-p['K4_100_eV_A4']))*p['zero_point_axis_variance_A2'])


def soft_basis():
    p=s.inputs()['lattice'];k=primitive_basis(p);result={'k':k}
    for alloc in ALLOCATIONS:
        w,v,D,meta=dispersion(p,k,alloc)
        result[alloc+'_omega']=w;result[alloc+'_vg']=v;result[alloc+'_D']=D
    result['old_omega'],result['old_vg'],result['old_D']=s.optical(p,k)
    np.savez_compressed(HERE/'frozen/native_optical_basis.npz',**result)
    return result,meta


def thermal_pole(p,bare,temperature):
    """Classical thermal covariance only: never add the frozen zero-point Q."""
    coefficient=2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4']
    theta=s.KB*np.asarray(temperature)
    if np.any(np.asarray(bare)<=0) or np.any(theta<=0):raise ValueError('Outside retained stable thermal optical branch')
    Q=2*theta/(bare+np.sqrt(bare*bare+4*coefficient*theta))
    evolved=bare+coefficient*Q;denominator=bare+2*coefficient*Q
    return evolved,Q,evolved/denominator,coefficient*s.KB/denominator


def radial_green(distance,age,velocity,mean_path,initial_radius):
    """Positive finite-front ballistic shell plus diffusive moment closure (nm/fs).

    Initial source is a uniform ball with the same peak density as the accepted
    2-nm Gaussian. No Gaussian tail is silently dropped; each kernel has unit norm.
    """
    age=np.asarray(age);R=velocity*age;radius=initial_radius
    rr=np.asarray(distance)
    denom=4*rr*R
    fraction=np.divide(radius*radius-(rr-R)**2,denom,out=np.zeros(np.broadcast_shapes(rr.shape,R.shape)),where=denom>0)
    fraction=np.where((rr==0)|(R==0),((rr+R)<=radius).astype(float),np.minimum(np.maximum(fraction,0.),1.))
    ballistic=fraction/(4*np.pi*radius**3/3)
    tau=mean_path/velocity;b=np.exp(-age/tau)
    diffusion=velocity*mean_path/3
    diffuse_radius=np.minimum(radius+velocity*age,np.sqrt(radius*radius+10*diffusion*age))
    diffuse=(rr<=diffuse_radius)/(4*np.pi*diffuse_radius**3/3)
    return b*ballistic+(1-b)*diffuse


def ball_transform(x):
    x=np.asarray(x);out=np.ones_like(x)
    small=abs(x)<1e-3
    out[small]=1-x[small]**2/10+x[small]**4/280
    out[~small]=3*(np.sin(x[~small])-x[~small]*np.cos(x[~small]))/x[~small]**3
    return out


def cloud(arm,track,sp,time_fs,positions_nm):
    """Actual causal primary electron front in the whole-foil source coordinates.

    Track intervals use two positive Gauss nodes on their elapsed portion.
    The subcell quadrature is a source/field projection, not another spatial domain.
    Thermal and secondary-electron components are added by the shared carrier
    solver; this routine returns their common primary transport Green.
    """
    times=np.r_[0.,[r['time_fs'] for r in track['timeline']]]
    depth=np.r_[0.,[r['depth_A']*.1 for r in track['timeline']]]
    electronic=np.diff(np.r_[0.,[r['electron_eV'] for r in track['timeline']]])
    center=(depth[:-1]+depth[1:])/2
    # Match accepted Gaussian peak, not its unbounded spatial support.
    initial=(3*(2*np.pi)**1.5*2.**3/(4*np.pi))**(1/3)
    out=np.zeros((len(time_fs),len(positions_nm)));ne=np.zeros(len(time_fs));minimum_radius=initial
    offsets=(0.,.74e-10/math.sqrt(2*track['initial_eV']*s.QE/s.MD)*1e15) if sp else (0.,)
    distance=np.sqrt(positions_nm[:,0,None]**2+(positions_nm[:,1,None]-center[None,:])**2)
    mirrored=np.sqrt(positions_nm[:,0,None]**2+(positions_nm[:,1,None]+center[None,:])**2)
    gx=np.array([-.5773502691896257,.5773502691896257])
    for it,t in enumerate(time_fs):
        for offset in offsets:
            elapsed=np.minimum(np.maximum(t-offset-times[:-1],0.),np.diff(times))
            for node in gx:
                emission=times[:-1]+elapsed*(1+node)/2+offset;age=t-emission;active=elapsed>0
                if not np.any(active):continue
                energy=electronic[active]*elapsed[active]/np.diff(times)[active]/2*np.exp(-age[active]/(arm['ne']*1e15))
                kernel=radial_green(distance[:,active],age[active][None,:],1.,arm['lambda_e']*1e9,initial)
                kernel+=radial_green(mirrored[:,active],age[active][None,:],1.,arm['lambda_e']*1e9,initial)
                out[it]+=kernel@energy;ne[it]+=energy.sum()
    return out*V0*1e27,ne,dict(initial_ball_radius_nm=initial,source='exact row energy, elapsed-interval positive GL2, causal ballistic and bounded diffusion; eV per fixed formula',
        front_speed_m_s=1e6,elastic_path_nm=arm['lambda_e']*1e9)


def quadrature(track,radius):
    # Resolves the hot core and retains its far front; nonlinear G16 pressure is
    # integrated on one fixed source quadrature; linear dilute tails stay causal.
    radial_edges=np.r_[0.,np.geomspace(.03,250.,28)]
    tdepth=np.r_[0.,[r['depth_A']*.1 for r in track['timeline']]]
    z=np.unique(np.r_[0.,5.12,np.arange(0.,tdepth[-1]+1.,1.),tdepth[-1]+np.array([5.,15.,40.,100.,250.])])
    zedges=np.r_[0.,(z[:-1]+z[1:])/2,z[-1]+(z[-1]-z[-2])/2]
    r=np.sqrt((radial_edges[:-1]**2+radial_edges[1:]**2)/2)
    rr,zz=np.meshgrid(r,z,indexing='ij');volume=np.pi*np.diff(radial_edges**2)[:,None]*np.diff(zedges)[None,:]
    positions=np.column_stack((rr.ravel(),zz.ravel()));volume=volume.ravel()
    # Optical packet focus is normalized to its actual physical volume/profile.
    focus=np.exp(-positions[:,0]**2/(radius*1e9)**2)
    focus*=((positions[:,1]>=0)&(positions[:,1]<=tdepth[-1]))
    weights=volume*focus;weights/=weights.sum()
    return positions,volume,weights


def gaussian_energy(q,p,w,pw,D,K,A,M):
    return (p*p+pw*pw)/(2*M)+M*D/(2*w*w)+K*(q*q+w*w)/2+A*(q**4+6*q*q*w*w+3*w**4)/4


def gaussian_step(q,p,w,pw,D,K,A,M,dt):
    """Energy-preserving discrete gradient of the Gaussian quartic Hamiltonian."""
    if not np.any(q) and not np.any(p):
        # Exact harmonic rotation at the self-consistent discrete-gradient
        # Hartree pole. Avoid nearly-double elliptic roots without phase lag.
        Q=w*w;C=w*pw/M;V=pw*pw/M**2+D/Q
        count=max(1,int(math.ceil(dt/2.)));step=dt/count
        for substep in range(count):
            lam=3*A*Q
            for iteration in range(16):
                omega=np.sqrt((K+lam)/M);co=np.cos(omega*step);si=np.sin(omega*step)/omega
                Qnew=co*co*Q+2*co*si*C+si*si*V
                target=1.5*A*(Qnew+Q)
                if np.max(abs(target-lam))<2e-14:break
                lam=target
            else:raise ValueError('Nonlinear covariance rotation did not converge')
            Cnew=-(omega*omega)*co*si*Q+(co*co-omega*omega*si*si)*C+co*si*V
            Q=Qnew;C=Cnew;V=(D+C*C)/Q
        return q,p,np.sqrt(Q),M*C/np.sqrt(Q)
    q0=q.copy();w0=w.copy();p0=p.copy();pw0=pw.copy();q=q.copy();w=w.copy();h=dt*dt/(2*M)
    for iteration in range(12):
        gq=K*(q+q0)/2+A*(q**3+q*q*q0+q*q0*q0+q0**3)/4+3*A*(w*w+w0*w0)*(q+q0)/4
        gw=K*(w+w0)/2+3*A*(w**3+w*w*w0+w*w0*w0+w0**3)/4+3*A*(q*q+q0*q0)*(w+w0)/4-M*D*(w+w0)/(2*w*w*w0*w0)
        fq=q-q0-dt*p0/M+h*gq;fw=w-w0-dt*pw0/M+h*gw
        if max(float(abs(fq).max()),float(abs(fw).max()))<2e-13:break
        a=1+h*(K/2+A*(3*q*q+2*q*q0+q0*q0)/4+3*A*(w*w+w0*w0)/4)
        b=h*3*A*w*(q+q0)/2;c=h*3*A*q*(w+w0)/2
        d=1+h*(K/2+3*A*(3*w*w+2*w*w0+w0*w0)/4+3*A*(q*q+q0*q0)/4+M*D*(w+2*w0)/(2*w**3*w0*w0))
        determinant=a*d-b*c;dq=(d*fq-b*fw)/determinant;dw=(a*fw-c*fq)/determinant
        fraction=1.
        while np.any(w-fraction*dw<=0):fraction*=.5
        q-=fraction*dq;w-=fraction*dw
    else:raise ValueError('Native Gaussian discrete-gradient solve did not converge')
    return q,2*M*(q-q0)/dt-p0,w,2*M*(w-w0)/dt-pw0


def add_population(q,p,w,pw,D,K,A,M,energy):
    """Positive covariance injection with exactly the specified per-site energy."""
    Q=w*w;C=w*pw/M;V=pw*pw/M**2+D/Q
    linear=1+3*A*(Q+q*q)/(2*K);quadratic=3*A/(4*K*K)
    amount=2*energy/(linear+np.sqrt(linear*linear+4*quadratic*energy))
    Qnew=Q+amount/K;Vnew=V+amount/M
    Dnew=Qnew*Vnew-C*C
    return np.sqrt(Qnew),M*C/np.sqrt(Qnew),Dnew


def damp_population(q,p,w,pw,D,K,A,M,dt,kappa,theta):
    """Positive phase-space relaxation: κ is energy loss in the harmonic limit.
    Nonlinear energy loss is measured and credited, never assumed or clipped.
    """
    Q=w*w;C=w*pw/M;V=pw*pw/M**2+D/Q
    Qquiet=2*theta/(K+np.sqrt(K*K+12*A*theta))
    decay=np.exp(-kappa*dt);mix=-np.expm1(-kappa*dt)
    Q=decay*Q+mix*Qquiet;C=decay*C;V=decay*V+mix*theta/M
    return q*np.sqrt(decay),p*np.sqrt(decay),np.sqrt(Q),M*C/np.sqrt(Q),Q*V-C*C



def source_basis(arm,ec,sp,track):
    """One reusable causal source projection for both allocations and meshes."""
    import run as model
    positions,volume,weights=quadrature(track,arm['radius'])
    stop=track['terminal_time_fs']+(.74e-10/math.sqrt(2*track['initial_eV']*s.QE/s.MD)*1e15 if sp else 0.)
    end=stop+14*max(arm['ne'],arm['opt'],arm['end'])*1e15
    first=np.arange(0.,stop+100.,.5);dense=np.r_[first,np.arange(first[-1]+2.,stop+6*arm['ne']*1e15+2.,2.)]
    times=np.unique(np.r_[dense,np.geomspace(max(dense[-1],1.),end,160)])
    primary,primary_energy,geometry=cloud(arm,track,sp,times,positions)
    tt,h,length,_=model.optical_histogram(track,sp,arm,times*1e-15)
    population=np.column_stack([np.interp(times,tt*1e15,h[:,j]) for j in range(9)])
    rr=positions[:,0];radius2=rr*rr+positions[:,1]**2
    secondary_radius=np.sqrt((arm['radius']*1e9)**2+2*1e6*arm['lambda_e']*max(arm['ne'],arm['te'])*1e18)
    row_time=np.array([r['time_fs'] for r in track['timeline']]);row_depth=np.array([r['depth_A']*.1 for r in track['timeline']])
    centers=np.interp(times,row_time,row_depth)/2
    secondary_shape=np.exp(-(rr[None,:]**2+(positions[None,:,1]-centers[:,None])**2)/(2*secondary_radius**2))
    secondary_shape*=radius2[None,:]<(geometry['initial_ball_radius_nm']+times[:,None])**2
    secondary_shape/=secondary_shape@volume[:,None]
    residual_ne=population[:,3]-primary_energy
    quadrature_error=float(np.max(abs(residual_ne-np.maximum(residual_ne,0.))))
    u=primary+(residual_ne+population[:,8])[:,None]*secondary_shape*V0*1e27
    return positions,volume,weights,times,geometry,tt,h,length,population,secondary_shape,quadrature_error,u


def packet(arm,ec,sp,track,mesh=16,basis=None):
    """Joined constitutive packet: actual source clock, moving electron Green,
    positive optical covariance, reciprocal native frequency work, causal ports.

    The packet is a retained source/phase-space basis of the whole foil. It has
    no independent beam, temperature boundary, exposure or spatial domain.
    """
    import run as model
    from scipy.linalg import expm
    from scipy.optimize import brentq
    p=s.inputs()['lattice'];allocation=arm['allocation'];material=Material(mesh)
    if basis is None:basis=source_basis(arm,ec,sp,track)
    positions,volume,weights,times,geometry,tt,h,length,population,secondary_shape,quadrature_error,u=basis
    weights=weights.copy()
    # The table is applied to actual local internal energy per formula, not
    # coarse energy divided by the occupied-D fraction.
    k=primitive_basis(p);oldw,oldv,oldD=s.optical(p,k)
    omega,vg,Dbase,bridge=dispersion(p,k,allocation,mesh=16)
    shell=np.asarray(p['shells'][0]['periodic_alias_offsets_A'])*1e-10
    F=np.cos(k@shell.T).mean(axis=1);sensitivity=np.ones(len(k)) if allocation=='onsite' else F
    profile=next(r for r in s.inputs()['profiles'] if r['case']['E_cm_keV']==ec)
    density=np.array([next((r[3] for r in profile['layers'] if r[0]<=z*1e-9<r[1]),0.) for z in positions[:,1]])
    weights*=density;weights/=weights.sum()
    column=sum(max(0.,min(length,r[1])-r[0])*r[3] for r in profile['layers'] if r[0]<length)
    sites=column*np.pi*arm['radius']**2
    mass=p['mass_eVc2']*s.QE/299792458.**2;M=mass*1e10/s.QE
    axial=p['K4_100_eV_A4'];cross=(3*p['K4_111_eV_A4']-axial)/2
    A=(3*axial+2*cross)/3;theta=s.KB*300;dofsites=3*sites
    Q=2*theta/(Dbase+np.sqrt(Dbase*Dbase+12*A*theta));w=np.sqrt(Q);q=np.zeros(len(k));mom=q.copy();pw=q.copy();det=Q*theta/M
    q0=q.copy();p0=mom.copy();w0=w.copy();pw0=pw.copy();det0=det.copy()
    wc=w.copy();pwc=pw.copy();detc=det.copy();off_correction=np.zeros(11);off_records=[];off_states=[];off_stress=[];off_pressures=[]
    thermal_work=[];off_thermal_work=[]
    _,_,_,thermal_slope=thermal_pole(p,Dbase,300.)
    Ab=model.network(arm)*1e-15;transport=np.zeros((11,11));transport[:9,:9]=Ab
    transport[:,4]=0.;transport[4,:]=0.;transport[9,:9]=Ab[4];transport[9,4]=0.;transport[10,:9]=-Ab.sum(axis=0)
    correction=np.zeros(11);previous_delta=np.zeros(len(k));propagators={}
    cumulative_incoming=np.interp(times,tt*1e15,h[:,9])
    records=[];coarse_maps=[];work_sum=0.;loss_sum=0.;maxdef=0.;maxglobal=0.;phase=np.zeros(len(k));phase_ref=np.zeros(len(k));last=0.;stress=[];states=[];pressures=[];dk_volume=[]
    for it,t in enumerate(times):
        if it==0:continue
        dt=t-last;last=t
        key=round(float(dt),10)
        if key not in propagators:propagators[key]=expm(transport*dt)
        correction=propagators[key]@correction
        off_correction=propagators[key]@off_correction
        off_incoming=off_correction[9];off_correction[9]=0.
        delta_incoming=correction[9];correction[9]=0.
        state=population[it]+correction[:9]
        electron=state[3]+state[8]
        if electron<=0:raise ValueError('Native work depleted total local electronic excess')
        norm=population[it,3]+population[it,8]
        actual_u=u[it]*(electron/norm) if norm>0 else np.zeros(len(volume))
        prior_energy=gaussian_energy(q,mom,w,pw,det,Dbase+previous_delta,A,M)
        prior_quiet=gaussian_energy(q0,p0,w0,pw0,det0,Dbase,A,M)
        # Nonlinear onsite/neighbor frequency work is credited oppositely to
        # the exact same electronic carrier through an implicit energy root.
        second=q*q+w*w
        def delta_at(scale):return float(weights@material.value(actual_u*scale))*sensitivity
        def balance(scale):
            change=delta_at(scale)-previous_delta
            return (scale-1)*electron+dofsites*.5*float(np.mean(second*change))
        upper=2.
        while balance(upper)<0:upper*=2.
        scale=brentq(balance,0.,upper,xtol=1e-13,rtol=1e-13)
        delta=delta_at(scale);work=dofsites*.5*float(np.mean(second*(delta-previous_delta)))
        correction[[3,8]]-=work*state[[3,8]]/electron
        actual_u*=scale;state[[3,8]]*=scale;work_sum+=work
        K=Dbase+delta
        q,mom,w,pw=gaussian_step(q,mom,w,pw,det,K,A,M,dt)
        q0,p0,w0,pw0=gaussian_step(q0,p0,w0,pw0,det0,Dbase,A,M,dt)
        _,_,wc,pwc=gaussian_step(q0,p0,wc,pwc,detc,Dbase,A,M,dt)
        # The incoherent source is a positive phase-space covariance increment;
        # no universal electronic phase is invented.
        base_incoming=cumulative_incoming[it]-cumulative_incoming[it-1]
        injection=base_incoming+delta_incoming
        if injection< -1e-10:raise ValueError('Negative optical source')
        w,pw,det=add_population(q,mom,w,pw,det,K,A,M,injection/dofsites)
        wc,pwc,detc=add_population(q0,p0,wc,pwc,detc,Dbase,A,M,(base_incoming+off_incoming)/dofsites)
        beforec=gaussian_energy(q0,p0,wc,pwc,detc,Dbase,A,M)
        before=gaussian_energy(q,mom,w,pw,det,K,A,M)
        before0=gaussian_energy(q0,p0,w0,pw0,det0,Dbase,A,M)
        q,mom,w,pw,det=damp_population(q,mom,w,pw,det,K,A,M,dt,1/(arm['opt']*1e15),theta)
        q0,p0,w0,pw0,det0=damp_population(q0,p0,w0,pw0,det0,Dbase,A,M,dt,1/(arm['opt']*1e15),theta)
        _,_,wc,pwc,detc=damp_population(q0,p0,wc,pwc,detc,Dbase,A,M,dt,1/(arm['opt']*1e15),theta)
        after=gaussian_energy(q,mom,w,pw,det,K,A,M);after0=gaussian_energy(q0,p0,w0,pw0,det0,Dbase,A,M)
        excess=dofsites*float(np.mean(after-after0));loss=dofsites*float(np.mean((before-after)-(before0-after0)))
        afterc=gaussian_energy(q0,p0,wc,pwc,detc,Dbase,A,M)
        off_excess=dofsites*float(np.mean(afterc-after0));off_loss=dofsites*float(np.mean((beforec-afterc)-(before0-after0)))
        # Replace, do not add to, the frozen optical damping ports.
        inherited_loss=base_incoming-(population[it,4]-population[it-1,4])
        correction[8]+=.2*(loss-inherited_loss);correction[6]+=.8*(loss-inherited_loss)
        off_correction[8]+=.2*(off_loss-inherited_loss);off_correction[6]+=.8*(off_loss-inherited_loss)
        loss_sum+=loss
        expected=dofsites*float(np.mean(prior_energy-prior_quiet))+injection+work-loss
        maxdef=max(maxdef,abs(excess-expected))
        global_defect=float(correction[:9].sum()+excess-population[it,4]+correction[10]);maxglobal=max(maxglobal,abs(global_defect))
        state=population[it]+correction[:9];state[4]=excess;states.append(state.copy())
        stress.append(dofsites*p['collective_strain_slope_eV_A2']*float(np.mean(q*q+w*w-q0*q0-w0*w0))/6)
        thermal_work.append(dofsites*.5*float(np.mean((q*q+w*w-q0*q0-w0*w0)*(p['collective_strain_slope_eV_A2']*s.ALPHA+thermal_slope))))
        # The analytic linear part includes the full causal carrier population;
        # only the nonlinear correction needs the finite near-front quadrature.
        pressures.append(material.p0*(state[3]+state[8])*V0+float((material.value(actual_u,'P')-material.p0*actual_u)@(volume*1e-27)))
        dk_volume.append(float(material.value(actual_u)@(volume*1e-27)))
        off_state=population[it]+off_correction[:9];off_state[4]=off_excess;off_states.append(off_state.copy())
        off_stress.append(dofsites*p['collective_strain_slope_eV_A2']*float(np.mean(wc*wc-w0*w0))/6)
        off_thermal_work.append(dofsites*.5*float(np.mean((wc*wc-w0*w0)*(p['collective_strain_slope_eV_A2']*s.ALPHA+thermal_slope))))
        off_u=u[it]*((off_state[3]+off_state[8])/norm) if norm>0 else np.zeros(len(volume))
        off_pressures.append(material.p0*(off_state[3]+off_state[8])*V0+float((material.value(off_u,'P')-material.p0*off_u)@(volume*1e-27)))
        variance=pw*pw/M**2+det/(w*w)
        theta_modes=M*variance;omega_eff=np.sqrt((K+3*A*(q*q+w*w))/M)*1e15
        gradbase=2*mass*omega[:,None]*vg/(s.QE/1e-20)
        gradF=-np.einsum('kb,bi->ki',np.sin(k@shell.T),shell)/len(shell)
        gradlive=gradbase+(float(weights@material.value(actual_u))*gradF if allocation=='coupling' else 0.)
        livevg=gradlive*s.QE/1e-20/(2*mass*omega_eff[:,None])
        phase+=omega_eff*dt*1e-15;phase_ref+=np.sqrt((Dbase+3*A*(q0*q0+w0*w0))/M)*dt
        if it%4==0 or it==len(times)-1:
            records.append((t*1e-15,excess,work_sum,loss_sum,float(actual_u.max()),float(state[3]),float(state[8]),
                theta_modes.copy(),q.copy(),mom.copy(),w*w,delta.copy(),phase.copy(),phase_ref.copy(),livevg.copy()))
            off_records.append((off_excess,M*(pwc*pwc/M**2+detc/(wc*wc))))
        if len(coarse_maps)<10 and t>=np.array([2,6,12,20,50,100,300,1000,3000,10000])[len(coarse_maps)]:
            coarse_maps.append((t*1e-15,actual_u.copy(),material.value(actual_u).copy(),material.value(actual_u,'P').copy()))
        previous_delta=delta
    retained=np.unique(np.r_[np.argmin(np.linalg.norm(k,axis=1)),np.linspace(0,len(k)-1,8).astype(int)])
    gamma_index=int(np.argmin(np.linalg.norm(k,axis=1)))
    saved_maps=coarse_maps if mesh==16 and ec==2.5 else []
    path=HERE/f"native_E{ec:g}_{arm['name']}_sp{sp}_m{mesh}.npz"
    np.savez_compressed(path,time_s=np.array([r[0] for r in records]),optical_excess_eV=np.array([r[1] for r in records]),
        native_work_eV=np.array([r[2] for r in records]),optical_loss_eV=np.array([r[3] for r in records]),
        maximum_U_eV_formula=np.array([r[4] for r in records]),hot_e_eV=np.array([r[5] for r in records]),thermal_e_excess_eV=np.array([r[6] for r in records]),
        velocity_theta_excess_eV=np.array([r[7]-theta for r in records]),
        variance_q_selected_A2=np.array([r[10][retained] for r in records]),projected_deltaK_eVA2=np.array([r[11][gamma_index] for r in records]),
        phase_selected_rad=np.array([r[12][retained] for r in records]),reference_phase_selected_rad=np.array([r[13][retained] for r in records]),
        group_velocity_selected_m_s=np.array([r[14][retained] for r in records]),selected_mode_indices=retained,k_m_inverse=k,base_K=Dbase,sites=sites,
        velocity_theta_no_e_excess_eV=np.array([r[1]-theta for r in off_records]),optical_no_e_coupling_eV=np.array([r[0] for r in off_records]),
        positions_r_z_nm=positions,quadrature_volume_nm3=volume,focus_weight=weights,
        map_time_s=np.array([r[0] for r in saved_maps]),map_U=np.array([r[1] for r in saved_maps]),
        map_deltaK=np.array([r[2] for r in saved_maps]),map_deltaP_Pa=np.array([r[3] for r in saved_maps]),
        actual_source_thermal_proxy=np.array('native table of total nonthermal+thermal excess is a constitutive hypothesis, not a measured electron temperature'))
    # Native pressure residence is a shot/cohort average of the nonlinear table,
    # not the equation of state applied to a dilute coarse mean.
    pressure_integral=np.trapezoid(pressures,times[1:]*1e-15)
    states=np.array(states);electron_integral=np.trapezoid(states[:,3]+states[:,8],times[1:]*1e-15)*s.QE
    optical_integral=np.trapezoid(states[:,4],times[1:]*1e-15)
    stress_integral=np.trapezoid(stress,times[1:]*1e-15)
    residence=np.trapezoid(states,times[1:]*1e-15,axis=0)
    baseline_residence=np.trapezoid(population[1:],times[1:]*1e-15,axis=0)
    off_states=np.array(off_states);off_residence=np.trapezoid(off_states,times[1:]*1e-15,axis=0)
    off_pressure=np.trapezoid(off_pressures,times[1:]*1e-15)
    off_optical=np.trapezoid(off_states[:,4],times[1:]*1e-15)
    return dict(arm=arm['name'],closure=arm['closure'],allocation=allocation,E_cm_keV=ec,species=sp+1,mesh=mesh,
        native_pressure_energy_coefficient=pressure_integral/electron_integral,native_pressure_residence_J_s=pressure_integral,
        maximum_U_eV_formula=float(u.max()),maximum_live_U_eV_formula=max(r[4] for r in records),old_predebit_reference_U_eV_formula=.03010 if ec==2.5 and sp else None,
        optical_native_work_eV=work_sum,maximum_packet_energy_defect_eV=maxdef,optical_sites=sites,
        maximum_complete_packet_defect_eV=maxglobal,carrier_residence_eV_s=residence.tolist(),
        carrier_residence_ratio=(residence/baseline_residence).tolist(),optical_strain_energy_ratio=stress_integral/optical_integral,
        thermal_energy_ratio_per_K=float(np.trapezoid(thermal_work,times[1:]*1e-15)/optical_integral),
        no_e_residence_ratio=(off_residence/baseline_residence).tolist(),no_e_carrier_residence_eV_s=off_residence.tolist(),
        no_e_optical_strain_energy_ratio=float(np.trapezoid(off_stress,times[1:]*1e-15)/off_optical),
        no_e_thermal_energy_ratio_per_K=float(np.trapezoid(off_thermal_work,times[1:]*1e-15)/off_optical),
        no_e_pressure_residence_J_s=float(off_pressure),
        dk_volume_residence_eVA2_m3_s=float(np.trapezoid(dk_volume,times[1:]*1e-15)),
        primary_quadrature_population_error_eV=quadrature_error,
        maximum_group_speed_m_s=float(np.linalg.norm(vg,axis=1).max()),bare_Gamma_K=Material(16).reference,cartesian_degrees_per_D=3,
        thermal_only_peak_U_eV_formula=float((population[:,8,None]*secondary_shape*V0*1e27).max()),
        primary_finite_front_source=geometry,bridge=bridge,path=path.name)


def prepare(cases=None,closures=None):
    """Build the single reusable source operator before whole-foil integration."""
    import run as model
    import fields
    import launch
    import inspect
    epoch=s.digest(HERE/'native.py')+s.digest(HERE/'launch.py')+s.digest(HERE/'fields.py')
    epoch+=inspect.getsource(model.network)+inspect.getsource(model.optical_histogram)+inspect.getsource(model.source_totals)
    epoch=__import__('hashlib').sha256(epoch.encode()).hexdigest()
    previous=json.loads((HERE/'native_summary.json').read_text()) if (HERE/'native_summary.json').exists() else {}
    valid=previous.get('epoch')==epoch
    rows=previous.get('rows',[]) if valid else [];source_rows=previous.get('sources',[]) if valid else []
    calibrations={}
    if valid and (HERE/'frozen/native_transfer.npz').exists():
        saved=np.load(HERE/'frozen/native_transfer.npz')['normalization'];calibrations={j:saved[j] for j in (0,1)}
    oldarms={a['name']:a['source_debits'] for a in previous.get('arms',[])} if valid else {}
    tracks=json.loads((HERE/'frozen/tracks.json').read_text())
    for physical in s.PHYSICAL_ARMS:
        if closures is not None and physical['name'] not in closures:continue
        variants=[a for a in s.ARMS if a['closure']==physical['name']]
        for arm in variants:
            arm['source_debits']=oldarms.get(arm['name'],{})
            model.carrier_kernel(arm)
        endpoints=np.load(HERE/f"paths_{variants[0]['name']}.npz")['heat_endpoints']
        for ec in reversed(s.ENERGIES):
            if cases is not None and ec not in cases:continue
            for sp in (0,1):
                if all(any(r['arm']==a['name'] and r['E_cm_keV']==ec and r['species']==sp+1 and r['mesh']==mesh for r in rows) for a in variants for mesh in (12,16)):continue
                tr=tracks[f'{ec:g}_{sp}'];basis=source_basis(variants[0],ec,sp,tr)
                spectrum,velocities,polarizations=launch.spectral_source(variants[0],ec,sp,tr,basis,endpoints)
                if sp not in calibrations:calibrations[sp]=launch.calibration(spectrum,sp)
                energy,wm,wh,we,defect=launch.project(spectrum,calibrations[sp])
                mechanical=float(wm.sum()/s.QE);electronic=float(we.sum()/s.QE)
                key=f'{ec:g}_{sp}'
                acoustic_path=f"source_E{ec:g}_{physical['name']}_sp{sp}.npz"
                np.savez_compressed(HERE/acoustic_path,energy_J=energy,mechanical_work_J=wm,heat_work_J=wh,electron_work_J=we,
                    source_amplitudes=spectrum,source_component_names=np.array(['mechanical_mean','sqrt_mechanical_variance','causal_lattice_heat','native_electron_pressure','frozen_instant_heat_reference']),
                    polarization=polarizations,group_speed_m_s=velocities,directions=launch.DIRECTIONS)
                raw_electron=(sp+1)*tr['components_eV'][0]
                source_rows.append(dict(closure=physical['name'],E_cm_keV=ec,species=sp+1,
                    acoustic_eV=float(energy.sum()/s.QE),mechanical_work_eV=mechanical,electron_work_eV=electronic,
                    heat_work_eV=float(wh.sum()/s.QE),maximum_signed_work_identity_defect_J=defect,
                    pressure_operator_split_relative_electron_energy=abs(electronic)/raw_electron,
                    acoustic_path=acoustic_path))
                for arm in variants:
                    arm['source_debits'][key]=dict(host=mechanical,electron=electronic)
                    # The radiation-work port debits the same parent before
                    # subsequent carrier propagation; spatial Green shape stays
                    # retarded, and its normalization is not a new source.
                    pos,vol,ww,times,geo,tt,hh,length,pop,secondary,error,u=basis
                    nt,nh,nlength,_=model.optical_histogram(tr,sp,arm,times*1e-15)
                    npop=np.column_stack([np.interp(times,nt*1e15,nh[:,j]) for j in range(9)])
                    old=pop[:,3]+pop[:,8];new=npop[:,3]+npop[:,8]
                    scale=np.divide(new,old,out=np.zeros_like(new),where=old>0)
                    adjusted=(pos,vol,ww,times,geo,nt,nh,nlength,npop,secondary,error,u*scale[:,None])
                    for mesh in (16,12):
                        result=packet(arm,ec,sp,tr,mesh,adjusted);rows.append(result)
                print(json.dumps(dict(native_source_complete=True,closure=physical['name'],E_cm_keV=ec,species=sp+1,
                    acoustic_eV=float(energy.sum()/s.QE),maximum_live_U_eV_formula=rows[-2]['maximum_live_U_eV_formula'])),flush=True)
    np.savez_compressed(HERE/'frozen/native_transfer.npz',normalization=np.stack([calibrations[j] for j in (0,1)]),
        reference=np.array('accepted G13 track band integrals; six-direction retained source projection, no fitted beam sound fraction'))
    model.table('native_source_work.csv',source_rows)
    (HERE/'native_summary.json').write_text(json.dumps(dict(epoch=epoch,rows=rows,sources=source_rows,
        arms=[dict(name=a['name'],source_debits=a.get('source_debits',oldarms.get(a['name'],{}))) for a in s.ARMS]),indent=2)+'\n')
    s.refresh_optical();load.cache_clear()
    return rows


@lru_cache(maxsize=1)
def load():
    data=json.loads((HERE/'native_summary.json').read_text())
    by_name={a['name']:a['source_debits'] for a in data['arms']}
    for arm in s.ARMS:arm['source_debits']=by_name[arm['name']]
    return data


def constitutive(arm,ec,mesh=16,control='joined'):
    data=load();rr=[r for r in data['rows'] if r['arm']==arm['name'] and r['E_cm_keV']==ec and r['mesh']==mesh]
    if len(rr)!=2:raise ValueError('Missing native species constitutive operator')
    rr.sort(key=lambda r:r['species']);fractions=np.array([.85,.15]);off=control=='native_e_off'
    residence_key='no_e_carrier_residence_eV_s' if off else 'carrier_residence_eV_s'
    ratio_key='no_e_residence_ratio' if off else 'carrier_residence_ratio'
    pressure_key='no_e_pressure_residence_J_s' if off else 'native_pressure_residence_J_s'
    stress_key='no_e_optical_strain_energy_ratio' if off else 'optical_strain_energy_ratio'
    thermal_key='no_e_thermal_energy_ratio_per_K' if off else 'thermal_energy_ratio_per_K'
    residence=np.array([r[residence_key] for r in rr]);ratio=np.array([r[ratio_key] for r in rr])
    raw=residence/ratio;aggregate=fractions@residence
    electron=float(aggregate[3]+aggregate[8])*s.QE
    pressure=sum(fr*r[pressure_key] for fr,r in zip(fractions,rr))
    stress=sum(fr*r[stress_key]*r[residence_key][4] for fr,r in zip(fractions,rr))
    thermal=sum(fr*r[thermal_key]*r[residence_key][4] for fr,r in zip(fractions,rr))
    dk=0. if off else sum(fr*r['dk_volume_residence_eVA2_m3_s'] for fr,r in zip(fractions,rr))
    return dict(pressure_ratio=pressure/electron,optical_strain_ratio=stress/aggregate[4],thermal_energy_ratio_per_K=thermal/aggregate[4],
        residence_ratio=aggregate/(fractions@raw),species_residence_ratio=ratio,dk_per_energy_density=dk/electron,
        pressure_residence_species_J_s=np.array([r[pressure_key] for r in rr]))
