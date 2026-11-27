"""Full-window electronic-gradient projection of the accepted bulk linear Green.
Azimuth is integrated analytically (J1 for signed radial force, J0 for axial).
The stripe result is a continuum loading/ray average, NOT an FCC registration
calculation. It does not replace the nonlinear motif energy ledger.
"""
import math
import numpy as np
from scipy.special import roots_legendre,j0,j1,jv
import whole


def rule(edges,order=8):
    x,w=roots_legendre(order);edges=np.asarray(edges)
    return (((edges[1:]+edges[:-1])[:,None]/2+np.diff(edges)[:,None]*x/2).ravel(),
        (np.diff(edges)[:,None]*w/2).ravel())


def radial(r,R):
    x=r/R;inside=x<1;ex=np.exp(-2*x*x)
    g=np.where(inside,(1-x*x)**2*ex,0.)
    d=np.where(inside,-4*x*(1-x*x)*(2-x*x)*ex/R,0.)
    return g,d


def union_area(R,d):
    return 2*np.pi*R*R-(2*R*R*np.arccos(d/(2*R))-.5*d*np.sqrt(4*R*R-d*d)) if d else np.pi*R*R


def union_transform(q,R,d):
    kx,ky=q[:2]
    if d==0:
        k=np.hypot(kx,ky);return np.pi*R*R*(2*j1(k*R)/(k*R) if k else 1.)
    stop=np.arcsin(d/(2*R));n=max(16,int(np.ceil((stop+np.pi/2)*R*max(abs(kx),abs(ky))*2/np.pi)))
    theta,w=rule(np.linspace(-np.pi/2,stop,n+1))
    x=R*np.sin(theta);y=R*np.cos(theta)
    a=np.sum(w*2*R*R*np.cos(theta)**2*np.sinc(ky*y/np.pi)*np.exp(-1j*kx*x))
    return a+np.exp(-1j*kx*d)*a.conjugate()


def axial(m,ec,sp,s):
    R=s['radius']*1e10;z,ts,mid,logu,normal=m.electron_profile(ec,sp,R)
    edges=np.unique(np.r_[z,mid,[r[1]*1e10 for r in m.p1.PROFILES[ec]['layers']]])
    edges=edges[(edges>=0)&(edges<=z[-1])]
    zz,w=rule(edges)
    U,grad,clock,gclock=m.electronic_coordinate(ec,sp,s,np.c_[zz*0,zz*0,zz])
    return zz,w*m.profile_nd(ec,zz)*1e-30,U,grad[:,2],clock,gclock[:,2],z[-1],ts[-1]


def moments(m,ec,sp,s,q,omega,geometry,order=8):
    """Finite Green pole/particular source integrals and free-ring q-squared tensor."""
    R=s['radius']*1e10;qr=np.linalg.norm(q[:2]);direction=q[:2]/qr if qr else np.array([1.,0.])
    transverse=np.array([-direction[1],direction[0]])
    z,wz,U,Uz,clock,clock_z,zmax,tmax=axial(m,ec,sp,s)
    kappa=1/(s['tau_opt']*1e15);beta=kappa;rho=1/(s['tau_ne']*1e15);alpha=kappa/2
    wd=math.sqrt(omega*omega-alpha*alpha);lam=alpha+1j*wd
    delay=.74/math.sqrt(2*m.p1.track(ec,sp)['initial_eV']/m.M) if sp and geometry=='axial' else 0.
    displacement=.74 if sp and geometry=='transverse' else 0.
    tref=tmax+R/10+delay
    n=max(24,int(np.ceil(R*max(qr,1/R)*2/np.pi)))
    r,wr=rule(np.linspace(0,R,n+1),order);g,gp=radial(r,R)
    N=union_area(R,displacement)*wz.sum()
    circle_N=np.pi*R*R*wz.sum()
    def source(l):
        radial0=np.sum(wr*r*j0(qr*r)*g*np.exp(l*r/10))
        radial1=np.sum(wr*r*j1(qr*r)*(gp+l*g/10)*np.exp(l*r/10))
        zz=np.exp(l*(clock-tref)-1j*q[2]*z)
        a0=np.sum(wz*U*zz);a1=np.sum(wz*(Uz+l*U*clock_z)*zz)
        out=np.r_[-2j*np.pi*direction*radial1*a0,2*np.pi*radial0*a1]
        if sp:out*=1+np.exp(l*delay-1j*q[0]*displacement)
        return out/N
    rates=[beta,beta+rho,lam,lam.conjugate()]
    I=np.array([source(l) for l in rates])
    # The tensor uses a single physically chosen pump pole (passed omega), not
    # mean phase or |sum A|². Its homogeneous part excludes the drifting offset.
    rr=(gp+lam*g/10)*np.exp(lam*r/10);rz=g*np.exp(lam*r/10)
    Z0=U*np.exp(lam*(clock-tref));Z1=(Uz+lam*U*clock_z)*np.exp(lam*(clock-tref))
    z00=np.sum(wz*Z0*Z0*np.exp(-2j*q[2]*z));z11=np.sum(wz*Z1*Z1*np.exp(-2j*q[2]*z))
    z01=np.sum(wz*Z0*Z1*np.exp(-2j*q[2]*z))
    tensor=np.zeros((3,3),complex)
    if not displacement:
        x=2*qr*r
        Rxy=np.pi*(np.outer(direction,direction)[:,:,None]*(j0(x)-jv(2,x))+
            np.outer(transverse,transverse)[:,:,None]*(j0(x)+jv(2,x)))
        tensor[:2,:2]=np.sum(Rxy*(wr*r*rr*rr)[None,None,:],axis=2)*z00
        tensor[2,2]=2*np.pi*np.sum(wr*r*j0(x)*rz*rz)*z11
        tensor[:2,2]=-2j*np.pi*direction*np.sum(wr*r*j1(x)*rr*rz)*z01
        tensor[2,:2]=tensor[:2,2]
        J2=2*np.pi*(np.sum(wr*r*abs(rr)**2)*np.sum(wz*abs(Z0)**2)+
            np.sum(wr*r*abs(rz)**2)*np.sum(wz*abs(Z1)**2))
        if sp:twin=1+np.exp(lam*delay);tensor*=twin*twin;J2*=abs(twin)**2
    else:
        # Integrate the smooth angular source interpolant, never the unresolved
        # exp(-2 i q.r) on a four-angle motif rule. Exact Bessel multipliers retain
        # every angular Fourier coefficient of this declared 32-angle interpolant.
        r,wr=rule(np.linspace(0,R+displacement,n+2),order)
        phi=2*np.pi*np.arange(32)/32;orders=np.fft.fftfreq(32,1/32).astype(int)
        xx=r[:,None]*np.cos(phi);yy=r[:,None]*np.sin(phi)
        sr=np.zeros((len(r),32,2),complex);sz=np.zeros((len(r),32),complex)
        for offset in (0.,displacement):
            dx=xx-offset;distance=np.hypot(dx,yy);gg,dd=radial(distance,R)
            amplitude=(dd+lam*gg/10)*np.exp(lam*distance/10)
            denom=np.maximum(distance,1e-100)
            sr[:,:,0]+=amplitude*dx/denom;sr[:,:,1]+=amplitude*yy/denom
            sz+=gg*np.exp(lam*distance/10)
        theta=np.arctan2(q[1],q[0]) if qr else 0.
        multipliers=((-1j)**orders)[None,:]*jv(orders[None,:],2*qr*r[:,None])*np.exp(1j*orders*theta)[None,:]
        def angular(value):return 2*np.pi*np.sum(np.fft.fft(value,axis=1)/32*multipliers,axis=1)
        for a in range(2):
            for b in range(2):tensor[a,b]=np.sum(wr*r*angular(sr[:,:,a]*sr[:,:,b]))*z00
            tensor[a,2]=np.sum(wr*r*angular(sr[:,:,a]*sz))*z01;tensor[2,a]=tensor[a,2]
        tensor[2,2]=np.sum(wr*r*angular(sz*sz))*z11
        J2=2*np.pi*(np.sum(wr*r*np.mean(np.sum(abs(sr)**2,axis=2),axis=1))*np.sum(wz*abs(Z0)**2)+
            np.sum(wr*r*np.mean(abs(sz)**2,axis=1))*np.sum(wz*abs(Z1)**2))
    L=rho/((beta-lam)*(beta+rho-lam));scale=1j*L/(m.M*wd)
    tensor*=scale*scale;J2*=abs(scale)**2
    gram_z=np.sum(wz*np.exp(-2j*q[2]*z))/wz.sum()
    gram=gram_z*union_transform(2*q,R,displacement)/union_area(R,displacement)
    G=.5*np.array([[1+gram.real,-gram.imag],[-gram.imag,1-gram.real]])
    inverse=np.linalg.pinv(G,rcond=1e-12)
    return dict(I=I,omega=omega,wd=wd,alpha=alpha,beta=beta,rho=rho,tref=tref,N=N,
        tensor=tensor,amplitude_norm_D_A2=float(J2),inverse_gram=inverse,gram=gram,
        radial_nodes=len(r),axial_nodes=len(z),circle_N=circle_N)


def states(m,d,delta):
    b=d['beta'];rho=d['rho'];a=d['alpha'];wd=d['wd'];I=d['I'];lam=a+1j*wd
    D1=(b-a)**2+wd*wd;D2=(b+rho-a)**2+wd*wd
    p1=I[0]*np.exp(-b*delta)/(m.M*D1);p2=-I[1]*np.exp(-(b+rho)*delta)/(m.M*D2)
    L=rho/((b-lam)*(b+rho-lam))
    h1=.5j*L*I[2]*np.exp(-lam*delta)/(m.M*wd)
    h2=-.5j*L.conjugate()*I[3]*np.exp(-lam.conjugate()*delta)/(m.M*wd)
    Q=p1+p2+h1+h2;V=-b*p1-(b+rho)*p2-lam*h1-lam.conjugate()*h2
    return Q,V,h1+h2


def execute(m,summary):
    rows=[];checks=[];critical=[]
    source={(r['E_cm_keV'],r['species']-1,r['scenario']):r for r in summary if r['environment']=='bulk' and r['branch']=='absent'}
    for ec in m.EC:
        for sp in (0,1):
            for s in m.p1.SCENARIOS:
                names,waves,frequencies=whole.mode_spec(m,ec,sp,s)
                for geometry in (['single'] if sp==0 else ['axial','transverse']):
                    for name,q,omega in zip(names,waves,frequencies):
                        d=moments(m,ec,sp,s,q,omega,geometry)
                        # h source is the near-Gamma electronic field, with the
                        # actual target frequency retained in its denominator.
                        pump=d if name=='Gamma' else moments(m,ec,sp,s,q,frequencies[0],geometry)
                        T=pump['tensor'];H=m.B*np.trace(T)*np.eye(3)+2*m.B*T+3*(m.A-m.B)*np.diag(np.diag(T))
                        h=float(np.linalg.svd(H,compute_uv=False)[0]/(2*m.M*pump['N']*omega**2))
                        kappa=2*d['alpha'];detuning=pump['wd']-omega
                        threshold=4*np.sqrt(detuning*detuning+(kappa/2)**2)/omega
                        budget=source[ec,sp,s['name']]['original_optical_positive_flow_eV']
                        critical.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],geometry=geometry,mode=name,
                            continuum_free_ring_h_per_chi_squared=h,pump_frequency_THz=pump['wd']*1000/(2*np.pi),
                            target_frequency_THz=omega*1000/(2*np.pi),detuning_fs_inverse=detuning,
                            frozen_detuned_threshold=threshold,linear_chi_at_frozen_threshold=np.sqrt(threshold/h) if h>0 else None,
                            source_energy_ceiling_eV=budget,meaning='linear free-ring extrapolation, not a nonlinear gain/yield fit or measured coefficient'))
                        for multiple in (0.,.25,.5,1.,2.,4.,8.):
                            delta=multiple/kappa;Q,V,ring=states(m,d,delta)
                            pQ=np.array([Q.real,-Q.imag]);pV=np.array([V.real,-V.imag])
                            E=float(.5*m.M*d['N']*(np.sum(pV*(d['inverse_gram']@pV))+omega**2*np.sum(pQ*(d['inverse_gram']@pQ))))
                            for branch,chi in m.BRANCHES:
                                rows.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment='bulk_linear_reference',
                                    geometry=geometry,branch=branch,mode=name,time_fs=d['tref']+delta,time_after_last_local_onset_fs=delta,
                                    qx_A_inverse=q[0],qy_A_inverse=q[1],qz_A_inverse=q[2],physical_occupied_D=d['N'],
                                    optical_frequency_THz=omega*1000/(2*np.pi),window_real_mode_probe_energy_eV=chi*chi*E,
                                    projected_normal_mode_number=chi*chi*E/(m.p1.HB*1e15*omega),
                                    source_energy_ceiling_eV=budget,linear_energy_below_source_ceiling=chi*chi*E<=budget,
                                    free_ring_2wd_stiffness_h=chi*chi*h*np.exp(-kappa*delta),
                                    curvature_model='long-wavelength Gamma-pole electronic field; nonlinear motif curvature fits are separate',
                                    electronic_component_only=True,
                                    mechanical_energy_ceiling_eV=source[ec,sp,s['name']]['original_direct_D_eV'],
                                    particular_coordinate_residual_A=float(abs(chi*(Q-ring)).max()),
                                    radial_nodes=d['radial_nodes'],axial_nodes=d['axial_nodes'],azimuth='analytic signed J1 radial / J0 axial',
                                    source='same energy-derived retarded displacive coordinate; beta=kappa, rho=1/tau_ne; no extra energy debit',
                                    support='post last geometric onset; linear-demand probe, not the capped nonlinear trajectory or a sum over independent modes',
                                    registration='continuum loading/ray-average projection; FCC registration remainder not fixed by profiles' if 'sheet' in name else 'continuous loading symmetry',
                                    remainder_energy_ceiling_eV=budget if 'sheet' in name else 0.))
                        if ec==2.5 and sp==0 and s==m.p1.SCENARIOS[-1] and name=='opposite_sheet_110':
                            R=s['radius']*1e10;k=np.linalg.norm(q[:2]);r,w=rule(np.linspace(0,R,max(24,int(np.ceil(2*k*R/np.pi)))+1))
                            actual=np.sum(w*r*np.pi*(1+j0(2*k*r)))
                            exact=.5*np.pi*R*R*(1+2*j1(2*k*R)/(2*k*R))
                            gain=actual/exact
                            radial_actual=-2j*np.pi*np.sum(w*r*j1(k*r)**2)
                            radial_exact=-1j*np.pi*R*R*(j1(k*R)**2-j0(k*R)*jv(2,k*R))
                            radial_gain=radial_actual/radial_exact
                            if max(abs(gain-1),abs(radial_gain-1))>1e-8:raise ValueError('Full-window mode-matched mass/source projection failed')
                            checks.append(dict(mode='accepted stripe: cosine axial source and matched radial J1 vector source',
                                mode_matched_gain=gain,radial_signed_mode_gain=float(radial_gain.real),
                                relative_difference=max(abs(gain-1),abs(radial_gain-1)),Gamma_radial_parity=0.,
                                normalization='same full disk, physical occupied column and real-mode mass metric; radial projection is negative imaginary',
                                limitation='continuum angular/source integration check, not an FCC registration certificate'))
    m.table('whole_patch_modes.csv',rows);m.table('critical_phase_drive.csv',critical)
    m.dump('full_window_checks.json',checks)
    return rows
