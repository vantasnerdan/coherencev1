"""G15 Pass 2. Run from the repository root: python g15_track/pass2/run.py.
Conditional signed forces, conservative optical work, physical FCC microblocks,
finite native pair cuts, and the same fixed-profile whole-foil consumers.
"""
from __future__ import annotations
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
import sys
sys.dont_write_bytecode=True
from functools import lru_cache
import argparse,csv,hashlib,importlib.util,json,math,time
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.interpolate import PchipInterpolator
from scipy.special import roots_legendre
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
P1=ROOT/'g15_track/pass1'
sys.path.insert(0,str(P1))
spec=importlib.util.spec_from_file_location('accepted_pass1',P1/'run.py')
p1=importlib.util.module_from_spec(spec);spec.loader.exec_module(p1)
spec=importlib.util.spec_from_file_location('accepted_lattice',ROOT/'g3_lattice/coherence/lattice.py')
g3=importlib.util.module_from_spec(spec);spec.loader.exec_module(g3)
EV=p1.QE;M=p1.MD*1e10/EV  # eV fs^2 / A^2
P=p1.DATA['lattice'];K,A,B,_=g3.environment_law(P,0,0.)
OMEGA=math.sqrt(K/M);KG=K-sum(s['J_sum_eV_A2'] for s in P['shells'])
EPOCH='Pass1 frozen lattice snapshot: provisional coherence/lattice_inputs.json epoch, not corrected rerun'
EC=[.25,.5,1.,1.5,2.,2.5]
BRANCHES=[('absent',0.),('gradient_plus',1.),('gradient_minus',-1.),
          ('gradient_large_plus',1e6),('gradient_large_minus',-1e6)]
SHAPE=(2,2,2)
JREAL=g3.shell_kernel(P,SHAPE);JFFT=np.fft.rfftn(JREAL)
# Sample the authoritative convolution once; eight-site dense action avoids tiny FFTs.
_basis=np.eye(int(np.prod(SHAPE))).reshape((int(np.prod(SHAPE)),)+SHAPE)
JMATRIX=np.ascontiguousarray(g3.coupling_force(_basis,JFFT,SHAPE).reshape(int(np.prod(SHAPE)),-1))

def neighbor_force(q):
    return (q.reshape(-1,int(np.prod(SHAPE)))@JMATRIX).reshape(q.shape)
VAC=json.loads((ROOT/'g10_quantum_site/pd-vacancy-r4-six-D-pair/pair.json').read_text())
VREF=json.loads((ROOT/'g10_quantum_site/pd-vacancy-r4-six-D-pair/native/reference.json').read_text())
BULK=json.loads((ROOT/'g10_quantum_site/pdd-nearest-pair/pair-wells.json').read_text())
AXIS={'bulk':np.array(BULK['closing_axis']),'vacancy':np.array(VAC['axis_A_to_B'])}
D0={'bulk':2.88782409436586,'vacancy':VAC['rest_separation_A']}
CUTS={env: p1.DATA['bulk_pair'] if env=='bulk' else p1.DATA['vacancy_pair']['rows'] for env in D0}

def dump(name,value):
    def conv(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    (HERE/name).write_text(json.dumps(value,indent=2,allow_nan=False,default=conv)+'\n')

def table(name,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def raw_table(name):
    with (P1/name).open() as f:return list(csv.DictReader(f))

def work_fraction(linear,quad,available):
    scale=np.ones_like(linear);over=linear+quad>available
    root=np.sqrt(linear*linear+4*quad*np.maximum(available,0))
    positive=over&(linear>=0)
    negative=over&(linear<0)
    scale[positive]=2*available[positive]/np.maximum(linear[positive]+root[positive],1e-100)
    scale[negative]=(-linear[negative]+root[negative])/np.maximum(2*quad[negative],1e-100)
    return np.clip(scale,0,1)

def stock_path(ec,sp,s,dt_fs=1.,extra=(),first_only=False):
    """Exact constant-source Pass1 ODE; sample, do not rerun its thermal producer."""
    tr=p1.track(ec,sp);ts,z,dep,src=p1.source(tr,sp);ts*=1e15
    delay=.74/math.sqrt(2*tr['initial_eV']/M) if sp else 0.
    stop=ts[-1]+8*max(s['tau_opt'],s['tau_ne'],s['tau_end'])*1e15
    edges=np.unique(np.r_[np.arange(0,stop+dt_fs,dt_fs),ts,ts+delay,np.asarray(extra)])
    edges=np.unique(np.r_[0.,edges[(edges>=0)&(edges<=stop)],stop])
    matrix=p1.network(s)*1e-15
    cache={};y=np.zeros(11);states=[];flows=[];times=[];inputs=[];integrals=[];keys=[];injected=0.
    offsets=[0.,delay] if sp and not first_only else [0.]
    for ta,tb in zip(edges[:-1],edges[1:]):
        dt=float(tb-ta);key=round(dt,11)
        if key not in cache:
            ext=np.zeros((33,33));ext[:11,:11]=matrix;ext[:11,11:22]=np.eye(11);ext[22:,:11]=np.eye(11)
            e=expm(ext*dt);cache[key]=(e[:11,:11],e[:11,11:22],e[22:,11:22])
        e,f,g=cache[key];u=np.zeros(11)
        for off in offsets:
            j=np.searchsorted(ts,(ta+tb)/2-off,side='right')-1
            if 0<=j<len(src):u+=src[j]/((sp+1)*(ts[j+1]-ts[j]))
        new=e@y+f@u
        integral=f@y+g@u
        inflow=float((matrix[4].copy()+np.eye(11)[4]*(-matrix[4,4]))@integral)
        for off in offsets:
            if abs(tb-ts[-1]-off)<1e-7:new[0]+=tr['terminal_unresolved_eV'];injected+=tr['terminal_unresolved_eV']
        injected+=float(u.sum()*dt)
        states.append(new.copy());flows.append(max(0.,inflow));times.append(tb);inputs.append(injected)
        integrals.append(integral);keys.append(key);y=new
    return dict(t=np.r_[0.,times],y=np.vstack([np.zeros(11),states]),flow=np.r_[0.,flows],
                tr=tr,z=z*1e10,ts=ts,dep=dep,total=injected,delay=delay,inputs=np.r_[0.,inputs],
                integral=np.vstack([np.zeros(11),integrals]),propagator=cache,step_keys=keys)

def profile_nd(ec,z_A):
    layers=p1.PROFILES[ec]['layers'];ends=np.array([x[1]*1e10 for x in layers])
    idx=np.minimum(np.searchsorted(ends,z_A),len(ends)-1)
    return np.array([x[3] for x in layers])[idx]

@lru_cache(maxsize=None)
def electron_profile(ec,sp,R):
    tr=p1.track(ec,sp);ts,z,dep,_=p1.source(tr,sp);z=z*1e10;ts=ts*1e15
    mid=(z[:-1]+z[1:])/2;nd=profile_nd(ec,mid)*1e-30
    atomic_volume=P['a_A']**3/4
    values=dep[:,0]*atomic_volume/((sp+1)*np.pi*R*R*np.diff(z))
    logu=np.log(np.maximum(values,1e-100))
    breaks=np.unique(np.r_[z,mid,[x[1]*1e10 for x in p1.PROFILES[ec]['layers']]])
    breaks=breaks[(breaks>=0)&(breaks<=z[-1])]
    gl,gw=roots_legendre(8);zz=(breaks[:-1,None]+breaks[1:,None])/2+np.diff(breaks)[:,None]*gl/2
    integral=float(np.sum(np.diff(breaks)[:,None]*gw/2*np.exp(np.interp(zz,mid,logu))))
    radial=float(np.sum(gw/2*(1-(gl+1)/2)**2*np.exp(-2*(gl+1)/2)))
    normal=(dep[:,0].sum()*atomic_volume/(sp+1))/(np.pi*R*R*integral*radial)
    return z,ts,mid,logu,normal

def electronic_coordinate(ec,sp,s,pos):
    R=s['radius']*1e10;z,ts,mid,logu,normal=electron_profile(ec,sp,R)
    zz=pos[:,2];rr=np.sum(pos[:,:2]**2,axis=1);taper=np.maximum(0.,1-rr/R**2)
    u=normal*np.exp(np.interp(zz,mid,logu))*np.exp(-2*rr/R**2)*taper*taper
    interval=np.clip(np.searchsorted(mid,zz,side='right')-1,0,max(0,len(mid)-2))
    slope=np.diff(logu)/np.diff(mid)
    dz=slope[interval] if len(slope) else np.zeros(len(zz))
    dz=np.where((zz>mid[0])&(zz<mid[-1]),dz,0.)
    radial=-4/R**2-4/np.maximum(R*R-rr,1e-100)
    grad=np.c_[u*radial*pos[:,0],u*radial*pos[:,1],u*dz]
    j=np.clip(np.searchsorted(z,zz,side='right')-1,0,len(z)-2)
    dt_dz=np.diff(ts)[j]/np.diff(z)[j]
    r=np.sqrt(rr);clock=np.interp(zz,z,ts)+r/10.
    gradclock=np.c_[pos[:,0]/np.maximum(10*r,1e-100),pos[:,1]/np.maximum(10*r,1e-100),dt_dz]
    return u,grad,clock,gradclock

def map_nodes(ec,sp,s,fine=False):
    """Annular/azimuth/depth quadrature; weights count physical occupied D."""
    tr=p1.track(ec,sp);ts,z,dep,_=p1.source(tr,sp);z*=1e10;ts*=1e15
    R=s['radius']*1e10
    re=np.unique(np.r_[0.,.12,.3,.65,1.4,3.,8.,R/4,R/2,R]) if fine else np.array([0.,.8,3.,R/2,R])
    re=np.sort(re[(re>=0)&(re<=R)])
    nz=12 if fine else 2;nphi=8 if fine else 4
    ze=np.linspace(0,z[-1],nz+1);zz=(ze[1:]+ze[:-1])/2
    rr=np.sqrt((re[1:]**2+re[:-1]**2)/2);phi=(np.arange(nphi)+.5)*2*np.pi/nphi
    r,ph,za=np.meshgrid(rr,phi,zz,indexing='ij');x=r*np.cos(ph);y=r*np.sin(ph)
    pos=np.array([x.ravel(),y.ravel(),za.ravel()]).T
    nd=profile_nd(ec,pos[:,2]);ann=np.repeat(np.pi*np.diff(re**2),nphi*nz)/nphi
    weights=ann*np.tile(np.diff(ze),len(rr)*nphi)*nd*1e-30
    N=p1.nd_column(ec,z[-1]*1e-10)*np.pi*s['radius']**2
    weights*=N/weights.sum()
    ee,gradient,t0,gradient_clock=electronic_coordinate(ec,sp,s,pos)
    en=np.interp(pos[:,2],z,np.r_[tr['initial_eV'],[v['energy_eV'] for v in tr['timeline']]])
    return pos,weights,ee,t0,en,N

def force_signature(pos,en,omega):
    """ZBL derivative integral with signed force and finite transit clock.
    CM-angle equal-mass recoil limits the close-collision impulse to E_projectile.
    Raw close kicks are exported and classified, never removed from display.
    """
    r=np.linalg.norm(pos[:,:2],axis=1);v=np.sqrt(2*en/M)
    gl,gw=roots_legendre(64);aeta=np.arcsinh(24.);u=np.sinh(gl*aeta);wt=gw*aeta*np.cosh(gl*aeta)
    b=np.maximum(r,.003);s=u[None,:]*b[:,None]
    dist=np.sqrt(b[:,None]**2+s*s);a=.8854*.529177/2
    F=14.39964*(p1.bca.zbl_phi(dist/a)/dist**2-p1.bca.zbl_dphi(dist/a)/(a*dist))
    dt=b[:,None]/v[:,None];osc=np.exp(1j*omega*s/v[:,None])
    pr=np.sum(F*b[:,None]/dist*osc*wt*dt,axis=1)
    pz=np.sum(F*(-s)/dist*osc*wt*dt,axis=1)
    # Derive reduced energy from the actual Coulomb screening length instead.
    eps=(en/2)*a/14.39964
    theta=p1.bca.cm_angle(eps,b/a)
    T=en*np.sin(theta/2)**2
    raw=np.c_[pr*pos[:,0]/b,pr*pos[:,1]/b,pz]
    unphased=np.sum(F*b[:,None]/dist*wt*dt,axis=1)
    scale=np.minimum(1.,np.sqrt(2*M*T)/np.maximum(abs(unphased),1e-100))
    # Retain the signed axial odd force, not an invented common recoil phase.
    return 1j*raw*scale[:,None]/(M*omega),T,raw,theta

def response_map(ec,sp,s,env,chi,geom='axial',fine=True):
    pos,w,ee,t0,en,N=map_nodes(ec,sp,s,fine)
    if env=='vacancy' and fine:
        native=np.array(VREF['positions_A'])[63:69];offset=native-native.mean(0)
        pos=(pos[:,None,:]+offset[None,:,:]).reshape(-1,3)
        w=np.repeat(w/6,6);en=np.repeat(en,6)
    z=pos[:,2];tr=p1.track(ec,sp)
    mech,T,raw,theta=force_signature(pos,en,OMEGA)
    qcap=.36 if env=='bulk' else .24
    cap=.5*K*qcap*qcap+.25*A*qcap**4
    strong=T>cap
    recoil=sum(p1.source(tr,sp)[2][:,2])/(sp+1)
    norm=min(1.,recoil/max(float(w@T),1e-100))
    bound=mech.copy()*math.sqrt(norm);bound[strong]=0
    ee,grad,t0,gradclock=electronic_coordinate(ec,sp,s,pos)
    F=chi*grad;G=-chi*ee[:,None]*gradclock
    if env=='vacancy':
        coords=np.array(VREF['positions_A'])[63:69];n=coords-coords.mean(0)
        n/=np.linalg.norm(n,axis=1)[:,None]
        F+=chi*ee[:,None]*n[np.arange(len(z))%6]/.1
    rho=1/(s['tau_ne']*1e15);beta=1/(s['tau_opt']*1e15)
    integral=(F-1j*OMEGA*G)*rho/((beta-1j*OMEGA)*(beta+rho-1j*OMEGA))
    elec=1j*integral/(M*OMEGA)
    mech_time=np.interp(pos[:,2],np.r_[0.,[r['depth_A'] for r in tr['timeline']]],np.r_[0.,[r['time_fs'] for r in tr['timeline']]])
    phase=np.exp(1j*OMEGA*t0)[:,None]
    mechanical_amp=bound*np.exp(1j*OMEGA*mech_time)[:,None]
    electronic_amp=elec*phase
    T2=np.zeros_like(T);strong2=np.zeros_like(strong);raw2=np.zeros_like(raw)
    if sp:
        delay=.74/math.sqrt(2*tr['initial_eV']/M) if geom=='axial' else 0.
        shift=np.array([.74,0.,0.]) if geom=='transverse' else np.zeros(3)
        second=pos.copy();second-=shift
        b2,T2,raw2,_=force_signature(second,en,OMEGA);strong2=T2>cap
        norm2=min(1.,recoil/max(float(w@T2),1e-100))
        b2*=math.sqrt(norm2);b2[strong2]=0
        ee2,grad2,t02,gclock2=electronic_coordinate(ec,sp,s,second)
        F2=chi*grad2;G2=-chi*ee2[:,None]*gclock2
        if env=='vacancy':F2+=chi*ee2[:,None]*n[np.arange(len(z))%6]/.1
        elec2=1j*(F2-1j*OMEGA*G2)*rho/((beta-1j*OMEGA)*(beta+rho-1j*OMEGA))/(M*OMEGA)
        mechanical_amp+=b2*np.exp(1j*OMEGA*(mech_time+delay))[:,None]
        electronic_amp+=elec2*np.exp(1j*OMEGA*(t02+delay))[:,None]
    amp=mechanical_amp+electronic_amp
    total_req=float(.5*K*np.sum(w[:,None]*abs(amp)**2))
    return dict(pos=pos,w=w,ee=ee,t0=t0,en=en,N=N,mech=bound,F=F,G=G,amp=amp,T=T,
                t_mech=mech_time,strong=strong,cap=cap,direct_scale=norm,raw=raw,theta=theta,
                T_second=T2,strong_second=strong2,raw_second=raw2,
                requested_eV=total_req,geom=geom,electronic_amplitude=electronic_amp)

def metrics(q,v,w,pos):
    """Mass metric is M I; q and v are three Cartesian physical coordinates."""
    z=q+1j*v/OMEGA
    total=float(np.sum(w[:,None]*abs(z)**2));N=float(w.sum())
    if total<=1e-100:return dict(gamma=0.,stripe=0.,q2=0.)
    gamma=float(np.sum(abs(np.sum(w[:,None]*z,axis=0))**2)/(N*total))
    # Existing opposite-sheet stripe: |q|=2pi/(5.3 A), along [110].
    wave=np.array([1.,1.,0.])/math.sqrt(2)*2*np.pi/5.3
    stripe=float(np.sum(abs(np.sum(w[:,None]*z*np.exp(-1j*(pos@wave))[:,None],axis=0))**2)/(N*total))
    tensor=np.einsum('i,ia,ib->ab',w,z,z)
    q2=float(np.linalg.norm(tensor)/max(total,1e-100))
    return dict(gamma=gamma,stripe=stripe,q2=q2)

def native_curves(env):
    """Native energies AND signed per-D cut derivatives; no numerical F0 bias."""
    if env=='bulk':
        rows=BULK['rows'];u=np.r_[0.,[r['u_A'] for r in rows]]
        conn=np.r_[0.,[r['pair_nonadditive_energy_eV'] for r in rows]]
        dc=np.r_[0.,[r['d_pair_energy_du_eV_A'] for r in rows]]
        ea=np.r_[0.,[r['one_A']['energy_free_eV']-BULK['zero_free_energy_eV'] for r in rows]]
        eb=np.r_[0.,[r['one_B']['energy_free_eV']-BULK['zero_free_energy_eV'] for r in rows]]
        da=np.r_[0.,[-np.dot(r['one_A']['forces_eV_A'][25],AXIS[env]) for r in rows]]
        db=np.r_[0.,[np.dot(r['one_B']['forces_eV_A'][63],AXIS[env]) for r in rows]]
    else:
        rows=VAC['rows'];u=np.r_[0.,[r['displacement_per_D_A'] for r in rows]]
        conn=np.r_[0.,[r['connected_eV'] for r in rows]]
        dc=np.r_[0.,[r['connected_dE_du_eV_A'] for r in rows]]
        ea=np.r_[0.,[r['A_eV'] for r in rows]];eb=np.r_[0.,[r['B_eV'] for r in rows]]
        slopes=[[],[]]
        for r in rows:
            for k,(kind,idx,sgn) in enumerate([('a',63,-1),('b',66,1)]):
                source=next(x['file'] for x in r['sources'] if f'-{kind}-' in x['file'])
                native=json.loads((ROOT/'g10_quantum_site/pd-vacancy-r4-six-D-pair/native'/source).read_text())
                slopes[k].append(sgn*np.dot(native['forces_eV_A'][idx],AXIS[env]))
        da=np.r_[0.,slopes[0]];db=np.r_[0.,slopes[1]]
    from scipy.interpolate import CubicHermiteSpline
    def curve(values,slope):
        ux=np.r_[u,u[-1]+.1];ex=np.r_[values,values[-1]+.05*slope[-1]]
        # Smooth even local quadratic+quartic matches BOTH first native data.
        # It avoids a cusp in W(u)/u² at the common-coordinate origin.
        quart=(slope[1]*u[1]-2*values[1])/(2*u[1]**4)
        quad=(4*values[1]-slope[1]*u[1])/(2*u[1]**2)
        return dict(spline=CubicHermiteSpline(ux,ex,np.r_[slope,0.],extrapolate=False),
                    u1=u[1],quad=quad,quart=quart)
    return [curve(conn,dc),curve(ea,da),curve(eb,db)],float(u[-1])
CURVES={env:native_curves(env) for env in D0}

def curve_value(c,u,limit,der=0):
    r=np.minimum(abs(u),limit+.1)
    a,b=c['quad'],c['quart']
    if der==0:local=a*r*r+b*r**4
    elif der==1:local=2*a*r+4*b*r**3
    else:local=2*a+12*b*r*r
    val=np.where(r<c['u1'],local,c['spline'](r,der))
    if der:val*=((abs(u)<limit+.1)*np.sign(u) if der==1 else abs(u)<limit+.1)
    return val

def connected_g(u,env,der=0):
    """One conditional smooth mirrored extension: Vc=-x*y*W((x-y)/2)/u²."""
    c,limit=CURVES[env][0][0],CURVES[env][1]
    r=abs(u);safe=np.maximum(r,c['u1'])
    w=curve_value(c,safe,limit)
    if der==0:
        outer=w/safe**2;inner=c['quad']+c['quart']*u*u
    elif der==1:
        d=curve_value(c,safe,limit,1)
        outer=(d/safe**2-2*w/safe**3)*np.sign(u);inner=2*c['quart']*u
    else:
        d=curve_value(c,safe,limit,1)
        outer=curve_value(c,safe,limit,2)/safe**2-4*d/safe**3+6*w/safe**4
        inner=np.full_like(u,2*c['quart'])
    return np.where(r<c['u1'],inner,outer)

def lattice_force(q,env,active=None):
    """Exact accepted G3 energy law in A/fs, with finite connected-cut replacement."""
    # q: (lane,block,3,nx,ny,nz), all local microblocks are physical FCC cells.
    linear=neighbor_force(q) if env=='bulk' else np.zeros_like(q)
    sq=q*q;r2=sq.sum(axis=2,keepdims=True)
    f=linear-(K+B*r2+(A-B)*sq)*q
    pot=.5*K*np.sum(q*q,axis=2)+.25*B*np.sum(q*q,axis=2)**2+.25*(A-B)*np.sum(q**4,axis=2)
    pot-=.5*np.sum(q*linear,axis=2)
    # One tagged native nearest pair in each physical microblock. 0 and +[110].
    flat=q.reshape(q.shape[:3]+(-1,));ff=f.reshape(f.shape[:3]+(-1,))
    ia,ib=(0,4) if env=='bulk' else (0,3)
    n=AXIS[env]
    qa=flat[:,:,:,ia];qb=flat[:,:,:,ib]
    x=qa@n;y=qb@n;u=.5*(x-y)
    curves,limit=CURVES[env]
    g=connected_g(u,env);dg=connected_g(u,env,1)
    pairmask=1.
    if active is not None:
        aa=active.reshape(active.shape[:2]+(-1,))
        pairmask=aa[:,:,ia]&aa[:,:,ib]
    # Replace only this bond's axial connected component, exactly once.
    j=float(JMATRIX[ia,ib]) if env=='bulk' else 0.
    fa=(y*(g-j)+.5*x*y*dg)*pairmask
    fb=(x*(g-j)-.5*x*y*dg)*pairmask
    ff[:,:,:,ia]+=fa[:,:,None]*n
    ff[:,:,:,ib]+=fb[:,:,None]*n
    pot.reshape(pot.shape[:2]+(-1,))[:,:,ia]+=(j-g)*x*y*pairmask
    if env=='vacancy':
        # Native cage has no calibrated mapping to this bulk microblock's bonds.
        # Replace tagged one-body *axial* potentials; rest is a declared G3 cage surrogate.
        for idx,c,sgn in [(ia,curves[1],1.),(ib,curves[2],-1.)]:
            x=np.sum(flat[:,:,:,idx]*n,axis=2)*sgn
            axquart=B+(A-B)*np.sum(n**4)
            old=.5*K*x*x+.25*axquart*x**4
            df=curve_value(c,x,limit,1)-(K*x+axquart*x**3)
            ff[:,:,:,idx]-=df[:,:,None]*n*sgn
            pot.reshape(pot.shape[:2]+(-1,))[:,:,idx]+=curve_value(c,x,limit)-old
        # Only the six original physical cage D exist in this coordinate chart.
        ff[:,:,:,6:]=0
        pot.reshape(pot.shape[:2]+(-1,))[:,:,6:]=0
    return f,pot,u

def hamiltonian_step(q,v,dt,env,active):
    """Energy-preserving discrete gradient of the *same* accepted finite energy.
    The local physical chart is kept bounded by promotion to the moving channel.
    """
    if dt<1e-5:
        force,_,_=lattice_force(q,env,active);force*=active[:,:,None,:,:,:]
        return q+dt*v+.5*dt*dt*force/M,v+dt*force/M
    x=q.copy();y=q+dt*v
    ia,ib=(0,4) if env=='bulk' else (0,3);n=AXIS[env]
    curves,limit=CURVES[env]
    aa=active.reshape(active.shape[:2]+(-1,));pairmask=aa[:,:,ia]&aa[:,:,ib]
    j=float(JMATRIX[ia,ib]) if env=='bulk' else 0.
    for iteration in range(12):
        mid=(x+y)/2
        linear=neighbor_force(mid) if env=='bulk' else np.zeros_like(mid)
        s0=np.sum(x*x,axis=2,keepdims=True);s1=np.sum(y*y,axis=2,keepdims=True)
        cubic=(x*x*x+x*x*y+x*y*y+y*y*y)/4
        force=linear-K*mid-B*(s0+s1)*(x+y)/4-(A-B)*cubic
        xf=x.reshape(x.shape[:3]+(-1,));yf=y.reshape(y.shape[:3]+(-1,))
        ff=force.reshape(force.shape[:3]+(-1,))
        xa=xf[:,:,:,ia]@n;xb=xf[:,:,:,ib]@n
        ya=yf[:,:,:,ia]@n;yb=yf[:,:,:,ib]@n
        u0=.5*(xa-xb);u1=.5*(ya-yb)
        def divided(c,a,b):
            diff=b-a
            out=curve_value(c,(a+b)/2,limit,1)
            np.divide(curve_value(c,b,limit)-curve_value(c,a,limit),diff,out=out,where=abs(diff)>1e-7)
            return out
        g0=connected_g(u0,env);g1=connected_g(u1,env)
        gd=connected_g((u0+u1)/2,env,1)
        np.divide(g1-g0,u1-u0,out=gd,where=abs(u1-u0)>1e-7)
        gb=(g0+g1)/2;pb=(xa*xb+ya*yb)/2
        ff[:,:,:,ia]+=((gb-j)*(xb+yb)/2+.5*pb*gd)[:,:,None]*n*pairmask[:,:,None]
        ff[:,:,:,ib]+=((gb-j)*(xa+ya)/2-.5*pb*gd)[:,:,None]*n*pairmask[:,:,None]
        if env=='vacancy':
            axquart=B+(A-B)*np.sum(n**4)
            for idx,c,sgn in [(ia,curves[1],1.),(ib,curves[2],-1.)]:
                a=np.sum(xf[:,:,:,idx]*n,axis=2)*sgn;b=np.sum(yf[:,:,:,idx]*n,axis=2)*sgn
                old=K*(a+b)/2+axquart*(a**3+a*a*b+a*b*b+b**3)/4
                ff[:,:,:,idx]-=(divided(c,a,b)-old)[:,:,None]*n*sgn
        force*=active[:,:,None,:,:,:]
        new=x+dt*v+dt*dt/(2*M)*force
        error=float(abs(new-y).max());y=new
        if error<2e-12:break
    else:raise ValueError('Physical-chart discrete gradient did not converge')
    return y,2*(y-x)/dt-v

def evolve(ec,sp,s,env):
    """Same conserved parent, physical local motifs, full-footprint receivers."""
    nodes=response_map(ec,sp,s,env,0.,'axial',False)
    nb=len(nodes['w']);micro=np.indices(SHAPE).reshape(3,-1).T
    micro=(micro-micro.mean(0))@g3.CELL.T*P['a_A']
    # Each block has the actual physical FCC spacing; centres cover the whole packet.
    pos=nodes['pos'][:,None,:]+micro[None,:,:]
    w=np.repeat((nodes['w']/np.prod(SHAPE))[:,None],len(micro),axis=1)
    if env=='vacancy':
        native=np.array(VREF['positions_A'])[63:69]
        pos[:,:6,:]=nodes['pos'][:,None,:]+native[None,:,:]-native.mean(0)
        w[:]=0;w[:,:6]=nodes['w'][:,None]/6
    wflat=w.ravel();position=pos.reshape(-1,3)
    consumer_take=np.searchsorted(np.cumsum(wflat)/wflat.sum(),(np.arange(64)+.5)/64)
    consumer_w=np.full(64,wflat.sum()/64)
    nbranch=len(BRANCHES);lanes=nbranch*3
    rng=np.random.default_rng(152027+int(ec*100)+sp)
    q=np.zeros((lanes,nb,3)+SHAPE);v=np.zeros_like(q)
    peak_q=np.zeros_like(q);peak_v=np.zeros_like(q);peak_t=np.zeros(lanes)
    mask=(w>0).reshape((nb,)+SHAPE)
    active=np.repeat(mask[None,...],lanes,axis=0)
    f,pot,u=lattice_force(q,env)
    Einit=np.sum((.5*M*np.sum(v*v,axis=2)+pot).reshape(lanes,nb,-1)*w[None,:,:],axis=(1,2))
    pool=np.zeros(lanes);out=np.zeros(lanes);work=np.zeros(lanes);direct=np.zeros(lanes);minimum_carrier=np.zeros(lanes)
    correction=np.zeros((lanes,11));correction_exit=np.zeros(lanes)
    carrier_integral=np.zeros((lanes,11));original_integral=np.zeros(11);carrier_records=[]
    limited=np.zeros(lanes);loss=np.zeros(lanes);maxdef=np.zeros(lanes)
    peak=np.zeros(lanes);minsep=np.full(lanes,D0[env]);maxrel=np.zeros(lanes)
    unsupported=np.zeros(lanes);summary=[];samples=[];detail=[]
    shells=np.array([float(r['D_D_separation_A'] if env=='bulk' else r['separation_A']) for r in CUTS[env]])
    pair_weight=nodes['w']/(6 if env=='vacancy' else np.prod(SHAPE))
    crossing=np.zeros((lanes,nb,len(shells)));residence=np.zeros_like(crossing)
    previous_vector=np.broadcast_to(D0[env]*AXIS[env],(lanes,nb,3)).copy()
    previous_common=np.zeros((lanes,nb));crossing_marks=[]
    from observables import ConsumerMeasures,inward_crossing
    measures=ConsumerMeasures(sys.modules[__name__],ec,sp,nodes['N'],wflat,pair_weight,shells,env)
    site_peak=np.zeros((lanes,nb,8));site_peak_q=np.zeros_like(q);site_peak_v=np.zeros_like(q)
    site_peak_t=np.zeros_like(site_peak);site_peak_pool=np.zeros_like(site_peak)
    last_record_time=0.
    import whole
    mode_names,qwave,omega=whole.mode_spec(sys.modules[__name__],ec,sp,s)
    plane=np.exp(-1j*(position@qwave.T))
    mode_q=np.zeros((lanes,len(omega),3),complex);mode_v=np.zeros_like(mode_q)
    mode_work=np.zeros((lanes,len(omega)));mode_loss=np.zeros_like(mode_work)
    mode_records=[];observer_records=[]
    escaped=np.zeros(lanes);escape_count=np.zeros(lanes);escape_marks=[]
    kappa=1/(s['tau_opt']*1e15)
    # F0 is evaluated at every physical site, including its signed radial parity.
    R=s['radius']*1e10
    ee_flat,grad,t0_flat,gclock=electronic_coordinate(ec,sp,s,position)
    ee=ee_flat.reshape(nb,-1);t0=t0_flat.reshape(nb,-1)
    polar=np.zeros_like(grad)
    if env=='vacancy':
        native=np.array(VREF['positions_A'])[63:69];n=native-native.mean(0);n/=np.linalg.norm(n,axis=1)[:,None]
        polar=np.tile(n[np.arange(len(micro))%6],(nb,1))*ee_flat[:,None]/.1
    F0=np.array([chi*(grad+polar) for _,chi in BRANCHES]).reshape(nbranch,nb,len(micro),3)
    G0=np.array([-chi*ee_flat[:,None]*gclock for _,chi in BRANCHES]).reshape(nbranch,nb,len(micro),3)
    F0=np.repeat(F0,3,axis=0).transpose(0,1,3,2).reshape(q.shape)
    G0=np.repeat(G0,3,axis=0).transpose(0,1,3,2).reshape(q.shape)
    # Causal uniform clock shifts supply the random-phase/common-resource control.
    F0[2::3]=0;G0[2::3]=0
    tarr=np.broadcast_to(t0.reshape((nb,)+SHAPE),(lanes,nb)+SHAPE).copy()
    shift=rng.uniform(0,2*np.pi/OMEGA,size=(nbranch,nb)+SHAPE)
    tarr[1::3]+=shift
    tr=p1.track(ec,sp);tz=np.r_[0.,[r['depth_A'] for r in tr['timeline']]]
    tt=np.r_[0.,[r['time_fs'] for r in tr['timeline']]]
    mechanical=np.interp(position[:,2],tz,tt).reshape((nb,)+SHAPE)
    marr=np.broadcast_to(mechanical,tarr.shape).copy();marr[1::3]+=shift
    from mechanics import ForceHistory
    from observer import CurvatureObserver
    mechanics=ForceHistory(sys.modules[__name__],ec,sp,position,env,nodes['direct_scale'])
    clock_shift=marr.reshape(lanes,-1)-mechanics.clock[None,:]
    observer=CurvatureObserver(sys.modules[__name__],env,position,wflat,plane,omega,
        np.mean(mechanical.reshape(nb,-1),axis=1))
    delay=.74/math.sqrt(2*tr['initial_eV']/M) if sp else 0.
    transit=mechanics.transit_fs.reshape((nb,)+SHAPE)
    path=stock_path(ec,sp,s,2.,np.r_[tarr.ravel(),marr.ravel(),tarr.ravel()+delay,marr.ravel()+delay,
        (marr-.25*transit).ravel(),(marr+.25*transit).ravel()])
    all_direct=float(path['dep'][:,2].sum())
    record_steps=np.unique(np.r_[1,np.linspace(1,len(path['t'])-1,40).astype(int)])
    y=path['y'];available_total=np.zeros(lanes);original_flow=0.
    optical_inflow_row=p1.network(s)[4].copy()*1e-15
    optical_inflow_row[4]=0.
    current_energy=Einit.copy()
    for step,(ta,tb) in enumerate(zip(path['t'][:-1],path['t'][1:]),1):
        dt=tb-ta;tm=(ta+tb)/2
        alpha=kappa/2;wd=np.sqrt(omega*omega-alpha*alpha)
        co=np.cos(wd*dt);si=np.sin(wd*dt);decay=np.exp(-alpha*dt)
        before_mode=.5*M*wflat.sum()*(np.sum(abs(mode_v)**2,axis=2)+omega[None,:]**2*np.sum(abs(mode_q)**2,axis=2))
        oldq=mode_q.copy()
        mode_q=decay*(mode_q*(co+alpha/wd*si)[None,:,None]+mode_v*(si/wd)[None,:,None])
        mode_v=decay*(mode_v*(co-alpha/wd*si)[None,:,None]-oldq*(omega*omega/wd*si)[None,:,None])
        after_mode=.5*M*wflat.sum()*(np.sum(abs(mode_v)**2,axis=2)+omega[None,:]**2*np.sum(abs(mode_q)**2,axis=2))
        mode_loss+=before_mode-after_mode
        key=path['step_keys'][step-1];eprop,fprop,_=path['propagator'][key]
        correction_integral=correction@fprop.T
        step_integral=path['integral'][step][None,:]+correction_integral
        carrier_integral+=step_integral;original_integral+=path['integral'][step]
        available_total+=np.maximum(step_integral@optical_inflow_row,0.)
        original_flow+=path['flow'][step]
        out+=kappa*step_integral[:,4]
        oldsum=correction.sum(axis=1)
        correction=correction@eprop.T
        correction_exit+=oldsum-correction.sum(axis=1)
        pool=y[step,4]+correction[:,4]
        if pool.min() < -1e-7:raise ValueError(f'Negative optical pool at {tb} fs: {pool.tolist()}; carrier minimum {(y[step][None,:]+correction).min(axis=1).tolist()}')
        pool=np.maximum(pool,0.)
        # Discrete-gradient Hamiltonian propagation, exact damping, signed source kick.
        q,v=hamiltonian_step(q,v,dt,env,active)
        before=.5*M*np.sum(v*v,axis=2)
        v*=math.exp(-kappa*dt)
        after=.5*M*np.sum(v*v,axis=2)
        dissip=np.sum((before-after).reshape(lanes,nb,-1)*w[None,:,:],axis=(1,2))
        loss+=dissip
        # Coherent loss is routed once, into the same accepted destinations.
        # The unforced prior's decay is subtracted: equilibrium is not beam energy.
        refd=np.repeat(dissip.reshape(-1,3)[:,2],3)
        correction[:,8]+=.2*(dissip-refd);correction[:,6]+=.8*(dissip-refd)
        # Integrate the actual signed ZBL derivative over this interval, including
        # its leading field. Each constituent sees the state left by the first.
        raw_dp=mechanics.interval(ta,tb,clock_shift)
        if sp:raw_dp+=mechanics.interval(ta,tb,clock_shift+delay)
        dp=raw_dp.reshape(lanes,nb,len(micro),3).transpose(0,1,3,2).reshape(q.shape)
        dp[2::3]=0.
        # Allocate direct bound work from the *existing* struck-D channel.
        dp*=active[:,:,None,:,:,:]
        oldkin=.5*M*np.sum(v*v,axis=2)
        dl=np.sum((np.sum(v*dp,axis=2)).reshape(lanes,nb,-1)*w[None,:,:],axis=(1,2))
        dq=np.sum((.5*np.sum(dp*dp,axis=2)/M).reshape(lanes,nb,-1)*w[None,:,:],axis=(1,2))
        available=np.maximum(0.,y[step,2]+correction[:,2])
        scale=work_fraction(dl,dq,available)
        v+=dp/M*scale[:,None,None,None,None,None]
        applied_impulse=dp*scale[:,None,None,None,None,None]
        signed=dl*scale+dq*scale*scale
        direct+=signed;correction[:,2]-=signed
        age=tm-tarr
        beta=1/(s['tau_opt']*1e15);rho=1/(s['tau_ne']*1e15)
        def pulse_value(a):
            b=np.maximum(a,0.)
            value=np.where(a>=0,np.exp(-beta*b)-np.exp(-(beta+rho)*b),0.)
            derivative=np.where(a>=0,-beta*np.exp(-beta*b)+(beta+rho)*np.exp(-(beta+rho)*b),0.)
            return value,derivative
        local,derivative=pulse_value(age)
        if sp:
            l2,d2=pulse_value(age-delay);local+=l2;derivative+=d2
        F=F0*local[:,:,None,:,:,:]+G0*derivative[:,:,None,:,:,:]
        F*=active[:,:,None,:,:,:]
        dv=F*dt/M
        site_linear=M*np.sum(v*dv,axis=2);site_quad=.5*M*np.sum(dv*dv,axis=2)
        requested=np.sum((site_linear+site_quad).reshape(lanes,nb,-1)*w[None,:,:],axis=(1,2))
        # A fractional quadrature weight cannot give one physical D more than
        # the entire incoming parent's available optical work.
        local_scale=work_fraction(site_linear,site_quad,np.broadcast_to(pool[:,None,None,None,None],site_linear.shape))
        dv*=local_scale[:,:,None,:,:,:]
        linear=np.sum((M*np.sum(v*dv,axis=2)).reshape(lanes,nb,-1)*w[None,:,:],axis=(1,2))
        quad=np.sum((.5*M*np.sum(dv*dv,axis=2)).reshape(lanes,nb,-1)*w[None,:,:],axis=(1,2))
        # Work limiter solves quadratic signed kick work, not an arbitrary phase share.
        scale=work_fraction(linear,quad,pool)
        v+=dv*scale[:,None,None,None,None,None]
        actual=linear*scale+quad*scale*scale
        pool-=actual;correction[:,4]-=actual;work+=actual;limited+=np.maximum(requested-actual,0)
        applied_impulse+=M*dv*scale[:,None,None,None,None,None]
        projected=np.einsum('lnc,n,nm->lmc',applied_impulse.transpose(0,1,3,4,5,2).reshape(lanes,-1,3),wflat,plane,optimize=True)/wflat.sum()
        before_kick=.5*M*wflat.sum()*np.sum(abs(mode_v)**2,axis=2)
        mode_v+=projected/M
        mode_work+=.5*M*wflat.sum()*np.sum(abs(mode_v)**2,axis=2)-before_kick
        if pool.min() < -1e-7:raise ValueError(f'Work limiter {tb} fs: pool={pool}, b={linear}, a={quad}, fraction={scale}')
        # Leave the known bound chart before applying a polynomial outside it.
        qcap=.36 if env=='bulk' else .24
        escaping=active&(np.sum(q*q+M/K*v*v,axis=2)>qcap*qcap)
        if np.any(escaping):
            _,prepot,_=lattice_force(q,env,active)
            pre=.5*M*np.sum(v*v,axis=2)+prepot
            for lane,block,i,j,k in np.argwhere(escaping):
                site=np.ravel_multi_index((i,j,k),SHAPE)
                escape_marks.append([tb,lane,block,site,*q[lane,block,:,i,j,k],*v[lane,block,:,i,j,k],w[block,site]])
            active&=~escaping
            q*=active[:,:,None,:,:,:];v*=active[:,:,None,:,:,:]
            _,postpot,_=lattice_force(q,env,active)
            post=.5*M*np.sum(v*v,axis=2)+postpot
            moved=np.sum((pre-post).reshape(lanes,nb,-1)*w[None,:,:],axis=(1,2))
            escaped+=moved;escape_count+=np.sum(escaping.reshape(lanes,nb,-1)*w[None,:,:],axis=(1,2))
            refm=np.repeat(moved.reshape(-1,3)[:,2],3)
            correction[:,2]+=moved-refm
        f,pot,u=lattice_force(q,env,active);f*=active[:,:,None,:,:,:]
        current_energy=np.sum((.5*M*np.sum(v*v,axis=2)+pot).reshape(lanes,nb,-1)*w[None,:,:],axis=(1,2))
        closure=current_energy-Einit+loss+escaped-direct-work
        maxdef=np.maximum(maxdef,abs(closure))
        reference_energy=np.repeat(current_energy.reshape(-1,3)[:,2],3)
        excess=current_energy-reference_energy
        rising=excess>peak
        peak_q[rising]=q[rising];peak_v[rising]=v[rising];peak_t[rising]=tb
        peak=np.maximum(peak,excess)
        flatq=q.reshape(lanes,nb,3,-1)
        flatv=v.reshape(lanes,nb,3,-1);ia,ib=(0,4) if env=='bulk' else (0,3)
        pairvec=D0[env]*AXIS[env]+flatq[:,:,:,ib]-flatq[:,:,:,ia]
        pairvel=flatv[:,:,:,ia]-flatv[:,:,:,ib]
        sep=np.linalg.norm(pairvec,axis=2)
        alive=active.reshape(lanes,nb,-1)
        pair_active=alive[:,:,ia]&alive[:,:,ib]
        bound_D=np.sum(alive*w[None,:,:],axis=(1,2))
        minsep=np.minimum(minsep,np.where(pair_active,sep,np.inf).min(axis=1))
        rel=.25*M*np.sum(pairvel*pairvel,axis=2)
        radial_v=np.sum(pairvel*pairvec,axis=2)/sep
        entered,root=inward_crossing(previous_vector,pairvec,shells)
        entered&=pair_active[:,:,None]
        common=.5*((flatq[:,:,:,ia]+flatq[:,:,:,ib])@AXIS[env])
        for lane,block,shell in np.argwhere(entered):
            fraction=root[lane,block,shell]
            vector=previous_vector[lane,block]+fraction*(pairvec[lane,block]-previous_vector[lane,block])
            velocity=-(pairvec[lane,block]-previous_vector[lane,block])/dt
            coordinate=.5*(D0[env]-vector@AXIS[env])
            centre=previous_common[lane,block]+fraction*(common[lane,block]-previous_common[lane,block])
            native=0.;supported=0<=coordinate<=CURVES[env][1]
            if supported:
                native=sum(float(curve_value(c,np.array(coordinate),CURVES[env][1])) for c in CURVES[env][0])
            connected=-(centre+coordinate)*(centre-coordinate)*float(connected_g(np.array(coordinate),env))
            crossing_marks.append([ta+fraction*dt,dt,lane,block,shell,*vector,*((M/2)*velocity),
                .25*M*float(velocity@velocity),native,connected,coordinate,centre,pair_weight[block],
                supported,env=='vacancy' or coordinate>.36])
        crossing+=entered
        residence+=((sep[:,:,None]<=shells)&(radial_v[:,:,None]>0)&pair_active[:,:,None])*dt
        previous_vector=pairvec.copy();previous_common=common.copy()
        maxrel=np.maximum(maxrel,np.where(pair_active,rel,0.).max(axis=1))
        unsupported+=np.mean((u<0)|(u>CURVES[env][1]),axis=1)*dt
        observer.add(tb,dt,q,active)
        minimum_carrier=np.minimum(minimum_carrier,(y[step][None,:]+correction).min(axis=1))
        if minimum_carrier.min() < -1e-6:raise ValueError('A physical carrier inventory became negative')
        v_all=v.transpose(0,1,3,4,5,2).reshape(lanes,-1,3)
        measures.add(dt,v_all,pool,pairvec,pairvel,pair_active)
        local_reference=.5*K*np.sum(q*q,axis=2)+.5*M*np.sum(v*v,axis=2)
        local_total=local_reference.reshape(lanes,nb,8)+pool[:,None,None]/nodes['N']
        rising_site=local_total>site_peak
        np.maximum(site_peak,local_total,out=site_peak)
        site_peak_t[rising_site]=tb
        site_peak_pool[rising_site]=np.broadcast_to(pool[:,None,None]/nodes['N'],site_peak.shape)[rising_site]
        np.copyto(site_peak_q,q,where=rising_site.reshape((lanes,nb)+SHAPE)[:,:,None])
        np.copyto(site_peak_v,v,where=rising_site.reshape((lanes,nb)+SHAPE)[:,:,None])
        if step in record_steps:
            qflat=q.transpose(0,1,3,4,5,2).reshape(lanes,-1,3)
            vflat=v.transpose(0,1,3,4,5,2).reshape(lanes,-1,3)
            Q=np.einsum('lnc,n,nm->lmc',qflat,wflat,plane,optimize=True)/wflat.sum()
            V=np.einsum('lnc,n,nm->lmc',vflat,wflat,plane,optimize=True)/wflat.sum()
            observed=observer.read(tb);observer_records.append(observed)
            h=observed['h'][:,:len(omega)]
            harmonic=.5*M*(np.einsum('lnc,n->l',vflat*vflat,wflat)[:,None]+omega[None,:]**2*np.einsum('lnc,n->l',qflat*qflat,wflat)[:,None])
            mode_records.append((tb,Q.copy(),V.copy(),mode_q.copy(),mode_v.copy(),mode_work.copy(),mode_loss.copy(),h,harmonic))
            carrier_records.append((tb,(y[step][None,:]+correction).copy(),y[step].copy(),
                carrier_integral.copy(),original_integral.copy()))
            for branch,(label,chi) in enumerate(BRANCHES):
                base=branch*3;ref=base+2
                for ctrl in (0,1):
                    lane=base+ctrl
                    qf=q[lane].transpose(0,2,3,4,1).reshape(-1,3)
                    vf=v[lane].transpose(0,2,3,4,1).reshape(-1,3)
                    qr=q[ref].transpose(0,2,3,4,1).reshape(-1,3)
                    vr=v[ref].transpose(0,2,3,4,1).reshape(-1,3)
                    mm=metrics(qf,vf,wflat,position)
                    coherent=current_energy[lane]-current_energy[ref]
                    samples.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,branch=label,
                        control='built' if ctrl==0 else 'matched_random',time_fs=tb,D_coordinate_excess_eV=coherent,
                        optical_population_eV=pool[lane],optical_inflow_eV=available_total[lane],electronic_work_eV=work[lane],
                        bound_chart_occupied_D=bound_D[lane],
                        mechanical_signed_work_eV=direct[lane],coherent_loss_eV=loss[lane]-loss[ref],population_loss_eV=out[lane],
                        energy_residual_eV=closure[lane]-closure[ref],
                        all_carrier_residual_eV=float(correction[lane].sum()+coherent+correction_exit[lane]),
                        canonical_Gamma_coordinate_fraction=mm['gamma'],canonical_stripe_coordinate_fraction=mm['stripe'],
                        canonical_q_squared_proxy=mm['q2'],
                        h_window_available=bool(observed['window_available'][lane,0]),
                        h_trend_robust=bool(observed['valid'][lane,0]),
                        h_2omega_fit=float(h[lane,0]) if observed['window_available'][lane,0] else None,
                        h_trend_leakage=float(observed['leak'][lane,0]) if observed['window_available'][lane,0] else None,
                        h_over_threshold=float(h[lane,0]/(2*kappa/omega[0])) if observed['valid'][lane,0] else None,
                        local_closing_oscillatory_reference_eV=float(observed['local_oscillatory_reference_eV'][lane]),
                        local_closing_fit_coverage=float(np.average(observed['local_fit_valid'][lane],weights=observer.block_w)),
                        local_closing_clock_mode_fraction=float(observed['local_clock_mode_fraction'][lane]),
                        minimum_pair_A=float(sep[lane].min()),maximum_relative_kinetic_eV=float(rel[lane].max()),
                        reference_energy_eV=current_energy[ref],reference_drift_eV=current_energy[ref]-Einit[ref]))
                    # Actual joint state quadrature drives all consumers, not independent marginal extrema.
                    detail.append(dict(key=(ec,sp,s['name'],env,label,ctrl),t=tb,dt=tb-last_record_time,
                        q=qf[consumer_take].astype(np.float32),v=vf[consumer_take].astype(np.float32),
                        reference_v=vr[consumer_take].astype(np.float32),w=consumer_w,pos=position[consumer_take],
                        pool=max(0.,pool[lane]),pair_u=u[lane].copy(),
                        pair_r=sep[lane].copy(),pair_v=np.sum(pairvel[lane]*pairvec[lane],axis=1)/sep[lane],
                        pair_kin=rel[lane].copy(),reference_pair_r=sep[ref].copy(),
                        reference_pair_v=np.sum(pairvel[ref]*pairvec[ref],axis=1)/sep[ref],
                        pair_vec=pairvec[lane].copy(),pair_velocity=pairvel[lane].copy(),
                        pair_w=pair_weight*pair_active[lane],energy=coherent,source_sites=nodes['N']))
            last_record_time=tb
    exposure_expected=nodes['N']*path['t'][-1]
    exposure_defect=float(abs(measures.moments[:,:,:,0].sum(axis=2)-exposure_expected).max())
    if exposure_defect>1e-9*max(1.,exposure_expected):raise ValueError('Consumer physical interval exposure failed')
    np.savez_compressed(HERE/f'peak_{ec:g}_{sp+1}_{s["name"]}_{env}.npz',
        position_A=position,occupied_weight=wflat,branch=np.array([x[0] for x in BRANCHES]),
        q_A=peak_q[::3].astype(np.float32),v_A_fs=peak_v[::3].astype(np.float32),time_fs=peak_t[::3],
        reference_q_A=peak_q[2::3].astype(np.float32),reference_v_A_fs=peak_v[2::3].astype(np.float32),
        escape_marks=np.array(escape_marks).reshape(-1,11),
        escape_columns=np.array(['time_fs','lane','block','site','qx_A','qy_A','qz_A','vx_A_fs','vy_A_fs','vz_A_fs','occupied_D_weight']),
        pair_shell_A=shells,coherent_inward_crossings=np.sum(crossing*pair_weight[None,:,None],axis=1),
        coherent_shell_residence_fs=np.sum(residence*pair_weight[None,:,None],axis=1),
        mode_names=np.array(mode_names),mode_q_A_inverse=qwave,mode_omega_fs_inverse=omega,
        mode_time_fs=np.array([x[0] for x in mode_records]),actual_Q_A=np.array([x[1] for x in mode_records]),
        actual_V_A_fs=np.array([x[2] for x in mode_records]),green_Q_A=np.array([x[3] for x in mode_records]),
        green_V_A_fs=np.array([x[4] for x in mode_records]),green_work_eV=np.array([x[5] for x in mode_records]),
        green_loss_eV=np.array([x[6] for x in mode_records]),mode_h=np.array([x[7] for x in mode_records]),
        harmonic_reference_eV=np.array([x[8] for x in mode_records]),
        mode_h_trend_leakage=np.array([x['leak'] for x in observer_records]),
        mode_h_window_available=np.array([x['window_available'] for x in observer_records]),
        mode_h_trend_robust=np.array([x['valid'] for x in observer_records]),
        local_modal_amplitude_A=np.array([x['local_amplitude'] for x in observer_records]).astype(np.complex64),
        local_modal_fit_valid=np.array([x['local_fit_valid'] for x in observer_records]),
        local_modal_fit_spread_A=np.array([x['local_fit_spread'] for x in observer_records]).astype(np.float32),
        local_modal_slow_offset_A=np.array([x['local_slow_offset'] for x in observer_records]).astype(np.float32),
        local_mode_vector=observer.vector,local_mode_omega_fs_inverse=observer.local_omega,
        local_mode_info=np.array(json.dumps(observer.local_info)),local_block_clock_phase=observer.local_phase,
        original_carriers_eV=np.array([x[2] for x in carrier_records]),
        corrected_carriers_eV=np.array([x[1] for x in carrier_records]),
        carrier_integral_eV_fs=np.array([x[3] for x in carrier_records]),
        original_carrier_integral_eV_fs=np.array([x[4] for x in carrier_records]),
        mechanical_candidate_transfer_eV=mechanics.T,mechanical_strong_mark=mechanics.strong,
        mechanical_raw_force_phasor=mechanics.raw,
        pair_crossings_by_block=crossing,pair_weights=pair_weight,
        pair_shell_crossing_marks=np.array(crossing_marks).reshape(-1,19),
        pair_crossing_columns=np.array(['time_linear_fs','bracket_fs','lane','block','shell',
            'rx_A','ry_A','rz_A','px_eV_fs_A','py_eV_fs_A','pz_eV_fs_A','relative_kinetic_eV',
            'native_opposed_cut_reference_eV','model_connected_eV','u_A','common_A','pair_weight',
            'within_native_positive_cut','PAW_overlap_uncertain']),
        site_peak_total_reference_eV=site_peak,site_peak_q_A=site_peak_q.astype(np.float32),
        site_peak_v_A_fs=site_peak_v.astype(np.float32),site_peak_time_fs=site_peak_t,
        site_peak_population_eV_D=site_peak_pool,beam_exposure_defect_D_fs=exposure_defect,
        **measures.arrays())
    for branch,(label,chi) in enumerate(BRANCHES):
        lane=branch*3;ref=lane+2
        summary.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,branch=label,
            chi_gradient=chi,chi_vacancy_per_A=chi/.1 if env=='vacancy' else 0.,
            occupied_D=nodes['N'],physical_microblocks=nb,physical_sites_per_microblock=8 if env=='bulk' else 6,
            numeric_slots_per_microblock=8,quadrature_nodes=nb*(8 if env=='bulk' else 6),
            final_bound_chart_D=bound_D[lane],chart_exit_meaning='conditional moving/unresolved continuation, not demonstrated ejection',
            peak_coherent_excess_eV=peak[lane],electronic_signed_work_eV=work[lane],mechanical_signed_work_eV=direct[lane],
            original_direct_D_eV=all_direct,strong_channel_eV=max(0.,all_direct-direct[lane]),
            escaped_from_bound_eV=escaped[lane]-escaped[ref],escaped_occupied_D=escape_count[lane],
            optical_positive_flow_eV=available_total[lane],original_optical_positive_flow_eV=original_flow,
            consumer_exposure_defect_D_fs=exposure_defect,
            requested_work_rejected_eV=limited[lane],
            minimum_carrier_inventory_eV=minimum_carrier[lane],
            max_integrator_energy_defect_eV=maxdef[lane],final_energy_residual_eV=float(current_energy[lane]-current_energy[ref]+loss[lane]-loss[ref]+escaped[lane]-escaped[ref]-direct[lane]-work[lane]),
            minimum_pair_A=minsep[lane],peak_pair_relative_kinetic_eV=maxrel[lane],time_beyond_native_support_fs=unsupported[lane],
            optical_lifetime_fs=s['tau_opt']*1e15,loss_branch='0.2 thermal electrons + 0.8 Pd high; original rates unchanged',
            coupling='exact G3 physical microblock; tagged axial bond replaced by native connected cut; off-cut extension conditional',epoch=EPOCH))
    return summary,samples,detail

def phase_maps():
    rows=[];marks=[];arrays={}
    for ec in EC:
        for sp in (0,1):
            for s in p1.SCENARIOS:
                for env in D0:
                    signatures={}
                    for label,chi in BRANCHES:
                        for geom in (['single'] if sp==0 else ['axial','transverse']):
                            if chi in (0.,1.):
                                d=response_map(ec,sp,s,env,chi,geom,True);signatures[geom,chi]=d
                            else:
                                base=signatures[geom,0.];unit=signatures[geom,1.]
                                amp=base['amp']+chi*unit['electronic_amplitude']
                                d=dict(base,amp=amp,F=chi*unit['F'],G=chi*unit['G'],
                                    electronic_amplitude=chi*unit['electronic_amplitude'],
                                    requested_eV=float(.5*K*np.sum(base['w'][:,None]*abs(amp)**2)))
                            amp=d['amp'];w=d['w'];pos=d['pos'];mm=metrics(amp.real,OMEGA*amp.imag,w,pos)
                            tr=p1.track(ec,sp);v0=math.sqrt(2*tr['initial_eV']/M)
                            ktime=OMEGA/v0
                            neighbor=float(np.cos(ktime*P['a_A']/math.sqrt(2)))
                            rows.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,branch=label,geometry=geom,
                                chi=chi,occupied_D=d['N'],quadrature_nodes=len(w),requested_free_harmonic_energy_eV=d['requested_eV'],
                                gamma_fraction=mm['gamma'],stripe_fraction=mm['stripe'],squared_phasor_fraction=mm['q2'],
                                local_neighbor_timing_correlation=neighbor,timing_wavelength_nm=2*np.pi/ktime/10,
                                direct_CSDA_allocation_scale=d['direct_scale'],bound_cap_eV=d['cap'],
                                strong_weighted_D=float(w@d['strong']),strong_raw_recoil_eV=float(w@(d['T']*d['strong'])),
                                second_strong_weighted_D=float(w@d['strong_second']),
                                second_strong_raw_recoil_eV=float(w@(d['T_second']*d['strong_second'])),
                                conditional='unbudgeted linear signature; actual capped nonlinear work is in patch_summary.csv'))
                            if chi==0:
                                for i in range(len(w)):
                                    marks.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,geometry=geom,node=i,
                                        second_raw_recoil_eV=d['T_second'][i],second_strong_moving_mark=bool(d['strong_second'][i]),
                                        x_A=pos[i,0],y_A=pos[i,1],z_A=pos[i,2],occupied_D_weight=w[i],
                                        mechanical_arrival_fs=d['t_mech'][i],electronic_arrival_fs=d['t0'][i],
                                        projectile_energy_eV=d['en'][i],raw_recoil_eV=d['T'][i],strong_moving_mark=bool(d['strong'][i]),
                                        force_impulse_x_eV_fs_A=d['raw'][i,0].real,force_impulse_y_eV_fs_A=d['raw'][i,1].real,
                                        force_impulse_z_real_eV_fs_A=d['raw'][i,2].real,force_impulse_z_imag_eV_fs_A=d['raw'][i,2].imag))
                            prefix=f'E{ec:g}_D{sp+1}_{s["name"]}_{env}_{geom}'
                            if chi==0:
                                arrays[prefix+'_mechanical_amplitude_A']=amp.astype(np.complex64)
                                arrays[prefix+'_position_A']=pos;arrays[prefix+'_occupied_weight']=w
                                arrays[prefix+'_mechanical_time_fs']=d['t_mech'];arrays[prefix+'_electronic_time_fs']=d['t0']
                                arrays[prefix+'_candidate_first_T_eV']=d['T'];arrays[prefix+'_candidate_second_T_eV']=d['T_second']
                                arrays[prefix+'_first_strong']=d['strong'];arrays[prefix+'_second_strong']=d['strong_second']
                            elif chi==1:
                                arrays[prefix+'_electronic_amplitude_per_chi_A']=d['electronic_amplitude'].astype(np.complex64)
                                arrays[prefix+'_F_per_chi_eV_A']=d['F'].astype(np.float32)
                                arrays[prefix+'_G_per_chi_eV_fs_A']=d['G'].astype(np.float32)
    table('phase_map.csv',rows);table('all_kick_marks.csv',marks);np.savez_compressed(HERE/'phase_fields.npz',**arrays)
    return rows

def beam_factors(beam):
    xx=np.linspace(-.005,.005,101);x,y=np.meshgrid(xx,xx,indexing='ij')
    w=p1.beam_weight(x,y,beam)
    spatial=float(np.mean(w*w)/np.mean(w)**2)
    temporal=(1+.7**2/2+18*.1+81*.1)/1.9**2 if beam=='hot_bursts' else 1.
    return spatial*temporal

def current(ec):
    return np.interp(2*ec,np.arange(.5,15.01,.5),[38,68,94,110,118,123,130,134,136,138,140,145,152,156,160,165,171,176,180,182,192,207,217,222,228,238,251,263,275,283])*1e-6

def signed_rate(lab,v,ch,theta,U):
    vb=math.sqrt(2*lab*1000/M)
    rel=v.copy();rel[:,2]-=vb
    er=.5*M*np.sum(rel*rel,axis=1)/1000
    return p1.moving_rate(er,theta,ch,U)*np.sqrt(er/lab)


def positive_loss_check():
    q=np.zeros((1,1,3)+SHAPE);v=np.zeros_like(q);active=np.ones((1,1)+SHAPE,bool)
    force=np.zeros_like(q);force[:,:,0]=.01
    dt=.5;kappa=1/(p1.SCENARIOS[1]['tau_opt']*1e15);work=0.;loss=0.
    for _ in range(200):
        q,v=hamiltonian_step(q,v,dt,'bulk',active)
        before=.5*M*np.sum(v*v);v*=math.exp(-kappa*dt);loss+=before-.5*M*np.sum(v*v)
        before=.5*M*np.sum(v*v);v+=force*dt/M;work+=.5*M*np.sum(v*v)-before
    _,pot,_=lattice_force(q,'bulk');energy=float(.5*M*np.sum(v*v)+np.sum(pot))
    amplitude=float(abs(q).max());times=[];energies=[]
    for step in range(6000):
        q,v=hamiltonian_step(q,v,dt,'bulk',active);v*=math.exp(-kappa*dt)
        if step%8==0:
            _,pot,_=lattice_force(q,'bulk')
            times.append((step+1)*dt);energies.append(float(.5*M*np.sum(v*v)+np.sum(pot)))
    fitted=-float(np.polyfit(times,np.log(energies),1)[0])
    result=dict(positive_forced_max_displacement_A=amplitude,positive_forced_work_eV=work,
        positive_forced_energy_eV=energy,positive_forced_loss_eV=loss,positive_forced_closure_eV=energy+loss-work,
        source_off_fitted_energy_kappa_fs_inverse=fitted,required_energy_kappa_fs_inverse=kappa,
        source_off_energy_lifetime_fs=1/fitted,source_off_amplitude_lifetime_from_energy_fs=2/fitted,
        source_off_fit_relative_error=abs(fitted/kappa-1),correct_energy_loss_checked=True,
        damping_convention='velocity-only drag exp(-kappa dt); underdamped mode amplitude exp(-kappa t/2); energy envelope exp(-kappa t)')
    if result['source_off_fit_relative_error']>.01 or abs(result['positive_forced_closure_eV'])>1e-8:
        raise ValueError('Positive/source-off energy-loss control failed')
    return result

def controls():
    shape=(1,1,3)+SHAPE;q=np.zeros(shape);f,pot,u=lattice_force(q,'bulk')
    pure_zero=float(abs(f).max())
    positive=positive_loss_check();dt=.5;active=np.ones((1,1)+SHAPE,bool)
    # Both the zero scalar state and a physical +/- prior use this same consumer.
    zq=np.zeros(shape);zv=np.zeros(shape)
    zq,zv=hamiltonian_step(zq,zv,dt,'bulk',active);zv-=dt*.1*K*zq/M
    zero_launch=float(abs(zq).max()+abs(zv).max())
    pq=np.zeros((2,1,3)+SHAPE);pv=np.zeros_like(pq);pq[0,:,0]=.05;pq[1,:,0]=-.05
    pa=np.ones((2,1)+SHAPE,bool)
    for _ in range(200):
        pq,pv=hamiltonian_step(pq,pv,dt,'bulk',pa)
        pv-=dt*.1*K*pq/M
    result=dict(zero_source_force_max_eV_A=pure_zero,uniform_scalar_zero_classical_launch_A=zero_launch,
        **positive,
        even_prior_mean_A=float(pq.mean()),even_prior_second_moment_A2=float(np.mean(pq**2)),
        energy_loss_lifetimes_fs=[s['tau_opt']*1e15 for s in p1.SCENARIOS],
        energy_loss_convention='kappa is energy envelope loss; velocity-only drag uses the full kappa',
        force_source='exact G3 environment and coupling functions; physical mass M',
        timing='each source has H(t-t0), D2 second inserted after actual axial delay; transverse is synchronous')
    if pure_zero!=0 or result['positive_forced_work_eV']<=0 or result['positive_forced_max_displacement_A']<=0:raise ValueError('Forced consumer control failed')
    dump('controls.json',result);return result

def figures(samples):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10,7))
    for s in p1.SCENARIOS:
        for branch,style in [('absent','--'),('gradient_large_plus','-')]:
            r=[x for x in samples if x['E_cm_keV']==2.5 and x['species']==1 and x['scenario']==s['name'] and x['environment']=='bulk' and x['branch']==branch and x['control']=='built']
            t=np.array([x['time_fs'] for x in r])/1000
            axes[0,0].plot(t,[x['D_coordinate_excess_eV'] for x in r],style,label=s['name']+' '+branch)
            axes[0,1].plot(t,[float(x['h_over_threshold']) if x['h_over_threshold'] not in (None,'') else np.nan for x in r],style)
            axes[1,0].plot(t,[x['local_closing_clock_mode_fraction'] for x in r],style)
            axes[1,1].plot(t,[x['minimum_pair_A'] for x in r],style)
    for ax,title in zip(axes.ravel(),['D-coordinate excess / parent (eV)','Resolved curvature pump / threshold','Local closing-mode clock fraction','Joint tagged-pair distance (A)']):ax.set(xlabel='Time (ps)',title=title);ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=6);fig.tight_layout();fig.savefig(HERE/'phase_evolution.png',dpi=130);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,7))
    for row,env in enumerate(D0):
        for ax,s in zip(axes[row],p1.SCENARIOS):
            with np.load(HERE/f'peak_2.5_1_{s["name"]}_{env}.npz',allow_pickle=False) as data:
                pos=data['position_A'];q=data['q_A'][1].transpose(0,2,3,4,1).reshape(-1,3)
                v=data['v_A_fs'][1].transpose(0,2,3,4,1).reshape(-1,3);keep=data['occupied_weight']>0
                z=(q+1j*v/OMEGA)@AXIS[env];valid=keep&(abs(z)>1e-15)
                sc=ax.scatter(pos[valid,2]/10,(pos[valid]@AXIS[env])/10,c=np.angle(z[valid]),s=12,cmap='twilight',vmin=-np.pi,vmax=np.pi)
                ax.scatter(pos[keep&~valid,2]/10,(pos[keep&~valid]@AXIS[env])/10,c='gray',s=5)
                strong=keep&data['mechanical_strong_mark']
                ax.scatter(pos[strong,2]/10,(pos[strong]@AXIS[env])/10,marker='x',c='red',s=12+5*np.log1p(data['mechanical_candidate_transfer_eV'][strong]))
            ax.set(title=env+' / '+s['name'],xlabel='Depth (nm)',ylabel='Native pair-axis projection (nm)')
    fig.colorbar(sc,ax=axes,label='Canonical pair-axis phase at common bound-energy peak; not a 2ω fit')
    fig.savefig(HERE/'phase_map.png',dpi=130);plt.close(fig)

def source_index():
    paths=['g15_track/pass1/inputs.json','g15_track/pass1/run.py','g15_track/pass1/bca.py','g15_track/pass1/yield_consumers.csv',
        'g15_track/pass1/joined_beam_history.npz','g15_track/pass1/packet_marks.csv','g15_track/pass1/whole_foil.csv',
        'g3_lattice/coherence/lattice.py','g10_quantum_site/pdd-nearest-pair/pair-wells.json',
        'g10_quantum_site/pd-vacancy-r4-six-D-pair/pair.json','g10_quantum_site/pd-vacancy-r4-six-D-pair/native/reference.json']
    paths=sorted(set(paths+['g12_residual/inputs_ledger.md']+
        ['g10_quantum_site/pd-vacancy-r4-six-D-pair/native/'+x['file'] for r in VAC['rows'] for x in r['sources']]))
    dump('source_index.json',dict(base='c037c8c479bc08f6ae3ed01916cf8582cafdcf38',sources=[dict(path=p,sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest(),bytes=(ROOT/p).stat().st_size) for p in paths],epoch=EPOCH))

def save_joint(details):
    arrays={k:np.array([d[k] for d in details]) for k in details[0] if k!='key'}
    arrays['keys']=np.array([list(map(str,d['key'])) for d in details])
    np.savez_compressed(HERE/'joint_states.npz',**arrays)

def saved_tables(name):
    with (HERE/name).open() as f:rows=list(csv.DictReader(f))
    for row in rows:
        for k,v in row.items():
            try:row[k]=float(v)
            except (ValueError,TypeError):pass
    return rows

def force_checks():
    state=np.sin(np.arange(3*np.prod(SHAPE))).reshape((3,)+SHAPE)*.037
    mask=np.ones(SHAPE,bool);mask[0,1,0]=False;mask[1,0,1]=False
    state*=mask
    exact=g3.coupling_force(state,JFFT,SHAPE)*mask
    sampled=neighbor_force(state)*mask
    defect=float(abs(exact-sampled).max())
    ec=2.5;sp=0;s=p1.SCENARIOS[0]
    z,ts,mid,logu,normal=electron_profile(ec,sp,s['radius']*1e10)
    j=len(z)//2
    pos=np.array([[s['radius']*1e10/3,s['radius']*1e10/7,z[j]+.37*(z[j+1]-z[j])]])
    u,grad,clock,gclock=electronic_coordinate(ec,sp,s,pos)
    age=s['tau_ne']*1e15*.4;at=float(clock[0]+age)
    beta=1/(s['tau_opt']*1e15);rho=1/(s['tau_ne']*1e15)
    f=math.exp(-beta*age)-math.exp(-(beta+rho)*age)
    fp=-beta*math.exp(-beta*age)+(beta+rho)*math.exp(-(beta+rho)*age)
    expected=float(grad[0,2]*f-u[0]*gclock[0,2]*fp)
    h=min(.001,(z[j+1]-z[j])*.001)
    def coordinate(point):
        uu,gg,tt,gtt=electronic_coordinate(ec,sp,s,point)
        a=at-tt
        return float((uu*np.where(a>=0,np.exp(-beta*np.maximum(a,0))-np.exp(-(beta+rho)*np.maximum(a,0)),0.))[0])
    plus=pos.copy();minus=pos.copy();plus[0,2]+=h;minus[0,2]-=h
    finite=(coordinate(plus)-coordinate(minus))/(2*h)
    relative=abs(finite-expected)/max(abs(expected),1e-15)
    null_force=np.zeros(3)*f-u[0]*np.zeros(3)*fp
    result=dict(masked_three_component_operator_defect_eV_A=defect,
        axial_force_over_chi_eV_A=expected,axial_fixed_time_difference_eV_A=finite,axial_relative_difference=relative,
        source_position_A=pos[0],observation_time_fs=at,finite_difference_step_A=h,
        actual_row_clock_derivative_fs_A=gclock[0,2],uniform_scalar_force_max_eV_A=float(abs(null_force).max()),
        causal_kernel_at_onset=math.exp(0.)-math.exp(0.),
        interface='continuous log-linear source coordinate; one-sided derivatives at row-centre knots; accepted piecewise-linear arrival clock',
        force='chi[grad(u0)f - u0 grad(t0)fprime]; f(0)=0, no onset impulse',
        boundary='same physical periodic eight-site G3 motif and mask; dense matrix sampled from authoritative function')
    if defect>1e-12 or relative>2e-5:raise ValueError('Actual force/operator witness failed')
    dump('force_checks.json',result)
    return result

def load_source():
    with np.load(HERE/'joint_states.npz',allow_pickle=False) as opened:
        arrays={k:opened[k] for k in opened.files}
    details=[]
    for i,key in enumerate(arrays['keys']):
        d={k:arrays[k][i] for k in arrays if k!='keys'}
        d['key']=(float(key[0]),int(key[1]),str(key[2]),str(key[3]),str(key[4]),int(key[5]));details.append(d)
    return saved_tables('patch_summary.csv'),saved_tables('patch_times.csv'),details

def producer_fingerprint():
    import inspect,whole
    functions=[electron_profile,electronic_coordinate,map_nodes,force_signature,response_map,stock_path,work_fraction,
        metrics,neighbor_force,native_curves,curve_value,connected_g,lattice_force,hamiltonian_step,evolve,whole.mode_spec]
    source=''.join(inspect.getsource(f) for f in functions)+json.dumps([P,BRANCHES,SHAPE],sort_keys=True)
    for name in ('mechanics.py','observer.py','observables.py','source_index.json'):source+=(HERE/name).read_text()
    return hashlib.sha256(source.encode()).hexdigest()

def force_fingerprint():
    import inspect
    return hashlib.sha256((''.join(inspect.getsource(f) for f in
        (electron_profile,electronic_coordinate,neighbor_force))+json.dumps([P,SHAPE],sort_keys=True)).encode()).hexdigest()

def hamiltonian_fingerprint():
    import inspect
    inputs=[r for r in json.loads((HERE/'source_index.json').read_text())['sources']
        if r['path'].startswith(('g10_quantum_site','g3_lattice'))]
    return hashlib.sha256((''.join(inspect.getsource(f) for f in
        (native_curves,curve_value,connected_g,lattice_force,hamiltonian_step))+json.dumps([P,inputs],sort_keys=True)).encode()).hexdigest()

def main(consume=False,resume=False,source_only=False,workers=3):
    start=time.monotonic()
    if consume:
        receipt=json.loads((HERE/'source_run.json').read_text())
        if not receipt['normal_exit'] or receipt['source_cases']!=72 or receipt['physics_fingerprint']!=producer_fingerprint():
            raise ValueError('Consumer requires the completed, matching 72-case source epoch; do not combine partial or stale archives')
        control=json.loads((HERE/'controls.json').read_text());phase=saved_tables('phase_map.csv')
        summary,samples,details=load_source()
    else:
        source_index();control=json.loads((HERE/'controls.json').read_text()) if (HERE/'controls.json').exists() else controls()
        if control.get('hamiltonian_fingerprint')!=hamiltonian_fingerprint():
            import checks as joining
            control=joining.run(sys.modules[__name__],control)
            control['hamiltonian_fingerprint']=hamiltonian_fingerprint();dump('controls.json',control)
        fingerprint=producer_fingerprint();contract=HERE/'producer_contract.json'
        if resume and (HERE/'joint_states.npz').exists():
            if not contract.exists() or json.loads(contract.read_text())['physics_fingerprint']!=fingerprint:
                raise ValueError('Saved trajectories belong to another force contract; do not reuse them')
            summary,samples,details=load_source()
        else:summary=[];samples=[];details=[]
        dump('producer_contract.json',dict(physics_fingerprint=fingerprint,source='signed retarded density gradient and exact sampled G3 operator'))
        checks=json.loads((HERE/'force_checks.json').read_text()) if (HERE/'force_checks.json').exists() else {}
        if checks.get('force_contract')!=force_fingerprint():
            checks=force_checks();checks['force_contract']=force_fingerprint();dump('force_checks.json',checks)
        phase=phase_maps()
        complete={(d['key'][0],d['key'][1],d['key'][2],d['key'][3]) for d in details}
        summary=[r for r in summary if (float(r['E_cm_keV']),int(r['species'])-1,r['scenario'],r['environment']) in complete]
        samples=[r for r in samples if (float(r['E_cm_keV']),int(r['species'])-1,r['scenario'],r['environment']) in complete]
        todo=[(ec,sp,s,env) for ec in EC for sp in (0,1) for s in p1.SCENARIOS for env in D0
            if (ec,sp,s['name'],env) not in complete]
        from concurrent.futures import ProcessPoolExecutor,as_completed
        with ProcessPoolExecutor(max_workers=workers) as executor:
            jobs={executor.submit(evolve,*case):case for case in todo}
            for job in as_completed(jobs):
                case=jobs[job];a,b,c=job.result();summary+=a;samples+=b;details+=c
                table('patch_summary.csv',summary);table('patch_times.csv',samples);save_joint(details)
                print(json.dumps(dict(done=[case[0],case[1]+1,case[2]['name'],case[3]])),flush=True)
        dump('source_run.json',dict(source_cases=len(summary)//len(BRANCHES),force_branches=len(BRANCHES),
            workers=workers,numerical_threads_per_worker=1,elapsed_s=time.monotonic()-start,
            physics_fingerprint=fingerprint,normal_exit=True))
        if source_only:return
    import whole,linear_join,consumers
    context=sys.modules[__name__]
    modes,transport=whole.global_modes(context,summary)
    linear_join.execute(context,summary)
    windows=whole.marked_windows(context)
    whole.patch_level_rates(context)
    y,p,r=consumers.run(context,details,summary,samples);figures(samples)
    # Accepted long-time foil closure is reused; signed phase-routing differences are explicit above.
    result=dict(model='conditional signed causal phase refinement; weighted physical microblock closure',epoch=EPOCH,
        full_source_cases=72,patch_cases=len(summary),phase_cases=len(phase),yield_rows=len(y),pair_rows=len(p),rate_rows=len(windows),
        maximum_patch_integrator_defect_eV=max(x['max_integrator_energy_defect_eV'] for x in summary),
        max_abs_signed_energy_residual_eV=max(abs(x['final_energy_residual_eV']) for x in summary),
        maximum_supported_beam_envelope_fraction_of_primary_residual=max(x['beam_envelope_mixture_fraction_of_residual'] for x in y),
        maximum_actual_h_over_threshold=max((float(x['h_over_threshold']) for x in samples if x['h_over_threshold'] not in (None,'')),default=None),
        mechanical_and_electronic_signed_work_separate=True,electron_susceptibility_unmeasured=True,
        prior='cold excess trajectories; preserved phase-symmetric scalar-prior control is separate; zero-point Hartree in K once',
        three_largest_omissions=['matched PdD transient force/curvature susceptibilities','interblock phase transport and vacancy cross-couplings','joint finite-shell to nuclear-contact transmission and physical vacancy prevalence'],
        controls=control,elapsed_s=time.monotonic()-start)
    dump('results.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--controls',action='store_true');parser.add_argument('--consume',action='store_true')
    parser.add_argument('--resume',action='store_true');parser.add_argument('--source-only',action='store_true')
    parser.add_argument('--workers',type=int,choices=range(1,6),default=3);args=parser.parse_args()
    HERE.mkdir(parents=True,exist_ok=True)
    if args.controls:controls()
    else:main(args.consume,args.resume,args.source_only,args.workers)
