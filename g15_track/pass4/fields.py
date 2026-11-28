"""Whole-foil finite-volume enthalpy, characteristic transport and reciprocal modes.
SI internally. Radiation is linearized about 300 K, not an unnamed heat sink.
"""
import math
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from sources import L,H,RHO,CV,KTH,ALPHA,BULK,QE,MD,HB,KB,optical
N=17; DX=L/N
ZE=np.array([0.,.25e-6,125e-6,H]); DZ=np.diff(ZE)
X=(np.arange(N)+.5)*DX-L/2; Z=(ZE[:-1]+ZE[1:])/2
XX,YY,ZZ=np.meshgrid(X,X,Z,indexing='ij'); VOL=np.broadcast_to(DX*DX*DZ,(N,N,3)).copy(); CAP=CV*VOL
BANDS=2.5e5*10.**np.arange(8)
CHANNELS=('end_D','Pd_recoil','moving_D','hot_e','D_optical','Pd_diffuse_low','Pd_ballistic_high','transient_structure','thermal_e','primary_D')
DIRECTIONS=np.r_[np.eye(3),-np.eye(3)]

class Calendar:
    """One physical Cox history, exactly integrated; Poisson marks are conditioned on it."""
    def __init__(self,seconds=2.,seed=441):
        rng=np.random.default_rng(seed);t=0.;on=bool(rng.random()<.1);edges=[0.];state=[]
        while t<seconds:
            t=min(seconds,t+rng.exponential(2e-6 if on else 18e-6));state.append(on);edges.append(t);on=not on
        self.edges=np.array(edges);self.state=np.array(state);self.seconds=seconds
        self.phase=float(rng.uniform(0,2*np.pi));self.omega=2*np.pi*1e5
        self.burst_integral=np.r_[0.,np.cumsum(np.diff(self.edges)*self.state)]
        self.norm=self.primitive(seconds)/seconds
    def primitive(self,t):
        a=np.asarray(t);c=np.minimum(np.maximum(a,0.),self.seconds);j=np.minimum(np.searchsorted(self.edges,c,side='right')-1,len(self.state)-1)
        burst=self.burst_integral[j]+(c-self.edges[j])*self.state[j]
        return c+.7/self.omega*(np.sin(self.omega*c+self.phase)-np.sin(self.phase))+9*burst
    def integral(self,a,b,name):
        if name!='hot_bursts':return np.maximum(np.asarray(b)-np.maximum(a,0.),0.)
        return (self.primitive(b)-self.primitive(a))/self.norm
    def value(self,t,name):
        t=np.asarray(t);j=np.minimum(np.maximum(np.searchsorted(self.edges,t,side='right')-1,0),len(self.state)-1)
        value=(1+.7*np.cos(self.omega*t+self.phase)+9*self.state[j])/self.norm if name=='hot_bursts' else np.ones_like(t)
        return value*(t>=0)*(t<=self.seconds)

class Heat:
    """Positive backward-Euler conservative FV, compatible thermoelastic enthalpy."""
    def __init__(self):
        shape=(N,N,3);n=np.prod(shape);rr=[];cc=[];vv=[];diag=np.zeros(n);self.exit=np.zeros((3,n))
        ids=np.arange(n).reshape(shape)
        for axis in range(3):
            for index in np.ndindex(shape):
                if index[axis]==shape[axis]-1:continue
                other=list(index);other[axis]+=1;other=tuple(other);i=ids[index];j=ids[other]
                area=DX*DZ[index[2]] if axis<2 else DX*DX
                distance=DX if axis<2 else (DZ[index[2]]+DZ[other[2]])/2
                g=KTH*area/distance;diag[i]+=g;diag[j]+=g;rr.extend((i,j));cc.extend((j,i));vv.extend((-g,-g))
        rad=4*.1*5.670374419e-8*300**3
        for index in np.ndindex(shape):
            i=ids[index]
            for axis in range(3):
                for side in (0,shape[axis]-1):
                    if index[axis]!=side:continue
                    if axis<2:
                        g=DX*DZ[index[2]]/(1/1e4+DX/(2*KTH));self.exit[0,i]+=g
                    else:
                        g=DX*DX/(1/(5+rad)+DZ[index[2]]/(2*KTH));self.exit[1,i]+=g*5/(5+rad);self.exit[2,i]+=g*rad/(5+rad)
        diag+=self.exit.sum(axis=0);rr.extend(range(n));cc.extend(range(n));vv.extend(diag)
        self.G=sparse.coo_matrix((vv,(rr,cc)),shape=(n,n)).tocsc();self.cached={}
    def step(self,u,amount,dt,capacity=None):
        cap=CAP.ravel() if capacity is None else np.asarray(capacity).ravel()
        if np.any(cap<=0):raise ValueError('Nonpositive physical enthalpy derivative')
        # Capacity changes are retained in the physical solve, not clipped away.
        if capacity is None:
            factor=self.cached.get(dt)
            if factor is None:factor=splu(sparse.diags(cap)+dt*self.G);self.cached[dt]=factor
        else:factor=splu(sparse.diags(cap)+dt*self.G)
        temp=factor.solve(np.asarray(u).ravel()+np.asarray(amount).ravel())
        exits=dt*(self.exit@temp);new=cap*temp
        defect=float(new.sum()+exits.sum()-np.asarray(u).sum()-np.asarray(amount).sum())
        if temp.min()<-1e-9:raise ValueError('Negative excess enthalpy; no numerical clipping')
        return new.reshape(VOL.shape),temp.reshape(VOL.shape),exits,defect


def deposit(points,weights):
    """Conservative point-to-cell map. Specular geometry is applied before this call."""
    p=np.asarray(points);ix=np.minimum(np.maximum(((p[:,0]+L/2)/DX).astype(int),0),N-1)
    iy=np.minimum(np.maximum(((p[:,1]+L/2)/DX).astype(int),0),N-1)
    iz=np.minimum(np.maximum(np.searchsorted(ZE,p[:,2],side='right')-1,0),2)
    index=(ix*N+iy)*3+iz;weights=np.asarray(weights)
    if np.iscomplexobj(weights):
        out=np.bincount(index,weights=weights.real,minlength=N*N*3)+1j*np.bincount(index,weights=weights.imag,minlength=N*N*3)
    else:out=np.bincount(index,weights=weights,minlength=N*N*3)
    return out.reshape(VOL.shape)


def shift(field,delta,reflect=.8,phase_sign=-1.):
    """Finite-volume uniform-cell characteristic. Exact translation overlap in x/y.

    Each face crossing attenuates by its stated reflection energy; removed
    energy is a boundary port, not inelastic absorption in the interior.
    The six-ray quadrature has at most one moving axis per call.
    """
    axis=int(np.argmax(abs(delta)));distance=float(delta[axis]);a=np.asarray(field)
    if np.ndim(reflect):reflect=reflect[axis]
    if distance==0:return a.copy(),0.,a.copy()
    centers=(XX,YY,ZZ)[axis].ravel();length=L if axis<2 else H;origin=-L/2 if axis<2 else 0.
    width=DX if axis<2 else DZ[np.indices(VOL.shape)[2]].ravel()
    # Two Gauss points preserve the in-cell launch distribution and boundary proximity.
    out=np.zeros_like(a);signed=np.zeros_like(a);escaped=0.
    for node in (-1/math.sqrt(3),1/math.sqrt(3)):
        start=centers+node*width/2-origin;raw=start+distance
        segment=np.floor(raw/length).astype(int);bounces=abs(segment)
        fold=np.mod(raw,2*length);where=np.where(fold<=length,fold,2*length-fold)+origin
        survival=reflect**bounces;sign=phase_sign**bounces
        points=np.column_stack((XX.ravel(),YY.ravel(),ZZ.ravel()));points[:,axis]=where
        weight=a.ravel()/2;out+=deposit(points,weight*survival);signed+=deposit(points,weight*np.sqrt(survival)*sign)
        escaped+=np.sum(weight*(1-survival))
    return out,escaped,signed


def continuous_shift(field,delta,reflect=.8):
    """Overlap-exact subcell transfer for displacements smaller than one cell."""
    axis=int(np.argmax(abs(delta)));distance=abs(float(delta[axis]));size=DX if axis<2 else DZ[0]
    if np.ndim(reflect):reflect=reflect[axis]
    if distance>=size:return shift(field,delta,reflect)[:2]
    fraction=distance/size if axis<2 else distance/DZ[None,None,:]
    step=1 if delta[axis]>=0 else -1
    transfer=field*fraction;out=field-transfer
    rolled=np.roll(transfer,step,axis=axis);boundary=[slice(None)]*3;boundary[axis]=0 if step>0 else -1
    lost=float(np.sum(rolled[tuple(boundary)]))*(1-reflect)
    # Reflected energy remains in its originating boundary cell.
    reflected=np.zeros_like(field);origin=[slice(None)]*3;origin[axis]=-1 if step>0 else 0
    reflected[tuple(origin)]=rolled[tuple(boundary)]*reflect
    rolled[tuple(boundary)]=0.;out=out+rolled+reflected
    return out,lost


def signed_shift(field,delta,reflect=.8):
    """First-moment normal-strain transfer; same geometry, signed free return."""
    axis=int(np.argmax(abs(delta)));distance=abs(float(delta[axis]));size=DX if axis<2 else DZ[0]
    if np.ndim(reflect):reflect=reflect[axis]
    if distance>=size:return shift(field,delta,reflect,-1.)[2]
    fraction=distance/size if axis<2 else distance/DZ[None,None,:]
    step=1 if delta[axis]>=0 else -1
    transfer=field*fraction;out=field-transfer;rolled=np.roll(transfer,step,axis=axis)
    boundary=[slice(None)]*3;boundary[axis]=0 if step>0 else -1
    origin=[slice(None)]*3;origin[axis]=-1 if step>0 else 0
    out[tuple(origin)]-=math.sqrt(reflect)*rolled[tuple(boundary)]
    rolled[tuple(boundary)]=0.
    return out+rolled


def vector_shift(field,delta,reflect=(.8,.8,1.),signed=False):
    moved=field;escaped=0.
    for axis in range(3):
        if delta[axis]==0:continue
        step=np.zeros(3);step[axis]=delta[axis]
        if signed:moved=signed_shift(moved,step,reflect)
        else:
            moved,lost=continuous_shift(moved,step,reflect);escaped+=lost
    return moved,escaped


def modal_transition(omega,gamma,dt):
    a=gamma/2;wd=np.sqrt(omega*omega-a*a);e=np.exp(-a*dt);s=np.sin(wd*dt)/wd;c=np.cos(wd*dt)
    return e*(c+a*s),e*s,-e*omega*omega*s,e*(c-a*s)


def oscillator(q,p,force,omega,gamma,dt):
    aa,ab,ba,bb=modal_transition(omega,gamma,dt);eq=force/(omega*omega);old=q
    q=aa*(q-eq)+ab*p+eq;p=ba*(old-eq)+bb*p
    work=float(force@(q-old))
    return q,p,work


def interaction_force(S,u0,pressure_volume,optical_stress_volume):
    """The derivative of the calibrated thermoelastic and optical Hamiltonian.
    Native electronic pressure is a signed work port, not ideal-gas 2E/3.
    """
    coefficient=3*BULK*ALPHA/CV*u0.ravel()+pressure_volume.ravel()-optical_stress_volume.ravel()
    return coefficient@S


def optical_factor(strain,temp,K0,dKdT,dKde):
    return 1+(dKde*(strain/3+ALPHA*temp)+dKdT*temp)/(2*K0)


def ray_fields(source,power,frequency,vg,kappa,scatter,calendar,time,beam,reflect=(.8,.8,1.),directions=None,diffusion_directions=None,wavevector=None,birth_delay=0.):
    """Positive retarded RTE age Green with elastic direction randomization.

    Scalar populations carry no invented phase. signed first/second coherent
    moments are sampled separately on actual physical parent windows.
    """
    tau=1/kappa;maxage=min(time,28*tau)
    if maxage<=0:return np.zeros_like(source),np.zeros_like(source),0.,0.,np.zeros_like(source,dtype=complex)
    edges=np.r_[0.,np.geomspace(max(1e-18,maxage*1e-7),maxage,32)]
    stored=np.zeros_like(source);loss=np.zeros_like(source);escape_power=0.;pathmoment=0.;phase_green=np.zeros_like(source,dtype=complex)
    rays=DIRECTIONS if directions is None else np.asarray(directions)
    diffuse_rays=DIRECTIONS if diffusion_directions is None else np.asarray(diffusion_directions)
    source_phase=source*np.exp(-2j*np.pi*frequency*birth_delay)
    for a,b in zip(edges[:-1],edges[1:]):
        age=(a+b)/2;weight=(math.exp(-kappa*a)-math.exp(-kappa*b))/kappa
        clock=calendar.value(time-age-birth_delay,beam);drive=clock*power
        ballistic=math.exp(-scatter*age);diffdistance=min(vg*age,math.sqrt(2*vg*vg/max(scatter,1e-30)*age))
        # Linear unit-impulse first-moment Green, in seconds. Its squared
        # modulus is not the RTE energy or a new source of power.
        phase_green_source=clock*source_phase
        for direction in rays:
            carrier=2*np.pi*frequency-(float(np.dot(wavevector,direction*vg)) if wavevector is not None else 0.)
            pole=kappa/2+scatter+1j*carrier
            phase_weight=(np.exp(-pole*a)-np.exp(-pole*b))/pole
            moved,_=vector_shift(phase_green_source,direction*vg*age,reflect,signed=True)
            phase_green+=phase_weight*moved/len(rays)
        for fraction,distance,angular in ((ballistic,vg*age,rays),(1-ballistic,diffdistance,diffuse_rays)):
            if fraction==0:continue
            for direction in angular:
                moved,escaped=vector_shift(drive*source,direction*distance,reflect)
                amount=weight*fraction/len(angular)
                stored+=amount*moved;loss+=kappa*amount*moved
                escape_power+=amount*kappa*escaped;pathmoment+=amount*float(np.sum(drive*source))*distance*distance
    return stored,loss,escape_power,pathmoment,phase_green


def path_green(arm,p,seed=446):
    """Conservative marked flight histories for every short-carrier loss port.

    A row is energy visiting one propagating segment. Its terminal position and
    clock are exactly the next row's source. No decay ends in a heat stock.
    Direction-randomizing collisions do not remove energy. Heat is an endpoint
    into the whole-foil FV enthalpy operator.
    """
    import native
    from sources import inputs
    rng=np.random.default_rng(seed);p_lattice=inputs()['lattice'];basis=native.primitive_basis(p_lattice);segments=[];endpoints=[]
    baseline_w,baseline_v,baseline_D,_=native.dispersion(p_lattice,basis,arm['allocation'])
    pole,_,_,_=native.thermal_pole(p_lattice,baseline_D,300.)
    pole_w=np.sqrt(pole*QE/1e-20/(p_lattice['mass_eVc2']*QE/299792458.**2))
    baseline_v*=baseline_w[:,None]/pole_w[:,None];baseline_w=pole_w
    # Origins: cutoff D, host recoil, struck D, nonthermal electrons.
    for origin in range(4):
        for sample in range(128):
            channel=(0,1,2,3)[origin];pos=np.zeros(3);age=0.;weight=1/128;phase=float(rng.uniform(0,2*np.pi))
            for visit in range(15):
                tau=[arm['end'],arm['recoil'],.73*arm['recoil'],arm['ne'],arm['opt'],2*arm['pd'],arm['pd'],arm['struct'],arm['te']][channel]
                duration=float(rng.exponential(tau));direction=rng.normal(size=3);direction/=np.linalg.norm(direction)
                omega=0.;wave=np.zeros(3);speed=0.;scatter=0.;displacement=None
                if channel in (0,1,2):
                    mass=MD if channel!=1 else 106.42*1.66053906660e-27
                    e=100. if channel==0 else (20. if channel==1 else 100.)
                    speed=math.sqrt(2*e*QE/mass);distance=2*speed*tau*(-math.expm1(-duration/(2*tau)))
                elif channel in (3,8):
                    speed=1e6;scatter=speed/arm['lambda_e'];flight=min(duration,1/scatter)
                    after=max(0.,duration-flight);random=rng.normal(size=3);random/=np.linalg.norm(random)
                    displacement=direction*speed*flight+random*min(speed*after,math.sqrt(2*speed*arm['lambda_e']*after))
                elif channel==4:
                    j=int(rng.integers(len(basis)));wave=basis[j].copy()
                    # The Γ representative is a finite-width curvature focus,
                    # not an immobile stock. ±q nodes sample its packet bandwidth.
                    if np.linalg.norm(wave)<1.:
                        wave=direction/arm['radius'];ww,vv,DD,_=native.dispersion(p_lattice,wave[None,:],arm['allocation'])
                        KK,_,_,_=native.thermal_pole(p_lattice,DD,300.)
                        omega=float(np.sqrt(KK[0]*QE/1e-20/(p_lattice['mass_eVc2']*QE/299792458.**2)))
                        velocity=vv[0]*ww[0]/omega
                    else:velocity=baseline_v[j];omega=float(baseline_w[j])
                    speed=float(np.linalg.norm(velocity));direction=velocity/speed if speed>0 else direction;distance=speed*duration
                else:
                    import launch
                    pol=int(rng.integers(3));speeds,_,groups=launch.eigenbasis(direction);cap=(3.4e12,4.7e12,6.8e12)[pol]
                    q=float(rng.random()**(1/3)*np.pi/2);wave=direction*(2*np.pi*cap/speeds[pol])*q
                    omega=2*np.pi*cap*np.sin(q);velocity=groups[pol]*np.cos(q)
                    speed=float(np.linalg.norm(velocity));scatter=speed/arm['lambda_pd'];flight=min(duration,1/scatter)
                    random=rng.normal(size=3);random/=np.linalg.norm(random);after=max(0.,duration-flight)
                    displacement=velocity*flight+random*min(speed*after,math.sqrt(2*speed*arm['lambda_pd']*after))
                finish=pos+(direction*distance if displacement is None else displacement)
                segments.append((origin,sample,channel,age,age+duration,*pos,*finish,weight,omega,phase,*wave,speed,scatter))
                pos=finish;age+=duration
                if channel==0:nextch=rng.choice([3,4,6,7],p=[.3,.25,.4,.05])
                elif channel==1:nextch=rng.choice([4,5,6],p=[.05,.35,.6])
                elif channel==2:nextch=rng.choice([3,4,6],p=[.25,.5,.25])
                elif channel==3:nextch=rng.choice([4,5,6],p=[arm['share'],(1-arm['share'])*.3,(1-arm['share'])*.7])
                elif channel==4:nextch=rng.choice([8,6],p=[.2,.8])
                else:nextch=-1
                if nextch<0:
                    endpoints.append((origin,sample,age,*pos,weight));break
                channel=int(nextch)
            else:raise ValueError('Unclosed carrier path')
    return np.array(segments),np.array(endpoints)
