"""Signed retarded receiver projection on the accepted eight source ports.

The reported mechanical observable is volumetric strain tr(epsilon), and K tr(epsilon)
is tensile-positive mean normal stress. The sampler returns its isotropic projection,
not a reconstructed shear tensor. Complex fields and energy are retained separately.
"""
from functools import lru_cache
from pathlib import Path
import math, time
import numpy as np
from scipy.special import j0
from numpy.polynomial.legendre import leggauss, legval, legder
from scipy.linalg import eigh
from common import HERE, QE, L, H, RHO, ENERGIES, BEAMS, inherited, emit, table

BULK = (227.1e9 + 2 * 176.04e9) / 3
FREQUENCY_COMPARATORS = (2.5054719912489504e-8, .00556844834632475)

def packet_duration(bandwidth,omega):
    """Named positive-bandwidth scenario, at least two complete carrier cycles."""
    return 2*np.pi*max(2,int(np.ceil(omega/bandwidth)))/omega


def wavelet(age, bandwidth, omega):
    """Causal Hann derivative; zero DC, orthogonal real quadratures, unit energy.

    Integer N>=2 makes integral(z²)=0 as well as integral(|z|²)=1 for
    z=wavelet*exp(-i omega t). Thus sqrt(2)*Re(A z) has energy |A|² for
    every launch phase, without a nonlinear receiver normalization.
    """
    duration=packet_duration(bandwidth,omega);t=np.clip(age,0.,duration)
    g=np.sin(np.pi*t/duration)**2
    derivative=np.pi/duration*np.sin(2*np.pi*t/duration)
    norm=np.sqrt(3*omega*omega*duration/8+np.pi*np.pi/(2*duration))
    return np.where((age>=0)&(age<=duration),(omega*g+1j*derivative)/norm,0.)

def directions(guided, order=4):
    az = 2 * np.pi * (np.arange(8) + .5) / 8
    if guided:
        return np.c_[np.cos(az), np.sin(az), np.zeros(8)], np.full(8, 1 / 8)
    mu, weight = leggauss(order)
    n = np.stack(np.broadcast_arrays(np.sqrt(1-mu[:, None]**2)*np.cos(az),
        np.sqrt(1-mu[:, None]**2)*np.sin(az), mu[:, None]), axis=-1).reshape(-1, 3)
    return n, np.repeat(weight / 16, 8)


def build_spectra(out, path4, energies=ENERGIES, closures=None):
    """Re-express exact frozen coefficients; do not repeat a source producer."""
    import json
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    s,native,fields,launch,response,acoustic,yielding=inherited(path4)
    tracks=json.loads((Path(path4)/'frozen/tracks.json').read_text())
    norm=np.load(Path(path4)/'frozen/native_transfer.npz')['normalization']
    arms=[a for a in s.ARMS if a['allocation']=='onsite' and (closures is None or a['closure'] in closures)]
    rows=[];start=time.monotonic()
    for arm in arms:
        for ec in energies:
            path=out/f"spectrum_E{ec:g}_{arm['closure']}.npz"
            if path.exists():continue
            save={}
            for sp in (0,1):
                tr=tracks[f'{ec:g}_{sp}'];marks=launch.recoil_marks(tr)
                q,lam,displ,duration,number=marks
                rng=np.random.default_rng(7721+sp)
                fluctuations=[rng.poisson(lam,size=(64,len(lam)))-lam for _ in range(sp+1)]
                accepted=np.load(Path(path4)/f"source_E{ec:g}_{arm['closure']}_sp{sp}.npz")
                z,t,dep,E,dt,dp=launch.track_arrays(tr)
                boundaries=np.r_[0.,[r['depth_A']*1e-10 for r in tr['timeline']]]
                delay=.74e-10/np.sqrt(2*tr['initial_eV']*s.QE/s.MD) if sp else 0.
                for pi,ch in enumerate(launch.ports()):
                    nn=accepted['directions'];ww=np.full(6,1/6)
                    if ch['guided']:ww[[2,5]]=0.;ww/=ww.sum()
                    values=accepted['source_amplitudes'][pi,:,:,:4].copy()
                    multiplicity=4 if ch['guided'] else 6
                    scale=np.sqrt(norm[sp,pi]*multiplicity)
                    values*=scale[:,None,None]
                    stochastic=np.zeros((3,6,64),complex);kval=np.zeros((3,6))
                    omega=2*np.pi*ch['frequency']
                    for pol in range(3):
                        for j,n in enumerate(nn):
                            speeds,ev,group=launch.eigenbasis(n);cap=(3.4e12,4.7e12,6.8e12)[pol]
                            k=2*np.pi*cap/speeds[pol]*np.arcsin(ch['frequency']/cap);kval[pol,j]=k
                            phase=np.exp(1j*omega*t-1j*k*n[2]*z)*np.sinc((omega*dt-k*n[2]*np.diff(boundaries))/(2*np.pi))
                            coefficient=(q@ev[:,pol])*np.expm1(1j*omega*duration-1j*k*(displ@n))*np.repeat(phase,number)*np.exp(-.5*(k*2e-9)**2)
                            random=fluctuations[0]@coefficient
                            if sp:random+=np.exp(1j*omega*delay)*(fluctuations[1]@coefficient)
                            stochastic[pol,j]=random*scale[pol]
                    for name,a in [('n',nn),('weight',ww),('a',values),('vg',accepted['group_speed_m_s'][pi]),('k',kval),('polar',accepted['polarization'][pi]),('stochastic',stochastic)]:
                        save[f'{name}_{sp}_{pi}']=a
                    total=values[...,0]+values[...,2]+values[...,3]
                    energy=np.sum(ww[None,:]*(abs(total)**2+abs(values[...,1])**2))
                    target=float(accepted['energy_J'][pi].sum())
                    rows.append(dict(E_cm_keV=ec,closure=arm['closure'],species=sp+1,port=pi,frequency_Hz=ch['frequency'],
                        energy_J=target,replayed_source_energy_J=float(energy),source_relative_energy_defect=float(energy/target-1) if target else 0.,
                        stochastic_variance_relative_sample_error=float(np.sum(ww[None,:]*(np.mean(abs(stochastic)**2,axis=-1)-abs(values[...,1])**2))/max(float(np.sum(ww[None,:]*abs(values[...,1])**2)),1e-300)),
                        semantics='accepted complex source and exact diagonal covariance; 64 physical marked recoil realizations common across bands, copies and D2 echoes',
                        angular_support='six accepted directions; off-axis interpolation is an explicit conditional model, not measured source support'))
            np.savez_compressed(path,**save)
            print(f'signed_source E={ec:g} closure={arm["closure"]} elapsed_s={time.monotonic()-start:.3f}',flush=True)
    if rows:table(out/'source_projection.csv',rows)
    return rows

@lru_cache(maxsize=64)
def guided_profile(path4, frequency, pol):
    s,native,fields,launch,response,acoustic,yielding=inherited(path4)
    guide=acoustic.guided(frequency,pol); k=guide['k']; degree=14; count=degree+1
    I=np.eye(3); R=launch.ROTATION
    C4=176.04e9*np.einsum('ij,kl->ijkl',I,I)+71.73e9*(np.einsum('ik,jl->ijkl',I,I)+np.einsum('il,jk->ijkl',I,I))
    C4+=(227.1e9-176.04e9-2*71.73e9)*np.einsum('ia,ja,ka,la->ijkl',R,R,R,R)
    pairs=((0,0),(1,1),(2,2),(1,2),(0,2),(0,1)); C=np.array([[C4[i,j,a,b] for a,b in pairs] for i,j in pairs])
    def basis(z):
        ph=np.array([legval(z,np.eye(count)[j])*np.sqrt((2*j+1)/(RHO*H)) for j in range(count)]).T
        dz=np.array([legval(z,legder(np.eye(count)[j]))*2/H*np.sqrt((2*j+1)/(RHO*H)) for j in range(count)]).T
        b=np.zeros((len(z),6,3*count),complex); kx=ky=k/np.sqrt(2)
        b[:,0,:count]=1j*kx*ph;b[:,1,count:2*count]=1j*ky*ph;b[:,2,2*count:]=dz
        b[:,3,count:2*count]=dz;b[:,3,2*count:]=1j*ky*ph;b[:,4,:count]=dz;b[:,4,2*count:]=1j*kx*ph
        b[:,5,:count]=1j*ky*ph;b[:,5,count:2*count]=1j*kx*ph
        return b
    z,w=leggauss(18);b=basis(z); matrix=np.einsum('zai,ab,zbj,z->ij',b.conj(),C,b,w*H/2,optimize=True)
    val,vec=eigh(matrix,subset_by_index=[0,17]);j=guide['selected_mode'];om=np.sqrt(max(val[j],0.))
    mode=vec[:,j].copy();pivot=int(np.argmax(abs(mode)));mode*=np.exp(-1j*np.angle(mode[pivot]))
    depths=np.unique(np.r_[0.,np.geomspace(1e-10,H,128),H]);coordinate=2*depths/H-1
    normal=np.sqrt((2*np.arange(count)+1)/(RHO*H))
    div_coeff=1j*k/np.sqrt(2)*(mode[:count]+mode[count:2*count])*normal
    div_coeff[:-1]+=2/H*legder(mode[2*count:]*normal)
    div_coeff/=om
    div=legval(coordinate,div_coeff)
    ph=np.array([legval(coordinate,np.eye(count)[m])*np.sqrt((2*m+1)/(RHO*H)) for m in range(count)]).T
    displacement=ph@mode.reshape(3,count).T
    return depths,div,displacement,guide['k'],guide['vp'],guide['vg'],div_coeff
 
 
def guided_loading(div_coeff,density):
    """Exact degree-28 modal energy overlap on every constant-D layer.

    Mass-normalized div(u)/omega has units s/sqrt(kg), hence h is Pa^-1.
    There is no extra foil-thickness divisor or second density multiplier.
    """
    edge=np.unique(np.r_[0.,density.edge[(density.edge>0)&(density.edge<H)],H])
    q,w=leggauss(15);half=np.diff(edge)/2
    z=(edge[:-1]+edge[1:])[:,None]/2+half[:,None]*q
    measure=half[:,None]*w*abs(legval(2*z/H-1,div_coeff))**2
    h=float(measure.sum())
    return h,float(np.sum(measure*density.local(z))/h) if h>0 else 0.

class Green:
    def __init__(self,bundle,path4,ec,closure='nominal_retained',allocation='onsite',temperature=300.,bandwidth_scale=1.,control='joined'):
        self.path4=str(Path(path4).resolve());self.ec=ec;self.closure=closure;self.allocation=allocation
        self.s,self.native,self.fields,self.launch,self.response,self.acoustic,self.yielding=inherited(path4)
        self.arm=next(a for a in self.s.ARMS if a['closure']==closure and a['allocation']==allocation)
        self.data=np.load(Path(bundle)/f'spectrum_E{ec:g}_{closure}.npz');self.density=self.acoustic.Density(ec)
        self.channels=self.launch.ports();self.temperature=temperature
        import csv
        shape=list(csv.DictReader((Path(bundle)/'spectral_shape.csv').open()))
        self.bandwidth=[]
        for ch in self.channels:
            rows=[r for r in shape if r['scenario']=='track' and int(r['species'])==0 and int(r['band'])==ch['band']]
            if rows and not ch['guided']:
                frequency=np.array([float(r['frequency_Hz']) for r in rows]);weight=np.array([float(r['energy_J']) for r in rows]);weight/=weight.sum()
                sigma=np.sqrt(float(weight@((frequency-float(weight@frequency))**2)))
            else:sigma=ch['frequency']/4
            self.bandwidth.append(2*np.pi*sigma*bandwidth_scale)
        if control not in ('joined','zero_real','zero_native'):raise ValueError('Unknown susceptibility control')
        self.material={};self.loading={};self.susceptibility={};self.guides={};self.guide_phase_scope={};self.guide_effective_density={}
        self.source_phase={(pi,sp):np.ones((3,6),complex) for pi in range(8) for sp in (0,1)}
        import json
        tracks=json.loads((Path(path4)/'frozen/tracks.json').read_text())
        self.tracks=tracks
        self.source_depth=[]
        self.source_extent=[]
        for sp in (0,1):
            tr=tracks[f'{ec:g}_{sp}'];z,t,dep,E,dt,dp=self.launch.track_arrays(tr)
            total=dep.sum(axis=1);terminal=tr['terminal_unresolved_eV']*QE
            # Historical field name: terminal_depth_A is already SI metres;
            # timeline depth_A rows are Angstroms and were converted above.
            self.source_depth.append(float((z@total+tr['terminal_depth_A']*terminal)/(total.sum()+terminal)))
            self.source_extent.append(float(tr['terminal_depth_A']+2e-9))
        for pi,ch in enumerate(self.channels):
            cov=self.response.covariance(self.arm,ch['frequency'],temperature)
            B=cov['bulk_per_D_J'];derivative=cov['bulk_derivative_per_D_J_s']
            if control=='zero_real':B=1j*B.imag;derivative=1j*derivative.imag
            if control=='zero_native':B=0j;derivative=0j
            self.susceptibility[pi]=(B,derivative,cov['availability_per_D_J'])
            for pol in range(3):
                if ch['guided']:
                    profile=guided_profile(self.path4,ch['frequency'],pol)
                    self.guides[pi,pol]=profile
                    zz,div,displacement,k,vp,vg,div_coeff=profile
                    for sp in (0,1):
                        anchor,scope=self.guide_source_phase(pi,pol,sp,profile)
                        self.source_phase[pi,sp][pol]=anchor;self.guide_phase_scope[pi,pol,sp]=scope
                    hr,loading_density=guided_loading(div_coeff,self.density)
                    self.guide_effective_density[pi,pol]=loading_density
                    with (Path(path4)/'frozen/guided_branches.csv').open() as stream:
                        branch=[r for r in csv.DictReader(stream) if int(r['branch'])==pol]
                    kk=np.array([float(r['k_parallel_m_1']) for r in branch]);ff=np.array([float(r['frequency_Hz']) for r in branch])
                    order=np.argsort(kk);kk=kk[order];ff=ff[order];vv=np.gradient(2*np.pi*ff,kk)
                    vglog=float(np.interp(ch['frequency'],ff,np.gradient(vv,ff)))*ch['frequency']/vg
                    hlog=0. # retained guided eigenprofile is fixed within this port
                else:
                    vp=float(np.mean(2*np.pi*ch['frequency']/self.data[f'k_0_{pi}'][pol]))
                    cap=(3.4e12,4.7e12,6.8e12)[pol];x=ch['frequency']/cap
                    speed=vp*np.arcsin(x)/x;vg=speed*np.sqrt(1-x*x)
                    vglog=-x*x/(1-x*x);hlog=-2*(1-vp/vg)
                    hr=float(np.mean(np.sum(self.data[f'n_0_{pi}']*self.data[f'polar_0_{pi}'][pol],axis=1)**2))/(RHO*vp*vp)
                    loading_density=float(self.density.local(0.))
                r=self.response.coefficients(self.arm,ch['frequency'],loading_density,vp,vg,hr,temperature)
                intrinsic=vg*np.log(10)/10*32*ch['frequency']/1e6 if self.arm['loss']=='Cu' else 1/50e-6+vg/.2
                self.material[pi,pol]=(vp,vg,intrinsic,r,hr)
                omega=2*np.pi*ch['frequency']
                phase=-.5*omega*hr*B.real
                group=-.5*hr*(B.real*(1+hlog-vglog)+omega*derivative.real)
                loss=-omega*hr*B.imag
                if loss<0:raise ValueError('Nonpassive native attenuation')
                self.loading[pi,pol]=(phase,group,loss)
    def guide_source_phase(self,pi,pol,sp,profile):
        """Gauge-covariant known normal impulse; residual force phase is a scenario."""
        zz,div,displacement,k,vp,vg,div_coeff=profile
        tr=self.tracks[f'{self.ec:g}_{sp}'];z,t,dep,E,dt,dp=self.launch.track_arrays(tr)
        omega=2*np.pi*self.channels[pi]['frequency']
        source_z=np.interp(z,zz,displacement[:,2].real)+1j*np.interp(z,zz,displacement[:,2].imag)
        impulse=np.sum(np.conj(source_z)*dp*np.exp(1j*omega*t)*np.sinc(omega*dt/(2*np.pi)))
        source_value=np.array([np.interp(self.source_depth[sp],zz,displacement[:,axis].real)+1j*np.interp(self.source_depth[sp],zz,displacement[:,axis].imag) for axis in range(3)])
        reference=np.conj(source_value[int(np.argmax(abs(source_value)))])
        if abs(reference)==0:raise ValueError('Zero guided source eigenvector')
        phase=np.full(6,reference/abs(reference),complex);known=0
        scale=float(np.max(abs(displacement))*np.sum(dp))
        nn=self.data[f'n_{sp}_{pi}'];polar=self.data[f'polar_{sp}_{pi}'][pol]
        boundaries=np.r_[0.,[r['depth_A']*1e-10 for r in tr['timeline']]]
        for j,n in enumerate(nn):
            source_k=self.data[f'k_{sp}_{pi}'][pol,j]
            bulk=polar[j,2]*np.sum(dp*np.exp(1j*omega*t-1j*source_k*n[2]*z)*np.sinc((omega*dt-source_k*n[2]*np.diff(boundaries))/(2*np.pi)))
            if abs(impulse)>1e-10*scale and abs(bulk)>1e-12*np.sum(dp):
                phase[j]=impulse/abs(impulse)*np.conj(bulk)/abs(bulk);known+=int(self.data[f'weight_{sp}_{pi}'][j]>0)
        return phase,dict(known_normal_impulse_angular_nodes=known,
            residual='same modal gauge, but residual pressure/recoil/LTE force phase remains an explicitly uncalibrated scenario')

    def _angular(self,pi,sp,n):
        nn=self.data[f'n_{sp}_{pi}'];support=np.flatnonzero(self.data[f'weight_{sp}_{pi}']>0);dots=n@nn[support].T
        # Barycentric nearest three angular directions; signed complex values,
        # variance remains a second moment and is NEVER assigned a fixed phase.
        local=np.argpartition(dots,-3,axis=-1)[...,-3:];dist=np.maximum(1-np.take_along_axis(dots,local,axis=-1),1e-12);ids=support[local]
        weight=1/dist;weight/=weight.sum(axis=-1,keepdims=True)
        a=self.data[f'a_{sp}_{pi}']
        coherent=(a[...,0]+a[...,2]+a[...,3])*self.source_phase[pi,sp]
        aa=np.sum(np.moveaxis(coherent[:,ids],0,-2)*weight[...,None,:],axis=-1)
        variance=np.sum(np.moveaxis(abs(a[:,ids,1])**2,0,-2)*weight[...,None,:],axis=-1)
        return aa,variance
    def components(self,receiver,source_xy,sp,recoil_mark=None,returns=True,max_port=7,min_port=0):
        """Retarded components; one reusable law for shots and Cox convolution."""
        receiver=np.asarray(receiver);source_xy=np.asarray(source_xy)
        shape=np.broadcast_shapes(source_xy.shape[:-1],receiver.shape[:-1])
        src=np.zeros(shape+(3,));src[...,:2]=source_xy
        src[...,2]=self.source_depth[sp]
        images=[(src,1.,0)]
        if returns:
            # Free-face normal volumetric strain reverses under mirror. Edge
            # amplitude sqrt(.8) is the accepted passive energy reflection.
            for axis,length,origin in ((0,L,-L/2),(1,L,-L/2)):
                for wall in (origin,origin+length):
                    image=src.copy();image[...,axis]=2*wall-src[...,axis]
                    images.append((image,-np.sqrt(.8),1))
            for wall in (0.,H):
                image=src.copy();image[...,2]=2*wall-src[...,2]
                images.append((image,-1.,2))
        for pi,ch in enumerate(self.channels):
            if pi>max_port or pi<min_port:continue
            frequency=ch['frequency'];omega=2*np.pi*frequency
            bandwidth=self.bandwidth[pi]
            for image,reflection,kind in images:
                if ch['guided'] and kind==2:continue # free faces already in exact guided mode
                delta=receiver-image
                if ch['guided']:delta=delta.copy();delta[...,2]=0.
                distance=np.linalg.norm(delta,axis=-1)
                core=2e-9;rr=np.sqrt(distance*distance+core*core)
                n=delta/rr[...,None]
                aa,var=self._angular(pi,sp,n)
                for pol in range(3):
                    vp,vg,intrinsic,r,hr=self.material[pi,pol]
                    if ch['guided']:
                        zz,div,displacement,k,vp,vg,div_coeff=self.guides[pi,pol]
                        profile=legval(2*receiver[...,2]/H-1,div_coeff)
                        metric=np.exp(.25j*np.pi)*profile/np.sqrt(2*np.pi*rr*vg)
                        launch_delay=receiver[...,2]/4460.
                    else:
                        k=omega/vp
                        nearest=np.argmax(n@self.data[f'n_{sp}_{pi}'].T,axis=-1)
                        vectors=self.data[f'polar_{sp}_{pi}'][pol,nearest]
                        projection=np.sum(n*vectors,axis=-1)
                        metric=1j*projection/(np.sqrt(4*np.pi*rr*rr*vg*RHO)*vp)
                        launch_delay=0.
                    host_time=distance/vg
                    if ch['guided']:
                        integrated_density=self.guide_effective_density[pi,pol]*host_time
                    else:
                        vz=np.divide(receiver[...,2]-image[...,2],host_time,out=np.zeros_like(host_time),where=host_time>0)
                        start=image[...,2]
                        # Vectorized primitive is the exact reflected depth-loading
                        # line integral. No homogeneous birth-density hard cutoff.
                        integrated_density=np.divide(self.density.primitive(start+vz*host_time)[0]-self.density.primitive(start)[0],vz,
                            out=np.asarray(self.density.local(start)*host_time),where=abs(vz)>1e-12)
                    phase_per_D,group_per_D,loss_per_D=self.loading[pi,pol]
                    arrival=host_time+launch_delay+group_per_D*integrated_density
                    if np.any(arrival<0):raise ValueError('Noncausal retained group-delay projection')
                    damping=np.exp(-np.minimum((intrinsic*host_time+loss_per_D*integrated_density)/2,745.))
                    phase=k*distance-omega*arrival+phase_per_D*integrated_density
                    amplitude=aa[...,pol]
                    if recoil_mark is not None:
                        mark=np.asarray(recoil_mark,dtype=int)%64
                        nearest=np.argmax(n@self.data[f'n_{sp}_{pi}'].T,axis=-1)
                        amplitude=amplitude+self.data[f'stochastic_{sp}_{pi}'][pol,nearest,mark]*self.source_phase[pi,sp][pol,nearest]
                    amplitude=amplitude*np.exp(1j*k*n[...,2]*self.source_depth[sp])
                    coefficient=reflection*np.sqrt(2.)*metric*damping*amplitude*np.exp(1j*phase)
                    yield pi,omega,bandwidth,arrival,coefficient
    def parent(self,receiver,source_xy,age,sp,recoil_mark,returns=True,max_port=7,min_port=0):
        """Total same-parent contribution, before any powers or cumulants."""
        age=np.asarray(age);shape=np.broadcast_shapes(np.shape(source_xy)[:-1],age.shape)
        per=np.zeros(shape+(8,),complex);intensity=np.zeros(shape)
        for pi,omega,bandwidth,arrival,coefficient in self.components(receiver,source_xy,sp,recoil_mark,returns,max_port,min_port):
            elapsed=age-arrival
            field=coefficient*wavelet(elapsed,bandwidth,omega)*np.exp(-1j*omega*elapsed)
            per[...,pi]+=field;intensity+=abs(field)**2
        return per.sum(axis=-1),intensity,per

    def mean(self,receiver,beam,times,calendar,returns=True,per_band=False):
        """Analytic causal Cox convolution at fixed or batched real receivers."""
        def primitive(T,a,omega,nu):
            duration=packet_duration(a,omega);t=np.clip(T,0.,duration);beta=2*np.pi/duration
            def integral(v):
                return t*np.exp(-.5j*v*t)*np.sinc(v*t/(2*np.pi))
            norm=np.sqrt(3*omega*omega*duration/8+np.pi*np.pi/(2*duration))
            value=(omega/2*integral(nu)-(omega-beta)/4*integral(nu-beta)-(omega+beta)/4*integral(nu+beta))/norm
            return np.where(T>=duration,0.,value) if nu==omega else value
        q,w=leggauss(12);x=.005*q
        xx,yy=np.meshgrid(x,x,indexing='ij');xy=np.c_[xx.ravel(),yy.ravel()]
        spatial=(.005**2*np.outer(w,w)).ravel()*self.s.beam(xy[:,0],xy[:,1],beam)/self.s.beam_moments(beam)[0]
        times=np.atleast_1d(times);answer=np.zeros((len(times),8));receiver=np.asarray(receiver)
        source=xy[None,:,:] if receiver.ndim>1 else xy
        target=receiver[:,None,:] if receiver.ndim>1 else receiver
        for sp,fraction in ((0,.85),(1,.15)):
            for pi,omega,a,arrival,c in self.components(target,source,sp,None,returns):
                T=times[:,None]-arrival;duration=packet_duration(a,omega)
                value=primitive(T,a,omega,omega)
                if beam=='hot_bursts':
                    angle=calendar.omega*T+calendar.phase
                    value+=.35*(np.exp(1j*angle)*primitive(T,a,omega,omega+calendar.omega)+np.exp(-1j*angle)*primitive(T,a,omega,omega-calendar.omega))
                    first=np.maximum(np.searchsorted(calendar.edges,np.maximum(0.,T-duration),side='right')-1,0)
                    last=np.minimum(np.searchsorted(calendar.edges,np.maximum(0.,T),side='right'),len(calendar.state))
                    for offset in range(int(np.max(last-first,initial=0))):
                        index=np.minimum(first+offset,len(calendar.state)-1)
                        on=(index<last)&calendar.state[index]
                        value+=9*on*(primitive(T-calendar.edges[index],a,omega,omega)-primitive(T-calendar.edges[index+1],a,omega,omega))
                    value/=calendar.norm
                answer[:,pi]+=fraction*self.s.current(self.ec)/QE*(np.real(c*value)@spatial)
        return answer if per_band else answer.sum(axis=1)


def receiver_nodes():
    base=np.array([(x,y,z) for x in np.linspace(-L/2,L/2,18) for y in np.linspace(-L/2,L/2,18) for z in (0.,H/2,H)])
    # Exact coordinates, not a propagation patch: all receivers call the same Green.
    nested=[]
    for center in ((.0025,-.0025,1e-8),(-.0025,.0025,1e-8),(0.,0.,1e-8)):
        for scale in (1.25e-3,1.25e-4,1.25e-5,1.25e-6,1.25e-7,1.25e-8,1.25e-9,1.25e-10):
            for sign in (-1,0,1):
                for axis in range(3):
                    point=np.array(center);point[axis]+=sign*scale
                    if 0<=point[2]<=H and np.max(abs(point[:2]))<=L/2:nested.append(point)
    return np.unique(np.r_[base,np.array(nested)],axis=0)
