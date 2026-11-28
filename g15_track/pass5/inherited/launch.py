"""Band source operator, not a G13 producer: same tracks, native pressure,
actual onward heat Green, and a finite signed mechanical-work decomposition.
The retained angular rule is six directions; accepted G13 band integrals fix
only its transfer normalization. No sound fraction is fitted to beam power.
"""
import csv,math
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.special import j0
import sources as s
import native
MPD=106.4*1.66053906660e-27;ND=s.RHO/MPD
DIRECTIONS=np.r_[np.eye(3),-np.eye(3)]
ROTATION=Rotation.random(64,random_state=173).as_matrix()[27]

def ports():
    """Band 1 explicitly straddles the accepted 10-MHz guided/bulk handover."""
    result=[]
    for band in range(7):
        lo=2.5e5*10**band;hi=10*lo
        if band<2:result.append(dict(band=band,guided=True,frequency=math.sqrt(lo*min(hi,1e7))))
        if hi>1e7:result.append(dict(band=band,guided=False,frequency=math.sqrt(max(lo,1e7)*hi)))
    return result


def reference_ports(sp):
    reference=np.load(s.HERE/'frozen/acoustic_source.npz')['energy_J'][0,sp]
    with (s.HERE/'frozen/acoustic_frequency.csv').open() as stream:
        rows=[r for r in csv.DictReader(stream) if r['scenario']=='track' and int(r['species'])==sp and int(r['band'])==1]
    bulk=sum(float(r['energy_J']) for r in rows);fraction=1-bulk/reference[1].sum()
    if not 0<=fraction<=1:raise ValueError('Accepted guided/bulk source partition is inconsistent')
    # The frozen file preserves the exact band-total partition, not its
    # per-polarization split. Retain each recorded bucket's relative share.
    return np.array([reference[r['band']]*(fraction if r['guided'] else 1-fraction)
        if r['band']==1 else reference[r['band']] for r in ports()])


def eigenbasis(n):
    lam=176.04e9;mu=71.73e9;anisot=227.1e9-lam-2*mu
    projections=n@ROTATION
    tensor=mu*np.eye(3)+(lam+mu)*np.outer(n,n)+anisot*(ROTATION*(projections*projections))@ROTATION.T
    values,vectors=np.linalg.eigh(tensor)
    for j in range(3):
        i=int(np.argmax(abs(vectors[:,j])))
        if vectors[i,j]<0:vectors[:,j]*=-1
    speeds=np.sqrt(values/s.RHO)
    # d omega/d k from the same Christoffel tensor, not the phase velocity.
    group=[]
    for c,e in zip(speeds,vectors.T):
        elastic=mu*n+(lam+mu)*float(e@n)*e+anisot*np.sum(ROTATION*((e@ROTATION)**2*projections)[None,:],axis=1)
        group.append(elastic/(s.RHO*c))
    return speeds,vectors,np.array(group)


def track_arrays(tr):
    z=np.r_[0.,[r['depth_A']*1e-10 for r in tr['timeline']]]
    t=np.r_[0.,[r['time_fs']*1e-15 for r in tr['timeline']]]
    dep=np.diff(np.array([[0.,0.,0.]]+[[r[k] for k in ('electron_eV','host_recoil_eV','D_recoil_eV')] for r in tr['timeline']]),axis=0)*s.QE
    E=np.r_[tr['initial_eV'],[r['energy_eV'] for r in tr['timeline']]]
    momentum=np.sqrt(2*s.MD*E*s.QE)
    return (z[:-1]+z[1:])/2,(t[:-1]+t[1:])/2,dep,(E[:-1]+E[1:])/2,np.diff(t),-np.diff(momentum)


def recoil_marks(tr):
    """Accepted ZBL impact-area/azimuth quadrature; moving kinetic energy once."""
    bca=s.load_file(s.HERE/'frozen/bca.py','launch_bca');z,t,dep,E,dt,dp=track_arrays(tr)
    x,w=np.polynomial.legendre.leggauss(8);cell=ND**(-1/3);pmax=1/np.sqrt(np.pi*ND*1e-27*cell*1e9)
    a=.8854*bca.A0/(1+46**.23);eps=a*E*106.4/(bca.M_D+106.4)/(46*bca.E2)
    theta=bca.cm_angle(np.repeat(eps,len(x)),np.tile(pmax*np.sqrt((x+1)/2)/a,len(E))).reshape(len(E),-1)
    Tmax=4*s.MD*MPD/(s.MD+MPD)**2*E;T=np.repeat(Tmax[:,None]*np.sin(theta/2)**2,8,axis=1)
    prob=np.tile(np.repeat(w/2,8)/8,(len(E),1));cz=np.sqrt(T/Tmax[:,None]);st=np.sqrt(1-cz*cz);az=2*np.pi*(np.arange(8)+.5)/8
    phi=np.tile(az,8)[None,:];direction=np.stack((st*np.cos(phi),st*np.sin(phi),cz),axis=-1)
    q=np.sqrt(2*MPD*T*s.QE)[...,None]*direction;intensity=prob*(dep[:,1]/s.QE/(prob*T).sum(axis=1))[:,None]
    xx,ww=np.polynomial.legendre.leggauss(24);lo=np.log(.1);hi=np.log(np.maximum(T.ravel(),.100001))
    energy=np.exp(lo+(xx+1)[:,None]/2*(hi-lo));stopping=bca.sn_zbl(energy,46,106.4,Z1=46,M1=106.4)
    distance=(ww[:,None]/2*(hi-lo)*energy/(ND*1e-27*stopping)).sum(axis=0)*1e-9
    speed=np.sqrt(2*T.ravel()*s.QE/MPD);duration=np.divide(2*distance,speed,out=np.zeros_like(distance),where=speed>0)
    return q.reshape(-1,3),intensity.ravel(),(direction.reshape(-1,3)*distance[:,None]),duration,T.shape[1]


def source_amplitudes(tr,sp,frequency,n,pol,heat_green,pressure_integral,marks):
    c,ev,group=eigenbasis(n);cap=(3.4e12,4.7e12,6.8e12)[pol];f=frequency
    k=2*np.pi*cap/c[pol]*math.asin(f/cap);omega=2*np.pi*f;polarization=ev[:,pol]
    z,t,dep,E,dt,dp=track_arrays(tr);delay=.74e-10/math.sqrt(2*tr['initial_eV']*s.QE/s.MD) if sp else 0.
    parent=(1+np.exp(1j*omega*delay)) if sp else 1.
    boundaries=np.r_[0.,[r['depth_A']*1e-10 for r in tr['timeline']]]
    phase=np.exp(1j*omega*t-1j*k*n[2]*z);smooth=np.sinc((omega*dt-k*n[2]*np.diff(boundaries))/(2*np.pi))
    phase*=smooth
    mechanical=polarization[2]*np.sum(phase*dp)*parent
    q,lam,displ,duration,number=marks;qproj=q@polarization
    pair=np.expm1(1j*omega*duration-1j*k*(displ@n));jump=qproj*pair*np.repeat(phase,number)
    mechanical+=np.sum(lam*jump)*parent;variance=np.sum(lam*abs(jump)**2)*(sp+1)
    totals=np.zeros(4,complex)
    totals[1:]=np.array([np.sum(phase*dep[:,j]) for j in (1,2,0)])*parent
    totals[0]=tr['terminal_unresolved_eV']*s.QE*np.exp(1j*omega*tr['terminal_time_fs']*1e-15-1j*k*n[2]*tr['terminal_depth_A'])*parent
    late=0j
    for origin in range(4):
        rows=heat_green[heat_green[:,0]==origin]
        transfer=np.sum(rows[:,6]*np.exp(1j*omega*rows[:,2]-1j*k*(rows[:,3:6]@n)))
        late+=totals[origin]*transfer
    gamma_l=3*s.BULK*s.ALPHA/s.CV
    thermal=-1j*k*float(n@polarization)*gamma_l*late/(s.KTH/s.CV*k*k-1j*omega)
    electronic=-1j*k*float(n@polarization)*pressure_integral(omega,k,n)
    # The 2-nm source form factor belongs to the impact/collision source;
    # native pressure already includes its explicit transported spatial shape.
    form=np.exp(-.5*(k*2e-9)**2)
    oldheat=np.sum(phase*(dep[:,0]+dep[:,2]))*parent+totals[0]
    nodes,node_weights=np.polynomial.legendre.leggauss(3);recoil_energy=np.sum(q*q,axis=1)/(2*MPD)
    for fraction,weight in zip((nodes+1)/2,node_weights/2):
        shifted=np.exp(1j*omega*duration*(1-np.sqrt(1-fraction))-1j*k*(displ@n)*fraction)
        oldheat+=weight*np.sum(np.repeat(phase,number)*shifted*lam*recoil_energy)*parent
    oldthermal=-1j*k*float(n@polarization)*gamma_l*oldheat/(s.KTH/s.CV*k*k-1j*omega)
    return mechanical*form,variance*form*form,thermal*form,electronic*form,oldthermal*form,k,float(np.linalg.norm(group[pol])*np.sqrt(1-(f/cap)**2)),polarization


def spectral_source(arm,ec,sp,tr,basis,endpoints):
    positions,volume,weights,times,geometry,tt,h,length,population,secondary_shape,error,u=basis
    material=native.Material(16);pressure=material.value(u,'P')
    nonlinear=pressure-material.p0*u
    # Near-front nonlinear correction plus an exactly normalized linear carrier
    # moment. The finite aperture error and missing q calibration remain visible.
    electron=population[:,3]+population[:,8]
    # The analytic linear Green uses the causal flight/scatter second moment;
    # its finite angular closure is the same one used by the global RTE.
    age=np.minimum(times*1e-15,arm['ne']);width=np.minimum(1e6*times*1e-15,np.sqrt(2*1e6*arm['lambda_e']*age))
    zcenter=np.interp(times,np.array([r['time_fs'] for r in tr['timeline']]),np.array([r['depth_A']*1e-10 for r in tr['timeline']]))/2
    def integral(omega,k,n):
        spatial=j0(k*np.sqrt(1-n[2]*n[2])*positions[:,0]*1e-9)*np.exp(-1j*k*n[2]*positions[:,1]*1e-9)
        local=nonlinear@(spatial*volume*1e-27)
        analytic=material.p0*native.V0*electron*native.ball_transform(k*width)*np.exp(-1j*k*n[2]*zcenter)
        return np.trapezoid((local+analytic)*np.exp(1j*omega*times*1e-15),times*1e-15)
    marks=recoil_marks(tr)
    channels=ports();result=np.zeros((len(channels),3,6,5),complex)
    velocities=np.zeros((len(channels),3,6));polarizations=np.zeros((len(channels),3,6,3))
    for port,channel in enumerate(channels):
        f=channel['frequency']
        for pol in range(3):
            for direction,n in enumerate(DIRECTIONS):
                a,v,thermal,electronic,oldthermal,k,vg,polarization=source_amplitudes(tr,sp,f,n,pol,endpoints,integral,marks)
                result[port,pol,direction]=a,math.sqrt(v),thermal,electronic,oldthermal
                velocities[port,pol,direction]=vg;polarizations[port,pol,direction]=polarization
    return result,velocities,polarizations


def calibration(spectrum,sp):
    reference=reference_ports(sp)
    norm=abs(spectrum[...,0]+spectrum[...,4])**2+abs(spectrum[...,1])**2
    # Guided modes propagate in the foil plane. No fictitious low-MHz vertical
    # branch is launched; its depth eigenprofile is a separate retained mode.
    guided=np.array([r['guided'] for r in ports()]);norm[guided,:,2]=0.;norm[guided,:,5]=0.
    denominator=norm.sum(axis=2)
    if np.any(denominator<=0):raise ValueError('Missing supported source normalization')
    return reference/denominator


def project(spectrum,normalization):
    a=spectrum[...,0];variance=abs(spectrum[...,1])**2;thermal=spectrum[...,2];electronic=spectrum[...,3]
    total=a+thermal+electronic;factor=np.broadcast_to(normalization[:,:,None],total.shape).copy()
    guided=np.array([r['guided'] for r in ports()]);factor[guided,:,2]=0.;factor[guided,:,5]=0.
    energy=factor*(abs(total)**2+variance)
    work_mechanical=factor*(np.real(np.conj(a)*total)+variance)
    work_heat=factor*np.real(np.conj(thermal)*total);work_electron=factor*np.real(np.conj(electronic)*total)
    defect=float(np.max(abs(energy-work_mechanical-work_heat-work_electron)))
    return energy,work_mechanical,work_heat,work_electron,defect
