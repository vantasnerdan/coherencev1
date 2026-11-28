"""Whole-foil retarded band characteristics with conservative timed loss ports.

A positive lifetime quadrature, not instantaneous stocks, is the computational
state. Each electrical parent drives every ray through one shared source clock.
The inherited attenuation laws are *intrinsic-floor assumptions* in the
explicitly decomposed scenario; the coupling-off control is their total-loss
reference. They are never added to an asserted measured PdD total.
"""
import csv,math
from functools import lru_cache
import numpy as np
from scipy.linalg import eigh
from numpy.polynomial.legendre import leggauss,legval,legder
from scipy.signal import fftconvolve
import sources as s
import fields as f
import launch,native,response
from retarded import DelayMeasure


class Density:
    def __init__(self,ec):
        profile=next(r for r in s.inputs()['profiles'] if r['case']['E_cm_keV']==ec)
        a=np.array(profile['layers']);self.edge=np.r_[0.,a[:,1]];self.nd=a[:,3]
        self.F=np.r_[0.,np.cumsum(self.nd*np.diff(self.edge))]
        self.G=np.r_[0.,np.cumsum(self.nd*np.diff(self.edge**2)/2)]
    def local(self,z):
        r=np.mod(z,2*s.H);r=np.where(r<=s.H,r,2*s.H-r)
        return self.nd[np.minimum(np.searchsorted(self.edge[1:],r,side='right'),len(self.nd)-1)]
    def primitive(self,z):
        z=np.asarray(z);period=2*s.H;q=np.floor(z/period);r=z-q*period;u=np.where(r<=s.H,r,period-r)
        j=np.minimum(np.searchsorted(self.edge[1:],u,side='right'),len(self.nd)-1)
        F=self.F[j]+self.nd[j]*(u-self.edge[j]);G=self.G[j]+self.nd[j]*(u*u-self.edge[j]**2)/2
        Fp=np.where(r<=s.H,F,2*self.F[-1]-F)
        Gp=np.where(r<=s.H,G,period*self.F[-1]-period*F+G)
        norm=2*self.F[-1];moment=s.H*norm
        return q*norm+Fp,q*moment+period*norm*q*(q-1)/2+q*period*Fp+Gp
    def flight(self,z,v,t):
        if abs(v)<1e-12:return self.local(z)*t
        return (self.primitive(z+v*t)[0]-self.primitive(z)[0])/v


@lru_cache(maxsize=32)
def guided(frequency,pol):
    with (s.HERE/'frozen/guided_branches.csv').open() as stream:
        rows=[r for r in csv.DictReader(stream) if int(r['branch'])==pol]
    kval=np.array([float(r['k_parallel_m_1']) for r in rows]);freq=np.array([float(r['frequency_Hz']) for r in rows])
    order=np.argsort(kval);kval=kval[order];freq=freq[order]
    k=float(np.interp(frequency,freq,kval));vg=float(np.interp(frequency,freq,np.gradient(2*np.pi*freq,kval)))
    # Exact accepted free-face Legendre Galerkin form, selected grain 27,
    # representative in-plane direction. This is not a G13 producer/sweep.
    degree=14;n=degree+1;theta=np.pi/4;R=launch.ROTATION
    I=np.eye(3);C4=176.04e9*np.einsum('ij,kl->ijkl',I,I)+71.73e9*(np.einsum('ik,jl->ijkl',I,I)+np.einsum('il,jk->ijkl',I,I))
    C4+=(227.1e9-176.04e9-2*71.73e9)*np.einsum('ia,ja,ka,la->ijkl',R,R,R,R)
    pairs=((0,0),(1,1),(2,2),(1,2),(0,2),(0,1));C6=np.array([[C4[i,j,a,b] for a,b in pairs] for i,j in pairs])
    def basis(z):
        ph=np.array([legval(z,np.eye(n)[j])*np.sqrt((2*j+1)/(s.RHO*s.H)) for j in range(n)]).T
        dz=np.array([legval(z,legder(np.eye(n)[j]))*2/s.H*np.sqrt((2*j+1)/(s.RHO*s.H)) for j in range(n)]).T
        b=np.zeros((len(z),6,3*n),complex);kx=k*np.cos(theta);ky=k*np.sin(theta)
        b[:,0,:n]=1j*kx*ph;b[:,1,n:2*n]=1j*ky*ph;b[:,2,2*n:]=dz
        b[:,3,n:2*n]=dz;b[:,3,2*n:]=1j*ky*ph;b[:,4,:n]=dz;b[:,4,2*n:]=1j*kx*ph
        b[:,5,:n]=1j*ky*ph;b[:,5,n:2*n]=1j*kx*ph
        return ph,b
    z,w=leggauss(degree+3);ph,b=basis(z)
    matrix=np.einsum('zai,ab,zbj,z->ij',b.conj(),C6,b,w*s.H/2,optimize=True)
    val,vec=eigh(matrix,subset_by_index=[0,17]);omegas=np.sqrt(np.maximum(val,0.))
    j=pol if pol<2 else 2+int(np.argmin(abs(omegas[2:]-2*np.pi*frequency)))
    vector=vec[:,j];om=omegas[j];energy=[];hydro=[];centers=[];vz=[];sign=[]
    for lo,hi in zip(f.ZE[:-1],f.ZE[1:]):
        zz=(lo+hi)/2+(hi-lo)*z/2;pp,bb=basis(2*zz/s.H-1)
        U=np.column_stack([pp@vector[a*n:(a+1)*n] for a in range(3)]);strain=bb@vector
        metric=s.RHO*om*om*np.sum(abs(U)**2,axis=1)+np.einsum('za,ab,zb->z',strain.conj(),C6,strain).real
        ww=w*(hi-lo)/2;div=strain[:,:3].sum(axis=1);den=float(metric@ww)
        energy.append(den);centers.append(float((zz*metric)@ww/den));hydro.append(float(2*(abs(div)**2@ww)/den))
        vz.append(float(2*om*om*(abs(U[:,2])**2@ww)/den));sign.append(float(np.angle(div@ww)))
    energy=np.array(energy);energy/=energy.sum()
    return dict(k=k,vp=2*np.pi*frequency/k,vg=vg,profile=energy,z=np.array(centers),
        hydro=np.array(hydro),velocity_z=np.array(vz),hydro_phase=np.array(sign),selected_mode=j,
        profile_pole_Hz=float(om/(2*np.pi)))


def source_points():
    cells=np.flatnonzero(s.beam_cells(np.linspace(-s.L/2,s.L/2,f.N+1),'even_steady').ravel())
    points=[];ids=[]
    for cell in cells:
        ix,iy=divmod(int(cell),f.N)
        lx=max(-.005,-s.L/2+ix*f.DX);hx=min(.005,-s.L/2+(ix+1)*f.DX)
        ly=max(-.005,-s.L/2+iy*f.DX);hy=min(.005,-s.L/2+(iy+1)*f.DX)
        for a in (-1/np.sqrt(3),1/np.sqrt(3)):
            for b in (-1/np.sqrt(3),1/np.sqrt(3)):
                points.append(((lx+hx+a*(hx-lx))/2,(ly+hy+b*(hy-ly))/2));ids.append(cell)
    return np.array(points),np.array(ids)


class Operator:
    def __init__(self,arm,ec,beam,tracks,control='joined',temperature=300.,linear_strain=0.):
        self.arm=arm;self.ec=ec;self.beam=beam;self.control=control;self.rows=[];self.volume=[];self.boundary=[];self.rays=[];self.heat=[];self.bath=[]
        self.points,self.cells=source_points();self.weights=s.beam_cells(np.linspace(-s.L/2,s.L/2,f.N+1),beam).ravel()[self.cells]/4
        self.density=Density(ec);self.coupled=control not in ('strain_off','frozen','inherited_total')
        source=[np.load(s.HERE/f"source_E{ec:g}_{arm['closure']}_sp{sp}.npz") for sp in (0,1)]
        self.energy=np.stack([a['energy_J'] for a in source]);self.total_per_parent=np.einsum('s,spqd->',np.array([.85,.15]),self.energy)
        p=s.inputs()['lattice'];self.temperature=temperature;self.linear_strain=linear_strain
        self.response={a['frequency']:response.covariance(arm,a['frequency'],temperature,linear_strain,
            frozen=control=='frozen_G3_reference') for a in launch.ports()}
        for port,ch in enumerate(launch.ports()):
            frequency=ch['frequency'];w=2*np.pi*frequency;r=self.response[frequency]
            for pol in range(3):
                guide=guided(frequency,pol) if ch['guided'] else None
                for di,n in enumerate(launch.DIRECTIONS):
                    energies=self.energy[:,port,pol,di]
                    if not np.any(energies):continue
                    if guide:
                        # Normal launch marks couple into guided in-plane motion,
                        # not an unsupported stationary axial guided ray.
                        if abs(n[2])>.5:n=np.array([1.,1.,0.])*np.sign(n[2])/np.sqrt(2)
                        k=guide['k'];vg=guide['vg'];vp=guide['vp'];velocity=vg*n
                        depths=guide['z'];fractions=guide['profile'];hydro=guide['hydro'];velocity_metric=guide['velocity_z']
                    else:
                        speeds,polarization,group=launch.eigenbasis(n);cap=(3.4e12,4.7e12,6.8e12)[pol]
                        k=2*np.pi*cap/speeds[pol]*math.asin(frequency/cap);vp=w/k
                        velocity=group[pol]*np.sqrt(1-(frequency/cap)**2);vg=float(np.linalg.norm(velocity))
                        # Actual source energy-weighted depth, not a uniform 250-nm stock.
                        zz=[]
                        for sp in (0,1):
                            tr=tracks[f'{ec:g}_{sp}'];z,t,dep,_,_,_=launch.track_arrays(tr)
                            zz.append(float((z@dep.sum(axis=1)+tr['terminal_depth_A']*tr['terminal_unresolved_eV']*s.QE)/(dep.sum()+tr['terminal_unresolved_eV']*s.QE)))
                        depths=np.array([float(np.dot([.85,.15],zz))]);fractions=np.ones(1)
                        hydro=np.array([float(n@polarization[:,pol])**2/(s.RHO*vp*vp)])
                        velocity_metric=np.array([float(polarization[2,pol])**2/s.RHO])
                    intrinsic=vg*np.log(10)/10*32*frequency/1e6 if arm['loss']=='Cu' else 1/50e-6+vg/.2
                    for depth,z0 in enumerate(depths):
                        hr=hydro[depth];beta=float(-w*r['bulk_per_D_J'].imag*hr) if self.coupled else 0.
                        if beta<0:raise ValueError('Nonpassive optical loading')
                        correction=float(r['bulk_per_D_J'].real*hr) if self.coupled else 0.
                        slope=0. if guide else 2*(1-vg/vp)
                        transport=hr*(r['bulk_per_D_J'].real*(1+.5*vp/vg*slope)-2*r['availability_per_D_J']) if self.coupled else 0.
                        delay=z0/4460. if guide else 0.
                        nd0=float(self.density.local(z0));minimum=min(float(self.density.nd.min()),nd0)
                        horizon=28/(intrinsic+beta*minimum)+delay
                        edges=np.unique(np.r_[0.,delay,np.geomspace(max(1e-18,horizon*1e-8),horizon,48)])
                        def geometry(age,velocity=velocity.copy(),delay=delay,z0=float(z0),guide=bool(guide)):
                            a=np.asarray(age);active=np.maximum(a-delay,0.)
                            xy=self.points[None,:,:]+active[...,None,None]*velocity[:2]
                            z=np.minimum(4460*a,z0) if guide else z0+velocity[2]*a
                            raw=(xy+s.L/2)/s.L;counts=np.abs(np.floor(raw)).astype(int)
                            fold=np.mod(xy+s.L/2,2*s.L);xy=np.where(fold<=s.L,fold,2*s.L-fold)-s.L/2
                            zfold=np.mod(z,2*s.H);zfold=np.where(zfold<=s.H,zfold,2*s.H-zfold)
                            counts=np.concatenate((counts,np.broadcast_to(np.abs(np.floor(z/s.H)).astype(int)[...,None,None],counts.shape[:-1]+(1,))),axis=-1)
                            return xy,zfold,counts
                        def integrated(age,velocity=velocity.copy(),delay=delay,z0=float(z0),nd0=nd0,guide=bool(guide)):
                            a=np.asarray(age)
                            if guide:return self.density.flight(0.,4460.,np.minimum(a,delay))+nd0*np.maximum(a-delay,0.)
                            return self.density.flight(z0,velocity[2],a)
                        def proper(age,integrated=integrated,transport=transport):
                            return np.asarray(age)-transport*integrated(age)
                        if 1-transport*float(self.density.nd.max())<=0:raise ValueError('Noncausal projected acoustic energy front')
                        def survival(age,geometry=geometry,integrated=integrated,intrinsic=intrinsic,beta=beta):
                            _,_,counts=geometry(age);a=np.asarray(age);bounce=counts[...,:2].sum(axis=-1)
                            return np.exp(-intrinsic*a-beta*integrated(a))[...,None]*.8**bounce
                        sa=survival(edges);exit_by_bin=np.zeros((len(edges)-1,len(self.points)))
                        # Boundary events are real first passages, with one shared
                        # parent clock. No kappa*escaped-occupancy surrogate.
                        for axis in (0,1):
                            speed=velocity[axis]
                            if speed==0:continue
                            starts=self.points[:,axis]+s.L/2;end=starts+speed*max(0.,horizon-delay)
                            qlo=int(np.floor(min(starts.min(),end.min())/s.L));qhi=int(np.ceil(max(starts.max(),end.max())/s.L))
                            for qplane in range(qlo,qhi+1):
                                time=delay+(qplane*s.L-starts)/speed;valid=(time>delay)&(time<horizon)
                                if not np.any(valid):continue
                                ids=np.flatnonzero(valid);tt=time[valid];other=1-axis
                                raw_other=(self.points[ids,other]+s.L/2+velocity[other]*(tt-delay))/s.L
                                near=np.abs(raw_other-np.rint(raw_other))<32*np.finfo(float).eps*np.maximum(1.,abs(raw_other))
                                previous=np.abs(np.floor(raw_other)).astype(int)
                                if velocity[other]>0:previous=np.where(near,np.rint(raw_other).astype(int)-1,previous)
                                elif velocity[other]<0:previous=np.where(near,-np.rint(raw_other).astype(int),previous)
                                previous+=near*(other<axis)
                                before=qplane-1 if speed>0 else -qplane
                                amount=.2*.8**(before+previous)*np.exp(-intrinsic*tt-beta*integrated(tt))
                                ii=np.searchsorted(edges,tt,side='right')-1;np.add.at(exit_by_bin,(ii,ids),amount)
                                self.boundary.append((proper(tt),amount*self.weights[ids]*fractions[depth],energies.copy(),port,pol,di,axis))
                        absorption=sa[:-1]-sa[1:]-exit_by_bin
                        error=float(absorption.min())
                        if error < -2e-13:raise ValueError('Negative conservative acoustic loss '+str(error))
                        # Keep the signed floating-point remainder in the record;
                        # it is measured, not clipped into a positive source.
                        termination=sa[-1].copy();tail=float(np.sum(self.weights*termination))
                        means=[];unitmaps=[]
                        for j,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
                            effective=(intrinsic*(b-a)+beta*float(integrated(b)-integrated(a)))/(b-a)
                            x=effective*(b-a)
                            mean=a+(b-a)*(.5-x/12+x**3/720) if x<1e-3 else a+1/effective-(b-a)/np.expm1(min(x,700.))
                            xy,zz,_=geometry(mean);point=np.column_stack((xy[0],np.full(len(self.points),float(zz))))
                            amount=absorption[j]*self.weights*fractions[depth]
                            spatial=f.deposit(point,amount)
                            native_hazard=beta*float(integrated(b)-integrated(a));native_fraction=native_hazard/(native_hazard+intrinsic*(b-a))
                            self.add_loss(float(proper(mean)),point,amount,energies,native_fraction,port,pol,di)
                        xy,zz,_=geometry(horizon);point=np.column_stack((xy[0],np.full(len(self.points),float(zz))))
                        self.add_loss(float(proper(horizon)),point,termination*self.weights*fractions[depth],energies,
                            beta*nd0/(intrinsic+beta*nd0),port,pol,di)
                        # Snapshot spatial Green retains the causal characteristic.
                        self.rays.append((port,pol,di,energies.copy(),fractions[depth],edges,geometry,survival,integrated,proper,k,n.copy(),velocity.copy(),intrinsic,beta,correction,hr,velocity_metric[depth],delay,
                            float(guide['hydro_phase'][depth]) if guide else 0.))
                        mean_hazard=beta*float(self.density.local(z0))
                        coupled=response.coefficients(arm,frequency,nd0,vp,vg,hr,temperature,linear_strain,
                            frozen=control=='frozen_G3_reference',hydro_logarithmic_k_slope=slope)
                        self.rows.append(dict(arm=arm['name'],closure=arm['closure'],allocation=arm['allocation'],beam=beam,E_cm_keV=ec,
                            control=control,port=port,band=ch['band'],polarization=pol,direction=di,guided=ch['guided'],frequency_Hz=frequency,
                            phase_velocity_m_s=vp,host_group_velocity_m_s=vg,group_x_m_s=float(velocity[0]),group_y_m_s=float(velocity[1]),group_z_m_s=float(velocity[2]),
                            coupled_phase_velocity_at_birth_m_s=coupled['phase_velocity'] if self.coupled else vp,
                            coupled_group_velocity_at_birth_m_s=coupled['group_velocity'] if self.coupled else vg,
                            coupled_energy_velocity_at_birth_m_s=coupled['energy_velocity'] if self.coupled else vg,
                            source_power_W=float(s.current(ec)/s.QE*np.dot([.85,.15],energies)*fractions[depth]),
                            intrinsic_energy_loss_rate_s_inverse=intrinsic,native_D_loss_rate_at_birth_s_inverse=mean_hazard,
                            elastic_scattering_rate_s_inverse=0.,intrinsic_inelastic_length_m=vg/intrinsic,
                            native_bulk_dynamic_real_at_birth_Pa=float(r['bulk_per_D_J'].real*nd0),native_bulk_dynamic_imag_at_birth_Pa=float(r['bulk_per_D_J'].imag*nd0),
                            age_quadrature_tail_fraction=tail,minimum_signed_roundoff_loss=error,
                            boundary='free faces R=1; edge R=.8, external mount T+A=.2; normal-incidence strain sign minus, otherwise stochastic phase',
                            loss_semantics='decomposed intrinsic floor + same retarded D susceptibility' if self.coupled else 'inherited total-loss reference'))
        self.delays=np.array([x[0] for x in self.volume]);self.maps=np.array([x[1] for x in self.volume])
        self.weights_energy=np.array([x[2] for x in self.volume]);self.native_fraction=np.array([x[3] for x in self.volume])
        self.exit_delays=np.concatenate([x[0] for x in self.boundary]) if self.boundary else np.zeros(0)
        self.exit_weights=np.concatenate([x[1][:,None]*x[2][None,:] for x in self.boundary]) if self.boundary else np.zeros((0,2))
        self.mass=np.sum(self.maps,axis=(1,2,3))[:,None]*self.weights_energy
        fraction=np.array([.85,.15])
        self._exit=DelayMeasure(self.exit_delays,self.exit_weights@fraction)
        self._inventory=DelayMeasure(np.r_[self.delays,self.exit_delays],np.r_[self.mass@fraction,self.exit_weights@fraction])
        self._band_inventory={}
        for band in range(7):
            for pol in range(3):
                take=np.array([launch.ports()[x[4]]['band']==band and x[5]==pol for x in self.volume])
                boundary=[x for x in self.boundary if launch.ports()[x[3]]['band']==band and x[4]==pol]
                ages=np.concatenate([x[0] for x in boundary]) if boundary else np.zeros(0)
                weights=np.concatenate([x[1]*float(x[2]@fraction) for x in boundary]) if boundary else np.zeros(0)
                self._band_inventory[band,pol]=DelayMeasure(np.r_[self.delays[take],ages],np.r_[(self.mass@fraction)[take],weights])
        closure=self.mass.sum(axis=0)+self.exit_weights.sum(axis=0)-self.energy.sum(axis=(1,2,3))
        if np.max(abs(closure))>2e-12*max(self.energy.sum(),1e-300):raise ValueError('Acoustic source ports do not close')
        self.maximum_delay=max(float(self.delays.max()),float(self.exit_delays.max(initial=0.)))
        self.heat_delays=np.array([x[0] for x in self.heat]);self.heat_maps=np.array([x[1] for x in self.heat])
        self.heat_energy=np.array([x[2] for x in self.heat]);self.heat_native=np.array([x[3] for x in self.heat])
        self.maximum_delay=max(self.maximum_delay,float(self.heat_delays.max()))
        self.volume.clear();self.heat.clear()

    def add_loss(self,age,points,amount,energies,native_fraction,port,pol,di):
        self.volume.append((age,f.deposit(points,amount),energies.copy(),native_fraction,port,pol,di))
        # The same Pd/electron bath receives the debit once. Vectorized angular
        # deposition is identical to the six individual conservative deposits.
        directions=launch.DIRECTIONS;weights=np.repeat(amount/6,6)
        for channel,part,tau,path in ((6,1-.2*native_fraction,self.arm['pd'],self.arm['lambda_pd']),
                                     (8,.2*native_fraction,self.arm['te'],self.arm['lambda_e'])):
            if part==0:continue
            if channel==8:speed=np.full(6,1e6)
            else:
                q=((np.arange(6)+.5)/6)**(1/3)*np.pi/2
                speed=np.array([2440.,3300.,4460.,2440.,3300.,4460.])*np.cos(q)
            destinations=[]
            for duration in (tau,tau/2):
                flight=np.minimum(duration,path/speed)
                distance=np.minimum(speed*duration,speed*flight+np.sqrt(2*speed*path*np.maximum(0.,duration-flight)))
                location=(points[:,None,:]+directions[None,:,:]*distance[None,:,None]).reshape(-1,3)
                for axis,length,origin in ((0,s.L,-s.L/2),(1,s.L,-s.L/2),(2,s.H,0.)):
                    raw=np.mod(location[:,axis]-origin,2*length);location[:,axis]=np.where(raw<=length,raw,2*length-raw)+origin
                destinations.append(f.deposit(location,weights))
            self.heat.append((age+tau,destinations[0],energies*part,native_fraction))
            self.bath.append((age,age+tau,destinations[1],energies*part,channel))

    def bath_snapshot(self,calendar,time):
        output=np.zeros((len(f.CHANNELS),)+f.VOL.shape);rate=s.current(self.ec)/s.QE
        for a,b,spatial,energy,channel in self.bath:
            number=rate*((b-a) if calendar is None else float(calendar.integral(time-b,time-a,self.beam)))
            output[channel]+=number*float(energy@np.array([.85,.15]))*spatial
        return output

    def inventory(self,calendar,time):
        return s.current(self.ec)/s.QE*self._inventory.inventory(calendar,time,self.beam)

    def exits(self,calendar,a,b):
        return s.current(self.ec)/s.QE*(self._exit.departed(calendar,b,self.beam)-self._exit.departed(calendar,a,self.beam))

    def heat_shape(self):
        total=np.einsum('l,lijk->ijk',self.heat_energy@np.array([.85,.15]),self.heat_maps)
        optical=np.einsum('l,lijk->ijk',(self.heat_energy@np.array([.85,.15]))*self.heat_native,self.heat_maps)
        return total,optical

    def snapshot(self,calendar,time):
        output=np.zeros((7,3)+f.VOL.shape);phase=np.zeros((8,3,6)+f.VOL.shape,complex)
        for port,pol,di,energies,fraction,edges,geometry,survival,integrated,proper,k,n,velocity,intrinsic,beta,correction,hr,vmetric,delay,hphase in self.rays:
            frequency=launch.ports()[port]['frequency'];power=s.current(self.ec)/s.QE*np.dot([.85,.15],energies)*fraction
            for a,b in zip(edges[:-1],edges[1:]):
                age=(a+b)/2;xy,z,_=geometry(age);survive=survival(age)[0]
                points=np.column_stack((xy[0],np.full(len(self.points),float(z))))
                duration=float(proper(b)-proper(a));clock=float(calendar.value(time-float(proper(age)),self.beam))
                weight=power*duration*clock*self.weights*survive
                output[launch.ports()[port]['band'],pol]+=f.deposit(points,weight)
                angle=-2*np.pi*frequency*float(proper(age))+k*float(n@velocity)*max(0.,age-delay)-.5*2*np.pi*frequency*correction*float(integrated(age))+hphase
                # Only normal-incidence free-return signs are determined. An
                # oblique unresolved boundary removes the first moment, not E.
                count=geometry(age)[2][0];known=int(np.argmax(abs(n)))
                coherent=np.sum(np.delete(count,known,axis=1),axis=1)==0;sign=(-1.)**count[:,known]
                phase[port,pol,di]+=f.deposit(points,duration*clock*self.weights*np.sqrt(survive)*sign*coherent*np.exp(1j*angle))
        normalizers=np.ones((7,3))
        for band in range(7):
            for pol in range(3):
                exact=s.current(self.ec)/s.QE*self._band_inventory[band,pol].inventory(calendar,time,self.beam)
                raw=float(output[band,pol].sum())
                if raw<=0 and exact>0:raise ValueError('Retarded quadrature misses populated acoustic band')
                if raw>0:normalizers[band,pol]=exact/raw;output[band,pol]*=exact/raw
        return output,phase,normalizers

    def delayed_heat(self,counts,dt,start=0):
        """Actual common species counts, with conditional spatial beam moments."""
        steps=counts.shape[1];lag=self.heat_delays/dt;integer=lag.astype(int);fraction=lag-integer
        length=int(integer.max())+2;result=np.zeros((steps-start,)+f.VOL.shape)
        for sp in (0,1):
            for lo in range(0,f.VOL.size,48):
                hi=min(f.VOL.size,lo+48);kernel=np.zeros((length,hi-lo))
                weighted=self.heat_maps.reshape(len(self.heat_maps),-1)[:,lo:hi]*self.heat_energy[:,sp,None]
                np.add.at(kernel,integer,weighted*(1-fraction[:,None]));np.add.at(kernel,integer+1,weighted*fraction[:,None])
                conv=fftconvolve(counts[sp,:,None],kernel,mode='full',axes=0)[start:steps]
                result.reshape(steps-start,-1)[:,lo:hi]+=conv
        return result
# Operational kernel reuse; the physics above is independently fingerprinted.
_family=None
_operators={}


def cached_operator(arm,ec,beam,tracks,control='joined',temperature=300.,linear_strain=0.):
    """Reuse only identical geometry, sources, density and linear material state."""
    global _family,_operators
    shape='even' if beam=='even_steady' else 'hot'
    family=(arm['closure'],ec,shape,temperature,linear_strain,arm.get('mesh',16))
    if family!=_family:
        _operators={};_family=family
        __import__('gc').collect()
    key=(control not in ('strain_off','frozen','inherited_total'),control=='frozen_G3_reference')
    if key not in _operators:_operators[key]=Operator(arm,ec,beam,tracks,control,temperature,linear_strain)
    result=_operators[key]
    # Onsite/coupling have the same U=0 native branch here. Hot steady/Cox
    # share spatial weights; only their calendar application differs.
    result.arm=arm;result.beam=beam;result.control=control
    result.rows=[dict(row,arm=arm['name'],allocation=arm['allocation'],beam=beam,control=control) for row in result.rows]
    return result
