"""Actual-coordinate conditional exterior field and globally named local parents.

This is an explicitly conditional far closure, NOT a measured atomic phase chart.
The finite spectral Gaussian exterior is assigned only its Green variance. Local
rare parents have physical times/counts, with one mark reused across all bands.
Source-interior elasticity and higher-return covariance remain separate omissions.
"""
from pathlib import Path
import math
import numpy as np
from numpy.polynomial.legendre import legval
from common import QE,RHO,L,H,ENERGIES,BEAMS,inherited,calendar
from receiver import Green


def _hash(x):
    x=np.asarray(x,dtype=np.uint64)
    with np.errstate(over='ignore'):
        x=(x ^ (x >> np.uint64(30)))*np.uint64(0xbf58476d1ce4e5b9)
        x=(x ^ (x >> np.uint64(27)))*np.uint64(0x94d049bb133111eb)
    return x ^ (x >> np.uint64(31))


def uniform(x,stream):
    with np.errstate(over='ignore'):
        value=_hash(np.asarray(x,dtype=np.uint64)+np.uint64(stream)*np.uint64(0x9e3779b97f4a7c15))
    return ((value >> np.uint64(11)).astype(float)+.5)/2**53


class Strain:
    def __init__(self,bundle,path4=None):
        self.bundle=Path(bundle);self.path4=Path(path4) if path4 else self.bundle/'inherited'
        self.s,self.native,self.fields,self.launch,self.response,self.acoustic,self.yielding=inherited(self.path4)
        self.calendar=calendar(self.path4);self.greens={}
        from low import Low
        self.low=Low(self.bundle,self.path4)
        self.near_radius=1e-6;self.cell=1e-6;self.bucket=1e-9
        rng=np.random.default_rng(517112)
        self.normal=rng.normal(size=(8,32,2));self.azimuth=rng.uniform(0,2*np.pi,size=(8,32))
        import csv
        shape=list(csv.DictReader((self.bundle/'spectral_shape.csv').open()))
        self.frequencies=[]
        for ch in self.launch.ports():
            rows=[r for r in shape if r['scenario']=='track' and int(r['species'])==0 and int(r['band'])==ch['band']]
            if rows and not ch['guided']:
                frequency=np.array([float(r['frequency_Hz']) for r in rows]);weight=np.array([float(r['energy_J']) for r in rows]);weight/=weight.sum()
                index=np.searchsorted(np.cumsum(weight),(np.arange(32)+.5)/32)
                self.frequencies.append(frequency[index])
            else:
                lo=2.5e5*10**ch['band'];hi=min(10*lo,1e7)
                self.frequencies.append(np.geomspace(lo,hi,32))
        from numpy.polynomial.legendre import leggauss
        q,w=leggauss(2);edge=np.linspace(np.log(2e-9),np.log(2*L),33)
        logarithm=((edge[:-1,None]+edge[1:,None])/2+(edge[1:,None]-edge[:-1,None])*q/2).ravel()
        radial_weight=((edge[1:,None]-edge[:-1,None])*w/2*np.exp(2*((edge[:-1,None]+edge[1:,None])/2+(edge[1:,None]-edge[:-1,None])*q/2))).ravel()
        azimuth=2*np.pi*(np.arange(8)+.5)/8
        self.source_offset=np.stack(np.broadcast_arrays(np.exp(logarithm[:,None])*np.cos(azimuth),np.exp(logarithm[:,None])*np.sin(azimuth)),axis=-1).reshape(-1,2)
        self.source_area=np.repeat(radial_weight*2*np.pi/8,8)
    def green(self,ec,allocation,closure='nominal_retained'):
        key=(ec,allocation,closure)
        if key not in self.greens:self.greens[key]=Green(self.bundle,self.path4,ec,closure=closure,allocation=allocation)
        return self.greens[key]
    def variance(self,g,beam,xyz,time_s):
        """Whole-beam direct Green variance, actual xyz and reflected D(z) loading.

        Exterior first-return interference is not represented by this Gaussian
        variance and must be priced with the independent total-parent quadrature.
        """
        xyz=np.asarray(xyz);xy=xyz[:,None,:2]+self.source_offset[None,:,:]
        weight=self.s.beam(xy[...,0],xy[...,1],beam)*self.source_area[None,:]/self.s.beam_moments(beam)[0]
        result=np.zeros((len(xyz),8))
        for sp,fraction in ((0,.85),(1,.15)):
            source=np.concatenate((xy,np.full(xy.shape[:-1]+(1,),g.source_depth[sp])),axis=-1)
            delta=xyz[:,None,:]-source
            distance=np.linalg.norm(delta,axis=-1);parallel=np.linalg.norm(delta[...,:2],axis=-1)
            for pi,ch in enumerate(g.channels):
                rr=np.maximum(parallel if ch['guided'] else distance,2e-9)
                direction=delta/np.maximum(distance,2e-9)[...,None]
                if ch['guided']:direction=np.concatenate((delta[...,:2]/np.maximum(parallel,2e-9)[...,None],np.zeros(rr.shape+(1,))),axis=-1)
                aa,var=g._angular(pi,sp,direction);energy=abs(aa)**2+var
                for pol in range(3):
                    vp,vg,intrinsic,response,hr=g.material[pi,pol]
                    if ch['guided']:
                        z,div,displacement,k,vp,vg,div_coeff=g.guides[pi,pol]
                        profile=legval(2*xyz[:,2]/H-1,div_coeff)
                        metric=abs(profile[:,None])**2/(2*np.pi*rr*vg)
                    else:
                        # Longitudinal projection of the retained Christoffel
                        # polarization; no transverse velocity is counted as z.
                        nearest=np.argmax(direction@g.data[f'n_{sp}_{pi}'].T,axis=-1)
                        polvec=g.data[f'polar_{sp}_{pi}'][pol,nearest]
                        projection=np.sum(direction*polvec,axis=-1)**2
                        metric=projection/(4*np.pi*rr*rr*vg*RHO*vp*vp)
                    age=rr/vg
                    if ch['guided']:
                        integral=g.guide_effective_density[pi,pol]*age
                    else:
                        vz=np.divide(delta[...,2],age,out=np.zeros_like(age),where=age>0)
                        integral=np.divide(g.density.primitive(g.source_depth[sp]+vz*age)[0]-g.density.primitive(g.source_depth[sp])[0],vz,
                            out=g.density.local(g.source_depth[sp])*age,where=abs(vz)>1e-12)
                    survive=np.exp(-np.minimum(intrinsic*age+g.loading[pi,pol][2]*integral,745.))
                    keep=(rr>=self.near_radius) if pi>=4 else 1.
                    arrival=age+g.loading[pi,pol][1]*integral+(xyz[:,2,None]/4460. if ch['guided'] else 0.)
                    clock=self.calendar.value(np.asarray(time_s)[:,None]-arrival,beam)
                    result[:,pi]+=fraction*self.s.current(g.ec)/QE*np.sum(energy[...,pol]*metric*survive*keep*weight*clock,axis=1)
        return result
    def _birth_bins(self,g,xyz,t,centers,band):
        """Union retarded finite-support windows BEFORE drawing parent counts."""
        from receiver import packet_duration
        rows=[set() for _ in t];spread=np.sqrt(2)*self.cell/2
        lo_port=4 if band is None else band;hi_port=7 if band is None else band
        for sp in (0,1):
            src=np.c_[centers,np.full(len(t),g.source_depth[sp])];images=[src]
            for axis in (0,1):
                for wall in (-L/2,L/2):
                    image=src.copy();image[:,axis]=2*wall-src[:,axis];images.append(image)
            for wall in (0.,H):
                image=src.copy();image[:,2]=2*wall-src[:,2];images.append(image)
            for index,(pi,omega,a,arrival,c) in enumerate(g.components(xyz,centers,sp,min_port=lo_port,max_port=hi_port)):
                image=images[(index//3)%7];pol=index%3
                distance=np.linalg.norm(xyz-image,axis=1);vg=g.material[pi,pol][1]
                limit=(1+g.loading[pi,pol][1]*g.density.local(image[:,2]))/vg
                inverse_speed=np.divide(arrival,distance,out=limit.copy(),where=distance>0)
                lower=np.floor((t-arrival-inverse_speed*spread-packet_duration(a,omega))/self.bucket).astype(np.int64)
                upper=np.floor((t-np.maximum(0.,arrival-inverse_speed*spread))/self.bucket).astype(np.int64)
                for j in range(len(t)):rows[j].update(range(max(0,int(lower[j])),int(upper[j])+1))
        result=np.full((len(t),max(map(len,rows),default=0)),-1,dtype=np.int64)
        for j,row in enumerate(rows):result[j,:len(row)]=sorted(row)
        return result

    def _near(self,g,beam,xyz,t,independent=False):
        eligible=(np.minimum(abs(xyz[:,2]-g.source_depth[0]),abs(xyz[:,2]-g.source_depth[1]))<=self.near_radius)
        eligible&=(np.max(abs(xyz[:,:2]),axis=1)<=.005+self.near_radius)
        if not np.all(eligible):
            result=np.zeros((len(t),8),complex)
            if np.any(eligible):result[eligible]=self._near(g,beam,xyz[eligible],t[eligible],independent)
            return result
        n=len(t);radius=self.near_radius;cell=self.cell;bucket=self.bucket
        # Enumerate each physical parent once across all band/return support.
        ix0=np.floor(xyz[:,0]/cell).astype(np.int64);iy0=np.floor(xyz[:,1]/cell).astype(np.int64)
        time_bin,time_remainder=np.divmod(t,bucket);time_bin=time_bin.astype(np.int64)
        ecindex=ENERGIES.index(float(g.ec));beamindex=BEAMS.index(beam)
        output=np.zeros((n,8),complex);upper=12. if beam=='hot_bursts' else 1.
        bands=range(4,8) if independent else (None,)
        for band in bands:
            for dx in (-1,0,1):
                for dy in (-1,0,1):
                    ix=ix0+dx;iy=iy0+dy
                    centers=np.c_[(ix+.5)*cell,(iy+.5)*cell]
                    # Gaussian beam maximum is bounded over each tiny voxel by
                    # its corner/center envelope; the factor 1.01 covers slope.
                    local=1.01*self.s.beam(centers[:,0],centers[:,1],beam)/self.s.beam_moments(beam)[0]
                    lam=self.s.current(g.ec)/QE*local*cell*cell*bucket*upper
                    bins=self._birth_bins(g,xyz,t,centers,band)
                    for slot in range(bins.shape[1]):
                        it=bins[:,slot]
                        with np.errstate(over='ignore'):
                            key=(np.asarray(ix+100000,dtype=np.uint64)*np.uint64(0x9e3779b185ebca87)) ^ (np.asarray(iy+100000,dtype=np.uint64)*np.uint64(0xc2b2ae3d27d4eb4f)) ^ np.asarray(it,dtype=np.uint64)
                            key^=np.uint64((ecindex+1)*31337+(beamindex+1)*11213+(0 if band is None else band+1)*773)
                        u=uniform(key,0);prob=np.exp(-lam);cdf=prob.copy()
                        count=np.zeros(n,dtype=int);active=(u>cdf)&(it>=0)
                        # Exact inverse Poisson CDF for these small physical
                        # voxel intensities. No aggregated superparent marks.
                        k=0
                        while np.any(active):
                            k+=1;count[active]=k;prob*=lam/k;cdf+=prob;active=(u>cdf)&(it>=0)
                        for event in range(int(count.max(initial=0))):
                            take=count>event
                            ids=np.flatnonzero(take)
                            if not len(ids):continue
                            kk=key[take];stream=10+event*10
                            xy=np.c_[(ix[take]+uniform(kk,stream))*cell,(iy[take]+uniform(kk,stream+1))*cell]
                            local_birth=uniform(kk,stream+2)*bucket;birth=it[take]*bucket+local_birth
                            physical=self.s.beam(xy[:,0],xy[:,1],beam)/self.s.beam_moments(beam)[0]
                            accept=(birth<=t[take])&(birth>=0)&(np.linalg.norm(xy-xyz[take,:2],axis=1)<=radius)
                            accept&=uniform(kk,stream+3)<self.calendar.value(birth,beam)*np.divide(physical,local[take]*upper,out=np.zeros_like(physical),where=local[take]>0)
                            if not np.any(accept):continue
                            ids=ids[accept];xy=xy[accept];kk=kk[accept]
                            age=(time_bin[ids]-it[ids])*bucket+time_remainder[ids]-local_birth[accept]
                            species=(uniform(kk,stream+4)>=.85).astype(int);recoil=np.floor(uniform(kk,stream+5)*64).astype(int)
                            for sp in (0,1):
                                select=(species==sp)&(np.sum((xy-xyz[ids,:2])**2,axis=1)+(xyz[ids,2]-g.source_depth[sp])**2<=radius*radius)
                                if not np.any(select):continue
                                h,_,per=g.parent(xyz[ids[select]],xy[select],age[select],sp,recoil[select],
                                    min_port=4 if band is None else band,max_port=7 if band is None else band)
                                np.add.at(output,ids[select],per)
        return output
    def bands(self,ec,beam,xyz_m,mount='nominal',time_s=2,allocation='onsite',control='joined',closure='nominal_retained'):
        xyz=np.asarray(xyz_m,float);shape=np.broadcast_shapes(xyz.shape[:-1],np.shape(time_s));xyz=np.broadcast_to(xyz,shape+(3,));tt=np.broadcast_to(time_s,shape)
        if control=='strain_off':return np.zeros(shape+(8,),complex)
        if np.any(~np.isfinite(tt)) or np.any((tt<0)|(tt>2)):raise ValueError('Signed realization requires finite physical time in [0,2] s')
        g=self.green(float(ec),allocation,closure);points=xyz.reshape(-1,3);times=tt.ravel();values=np.empty((len(times),8),complex)
        for start in range(0,len(times),256):
            end=min(start+256,len(times));p=points[start:end];t=times[start:end]
            variance=self.variance(g,beam,p,t);value=np.zeros((len(t),8),complex)
            for pi,ch in enumerate(g.channels):
                vp=g.material[pi,2][0];frequency=self.frequencies[pi]
                clock=np.remainder(t[:,None].astype(np.longdouble)*frequency[None,:].astype(np.longdouble),1.)
                path=frequency[None,:]*(np.cos(self.azimuth[pi])[None,:]*p[:,0,None]+np.sin(self.azimuth[pi])[None,:]*p[:,1,None])/vp
                phase=2*np.pi*(np.asarray(clock,float)-path)
                coefficient=self.normal[pi]
                value[:,pi]=np.sqrt(variance[:,pi]/32)*np.sum((coefficient[None,:,0]+1j*coefficient[None,:,1])*np.exp(-1j*phase),axis=1)
            value+=self._near(g,beam,p,t,independent=control=='parent_resampled')
            value+=g.mean(p,beam,t,self.calendar,per_band=True)
            values[start:end]=value
        return values.reshape(shape+(8,))

    def sample(self,ec,beam,xyz_m,mount='nominal',time_s=2,allocation='onsite',control='joined',closure='nominal_retained'):
        field=self.bands(ec,beam,xyz_m,mount,time_s,allocation,control,closure)
        low=0. if control=='strain_off' else self.low.mean(float(ec),beam,xyz_m,time_s,closure,allocation)
        return (field.sum(axis=-1).real+low)[...,None,None]*np.eye(3)/3

    def covariance(self,ec,beam,xyz_m,temperature_K,time_s=2,allocation='onsite',closure='nominal_retained',band_field=None):
        field=self.bands(ec,beam,xyz_m,time_s=time_s,allocation=allocation,closure=closure) if band_field is None else np.asarray(band_field)
        temperature=np.broadcast_to(temperature_K,field.shape[:-1])
        result=np.zeros(field.shape[:-1]);g=self.green(float(ec),allocation,closure)
        flat=field.reshape(-1,8)
        for j,T in enumerate(temperature.ravel()):
            coefficients=np.array([self.response.covariance(g.arm,ch['frequency'],float(T))['theta_per_linear_strain_eV'] for ch in g.channels])
            result.ravel()[j]=np.real(coefficients@flat[j])
        return result
