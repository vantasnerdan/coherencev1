"""One source-to-foil carrier Green. Physical birth depth is retained below the
FV grid; a displacement is never applied to a fictitious uniform depth stock.
The lateral source is the declared finite-volume beam-cell distribution.
"""
import math
import numpy as np
from scipy import sparse
import sources as s
import fields as f


def birth_distribution(track,sp,origin):
    z=np.r_[0.,[r['depth_A']*1e-10 for r in track['timeline']]]
    t=np.r_[0.,[r['time_fs']*1e-15 for r in track['timeline']]]
    delay=.74e-10/math.sqrt(2*track['initial_eV']*s.QE/s.MD)/2 if sp else 0.
    if origin==0:return np.array([z[-1]]),np.array([1.]),float(t[-1]+delay)
    key={1:'host_recoil_eV',2:'D_recoil_eV',3:'electron_eV'}[origin]
    energy=np.diff(np.r_[0.,[r[key] for r in track['timeline']]])
    meanz=(z[:-1]+z[1:])/2;half=np.diff(z)/(2*math.sqrt(3))
    positions=np.r_[meanz-half,meanz+half];weights=np.r_[energy,energy]/(2*energy.sum())
    return positions,weights,float(((t[:-1]+t[1:])/2)@energy/energy.sum()+delay)


def depth_probability(birth,offset):
    positions,weights,age=birth;raw=np.mod(positions+offset,2*s.H)
    folded=np.where(raw<=s.H,raw,2*s.H-raw)
    index=np.minimum(np.searchsorted(f.ZE[1:],folded,side='right'),2)
    return np.bincount(index,weights=weights,minlength=3)


def axis_transfer(indices,offset):
    left=np.maximum(-.005,-s.L/2+indices*f.DX)+s.L/2+offset
    right=np.minimum(.005,-s.L/2+(indices+1)*f.DX)+s.L/2+offset
    width=right-left
    if np.any(width<=0):raise ValueError('Source cell outside physical beam')
    first=np.floor(left/f.DX).astype(int);unfolded=np.stack((first,first+1))
    amount=np.maximum(0.,np.minimum(right[None,:],(unfolded+1)*f.DX)-np.maximum(left[None,:],unfolded*f.DX))/width
    fold=np.mod(unfolded,2*f.N);destination=np.where(fold<f.N,fold,2*f.N-1-fold)
    return destination,amount


class Operator:
    def __init__(self,arm,ec,kernels,tracks,ratios):
        import run as model
        self.cells=np.flatnonzero(model.spatial_source('even_steady')[:,:,0].ravel());ns=len(self.cells);ng=f.VOL.size
        ix=self.cells//f.N;iy=self.cells%f.N;cols=np.arange(ns)
        self.fast=np.zeros((2,len(f.CHANNELS),ng,ns));self.entries=[];self.heat=[{},{}];self.delay_s=0.
        endpoints=np.load(s.HERE/f"paths_{arm['name']}.npz")['heat_endpoints']
        for sp in (0,1):
            track=tracks[f'{ec:g}_{sp}'];totals=model.source_totals(track,sp,arm)
            birth=[birth_distribution(track,sp,origin) for origin in range(4)]
            acoustic=float(np.load(s.HERE/f"source_E{ec:g}_{arm['closure']}_sp{sp}.npz")['energy_J'].sum())
            heat_energy=2000*ec*s.QE-acoustic
            if heat_energy<=0:raise ValueError('Acoustic source exceeds actual parent energy')
            for (channel,origin),(offset,age,weights,residence) in kernels.items():
                for delta,delay,weight in zip(offset,age,weights):
                    x,wx=axis_transfer(ix,delta[0]);y,wy=axis_transfer(iy,delta[1]);pz=depth_probability(birth[origin],delta[2])
                    rows=[];values=[]
                    for a in (0,1):
                        for b in (0,1):
                            for z in range(3):
                                rows.append((x[a]*f.N+y[b])*3+z);values.append(wx[a]*wy[b]*pz[z])
                    rows=np.array(rows);values=np.array(values)
                    coefficient=totals[origin]*s.QE*residence*weight*ratios[sp,channel]
                    for row,value in zip(rows,values):self.fast[sp,channel,row,cols]+=coefficient*value
                    self.entries.append((sp,channel,float(delay+birth[origin][2]),coefficient,rows,values))
            # In-flight primary kinetic energy precedes every deposited carrier.
            # It travels on the saved track; it is not added to the end kinetic.
            times=np.r_[0.,[r['time_fs']*1e-15 for r in track['timeline']]]
            energy=np.r_[track['initial_eV'],[r['energy_eV'] for r in track['timeline']]]
            for age,dt,e in zip((times[:-1]+times[1:])/2,np.diff(times),(energy[:-1]+energy[1:])/2):
                coefficient=(sp+1)*e*dt*s.QE;rows=np.array([self.cells*3]);values=np.ones((1,ns))
                self.fast[sp,9,rows[0],cols]+=coefficient
                self.entries.append((sp,9,float(age+birth[0][2]-times[-1]),coefficient,rows,values))
            for row in endpoints:
                origin=int(row[0]);delay=float(row[2]+birth[origin][2]);self.delay_s=max(self.delay_s,delay)
                self.heat[sp].setdefault(delay,[]).append((row[3:6].copy(),depth_probability(birth[origin],row[5]),
                    heat_energy*totals[origin]/totals[:4].sum()*row[6]))
        self.ns=ns;self.ng=ng;self.ix=ix;self.iy=iy
        self.steady_heat=np.zeros((2,ng,ns))
        for sp in (0,1):
            for delay,records in self.heat[sp].items():
                for delta,pz,coefficient in records:
                    x,wx=axis_transfer(ix,delta[0]);y,wy=axis_transfer(iy,delta[1])
                    for a in (0,1):
                        for b in (0,1):
                            for z in range(3):self.steady_heat[sp,(x[a]*f.N+y[b])*3+z,cols]+=coefficient*wx[a]*wy[b]*pz[z]

    def field(self,source,rate,calendar=None,at=None,beam=None):
        probability=source[:,:,0].ravel()[self.cells]
        if calendar is None:
            result=np.einsum('s,scri,i->cr',np.array([.85,.15]),self.fast,rate*probability)
        else:
            result=np.zeros((len(f.CHANNELS),self.ng))
            for sp,channel,age,coefficient,rows,values in self.entries:
                scale=(.85 if sp==0 else .15)*rate*float(calendar.value(at-age,beam))*coefficient
                for row,value in zip(rows,values):np.add.at(result[channel],row,scale*probability*value)
        return result.reshape((len(f.CHANNELS),)+f.VOL.shape)

    def heat_rate(self,source,rate):
        probability=source[:,:,0].ravel()[self.cells]
        return np.einsum('s,sri,i->r',np.array([.85,.15]),self.steady_heat,rate*probability).reshape(f.VOL.shape)

    def delayed_heat(self,counts,dt):
        """Conservative fractional-delay handoff of the very same beam counts.
        counts[species,time,source-cell]; source-parent phases are not redrawn.
        """
        steps=counts.shape[1];out=np.zeros((steps,self.ng));cols=np.arange(self.ns);operators=[{},{}]
        for sp in (0,1):
            for delay,records in self.heat[sp].items():
                position=delay/dt;integer=int(position);fraction=position-integer
                for delta,pz,coefficient in records:
                    x,wx=axis_transfer(self.ix,delta[0]);y,wy=axis_transfer(self.iy,delta[1])
                    for lag,part in ((integer,1-fraction),(integer+1,fraction)):
                        if part==0:continue
                        rr,cc,vv=operators[sp].setdefault(lag,([],[],[]))
                        for a in (0,1):
                            for b in (0,1):
                                for z in range(3):
                                    rr.extend(((x[a]*f.N+y[b])*3+z).tolist());cc.extend(cols.tolist());vv.extend((coefficient*part*wx[a]*wy[b]*pz[z]).tolist())
            for lag,(rr,cc,vv) in operators[sp].items():
                matrix=sparse.coo_matrix((vv,(rr,cc)),shape=(self.ng,self.ns)).tocsr()
                if lag<steps:out[lag:]+=(matrix@counts[sp,:steps-lag].T).T
        return out.reshape((steps,)+f.VOL.shape)
