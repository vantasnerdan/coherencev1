"""Accepted 303-mode prompt/pressure receiver without replaying its producer.

Mean forcing follows the original complete Cox history. Prompt/fast-pressure
innovations are carried as a correlation-independent Campbell covariance bound,
not a fabricated microscopic phase record. Full inherited thermoelastic windows
are retained as finite controls; their out-of-window continuation is unresolved.
"""
from functools import lru_cache
from pathlib import Path
import csv,json
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from common import QE,MD,L,ENERGIES,BEAMS,verify_selected,calendar,digest,emit,table,inherited


def build(out,path4):
    out=Path(out);path4=Path(path4)
    names=['frozen/modes.npz','frozen/source_modes.npz']
    names.extend(f'phase_{beam}_{closure}_{allocation}.npz' for beam in BEAMS for closure in ('fast_Cu','nominal_retained','long_retained') for allocation in ('onsite','coupling'))
    binding=verify_selected(path4,names)
    data=np.load(path4/'frozen/modes.npz');source=np.load(path4/'frozen/source_modes.npz')
    points=data['points'].copy();points[:,:2]-=L/2
    source_points=source['points'].copy();source_points[:,:2]-=L/2
    payload=dict(omega=data['omega'],speed=data['speed'],pressure=data['pressure'],points_m=points,
        trace_strain=data['strain'][:,:3].sum(axis=1),source_points_m=source_points,
        source_trace_strain=source['strain'][:,:3].sum(axis=1),source_xy_m=data['source_xy']-L/2)
    rows=[]
    for beam in BEAMS:
        for closure in ('fast_Cu','nominal_retained','long_retained'):
            for allocation in ('onsite','coupling'):
                name=f'phase_{beam}_{closure}_{allocation}.npz';a=np.load(path4/name)
                key=beam+'_'+closure+'_'+allocation
                q=a['modal_phase_space'][:,0,:len(data['omega'])]-a['static_reference_q']
                payload['reference_q_'+key]=q;payload['reference_time_'+key]=a['time_s']
                signed=q@payload['trace_strain'].T
                rows.append(dict(beam=beam,closure=closure,allocation=allocation,
                    physical_window_start_s=float(a['physical_source_window_start_s']),physical_exposure_s=float(a['exposure_s']),
                    sample_time_min_s=float(a['time_s'].min()),sample_time_max_s=float(a['time_s'].max()),max_abs_dynamic_volumetric_strain=float(abs(signed).max()),
                    FV_receiver_max_time_rms=float(np.sqrt(np.mean(signed*signed,axis=0)).max()),
                    source_parent_hash=str(a['source_parent_hash']),scope='exact saved finite full-force low-mode window; not a seconds-long tail/rate bound'))
    np.savez_compressed(out/'low_modes.npz',**payload)
    table(out/'low_accepted_window.csv',rows)
    emit(out/'low_inputs.json',dict(base='f6010129f11ba8010815b78a67b3b6bed1db7194',selected_source_files=binding,
        cache_sha256=digest(out/'low_modes.npz'),projection='accepted mass-normalized FV rows plus exact source-surface rows; linear declared interpolation, with separate 0-to-125nm surface layer',
        mean_scope='prompt stopping and fast native-pressure/optical-stress impulses, exact original Cox convolution; no forced zero mean',
        unresolved='thermoelastic/slow-force innovations and continuation outside exact saved low windows; no physical exclusion inferred'))
    return rows


@lru_cache(maxsize=2)
def _basis(bundle):
    a=np.load(Path(bundle)/'low_modes.npz');data={k:a[k] for k in ('omega','speed','pressure','points_m','trace_strain','source_points_m','source_trace_strain','source_xy_m')}
    points=data['points_m'];axes=[np.unique(points[:,j]) for j in range(3)]
    data['axes']=axes;data['interpolate']=RegularGridInterpolator(axes,data['trace_strain'].reshape(*(len(x) for x in axes),-1),bounds_error=False,fill_value=None)
    surface=data['source_points_m'];sx=np.unique(surface[:,0]);sy=np.unique(surface[:,1])
    data['surface']=RegularGridInterpolator((sx,sy),data['source_trace_strain'].reshape(len(sx),len(sy),-1),bounds_error=False,fill_value=None)
    xy=data['source_xy_m'];axis=np.unique(xy[:,0]);dx=L/17;cell=((axis+L/2)/dx).astype(int)
    edges=np.r_[-.005,(axis[:-1]+axis[1:])/2,.005]
    for j in range(len(axis)-1):
        if cell[j]!=cell[j+1]:edges[j+1]=-L/2+(cell[j]+1)*dx
    data['tile_edges']=edges;data['tile_ix']=np.searchsorted(axis,xy[:,0]);data['tile_iy']=np.searchsorted(axis,xy[:,1])
    return data


class Low:
    def __init__(self,bundle,path4):
        self.bundle=str(Path(bundle).resolve());self.data=_basis(self.bundle);self.path4=Path(path4)
        self.s=inherited(path4)[0];self.clock=calendar(path4)
        self.rows=json.loads((self.path4/'native_summary.json').read_text())['rows'];self.cache={}
    def projection(self,xyz):
        xyz=np.asarray(xyz,float);shape=xyz.shape[:-1];point=xyz.reshape(-1,3);result=self.data['interpolate'](point)
        first=float(self.data['axes'][2][0]);mask=point[:,2]<first
        if np.any(mask):
            surface=self.data['surface'](point[mask,:2]);at=point[mask].copy();at[:,2]=first
            near=self.data['interpolate'](at)
            result[mask]=surface+(near-surface)*(point[mask,2,None]/first)
        return result.reshape(shape+(len(self.data['omega']),))
    def case(self,ec,beam,closure,allocation):
        key=(ec,beam,closure,allocation)
        if key in self.cache:return self.cache[key]
        omega=self.data['omega'];speed=self.data['speed']
        gamma=(4460*np.log(10)/10*32*omega/(2*np.pi*1e6) if closure=='fast_Cu' else 1/50e-6+speed/.2)
        wd=np.sqrt(omega*omega-gamma*gamma/4);poles=-gamma/2+1j*wd
        probability=self.s.beam_cells(self.data['tile_edges'],beam)[self.data['tile_ix'],self.data['tile_iy']]
        J=np.sqrt(2*MD*2000*ec*QE);matrices=[]
        for sp in (0,1):
            row=next(r for r in self.rows if r['arm']==closure+'_'+allocation and r['E_cm_keV']==ec and r['species']==sp+1 and r['mesh']==16)
            pressure=row['native_pressure_residence_J_s']-row['optical_strain_energy_ratio']*row['carrier_residence_eV_s'][4]*QE
            matrices.append(J*(np.sqrt(2) if sp else 1)*self.data['pressure']+pressure*self.data['source_trace_strain'].T)
        rate=self.s.current(ec)/QE;mean=rate*(.85*matrices[0]+.15*matrices[1])@probability
        D=rate*sum(fr*(a*probability[None,:])@a.T for fr,a in zip((.85,.15),matrices))
        decay=(gamma[:,None]+gamma[None,:])/2
        K=2*decay/((decay*decay+(wd[:,None]-wd[None,:])**2)*(decay*decay+(wd[:,None]+wd[None,:])**2))
        covariance=D*K
        result=dict(omega=omega,gamma=gamma,wd=wd,poles=poles,mean_force=mean,covariance=covariance,impulses=matrices,
            source_probabilities=probability,source_rate_s=rate)
        self.cache[key]=result;return result
    @lru_cache(maxsize=16384)
    def filtered_clock(self,time,loss):
        omega=self.data['omega'];speed=self.data['speed']
        gamma=(4460*np.log(10)/10*32*omega/(2*np.pi*1e6) if loss=='fast_Cu' else 1/50e-6+speed/.2)
        pole=-gamma/2+1j*np.sqrt(omega*omega-gamma*gamma/4);decay=np.exp(pole*time)
        constant=np.expm1(pole*time)/pole
        w=self.clock.omega;phase=self.clock.phase
        ripple=.35*(np.exp(1j*phase)*(np.exp(1j*w*time)-decay)/(1j*w-pole)+np.exp(-1j*phase)*(np.exp(-1j*w*time)-decay)/(-1j*w-pole))
        burst=np.zeros(len(pole),complex)
        for start in range(0,len(pole),12):
            p=pole[start:start+12];horizon=40/float(np.min(-p.real))
            first=max(0,int(np.searchsorted(self.clock.edges,max(0.,time-horizon),side='right')-1));last=min(len(self.clock.state),int(np.searchsorted(self.clock.edges,time,side='right')))
            a=self.clock.edges[first:last];b=np.minimum(self.clock.edges[first+1:last+1],time);mask=self.clock.state[first:last]&(a<time)
            a=a[mask];b=b[mask]
            burst[start:start+12]=np.sum(np.exp(p[:,None]*(time-b))*np.expm1(p[:,None]*(b-a))/p[:,None],axis=1)
        return (constant+ripple+9*burst)/self.clock.norm
    def mean(self,ec,beam,xyz,time,closure='nominal_retained',allocation='onsite',dynamic=False):
        xyz=np.asarray(xyz,float);time=np.asarray(time,float);shape=np.broadcast_shapes(xyz.shape[:-1],time.shape)
        point=np.broadcast_to(xyz,shape+(3,)).reshape(-1,3);times=np.broadcast_to(time,shape).ravel();S=self.projection(point)
        case=self.case(float(ec),beam,closure,allocation);q=[]
        for t in times:
            if beam=='hot_bursts':z=self.filtered_clock(float(t),'fast_Cu' if closure=='fast_Cu' else 'retained')
            else:z=np.expm1(case['poles']*t)/case['poles']
            state=case['mean_force']*z.imag/case['wd']
            if dynamic:state=state-case['mean_force']/(case['omega']**2)
            q.append(state)
        return np.einsum('ij,ij->i',S,np.array(q)).reshape(shape)
    def innovation_rms_bound(self,ec,beam,xyz,closure='nominal_retained',allocation='onsite'):
        S=self.projection(xyz);case=self.case(float(ec),beam,closure,allocation)
        variance=np.einsum('...i,ij,...j->...',S,case['covariance'],S)
        if np.any(variance<0):raise ValueError('Low-mode covariance lost positivity')
        # Conservative intensity ceiling; it is not an empirical phase peak.
        return np.sqrt(variance*(11.7/self.clock.norm if beam=='hot_bursts' else 1.))
    def parent(self,ec,beam,xyz,source_xy,age,sp,closure='nominal_retained',allocation='onsite'):
        case=self.case(float(ec),beam,closure,allocation);S=self.projection(xyz);edges=self.data['tile_edges'];n=len(edges)-1
        source_xy=np.asarray(source_xy);ix=np.clip(np.searchsorted(edges,source_xy[...,0],side='right')-1,0,n-1);iy=np.clip(np.searchsorted(edges,source_xy[...,1],side='right')-1,0,n-1)
        lookup=np.empty((n,n),int);lookup[self.data['tile_ix'],self.data['tile_iy']]=np.arange(len(self.data['tile_ix']))
        impulse=np.moveaxis(case['impulses'][sp][:,lookup[ix,iy]],0,-1)
        age=np.asarray(age);response=np.exp(-case['gamma']*age[...,None]/2)*np.sin(case['wd']*age[...,None])/case['wd']
        return np.sum(S*impulse*response,axis=-1)*(age>=0)
