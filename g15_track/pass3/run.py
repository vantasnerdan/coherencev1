"""Replay Pass 3 from the accepted inputs: python g15_track/pass3/run.py.

Only this directory is written. The model is a dispersion-resolved Gaussian
Wigner closure, not atom-by-atom dynamics or an entanglement calculation.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse,csv,hashlib,importlib.util,json,math,sys,time
from pathlib import Path
from functools import lru_cache
sys.dont_write_bytecode=True
import numpy as np
from scipy.linalg import expm
from scipy.special import roots_genlaguerre,roots_legendre,gamma
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'g15_track/pass2'))
spec=importlib.util.spec_from_file_location('frozen_pass2',ROOT/'g15_track/pass2/run.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
p1=old.p1
EV=p1.QE;M=old.M;HB=p1.HB*1e15;P=old.P
EC=old.EC;COEF=np.array([0.,.001,-.001,.01,-.01,.1,-.1])
DT=1.;SAVE=4
VOLUME=P['a_A']**3/4
BASE='2d85230f2580efbee4d7de5e8aece2c9540b3c1c'

def dump(name,value):
    def convert(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    (HERE/name).write_text(json.dumps(value,indent=2,allow_nan=False,default=convert)+'\n')

def table(name,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (HERE/name).open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=fields);w.writeheader();w.writerows(rows)

def rows(path):
    with path.open() as stream:return list(csv.DictReader(stream))

_GX,_GW=roots_genlaguerre(8,.5)
_MU,_MW=roots_legendre(8)
def moving_rate(lab,theta,ch,U):
    """Same accepted relative-flux kernel with explicit broadcast semantics."""
    lab,theta=np.broadcast_arrays(np.asarray(lab,float),np.asarray(theta,float))
    e=lab[...,None,None]/2
    thermal=theta[...,None,None]*_GX[:,None]/2000
    relative=e+thermal-2*np.sqrt(e*thermal)*_MU[None,:]
    return np.sum(p1.sigma(relative,ch,U)*np.sqrt(relative/e)*_GW[:,None]/gamma(1.5)*_MW[None,:]/2,axis=(-2,-1))

def quantum():
    summaries=json.loads((ROOT/'g10_quantum_site/output/exact/exact_summary.json').read_text())
    s=next(r for r in summaries if r['case']=='PdD_300K_joint')
    q=s['position']['mean_square_radius']/3
    v=HB**2*s['momentum']['mean_square_radius']/(3*M*M)
    coupled=next(r for r in rows(ROOT/'g10_quantum_site/output/tables/all_quantities.csv') if r['case']=='PdD_300K')
    return q,v,dict(position_axis_A2=q,momentum_axis_hbar2_Aminus2=s['momentum']['mean_square_radius']/3,
        velocity_axis_A2_fs2=v,symmetrized_qp_eV_fs=0.,uncertainty_over_hbar2=q*M*M*v/HB**2,
        exact_tail_source='g10_quantum_site/output/exact/PdD_300K_joint.npz',
        coupled_axis_A2=float(coupled['r2_A2'])/3,coupled_axis_SE_A2=float(coupled['r2_block_SE_A2'])/3,
        coupled_kinetic_eV=float(coupled['mean_kinetic_eV']),
        no_extra_zero_point_quantum=True)
Q0,V0,QINFO=quantum()

def pair_input(env):
    if env=='bulk':
        r=rows(ROOT/'g10_quantum_site/pdd-nearest-pair/pair-cut.csv')
        u=np.array([0.]+[float(x['displacement_per_D_A']) for x in r])
        energy=np.array([0.]+[float(x['total_approach_free_energy_eV']) for x in r])
        one=np.array([0.]+[float(x['matched_one_body_sum_eV']) for x in r])
        d=old.D0['bulk']
    else:
        r=old.VAC['rows'];u=np.array([0.]+[x['displacement_per_D_A'] for x in r])
        energy=np.array([0.]+[x['AB_eV'] for x in r]);one=np.array([0.]+[x['A_eV']+x['B_eV'] for x in r]);d=old.D0['vacancy']
    # Even directional normal form, fixed by the FIRST supported native point.
    # The higher native points remain actual energy bands, not extrapolated contact.
    ksite=one[1]/u[1]**2
    krel=energy[1]/(2*u[1]**2)
    return dict(u=u,energy=energy,one=one,d0=d,ksite=ksite,krel=krel,
        minimum_native_A=d-2*u[-1],native_lowpoint_A=d-2*u[1],
        connected_curvature_eV_A2=(energy[1]-one[1])/u[1]**2)
PAIR={e:pair_input(e) for e in ('bulk','vacancy')}

@lru_cache(maxsize=None)
def modes(env,coupled):
    # Whole reciprocal-cell symbol. A cage's internal COM/relative coordinates
    # come from its native two-body cut, not a bulk cos(k.r) at a nonlattice r.
    grid=np.array(np.meshgrid(*([np.arange(4)/4]*3),indexing='ij')).reshape(3,-1).T
    transform=np.linalg.inv(old.g3.CELL).T*2*np.pi/P['a_A'];wave=grid@transform.T
    J=np.zeros(len(wave))
    if coupled:
        for shell in P['shells']:
            delta=np.array(shell['periodic_alias_offsets_A'])
            J+=shell['J_sum_eV_A2']*np.cos(wave@delta.T).mean(axis=1)
    local=old.K if env=='bulk' else PAIR[env]['ksite']
    values=[];weights=[];characters=[];priors=[]
    sectors=(0,) if env=='bulk' else (1,-1)
    for sector in sectors:
        native=local-sector*PAIR[env]['connected_curvature_eV_A2'] if env!='bulk' else local
        eig=(native if coupled or env=='bulk' else local)-J
        if np.min(eig)<=0:raise ValueError('Quiet operator is not positive')
        unique,idx=np.unique(np.round(eig,12),return_inverse=True)
        count=np.bincount(idx);weight=count/(len(wave)*len(sectors))
        neighbor=np.cos(wave@(old.AXIS[env]*old.D0[env])) if env=='bulk' else np.full(len(wave),sector)
        values.extend(unique);weights.extend(weight)
        characters.extend(np.bincount(idx,weights=neighbor)/count)
        priors.extend(np.full(len(unique),math.sqrt(native/M)))
    values=np.array(values)
    return values,np.sqrt(values/M),np.array(weights),np.array(characters),local,wave,np.array(priors)

def spatial(ec,sp,s):
    tr=p1.track(ec,sp);ts,z,dep,src=p1.source(tr,sp);ts=ts*1e15;z=z*1e10
    R=s['radius']*1e10
    radial_edges=np.unique([0.,.8,3.,R/4,R/2,.75*R,R])
    radial_edges=radial_edges[(radial_edges>=0)&(radial_edges<=R)]
    rr=np.sqrt((radial_edges[1:]**2+radial_edges[:-1]**2)/2)
    area=np.pi*np.diff(radial_edges**2)
    shape=np.exp(-2*rr**2/R**2)*(1-rr**2/R**2)**2
    shape*=np.pi*R*R/(shape@area)
    pos=[];vol=[];start=[];duration=[];deposit=[];parent_segment=[];occupied=[]
    ze=np.linspace(0,z[-1],9)
    for j,(za,zb) in enumerate(zip(ze[:-1],ze[1:])):
        fraction=np.maximum(0.,np.minimum(z[1:],zb)-np.maximum(z[:-1],za))/np.diff(z)
        cell_source=fraction@src
        midpoint=(za+zb)/2
        native_row=int(np.clip(np.searchsorted(z,midpoint,side='right')-1,0,len(ts)-2))
        local_energy=float(np.interp(midpoint,z,np.r_[tr['initial_eV'],[row['energy_eV'] for row in tr['timeline']]]))
        pulse_duration=min(ts[native_row+1]-ts[native_row],20/math.sqrt(2*local_energy/M))
        arrival=float(np.interp(midpoint,z,ts))
        occupied_column=p1.nd_column(ec,zb*1e-10)-p1.nd_column(ec,za*1e-10)
        for ir,r in enumerate(rr):
            pos.append([r,0.,(za+zb)/2]);vol.append(area[ir]*(zb-za))
            start.append(arrival+r/10.);duration.append(pulse_duration)
            deposit.append(cell_source/(sp+1)*VOLUME/(np.pi*R*R*(zb-za))*shape[ir])
            parent_segment.append(j)
            occupied.append(area[ir]*1e-20*occupied_column)
    pos=np.array(pos);vol=np.array(vol);formula=vol/VOLUME
    weights=np.array(occupied);occupancy=weights/formula
    # Each coarse depth cell integrates every overlapping original source row.
    deposit=np.array(deposit)
    offsets=np.array([0.,.74/math.sqrt(2*tr['initial_eV']/M)]) if sp else np.array([0.])
    total_expected=src.sum(axis=0)
    source_error=np.max(np.abs((formula[:,None]*deposit).sum(axis=0)*(sp+1)-total_expected))
    terminal_mask=np.array(parent_segment)==7
    terminal_profile=np.zeros(len(pos))
    radial_shape=np.exp(-2*pos[:,0]**2/R**2)*(1-pos[:,0]**2/R**2)**2
    terminal_profile[terminal_mask]=radial_shape[terminal_mask]
    terminal_profile*=tr['terminal_unresolved_eV']/np.dot(formula,terminal_profile)
    return dict(position=pos,formula=formula,occupied=weights,occupancy=occupancy,
        start=np.array(start),duration=np.array(duration),deposit=deposit,offsets=offsets,
        parent_segment=np.array(parent_segment),tr=tr,source_error=source_error,
        total_eV=tr['initial_eV']*(sp+1),radius_A=R,terminal_profile=terminal_profile,
        terminal_start=ts[-1]+pos[:,0]/10.)

def source_rate(d,t,first_only=False):
    f=np.zeros(len(d['start']))
    for off in (d['offsets'][:1] if first_only else d['offsets']):
        left=np.maximum(d['start']+off,t-DT/2)
        right=np.minimum(d['start']+off+d['duration'],t+DT/2)
        f+=np.maximum(right-left,0.)/(DT*d['duration'])
    terminal=np.zeros(len(f))
    for off in (d['offsets'][:1] if first_only else d['offsets']):
        terminal+=((d['terminal_start']+off>=t-DT/2)&(d['terminal_start']+off<t+DT/2))/DT
    result=d['deposit']*f[:,None]
    result[:,0]+=d['terminal_profile']*terminal
    return result

def scalar_transfer(q,c,v,k,dt):
    # Exact constant-curvature symplectic map, physical v=p/M.
    omega=np.sqrt(k/M);a=np.cos(omega*dt);b=np.sin(omega*dt)/omega;d=-omega*np.sin(omega*dt)
    return a*a*q+2*a*b*c+b*b*v,a*d*q+(a*a+b*d)*c+a*b*v,d*d*q+2*d*a*c+a*a*v

def covariance_case(ec,sp,s,env,coupled):

    d=spatial(ec,sp,s);n=len(d['formula']);eig,omega,mw,cosine,local,wave,prior_omega=modes(env,coupled)
    nm=len(eig);nb=len(COEF);shape=(nb,n,nm)
    if env=='bulk':qinitial,vinitial=Q0,V0
    else:
        qinitial=HB/(2*M*prior_omega)/np.tanh(HB*prior_omega/(2*p1.KB*300))
        vinitial=prior_omega**2*qinitial
    initial_site_q=float(np.dot(np.broadcast_to(qinitial,len(mw)),mw))
    q=np.full(shape,qinitial);c=np.zeros(shape);v=np.full(shape,vinitial)
    qr=np.full((n,nm),qinitial);cr=np.zeros_like(qr);vr=np.full_like(qr,vinitial)
    qb=qr.copy();cb=cr.copy();vb=vr.copy()
    y=np.zeros((n,11));dy=np.zeros((nb,n,11));kprev=np.zeros((nb,n))
    A=p1.network(s)*1e-15;E=expm(A*DT)
    ext=np.zeros((22,22));ext[:11,:11]=A;ext[:11,11:]=np.eye(11);F=expm(ext*DT)[:11,11:]
    stop=max(d['start']+d['duration'])+d['offsets'][-1]+8*max(s['tau_ne'],s['tau_opt'],s['tau_end'])*1e15
    times=np.arange(0,math.ceil(stop/DT)*DT+DT,DT)
    loss=math.exp(-DT/(s['tau_opt']*1e15));wocc=d['occupied'];formula=d['formula'];occ=d['occupancy']
    work=np.zeros(nb);returned=np.zeros(nb);debit=np.zeros(nb);diss=np.zeros(nb);foil=np.zeros(nb)
    foil_local=np.zeros((nb,n));local_defect=0.
    peak=np.zeros((nb,n));peakwidth=np.zeros((nb,n));minwidth=np.full((nb,n),np.inf)
    minactual=np.inf;minK=np.inf;closure=0.;maxworkfrac=0.;maxunc=0.;max_quiet_mismatch=0.
    archive={k:[] for k in ('t','q','v','c','qr','vr','cr','relq','relv','relc','electron','optical','dyopt','anomalous','energy','work','diss','foil','carrier_delta')}
    mode_times=[];mode_values=[]
    # Finite records retain dispersion-resolved local covariance, not summed energies as phase.
    for step,t in enumerate(times):
        if step:
            previous_k=eig[None,None,:]+kprev[:,:,None]
            active_dt=np.clip(t-d['start'],0.,DT)[:,None]
            local_loss=np.exp(-active_dt/(s['tau_opt']*1e15))
            q,c,v=scalar_transfer(q,c,v,previous_k,active_dt)
            qr,cr,vr=scalar_transfer(qr,cr,vr,eig[None,:],active_dt)
            qb,cb,vb=scalar_transfer(qb,cb,vb,eig[None,:],active_dt)
            delta_loss=1.5*(M*(v-vr[None,:,:])+previous_k*q-eig[None,None,:]*qr[None,:,:]-kprev[:,:,None]*qb[None,:,:])
            lost=np.einsum('bnm,m->bn',(1-local_loss)*delta_loss,mw)*occ[None,:]
            q=local_loss*q+(1-local_loss)*qb;c=local_loss*c+(1-local_loss)*cb;v=local_loss*v+(1-local_loss)*vb
            qr=local_loss*qr+(1-local_loss)*qb;cr=local_loss*cr+(1-local_loss)*cb;vr=local_loss*vr+(1-local_loss)*vb
            u=source_rate(d,t-DT/2)
            old_opt=y[:,4].copy()
            y=y@E.T+u@F.T
            before=dy.sum(axis=2)
            dy=dy@E.T
            port=before-dy.sum(axis=2);foil_local+=port
            foil+=(port*formula).sum(axis=1)
            # Optical population is represented INSIDE the covariance. Convert
            # its actual carrier inflow, then clear the duplicate delta stock.
            reference_feed=np.maximum(y[:,4]-loss*old_opt,0.)/np.maximum(occ,1e-30)
            actual_feed=reference_feed[None,:]+dy[:,:,4]/np.maximum(occ,1e-30)
            if np.min(actual_feed)<-1e-12:raise ValueError('Negative actual optical covariance source')
            actual_feed=np.maximum(actual_feed,0.)
            q+=actual_feed[:,:,None]/(3*previous_k);v+=actual_feed[:,:,None]/(3*M)
            qr+=reference_feed[:,None]/(3*eig[None,:]);vr+=reference_feed[:,None]/(3*M)
            dy[:,:,4]=0.
            dy[:,:,8]+=.2*lost;dy[:,:,6]+=.8*lost
            diss+=(lost*formula).sum(axis=1)
        electron=y[:,3]+y[:,8]
        shift=COEF[:,None]*electron[None,:]*local
        dqk=shift-kprev
        qmean=np.einsum('bnm,m->bn',q,mw)
        event_work=1.5*dqk*qmean*occ[None,:]
        available=y[:,3]+y[:,8]+dy[:,:,3]+dy[:,:,8]
        positive=event_work>0
        if np.any(positive & (event_work>available+1e-12)):
            raise ValueError('Specified curvature exceeds local electronic work; requires a physical work-limited branch')
        frac=np.divide(event_work,available,out=np.zeros_like(event_work),where=available>1e-20)
        maxworkfrac=max(maxworkfrac,float(np.max(frac)))
        share=np.divide(y[:,3],electron,out=np.ones(n),where=electron>1e-30)
        dy[:,:,3]-=event_work*share;dy[:,:,8]-=event_work*(1-share)
        w=(event_work*formula).sum(axis=1);work+=w;debit+=np.maximum(w,0);returned-=np.minimum(w,0)
        k=eig[None,None,:]+shift[:,:,None]
        minK=min(minK,float(k.min()))
        if minK<=0:raise ValueError('Curvature crosses the bound-well chart; inverted harmonic equilibrium is not used')
        kprev=shift
        qm=np.einsum('bnm,m->bn',q,mw);vm=np.einsum('bnm,m->bn',v,mw);cm=np.einsum('bnm,m->bn',c,mw)
        qrm=qr@mw;vrm=vr@mw;crm=cr@mw
        de=1.5*(M*(vm-vrm)+np.einsum('bnm,m->bn',k*q-eig[None,None,:]*qr[None,:,:],mw))
        peak=np.maximum(peak,de);peakwidth=np.maximum(peakwidth,qm/qrm);minwidth=np.minimum(minwidth,qm/qrm)
        total_excess=(de*wocc).sum(axis=1)
        err=total_excess+(dy.sum(axis=2)*formula).sum(axis=1)+foil
        closure=max(closure,float(np.max(abs(err))))
        local_defect=max(local_defect,float(np.max(abs(de*occ[None,:]+dy.sum(axis=2)+foil_local))))
        minactual=min(minactual,float(np.min(y[None,:,:]+dy)))
        maxunc=max(maxunc,float(np.max(np.maximum(0.,HB**2/4-M*M*(q*v-c*c)))))
        max_quiet_mismatch=max(max_quiet_mismatch,float(np.max(abs((qb@mw)/initial_site_q-1))))
        if step%SAVE==0 or step==len(times)-1:
            # Physical complex anomalous covariance in each actual quiet k mode.
            anomaly=(M*omega[None,None,:]*(q-qr[None,:,:])-(v-vr[None,:,:])*M/omega[None,None,:]+2j*M*(c-cr[None,:,:]))/(2*HB)
            # Pair covariances are reconstructed from the SAME k field; translation cancels by 1-cos(k.r).
            filt=2*mw*(1-cosine)
            relq=np.einsum('bnm,m->bn',q,filt);relv=np.einsum('bnm,m->bn',v,filt);relc=np.einsum('bnm,m->bn',c,filt)
            for key,value in dict(t=t,q=qm,v=vm,c=cm,qr=qrm,vr=vrm,cr=crm,relq=relq,relv=relv,relc=relc,
                electron=electron,optical=y[:,4]/np.maximum(occ,1e-30),dyopt=dy[:,:,4]/np.maximum(occ,1e-30),
                anomalous=np.einsum('bnm,m->bn',anomaly,mw),energy=total_excess,work=work.copy(),diss=diss.copy(),foil=foil.copy(),
                carrier_delta=(dy*formula[None,:,None]).sum(axis=1)).items():archive[key].append(value.copy() if hasattr(value,'copy') else value)
            if len(mode_times)<64 and (step==0 or step>=len(times)*len(mode_times)/64):
                mode_times.append(t);mode_values.append(anomaly.astype(np.complex64))
    out={k:np.asarray(value) for k,value in archive.items()}
    out.update(position=d['position'],occupied=wocc,formula=formula,coef=COEF,omega=omega,mode_weight=mw,neighbor_cosine=cosine,
        localK=np.array(local),initial_q=np.array(initial_site_q),initial_v=np.array(vinitial),peak=peak,variance_ratio_max=peakwidth,variance_ratio_min=minwidth,
        source_clock=d['start'],source_duration=d['duration'],source_deposit=d['deposit'],source_offsets=d['offsets'],
        source_segment=d['parent_segment'],native_support_A=np.array(PAIR[env]['minimum_native_A']))
    out['mode_snapshot_times']=np.array(mode_times);out['mode_anomalous']=np.array(mode_values)
    out['anomalous']=out['anomalous'].astype(np.complex64)
    for field,reference in (('q','qr'),('v','vr'),('c','cr')):
        out['delta_'+field]=(out.pop(field)-out[reference][:,None,:]).astype(np.float32)
    for field in ('relq','relv','relc'):
        values=out.pop(field);out['reference_'+field]=values[:,0,:]
        out['delta_'+field]=(values-values[:,0:1,:]).astype(np.float32)
    name=f'cov_{ec:g}_{sp+1}_{s["name"]}_{env}_{int(coupled)}.npz';np.savez_compressed(HERE/name,**out)
    summary=[]
    for b,coef in enumerate(COEF):
        hot=int(np.argmax(peak[b]));last=out['anomalous'][-1,b]
        modepair=last@wocc
        denom=abs(last)@wocc
        summary.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,neighbor_on=coupled,coefficient_per_eV=coef,
            parent_initial_eV=d['total_eV'],occupied_D=wocc.sum(),formula_sites=formula.sum(),source_integral_error_eV=d['source_error'],
            occupied_column_relative_error=wocc.sum()/(p1.nd_column(ec,d['tr']['timeline'][-1]['depth_A']*1e-10)*np.pi*s['radius']**2)-1,
            min_curvature_eV_A2=minK,max_electron_work_fraction=maxworkfrac,signed_work_eV=work[b],positive_work_eV=debit[b],returned_work_eV=returned[b],
            excess_stored_final_eV=out['energy'][-1,b],optical_dissipated_delta_eV=diss[b],foil_port_delta_eV=foil[b],
            conservation_error_eV=closure,min_actual_carrier_eV_formula=minactual,uncertainty_violation_eV2_fs2=maxunc,
            local_conservation_error_eV_formula=local_defect,retarded_whole_conservation_bound_eV=local_defect*formula.sum(),
            maximum_site_excess_eV=peak[b,hot],hot_radius_A=d['position'][hot,0],hot_depth_A=d['position'][hot,2],
            maximum_variance_ratio=float(peakwidth[b].max()),minimum_variance_ratio=float(minwidth[b].min()),
            native_quiet_variance_mismatch=max_quiet_mismatch,
            midpoint_probe_final_anomalous_over_local=float(abs(modepair)/denom) if denom else 0.,
            source_archive=name))
    return summary

def analytic_controls():
    from scipy.integrate import solve_ivp
    w=old.OMEGA;h=1e-4;decay=100.;pulse=[]
    for rise in (.01,15.,100.,500.):
        def rhs(t,x):
            f=h*math.exp(-t/decay)*(-math.expm1(-t/rise))
            return (np.array([[0,1],[-w*w*(1+f),0]])@x.reshape(2,2)).ravel()
        S=solve_ivp(rhs,(0,2000),np.eye(2).ravel(),rtol=2e-10,atol=2e-12,max_step=.5).y[:,-1].reshape(2,2)
        beta=.5*(S[0,0]-S[1,1]+1j*(w*S[0,1]+S[1,0]/w))
        weak=w*h/2*(1/rise)/abs((1/decay+2j*w)*(1/decay+1/rise+2j*w))
        pulse.append(dict(rise_fs=rise,decay_fs=decay,exact_beta=abs(beta),weak_beta=weak,sudden_step_beta=h/4,
            relative_weak_error=abs(abs(beta)/weak-1),symplectic_det=np.linalg.det(S)))
    return dict(finite_pulse=pulse,quantum=QINFO,clock_period_fs=2*np.pi/old.OMEGA,
        reference='matched same initial quantum covariance under exactly the same quiet operator; no extra half quantum',
        bath='C -> exp(-dt/T1) C + (1-exp(-dt/T1)) Cquiet(t); Gaussian CP noise is explicit; pure dephasing is zero unless separately displayed',
        first_moment='identically zero under even K and zero-mean Gaussian bath; not coherent displacement or evidence of entanglement')

def clock_witness():
    """One native-row/coarse-clock comparison on the same frozen spatial area."""
    ec=2.5;sp=1;s=p1.SCENARIOS[0];d=spatial(ec,sp,s)
    tr=d['tr'];ts,z,dep,src=p1.source(tr,sp);ts*=1e15;z*=1e10
    nu=2*old.OMEGA;A=p1.network(s)*1e-15
    response=np.linalg.solve(1j*nu*np.eye(11)-A,d['deposit'].T).T
    electron=response[:,3]+response[:,8]
    ends=np.linspace(0,z[-1],9);segment=d['parent_segment']
    ta=np.interp(ends[segment],z,ts)+d['position'][:,0]/10.
    tb=np.interp(ends[segment+1],z,ts)+d['position'][:,0]/10.
    coarse=electron*np.exp(-1j*nu*(ta+tb)/2)*np.sinc(nu*(tb-ta)/(2*np.pi))
    fine=electron*np.exp(-1j*nu*(d['start']+d['duration']/2))*np.sinc(nu*d['duration']/(2*np.pi))
    terminal=np.linalg.solve(1j*nu*np.eye(11)-A,np.eye(11)[:,0])
    terminal=(terminal[3]+terminal[8])*d['terminal_profile']*np.exp(-1j*nu*d['terminal_start'])
    sibling=1+np.exp(-1j*nu*d['offsets'][1])
    coarse=(coarse+terminal)*sibling;fine=(fine+terminal)*sibling
    j=int(np.argmax(abs(fine)));w=d['occupied']/d['occupied'].sum()
    native=int(np.clip(np.searchsorted(z,d['position'][j,2])-1,0,len(ts)-2))
    return dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],radius_A=d['radius_A'],
        selected_position_A=d['position'][j],cohort_formula_sites=d['formula'][j],
        electron_deposit_eV_formula_unchanged=d['deposit'][j,3],
        coarse_local_duration_fs=tb[j]-ta[j],native_record_duration_fs=ts[native+1]-ts[native],
        actual_fine_local_duration_fs=d['duration'][j],actual_retarded_start_fs=d['start'][j],
        coarse_2omega_integral_eV_fs=abs(coarse[j]),fine_2omega_integral_eV_fs=abs(fine[j]),
        coarse_weak_beta_at_0p1_per_eV=old.OMEGA*.1/2*abs(coarse[j]),
        fine_weak_beta_at_0p1_per_eV=old.OMEGA*.1/2*abs(fine[j]),
        coarse_whole_footprint_beta=old.OMEGA*.1/2*abs(w@coarse),
        fine_whole_footprint_beta=old.OMEGA*.1/2*abs(w@fine),
        per_parent_source_integral_error_eV=d['source_error'],
        rule='same cohort energy/formula and cylinder area; original local row clock capped by accepted2nm local kernel/speed; no coarse energy divided by native-row volume')

def normalization_witness():
    result=[]
    for ec in (.25,2.5):
        for sp in (0,1):
            for s in p1.SCENARIOS:
                d=spatial(ec,sp,s)
                expected=p1.nd_column(ec,d['tr']['timeline'][-1]['depth_A']*1e-10)*np.pi*s['radius']**2
                relative=d['occupied'].sum()/expected-1
                if abs(relative)>1e-12 or d['start'].min()<0:raise ValueError('Source normalization or causal support')
                result.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],weighted_D=d['occupied'].sum(),
                    accepted_profile_D=expected,relative_error=relative,formula_sites=d['formula'].sum(),
                    source_integral_error_eV=d['source_error'],earliest_local_source_fs=d['start'].min()))
    return result

def bindings():
    names=['g15_track/pass1/inputs.json','g15_track/pass1/run.py','g15_track/pass1/bca.py','g15_track/pass1/yield_consumers.csv',
        'g15_track/pass1/whole_foil.csv','g15_track/pass1/joined_beam_history.npz','g3_lattice/coherence/lattice.py',
        'g10_quantum_site/output/exact/PdD_300K_joint.npz','g10_quantum_site/output/exact/PdD_300K_joint_g8_joint.json',
        'g10_quantum_site/output/exact/PdD_300K_basis8.npz','g10_quantum_site/output/exact/PdD_77K_joint.npz',
        'g10_quantum_site/output/exact/exact_summary.json','g10_quantum_site/output/tables/all_quantities.csv',
        'g10_quantum_site/pdd-nearest-pair/pair-cut.csv','g10_quantum_site/pd-vacancy-r4-six-D-pair/pair.json',
        'g15_track/pass2/run.py','g10_quantum_site/pdd-nearest-pair/pair-wells.json',
        'g13_rogue/pass5/readout.csv','g13_rogue/pass5/transport_ledger.csv','g13_rogue/pass5/inputs_ledger.md',
        'g10_quantum_site/pd-vacancy-r4-six-D-pair/native/reference.json','g12_residual/inputs_ledger.md']
    return dict(base=BASE,files=[dict(path=n,bytes=(ROOT/n).stat().st_size,sha256=hashlib.sha256((ROOT/n).read_bytes()).hexdigest()) for n in names])

def energy_cases(ec):
    summary=[]
    for sp in (0,1):
        for s in p1.SCENARIOS:
            for env in PAIR:
                for coupled in (False,True):summary+=covariance_case(ec,sp,s,env,coupled)
    return summary

def main(consume=False):
    if not consume:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=3) as pool:
            summary=[r for case in pool.map(energy_cases,EC) for r in case]
        table('phase_summary.csv',summary)
        controls=json.loads((HERE/'controls.json').read_text()) if (HERE/'controls.json').exists() else analytic_controls()
        controls['source_clock_comparison']=clock_witness()
        controls['normalization_endpoints']=normalization_witness()
        dump('controls.json',controls);dump('source_index.json',bindings())
    cspec=importlib.util.spec_from_file_location('pass3_consumers',HERE/'consumers.py')
    consumers=importlib.util.module_from_spec(cspec);cspec.loader.exec_module(consumers)
    result=consumers.run(sys.modules[__name__])
    dump('results.json',result);print(json.dumps(result,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--consume',action='store_true')
    args=parser.parse_args();main(args.consume)
