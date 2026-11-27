"""Four genuine global Gaussian modes, conditional on a frozen footprint.

This is a retained-band consumer, not a reconstruction or bound on omitted
modes. Its diagnostic work is not a new source or reservoir and is never
debited from electrons or added to the full source ledger a second time.
"""
import importlib.util
import math
from pathlib import Path
import numpy as np


def trial(z,L,count):
    q=np.pi*z/L
    columns=[np.cos(q),np.sin(q),np.cos(2*q),np.sin(2*q)]
    return np.column_stack(columns[:count])


def basis(m,d,ec,sp,s,env,coupled,receiver):
    owner,pos,w,delay,start=receiver.rule(m,d,ec,sp,s)
    N=float(w.sum());L=float(d['position'][:,2].max()*16/15)
    count=4 if env=='bulk' else 2
    phi=trial(pos[:,2],L,count);G=phi.T@(w[:,None]*phi)/N
    eig,U=np.linalg.eigh(G)
    if eig.min()<=1e-10:raise ValueError('Retained physical Gram is singular')
    white=U/np.sqrt(eig)
    sectors=(0,) if env=='bulk' else (1,-1)
    matrices=[];functions=[];qinitial=[];vinitial=[];mass=[];volume=[]
    nd=m.old.profile_nd(ec,pos[:,2])*1e-30
    volume_weight=w/nd
    local=float(d['localK'])
    for sector in sectors:
        native=local if env=='bulk' else local-sector*m.PAIR[env]['connected_curvature_eV_A2']
        cold=native if coupled else local
        H=cold*G
        if coupled:
            for shell in m.P['shells']:
                J=shell['J_sum_eV_A2'];offset=np.asarray(shell['periodic_alias_offsets_A'])
                H-=J*G
                for dz in offset[:,2]:
                    difference=trial(pos[:,2]+dz,L,count)-phi
                    H+=J/(2*len(offset)*N)*(difference.T@(w[:,None]*difference))
        values,rotation=np.linalg.eigh(white.T@H@white)
        if values.min()<=0:raise ValueError('Retained quiet Hamiltonian is unstable')
        functions.append(phi@white@rotation);matrices.extend(values)
        if env=='bulk':q,v=m.Q0,m.V0
        else:
            omega=math.sqrt(native/m.M)
            q=m.HB/(2*m.M*omega)/math.tanh(m.HB*omega/(2*m.p1.KB*300));v=omega*omega*q
        qinitial.extend([q]*count);vinitial.extend([v]*count)
        sector_count=N if env=='bulk' else N/2
        mass.append(m.M*sector_count*G)
        volume.append(phi.T@(volume_weight[:,None]*phi))
    return dict(owner=owner,position=pos,weight=w,delay=delay,start=start,N=N,L=L,
        functions=functions,K0=np.asarray(matrices),q0=np.asarray(qinitial),v0=np.asarray(vinitial),
        physical_mass_Gram=np.asarray(mass),geometric_volume_Gram=np.asarray(volume),
        raw_Gram_eigenvalues=eig,sector_counts=np.array([N if env=='bulk' else N/2]*len(sectors)))


def drivers(m,d,b,s,end):
    """Positive causal quadrature, with a CP-consistent population interpolant."""
    t=np.arange(0.,end+0.5,1.);T=len(t);K=np.zeros((T,4,4));population=np.zeros_like(K)
    tau=s['tau_opt']*1e15;h=float(d['t'][1]-d['t'][0])
    mask=np.isclose(d['t']/h,np.rint(d['t']/h),atol=1e-10)
    oldt=d['t'][mask];electron=d['electron'][mask];optical=d['optical'][mask]
    decay=math.exp(-h/tau);minimum_feed=0.
    size=b['functions'][0].shape[1]
    for node in range(len(d['occupied'])):
        chosen=np.flatnonzero(b['owner']==node)
        if not len(chosen):continue
        query=t[:,None]-b['delay'][chosen];causal=t[:,None]>=b['start'][chosen]
        e=np.interp(query.ravel(),oldt,electron[:,node],left=0.,right=0.).reshape(query.shape)*causal
        index=np.clip(np.floor(query/h).astype(int),0,len(oldt)-2)
        age=np.clip(query-oldt[index],0.,h)
        feed=optical[index+1,node]-decay*optical[index,node]
        minimum_feed=min(minimum_feed,float(feed.min()))
        if feed.min()<-1e-11:raise ValueError('Frozen optical reference has a negative inflow')
        # Only measured floating-point roundoff in a positive inflow is removed.
        feed=np.maximum(feed,0.)
        p=np.exp(-age/tau)*optical[index,node]+feed*(-np.expm1(-age/tau))/(-math.expm1(-h/tau))
        p*=causal*(query>=0)*(query<=oldt[-1])
        for sector,fn in enumerate(b['functions']):
            f=fn[chosen];outer=f[:,:,None]*f[:,None,:]*b['weight'][chosen,None,None]/b['N']
            sl=slice(sector*size,(sector+1)*size)
            K[:,sl,sl]+=np.einsum('tn,nij->tij',e,outer,optimize=True)
            population[:,sl,sl]+=np.einsum('tn,nij->tij',p,outer,optimize=True)
    K[t<b['start'].min()]=0.;population[t<b['start'].min()]=0.
    return t,K,population,minimum_feed


def transfer(K,M):
    values,U=np.linalg.eigh(K)
    if values.min()<=0:raise ValueError('Driven retained Hamiltonian crosses zero')
    omega=np.sqrt(values/M);cs=np.cos(omega);sn=np.sin(omega)
    def rotate(diag):return (U*diag[...,None,:])@U.swapaxes(-1,-2)
    c=rotate(cs);a=rotate(sn/omega);b=rotate(-omega*sn)
    S=np.concatenate((np.concatenate((c,a),axis=-1),np.concatenate((b,c),axis=-1)),axis=-2)
    inverse_root=rotate(1/np.sqrt(values))
    return S,inverse_root,float(values.min())


def anomalous(m,delta,omega):
    a=np.sqrt(m.M*omega/(2*m.HB));b=np.sqrt(m.M/(2*m.HB*omega))
    Q=delta[...,:4,:4];V=delta[...,4:,4:];C=delta[...,:4,4:]
    return Q*a[:,None]*a[None,:]-V*b[:,None]*b[None,:]+1j*(C*a[:,None]*b[None,:]+C.swapaxes(-1,-2)*b[:,None]*a[None,:])


def energy(C,K,M):
    return .5*(M*np.trace(C[...,4:,4:],axis1=-2,axis2=-1)+np.einsum('...ij,...ji->...',K,C[...,:4,:4]))


def persistence(t,signal,start):
    indices=np.flatnonzero(t>=start)
    if not len(indices) or signal[indices].max()==0:return 0.,False,0
    peak=indices[np.argmax(signal[indices])];later=np.flatnonzero((t>=t[peak])&(signal<=signal[peak]/math.e))
    return (float(t[later[0]]-t[peak]),False,int(peak)) if len(later) else (float(t[-1]-t[peak]),True,int(peak))


def case(m,d,ec,sp,s,env,coupled,receiver,local_phases):
    b=basis(m,d,ec,sp,s,env,coupled,receiver)
    end=float(np.load(m.HERE/f'field_pairs_{ec:g}_{sp+1}_{s["name"]}_{env}_{int(coupled)}.npz')['t_fs'][-1])
    t,E,P,minfeed=drivers(m,d,b,s,end)
    coefficients=np.r_[m.COEF,m.COEF,0.];mixed=np.r_[np.ones(7,bool),np.zeros(7,bool),True]
    K0=np.diag(b['K0']);omega=np.sqrt(b['K0']/m.M);nbranch=len(coefficients)
    C0=np.diag(np.r_[b['q0'],b['v0']]);C=np.broadcast_to(C0,(nbranch,8,8)).copy();bare=C0.copy()
    S0,_,_=transfer(K0[None],m.M);S0=S0[0]
    work=np.zeros(nbranch);loss=np.zeros(nbranch);previous=np.broadcast_to(K0,(nbranch,4,4)).copy()
    record=[C.copy()];works=[work.copy()];losses=[loss.copy()];defect=0.;minK=float(b['K0'].min())
    r=math.exp(-1/(s['tau_opt']*1e15));PSDcorrection=0.;uncertainty=math.inf;deterror=0.
    J=np.block([[np.zeros((4,4)),np.eye(4)],[-np.eye(4),np.zeros((4,4))]])
    scale=np.r_[np.sqrt(m.M*omega/m.HB),np.sqrt(m.M/(m.HB*omega))]
    mismatch=0.
    for j in range(1,len(t)):
        midpoint=(E[j-1]+E[j])/2
        drive=coefficients[:,None,None]*float(d['localK'])*midpoint
        drive[~mixed]*=np.eye(4)
        K=K0+drive
        work+=.5*np.einsum('bij,bji->b',K-previous,C[:,:4,:4])
        S,root,minimum=transfer(K,m.M);minK=min(minK,minimum)
        deterror=max(deterror,float(abs(np.linalg.det(S)-1).max()))
        C=S@C@S.swapaxes(-1,-2);bare=S0@bare@S0.T
        before=energy(C,K,m.M)
        C=r*C+(1-r)*bare
        loss+=before-energy(C,K,m.M)
        inflow=P[j]-r*P[j-1];inflow=(inflow+inflow.T)/2
        eig,U=np.linalg.eigh(inflow)
        tolerance=1e-11*max(float(np.linalg.norm(P[j])),1.)
        if eig.min()<-tolerance:raise ValueError('Projected reference noise is not positive')
        PSDcorrection+=float(-eig[eig<0].sum());inflow=(U*np.maximum(eig,0))@U.T
        C[:,:4,:4]+=root@inflow@root/3
        C[:,4:,4:]+=inflow/(3*m.M)
        newdrive=coefficients[:,None,None]*float(d['localK'])*E[j]
        newdrive[~mixed]*=np.eye(4);previous=K0+newdrive
        work+=.5*np.einsum('bij,bji->b',previous-K,C[:,:4,:4])
        balance=energy(C,previous,m.M)-energy(C[-1],K0,m.M)-work+work[-1]+loss-loss[-1]
        defect=max(defect,float(abs(balance).max()))
        physical=C*scale[None,:,None]*scale[None,None,:]
        uncertainty=min(uncertainty,float(np.linalg.eigvalsh(physical+.5j*J).min()))
        mismatch=max(mismatch,float(abs(np.diag(bare[:4,:4])/b['q0']-1).max()))
        record.append(C.copy());works.append(work.copy());losses.append(loss.copy())
    C=np.asarray(record);work=np.asarray(works);loss=np.asarray(losses)
    delta=C-C[:,-1:];A=anomalous(m,delta,omega);norm=np.linalg.norm(A,axis=(-1,-2))
    post=float(max(d['source_clock']+d['source_duration'])+d['source_offsets'][-1]+b['delay'].max())
    filename=f'band_{ec:g}_{sp+1}_{s["name"]}_{env}_{int(coupled)}.npz'
    np.savez_compressed(m.HERE/filename,t_fs=t,covariance_q_velocity=C,quiet_branch=nbranch-1,
        coefficients=coefficients,mixing_on=mixed,reference_omega_per_fs=omega,
        mode_functions=np.stack(b['functions']),quadrature_position_A=b['position'],occupied_weight=b['weight'],
        physical_mass_Gram_eV_fs2_Aminus2=b['physical_mass_Gram'],
        geometric_volume_Gram_A3=b['geometric_volume_Gram'],sector_oscillator_count=b['sector_counts'],
        signed_work_eV_per_retained_axis=work,changed_bath_outflow_eV_per_retained_axis=loss-loss[:,-1:],
        frozen_projected_electron_eV_formula=E,ordinary_population_eV_D=P)
    rows=[];pairs=[]
    for branch,(coefficient,mix) in enumerate(zip(coefficients[:-1],mixed[:-1])):
        time,censored,peak=persistence(t,norm[:,branch],post)
        match=next(x for x in local_phases if float(x['coefficient_per_eV'])==coefficient)
        qend=energy(delta[-1,branch],K0,m.M)
        row=dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],environment=env,neighbor_on=coupled,
            coefficient_per_eV=float(coefficient),mode_mixing_on=bool(mix),retained_real_modes=4,
            physical_D_in_footprint=b['N'],retained_fraction_of_axial_modes=4/b['N'],
            first_arrival_fs=float(b['start'].min()),complete_window_end_fs=float(t[-1]),
            post_deposition_start_fs=post,persistence_1e_fs=time,persistence_right_censored=censored,
            local_field_persistence_1e_fs=float(match['post_source_1e_persistence_fs']),
            peak_anomalous_Frobenius=float(norm[:,branch].max()),post_deposition_peak_anomalous=float(norm[peak,branch]),
            final_anomalous_Frobenius=float(norm[-1,branch]),final_H0_excitation_delta_eV=qend,
            signed_work_eV_per_retained_axis=float(work[-1,branch]),
            changed_bath_outflow_eV_per_retained_axis=float(loss[-1,branch]-loss[-1,-1]),
            K_final_relative_to_peak=float(np.linalg.norm(E[-1])/max(np.linalg.norm(E,axis=(1,2)).max(),1e-300)),
            minimum_retained_K_eV_A2=minK,maximum_energy_defect_eV=defect,
            minimum_quantum_uncertainty_eigenvalue=uncertainty,maximum_symplectic_det_error=deterror,
            zero_drive_covariance_difference=float(abs(delta[:,0]).max()),
            quiet_minimum_position_trace_A2=float(np.trace(C[:,-1,:4,:4],axis1=-2,axis2=-1).min()),
            native_quiet_variance_mismatch=mismatch,positive_noise_roundoff_repair_eV=PSDcorrection,
            frozen_reference_feed_minimum_eV=minfeed,
            support='four global modes only; axial cos/sin includes +/-q; no omitted-spectrum bound',
            state='conditional joint Gaussian from same native moments; global quiet clock, not a measured G10 intersite state',
            ledger='retained-band diagnostic only; no second electron debit or added energy reservoir')
        rows.append(row)
        for snapshot in np.unique([peak,len(t)-1]):
            for a in range(4):
                for c in range(4):pairs.append(dict(**{k:row[k] for k in ('E_cm_keV','species','scenario','environment','neighbor_on','coefficient_per_eV','mode_mixing_on')},
                    time_fs=float(t[snapshot]),mode_a=a,mode_b=c,omega_a_per_fs=float(omega[a]),omega_b_per_fs=float(omega[c]),
                    anomalous_real=float(A[snapshot,branch,a,c].real),anomalous_imag=float(A[snapshot,branch,a,c].imag)))
    return rows,pairs


def run(m):
    spec=importlib.util.spec_from_file_location('pass3_consumers',m.HERE/'consumers.py')
    consumer=importlib.util.module_from_spec(spec);spec.loader.exec_module(consumer)
    _,receiver=consumer.scientific_modules(m);rows=[];pairs=[];phases=m.rows(m.HERE/'phase_evolution.csv')
    ec=2.5;sp=1;s=next(x for x in m.p1.SCENARIOS if x['name']=='fast_acoustic')
    for env in m.PAIR:
        for coupled in (False,True):
            d=consumer.load_case(m,f'cov_{ec:g}_{sp+1}_{s["name"]}_{env}_{int(coupled)}.npz')
            local=[x for x in phases if float(x['E_cm_keV'])==ec and int(x['species'])==sp+1 and x['scenario']==s['name'] and x['environment']==env and x['neighbor_on']==str(coupled)]
            a,b=case(m,d,ec,sp,s,env,coupled,receiver,local);rows+=a;pairs+=b
    m.table('whole_mode_summary.csv',rows);m.table('whole_mode_pairs.csv',pairs)
    return dict(retained_real_modes=4,frozen_footprints=4,
        maximum_retained_energy_defect_eV=max(x['maximum_energy_defect_eV'] for x in rows),
        minimum_retained_uncertainty_eigenvalue=min(x['minimum_quantum_uncertainty_eigenvalue'] for x in rows),
        maximum_zero_drive_covariance_difference=max(x['zero_drive_covariance_difference'] for x in rows))


if __name__=='__main__':
    import json
    spec=importlib.util.spec_from_file_location('pass3_frozen',Path(__file__).with_name('run.py'))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    result=run(m);existing=json.loads((m.HERE/'results.json').read_text());existing['retained_global_modes']=result
    m.dump('results.json',existing);print(json.dumps(result,indent=2))
