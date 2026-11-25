"""Item 5 only: fixed-site random streams and paired finite-ensemble readouts.

Keeps the accepted loaded-surface counterfactual and all deterministic inputs.
No seed, stopping, transport or contact producer is invoked.
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import importlib.util
import json
import math
import multiprocessing
from pathlib import Path
import sys
import time
sys.dont_write_bytecode=True
import numpy as np
from numba import njit
from scipy.stats import beta, t as student_t

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'g13_rogue'
OUT=HERE/'output/statistics_addendum'
BASE='09b1b0779ed83f014d3e9ab21328fcd106f2a1bb'
TARGET_ENERGY=2.5
TARGET_ETA=.1
TARGET_THRESHOLD=.09328814888270068
ALPHA=.05
DT=2.
WINDOW=32800.
BURN_STEPS=8200
ALL_PHASES=64
FIRST_FOCUS_PHASES=128
MAX_FOCUS_PHASES=256
BATCH=16
AXIS=np.array([1.,1.,0.])/np.sqrt(2)
INPUTS={}


def source(path):
    path=ROOT/path
    data=path.read_bytes()
    INPUTS[str(path.relative_to(ROOT))]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
    return data


def emit(name,data):
    (OUT/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def table(name,data):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]))
        w.writeheader();w.writerows(data)


def load_model():
    spec=importlib.util.spec_from_file_location('g13_fixed_site_graph',ROOT/'g10_quantum_site/model.py')
    model_module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model_module)
    model=model_module.full()
    source('g10_quantum_site/model.py');source(model['source'])
    coefficient=json.loads(source('g13_rogue/output/inputs.json'))['coefficients']
    e0=coefficient['discrete_gap_eV'];chi=coefficient['Kerr_eV_per_quantum']
    kinetic=coefficient['kinetic_eV_A2'];mass=coefficient['mass_eV_fs2_A2']
    hbar=coefficient['hbar_eV_fs'];tau=coefficient['amplitude_lifetime_fs']
    L=e0*np.eye(32)-kinetic/e0*model['J']
    n=32
    coords=(model['positions']-model['positions'][0])/2.042
    if np.max(abs(coords-np.rint(coords)))>1e-8:
        raise ValueError('The fixed Fourier basis requires the accepted unrelaxed bulk graph')
    coords=np.rint(coords).astype(int)%4
    q=np.array([(a,b,c) for a in range(2) for b in range(4) for c in range(4)],int)
    # Exact fourth roots of unity fix the degenerate-mode gauge by construction.
    roots=np.array([1.+0j,1j,-1.+0j,-1j])
    W=roots[(coords@q.T)%4]/np.sqrt(n)
    diagonal=W.conj().T@L@W
    energy=np.diag(diagonal).real
    if np.max(abs(diagonal-np.diag(energy)))>1e-12 or energy.min()<=0:
        raise ValueError('Fixed Fourier functions do not diagonalize this operator')
    thermal=1/np.expm1(energy/(8.617333262145e-5*300.))
    covariance=(W*thermal)@W.conj().T
    root=((W*np.sqrt(thermal))@W.conj().T).real
    linear=(W*np.exp(-1j*energy*DT/hbar))@W.conj().T
    decay=np.exp(-DT/tau)
    positions=model['positions']
    pairs=np.array([(i,j) for i in range(n) for j in range(i+1,n)
        if np.isclose(np.linalg.norm((positions[j]-positions[i])-np.rint((positions[j]-positions[i])/8.168)*8.168),2.887824,atol=2e-5)])
    dr=positions[pairs[:,1]]-positions[pairs[:,0]]
    dr-=np.rint(dr/8.168)*8.168
    projection=(dr@AXIS)/np.linalg.norm(dr,axis=1)
    with (HERE/'output/packet_geometry.csv').open() as f:
        configs=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    source('g13_rogue/output/packet_geometry.csv')
    distance=positions-positions[model['site']]
    distance-=np.rint(distance/8.168)*8.168
    local_shape=np.exp(-np.sum(distance**2,axis=1)/(2*2.887824**2)).astype(complex)
    wave=np.exp(1j*(distance@AXIS)*.25)
    shapes=[]
    for c in configs:
        opt=c['optical_received_eV'];ac=c['acoustic_received_eV']
        shape=np.sqrt(opt)*local_shape/np.sqrt(np.vdot(local_shape,L@local_shape).real)
        shape+=np.sqrt(ac)*wave/np.sqrt(np.vdot(wave,L@wave).real)
        shape*=np.sqrt((opt+ac)/np.vdot(shape,L@shape).real)
        shapes.append(shape)
    return dict(L=L,W=W,q=q,energy=energy,thermal=thermal,covariance=covariance,
        root=root,linear=linear,decay=decay,bath_root=root*np.sqrt(1-decay**2),
        e0=e0,chi=chi,hbar=hbar,tau=tau,mass=mass,scale=2*np.sqrt(kinetic/e0),
        pairs=pairs,projection=projection,original_indices=model['original_indices'],
        configs=configs,shapes=np.asarray(shapes),thresholds=np.r_[0.,np.geomspace(1e-6,2.,72)])


def basis_witness(m):
    eigen,V=np.linalg.eigh(m['L'])
    rng=np.random.default_rng(8721)
    rotated=V.copy()
    for e in np.unique(np.round(eigen,12)):
        ix=np.flatnonzero(np.isclose(eigen,e,atol=5e-13,rtol=0))
        if len(ix)>1:
            rotation=np.linalg.qr(rng.normal(size=(len(ix),len(ix))))[0]
            rotated[:,ix]=V[:,ix]@rotation
    weights=np.sqrt(1/np.expm1(eigen/(8.617333262145e-5*300.)))
    ref=(V*weights)@V.T
    rotation_root=(rotated*weights)@rotated.T
    perm=np.arange(31,-1,-1)
    permutation_root=(V[:,perm]*weights[perm])@V[:,perm].T
    g=rng.normal(size=32)+1j*rng.normal(size=32)
    result=dict(fixed_basis='32 distinct FCC plane waves; qx=0,1 and qy,qz=0,1,2,3; exact fourth roots of unity',
        variates='complex Gaussian in source original site-index order, then positive covariance root',
        max_fourier_orthogonality_defect=float(abs(m['W'].conj().T@m['W']-np.eye(32)).max()),
        max_covariance_root_defect=float(abs(m['root']@m['root'].T-m['covariance']).max()),
        max_root_vs_symmetric_spectral_function=float(abs(m['root']-ref).max()),
        max_root_change_under_degenerate_rotation=float(abs(ref-rotation_root).max()),
        max_root_change_under_eigenpair_permutation=float(abs(ref-permutation_root).max()),
        fixed_site_draw_rotation_difference=float(abs(g@ref.T-g@rotation_root.T).max()),
        old_mode_assigned_draw_rotation_difference=float(abs((g*weights)@V.T-(g*weights)@rotated.T).max()),
        mode_energy_eV=m['energy'].tolist(),source_original_site_order=m['original_indices'].tolist(),
        random_stream='NumPy Philox raw uint64; per-phase SeedSequence([13991, phase_id, 0]); explicitly ordered Box-Muller transform',
        packet_phase_stream='independent Philox SeedSequence([13991, phase_id, stable_case_key, 1])',
        floating_point='Random variates no longer depend on eigensolver eigenvector signs/rotations/permutations. sin/cos/log and floating-point propagation can differ across libraries; no cross-host bitwise trajectory identity has been observed.')
    for key in ('max_fourier_orthogonality_defect','max_covariance_root_defect','max_root_vs_symmetric_spectral_function',
                'max_root_change_under_degenerate_rotation','max_root_change_under_eigenpair_permutation','fixed_site_draw_rotation_difference'):
        if result[key]>=1e-12:
            raise ArithmeticError((key,result[key]))
    if result['old_mode_assigned_draw_rotation_difference']<.01:
        raise ArithmeticError('Wrong-basis positive control did not discriminate')
    emit('basis.json',result)
    return result


def gaussian(generator,n=32):
    raw=generator.bit_generator.random_raw(2*n).reshape(2,n)
    u=((raw[0]>>11).astype(np.float64)+1.)*2.**-53
    theta=((raw[1]>>11).astype(np.float64))*(2*np.pi*2.**-53)
    radius=np.sqrt(-np.log(u))
    return radius*(np.cos(theta)+1j*np.sin(theta))


@njit(cache=False)
def local_bins(action,thresholds,e0,chi,scale,previous,hist,cross,maxima,support):
    for b in range(2):
        for p in range(action.shape[1]):
            for c in range(action.shape[2]):
                for i in range(32):
                    a=action[b,p,c,i]
                    harmonic=e0*a
                    own=harmonic+(.5*chi*a*a if b else 0.)
                    support[b,p,c,0]+=int(scale*math.sqrt(a)>.8)
                    support[b,p,c,1]+=int(a>1.)
                    for kind in range(2):
                        value=own if kind==0 else harmonic
                        index=np.searchsorted(thresholds,value,side='left')
                        hist[kind,b,p,c,index]+=1
                        old=previous[kind,b,p,c,i]
                        if old>=0 and index>old:
                            cross[kind,b,p,c,old]+=1
                            cross[kind,b,p,c,index]-=1
                        previous[kind,b,p,c,i]=index
                        maxima[kind,b,p,c]=max(maxima[kind,b,p,c],value)


@njit(cache=False)
def closing_bins(fields,action,velocity,thresholds,pairs,projection,mass,scale,chi,hbar,time_fs,
                 previous,hist,cross,maxima,peak):
    for b in range(2):
        for p in range(fields.shape[1]):
            for c in range(fields.shape[2]):
                for k in range(len(pairs)):
                    i,j=pairs[k,0],pairs[k,1]
                    common=(velocity[b,p,c,i]-velocity[b,p,c,j])*projection[k]
                    shift=0.
                    if b:
                        shift=scale*chi/hbar*(action[b,p,c,i]*fields[b,p,c,i].imag-
                                              action[b,p,c,j]*fields[b,p,c,j].imag)*projection[k]
                    closing=common+shift
                    actual=mass/4*max(closing,0.)**2
                    mapped=mass/4*max(common,0.)**2
                    if actual>maxima[0,b,p,c]:
                        peak[b,p,c,0]=time_fs
                        peak[b,p,c,1]=k
                        peak[b,p,c,2]=closing
                        peak[b,p,c,3]=common
                        peak[b,p,c,4]=shift
                        peak[b,p,c,5]=mass/4*common*common
                        peak[b,p,c,6]=mass/2*common*shift
                        peak[b,p,c,7]=mass/4*shift*shift
                        peak[b,p,c,8]=(action[b,p,c,i]+action[b,p,c,j])/2
                        peak[b,p,c,9]=max(action[b,p,c,i],action[b,p,c,j])
                    for kind in range(2):
                        value=actual if kind==0 else mapped
                        index=np.searchsorted(thresholds,value,side='left')
                        hist[kind,b,p,c,index]+=1
                        old=previous[kind,b,p,c,k]
                        if old>=0 and index>old:
                            cross[kind,b,p,c,old]+=1
                            cross[kind,b,p,c,index]-=1
                        previous[kind,b,p,c,k]=index
                        maxima[kind,b,p,c]=max(maxima[kind,b,p,c],value)


def worker(job):
    phase_ids,config_ids,m=job
    ns,nc=len(phase_ids),len(config_ids)
    rng=[np.random.Generator(np.random.Philox(np.random.SeedSequence([13991,int(p),0]))) for p in phase_ids]
    noise=np.empty((ns,32),complex)
    for p,g in enumerate(rng):noise[p]=gaussian(g)
    initial=noise@m['root'].T
    a=np.stack((initial.copy(),initial.copy()))
    for _ in range(BURN_STEPS):
        for p,g in enumerate(rng):noise[p]=gaussian(g)
        innovation=noise@m['bath_root'].T
        a[1]*=np.exp(-.5j*m['chi']*abs(a[1])**2*DT/m['hbar'])
        a=(a@m['linear'].T)*m['decay']+innovation
        a[1]*=np.exp(-.5j*m['chi']*abs(a[1])**2*DT/m['hbar'])
    fields=np.repeat(a[:,:,None,:],nc,axis=2)
    for c,global_id in enumerate(config_ids):
        if global_id==0:continue
        spec=m['configs'][global_id-1]
        encoded=json.dumps([spec['E_cm_keV'],spec['eta'],int(spec['zone'])],separators=(',',':')).encode()
        case_key=int.from_bytes(hashlib.sha256(encoded).digest()[:4],'little')
        for p,phase in enumerate(phase_ids):
            kick_rng=np.random.Generator(np.random.Philox(np.random.SeedSequence([13991,int(phase),case_key,1])))
            angle=float(kick_rng.random())*2*np.pi
            kick=(np.cos(angle)+1j*np.sin(angle))*m['shapes'][global_id-1]
            fields[:,p,c]+=kick
    bins=len(m['thresholds'])+1
    shape=(2,2,ns,nc,bins)
    hist_local=np.zeros(shape,np.int32);hist_close=np.zeros(shape,np.int32)
    cross_local=np.zeros(shape,np.int32);cross_close=np.zeros(shape,np.int32)
    previous_local=np.full((2,2,ns,nc,32),-1,np.int16)
    previous_close=np.full((2,2,ns,nc,len(m['pairs'])),-1,np.int16)
    max_local=np.zeros((2,2,ns,nc));max_close=np.zeros_like(max_local)
    peak=np.zeros((2,ns,nc,10));support=np.zeros((2,ns,nc,2),np.int32)
    for step in range(int(WINDOW/DT)+1):
        action=abs(fields)**2
        # Real-coordinate derivative of the shared L and loss. The nonlinear
        # clock term is separately inserted by closing_bins, never local K4 energy.
        velocity=m['scale']*((fields.imag@m['L'].T)/m['hbar']-fields.real/m['tau'])
        local_bins(action,m['thresholds'],m['e0'],m['chi'],m['scale'],previous_local,
                   hist_local,cross_local,max_local,support)
        closing_bins(fields,action,velocity,m['thresholds'],m['pairs'],m['projection'],m['mass'],m['scale'],m['chi'],m['hbar'],
                     step*DT,previous_close,hist_close,cross_close,max_close,peak)
        if step==int(WINDOW/DT):break
        for p,g in enumerate(rng):noise[p]=gaussian(g)
        innovation=noise@m['bath_root'].T
        fields[1]*=np.exp(-.5j*m['chi']*action[1]*DT/m['hbar'])
        fields=(fields@m['linear'].T)*m['decay']+innovation[None,:,None,:]
        fields[1]*=np.exp(-.5j*m['chi']*abs(fields[1])**2*DT/m['hbar'])
    return dict(phase_ids=np.asarray(phase_ids),config_ids=np.asarray(config_ids),hist_local=hist_local,hist_close=hist_close,
        cross_local=cross_local,cross_close=cross_close,max_local=max_local,max_close=max_close,peak=peak,support=support)


PHASE_AXES=dict(hist_local=2,hist_close=2,cross_local=2,cross_close=2,max_local=2,max_close=2,peak=1,support=1)


def combine(parts):
    return dict(phase_ids=np.concatenate([p['phase_ids'] for p in parts]),config_ids=parts[0]['config_ids'],
                **{key:np.concatenate([p[key] for p in parts],axis=axis) for key,axis in PHASE_AXES.items()})


def subset(data,ids):
    indexes=[int(np.flatnonzero(data['config_ids']==i)[0]) for i in ids]
    return dict(phase_ids=data['phase_ids'],config_ids=np.asarray(ids),
                **{key:np.take(data[key],indexes,axis=axis+1) for key,axis in PHASE_AXES.items()})


def interval_category(values,weights,alpha=ALPHA):
    """Simultaneous CP intervals for each nonzero signed category in each stratum."""
    n=len(values)
    categories=(-2,-1,1,2)
    individual_alpha=alpha/(len(weights)*len(categories))
    lower=upper=0.
    nonzero=[]
    for zone,weight in enumerate(weights):
        nonzero.append(int(np.count_nonzero(values[:,zone])))
        for category in categories:
            count=int(np.count_nonzero(values[:,zone]==category))
            lo=float(beta.ppf(individual_alpha/2,count,n-count+1)) if count else 0.
            hi=float(beta.ppf(1-individual_alpha/2,count+1,n-count)) if count<n else 1.
            lower+=weight*category*(lo if category>0 else hi)
            upper+=weight*category*(hi if category>0 else lo)
    return float(lower),float(upper),nonzero


def response(data,m,case,threshold,kind=0):
    ids=[i for i,g in enumerate(data['config_ids']) if g and
         m['configs'][g-1]['E_cm_keV']==case[0] and m['configs'][g-1]['eta']==case[1]]
    weights=np.array([m['configs'][data['config_ids'][i]-1]['flux_nm2_s']*
                      m['configs'][data['config_ids'][i]-1]['area_stratum_nm2']/32 for i in ids])
    hits=(data['max_close'][kind]>threshold).astype(np.int8)
    differences=hits[:,:,ids]-hits[:,:,[0]]
    paired=differences[1]-differences[0]
    return differences@weights,paired,weights


def target_readout(data,m):
    rates,categories,weights=response(data,m,(TARGET_ENERGY,TARGET_ETA),TARGET_THRESHOLD)
    difference=rates[1]-rates[0]
    lo,hi,nonzero=interval_category(categories,weights,alpha=ALPHA/2)
    return dict(phases=len(data['phase_ids']),threshold_eV=TARGET_THRESHOLD,
        nonlinear_minus_linear_per_D_s=float(difference.mean()),
        empirical_cluster_SE_per_D_s=float(difference.std(ddof=1)/np.sqrt(len(difference))),
        simultaneous_95_category_interval_per_D_s=[lo,hi],
        interval_error_budget='0.025 at each of at most two focused looks; overall simultaneous coverage at least 95%',
        nonzero_phase_counts_by_radial_stratum=nonzero,
        largest_stratum_single_phase_resolution_per_D_s=float(weights.max()/len(difference)),
        all_responses_zero=bool(not np.any(categories)),
        stop_interval_excludes_zero=bool(lo>0 or hi<0),stop_half_width_at_most_1000=bool((hi-lo)/2<=1000.))


def distributions(data,m,ensemble):
    phases=len(data['phase_ids']);times=int(WINDOW/DT)+1
    out=[]
    # Per-phase counts keep time/site correlations inside the independent cluster.
    occ_local=np.cumsum(data['hist_local'][...,::-1],axis=-1)[...,::-1][...,1:]/(32*times)
    occ_close=np.cumsum(data['hist_close'][...,::-1],axis=-1)[...,::-1][...,1:]/(len(m['pairs'])*times)
    up_local=np.cumsum(data['cross_local'],axis=-1)[...,:-1]/32
    up_close=np.cumsum(data['cross_close'],axis=-1)[...,:-1]/32
    for b,label in enumerate(('linear','nonlinear')):
        for c,g in enumerate(data['config_ids']):
            config=m['configs'][g-1] if g else dict(E_cm_keV=0.,eta=0.,zone=-1.)
            for k,threshold in enumerate(m['thresholds']):
                row=dict(ensemble=ensemble,model=label,E_cm_keV=config['E_cm_keV'],eta=config['eta'],zone=int(config['zone']),
                         threshold_eV=float(threshold),phases=phases,phase_first=int(data['phase_ids'][0]),phase_last=int(data['phase_ids'][-1]))
                for kind,name in enumerate(('own','common')):
                    local_hits=data['max_local'][kind,b,:,c]>threshold
                    close_hits=data['max_close'][kind,b,:,c]>threshold
                    row.update({name+'_local_occupancy':float(occ_local[kind,b,:,c,k].mean()),
                        name+'_closing_occupancy':float(occ_close[kind,b,:,c,k].mean()),
                        name+'_local_max_survival':float(local_hits.mean()),name+'_closing_max_survival':float(close_hits.mean()),
                        name+'_local_upcrossings_per_D_window':float(up_local[kind,b,:,c,k].mean()),
                        name+'_pair_upcrossings_per_D_window':float(up_close[kind,b,:,c,k].mean()),
                        name+'_closing_zero_exceedances':bool(not np.any(close_hits))})
                row['maximum_survival_one_count_resolution']=1/phases
                row['zero_hit_one_sided_95_probability_upper']=1-(ALPHA/2 if phases>ALL_PHASES else ALPHA)**(1/phases)
                out.append(row)
    return out,(occ_local,occ_close,up_local,up_close)


def rates_readout(data,m,ensemble,measures):
    occ_local,occ_close,up_local,up_close=measures
    output=[];maxout=[]
    cases=list(dict.fromkeys((m['configs'][g-1]['E_cm_keV'],m['configs'][g-1]['eta']) for g in data['config_ids'] if g))
    n=len(data['phase_ids']);window_s=WINDOW*1e-15
    for energy,eta in cases:
        ids=[i for i,g in enumerate(data['config_ids']) if g and m['configs'][g-1]['E_cm_keV']==energy and m['configs'][g-1]['eta']==eta]
        lambdas=np.array([m['configs'][data['config_ids'][i]-1]['flux_nm2_s']*m['configs'][data['config_ids'][i]-1]['area_stratum_nm2'] for i in ids])
        for k,threshold in enumerate(m['thresholds']):
            per_model,category,weights=response(data,m,(energy,eta),threshold)
            paired=per_model[1]-per_model[0]
            lo,hi,nonzero=interval_category(category,weights,alpha=ALPHA/2 if n>ALL_PHASES else ALPHA)
            for current in (.1,1.,10.):
                active=float(current*lambdas.sum()*window_s)
                for b,label in enumerate(('linear','nonlinear')):
                    maxout.append(dict(ensemble=ensemble,model=label,E_cm_keV=energy,eta=eta,current_scale=current,threshold_eV=float(threshold),phases=n,
                        ambient_closing_max_survival=float(np.mean(data['max_close'][0,b,:,0]>threshold)),
                        single_packet_closing_max_excess_per_D_s=float(current*per_model[b].mean()),
                        excess_cluster_SE_per_D_s=float(current*per_model[b].std(ddof=1)/np.sqrt(n)),
                        paired_nonlinear_minus_linear_per_D_s=float(current*paired.mean()),
                        paired_cluster_SE_per_D_s=float(current*paired.std(ddof=1)/np.sqrt(n)),
                        paired_95_category_lower_per_D_s=current*lo,paired_95_category_upper_per_D_s=current*hi,
                        nonzero_phase_counts_by_stratum=';'.join(map(str,nonzero)),
                        largest_stratum_single_phase_resolution_per_D_s=float(current*weights.max()/n),
                        all_paired_responses_zero=bool(not np.any(category))))
                    row=dict(ensemble=ensemble,model=label,E_cm_keV=energy,eta=eta,current_scale=current,threshold_eV=float(threshold),phases=n,
                        ambient_site_occupancy=float(occ_local[0,b,:,0,k].mean()),
                        ambient_pair_closing_occupancy=float(occ_close[0,b,:,0,k].mean()),
                        ambient_local_initiations_per_D_s=float(up_local[0,b,:,0,k].mean()/window_s),
                        ambient_pair_initiations_per_D_s=float(up_close[0,b,:,0,k].mean()/window_s),
                        finite_memory_active_mean=active,two_or_more_arrivals_probability=float(-np.expm1(-active)-active*np.exp(-active)),
                        overlap_tail_coefficient=None)
                    # Integer indexing above can transpose the phase/stratum axes;
                    # take() fixes the declared cluster-first layout explicitly.
                    for name,array,multiplier in (('site_occupancy',occ_local,window_s),('pair_occupancy',occ_close,window_s),
                        ('local_initiations',up_local,1.),('pair_initiations',up_close,1.)):
                        delta=np.take(array[0,b,:,:,k],ids,axis=1)-array[0,b,:,0,k,None]
                        clustered=delta@(current*lambdas)*multiplier
                        suffix='_per_D_s' if 'initiations' in name else ''
                        row['single_packet_'+name+'_excess'+suffix]=float(clustered.mean())
                        row['single_packet_'+name+'_cluster_SE'+suffix]=float(clustered.std(ddof=1)/np.sqrt(n))
                    output.append(row)
    return output,maxout


def compare(data,m,ensemble):
    records=[]
    for c,g in enumerate(data['config_ids']):
        spec=m['configs'][g-1] if g else dict(E_cm_keV=0.,eta=0.,zone=-1.)
        L=data['max_close'][0,0,:,c];N=data['max_close'][0,1,:,c]
        common=data['max_close'][1,1,:,c]
        paired=N-L
        n=len(L);se=float(paired.std(ddof=1)/np.sqrt(n))
        t_quantile=student_t.ppf(.9875 if n>ALL_PHASES else .975,n-1)
        interval=t_quantile*se
        common_delta=common-L
        common_se=float(common_delta.std(ddof=1)/np.sqrt(n))
        logratio=np.log(N/L)
        logratio_error=float(t_quantile*logratio.std(ddof=1)/np.sqrt(n))
        peak=data['peak'][1,:,c]
        own_peak=N
        if np.max(abs(peak[:,5:8].sum(axis=1)-own_peak))>1e-12:
            raise ArithmeticError('The same-time velocity decomposition failed')
        records.append(dict(ensemble=ensemble,E_cm_keV=spec['E_cm_keV'],eta=spec['eta'],zone=int(spec['zone']),phases=n,
            linear_largest_eV=float(L.max()),nonlinear_largest_eV=float(N.max()),largest_ratio=float(N.max()/L.max()),
            median_linear_window_max_eV=float(np.median(L)),median_nonlinear_window_max_eV=float(np.median(N)),
            median_common_clock_nonlinear_window_max_eV=float(np.median(common)),
            mean_linear_window_max_eV=float(L.mean()),mean_nonlinear_window_max_eV=float(N.mean()),
            paired_mean_maximum_difference_eV=float(paired.mean()),paired_mean_difference_SE_eV=se,
            paired_mean_difference_t95_eV=[float(paired.mean()-interval),float(paired.mean()+interval)],
            paired_mean_common_clock_minus_linear_eV=float((common-L).mean()),
            paired_common_clock_difference_SE_eV=common_se,
            paired_common_clock_difference_t95_eV=[float(common_delta.mean()-t_quantile*common_se),float(common_delta.mean()+t_quantile*common_se)],
            paired_geometric_mean_ratio=float(np.exp(logratio.mean())),
            paired_logratio_t95=[float(np.exp(logratio.mean()-logratio_error)),float(np.exp(logratio.mean()+logratio_error))],
            median_actual_to_common_clock_same_time_ratio=float(np.median(own_peak/np.maximum(peak[:,5],1e-300))),
            mean_common_kinetic_at_actual_peak_eV=float(peak[:,5].mean()),
            mean_clock_cross_term_at_actual_peak_eV=float(peak[:,6].mean()),mean_clock_square_at_actual_peak_eV=float(peak[:,7].mean()),
            mean_pair_action_at_actual_peak=float(peak[:,8].mean()),
            mean_common_harmonic_maximum_linear_eV=float(data['max_local'][1,0,:,c].mean()),
            mean_common_harmonic_maximum_nonlinear_eV=float(data['max_local'][1,1,:,c].mean()),
            mean_own_local_maximum_nonlinear_eV=float(data['max_local'][0,1,:,c].mean())))
    return records


def export(all_data,extra,m,producer_hash):
    focus_ids=[0]+[i+1 for i,c in enumerate(m['configs']) if c['E_cm_keV']==TARGET_ENERGY and c['eta']==TARGET_ETA]
    focus=combine([subset(all_data,focus_ids),extra])
    survival=[];rates=[];maxima=[];comparisons=[];phase_rows=[]
    for data,ensemble in ((all_data,'all_source_64'),(focus,'focused_'+str(len(focus['phase_ids'])))):
        rows,measures=distributions(data,m,ensemble);survival.extend(rows)
        r,mm=rates_readout(data,m,ensemble,measures);rates.extend(r);maxima.extend(mm)
        comparisons.extend(compare(data,m,ensemble))
        for b,label in enumerate(('linear','nonlinear')):
            for p,phase_id in enumerate(data['phase_ids']):
                for c,g in enumerate(data['config_ids']):
                    cfg=m['configs'][g-1] if g else dict(E_cm_keV=0.,eta=0.,zone=-1.)
                    peak=data['peak'][b,p,c];pair=m['pairs'][int(peak[1])]
                    phase_rows.append(dict(ensemble=ensemble,phase_id=int(phase_id),model=label,E_cm_keV=cfg['E_cm_keV'],eta=cfg['eta'],zone=int(cfg['zone']),
                        own_local_max_eV=float(data['max_local'][0,b,p,c]),common_harmonic_local_max_eV=float(data['max_local'][1,b,p,c]),
                        inward_kinetic_max_eV=float(data['max_close'][0,b,p,c]),common_clock_inward_max_eV=float(data['max_close'][1,b,p,c]),
                        peak_time_fs=float(peak[0]),pair_i_original=int(m['original_indices'][pair[0]]),pair_j_original=int(m['original_indices'][pair[1]]),
                        signed_closing_A_fs=float(peak[2]),common_closing_A_fs=float(peak[3]),nonlinear_clock_closing_A_fs=float(peak[4]),
                        same_time_common_kinetic_eV=float(peak[5]),same_time_clock_cross_eV=float(peak[6]),same_time_clock_square_eV=float(peak[7]),
                        pair_mean_action_at_peak=float(peak[8]),pair_max_action_at_peak=float(peak[9])))
    table('survivals.csv',survival);table('rates.csv',rates);table('current_maxima.csv',maxima);table('phase_maxima.csv',phase_rows)
    emit('comparison.json',comparisons)
    target=target_readout(focus,m)
    summary=dict(schema='g13-fixed-site-statistics-v1',base=BASE,scope='retained loaded-surface counterfactual only; not the fixed-profile track or deep receiver',
        all_source_phases=ALL_PHASES,focused_phases=len(focus['phase_ids']),target=target,
        interval='Pointwise in threshold. Simultaneous exact Clopper-Pearson intervals for all four nonzero signed categories in each radial stratum, Bonferroni family 95%; focused cases reserve alpha=0.025 per look for at most two looks, preserving at least 95% overall coverage. Correlation between strata is retained within each independent phase cluster.',
        zero_sample_semantics='A zero empirical SE is not a zero physical response; category intervals and one-phase resolution remain nonzero.',
        stochastic_model='Harmonic Bose OU bath in fixed site order; the nonlinear branch uses the same bath, not an exact nonlinear quantum canonical ensemble.',
        closing_score='Projected relative kinetic energy m max(c,0)^2/4 only; local quartic energy is not added to that score.',
        common_clock='Use the same nonlinear fields with only the shared linear L/loss coordinate derivative; removes the explicit chi*n phase-frequency term without changing field evolution.',
        nonlinear_peak_decomposition='At the same actual peak: K=m*c0^2/4 + m*c0*cchi/2 + m*cchi^2/4. These are not independently selected extrema.',
        local_scores='Both branches retain the common harmonic score Delta01*|a|^2; own nonlinear score adds chi*|a|^4/2.',
        window_maximum_t_intervals='Approximate independent phase-cluster t intervals; focused comparisons reserve half the error budget at each of at most two looks, not independent time/site samples.',
        physical_rate=None,remaining_uncertainty='Finite independent phase count, short supported envelope closure and finite time/radial support; no material enhancement is inferred.',
        support={})
    for data,name in ((all_data,'all_source_64'),(focus,'focused')):
        samples=len(data['phase_ids'])*len(data['config_ids'])*32*(int(WINDOW/DT)+1)
        summary['support'][name]={label:dict(coordinate_support_fraction=1-float(data['support'][b,:,:,0].sum()/samples),
            low_action_support_fraction=1-float(data['support'][b,:,:,1].sum()/samples)) for b,label in enumerate(('linear','nonlinear'))}
    emit('summary.json',summary)
    physical_fields=('E_cm_keV','eta','zone','radius_nm','annulus_area_nm2','area_stratum_nm2','arrival_fs',
                     'packet_harmonic_energy_eV','optical_received_eV','acoustic_received_eV','flux_nm2_s')
    emit('inputs.json',dict(base=BASE,files=INPUTS,producer_script_sha256=producer_hash,
        consumer_script_sha256=INPUTS['g13_rogue/statistics_addendum.py']['sha256'],
        coefficient_semantics='Unmodified accepted deterministic coefficients; no other producer replayed.',
        configurations=[{key:c[key] for key in physical_fields} for c in m['configs']],
        ignored_source_columns='Old packet_geometry.csv mean/min/max injection work belongs to the accepted 12-phase run, not this ensemble; it is not used or carried as a new statistic.',
        dt_fs=DT,window_fs=WINDOW,burn_steps=BURN_STEPS,
        independent_phase_ids_all=all_data['phase_ids'].tolist(),independent_phase_ids_focus=focus['phase_ids'].tolist()))
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=4)
    parser.add_argument('--basis-only',action='store_true')
    parser.add_argument('--readout-only',action='store_true')
    args=parser.parse_args()
    if not 1<=args.workers<=5:parser.error('Use one through five explicitly single-threaded workers')
    OUT.mkdir(exist_ok=True)
    source('g13_rogue/statistics_addendum.py');source('g13_rogue/output/statistics_addendum/contract.json')
    m=load_model()
    if args.basis_only:
        print(json.dumps(basis_witness(m)));return
    if not (OUT/'basis.json').exists():basis_witness(m)
    started=time.monotonic()
    producer_hash=INPUTS['g13_rogue/statistics_addendum.py']['sha256']
    if args.readout_only:
        previous=json.loads((OUT/'inputs.json').read_text())
        producer_hash=previous.get('producer_script_sha256',previous['files']['g13_rogue/statistics_addendum.py']['sha256'])
        with np.load(OUT/'ensemble_all.npz') as z:all_data={k:z[k] for k in z.files}
        with np.load(OUT/'ensemble_focus_extra.npz') as z:extra={k:z[k] for k in z.files}
    else:
        full_ids=list(range(len(m['configs'])+1))
        focus_ids=[0]+[i+1 for i,c in enumerate(m['configs']) if c['E_cm_keV']==TARGET_ENERGY and c['eta']==TARGET_ETA]
        jobs=[(list(range(i,i+BATCH)),full_ids,m) for i in range(0,ALL_PHASES,BATCH)]
        jobs += [(list(range(i,i+BATCH)),focus_ids,m) for i in range(ALL_PHASES,FIRST_FOCUS_PHASES,BATCH)]
        with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('fork')) as pool:
            result=list(pool.map(worker,jobs))
            all_data=combine(result[:ALL_PHASES//BATCH]);extra=combine(result[ALL_PHASES//BATCH:])
            first=target_readout(combine([subset(all_data,focus_ids),extra]),m)
            emit('sampling_decision.json',dict(first_focused_result=first,criterion='interval excludes zero or half-width <=1000 per D/s; otherwise extend only focused case to 256'))
            if not (first['stop_interval_excludes_zero'] or first['stop_half_width_at_most_1000']):
                jobs=[(list(range(i,i+BATCH)),focus_ids,m) for i in range(FIRST_FOCUS_PHASES,MAX_FOCUS_PHASES,BATCH)]
                extra=combine([extra,*list(pool.map(worker,jobs))])
        np.savez_compressed(OUT/'ensemble_all.npz',**all_data)
        np.savez_compressed(OUT/'ensemble_focus_extra.npz',**extra)
    summary=export(all_data,extra,m,producer_hash)
    print(json.dumps(dict(summary='g13_rogue/output/statistics_addendum/summary.json',target=summary['target'],
        all_source_phases=ALL_PHASES,focused_phases=summary['focused_phases'],elapsed_seconds=time.monotonic()-started,
        explicitly_single_threaded_workers=args.workers)))


if __name__=='__main__':main()
