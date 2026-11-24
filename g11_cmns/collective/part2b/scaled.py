"""Scaled shared-mode family with source-style intermediate-state loss.

hbar=omega=1. N identical spins use the exact J=N/2 permutation-symmetric
sector. A shifted oscillator window retains exact sqrt(n) factors. Loss is a
specified Feshbach-like absorber, not ordinary optical damping. Lost norm is
never renormalized into a successful transition.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys
import time

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.optimize import root, brentq
import quimb as qu
import qutip as qt

HERE=Path(__file__).resolve().parent
VERSION='scaled-spin-boson-v2'


def save(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')


def make(q,N,n,g,delta,loss,gamma,energy=0.,absorber_reference=None):
    radius=N*q+max(12,q)
    ks=np.arange(max(-n,-radius),radius+1)
    ms=np.arange(N+1)
    spin=np.sqrt((N-ms[:-1])*(ms[:-1]+1))
    sx=sp.diags([spin,spin],[-1,1],shape=(N+1,N+1),format='csr')
    raising=np.sqrt(n+ks[1:])
    ax=sp.diags([raising,raising],[-1,1],shape=(len(ks),len(ks)),format='csr')
    mm,kk=np.meshgrid(ms,ks,indexing='ij')
    # Combined parity is exact. The even-q destination is in the other block.
    keep=((mm+kk)%2==0).ravel()
    m=mm.ravel()[keep];k=kk.ravel()[keep]
    diagonal=delta*m+k
    h=sp.diags(diagonal,format='csr')+g/np.sqrt(N)*sp.kron(sx,ax,format='csr')[keep][:,keep]
    endpoints=(k==-q*m)
    reference_gap,reference_energy=absorber_reference if absorber_reference is not None else (delta,energy)
    distance=reference_gap*m+k-reference_energy
    boundary=(abs(distance)<=1e-10)&(~endpoints)
    absorbing=((~endpoints)&(~boundary)) if loss=='symmetric' else ((distance< -1e-10)&(~endpoints))
    boundary_gap=float(np.min(abs(distance[~endpoints])))
    boundary_count=int(boundary.sum())
    if loss=='none':absorbing[:]=False
    if loss=='lower_limit':
        # Exact infinite-width limit: remove the lossy intermediate states.
        retain=~absorbing;h=h[retain][:,retain]
        m=m[retain];k=k[retain];diagonal=diagonal[retain];endpoints=endpoints[retain]
        absorbing=np.zeros(len(m),bool)
    widths=float(gamma)*absorbing
    hc=h.astype(complex)-.5j*sp.diags(widths,format='csr')
    endpoint_indices=np.flatnonzero(endpoints)
    ai=int(np.flatnonzero((m==0)&(k==0))[0])
    where_b=np.flatnonzero((m==1)&(k==-q))
    bi=int(where_b[0]) if len(where_b) else None
    return dict(h=h,hc=hc,m=m,k=k,endpoints=endpoint_indices,A=ai,B=bi,widths=widths,
                edge=(abs(k)==radius),oscillator_radius=radius,dimension=len(m),diagonal=diagonal,
                mask_boundary_gap=boundary_gap,mask_boundary_count=boundary_count,
                absorber_reference=[float(reference_gap),float(reference_energy)])


def schur(q,N,n,g,delta,loss,gamma,energy,absorber_reference=None):
    model=make(q,N,n,g,delta,loss,gamma,energy,absorber_reference)
    p=model['endpoints'];mask=np.ones(model['dimension'],bool);mask[p]=False
    h=model['hc'];hpp=h[p][:,p].toarray();hpq=h[p][:,mask]
    # Every q-grid endpoint is treated equally and kept out of the absorber.
    rhs=h[mask][:,p].toarray()
    solved=spsolve(energy*sp.eye(mask.sum(),format='csc')-h[mask][:,mask].tocsc(),rhs)
    eff=hpp+hpq@solved
    a=int(np.flatnonzero(p==model['A'])[0]);b=int(np.flatnonzero(p==model['B'])[0])
    return eff,a,b,model


def resonance(q,N,n,g,loss,gamma,fixed_reference=False):
    if q%2==0 or g==0:
        return float(q),0.,0j,dict(optimizer='exact zero/opposite-parity control',residual=0.)
    phase=2*np.pi*np.arange(512)/512
    dressed=lambda d:float(np.sqrt(d*d+16*g*g*n/N*np.cos(phase)**2).mean())
    if dressed(0)>=q:raise ValueError('Requested dressed order below this coupling branch')
    guide=brentq(lambda d:dressed(d)-q,0.,float(q))
    x=np.array([guide,N/2*(guide-q)])
    absorber_reference=None
    def solve_at(kind,width,guess):
        def equations(z):
            eff,a,b,_=schur(q,N,n,g,z[0],kind,width,z[1],absorber_reference)
            return [eff[a,a].real-z[1],eff[b,b].real-z[1]]
        answer=root(equations,guess,tol=1e-10)
        residual=float(np.linalg.norm(equations(answer.x)))
        if residual>2e-7 or not 0<answer.x[0]<1.5*q:
            raise ValueError(f'Joint resonance unresolved: {kind} width{width}, {answer.x}, residual{residual}')
        return answer.x,residual
    x,residual=solve_at('none',0.,x)
    anchor=x.copy()
    if fixed_reference:absorber_reference=anchor.tolist()
    if abs(dressed(x[0])-q)>.15:
        raise ValueError('Wrong dressed-order branch')
    track=[dict(width=0.,gap=float(x[0]),energy=float(x[1]))]
    def advance(kind,previous_width,width,guess,depth=0):
        try:
            value,error=solve_at(kind,width,guess)
        except ValueError:
            if depth>=8 or width-previous_width<1e-3*max(1.,previous_width):
                raise
            middle=(previous_width+width)/2
            halfway,_=advance(kind,previous_width,middle,guess,depth+1)
            return advance(kind,middle,width,halfway,depth+1)
        track.append(dict(width=float(width),gap=float(value[0]),energy=float(value[1])))
        return value,error
    if loss!='none':
        target=1e5 if loss=='lower_limit' else gamma
        previous_width=0.
        for width in [w for w in (.1,.3,1.,3.,10.,30.,100.,300.,1000.,3000.,1e4,3e4,1e5) if w<target]+[target]:
            kind='lower' if loss=='lower_limit' else loss
            x,residual=advance(kind,previous_width,width,x)
            previous_width=width
        if loss=='lower_limit':x,residual=solve_at(loss,0.,x)
    eff,a,b,final_model=schur(q,N,n,g,x[0],loss,gamma,x[1],absorber_reference)
    return float(x[0]),float(x[1]),complex(eff[a,b]),dict(
        optimizer='classical dressed-order guide, finite BW correction, continuous loss-width branch',
        residual=residual,lossless_dressed_order=float(dressed(anchor[0])),
        classical_gap_guide=guide,lossless_anchor=anchor.tolist(),loss_branch=track,
        absorber_reference=absorber_reference,
        mask_boundary_rule='within1e-10 quantum of the reference energy: stable; all q-grid endpoints stable',
        mask_boundary_gap=final_model['mask_boundary_gap'],mask_boundary_count=final_model['mask_boundary_count'],
        all_q_grid_endpoint_diagonals_real=np.real(np.diag(eff)).tolist(),
        all_q_grid_endpoint_widths=(-2*np.imag(np.diag(eff))).tolist())


def spectrum(model):
    h=model['hc'].toarray()
    if not np.any(model['widths']):
        values,vectors=qu.eigh(h.real)
        left=vectors.T.astype(complex)
        conditioning=1.
    else:
        # H is complex symmetric. Its transpose-bilinear dual is exact away
        # from defective points. Measure the biorthogonality, do not assume it.
        values,vectors=la.eig(h,check_finite=False)
        overlaps=(vectors*vectors).sum(0)
        conditioning=float(np.max(1/abs(overlaps)))
        if np.min(abs(overlaps))<1e-8:raise ValueError('Near-defective absorber spectrum')
        left=vectors.T/overlaps[:,None]
        defect=float(np.linalg.norm(left@vectors-np.eye(len(values)),ord=np.inf))
        if defect>2e-5:raise ValueError('Biorthogonal spectral defect '+str(defect))
    residual=float(np.max(abs(h@vectors-vectors*values)))
    return values,vectors,left,conditioning,residual


def trajectory(model,values,vectors,left,source,target,times):
    coeff=left[:,source]
    phases=np.exp(-1j*values[:,None]*times[None,:])
    psi=vectors@(coeff[:,None]*phases)
    prob=abs(psi)**2
    surviving=prob.sum(0)
    if np.min(surviving)<-1e-8 or np.max(surviving)>1+5e-6:raise ValueError('Unphysical norm')
    target_prob=prob[target] if target is not None else np.zeros(len(times))
    end=model['endpoints'];endprob=prob[end]
    peak=float(target_prob.max())
    threshold=peak/2
    crossed=np.flatnonzero(target_prob>=threshold) if peak>1e-24 else np.array([],int)
    if len(crossed):
        j=int(crossed[0])
        t50=float(times[j])
        if j and target_prob[j]>target_prob[j-1]:
            t50=float(times[j-1]+(times[j]-times[j-1])*(threshold-target_prob[j-1])/(target_prob[j]-target_prob[j-1]))
    else:t50=None
    m_initial=model['m'][source];k_initial=model['k'][source]
    direction=1 if target is not None and model['m'][target]>m_initial else -1
    candidates=[]
    for ix,idx in enumerate(end):
        moved=int((model['m'][idx]-m_initial)*direction)
        if moved>0:candidates.append(dict(spin_steps=moved,phonons=abs(int(model['k'][idx]-k_initial)),
            peak_unconditioned_probability=float(endprob[ix].max())))
    largest=max((r['phonons'] for r in candidates if r['peak_unconditioned_probability']>=.01),default=0)
    rate=(threshold/t50) if t50 is not None and t50>0 else None
    return dict(target_peak=peak,time_to_half_peak=t50,event_metric_per_dimensionless_time=rate,
        largest_endpoint_quanta_at_one_percent=largest,endpoint_peaks=candidates,
        max_boundary_probability=float(prob[model['edge']].sum(0).max()),
        initial_reconstruction_error=float(np.max(abs(psi[:,0]-np.eye(model['dimension'])[:,source]))),
        minimum_survival=float(surviving.min()),final_sink_probability=float(max(0.,1-surviving[-1])),
        maximum_conditional_target=float(np.max(np.divide(target_prob,surviving,out=np.zeros_like(target_prob),where=surviving>1e-14))),
        times=times.tolist(),target_unconditioned=target_prob.tolist(),survival=surviving.tolist(),
        sink=(1-surviving).tolist(),mean_spin=(model['m']@prob).tolist(),mean_oscillator_offset=(model['k']@prob).tolist())


def key(config):
    return hashlib.sha256((VERSION+json.dumps(config,sort_keys=True)).encode()).hexdigest()


def run(config):
    started=time.perf_counter()
    q,N,n,g,loss,gamma=[config[k] for k in ('q','N','n','g_over_quantum','loss','gamma_over_omega')]
    delta,energy,v,match=resonance(q,N,n,g,loss,gamma,config.get('fixed_loss_reference',False))
    if config.get('detuning'):delta+=config['detuning']
    model=make(q,N,n,g,delta,loss,gamma,energy,match.get('absorber_reference'))
    values,vectors,left,cond,residual=spectrum(model)
    effective=abs(v)
    # A finite observation limit is explicit. A censored weak case is not zero.
    tmax=min(2*np.pi/max(effective,1e-12),1e7)
    times=np.unique(np.r_[0.,np.geomspace(1e-4,tmax,260),np.linspace(0,tmax,641)])
    up=trajectory(model,values,vectors,left,model['A'],model['B'],times)
    if model['B'] is not None:
        down=trajectory(model,values,vectors,left,model['B'],model['A'],times)
    else:
        down=dict(up);down['target_peak']=0.;down['event_metric_per_dimensionless_time']=None
        down['control_note']='opposite-parity target; reverse is identical by block symmetry'
    initial_up=float(model['diagonal'][model['A']]+n)
    initial_down=float(model['diagonal'][model['B']]+n) if model['B'] is not None else None
    return dict(config=config,key=key(config),model_version=VERSION,bare_gap_over_quantum=delta,
        resonant_exchange_quanta=q,BW_reference_energy=energy,indirect_BW_coupling_real_imag=[v.real,v.imag],
        indirect_splitting_scale=2*abs(v),rate_semantics='half of peak unconditioned endpoint probability divided by first time reaching that level; a finite-preparation metric, not deposited power or fusion rate',
        observation_limit=tmax,time_censored=bool(tmax>=1e7),
        initial_bare_total_energy_quanta=dict(up=initial_up,down=initial_down),
        preparation='permutation-symmetric |J=N/2,m=-J,n> and reverse one-bright-excitation |J,m=-J+1,n-q>; zero nuclear excitation energy in up, one in down',
        loss_semantics='same absorber width on both intermediate sides, or only bare-product H_j<E_BW; all q-grid endpoints stable; sink receives lost norm, no renormalization',
        resonance=match,dimension=model['dimension'],oscillator_radius=model['oscillator_radius'],
        spectral_residual=residual,biorthogonal_condition=cond,up=up,down=down,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),runtime_seconds=time.perf_counter()-started)


def configs(suite):
    output=[]
    orders=(3,5,15) if suite=='probe' else (3,5,9,15,21)
    for q in orders:
        base=dict(q=q,N=2,n=2000*q,g_over_quantum=.085*np.sqrt(q/150),loss='none',gamma_over_omega=0.)
        if suite=='probe':
            for loss,width in (('none',0.),('symmetric',100.),('lower',100.),('lower_limit',0.)):
                output.append(dict(base,loss=loss,gamma_over_omega=width))
        else:
            family=[base]
            family += [dict(base,N=N) for N in (4,8)]
            family += [dict(base,n=f*q) for f in (500,1000,4000)]
            family += [dict(base,g_over_quantum=base['g_over_quantum']*f) for f in (.5,.75,1.25)]
            for x in family:
                for loss,width in (('none',0.),('symmetric',100.),('lower',100.),('lower_limit',0.)):
                    output.append(dict(x,loss=loss,gamma_over_omega=width))
            for width in (1.,10.,1000.):
                for loss in ('symmetric','lower'):output.append(dict(base,loss=loss,gamma_over_omega=width))
            output.append(dict(base,detuning=.2))
    # No sweep is needed for exact zero/opposite-parity corners.
    if suite=='full':
        output += [dict(q=4,N=2,n=8000,g_over_quantum=.04,loss='none',gamma_over_omega=0.),
                   dict(q=3,N=2,n=6000,g_over_quantum=0.,loss='none',gamma_over_omega=0.)]
    return [dict(c,fixed_loss_reference=True) if c['loss'] in ('symmetric','lower','lower_limit') else c for c in output]


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--suite',choices=('probe','full'),default='probe')
    parser.add_argument('--repair',action='store_true',help='Retry only previously unresolved cases with adaptive branch continuation')
    args=parser.parse_args();HERE.mkdir(parents=True,exist_ok=True)
    path=HERE/'runs.jsonl';rows=[]
    if path.exists():
        rows=[json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    latest={row['key']:row for row in rows}
    done={identity for identity,row in latest.items() if not (args.repair and 'failure' in row)}
    for config in configs(args.suite):
        if key(config) in done:continue
        try:
            result=run(config)
        except ValueError as error:
            result=dict(config=config,key=key(config),model_version=VERSION,failure=str(error),
                supported=False,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        with path.open('a') as f:f.write(json.dumps(result,allow_nan=False)+'\n')
        rows.append(result);done.add(result['key'])
        if 'failure' in result:
            print(json.dumps({'unsupported_case':config,'reason':result['failure']}),flush=True)
            continue
        print(json.dumps({**config,'delta':result['bare_gap_over_quantum'],'BW':result['indirect_splitting_scale'],
            'up_peak':result['up']['target_peak'],'up_rate':result['up']['event_metric_per_dimensionless_time'],
            'sink':result['up']['final_sink_probability'],'edge':result['up']['max_boundary_probability'],
            'dimension':result['dimension']}),flush=True)
    save(HERE/'environment.json',dict(model_version=VERSION,python=sys.version.split()[0],numpy=np.__version__,
        quimb=qu.__version__,qutip=qt.__version__,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sources=['https://jcmns.org/api/v1/articles/72149-energy-exchange-in-the-lossy-spin-boson-model.pdf',
                 'https://arxiv.org/html/0709.3557'],completed_cases=len(rows)))


if __name__=='__main__':main()
