"""Fit only the declared finite-family coordinates; retain log extrapolation.

Decimal path sums protect the weak high-order cancellation control. They are
leading amplitudes, not unconditioned population rates. Full-H propagation is
read from scaled.py outputs, not replaced by this recurrence.
"""
from pathlib import Path
from decimal import Decimal as D, localcontext
import csv
import json
import math
import hashlib
import sys
import numpy as np
from scipy.special import gammaln
import scipy.sparse as sp
import qutip as qt
import quimb as qu
import scaled

HERE=Path(__file__).resolve().parent


def precise_coefficient(q,N,lower,precision=100):
    with localcontext() as context:
        context.prec=precision
        levels={0:D(1)}
        for step in range(1,q+1):
            new={}
            for j,amplitude in levels.items():
                transitions=[]
                if j<N:transitions.append((j+1,D((N-j)*(j+1)).sqrt()))
                if j>0:transitions.append((j-1,D(j*(N-j+1)).sqrt()))
                for target,factor in transitions:
                    if step==q:
                        if target!=1:continue
                        denominator=D(1)
                    else:
                        if lower and q*target-step<0:continue
                        denominator=D(step-q*target)
                    new[target]=new.get(target,D(0))+amplitude*factor/denominator
            levels=new
        result=levels.get(1,D(0))/(D(N)**(D(q)/2))
        return str(result),float(abs(result).ln()) if result else None


def line_fit(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    coefficients=np.polyfit(x,y,1);pred=np.polyval(coefficients,x)
    return dict(slope=float(coefficients[0]),intercept=float(coefficients[1]),
        max_abs_log_residual=float(abs(pred-y).max()),rms_log_residual=float(np.sqrt(np.mean((pred-y)**2))),
        data_points=len(x),x_range=[float(x.min()),float(x.max())])


def family_tables(rows):
    table=[]
    for row in rows:
        if 'failure' in row:continue
        for direction in ('up','down'):
            r=row[direction];cfg=row['config']
            table.append(dict(**cfg,direction=direction,bare_gap=row['bare_gap_over_quantum'],
                endpoint_q=row['resonant_exchange_quanta'],peak_probability=r['target_peak'],
                metric_rate=r['event_metric_per_dimensionless_time'] if row.get('order_supported',True) else None,
                raw_bare_bin_metric=r['event_metric_per_dimensionless_time'],order_supported=row.get('order_supported',True),
                time_half=r['time_to_half_peak'],
                sink=r['final_sink_probability'],minimum_survival=r['minimum_survival'],
                largest_endpoint_quanta=r['largest_endpoint_quanta_at_one_percent'],
                edge_probability=r['max_boundary_probability'],time_censored=row['time_censored'],
                spectral_residual=row['spectral_residual'],dimension=row['dimension']))
    names=list(dict.fromkeys(k for row in table for k in row))
    with (HERE/'transfer.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=names);writer.writeheader();writer.writerows(table)
    fits=[]
    for q in (3,5,9,15,21):
        base_n=2000*q;base_g=.085*np.sqrt(q/150)
        for loss in ('none','symmetric','lower','lower_limit'):
            for direction in ('up','down'):
                for axis in ('N','n','g_over_quantum'):
                    data=[]
                    for row in table:
                        if row['q']!=q or row['loss']!=loss or row['direction']!=direction:continue
                        if row.get('detuning'):continue
                        if loss in ('symmetric','lower') and row['gamma_over_omega']!=100.:continue
                        if axis!='N' and row['N']!=2:continue
                        if axis!='n' and row['n']!=base_n:continue
                        if axis!='g_over_quantum' and abs(row['g_over_quantum']/base_g-1)>1e-10:continue
                        if row['metric_rate'] is None or row['peak_probability']<1e-12:continue
                        if row['edge_probability']>1e-5:continue
                        if row['time_censored']:continue
                        data.append(row)
                    unique={row[axis]:row for row in data};data=list(unique.values())
                    if len(data)<3:continue
                    fit=line_fit([math.log(row[axis]) for row in data],[math.log(row['metric_rate']) for row in data])
                    fits.append(dict(q=q,loss=loss,direction=direction,axis=axis,**fit,
                        axis_values=[row[axis] for row in data],rate_values=[row['metric_rate'] for row in data],
                        interpretation='local finite-preparation occupancy metric; no physical foil/fusion rate',
                        fit_is_locally_powerlike=fit['max_abs_log_residual']<.35))
    return table,fits


def consumer_proof(rows):
    row=next(r for r in rows if r['config']['q']==3 and r['config']['N']==2 and r['config']['loss']=='lower' and r['config']['gamma_over_omega']==100.)
    cfg=row['config'];m=scaled.make(cfg['q'],cfg['N'],cfg['n'],cfg['g_over_quantum'],row['bare_gap_over_quantum'],cfg['loss'],cfg['gamma_over_omega'],row['BW_reference_energy'],row['resonance'].get('absorber_reference'))
    dim=m['dimension'];sink=dim;h=sp.block_diag([m['h'],sp.csr_matrix((1,1))],format='csr')
    jumps=[qt.Qobj(sp.csr_matrix(([math.sqrt(width)],([sink],[j])),shape=(dim+1,dim+1))) for j,width in enumerate(m['widths']) if width]
    source=np.eye(dim+1)[m['A']];times=np.linspace(0,30,81)
    target_op=qt.Qobj(sp.csr_matrix(([1.],([m['B']],[m['B']])),shape=h.shape))
    active_op=qt.Qobj(sp.diags(np.r_[np.ones(dim),0.],format='csr'))
    solution=qt.mesolve(qt.Qobj(h),qt.Qobj(source),times,jumps,e_ops=[target_op,active_op],
        options={'rtol':2e-10,'atol':1e-13,'nsteps':100000})
    states=qu.expm_multiply(-1j*m['hc'],source[:dim],start=0,stop=float(times[-1]),num=len(times))
    probability=abs(np.asarray(states))**2;target=probability[:,m['B']];survival=probability.sum(1)
    discrepancy=max(float(abs(target-solution.expect[0]).max()),float(abs(survival-solution.expect[1]).max()))
    wrong=target/np.maximum(survival,1e-300)
    if discrepancy>3e-8:raise AssertionError(discrepancy)
    if np.max(abs(wrong-target))<.01:raise AssertionError('Control did not discriminate conditioning')
    return dict(qutip_full_sink_vs_quimb_sector_max_abs=discrepancy,
        wrong_survival_conditioning_max_abs=float(np.max(abs(wrong-target))),
        maximum_unconditioned_destination=float(target.max()),minimum_survival=float(survival.min()),
        norm_plus_sink_identity=True,case=cfg,proof_time_window=times.tolist())


def branch_witness(rows):
    selected=[r for r in rows if 'failure' not in r and r['config']['q']%2 and r['config']['g_over_quantum']>0
              and (r['config']['N']==2 or r['config']['q'] in (3,5,9,15))]
    selected=[r for r in selected if r['config']['gamma_over_omega'] in (0.,100.) and not r['config'].get('detuning')]
    witnesses=[]
    group=lambda c:tuple(c[k] for k in ('q','N','n','g_over_quantum'))
    common_by_group={group(r['config']):min(s['observation_limit'] for s in selected if group(s['config'])==group(r['config'])) for r in selected}
    old_path=HERE/'branch_proof.json'
    old=json.loads(old_path.read_text()) if old_path.exists() else {}
    prior=old.get('cases',old.get('q15',[]))
    cache_path=HERE/'cluster_cache.json'
    if cache_path.exists():prior+=json.loads(cache_path.read_text())['cases']
    old_by_key={scaled.key(x['config']):x for x in prior if not x['config'].get('detuning')}
    for row in selected:
        cfg=row['config'];track=list(row['resonance']['loss_branch']);reference=row['resonance'].get('absorber_reference')
        q=cfg['q'];N=cfg['N'];n=cfg['n'];common=common_by_group[group(cfg)]
        previous_result=old_by_key.get(scaled.key(cfg))
        if previous_result is not None:
            witnesses.append(previous_result)
            continue
        if cfg['loss']=='lower_limit':
            track.append(dict(width='infinite',gap=row['bare_gap_over_quantum'],energy=row['BW_reference_energy']))
        first=track[0]
        anchor=scaled.make(q,N,n,cfg['g_over_quantum'],first['gap'],'none',0.,first['energy'],reference)
        coordinates=list(zip(anchor['m'],anchor['k']));coordinate_index={x:i for i,x in enumerate(coordinates)}
        previous=None;steps=[]
        for point in track:
            kind='lower_limit' if point['width']=='infinite' else ('none' if point['width']==0 else ('lower' if cfg['loss']=='lower_limit' else cfg['loss']))
            width=0. if point['width']=='infinite' else point['width']
            model=scaled.make(q,N,n,cfg['g_over_quantum'],point['gap'],kind,width,point['energy'],reference)
            count=N+1
            if model['dimension']<=256:
                values,vectors,_,_,_=scaled.spectrum(model)
            else:
                values,vectors=scaled.la.eig(model['hc'].toarray(),check_finite=False)
                defect=float(np.max(abs(model['hc']@vectors-vectors*values)))
                if defect>1e-7:raise ValueError('Cluster eigensolve residual '+str(defect))
            embedded=np.zeros((len(coordinates),vectors.shape[1]),complex)
            ids=[coordinate_index[x] for x in zip(model['m'],model['k'])]
            embedded[ids]=vectors
            if previous is None:
                choice=np.argsort(abs(values-point['energy']))[:count]
            else:
                score=np.sum(abs(previous.conj().T@embedded)**2,axis=0)
                choice=np.argsort(score)[-count:]
            current=np.linalg.qr(embedded[:,choice])[0]
            overlap=float(np.linalg.svd(previous.conj().T@current,compute_uv=False).min()) if previous is not None else 1.
            steps.append(dict(width=point['width'],bare_gap=point['gap'],reference_energy=point['energy'],
                cluster_poles_real_imag=[[float(values[j].real),float(values[j].imag)] for j in choice],
                prior_cluster_min_cosine=overlap,cluster_center_error=float(abs(values[choice].real.mean()-point['energy']))))
            previous=current
            if overlap<1e-4:break
        same_window={}
        for direction in ('up','down'):
            r=row[direction];time=np.asarray(r['times']);prob=np.asarray(r['target_unconditioned'])
            inside=time<=common;peak_index=int(np.argmax(np.where(inside,prob,-1.)))
            same_window[direction]=dict(window=common,peak_probability=float(prob[inside].max()),
                sink_at_window=float(np.interp(common,time,r['sink'])),
                mean_physical_oscillator_offset_at_peak=float(r['mean_oscillator_offset'][peak_index]),
                mean_nuclear_steps_at_peak=float(r['mean_spin'][peak_index]),
                target_joint_photon_change=q,initial_physical_Fock_n=n if direction=='up' else n-q,
                gauge='original Fock basis; no displaced occupation used as source energy')
        witnesses.append(dict(config=cfg,lossless_dressed_order=row['resonance']['lossless_dressed_order'],
            mask_boundary_gap=row['resonance'].get('mask_boundary_gap'),
            mask_boundary_count=row['resonance'].get('mask_boundary_count'),
            same_physical_dimensionless_window=same_window,individually_scaled_window=row['observation_limit'],
            full_spectrum_cluster_track=steps))
        scaled.save(cache_path,dict(cases=witnesses))
    return dict(cases=witnesses,rejected_branch=dict(bare_gap=14.325519,dressed_order=17,
        reason='nonzero bare15-photon bin did not identify the dressed exchange order; excluded from fits and selected output',
        retained_diagnostic='rejected17-quantum probe retained locally; excluded from the draft'),
        source_benchmark='0709.3557 uses bare11,dressed15 and ncrit about27000; this is a convention check, not a claim to reproduce that exact parameter point')


def rate_extensions(table,fits):
    laws=[]
    for loss in ('none','symmetric','lower','lower_limit'):
        for direction in ('up','down'):
            base=[]
            for q in (3,5,9,15,21):
                gg=.085*np.sqrt(q/150)
                candidate=[r for r in table if r['q']==q and r['N']==2 and r['n']==2000*q and r['loss']==loss
                           and r['direction']==direction and abs(r['g_over_quantum']/gg-1)<1e-10
                           and not r.get('detuning') and not r['time_censored'] and r['metric_rate'] is not None
                           and r['peak_probability']>=1e-12
                           and (loss not in ('lower','symmetric') or r['gamma_over_omega']==100.)]
                if candidate:base.append(candidate[0])
            if len(base)<3:continue
            qs=np.array([r['q'] for r in base],float)
            matrix=np.column_stack([np.ones(len(qs)),qs,qs*np.log(qs)])
            y=np.log([r['metric_rate'] for r in base]);coefficient=np.linalg.lstsq(matrix,y,rcond=None)[0]
            axes={}
            for axis in ('N','n','g_over_quantum'):
                possible=[f for f in fits if f['loss']==loss and f['direction']==direction and f['axis']==axis and f['fit_is_locally_powerlike']]
                if len(possible)>=2:
                    qs_axis=np.array([f['q'] for f in possible],float)
                    slopes=np.array([f['slope'] for f in possible])
                    linear=np.polyfit(qs_axis,slopes,1)
                    axes[axis]=dict(form='p(q)=a+bq',a=float(linear[1]),b=float(linear[0]),
                        q_range=[float(qs_axis.min()),float(qs_axis.max())],
                        max_abs_exponent_residual=float(abs(np.polyval(linear,qs_axis)-slopes).max()),
                        observed_axis_fits=possible)
            laws.append(dict(loss=loss,direction=direction,q_form='ln r=a+bq+cq ln q at N2,n2000q,g=.085sqrt(q/150)',
                coefficients=coefficient.tolist(),q_range=[float(qs.min()),float(qs.max())],
                max_abs_log_residual=float(abs(matrix@coefficient-y).max()),
                axis_exponent_laws=axes,joint_extension_assumption='log-separable one-at-a-time laws; combined foil coordinates were not sampled',
                observed_destination_peak_range=[min(r['peak_probability'] for r in base),max(r['peak_probability'] for r in base)]))
    return laws




def main():
    rows=[json.loads(line) for line in (HERE/'runs.jsonl').read_text().splitlines() if line.strip()]
    rows=list({r['key']:r for r in rows}.values())
    wanted={scaled.key(config) for config in scaled.configs('full')}
    rows=[r for r in rows if r.get('model_version')==scaled.VERSION and r['key'] in wanted]
    branches=branch_witness(rows)
    branch_by_key={scaled.key(w['config']):w for w in branches['cases']}
    for row in rows:
        if 'failure' in row:continue
        cfg=row['config'];witness=branch_by_key.get(row['key'])
        retained=witness is None or min(s['prior_cluster_min_cosine'] for s in witness['full_spectrum_cluster_track'])>=1e-4
        row['order_supported']=retained and not (cfg['q']>=15 and cfg['N']>=4)
        if not row['order_supported']:
            row['order_support_note']='order-resolved high-N/high-q region unsupported; raw bare-bin populations are not fitted as transfer rates'
    table,fits=family_tables(rows)
    scaled.save(HERE/'results.json',dict(runs=rows,selected_configurations=len(rows),expected_configurations=len(wanted),
        unresolved=[r for r in rows if 'failure' in r],
        unsupported_order_labels=[r['config'] for r in rows if 'failure' not in r and not r.get('order_supported',True)],
        source_frame='original Fock basis; no displacement-gauge occupation relabelled'))
    metric_laws=rate_extensions(table,fits)
    leading=[]
    for q in (3,5,9,15,21):
        for N in (1,2,4,8,16,64,256,1024):
            for loss in ('none','lower_limit'):
                coefficient,logc=precise_coefficient(q,N,loss=='lower_limit')
                leading.append(dict(q=q,N=N,loss=loss,coefficient=coefficient,ln_abs_coefficient=logc))
    control=[]
    for N in (1,2,4,16,64):
        none=next(r for r in leading if r['q']==5 and r['N']==N and r['loss']=='none')
        lower=next(r for r in leading if r['q']==5 and r['N']==N and r['loss']=='lower_limit')
        ratio=float(abs(D(lower['coefficient'])/D(none['coefficient'])))
        expected=(N-1)*(15*N-18)/18
        if abs(ratio-expected)>1e-10*max(1,expected):raise AssertionError((N,ratio,expected))
        control.append(dict(N=N,measured_ratio=ratio,source_closed_ratio=expected))
    leading_fits=[]
    for q in (3,5,9,15,21):
        for loss in ('none','lower_limit'):
            subset=[r for r in leading if r['q']==q and r['loss']==loss and r['ln_abs_coefficient'] is not None and r['N']>=16]
            leading_fits.append(dict(q=q,loss=loss,**line_fit([math.log(r['N']) for r in subset],[r['ln_abs_coefficient'] for r in subset])))
    qfit=[]
    for loss in ('none','lower_limit'):
        subset=[r for r in leading if r['N']>=16 and r['loss']==loss and r['ln_abs_coefficient'] is not None]
        qs=np.array([r['q'] for r in subset],float)
        ns=np.log(np.array([r['N'] for r in subset],float)/1024)
        matrix=np.column_stack([np.ones(len(qs)),qs,qs*np.log(qs),ns,qs*ns])
        values=np.array([r['ln_abs_coefficient'] for r in subset])
        coef=np.linalg.lstsq(matrix,values,rcond=None)[0]
        qfit.append(dict(loss=loss,form='a+bq+cq ln(q)+(d+eq)ln(N/1024)',coefficients=coef.tolist(),
            q_range=[3,21],N_range=[16,1024],N_reference=1024,
            max_abs_log_residual=float(abs(matrix@coef-values).max())))
    extrap=[]
    # This is a formal log extension, not a prediction. Physical odd g and the
    # required mode population are not supplied by the foil data.
    for label,gap,quantum,g_single in [('one_keV',1000.,.04162028883288571,None),
                                      ('real_DD_gap',23846530.,.04162028883288571,None),
                                      ('standalone_odd_gap',33980312.65027593,.037,6.055943351805947)]:
        q=int(math.ceil(gap/quantum));q+=1-q%2
        for N in (32768,10**6,10**11):
            n=.0006*N/quantum
            for fit in qfit:
                a,b,cc,dd,ee=fit['coefficients']
                coeff_log=a+b*q+cc*q*math.log(q)+(dd+ee*q)*(math.log(N)-math.log(1024))
                admissible=n>=q
                oscillator_log=.5*(gammaln(n+1)-gammaln(n-q+1)) if admissible else None
                glog=math.log(g_single/quantum) if g_single is not None else None
                log10= (coeff_log+oscillator_log+q*glog)/math.log(10) if admissible and glog is not None else None
                extrap.append(dict(label=label,gap_eV=gap,quantum_eV=quantum,q=q,N=N,n_from_archived_energy=n,loss=fit['loss'],
                    formal_ln_coefficient=coeff_log,ln_oscillator_factor=oscillator_log,conditional_log10_amplitude_over_quantum=log10,
                    formal_log10_amplitude_at_g_over_quantum_one=(coeff_log+oscillator_log)/math.log(10) if admissible else None,
                    coefficient_of_log10_g_over_quantum=q,
                    physical_g_over_quantum=None,standalone_g_over_quantum=g_single/quantum if g_single is not None else None,
                    occupied_source_has_q_quanta=admissible,fit_q_range=[3,21],fit_N_range=[16,1024],
                    extrapolated_population_support='none at these coordinates; fitted destination probabilities exist only inside the numerical family',
                    status='formal extension outside fitted range; no physical rate; n<q forbids this prepared-state transition' if not admissible else 'formal extension outside fitted range; physical coupling/participation unknown'))
    scaled.save(HERE/'fits.json',dict(transfer_metric_fits=fits,transfer_metric_gap_laws=metric_laws,precise_leading_coefficients=leading,leading_N_fits=leading_fits,
        leading_q_fits=qfit,q5_source_control=control,extrapolations=extrap,
        limits='Leading amplitudes and unconditioned metric rates are fitted separately. No sink loss or survival normalization is counted as successful transfer.',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    if '--reuse-proof' not in sys.argv:
        scaled.save(HERE/'consumer_proof.json',consumer_proof(rows))
    else:
        proof=json.loads((HERE/'consumer_proof.json').read_text())
        if proof['case'] not in [r['config'] for r in rows if 'failure' not in r]:
            raise ValueError('Saved consumer proof has no unchanged case')
    scaled.save(HERE/'branch_proof.json',branches)
    print(json.dumps({'cases':len(rows),'metric_fits':len(fits),'leading_fits':len(leading_fits)}))


if __name__=='__main__':main()
