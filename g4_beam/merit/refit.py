"""Refit three causes to actual energy/run-bound transported profiles.

No free per-material amplitude. Same fixed proton nuclear background, BH neutron,
85/15 current mixture/spread, common stitch and high-voltage focus as pass4.
Three matched extra-current powers0/1/2; no change to supplied profiles during fit.
"""
from pathlib import Path
import json
import sys
import time
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from g4_beam.forward import model as M
from g4_beam.forward.consumer import forward_yields
from g4_beam.forward.profile_thermal import transport_stopping
from g4_beam.merit.api import ENERGIES, write_csv
from g4_beam.merit.responses import Response
HERE=Path(__file__).resolve().parent


def fit(cause,power,profiles,old,responses,*,initial_params=None):
    start=time.monotonic();matlist=M.MATS
    initial=dict(old['params'])
    names=['s']+['Ue_'+m for m in matlist]
    lo=[.3]+[0.]*4;hi=[3.]+[15.]*4;lognames={'s'}
    if cause=='screeningrise':
        names+=['Ues_'+m for m in matlist]+['E0'];lo += [0.]*4+[.05];hi += [30.]*4+[3.]
    elif cause=='surface':
        names+=['Ues_'+m for m in matlist]+['d'];lo += [0.]*4+[1e-7];hi += [200.]*4+[100.];lognames.add('d')
        # The transported sheet already supplies surface nuclei; the screening
        # enhancement starts modestly instead of inheriting a uniform-density amplitude.
        for m in matlist:initial['Ues_'+m]=initial['Ue_'+m]+.3
    else:
        names+=['f_n','E_n'];lo += [1e-5,.5];hi += [5.,26.];lognames.add('f_n')
    if cause=='surface' and '--surface-alternative' in sys.argv:
        initial['d']=5.
        for m in matlist:initial['Ues_'+m]=2.5
    if initial_params is not None:initial.update(initial_params)
    encode=lambda vals:np.array([np.log(v) if n in lognames else v for n,v in zip(names,vals)])
    lower=encode(lo);upper=encode(hi)
    x0=encode([max(initial[n],lo[i]+1e-10) for i,n in enumerate(names)])
    data=M.load_data(str(HERE.parent/'forward/errbars.csv'))
    fraction=lambda e:M.collected_fraction(e,1,0)
    base_meas=dict(current_A=lambda e:M.I_cup(e)*1e-6/fraction(e),collected_fraction=fraction,
                   neutral_reference_current_A=M.I_INT*1e-6,extra_current_power=power)
    calls=0
    def predictions(d,exact=False):
        out={}
        for m in matlist:
            pars={}
            cand={'screeningrise':'sigmoid_tt','surface':'layer','neutrals':'fast_neutrals'}[cause]
            if cause=='screeningrise':pars=dict(Ue_s=d['Ues_'+m],E0=d['E0'])
            if cause=='surface':pars=dict(Ue_s=d['Ues_'+m],d_nm=d['d'])
            if cause=='neutrals':pars=dict(f_n=d['f_n'],E_n_lab_keV=d['E_n'],delta=1)
            kw=dict(Ue_keV=d['Ue_'+m],candidate=cand,parameters=pars,stopping_fn=transport_stopping)
            ne=ENERGIES[ENERGIES>=4.5]
            if exact:
                p=forward_yields(ENERGIES,profiles[m],measurement=dict(base_meas,channel='p'),**kw)
                n=forward_yields(ne,profiles[m],measurement=dict(base_meas,channel='n'),**kw)
            else:
                signal={}
                for ch in ('p','n'):
                    values=[]
                    for e in ENERGIES:
                        response=responses[m][float(e)]
                        yi=response.ions(d['Ue_'+m],d.get('Ues_'+m,d['Ue_'+m]),cause,d,ch)
                        yn=0.
                        if cause=='neutrals':
                            true_current=M.I_cup(e)*1e-6/fraction(e)
                            yn=d['f_n']*M.I_INT*1e-6/true_current*response.neutral(d['E_n'],d['Ue_'+m],ch)
                        values.append((yi+yn)/fraction(e)/(M.I_cup(e)*1e-6)**power)
                    signal[ch]=np.asarray(values)
                p={'reported_signal':signal['p']}
                n={'reported_signal':signal['n'][ENERGIES>=4.5]}
            stitch=n['reported_signal'][0]/p['reported_signal'][np.flatnonzero(ENERGIES==4.5)[0]]
            pred=np.empty(len(ENERGIES));low=ENERGIES<=4
            pred[low]=p['reported_signal'][low]*stitch*d['s']/n['reported_signal'][-1]
            pred[~low]=n['reported_signal']/n['reported_signal'][-1]
            out[m]=pred
        return out
    def unpack(x):return {n:float(np.exp(v) if n in lognames else v) for n,v in zip(names,x)}
    def residual(x):
        nonlocal calls
        calls+=1;pred=predictions(unpack(x));res=[]
        for m in matlist:
            e,y,s=data[m];v=np.interp(e,ENERGIES,pred[m])
            res.extend((np.log(v)-np.log(y))/s)
        if calls==1:print('initial',cause,power,'chi2',float(np.dot(res,res)),'warm_s',time.monotonic()-start,flush=True)
        return np.asarray(res)
    if '--pilot' in sys.argv:
        a=time.monotonic();residual(x0);warm=time.monotonic()-a
        a=time.monotonic();residual(x0);hot=time.monotonic()-a
        print(json.dumps(dict(cause=cause,power=power,warm_s=warm,cached_evaluation_s=hot)),flush=True)
        return None
    result=least_squares(residual,np.clip(x0,lower+1e-9,upper-1e-9),bounds=(lower,upper),x_scale='jac',
                         max_nfev=250,ftol=2e-6,xtol=2e-6,gtol=2e-6)
    pars=unpack(result.x);surrogate=predictions(pars);pred=predictions(pars,exact=True)
    exact_res=np.concatenate([(np.log(np.interp(data[m][0],ENERGIES,pred[m]))-np.log(data[m][1]))/data[m][2] for m in matlist])
    discrepancy=max(float(np.max(np.abs(pred[m]/surrogate[m]-1))) for m in matlist)
    row=dict(cause=cause,cand={'screeningrise':'sigmoid_tt','surface':'layer','neutrals':'neutral d1'}[cause],
             meas='M2 focus; explicit current reading',current_power=power,b_hi=1,b_lo=0,
             params=pars,chi2=float(exact_res@exact_res),surrogate_chi2=float(result.fun@result.fun),
             surrogate_max_relative_error=discrepancy,n=48,k=len(names),dof=48-len(names),
             pass4_best_chi2=old['chi2'],pass4_best_meas=old['meas'],
             success=bool(result.success),message=result.message,nfev=result.nfev,calls=calls,
             optimality=float(result.optimality),wall_s=time.monotonic()-start,
             predicted_R={m:pred[m].tolist() for m in matlist})
    print('FINISHED',cause,power,'chi2',row['chi2'],'success',row['success'],'wall_s',row['wall_s'],flush=True)
    return row


def main():
    records=json.loads((HERE/'profiles.json').read_text())['records']
    profiles={m:{float(r['case']['E_cm_keV']):(r['layers'],r['nads_m2'],r['state']) for r in records if r['case']['condition']==m} for m in M.MATS}
    source=json.loads((HERE.parent/'forward/fits.json').read_text())
    candidates={'screeningrise':('sigmoid_tt',),'surface':('layer',),'neutrals':('neutral d0','neutral d1')}
    results=[]
    responses={m:{float(e):Response(float(e),profiles[m][float(e)]) for e in ENERGIES} for m in M.MATS}
    if '--surface-alternative' in sys.argv:
        old=min((r for r in source if r['cand']=='layer'),key=lambda r:r['chi2'])
        alternative=fit('surface',0,profiles,old,responses)
        results=json.loads((HERE/'fits_transported.json').read_text())
        i=next(i for i,r in enumerate(results) if r['cause']=='surface' and r['current_power']==0)
        previous=results[i]
        selected=min((previous,alternative),key=lambda r:r['chi2'])
        selected['initial_condition_comparison']={'original_chi2':previous['chi2'],
                'finite_5nm_start_chi2':alternative['chi2'],'finite_5nm_start_success':alternative['success']}
        results[i]=selected
        (HERE/'fits_transported.json').write_text(json.dumps(results,indent=2))
        write_csv(HERE/'refit_summary.csv',[{k:r[k] for k in ('cause','current_power','chi2','pass4_best_chi2','pass4_best_meas','success','nfev','calls','optimality','wall_s')} for r in results])
        return
    for cause,cands in candidates.items():
        old=min((r for r in source if r['cand'] in cands),key=lambda r:r['chi2'])
        if cause=='neutrals':
            for m in M.MATS:
                for response in responses[m].values():response.build_neutral_table()
                print('compiled neutral response',m,flush=True)
        for power in range(3):
            row=fit(cause,power,profiles,old,responses)
            if row is None:return
            results.append(row)
            (HERE/'fits_transported.json').write_text(json.dumps(results,indent=2))
            write_csv(HERE/'refit_summary.csv',[{k:r[k] for k in ('cause','current_power','chi2','pass4_best_chi2','pass4_best_meas','success','nfev','calls','optimality','wall_s')} for r in results])
    (HERE/'refit_status.json').write_text(json.dumps({'all_completed':True,'all_solver_success':all(r['success'] for r in results),'readings':[0,1,2],'focus':{'b_hi':1,'b_lo':0}},indent=2))

if __name__=='__main__':main()
