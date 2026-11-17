"""One joint shape/SI1 likelihood for the three existing physical causes.
Geometry and all nuclear quadrature come from merit.responses/forward.consumer.
"""
import time
import numpy as np
from scipy.optimize import least_squares
from g4_beam.forward import model as M
from g4_beam.merit.api import ENERGIES,evaluate_fit,_reported,paper_rows



def shape_residual(predictions,data):
    return np.concatenate([(np.log(np.interp(data[m][0],ENERGIES,predictions[m]))-np.log(data[m][1]))/data[m][2] for m in M.MATS])


def ratio_errors():
    table={(r['material'],float(r['E_cm'])):float(r['sigma_lnR']) for r in paper_rows()}
    return np.array([np.hypot(table['Pd',float(e)],table['PdD',float(e)]) for e in M.SI1_E])


def joint_fit(cause,targets,responses,baseline,data,*,initial_params=None):
    start=time.monotonic();names=['s']+['Ue_'+m for m in M.MATS]
    lo=[.3]+[0.]*4;hi=[3.]+[15.]*4;logs={'s'}
    if cause=='screeningrise':
        names+=['Ues_'+m for m in M.MATS]+['E0'];lo +=[0.]*4+[.05];hi +=[30.]*4+[3.]
    elif cause=='surface':
        names+=['Ues_'+m for m in M.MATS]+['d'];lo +=[0.]*4+[1e-7];hi +=[200.]*4+[100.];logs.add('d')
    elif cause=='neutrals':
        names+=['f_n','E_n'];lo +=[1e-5,.5];hi +=[5.,26.];logs.add('f_n')
        for material_responses in responses.values():
            for response in material_responses.values():
                if response.neutral_tables is None:response.build_neutral_table()
    else:raise ValueError('Existing three causes only')
    initial=baseline['params'] if initial_params is None else initial_params
    encode=lambda values:np.array([np.log(v) if n in logs else v for n,v in zip(names,values)])
    lower=encode(lo);upper=encode(hi);x0=encode([max(initial[n],lo[i]+1e-10) for i,n in enumerate(names)])
    observed=M.counts_per_C_table('PdD')/M.counts_per_C_table('Pd');errors=ratio_errors()
    positions=np.array([np.flatnonzero(ENERGIES==e)[0] for e in M.SI1_E])
    fraction=M.collected_fraction(ENERGIES,baseline['b_hi'],baseline['b_lo'])
    cup=M.I_cup(ENERGIES)*1e-6;current=cup/fraction;power=baseline['current_power'];low=ENERGIES<=4
    calls=0
    def unpack(x):return {n:float(np.exp(v) if n in logs else v) for n,v in zip(names,x)}
    def fast(d):
        physical={};pred={}
        for material in M.MATS:
            signals={}
            for ch in ('p','n'):
                yy=np.array([responses[material][float(e)].ions(d['Ue_'+material],d.get('Ues_'+material,d['Ue_'+material]),cause,d,ch) for e in ENERGIES])
                if cause=='neutrals':
                    yy += np.array([d['f_n']*M.I_INT*1e-6/current[i]*responses[material][float(e)].neutral(d['E_n'],d['Ue_'+material],ch) for i,e in enumerate(ENERGIES)])
                if ch=='p':physical[material]=yy
                signals[ch]=yy/fraction/cup**power
            stitch=signals['n'][ENERGIES==4.5][0]/signals['p'][ENERGIES==4.5][0]
            pred[material]=np.where(low,signals['p']*stitch*d['s'],signals['n'])/signals['n'][-1]
        return pred,physical['PdD'][positions]/physical['Pd'][positions]
    def residual(x):
        nonlocal calls
        calls+=1;pred,ratio=fast(unpack(x))
        return np.r_[shape_residual(pred,data),np.log(ratio/observed)/errors]
    result=least_squares(residual,np.clip(x0,lower+1e-9,upper-1e-9),bounds=(lower,upper),x_scale='jac',
                         max_nfev=250,ftol=2e-6,xtol=2e-6,gtol=2e-6)
    fit=dict(baseline,params=unpack(result.x),n=55,k=len(names),dof=55-len(names),success=bool(result.success),nfev=result.nfev,message=result.message)
    exact={m:_reported(targets[m],ENERGIES,cause,fit) for m in M.MATS};sr=shape_residual(exact,data)
    physical={m:evaluate_fit(targets[m],M.SI1_E,cause,fit)['physical_yield_per_C'] for m in ('Pd','PdD')}
    ratio=physical['PdD']/physical['Pd'];rr=np.log(ratio/observed)/errors
    shape=float(sr@sr);ratio_chi=float(rr@rr);joint=shape+ratio_chi;fit['chi2']=joint
    discrepancy=abs(joint-float(result.fun@result.fun))
    assert discrepancy<.05,('Compiled/direct joint likelihood discrepancy',cause,discrepancy)
    return dict(shape_chi2=shape,ratio_chi2=ratio_chi,joint_chi2=joint,
                mean_log_model_over_SI1=float(np.mean(np.log(ratio/observed))),geomean_model_PdD_over_Pd=float(np.exp(np.mean(np.log(ratio)))),
                min_model_PdD_over_Pd=float(np.min(ratio)),max_model_PdD_over_Pd=float(np.max(ratio)),
                ratio_all_energies_in_measured_band=bool(np.all((ratio>=min(observed))&(ratio<=max(observed)))),
                ratio_geomean_in_measured_band=bool(min(observed)<=np.exp(np.mean(np.log(ratio)))<=max(observed)),
                model_ratio_by_energy=ratio.tolist(),SI1_E_cm_keV=M.SI1_E.tolist(),fit=fit,success=bool(result.success),
                compiled_joint_chi2=float(result.fun@result.fun),nfev=result.nfev,calls=calls,wall_s=time.monotonic()-start)
