"""Settled pass5 merit outputs, sensitivities and finite consumer identities."""
from pathlib import Path
import json
import sys
import time
import platform
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from g4_beam.merit.api import *
from g4_beam.forward import model as M
from g4_beam.forward.profile_thermal import transport_stopping
HERE=Path(__file__).resolve().parent


def scaled_profile(base,factor,depth=None,scale_ads=False):
    def profile(e):
        layers,nads,state=base(e);out=[]
        for row in layers:
            a=np.array(row,float)
            if depth is not None and a[0]<depth<a[1]:
                left=a.copy();left[1]=depth;left[3]*=factor
                right=a.copy();right[0]=depth;out.extend([left.tolist(),right.tolist()])
            else:
                if depth is None or a[0]<depth:a[3]*=factor
                out.append(a.tolist())
        return out,nads*(factor if scale_ads else 1),dict(state,diagnostic_D_factor=factor)
    return profile


def frozen_stopping(factor,depth=None):
    def stop(e,layer):
        a=list(layer)
        if depth is None or a[0]<depth:a[3]/=factor
        return transport_stopping(e,a)
    return stop


def sealed_target():
    def profile(e):
        layers,nads,state=reference().profiles(e);out=[]
        nh=6.8e28/(1+.16*.7)
        for row in layers:
            a=np.array(row,float)
            if a[0]<400e-9:a[2]=nh;a[3]=.7*nh;a[6]=nh
            out.append(a.tolist())
        return out,nads,dict(state,design='specified0.7D/Pd beam-side hydride; ideal no-stopping seal')
    return Target('ideal_sealed_0.7D_per_Pd',profile,stopping_fn=transport_stopping)


def main():
    start=time.monotonic();base=reference();tables=[];controls=[];sensitivity=[]
    # Primary self-reference is the exact same target object and denominator.
    for material in M.MATS:tables.extend(compare(reference(material)))
    uniform=[]
    for material in M.MATS:
        uniform.extend(compare(uniform_reference(material)))
    write_csv(HERE/'uniform_closure_compare.csv',uniform)
    # Global homogeneity identity is deliberately tested with matched stopping;
    # the physical D-stopping effect is a separate diagnostic target.
    for factor in (.1,3.,50.):
        profile=scaled_profile(base.profiles,factor,scale_ads=True)
        frozen=Target(f'global_D_x{factor}_matched_stopping',profile,stopping_fn=frozen_stopping(factor))
        physical=Target(f'global_D_x{factor}_physical_stopping',profile,stopping_fn=transport_stopping)
        local=Target(f'top5nm_D_x{factor}',scaled_profile(base.profiles,factor,depth=5e-9),stopping_fn=transport_stopping)
        for target in (frozen,physical,local):tables.extend(compare(target))
        for cause in CAUSES:
            a=evaluate(base,ENERGIES,cause)['physical_yield_per_C'];b=evaluate(frozen,ENERGIES,cause)['physical_yield_per_C']
            scale_error=float(np.max(np.abs(b/a/factor-1)))
            shape_error=float(np.max(np.abs((b/b[-1])/(a/a[-1])-1)))
            controls.append(dict(cause=cause,factor=factor,max_homogeneity_relative_error=scale_error,max_self_R_relative_error=shape_error))
            assert scale_error<3e-5 and shape_error<3e-5,(cause,factor,scale_error,shape_error)
            values=evaluate(local,FOUR,cause)['physical_yield_per_C']/evaluate(base,FOUR,cause)['physical_yield_per_C']
            for e,v in zip(FOUR,values):sensitivity.append(dict(cause=cause,E_cm_keV=float(e),density_factor=factor,near_surface_depth_nm=5.,yield_factor=float(v),log_elasticity=float(np.log(v)/np.log(factor)),stopping='physical D stopping; adsorbate unchanged'))
    for cause in CAUSES:
        for e in FOUR:
            initial=float(evaluate(base,[e],cause)['physical_yield_per_C'][0])
            def gain(shift):
                target=Target('screening_root',base.profiles,stopping_fn=transport_stopping,Ue_shift_keV=shift)
                return float(evaluate(target,[e],cause)['physical_yield_per_C'][0]/initial)
            extra=brentq(lambda x:np.log(gain(x))-np.log(10.),0.,15.,xtol=2e-7)
            sensitivity.append(dict(cause=cause,E_cm_keV=float(e),additional_Ue_for_10x_keV=extra,attained_yield_factor=gain(extra),screening_shift='uniform added shift of fitted cause potential; neutral reactions also change'))
    sealed=sealed_target();sealed_rows=compare(sealed)
    for row in sealed_rows:
        denominator=evaluate(uniform_reference(),[row['E_cm_keV']],row['cause'])['physical_yield_per_C'][0]
        row['merit_uniform_closure_same_refit']=row['physical_reactions_per_true_C']/denominator
    tables.extend(sealed_rows)
    write_csv(HERE/'compare.csv',tables);write_csv(HERE/'density_normalization.csv',controls)
    write_csv(HERE/'sensitivities.csv',sensitivity);write_csv(HERE/'sealed_reference_spread.csv',sealed_rows)
    # Actual SI1 ratios; no extrapolation or6.5keV experimental normalizer.
    ratios=[]
    for cause in CAUSES:
        y={m:evaluate(reference(m),M.SI1_E,cause)['physical_yield_per_C'] for m in M.MATS}
        for a,b in [('PdD','Pd'),('TiD','Ti'),('TiD','PdD')]:
            obs=M.counts_per_C_table(a)/M.counts_per_C_table(b)
            for i,e in enumerate(M.SI1_E):ratios.append(dict(cause=cause,E_cm_keV=float(e),numerator=a,denominator=b,model_same_E_ratio=float(y[a][i]/y[b][i]),measured_SI1_counts_ratio=float(obs[i]),model_over_observed=float(y[a][i]/y[b][i]/obs[i]),normalization='same-E shared cup collection and detector acceptance assumed'))
    write_csv(HERE/'transported_material_ratios.csv',ratios)
    closure=[]
    old=json.loads((HERE.parent/'forward/fits.json').read_text())
    for cause,names in {'screeningrise':('sigmoid_tt',),'surface':('layer',),'neutrals':('neutral d0','neutral d1')}.items():
        fit=min((r for r in old if r['cand'] in names),key=lambda r:r['chi2']);p=fit['params'];yield_old={};yield_new={}
        for m in M.MATS:
            t=M.Target(m);par={'Ue':p['Ue_'+m]}
            if cause!='neutrals':par['Ue_s']=p['Ues_'+m]
            if cause=='screeningrise':par['E0']=p['E0']
            if cause=='surface':par['d']=p['d']
            yf=M.yield_function(t,{'screeningrise':'sigmoid_tt','surface':'layer','neutrals':'const'}[cause],par)
            yi=M.beam_average(lambda en:yf(en,'p'),2*FOUR)
            cand={'screeningrise':'sigmoid_tt','surface':'layer','neutrals':'fast_neutrals'}[cause]
            pars={}
            if cause=='screeningrise':pars=dict(Ue_s=par['Ue_s'],E0=p['E0'])
            if cause=='surface':pars=dict(Ue_s=par['Ue_s'],d_nm=p['d'])
            measurement={}
            if cause=='neutrals':
                measurement=dict(current_A=lambda e:M.I_cup(e)*1e-6/M.collected_fraction(e,1,0),neutral_reference_current_A=M.I_INT*1e-6)
                yi+=p['f_n']*M.I_INT/(M.I_cup(FOUR)/M.collected_fraction(FOUR,1,0))*float(yf(np.array([p['E_n']]),'p')[0])
                pars=dict(f_n=p['f_n'],E_n_lab_keV=p['E_n'],delta=1)
            yield_old[m]=yi/QE
            yield_new[m]=forward_yields(FOUR,uniform_reference(m).profiles,Ue_keV=p['Ue_'+m],candidate=cand,parameters=pars,measurement=measurement)['physical_yield_per_C']
        for m in M.MATS:
            for i,e in enumerate(FOUR):closure.append(dict(cause=cause,material=m,E_cm_keV=float(e),pass4_uniform_ratio=float(yield_old[m][i]/yield_old['PdD'][i]),consumer_uniform_ratio=float(yield_new[m][i]/yield_new['PdD'][i]),relative_difference=float((yield_new[m][i]/yield_new['PdD'][i])/(yield_old[m][i]/yield_old['PdD'][i])-1),geometry='pass4 angular-first-moment layer versus consumer straight-CSDA; other causes matched'))
    write_csv(HERE/'pass4_ratio_closure.csv',closure)
    # One actual neutral target response and no detector contamination in merit.
    n0=evaluate(base,FOUR,'neutrals')['neutral_yield_per_C']
    n1=evaluate(sealed,FOUR,'neutrals')['neutral_yield_per_C']
    assert np.any(np.abs(n1/n0-1)>1e-3),'Changed target neutral response was frozen'
    for c in CAUSES:assert np.all(np.asarray(merit(base,FOUR,cause=c))==1.)
    import scipy
    summary=dict(runtime={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
                 wall_s=time.monotonic()-start,self_reference_exact=True,neutral_sealed_to_reference=(n1/n0).tolist(),
                 max_controlled_density_self_R_error=max(r['max_self_R_relative_error'] for r in controls),
                 source='853c37738caca76ba70cefd8f4df7350bc204937 plus local listed shared consumer refinement',
                 reference='actual energy/run-end transported PdD + adsorbate, refitted causes',
                 anchor='SI1 same-energy counts/C only;3.5–6.5keV null')
    (HERE/'verification.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
    plot()


def plot():
    import csv
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=list(csv.DictReader((HERE/'compare.csv').open()))
    fig,axes=plt.subplots(2,3,figsize=(14,8),sharex=True)
    for j,cause in enumerate(CAUSES):
        for name,label in [('UC_Davis_transported_PdD','PdD reference'),('ideal_sealed_0.7D_per_Pd','0.7D/Pd + retained adsorbate'),('global_D_x3.0_matched_stopping','global3x, fixed stopping'),('top5nm_D_x50.0','top5nm50x')]:
            rr=[r for r in rows if r['target']==name and r['cause']==cause]
            if not rr:continue
            e=[float(r['E_cm_keV']) for r in rr]
            axes[0,j].semilogy(e,[float(r['fixed_PdD_reference_curve']) for r in rr],label=label)
            axes[1,j].semilogy(e,[float(r['self_R']) for r in rr],label=label)
        axes[0,j].set_title(cause)
        axes[1,j].set_xlabel('CM energy (keV)')
        for ax in axes[:,j]:ax.grid(alpha=.2)
    axes[0,0].set_ylabel('Physical yield / fixed PdD6.5 yield')
    axes[1,0].set_ylabel('Physical self-R=Y(E)/Y6.5')
    axes[0,2].legend(fontsize=8);fig.tight_layout();fig.savefig(HERE/'merit.png',dpi=160);plt.close(fig)
    # Actual SI1 markers need an explicit model bridge to the fixed modeled
    # PdD6.5 scale. No experimental high-energy absolute anchors are invented.
    fig,axes=plt.subplots(1,3,figsize=(15,5),sharex=True)
    colors={'Pd':'tab:blue','PdD':'tab:orange','Ti':'tab:green','TiD':'tab:purple'}
    for ax,cause in zip(axes,CAUSES):
        base=[r for r in rows if r['target']=='UC_Davis_transported_PdD' and r['cause']==cause]
        for material in M.MATS:
            model=[r for r in rows if r['target']=='UC_Davis_transported_'+material and r['cause']==cause]
            ax.semilogy([float(r['E_cm_keV']) for r in model],[float(r['fixed_PdD_reference_curve']) for r in model],
                        color=colors[material],lw=1,label=material+' model')
            points=[r for r in base if r['experiment_'+material+'_fixed_PdD_anchor']]
            ax.scatter([float(r['E_cm_keV']) for r in points],[float(r['experiment_'+material+'_fixed_PdD_anchor']) for r in points],
                       color=colors[material],marker='x',s=30,label=material+' SI1 bridge')
        design=[r for r in rows if r['target']=='ideal_sealed_0.7D_per_Pd' and r['cause']==cause]
        ax.semilogy([float(r['E_cm_keV']) for r in design],[float(r['fixed_PdD_reference_curve']) for r in design],
                    color='black',lw=2,label='0.7D/Pd + retained adsorbate')
        ax.set_title(cause+(' (poor fit)' if cause=='surface' else ''))
        ax.set_xlabel('CM energy (keV)');ax.grid(alpha=.2)
    axes[0].set_ylabel('Physical yield / fixed modeled PdD6.5')
    handles,labels=axes[-1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=5,fontsize=8)
    fig.suptitle('SI1 crosses = measured same-energy ratios × modeled PdD(E)/PdD6.5; no3.5–6.5keV absolute data')
    fig.tight_layout(rect=(0,.15,1,.94));fig.savefig(HERE/'fixed_reference_comparison.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(3,3,figsize=(13,11),sharex=True)
    for i,cause in enumerate(CAUSES):
        base=[r for r in rows if r['target']=='UC_Davis_transported_PdD' and r['cause']==cause]
        design=[r for r in rows if r['target']=='ideal_sealed_0.7D_per_Pd' and r['cause']==cause]
        for power in range(3):
            ax=axes[i,power]
            for rr,label in ((base,'refitted PdD'),(design,'0.7D/Pd + retained adsorbate')):
                ee=np.array([float(r['E_cm_keV']) for r in rr])
                ax.semilogy(ee,[float(r[f'self_R_current_power_{power}']) for r in rr],label=label)
                u=float(rr[0][f'Ue_effective_{power}_keV']);high=ee>=2.5
                ax.semilogy(ee[high],M.sigma(ee[high],u)/M.sigma(6.5,u),':',alpha=.6)
            observed=[r for r in base if r['experiment_PdD_self_R']]
            ee=np.array([float(r['E_cm_keV']) for r in observed])
            rr=np.array([float(r['experiment_PdD_self_R']) for r in observed])*(M.I_cup(6.5)/M.I_cup(ee))**power
            ax.scatter(ee,rr,color='black',s=12,label='actual PdD Fig2')
            ax.set_title(cause+f'; extra current power{power}');ax.grid(alpha=.2)
            if i==2:ax.set_xlabel('CM energy (keV)')
            if power==0:ax.set_ylabel('Stitched reported self-R')
    axes[0,2].legend(fontsize=8);fig.tight_layout();fig.savefig(HERE/'paper_style_comparison.png',dpi=160);plt.close(fig)

if __name__=='__main__':
    if '--plot-only' in sys.argv:plot()
    else:main()
