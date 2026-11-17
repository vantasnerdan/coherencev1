"""Synthesize finite density grids, actual yield quotients and lab responses.
Replot mode does not fit or evolve a target. No Zenodo data are accessed.
"""
from pathlib import Path
import csv
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from g4_beam.merit.api import reference,bestfits,CAUSES,ENERGIES,FOUR,evaluate_fit,Target,write_csv
from g4_beam.merit.report import sealed_target
from g4_beam.forward import model as M
from g4_beam.forward.profile_thermal import transport_stopping
from g4_beam.causes.density import changed_target
from g4_beam.causes.sheet import occupied_target
from g4_beam.causes.references import (ZERO_SHEET_SURFACE,zero_sheet_fit,zero_sheet_reference,
                                      without_exposed_pd_sheet,pd_sheet_identity)
HERE=Path(__file__).resolve().parent


def physical_sealed(keep):
    source=sealed_target()
    def profile(e):
        layers,nads,state=source.profiles(e)
        return layers,nads if keep else 0.,dict(state,sheet_kept=keep)
    return Target('0.7D_Pd_'+('retained_sheet' if keep else 'no_sheet'),profile,stopping_fn=transport_stopping)


def main():
    one=json.loads((HERE/'density1d_fits.json').read_text());joint=json.loads((HERE/'density2d_fits.json').read_text())
    best=min(joint,key=lambda r:r['joint_chi2'])
    region=[r for r in joint if r['requested_region']]
    rows=[]
    # Recompute both physical numerator and denominator at each actual inferred
    # fit. Never divide an old merit by the density multiplier alone.
    for entry in [best]+[r for r in region if r is not best]:
        fit=entry['fit'];ref=changed_target('PdD',entry['PdD_multiplier'])
        den=evaluate_fit(ref,FOUR,'screeningrise',fit)['physical_yield_per_C']
        for keep in (True,False):
            target=physical_sealed(keep);num=evaluate_fit(target,FOUR,'screeningrise',fit)['physical_yield_per_C']
            unscaled_den=evaluate_fit(reference(),FOUR,'screeningrise',fit)['physical_yield_per_C']
            for i,e in enumerate(FOUR):rows.append(dict(case='joint_minimum' if entry is best else 'accepted_grid_region',sheet_retained=keep,E_cm_keV=float(e),
                        Pd_multiplier=entry['Pd_multiplier'],PdD_multiplier=entry['PdD_multiplier'],
                        physical_design_per_true_C=float(num[i]),physical_inferred_reference_per_true_C=float(den[i]),
                        merit_inferred_reference=float(num[i]/den[i]),merit_m1_profile_at_same_joint_fit=float(num[i]/unscaled_den[i]),
                        shape_chi2=entry['shape_chi2'],joint_chi2=entry['joint_chi2'],instantaneous_density_control=True))
    write_csv(HERE/'inferred_baseline_sealed.csv',rows)
    summary=dict(joint_minimum_multipliers={'Pd':best['Pd_multiplier'],'PdD':best['PdD_multiplier']},
                 requested_region_count=len(region),reference_policy='m1 baseline remains canonical; inferred-reference quotient is an explicit separate diagnostic',
                 joint_minimum_D_per_Pd={'Pd':best['Pd'],'PdD':best['PdD']},four_energy_no_sheet={},four_energy_retained_sheet={})
    for keep,name in [(False,'four_energy_no_sheet'),(True,'four_energy_retained_sheet')]:
        for e in FOUR:
            selected=[r for r in rows if r['sheet_retained']==keep and r['E_cm_keV']==e]
            center=next(r for r in selected if r['case']=='joint_minimum')
            allowed=[r['merit_inferred_reference'] for r in selected if (r['Pd_multiplier'],r['PdD_multiplier']) in {(x['Pd_multiplier'],x['PdD_multiplier']) for x in region}]
            summary[name][str(float(e))]=dict(joint_minimum_merit=center['merit_inferred_reference'],sampled_region_range=[min(allowed),max(allowed)] if allowed else None)
    (HERE/'inferred_baseline_summary.json').write_text(json.dumps(summary,indent=2))
    ratios=[]
    def append_ratios(label,cause,fit,targets):
        y={m:evaluate_fit(targets[m],M.SI1_E,cause,fit)['physical_yield_per_C'] for m in M.MATS}
        for a,b in [('PdD','Pd'),('TiD','Ti'),('TiD','PdD')]:
            observed=M.counts_per_C_table(a)/M.counts_per_C_table(b)
            for i,e in enumerate(M.SI1_E):ratios.append(dict(case=label,cause=cause,E_cm_keV=float(e),numerator=a,denominator=b,
                        model_ratio=float(y[a][i]/y[b][i]),measured_SI1_counts_ratio=float(observed[i]),model_over_measured=float(y[a][i]/y[b][i]/observed[i])))
    for cause in CAUSES:
        append_ratios('canonical_m1',cause,bestfits()[cause],{m:reference(m) for m in M.MATS})
        preferred=min((r for r in one if r['cause']==cause and r['kind']=='density'),key=lambda r:r['profile_chi2'])
        append_ratios('shared1d_profile_minimum',cause,preferred['fit'],{m:changed_target(m,preferred['multiplier']) for m in M.MATS})
    append_ratios('independent2d_joint_minimum','screeningrise',best['fit'],
                  {m:changed_target(m,best['Pd_multiplier'] if m=='Pd' else best['PdD_multiplier'] if m=='PdD' else 1.) for m in M.MATS})
    append_ratios(ZERO_SHEET_SURFACE,'surface',zero_sheet_fit(),
                  {m:without_exposed_pd_sheet(reference(m)) for m in M.MATS})
    sheet_rows=json.loads((HERE/'sheet_fits.json').read_text())
    for cause in CAUSES:
        optimum=min((r for r in sheet_rows if r['cause']==cause),key=lambda r:r['joint_chi2'])
        append_ratios('sheet_joint_minimum',cause,optimum['fit'],
                      {m:occupied_target(m,optimum['theta_Pd'] if m=='Pd' else optimum['theta_PdD'] if m=='PdD' else 1.) for m in M.MATS})
    write_zero_sheet_controls()
    verify_reference_projection()
    write_csv(HERE/'material_ratios.csv',ratios)
    plot()
    print(json.dumps(summary,indent=2),flush=True)


def write_zero_sheet_controls():
    fit=zero_sheet_fit();rows=[]
    denominator=evaluate_fit(zero_sheet_reference(),ENERGIES,'surface',fit)['physical_yield_per_C']
    for keep in (True,False):
        original=physical_sealed(keep);projected=without_exposed_pd_sheet(original)
        numerator=evaluate_fit(projected,ENERGIES,'surface',fit)['physical_yield_per_C']
        for i,e in enumerate(ENERGIES):
            rows.append(dict(target=original.name,reference_scenario=ZERO_SHEET_SURFACE,cause='surface',
                    E_cm_keV=float(e),shape_fit_chi2=fit['chi2'],original_density_control_sheet_retained=keep,
                    physical_design_per_true_C=float(numerator[i]),physical_reference_per_true_C=float(denominator[i]),
                    merit_zero_sheet_reference=float(numerator[i]/denominator[i]),
                    exposed_Pd_sheet_removed_on_both_sides=True,instantaneous_density_control=True,
                    **pd_sheet_identity(projected,float(e))))
    write_csv(HERE/'zero_sheet_density_control.csv',rows)


def verify_reference_projection():
    """One consumer witness for the new scenario's consequential area branch."""
    fit=zero_sheet_fit();pd=reference('PdD');other=reference('TiD')
    before=[pd.profiles(float(e))[1] for e in FOUR]
    zero_ref=evaluate_fit(zero_sheet_reference(),FOUR,'surface',fit)['physical_yield_per_C']
    projected=without_exposed_pd_sheet(pd)
    self_y=evaluate_fit(projected,FOUR,'surface',fit)['physical_yield_per_C']
    # Actual supplied profiles, used only as an algebraic parallel-column
    # witness, NOT as an additional ALD candidate or physical prediction.
    mixed=Target('algebraic_column_witness_not_a_design',pd.profiles,stopping_fn=transport_stopping,
                 columns=((.25,other.profiles),(.75,pd.profiles)))
    other_same_fit=Target('witness_nonPd_column',other.profiles,stopping_fn=transport_stopping)
    expected=.25*evaluate_fit(other_same_fit,FOUR,'surface',fit)['physical_yield_per_C']+.75*zero_ref
    transformed=without_exposed_pd_sheet(mixed)
    observed=evaluate_fit(transformed,FOUR,'surface',fit)['physical_yield_per_C']
    area_error=float(np.max(np.abs(observed/expected-1)))
    self_error=float(np.max(np.abs(self_y/zero_ref-1)))
    assert area_error<2e-12 and self_error<2e-12
    assert before==[pd.profiles(float(e))[1] for e in FOUR]
    for e in FOUR:
        assert pd_sheet_identity(transformed,float(e))['retains_beam_side_Pd_sheet_area_fraction']==0.
        assert transformed.columns[0][1](float(e))[1]==other.profiles(float(e))[1]
    one=json.loads((HERE/'density1d_fits.json').read_text())
    density=json.loads((HERE/'density2d_fits.json').read_text())
    sheet=json.loads((HERE/'sheet_fits.json').read_text())
    result=dict(zero_sheet_self_reference_max_absolute_error=self_error,
                zero_sheet_area_composition_max_relative_error=area_error,
                original_Pd_sheet_unchanged=True,nonPd_adsorbate_preserved=True,
                both_sides_use_same_explicit_surface_fit=True,
                density1d_rows=len(one),density2d_rows=len(density),sheet_rows=len(sheet),
                all_sheet_solver_success=all(r['success'] for r in sheet),
                max_sheet_compiled_direct_chi2_error=max(abs(r['joint_chi2']-r['compiled_joint_chi2']) for r in sheet),
                sampled_sheet_region_counts={c:sum(r['requested_region'] for r in sheet if r['cause']==c) for c in CAUSES},
                continuous_confidence_enclosure=False,joined_transport_triggered=any(r['requested_region'] for r in sheet))
    (HERE/'verification.json').write_text(json.dumps(result,indent=2))


def plot_sheet():
    import matplotlib.pyplot as plt
    rows=json.loads((HERE/'sheet_fits.json').read_text())
    fig,axes=plt.subplots(3,2,figsize=(13,14),sharex=True,sharey=True)
    ticks=sorted({r['theta_Pd'] for r in rows});xx,yy=np.meshgrid(ticks,ticks)
    for axisrow,cause in zip(axes,CAUSES):
        rr=[r for r in rows if r['cause']==cause];grid={(r['theta_Pd'],r['theta_PdD']):r for r in rr}
        best=min(rr,key=lambda r:r['joint_chi2']);accepted=[r for r in rr if r['requested_region']]
        for ax,key,label,cmap in [(axisrow[0],'shape_chi2','Shape chi2','viridis'),
                                  (axisrow[1],'mean_log_model_over_SI1','Mean log(model / SI1 PdD:Pd)','coolwarm')]:
            values=np.array([[grid[a,b][key] for a in ticks] for b in ticks])
            image=ax.pcolormesh(xx,yy,values,shading='nearest',cmap=cmap)
            fig.colorbar(image,ax=ax,label=label)
            ax.plot([0,1],[0,1],'k--',lw=.8,label='equal occupation')
            if accepted:ax.scatter([r['theta_Pd'] for r in accepted],[r['theta_PdD'] for r in accepted],
                                  marker='s',s=140,facecolors='none',edgecolors='white',lw=1.8,label='shape delta<=4 + all7 ratios in band')
            ax.scatter([best['theta_Pd']],[best['theta_PdD']],marker='*',s=170,color='red',edgecolors='black',label='joint minimum')
            ax.set_title(f'{cause}: {label}; accepted cells={len(accepted)}')
            ax.set_xticks(ticks);ax.set_yticks(ticks);ax.set_xlabel('Pd sheet occupation (ML)');ax.legend(fontsize=7)
        axisrow[0].set_ylabel('PdD sheet occupation (ML)')
    fig.suptitle('Fixed original m1 bulk: diagnostic sheet scan, not a joined transport state')
    fig.tight_layout();fig.savefig(HERE/'sheet_maps.png',dpi=160);plt.close(fig)


def plot():
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    one=json.loads((HERE/'density1d_fits.json').read_text());joint=json.loads((HERE/'density2d_fits.json').read_text())
    fig,axes=plt.subplots(1,3,figsize=(14,4.8))
    for ax,cause in zip(axes,CAUSES):
        rows=sorted((r for r in one if r['cause']==cause and r['kind']=='density'),key=lambda r:r['multiplier'])
        x=[r['multiplier'] for r in rows];y=[r['profile_chi2'] for r in rows];f=[r['frozen_chi2'] for r in rows]
        ax.semilogx(x,y,'o-',label='profiled cause parameters');ax.semilogx(x,f,'--',label='frozen m1 parameters')
        low=min(y);ax.axhline(low+1,color='gray',lw=.8,ls=':');ax.axhline(low+4,color='gray',lw=.8,ls='--')
        sheet=next(r for r in one if r['cause']==cause and r['kind']=='sheet0')
        ax.scatter([1],[sheet['profile_chi2']],marker='x',s=60,color='black',label='m1, sheet removed (profiled)')
        ax.set_title(cause);ax.set_xlabel('Shared Pd/PdD top10nm D multiplier');ax.grid(alpha=.2)
    axes[0].set_ylabel('48-point shape chi2');axes[-1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(HERE/'density1d.png',dpi=160);plt.close(fig)
    x=sorted({r['Pd_multiplier'] for r in joint});y=sorted({r['PdD_multiplier'] for r in joint})
    grid={(r['Pd_multiplier'],r['PdD_multiplier']):r for r in joint}
    xx,yy=np.meshgrid(x,y);shape=np.array([[grid[a,b]['shape_chi2'] for a in x] for b in y])
    ratio=np.array([[grid[a,b]['mean_log_model_over_SI1'] for a in x] for b in y])
    fig,axes=plt.subplots(1,2,figsize=(13,5.5),sharex=True,sharey=True)
    for ax,values,label,cmap in [(axes[0],shape,'Shape chi2 at joint optimum','viridis'),(axes[1],ratio,'Mean log(model / SI1 PdD:Pd)','coolwarm')]:
        image=ax.pcolormesh(xx,yy,values,shading='nearest',cmap=cmap);fig.colorbar(image,ax=ax,label=label)
        accepted=[r for r in joint if r['requested_region']]
        if accepted:ax.scatter([r['Pd_multiplier'] for r in accepted],[r['PdD_multiplier'] for r in accepted],s=120,marker='s',facecolors='none',edgecolors='white',linewidths=1.8,label='shape delta<=4 + all7 ratios in observed band')
        best=min(joint,key=lambda r:r['joint_chi2']);ax.scatter([best['Pd_multiplier']],[best['PdD_multiplier']],marker='*',s=180,color='red',edgecolors='black',label='joint minimum')
        ax.set_xscale('log');ax.set_yscale('log');ax.set_xlabel('Pd top10nm D multiplier');ax.set_title(label)
        ax.set_xticks(x,[f'{v:.2g}' for v in x]);ax.set_yticks(y,[f'{v:.2g}' for v in y]);ax.legend(fontsize=7,loc='best')
    axes[0].set_ylabel('PdD top10nm D multiplier')
    fig.suptitle('Finite7×7 grid; white squares are sampled cells, not a continuous enclosure')
    fig.tight_layout();fig.savefig(HERE/'density2d_maps.png',dpi=160);plt.close(fig)
    lab=list(csv.DictReader((HERE/'laboratory_ratios.csv').open()))
    fig,ax=plt.subplots(figsize=(10,4.5));indices=np.arange(len(lab));width=.23
    for i,cause in enumerate(CAUSES):ax.bar(indices+(i-1)*width,[float(r[cause+'_ratio']) for r in lab],width,label=cause)
    ax.set_yscale('log');ax.set_ylabel('Changed / baseline proton counts per incident charged C')
    ax.set_xticks(indices,[r['knob']+'='+r['value']+'\nE='+r['E_cm_keV']+'keV' for r in lab]);ax.legend();ax.grid(axis='y',alpha=.2)
    fig.tight_layout();fig.savefig(HERE/'laboratory_ratios.png',dpi=160);plt.close(fig)
    rows=list(csv.DictReader((HERE/'sealed_compare.csv').open()));inferred=list(csv.DictReader((HERE/'inferred_baseline_sealed.csv').open()))
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),sharex=True)
    for ax,cause in zip(axes,CAUSES):
        for name,label in [('sealed0.7_retained_sheet','0.7D/Pd + retained sheet'),('sealed0.7_no_sheet','0.7D/Pd, sheet removed')]:
            rr=[r for r in rows if r['target']==name and r['cause']==cause]
            ax.semilogy([float(r['E_cm_keV']) for r in rr],[float(r['merit_PdD']) for r in rr],label=label)
        if cause=='screeningrise':
            rr=[r for r in inferred if r['case']=='joint_minimum' and r['sheet_retained']=='False']
            ax.semilogy([float(r['E_cm_keV']) for r in rr],[float(r['merit_inferred_reference']) for r in rr],'k*--',label='sheet removed / inferred baseline')
        if cause=='surface':
            rr=[r for r in csv.DictReader((HERE/'zero_sheet_density_control.csv').open()) if r['original_density_control_sheet_retained']=='False']
            ax.semilogy([float(r['E_cm_keV']) for r in rr],[float(r['merit_zero_sheet_reference']) for r in rr],
                        'k*--',label='sheet removed / refitted zero-sheet')
        ax.axhline(10,color='gray',ls=':');ax.set_title(cause);ax.set_xlabel('CM energy (keV)');ax.grid(alpha=.2)
    axes[0].set_ylabel('Instantaneous density-control merit')
    for ax in axes:ax.legend(fontsize=7)
    fig.suptitle('Density controls, not an integrated ALD-device headline')
    fig.tight_layout();fig.savefig(HERE/'sealed_density_control.png',dpi=160);plt.close(fig)
    plot_sheet()

if __name__=='__main__':
    if '--plot-only' in sys.argv:plot()
    else:main()
