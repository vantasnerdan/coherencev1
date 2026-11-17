"""Replot delivered design/erosion CSVs; no transport or refitting."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
CAUSES=('screeningrise','surface','neutrals')
ENERGY=(.25,1.,2.5,6.5)
LABELS={1:'seal',2:'deuteride cap',3:'D-bearing oxide',4:'donor promoter',5:'BZY front/back',6:'Pd/oxide laminate'}


def main():
    rows=list(csv.DictReader((HERE/'all_designs.csv').open()))
    fig,axes=plt.subplots(1,3,figsize=(14,5),sharex=True)
    for ax,cause in zip(axes,CAUSES):
        for family in range(1,7):
            best=[]
            for e in ENERGY:
                r=max((r for r in rows if r['cause']==cause and int(r['family'])==family and float(r['E_cm_keV'])==e),key=lambda r:float(r['merit_PdD']))
                best.append(float(r['merit_PdD']))
            ax.semilogy(ENERGY,best,'o-',label=LABELS[family])
        ax.axhline(10,color='black',ls=':',lw=1);ax.set_xlabel('CM energy (keV)');ax.grid(alpha=.2)
        ax.set_title(cause+(' (canonical full-sheet fit chi2=262)' if cause=='surface' else ''))
    axes[0].set_ylabel('Best scenario end-run merit / transported PdD')
    axes[-1].legend(fontsize=8)
    fig.suptitle('Scenario envelopes, not confidence intervals or demonstrated ALD performance')
    fig.tight_layout();fig.savefig(HERE/'designs.png',dpi=160);plt.close(fig)
    erows=list(csv.DictReader((HERE/'early_erosion.csv').open()))
    fig,axes=plt.subplots(2,3,figsize=(14,8),sharex=True)
    for j,cause in enumerate(CAUSES):
        for i,e in enumerate((.25,2.5)):
            ax=axes[i,j]
            actual=[r for r in erows if r['cause']==cause and float(r['E_cm_keV'])==e and r.get('reference_age_s')]
            names=sorted({r['target'] for r in actual},key=lambda name:max(float(r['merit_PdD']) for r in actual if r['target']==name),reverse=True)[:3]
            for name in names:
                rr=sorted((r for r in actual if r['target']==name),key=lambda r:float(r['reference_age_s']))
                ax.semilogy([float(r['reference_age_s'])/3600 for r in rr],[float(r['merit_PdD']) for r in rr],'o-',label=name.replace('_coverage1.0',''))
            ax.axhline(10,color='black',ls=':',lw=1);ax.grid(alpha=.2)
            ax.set_title(cause+f'; E={e}keV');ax.legend(fontsize=6)
            if i==1:ax.set_xlabel('Beam time (h)')
            if j==0:ax.set_ylabel('Merit / same-age PdD')
    fig.tight_layout();fig.savefig(HERE/'seal_erosion.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5))
    for family in range(1,7):
        values=[max(float(r['refitted_zero_sheet_surface_merit']) for r in rows if int(r['family'])==family and float(r['E_cm_keV'])==e) for e in ENERGY]
        ax.semilogy(ENERGY,values,'o-',label=LABELS[family])
    ax.axhline(10,color='gray',ls=':');ax.set_xlabel('CM energy (keV)')
    ax.set_ylabel('End-run merit / refitted zero-sheet surface reference')
    ax.set_title('Shape-only reference chi2=44.206; fixed bulk, not joined transport')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(HERE/'secondary_reference.png',dpi=160);plt.close(fig)
    resolved=list(csv.DictReader((HERE/'best_resolved_integrals.csv').open()))
    fig,axes=plt.subplots(2,2,figsize=(12,9))
    scenarios=[('canonical_full_sheet',cause) for cause in CAUSES]+[('refitted_zero_sheet_surface','surface')]
    for ax,(scenario,cause) in zip(axes.flat,scenarios):
        selected=sorted((r for r in resolved if r['reference_scenario']==scenario and r['cause']==cause),key=lambda r:float(r['E_cm_keV']))
        for field,label,style in [('integrated_merit','charge-integrated','ko-'),('end_merit','end-state','b.--'),('peak_merit','resolved peak','r.:')]:
            ax.semilogy([float(r['E_cm_keV']) for r in selected],[float(r[field]) for r in selected],style,label=label)
        ax.axhline(10,color='gray',ls=':');ax.set_title(cause+' / '+scenario.replace('_',' '))
        ax.set_xlabel('CM energy (keV)');ax.set_ylabel('Design / same-age PdD');ax.legend();ax.grid(alpha=.2)
    fig.suptitle('Best charge-integrated candidate within the resolved set; conditional material endpoints')
    fig.tight_layout();fig.savefig(HERE/'resolved_integrals.png',dpi=160);plt.close(fig)
    profiles=list(csv.DictReader((HERE/'end_run_depth_profiles.csv').open()))
    best=list(csv.DictReader((HERE/'best_per_family.csv').open()))
    fig,axes=plt.subplots(2,3,figsize=(15,9))
    for ax,family in zip(axes.flat,range(1,7)):
        target=next(r['best_endpoint_target'] for r in best if int(r['family'])==family and r['cause']=='screeningrise' and float(r['E_cm_keV'])==6.5)
        for e in ENERGY:
            rr=[r for r in profiles if r['target']==target and float(r['E_cm_keV'])==e and r['column']=='covered' and float(r['zt_nm'])<=200]
            depth=[(float(r['zt_nm'])+float(r['zb_nm']))/2 for r in rr]
            concentration=np.array([float(r['D_per_metal']) for r in rr])
            ax.step(depth,np.ma.masked_less_equal(concentration,0),where='mid',label=f'{e:g}keV')
        ax.set_yscale('log');ax.set_xlim(0,200);ax.set_xlabel('Depth (nm)');ax.set_ylabel('Volumetric D / metal')
        ax.set_title(LABELS[family]+'; '+target[:42],fontsize=9);ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.suptitle('Reported end-run depth profiles; zero-D gaps masked, separate adsorbed sheets in CSV')
    fig.tight_layout();fig.savefig(HERE/'depth_profiles.png',dpi=160);plt.close(fig)
    spread=list(csv.DictReader((HERE/'cap_endpoint_spread.csv').open()))
    nominal=list(csv.DictReader((HERE/'nominal_cap_integrals.csv').open()))
    fig,axes=plt.subplots(2,2,figsize=(13,9))
    for ax,(scenario,cause) in zip(axes.flat,scenarios):
        rr=sorted((r for r in spread if r['material']=='TiD2_cap' and float(r['thickness_nm'])==20 and r['reference_scenario']==scenario and r['cause']==cause),key=lambda r:float(r['E_cm_keV']))
        for label,color in [('retained_slow','C0'),('retained_fast','C1'),('lost_slow','C2'),('lost_fast','C3')]:
            ax.semilogy([float(r['E_cm_keV']) for r in rr],[float(r[label+'_integrated_merit']) for r in rr],'o-',color=color,label=label.replace('_',' '))
        nn=sorted((r for r in nominal if r['reference_scenario']==scenario and r['cause']==cause),key=lambda r:float(r['E_cm_keV']))
        ax.semilogy([float(r['E_cm_keV']) for r in nn],[float(r['integrated_merit']) for r in nn],'k*--',label='retained, nominal sputter')
        anchored=[r for r in rr if r['retained_slow_TiD_counts_calibrated_integrated_merit']]
        ax.semilogy([float(r['E_cm_keV']) for r in anchored],[float(r['retained_slow_TiD_counts_calibrated_integrated_merit']) for r in anchored],'s:',color='C0',label='slow retained, SI1 TiD-counts proxy')
        ax.axhline(10,color='gray',ls=':');ax.set_title(cause+' / '+scenario.replace('_',' '),fontsize=9)
        ax.set_xlabel('CM energy (keV)');ax.set_ylabel('Charge-integrated merit');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('20nm TiD2: retention and sputtering decide reach; calibrated points are counts proxies')
    fig.tight_layout();fig.savefig(HERE/'cap_endpoints.png',dpi=160);plt.close(fig)
    thresholds=list(csv.DictReader((HERE/'cap_sputter_thresholds.csv').open()))
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for ax,(scenario,cause) in zip(axes.flat,scenarios):
        for reading,style in [('raw_model','o--'),('SI1_TiD_counts_calibrated','s-')]:
            rr=sorted((r for r in thresholds if r['reference_scenario']==scenario and r['cause']==cause and r['criterion']==reading and r.get('largest_verified_sputter_factor')),key=lambda r:float(r['E_cm_keV']))
            if rr:
                x=np.array([float(r['E_cm_keV']) for r in rr]);lo=np.array([float(r['largest_verified_sputter_factor']) for r in rr])
                hi=np.array([float(r['first_failing_factor']) if r['first_failing_factor'] else float(r['largest_verified_sputter_factor']) for r in rr])
                ax.errorbar(x,lo,yerr=[np.zeros(len(lo)),hi-lo],fmt=style,capsize=3,label=reading.replace('_',' '))
        ax.axhline(1,color='gray',ls=':');ax.set_yscale('log');ax.set_xlabel('CM energy (keV)')
        ax.set_ylabel('Sputter factor for integrated merit >=10');ax.set_title(cause+' / '+scenario.replace('_',' '),fontsize=9);ax.grid(alpha=.2)
        handles,_=ax.get_legend_handles_labels()
        if handles:ax.legend(fontsize=7)
        else:ax.text(.5,.5,'No >=10 factor found in scan',ha='center',transform=ax.transAxes)
    fig.suptitle('Finite rightmost-crossing brackets; SI1-calibrated headline only where measured')
    fig.tight_layout();fig.savefig(HERE/'sputter_thresholds.png',dpi=160);plt.close(fig)

if __name__=='__main__':main()
