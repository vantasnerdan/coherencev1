"""Plot only the confirmed candidates and report physical absolute scale."""
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import model as M
import fits as F


def model_curve(result,material,energies):
    c=F.CANDS[result['cand']];d=result['params']
    p={'Ue':d['Ue_'+material]}
    for key,out in [('Ues_'+material,'Ue_s'),('E0','E0'),('d','d')]:
        if key in d:p[out]=d[key]
    if c['phys']=='thermal':
        p['thermal_per_charge']=lambda ee:F.PT.charge_probability(material,ee,d['A'])
    meas=dict(F.MEAS[result['meas']],s=d['s'])
    if 'neutral' in c:
        meas.update(f_n=d['f_n'],E_n=c.get('En_fixed',d.get('E_n')),delta=c['neutral'])
    return M.relative_yield(F.TARGETS[material],energies,c['phys'],p,meas)


def main():
    results=json.load(open(HERE/'fits.json'))
    candidates=[('no plateau',['none']),('fixed-energy fast neutrals',['neutral d0','neutral d1']),
                ('wide-bounds layer',['layer']),('sigmoid_tt',['sigmoid_tt']),('thermal spike',['thermal'])]
    best=[(label,min((r for r in results if r['cand'] in names),key=lambda r:r['aic'])) for label,names in candidates]
    e=np.linspace(.25,6.5,180)
    fig,axs=plt.subplots(2,2,figsize=(11,8))
    replot=[]
    for ax,mat in zip(axs.flat,M.MATS):
        ee,y,s=F.DATA[mat]
        for label,r in best:
            pred=model_curve(r,mat,e)
            ax.plot(e,pred,label=label+' ['+r['meas'].split()[0]+']',lw=1.6)
            replot.extend((mat,label,float(a),float(b)) for a,b in zip(e,pred))
        ax.errorbar(ee,y,yerr=y*s,fmt='o',color='black',ms=4,capsize=2,label='Fig2 digitized')
        ax.set(title=mat,yscale='log',ylim=(1e-4,1.5),xlabel='E_cm (keV)',ylabel='R(E) = signal / signal(6.5 keV)')
        ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Karahadian: confirmed mechanisms, explicit current/collection hypotheses\nBest AIC within each mechanism; yield-curve similarity is not mechanism identification')
    fig.tight_layout();fig.savefig(HERE/'fig_fits.png',dpi=170)
    with open(HERE/'replot.csv','w') as f:
        f.write('material,candidate,E_cm_keV,R\n')
        for row in replot:f.write('%s,%s,%.8g,%.8g\n'%row)
    bulk={r['meas']:r for r in results if r['cand']=='bulk 2.5-6.5'}
    with open(HERE/'abs_yield.csv','w') as f:
        f.write('material,normalization,Ue_keV,E_cm_keV,observed_counts_per_C,physical_reactions_per_C,required_efficiency_at_assumed_x,x_assumed\n')
        for mat in M.MATS:
            t=F.TARGETS[mat]
            for mn,r in bulk.items():
                u=r['params']['Ue_'+mat]
                for ec in (2.5,3.):
                    probability=M.beam_average(lambda a:t.Y_const(a,u),np.array([2*ec]))[0]
                    physical=probability/1.602176634e-19
                    count=M.counts_per_C_table(mat)[list(M.SI1_E).index(ec)]
                    f.write('%s,%s,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n'%(mat,mn,u,ec,count,physical,count/physical,t.x))
    for label,r in best:print(label,r['meas'],'chi2=',r['chi2'],'AIC=',r['aic'],'BIC=',r['bic'],r['params'])


if __name__=='__main__':main()
