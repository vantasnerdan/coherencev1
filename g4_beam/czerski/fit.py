"""Primary Czerski v3/EXFOR fits and shared-material comparisons.

No invented statistical errors: objectives are explicitly balanced log-residual
scores, NOT chi2/AIC. Only the Karahadian block has the existing chi2 convention.
Czerski source Table2 M is held fixed; it is not renamed a density measurement.
"""
from pathlib import Path
import csv
import json
import sys
import numpy as np
_trapz = getattr(np, 'trapezoid', None) or np.trapz  # NumPy 1.x and 2.x
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from g4_beam.forward import model as M
from g4_beam.forward import thermal
from g4_beam.forward import profile_thermal as PT
HERE=Path(__file__).resolve().parent
CANDIDATES=('none','layer','sigmoid_tt','fast_neutrals','thermal')
TABLE={r['target']:r for r in csv.DictReader(open(HERE/'data/table2.csv'))}
# Homogeneous densities derived from Table1, NOT a new measurement.
COMPOSITION={'Zr':(40,91.224,1.74,5.61),'Ti':(22,47.867,1.93,3.97),'Pd':(46,106.42,.054,12.02)}
EFF={'Zr':.0032,'Ti':.0032,'Pd':150e-6/(4*np.pi*.06**2)}
HOT={'Zr':.340,'Ti':.370,'Pd':.400}


class SourceTarget(M.Target):
    """Source electronic sqrt(E) stopping, absolute coefficient from Table2 M.

    M=.001/e * efficiency * nD/C_lab * 2sqrt(2)S0/BG * 1e-31.
    S0=57 keV barn is the source convention; retained BH variation is explicit.
    This calibration uses the source's high-energy asymptotic approximation,
    not a guessed ICRU Zr coefficient. Angular first moment is as in Target.
    """
    def __init__(self,mat):
        self.mat=mat
        z,m,x,rho=COMPOSITION[mat]
        self.x=x
        self.reaction='proton'
        self.n=rho*6.02214076e23/(m+2.0141*x) # host cm^-3
        self.E=np.geomspace(.1,27.,1000)
        self._lnE=np.log(self.E)
        nd=self.x*self.n*1e6
        c=(.001/thermal.QE)*EFF[mat]*nd*1e-31*2*np.sqrt(2)*57000/M.BG/float(TABLE[mat]['yield_normalization_fit'])
        self.dEdx=c*np.sqrt(self.E)*1e-6 # eV/nm
        self.S=self.dEdx/(self.n*1e-22)
        self.Sn_host=M.Sn_atom(self.E,z,m)
        gam=4*M.M_D*m/(M.M_D+m)**2
        inv_ltr=self.n*2*self.Sn_host*1e-15/(gam*self.E*1000)*1e-7
        de=np.diff(self.E)*1000
        w=1/self.dEdx
        self.path=np.r_[0,np.cumsum((w[1:]+w[:-1])*.5*de)]
        g=inv_ltr*w
        self.Lam=np.r_[0,np.cumsum((g[1:]+g[:-1])*.5*de)]
        q=np.exp(self.Lam)*w
        self.Q=np.r_[0,np.cumsum((q[1:]+q[:-1])*.5*de)]


TARGETS={m:SourceTarget(m) for m in TABLE}
POINTS=list(csv.DictReader(open(HERE/'data/figure_points.csv')))
EXFOR=list(csv.DictReader(open(HERE/'data/exfor_o1084003_ti.csv')))
KARA=M.load_data(ROOT/'g4_beam/forward/errbars.csv')
KT={m:M.Target(m) for m in M.MATS}


def specs(c,mat,blocks):
    u=float(TABLE[mat]['screening_energy_cold_eV'])/1000
    sp=[('Ue',u,0.,3.,False)]
    if c in ('layer','sigmoid_tt'):
        sp += [('Ue_s',3.,0.,200. if c=='layer' else 30.,False)]
        sp += [('d',.003,1e-7,100.,True)] if c=='layer' else [('E0',1.,.05,3.,False)]
    if c=='fast_neutrals' and 'karahadian' in blocks:
        sp += [('f_n',.2,1e-8,5.,True),('E_n',5.,.5,26.,False)]
    if 'karahadian' in blocks:sp += [('s',1.,.3,3.,True)]
    if 'czerski' in blocks and mat=='Zr':
        sp += [('neutral_2a_D2',3e-4,1e-7,.003,True),('neutral_7a_D1',3e-4,1e-7,.003,True)]
    return sp


def unpack(sp,v):return {s[0]:float(np.exp(a) if s[4] else a) for s,a in zip(sp,v)}


def physical_yfun(t,c,p,mat):
    phys={'none':'const','fast_neutrals':'const','thermal':'const'}.get(c,c)
    yf=M.yield_function(t,phys,p)
    if c!='thermal':return yf
    nd=t.x*t.n*1e6
    extra=float(thermal.probability(nd,Ue_hot_keV=HOT[mat],E_A_eV=thermal.BARRIERS[mat],exposure_m3_s=p['A']))
    return lambda e,ch:yf(e,ch)+(extra if ch=='p' else 0.)


def source_prediction(mat,c,p,points):
    t=TARGETS[mat];yf=physical_yfun(t,c,p,mat)
    e=np.array([float(r['energy_cm_keV']) for r in points])
    y=yf(2*e,'p')
    # Known neutral contamination is a beam-state measurement term in EVERY
    # candidate, never a new material candidate or a Karahadian cup correction.
    for i,r in enumerate(points):
        en=None;fn=0.
        if mat=='Zr' and r['figure_panel']=='2a' and r['beam_species']=='D2+' and e[i]<5:
            en=10.;fn=p.get('neutral_2a_D2',3e-4)
        elif mat=='Zr' and r['figure_panel']=='7a' and e[i]<10:
            en=20.;fn=p.get('neutral_7a_D1',3e-4)
        elif r['beam_species']=='D3+' and mat=='Zr':
            en=4.;fn=3e-4 # source estimate; reported explicitly, not silently subtracted
        if en is not None:y[i]+=fn*float(yf(np.array([en]),'p')[0])
    return y*EFF[mat]/thermal.QE*.001


def blocks_for(mat):
    # Fig7 reuses the high-energy D1+ data of Fig2. Keep its distinct
    # below-3keV CM neutral-dominated extension, not a second high-E likelihood.
    out={'czerski':[r for r in POINTS if r['observable']=='yield' and r['target']==mat
                    and not (r['figure_panel']=='7a' and float(r['energy_cm_keV'])>=3.)]}
    if mat=='Ti':out['exfor']=EXFOR
    if mat in KT:out['karahadian']=KARA[mat]
    return out


def predict(mat,c,p,block,data,kpower=0):
    if block=='czerski':
        pp=dict(p)
        if c=='thermal' and 'A' not in pp:pp['A']=p['A_czerski']
        return source_prediction(mat,c,pp,data)
    if block=='exfor':
        e=np.array([float(r['energy_cm_keV']) for r in data])
        pp=dict(p)
        if c=='thermal' and 'A' not in pp:pp['A']=p['A_exfor']
        yf=physical_yfun(TARGETS[mat],c,pp,mat)
        # No above-nominal neutral population is supplied by EXFOR; there is
        # no evidence to manufacture one. The normalization is 3.225keV CM.
        return yf(2*e,'p')/yf(np.array([6.45]),'p')[0]
    e=data[0]
    if c=='thermal':
        amplitude=p['A'] if 'A' in p else p['A_karahadian']
        pp=dict(p,thermal_per_charge=lambda ee:PT.charge_probability(mat,ee,amplitude))
        return M.relative_yield(KT[mat],e,'thermal',pp,{'k':kpower,'s':p['s']})
    meas={'k':kpower,'s':p['s']}
    if c=='fast_neutrals':meas.update(f_n=p['f_n'],E_n=p['E_n'],delta=1)
    return M.relative_yield(KT[mat],e,{'none':'const','fast_neutrals':'const'}.get(c,c),p,meas)


def observed(block,data):
    if block=='czerski':return np.array([float(r['value']) for r in data])
    if block=='exfor':return np.array([float(r['yield_relative_to_lab_6p45keV']) for r in data])
    return data[1]


def fit_one(mat,c,blocks,kpower=0):
    sp=specs(c,mat,blocks)
    lo=np.array([np.log(s[2]) if s[4] else s[2] for s in sp]);hi=np.array([np.log(s[3]) if s[4] else s[3] for s in sp])
    x0=np.array([np.log(s[1]) if s[4] else s[1] for s in sp])
    def residual(v):
        p=unpack(sp,v)
        return np.concatenate([(np.log(predict(mat,c,p,b,d,kpower))-np.log(observed(b,d)))/np.sqrt(len(observed(b,d))) for b,d in blocks.items()])
    best=None
    rng=np.random.default_rng(43)
    for j in range(3):
        start=x0 if j==0 else np.clip(x0+rng.normal(0,.1,len(x0))*(hi-lo),lo+1e-9,hi-1e-9)
        r=least_squares(residual,start,bounds=(lo,hi),max_nfev=700,ftol=2e-8,xtol=2e-8,gtol=2e-7,x_scale='jac')
        if best is None or r.cost<best.cost:best=r
    p=unpack(sp,best.x)
    per={}
    for b,d in blocks.items():
        lr=np.log(predict(mat,c,p,b,d,kpower)/observed(b,d))
        per[b]=dict(n=len(lr),log_sse=float(lr@lr),log_rmse=float(np.sqrt(np.mean(lr**2))),rms_factor=float(np.exp(np.sqrt(np.mean(lr**2)))))
        if b=='karahadian':per[b]['chi2_existing_error_convention']=float(np.sum((lr/d[2])**2))
    return dict(material=mat,candidate=c,blocks=list(blocks),karahadian_current_power=kpower,k=len(sp),balanced_log_score=float(2*best.cost),parameters=p,per_dataset=per,optimizer_success=bool(best.success),optimizer_message=best.message)


def proton_energy_lab(e_lab_keV):
    md,mp,mt=2.0141018,1.007825,3.0160493
    q=4032.66
    cs=np.cos(np.deg2rad(135))
    return ((np.sqrt(md*mp*np.asarray(e_lab_keV))*cs+np.sqrt(md*mp*np.asarray(e_lab_keV)*cs**2+(mt+mp)*((mt-md)*np.asarray(e_lab_keV)+mt*q)))/(mt+mp))**2


def spectral_comparison(best_source):
    result=[]
    for panel,mat,ypanel in [('2b','Zr','2a'),('4b','Pd','4a'),('7b','Zr','7a')]:
        r=best_source[mat];c=r['candidate'];p=r['parameters'];t=TARGETS[mat]
        yf=physical_yfun(t,c,p,mat)
        pts=[a for a in POINTS if a['figure_panel']==panel]
        predicted=[]
        for a in pts:
            e=float(a['energy_cm_keV']);el=2*e
            fake=dict(a,figure_panel=ypanel,observable='yield')
            total=source_prediction(mat,c,p,[fake])[0]*thermal.QE/(.001*EFF[mat])
            base=float(M.yield_function(t,{'none':'const','fast_neutrals':'const','thermal':'const'}.get(c,c),p)(np.array([el]),'p')[0])
            # Direct reaction-energy average along slowing down (not at E0).
            grid=t.E[t.E<=el]
            if len(grid)<2:direct_ep=float(proton_energy_lab(el))
            else:
                kernel=M.reaction_sigma(grid/2,p['Ue'],'p')/np.interp(grid,t.E,t.S)
                direct_ep=float(_trapz(proton_energy_lab(grid)*kernel,grid)/_trapz(kernel,grid))
            extra=max(total-base,0.)
            extra_ep=float(proton_energy_lab(0))
            if mat=='Zr' and a['beam_species']=='D2+':extra_ep=float(proton_energy_lab(10))
            if panel=='7b':extra_ep=float(proton_energy_lab(20))
            predicted.append((base*direct_ep+extra*extra_ep)/max(total,1e-300))
        pred=np.array(predicted);obs=np.array([float(a['value']) for a in pts]);ener=np.array([float(a['energy_cm_keV']) for a in pts])
        high=ener>=3
        # Calibrate detector channels only on high-energy line positions, then
        # keep the low-energy spectral distinction as a dependent prediction.
        gain,offset=np.polyfit(pred[high],obs[high],1)
        for a,pe,ob in zip(pts,pred,obs):
            result.append(dict(panel=panel,target=mat,species=a['beam_species'],E_cm_keV=float(a['energy_cm_keV']),observed_channel=ob,predicted_channel=float(gain*pe+offset),calibration_high_energy_only=True,candidate=c))
    return result


def fit_thermal_global(datasets, kpower=0, common_A=False):
    """One amplitude per DATASET, never per material; optionally one for all."""
    datasets=tuple(datasets)
    material_data={}
    if 'czerski' in datasets:
        for mat in TABLE:material_data.setdefault(mat,{})['czerski']=blocks_for(mat)['czerski']
    if 'exfor' in datasets:material_data.setdefault('Ti',{})['exfor']=EXFOR
    if 'karahadian' in datasets:
        for mat in M.MATS:material_data.setdefault(mat,{})['karahadian']=KARA[mat]
    sp=[('Ue_'+mat,float(TABLE[mat]['screening_energy_cold_eV'])/1000 if mat in TABLE else 1.,0.,3.,False) for mat in material_data]
    ampnames=['A'] if common_A or len(datasets)==1 else ['A_'+d for d in datasets]
    sp += [(name,1.87e-34,1e-50,1e-20,True) for name in ampnames]
    if 'karahadian' in datasets:sp += [('s',1.,.3,3.,True)]
    if 'czerski' in datasets:
        sp += [('neutral_2a_D2',3e-4,1e-7,.003,True),('neutral_7a_D1',3e-4,1e-7,.003,True)]
    lo=np.array([np.log(s[2]) if s[4] else s[2] for s in sp]);hi=np.array([np.log(s[3]) if s[4] else s[3] for s in sp])
    x0=np.array([np.log(s[1]) if s[4] else s[1] for s in sp])
    counts={d:sum(len(observed(d,blocks[d])) for blocks in material_data.values() if d in blocks) for d in datasets}
    def local(global_p,mat):
        return dict(global_p,Ue=global_p['Ue_'+mat])
    def residual(v):
        pp=unpack(sp,v);out=[]
        for mat,blocks in material_data.items():
            p=local(pp,mat)
            for block,data in blocks.items():
                out.extend(np.log(predict(mat,'thermal',p,block,data,kpower)/observed(block,data))/np.sqrt(counts[block]))
        return np.array(out)
    best=None;rng=np.random.default_rng(71)
    for j in range(3):
        start=x0 if j==0 else np.clip(x0+rng.normal(0,.06,len(x0))*(hi-lo),lo+1e-9,hi-1e-9)
        r=least_squares(residual,start,bounds=(lo,hi),max_nfev=900,ftol=2e-8,xtol=2e-8,gtol=2e-7,x_scale='jac')
        if best is None or r.cost<best.cost:best=r
    gp=unpack(sp,best.x)
    group='+'.join(datasets)+('_one_A' if common_A else '_separate')
    entries=[];aggregate={d:[] for d in datasets}
    for mat,blocks in material_data.items():
        p=local(gp,mat);per={}
        for block,data in blocks.items():
            lr=np.log(predict(mat,'thermal',p,block,data,kpower)/observed(block,data));aggregate[block].extend(lr)
            per[block]=dict(n=len(lr),log_sse=float(lr@lr),log_rmse=float(np.sqrt(np.mean(lr**2))),rms_factor=float(np.exp(np.sqrt(np.mean(lr**2)))))
            if block=='karahadian':per[block]['chi2_existing_error_convention']=float(np.sum((lr/data[2])**2))
        entries.append(dict(material=mat,candidate='thermal',blocks=list(blocks),global_scope=group,karahadian_current_power=kpower,k=len(sp),balanced_log_score=float(2*best.cost),parameters=p,per_dataset=per,optimizer_success=bool(best.success),optimizer_message=best.message))
    summary=dict(scope=group,current_power=kpower,parameters=gp,balanced_log_score=float(2*best.cost),optimizer_success=bool(best.success),per_dataset={d:dict(n=len(a),rms_factor=float(np.exp(np.sqrt(np.mean(np.asarray(a)**2))))) for d,a in aggregate.items()})
    print('THERMAL SHARED',summary,flush=True)
    return entries,summary


def main():
    results=[]
    membership=[]
    for i,r in enumerate(POINTS,2):
        selected=r['observable']=='yield' and not (r['figure_panel']=='7a' and float(r['energy_cm_keV'])>=3.)
        membership.append(dict(source_csv_line=i,figure_panel=r['figure_panel'],target=r['target'],species=r['beam_species'],used_in_yield_objective=selected,reason='measured yield' if selected else 'dependent observable or repeated high-energy Fig7 D1 data'))
    with open(HERE/'objective_membership.csv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(membership[0]));w.writeheader();w.writerows(membership)
    normalizations=[]
    for mat,t in TARGETS.items():
        normalizations.append(dict(material=mat,source_M=float(TABLE[mat]['yield_normalization_fit']),D_per_host=t.x,nD_m3=t.x*t.n*1e6,detector_fraction=EFF[mat],C_lab_keV_sqrt_per_m=float(t.dEdx[0]*1e6/np.sqrt(t.E[0])),meaning='Table1-derived homogeneous density; source-M-calibrated effective sqrt stopping; not independent density or efficiency determination'))
    with open(HERE/'normalization.csv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(normalizations[0]));w.writeheader();w.writerows(normalizations)
    for mat in TABLE:
        allblocks=blocks_for(mat)
        for c in CANDIDATES:
            if c=='thermal':continue
            source=fit_one(mat,c,{'czerski':allblocks['czerski']});results.append(source)
            print(mat,c,'czerski',source['per_dataset'],flush=True)
            if mat=='Ti':
                results.append(fit_one(mat,c,{'exfor':EXFOR}))
                results.append(fit_one(mat,c,{'czerski':allblocks['czerski'],'exfor':EXFOR}))
            if mat in KT:
                for k in (0,1):
                    r=fit_one(mat,c,allblocks,k);results.append(r)
                    print(mat,c,'shared M'+str(k),r['per_dataset'],flush=True)
    thermal_groups=[]
    for ds in [('czerski',),('exfor',),('czerski','exfor')]:
        entries,summary=fit_thermal_global(ds,common_A=len(ds)>1);results+=entries;thermal_groups.append(summary)
    for k in (0,1):
        for ds,shared in [(('czerski','karahadian'),False),(('czerski','karahadian'),True),(('czerski','exfor','karahadian'),True)]:
            entries,summary=fit_thermal_global(ds,k,shared);results+=entries;thermal_groups.append(summary)
    json.dump(thermal_groups,open(HERE/'thermal_shared.json','w'),indent=2)
    json.dump(results,open(HERE/'fits.json','w'),indent=2)
    rows=[]
    for r in results:
        for block,metrics in r['per_dataset'].items():
            rows.append(dict(material=r['material'],candidate=r['candidate'],fit_scope='+'.join(r['blocks']),global_scope=r.get('global_scope','per_material'),current_power=r['karahadian_current_power'],dataset=block,**metrics,parameters=json.dumps(r['parameters'],sort_keys=True)))
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with open(HERE/'cross_dataset.csv','w') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    render(results)


def source_fit(results,mat,c):
    return next(r for r in results if r['material']==mat and r['candidate']==c and r['blocks']==['czerski'] and r.get('global_scope','czerski_separate')=='czerski_separate')


def render(results):
    source={m:min((source_fit(results,m,c) for c in CANDIDATES),key=lambda r:r['per_dataset']['czerski']['log_sse']) for m in TABLE}
    spectral=[]
    for c in CANDIDATES:
        spectral+=spectral_comparison({m:source_fit(results,m,c) for m in TABLE})
    enhancement=[]
    for c in CANDIDATES:
        r=source_fit(results,'Ti',c)
        for a in POINTS:
            if a['figure_panel']!='3b':continue
            e=float(a['energy_cm_keV'])
            fake=dict(a,figure_panel='3a',beam_species='D3+')
            pred=source_prediction('Ti',c,r['parameters'],[fake])[0]
            denom=float(TABLE['Ti']['yield_normalization_fit'])*np.exp(-np.sqrt(986/(e+.025)))
            enhancement.append(dict(candidate=c,E_cm_keV=e,observed_enhancement=float(a['value']),predicted_enhancement=pred/denom,independent_yield_sample=False))
    with open(HERE/'enhancement_comparison.csv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(enhancement[0]));w.writeheader();w.writerows(enhancement)
    with open(HERE/'spectral_comparison.csv','w') as f:w=csv.DictWriter(f,fieldnames=list(spectral[0]));w.writeheader();w.writerows(spectral)
    fig,axs=plt.subplots(2,2,figsize=(11,8.5))
    colors=dict(zip(CANDIDATES,['#222222','#cd6600','#1b9954','#355ab5','#b32265']))
    replots=[]
    for ax,mat in zip(axs.flat,('Zr','Ti','Pd')):
        pts=blocks_for(mat)['czerski']
        for species,marker in [('D1+','s'),('D2+','o'),('D3+','^')]:
            a=[r for r in pts if r['beam_species']==species]
            ax.scatter([float(r['energy_cm_keV']) for r in a],[float(r['value']) for r in a],marker=marker,s=24,label=species,facecolors='none',edgecolors='#111111')
        for c in CANDIDATES:
            r=source_fit(results,mat,c)
            e=np.geomspace(.625,10,130);fake=[dict(target=mat,figure_panel='2a' if mat=='Zr' else '3a' if mat=='Ti' else '4a',beam_species='D3+' if ee<2.5 else 'D1+',energy_cm_keV=ee) for ee in e]
            y=source_prediction(mat,c,r['parameters'],fake)
            ax.plot(e,y,color=colors[c],label=c,lw=1.5)
            replots += [dict(panel=mat,candidate=c,E_cm_keV=float(ee),prediction=float(yy)) for ee,yy in zip(e,y)]
        if mat=='Zr':
            fixed=source_fit(results,mat,'thermal')['parameters']
            for panel,species,label in [('2a','D2+','thermal fit + known D2 neutrals'),('7a','D1+','thermal fit + known D1 neutrals')]:
                en=np.linspace(.65,3.5,90)
                state=[dict(target=mat,figure_panel=panel,beam_species=species,energy_cm_keV=ee) for ee in en]
                yp=source_prediction(mat,'thermal',fixed,state)
                ax.plot(en,yp,'--',lw=1.,label=label)
                replots += [dict(panel=mat,candidate=label,E_cm_keV=float(ee),prediction=float(yy)) for ee,yy in zip(en,yp)]
        ax.set(title=mat,yscale='log',xlabel='E_cm (keV)',ylabel='counts / (microampere kilosecond)');ax.grid(alpha=.2)
    ax=axs.flat[3]
    e=np.array([float(r['energy_cm_keV']) for r in EXFOR]);y=observed('exfor',EXFOR);ax.scatter(e,y,color='black',s=22,label='O1084003 Ti')
    for c in CANDIDATES:
        r=next(r for r in results if r['material']=='Ti' and r['candidate']==c and r['blocks']==['exfor'])
        yp=predict('Ti',c,r['parameters'],'exfor',EXFOR)
        ax.plot(e,yp,color=colors[c],label=c)
        replots += [dict(panel='EXFOR',candidate=c,E_cm_keV=float(ee),prediction=float(yy)) for ee,yy in zip(e,yp)]
    ax.set(title='EXFOR: reference 3.225 keV CM',yscale='log',xlabel='E_cm (keV)',ylabel='Y / Y(3.225 keV)');ax.grid(alpha=.2)
    for ax in axs.flat:ax.legend(fontsize=7,ncol=2)
    fig.suptitle('Source-normalized primary data; thermal A shared across all three materials\nZr includes separate source-identified neutral beam states; EXFOR reference is 3.225 keV CM')
    fig.tight_layout();fig.savefig(HERE/'primary_fits.png',dpi=170);plt.close(fig)
    with open(HERE/'replot.csv','w') as f:w=csv.DictWriter(f,fieldnames=list(replots[0]));w.writeheader();w.writerows(replots)
    fig,axs=plt.subplots(1,3,figsize=(12,3.8))
    for ax,panel in zip(axs,('2b','4b','7b')):
        a=[r for r in spectral if r['panel']==panel and r['candidate']=='thermal']
        for species in sorted(set(r['species'] for r in a)):
            b=[r for r in a if r['species']==species];e=[r['E_cm_keV'] for r in b]
            ax.scatter(e,[r['observed_channel'] for r in b],s=20,label=species)
        for c in ('none','thermal','layer'):
            b=sorted([r for r in spectral if r['panel']==panel and r['candidate']==c],key=lambda r:(r['species'],r['E_cm_keV']))
            for species in sorted(set(r['species'] for r in b)):
                bb=[r for r in b if r['species']==species]
                ax.plot([r['E_cm_keV'] for r in bb],[r['predicted_channel'] for r in bb],color=colors[c],lw=1,label=c if species==b[0]['species'] else None)
        ax.set(title='Fig '+panel,xlabel='E_cm (keV)',ylabel='proton peak channel');ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.suptitle('Dependent spectral check: high-energy-only channel calibration; no duplicated yield likelihood')
    fig.tight_layout();fig.savefig(HERE/'spectral_check.png',dpi=170)


if __name__=='__main__':
    render(json.load(open(HERE/'fits.json'))) if '--plot-only' in sys.argv else main()
