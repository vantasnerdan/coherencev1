"""Build the joined bulk table/page from actual isolated and tagged results."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import well

HERE, OUT = well.HERE, well.HERE/'results'
TABLES = HERE/'tables'
VARS = ('old','quadratic','quartic',*well.VARIANTS)

def dump(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def table(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def archive(variant,temperature,tag):
    root = OUT if variant in well.VARIANTS else well.PRIOR/'results'
    return np.load(root/f'D_{temperature}_{variant}_{tag}_distributions.npz')

def source_rows():
    old=json.loads((well.PRIOR/'results/summary_main.json').read_text())
    checks=json.loads((well.PRIOR/'results/summary_check.json').read_text())
    for variant in well.VARIANTS:
        old+=json.loads((OUT/f'summary_{variant}_main.json').read_text())
        checks+=json.loads((OUT/f'summary_{variant}_check.json').read_text())
    return old,checks

def potential():
    anchors=[]
    for j,axis in enumerate(well.RAYS):
        for p in well.ANCHORS[j]:
            row=dict(ray=axis,**p)
            row['predictions']={v:dict(energy_eV=float(well.radial(p['r_A'],v)[0][j]),
                outward_slope_eV_A=float(well.radial(p['r_A'],v)[1][j])) for v in VARS}
            anchors.append(row)
    signed=[]
    for p in well.WIDE['rows']:
        if p['material']!='PdD':continue
        u=np.array(p['ray_unit']);e,f=well.energy_force(p['displacement_A']*u,'wide_quadratic')
        signed.append(dict(case=p['case'],sample_kind=p['sample_kind'],parent_native_case=p.get('parent_native_case'),
            displacement_A=p['displacement_A'],energy_residual_eV=float(e-p['delta_energy_free_eV']),
            fixed_axis_force_residual_eV_A=float(f@u-p['native_force_parallel_eV_A'])))
    monotonic={}
    for axis,spline in zip(well.RAYS,well.SPLINES):
        candidates=list(spline.x)
        for j in range(len(spline.x)-1):
            a,b=spline.c[0,j],spline.c[1,j]
            if a!=0:
                t=-b/(3*a)
                if 0<t<spline.x[j+1]-spline.x[j]:candidates.append(spline.x[j]+t)
        monotonic[axis]=dict(min_slope_eV_A=float(spline(candidates,1).min()),
                            min_energy_eV=float(spline(candidates).min()))
        if monotonic[axis]['min_slope_eV_A']<0:
            raise ValueError('Anchor interpolation needs explicit monotonicity resolution')
    rng=np.random.default_rng(8025);x=rng.normal(size=(200,3))
    x*=rng.uniform(.001,1.2,len(x))[:,None]/np.linalg.norm(x,axis=1)[:,None]
    u=np.array([.71,.43,.56]);u/=np.linalg.norm(u)
    x=np.concatenate((x,np.array([r*u for r in (.03,.12,.24,.48,.80)])))
    gradients={}
    for v in well.VARIANTS:
        h=1e-6
        fd=np.column_stack([(well.energy_force(x+np.eye(3)[i]*h,v)[0]-well.energy_force(x-np.eye(3)[i]*h,v)[0])/(2*h) for i in range(3)])
        gradients[v]=float(np.max(abs(fd+well.energy_force(x,v)[1])))
    result=dict(anchors=anchors,signed_residuals=signed,radial_positivity=monotonic,
        max_abs_force_plus_gradient_eV_A=gradients,finite_difference_step_A=1e-6,
        K_eV_A2=well.K,F0_eV_A=[0.,0.,0.],angular_closure='unchanged accepted cubic nonnegative three-ray weights; not measured angular volume')
    dump(OUT/'potential.json',result)
    r=np.linspace(0,1.05,421);fig,axes=plt.subplots(2,3,figsize=(11,6))
    rows=[]
    for j,axis in enumerate(well.RAYS):
        for v in VARS:
            E,S=well.radial(r,v)
            for rr,e,s in zip(r,E[:,j],S[:,j]):rows.append(dict(ray=axis,variant=v,r_A=rr,energy_eV=e,outward_slope_eV_A=s))
            style=':' if v=='old' else ('--' if v in ('quartic','wide_quartic') else '-')
            color='black' if v=='old' else ('#dc8c27' if v in ('quadratic','quartic') else '#256aa0')
            for i,a in enumerate(axes[:,j]):a.plot(r,E[:,j] if i==0 else S[:,j],style,color=color,label=v)
        pts=well.ANCHORS[j]
        axes[0,j].scatter([p['r_A'] for p in pts],[p['energy_eV'] for p in pts],c='black',s=13)
        axes[1,j].scatter([p['r_A'] for p in pts],[p['slope_eV_A'] for p in pts],c='black',s=13)
        axes[0,j].set_title(f'PdD [{axis}]')
        for a in axes[:,j]:a.set_xlim(0,1.);a.axvline(.48,c='gray',lw=.6);a.axvline(.8,c='gray',lw=.6);a.set_xlabel('r (Å)');a.set_ylim(bottom=0)
    axes[0,0].set_ylabel('E (eV)');axes[1,0].set_ylabel('outward −F (eV/Å)');axes[0,0].legend(fontsize=6)
    fig.tight_layout();fig.savefig(HERE/'wells.png',dpi=160);plt.close(fig)
    table(TABLES/'energy_force.csv',rows)
    return result

def pdd():
    rows,checks=source_rows();common=[];numeric=[]
    fig,axes=plt.subplots(2,4,figsize=(14,6.5))
    for ti,T in enumerate((300,77)):
        zs={v:archive(v,T,'main') for v in VARS}
        for v in VARS:
            row=next(p for p in rows if p['variant']==v and p['temperature_K']==T)
            check=next(p for p in checks if p['variant']==v and p['temperature_K']==T)
            z,zf=zs[v],archive(v,T,'check')
            rr=z['position_radius'];tail=float(np.interp(.8,rr,z['position_tail']))
            fine_tail=float(np.interp(.8,zf['position_radius'],zf['position_tail']))
            resolved=tail>0 and fine_tail>0 and abs(tail-fine_tail)<.2*abs(fine_tail)
            num=dict(variant=v,temperature_K=T,delta_msd_A2=check['mean_r2_A2']-row['mean_r2_A2'],
                delta_kinetic_eV=check['mean_kinetic_eV']-row['mean_kinetic_eV'],
                delta_probability_gt_048=check['exterior_048']-row['exterior_048'],
                probability_gt_080_main=tail,probability_gt_080_check=fine_tail,
                probability_gt_080_empirically_resolved=bool(resolved))
            for kind in ('position','momentum'):
                a=z[f'{kind}_radius'];d=z[f'{kind}_density'];df=zf[f'{kind}_density']
                num[f'{kind}_density_integrated_L1_difference']=float(np.trapezoid(4*np.pi*a*a*abs(d-df),a))
            numeric.append(num)
            lines=np.array([np.interp(.48,rr,z['position_lines'][:,j]) for j in range(3)])
            common.append(dict(material='PdD',isotope='D',temperature_K=T,variant=v,
                msd_A2=row['mean_r2_A2'],kinetic_eV=row['mean_kinetic_eV'],probability_r_gt_048=row['exterior_048'],
                probability_r_gt_080=tail if resolved else None,
                support=dict(consumer='isolated 3D eigenstate/Fourier canonical',npar=row['npar'],quadrature_order=row['quadrature_order'],
                    radial_data_extent_A=.24 if v=='old' else (.48 if v in ('quadratic','quartic') else .8),
                    angular='same cubic closure; no angular data volume',line_density_anisotropy_at_048=float(lines.max()/lines.min()),
                    numerical_comparison=num)))
        for j,(kind,key,unit) in enumerate([('position','density','Å⁻³'),('position','measure','Å⁻¹'),('momentum','density','Å³'),('momentum','measure','Å')]):
            ax=axes[ti,j];a=zs['old'][f'{kind}_radius']
            ax.semilogy(a,zs['old'][f'{kind}_{key}'],':',c='black',label='old')
            for va,vb,color,label in [('quadratic','quartic','#dc8c27','.48 exterior pair'),('wide_quadratic','wide_quartic','#256aa0','.80 exterior pair')]:
                lo,hi=zs[va][f'{kind}_{key}'],zs[vb][f'{kind}_{key}']
                ax.plot(a,lo,c=color,label=label);ax.fill_between(a,np.minimum(lo,hi),np.maximum(lo,hi),color=color,alpha=.25)
            ax.set_title(f'D {T} K: {kind} {key}');ax.set_ylabel(unit)
            ax.set_xlabel('r (Å)' if kind=='position' else 'k=p/ℏ (Å⁻¹)')
            ax.set_xlim(0,1.05 if kind=='position' else 35);ax.set_ylim(bottom=max(zs['old'][f'{kind}_{key}'].max()*1e-6,1e-12))
            if kind=='position':ax.axvline(.48,c='gray',lw=.5);ax.axvline(.8,c='gray',lw=.5)
            gridrows=[]
            for i,x in enumerate(a):
                d=dict(coordinate=float(x))
                for v in VARS:
                    z=zs[v]
                    for k in ('density','measure','tail'):d[f'{v}_{k}']=float(z[f'{kind}_{k}'][i])
                    for ray,l in zip(well.RAYS,z[f'{kind}_lines'][i]):d[f'{v}_line_{ray}']=float(l)
                gridrows.append(d)
            table(TABLES/f'D_{T}_{kind}.csv',gridrows)
        signed=[]
        for i,s in enumerate(zs['old']['s111_A']):signed.append(dict(s_A=float(s),**{v:float(zs[v]['marginal111_A_inverse'][i]) for v in VARS}))
        table(TABLES/f'D_{T}_signed111.csv',signed)
    axes[0,0].legend(fontsize=7);fig.tight_layout();fig.savefig(HERE/'distributions.png',dpi=160);plt.close(fig)
    dump(OUT/'numerical.json',numeric)
    return common

def nih_figure():
    with np.load(HERE/'nih/distributions.npz') as z:
        fig,axes=plt.subplots(2,4,figsize=(14,6.5))
        for i,T in enumerate((300,77)):
            for j,(kind,coordinate,key,unit) in enumerate([
                ('position','position_r_A','position_spherical_density_A_inverse3','Å⁻³'),
                ('position','position_r_A','position_radial_measure_A_inverse','Å⁻¹'),
                ('momentum','momentum_k_A_inverse','momentum_spherical_density_A3','Å³'),
                ('momentum','momentum_k_A_inverse','momentum_radial_measure_A','Å')]):
                ax=axes[i,j];a=z[f'old_{T}K_{coordinate}']
                old=z[f'old_{T}K_{key}'];w=z[f'w_{T}K_{key}'];w2=z[f'w2_{T}K_{key}']
                ax.semilogy(a,old,':',c='black',label='old')
                ax.plot(a,w,c='#256aa0',label='111 correction w')
                ax.plot(a,w2,'--',c='#256aa0',label='111 correction w²')
                ax.fill_between(a,np.minimum(w,w2),np.maximum(w,w2),color='#256aa0',alpha=.25)
                ax.plot(a,z[f'w_stiff_exterior_{T}K_{key}'],'--',c='#ad5834',lw=.8,label='w, stiff exterior')
                ax.set_title(f'NiH H {T} K: {kind} '+('density' if 'density' in key else 'measure'))
                ax.set_ylabel(unit);ax.set_xlabel('r (Å)' if kind=='position' else 'k=p/ℏ (Å⁻¹)')
                ax.set_xlim(0,1. if kind=='position' else 40);ax.set_ylim(bottom=max(old.max()*1e-6,1e-12))
        axes[0,0].legend(fontsize=6);fig.tight_layout();fig.savefig(HERE/'nih-distributions.png',dpi=160);plt.close(fig)


def main():
    TABLES.mkdir(exist_ok=True)
    potential()
    isolated=pdd()
    coupled=json.loads((OUT/'coupled.json').read_text())
    rows=isolated+coupled['observable_rows']
    nih=json.loads((HERE/'nih/observable_table.json').read_text())
    for original in nih:
        row=dict(original)
        isolated_nih=row['variant'].startswith('isolated_')
        row['support']=dict(
            consumer='isolated 3D eigenstate/Fourier canonical' if isolated_nih else 'one-site finite-P reweight; centroid-virial kinetic',
            model_and_numeric='nih/support.json',
            evidence_key=(row['variant'].removeprefix('isolated_')+f'_{row["temperature_K"]}K') if isolated_nih else str(row['temperature_K']),
            angular='111 correction weighted by w or w²; 100/110 remain old-only, no measured angular volume')
        rows.append(row)
    nih_figure()
    dump(OUT/'observables.json',dict(rows=rows,units=dict(msd_A2='Å²',kinetic_eV='eV',probability='dimensionless'),
        separation='support.consumer distinguishes isolated and one-site finite-P rows; null is unsupported, not zero'))
    flat=[]
    for row in rows:
        flat.append({k:row.get(k) for k in ('material','isotope','temperature_K','variant','msd_A2','kinetic_eV','probability_r_gt_048','probability_r_gt_080')} | {'consumer':row.get('support',{}).get('consumer','NiH: see nih/README.md')})
    table(TABLES/'observables.csv',flat)
    table(TABLES/'coupled_characteristic_111.csv',coupled['characteristic_rows'])
    print(json.dumps(dict(isolated=isolated,rows=len(rows)),indent=2))

if __name__=='__main__':
    main()
