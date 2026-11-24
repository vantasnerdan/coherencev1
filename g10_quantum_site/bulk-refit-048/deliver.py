"""Assemble the scoped old/new report from actual solver and coupled outputs."""
from pathlib import Path
import csv
import hashlib
import json
import numpy as np
from scipy.constants import physical_constants
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from well import HERE,OLD,WIDE,ANCHORS,RAYS,VECTORS,K,K4,R0,RMAX,VARIANTS,energy_force,radial,SPLINES

OUT=HERE/'results';TABLES=HERE/'tables'


def dump(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n')


def table(path,headers,arrays):
    with path.open('w') as f:
        w=csv.writer(f);w.writerow(headers);w.writerows(zip(*arrays))


def potential_evidence():
    evidence=[]
    for axis,vec,points in zip(RAYS,VECTORS,ANCHORS):
        for point in points:
            r=point['r_A'];e,f=energy_force(r*vec,'quadratic')
            olde,oldf=energy_force(r*vec,'old')
            evidence.append(dict(ray=axis,**point,fitted_energy_eV=float(e),fitted_force_parallel_eV_A=float(f@vec),old_energy_error_eV=float(olde-point['energy_eV']),old_force_error_eV_A=float(oldf@vec+point['slope_eV_A'])))
    signed=[]
    for row in WIDE['rows']:
        if row['material']!='PdD':continue
        vec=np.array(row['ray_unit']);e,f=energy_force(row['displacement_A']*vec,'quadratic')
        signed.append(dict(case=row['case'],ray=row['direction'],displacement_A=row['displacement_A'],sample_kind=row['sample_kind'],parent_native_case=row.get('parent_native_case'),energy_residual_eV=float(e-row['delta_energy_free_eV']),force_parallel_residual_eV_A=float(f@vec-row['native_force_parallel_eV_A'])))
    rng=np.random.default_rng(48025);xyz=rng.normal(size=(500,3));xyz*=rng.uniform(.001,1.2,500)[:,None]/np.linalg.norm(xyz,axis=1)[:,None]
    # Include nonspecial directions at every splice, not just symmetry rays.
    u=np.array([.71,.43,.56]);u/=np.linalg.norm(u)
    xyz=np.concatenate([xyz,np.array([r*u for r in (.03,.12,.24,.48)])])
    grad={};splices={}
    for variant in VARIANTS:
        h=1e-6;E,F=energy_force(xyz,variant)
        finite=np.column_stack([(energy_force(xyz+np.eye(3)[j]*h,variant)[0]-energy_force(xyz-np.eye(3)[j]*h,variant)[0])/(2*h) for j in range(3)])
        grad[variant]=dict(max_abs_force_plus_energy_gradient_eV_A=float(np.max(abs(F+finite))),step_A=h)
        splices[variant]=[]
        for r in (.03,.12,.24,.48):
            eps=1e-8;em,fm=energy_force((r-eps)*u,variant);ep,fp=energy_force((r+eps)*u,variant)
            splices[variant].append(dict(r_A=r,energy_two_sided_difference_eV=float(abs(ep-em)),force_two_sided_difference_eV_A=float(np.max(abs(fp-fm))),separation_A=2*eps))
    monotone={}
    for j,sp in enumerate(SPLINES):
        candidates=list(sp.x)
        for k in range(len(sp.x)-1):
            a,b=sp.c[0,k],sp.c[1,k]
            if abs(a)>0:
                t=-b/(3*a)
                if 0<t<sp.x[k+1]-sp.x[k]:candidates.append(sp.x[k]+t)
        monotone[RAYS[j]]=dict(min_slope_eV_A=float(sp(np.array(candidates),1).min()),min_energy_eV=float(sp(np.array(candidates)).min()))
    result=dict(anchors=evidence,signed_native_and_derived_residuals=signed,gradient=grad,splices=splices,radial_segment_positivity=monotone,near_K_eV_A2=K,origin='exact old cubic quartic polynomial for r <= .03 A; common harmonic K, no p/r6 angular residue at origin',closure='w100=1-4q+9p, w110=4q-36p, w111=27p; nonnegative by Schur/AM-GM, sum=1; C1 radial joins',support='native rays only; neither complete ball nor isotope-specific electronic surface')
    dump(OUT/'potential_evidence.json',result)
    table(TABLES/'energy_force_anchors.csv',['ray','r_A','energy_eV','outward_slope_eV_A','old_energy_error_eV','old_force_error_eV_A'],[[v[key] for v in evidence] for key in ('ray','r_A','energy_eV','slope_eV_A','old_energy_error_eV','old_force_error_eV_A')])
    r=np.linspace(0,.95,381)
    arrays=[r];headers=['r_A']
    fig,ax=plt.subplots(1,2,figsize=(10,3.6))
    for j,ray in enumerate(RAYS):
        for var,ls in [('old',':'),('quadratic','-'),('quartic','--')]:
            E,F=energy_force(r[:,None]*VECTORS[j],var)
            arrays.extend([E,F@VECTORS[j]]);headers.extend([f'{ray}_{var}_energy_eV',f'{ray}_{var}_force_eV_A'])
            ax[0].plot(r,E,ls,label=f'{ray} {var}');ax[1].plot(r,-F@VECTORS[j],ls)
        ax[0].scatter([p['r_A'] for p in ANCHORS[j]],[p['energy_eV'] for p in ANCHORS[j]],s=18)
        ax[1].scatter([p['r_A'] for p in ANCHORS[j]],[p['slope_eV_A'] for p in ANCHORS[j]],s=18)
    for a in ax:a.axvline(.48,color='k',lw=.6);a.set_xlabel('r (Å)');a.set_xlim(0,.85)
    ax[0].set_ylim(0,1.);ax[1].set_ylim(0,5);ax[0].set_ylabel('E (eV)');ax[1].set_ylabel('−F · outward ray (eV/Å)');ax[0].legend(fontsize=6,ncol=2)
    fig.tight_layout();fig.savefig(HERE/'wells.png',dpi=170);plt.close(fig)
    table(TABLES/'energy_force_curves.csv',headers,arrays)
    return result


def distribution_tables(rows):
    for iso in ('D','H','T'):
        for T in (300,77):
            tag='main' if iso=='D' else 'isotopes'
            available={v:OUT/f'{iso}_{T}_{v}_{tag}_distributions.npz' for v in VARIANTS}
            if not all(p.exists() for p in available.values()):continue
            z={v:np.load(p) for v,p in available.items()}
            for kind,units in [('position','A'),('momentum','A_inverse')]:
                r=z['old'][f'{kind}_radius'];arrays=[r];heads=[f'radius_{units}']
                for var in VARIANTS:
                    for key in ('density','measure','tail'):
                        arrays.append(z[var][f'{kind}_{key}']);heads.append(f'{var}_{key}')
                    for j,ray in enumerate(RAYS):arrays.append(z[var][f'{kind}_lines'][:,j]);heads.append(f'{var}_line{ray}_density')
                for key in ('density','measure','tail'):
                    stack=np.array([z[v][f'{kind}_{key}'] for v in ('quadratic','quartic')])
                    arrays.extend([stack.min(axis=0),stack.max(axis=0)]);heads.extend([f'new_continuation_{key}_min',f'new_continuation_{key}_max'])
                table(TABLES/f'{iso}_{T}_{kind}.csv',heads,arrays)
            if 's111_A' in z['old']:
                table(TABLES/f'{iso}_{T}_signed111.csv',['s_A','old_P_A_inverse','quadratic_P_A_inverse','quartic_P_A_inverse'],[z['old']['s111_A']]+[z[v]['marginal111_A_inverse'] for v in VARIANTS])


def archived_sampling_control(T,z):
    from scipy.interpolate import RegularGridInterpolator
    from solve import density,load_spectrum,kinetic
    s=load_spectrum(OUT/'D_old_main_spectrum.npz')
    positive=z['rho_positive_octant_A_minus3'][:51,:51,:51]
    indices=abs(np.arange(-50,51))
    rho=positive[np.ix_(indices,indices,indices)]
    x=np.linspace(-1.,1.,101)
    kp=np.linspace(0,42,61);nk,_=density(s,T,kp,True)
    ki=abs(np.arange(-60,61));nk=nk[np.ix_(ki,ki,ki)]
    mu,muw=np.polynomial.legendre.leggauss(16);phi=np.arange(32)*np.pi/16
    u=np.stack(np.broadcast_arrays(np.sqrt(1-mu[:,None]**2)*np.cos(phi),np.sqrt(1-mu[:,None]**2)*np.sin(phi),mu[:,None]+phi*0),axis=-1).reshape(-1,3)
    w=np.repeat(muw,32)/64
    result={}
    for kind,grid,values,r in [('position',x,rho,np.linspace(0,1,181)),('momentum',np.linspace(-42,42,121),nk,np.linspace(0,42,211))]:
        interp=RegularGridInterpolator((grid,)*3,values,bounds_error=False,fill_value=0.)
        d=interp((r[:,None,None]*u[None,:,:]).reshape(-1,3)).reshape(len(r),-1)@w
        result[kind]=float(np.trapezoid(4*np.pi*r**4*d,r))
    return dict(recreated_archived_reporting_r2_A2=result['position'],recreated_archived_reporting_kinetic_eV=result['momentum']*kinetic(s['mass']))


def numerical_evidence(rows):
    checks=json.loads((OUT/'summary_check.json').read_text())
    archive=json.loads((HERE/'inputs/archived-exact-summary.json').read_text())
    result=[]
    for row in rows:
        if row['isotope']!='D':continue
        T,var=row['temperature_K'],row['variant'];fine=next(v for v in checks if v['temperature_K']==T and v['variant']==var)
        z=np.load(OUT/f'D_{T}_{var}_main_distributions.npz');zf=np.load(OUT/f'D_{T}_{var}_check_distributions.npz')
        item=dict(temperature_K=T,variant=var,main_npar=row['npar'],check_npar=fine['npar'],main_order=row['quadrature_order'],check_order=fine['quadrature_order'],delta_r2_A2=fine['mean_r2_A2']-row['mean_r2_A2'],delta_kinetic_eV=fine['mean_kinetic_eV']-row['mean_kinetic_eV'],delta_exterior_048=fine['exterior_048']-row['exterior_048'])
        for kind in ('position','momentum'):
            r=z[f'{kind}_radius'];d=z[f'{kind}_density'];df=zf[f'{kind}_density'];m=z[f'{kind}_measure'];mf=zf[f'{kind}_measure']
            item[f'{kind}_integrated_density_L1_difference']=float(np.trapezoid(4*np.pi*r*r*abs(d-df),r))
            for key,a,af in [('density',d,df),('measure',m,mf)]:
                crossing=row[kind][f'candidate_{key}_1e6'];cf=fine[kind][f'candidate_{key}_1e6']
                valid=False;delta=None
                if crossing is not None and cf is not None:
                    i=int(np.argmin(abs(r-crossing)));delta=float(abs(a[i]-af[i])/max(a.max(),af.max()))
                    valid=bool(delta<2e-7 and abs(crossing-cf)<=2.01*(r[1]-r[0]))
                item[f'{kind}_{key}_1e6']=dict(resolved_by_comparison=valid,main_candidate=crossing,check_candidate=cf,absolute_difference_relative_peak=delta,reported_crossing=crossing if valid else None)
        if var=='old':
            old=next(v for v in archive if v['case']==f'PdD_{T}K_joint')
            item['archived_baseline_r2_A2']=old['position']['mean_square_radius']
            item['archived_baseline_kinetic_eV']=old['mean_kinetic_eV']
            item['new_solver_minus_archive_r2_A2']=row['mean_r2_A2']-old['position']['mean_square_radius']
            item['new_solver_minus_archive_kinetic_eV']=row['mean_kinetic_eV']-old['mean_kinetic_eV']
            item.update(archived_sampling_control(T,z))
            item['recreated_reporting_minus_archive_r2_A2']=item['recreated_archived_reporting_r2_A2']-old['position']['mean_square_radius']
            item['recreated_reporting_minus_archive_kinetic_eV']=item['recreated_archived_reporting_kinetic_eV']-old['mean_kinetic_eV']
        result.append(item)
    dump(OUT/'numerical_evidence.json',result)
    return result


def plots():
    fig,axs=plt.subplots(2,3,figsize=(12,6.5))
    for i,T in enumerate((300,77)):
        zs={v:np.load(OUT/f'D_{T}_{v}_main_distributions.npz') for v in VARIANTS}
        for j,(kind,key) in enumerate([('position','density'),('position','measure'),('momentum','density')]):
            ax=axs[i,j];r=zs['old'][f'{kind}_radius']
            old=zs['old'][f'{kind}_{key}'];a=zs['quadratic'][f'{kind}_{key}'];b=zs['quartic'][f'{kind}_{key}']
            scale=max(old.max(),a.max(),b.max())
            ax.semilogy(r,old/scale,':',c='k',label='old quartic');ax.semilogy(r,a/scale,c='#2468a2',label='new tangent quadratic');ax.fill_between(r,np.minimum(a,b)/scale,np.maximum(a,b)/scale,color='#2468a2',alpha=.25,label='exterior model band')
            if kind=='position':ax.axvline(.48,c='k',lw=.5);ax.set_xlim(0,1.15);ax.set_xlabel('r (Å)')
            else:ax.set_xlim(0,35);ax.set_xlabel('k = p/ℏ (Å⁻¹)')
            ax.set_ylim(1e-7,2);ax.set_title(f'D {T} K: {kind} {key}');ax.set_ylabel('relative to shared peak')
    axs[0,0].legend(fontsize=6);fig.tight_layout();fig.savefig(HERE/'distributions.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(9,3.5))
    for j,T in enumerate((300,77)):
        for var,ls in [('old',':'),('quadratic','-'),('quartic','--')]:
            z=np.load(OUT/f'D_{T}_{var}_main_distributions.npz');ax[j].plot(z['s111_A'],z['marginal111_A_inverse'],ls,label=var)
        ax[j].set_title(f'D {T} K');ax[j].set_xlabel('signed s = Δr · +[111]/√3 (Å)');ax[j].set_ylabel('plane marginal P(s) (Å⁻¹)');ax[j].set_xlim(-.9,.9)
    ax[0].legend(fontsize=7);fig.tight_layout();fig.savefig(HERE/'signed111.png',dpi=170);plt.close(fig)


def page(rows,numerics,potential):
    def select(iso,T,var):
        return next(v for v in rows if v['isotope']==iso and v['temperature_K']==T and v['variant']==var)
    lines=['# PdD: what the signed 0.48 Å cuts change','']
    for T in (300,77):
        old=select('D',T,'old');new=select('D',T,'quadratic');other=select('D',T,'quartic')
        lines.append(f"At **{T} K**, the isolated D mean-square displacement changes from **{old['mean_r2_A2']:.6f} to {new['mean_r2_A2']:.6f} Å²** with the softer exterior, or **{other['mean_r2_A2']:.6f} Å²** with old-quartic curvature outside the last cut. The primary change is {100*(new['mean_r2_A2']/old['mean_r2_A2']-1):+.2f}%.")
    lines+=['','These are predictions of a frozen-neighbour, single-site Born–Oppenheimer well, not measured nuclear tails. The actual cuts constrain three rays, **not a whole 0.48 Å ball**. Continuation spread remains a model uncertainty even where the quantum calculation converges. No fit to yield, profile or S factor is made.','','![Anchored energy and computed force](wells.png)','','## Anchors, regularity and exterior choices','',
    'For each ray, a piecewise cubic Hermite curve in **r** matches both electronic energy E and outward slope −F at 0.12, 0.24 and 0.48 Å. This is genuinely nonquartic; no energy/force compromise is hidden in a PCHIP fit. Independently native ±[110] rows are averaged after converting the fixed-axis force to its outward sign. Native +[100]/+[111] at 0.48 have inversion-derived negative rows, **not second independent calculations**. Old +0.12/+0.24 [100]/[111] and the signed near-K input are retained.',
    '',
    f"The largest signed-native averaging residual is {max(abs(v['energy_residual_eV']) for v in potential['signed_native_and_derived_residuals'] if v['sample_kind']=='native_scf'):.3g} eV in E and {max(abs(v['force_parallel_residual_eV_A']) for v in potential['signed_native_and_derived_residuals'] if v['sample_kind']=='native_scf'):.3g} eV/Å in projected F. This parity discrepancy is not a DFT-error bound. Every ray's original/averaged anchors, old-model residuals and signed residuals are in `results/potential_evidence.json` and `tables/energy_force_anchors.csv`.",
    '',
    f"Below 0.03 Å the exact old cubic quartic polynomial is retained, with common harmonic K={K:.9f} eV/Å². Thus the angular closure is polynomial at the origin, rather than attaching an O(r⁴) coefficient to a singular sextic angular invariant. At larger r set xᵢ=uᵢ²/r², q=x₁x₂+x₂x₃+x₃x₁ and p=x₁x₂x₃, and use V=w₁₀₀V₁₀₀+w₁₁₀V₁₁₀+w₁₁₁V₁₁₁ with w=(1−4q+9p, 4q−36p, 27p). These nonnegative weights sum to one (Schur/AM–GM), match all three cuts, and give a C¹ cubic-symmetric closure. Monotone positive radial segments and positive confining extensions exclude negative angular wells in this declared closure; they do not establish unmeasured angular physics.",
    '',
    'The two new models are identical through 0.48 Å. For t=r−0.48>0, **quadratic** (primary) is E_R+s_Rt+Kt²/2; **quartic** (alternative) is E_R+s_Rt+(K+3a_ray R²)t²/2+a_ray R t³+a_ray t⁴/4 using the positive old ray quartics. Both match E and slope, and confine positively. Their structural spread is not a tiny coefficient-fit error bar.',
    '',
    f"Computed force is the analytic −∇V, including angular derivatives. Nonspecial-direction finite differences, including every splice, give max absolute force/gradient discrepancy {max(v['max_abs_force_plus_energy_gradient_eV_A'] for v in potential['gradient'].values()):.3g} eV/Å at 10⁻⁶ Å step; two-sided C¹ splice diagnostics are recorded separately.",
    '',
    '## Matched 3D quantum result','',
    'The same 3D Cartesian Hermite Galerkin solver treats old and new wells, with all parity sectors, actual isotope masses and canonical weights. The distributions are **genuine eigenstate mixtures and Fourier transforms**, not Gaussian fits or classical momenta. Geometry, source F0=0 and old neighbour graph are unchanged. Positive-octant 3D arrays plus exact reflection symmetry and the saved spectra define the full distributions.',
    '',
    '| D temperature | Well | ⟨r²⟩ (Å²) | ⟨Kinetic⟩ (eV) | P(r>.24 Å) | P(r>.48 Å) | line anisotropy at .48 |',
    '|---|---|---:|---:|---:|---:|---:|']
    for T in (300,77):
        for v in ('old','quadratic','quartic'):
            row=select('D',T,v)
            lines.append(f"| {T} K | {v} | {row['mean_r2_A2']:.6f} | {row['mean_kinetic_eV']:.6f} | {row['exterior_024']:.5f} | {row['exterior_048']:.5f} | {row['position_line_anisotropy_at_048_max_over_min']:.3f} |")
    lines+=['','Anisotropy is max/min of the [100]/[110]/[111] spatial line densities at the stated radius, not a covariance symmetry breaking. Cubic symmetry fixes ⟨x²⟩=⟨y²⟩=⟨z²⟩ and zero signed mean. Full directional curves are retained, including momentum anisotropy.','',
    '![Actual position and momentum distributions](distributions.png)','',
    '**Density is not radial measure:** spatial ρ̄(r) has units Å⁻³, whereas 4πr²ρ̄(r)=dP/dr has units Å⁻¹. Momentum uses k=p/ℏ (Å⁻¹), n̄(k) (Å³) and 4πk²n̄(k) (Å). Tables retain both, their peaks and tails. The plot scales each observable by a shared old/new peak; the band is exterior-model spread, not numerical error. There is no .48 Å→k cutoff.','',
    '![Signed incoming-axis plane marginal](signed111.png)','',
    'The signed [111] table integrates a plane transverse to +[111]/√3: P(s) ds, s=Δr·n, is a location probability (Å⁻¹). It is not the density sampled on the ray. No physical slab offset/orientation or arrival dynamics is inferred from the bulk reference.','',
    '### Numerical sensitivity, separately from continuation','',
    '| D temperature | Well | Δ⟨r²⟩ (Å²) | Δ⟨Kinetic⟩ (eV) | ΔP(r>.48) | position / momentum density 10⁻⁶ |',
    '|---|---|---:|---:|---:|---|']
    for n in numerics:
        def cross(kind):
            x=n[f'{kind}_density_1e6']['reported_crossing']
            return 'unresolved' if x is None else f'{x:.3g}'
        lines.append(f"| {n['temperature_K']} K | {n['variant']} | {n['delta_r2_A2']:+.3g} | {n['delta_kinetic_eV']:+.3g} | {n['delta_exterior_048']:+.3g} | {cross('position')} Å / {cross('momentum')} Å⁻¹ |")
    for n in numerics:
        if n['variant']=='old':
            lines.append(f"The matched old {n['temperature_K']} K solution differs from the archived known positive baseline by {n['new_solver_minus_archive_r2_A2']:+.3g} Å² in ⟨r²⟩ and {n['new_solver_minus_archive_kinetic_eV']:+.3g} eV in kinetic energy because the archive used linearly interpolated radial integration. Reapplying that archived reporting convention to this same old solution reproduces it within {abs(n['recreated_reporting_minus_archive_r2_A2']):.2g} Å² and {abs(n['recreated_reporting_minus_archive_kinetic_eV']):.2g} eV; the old/new table instead uses direct matrix moments.")
    lines+=['','Δ compares 10→12 parity states per Cartesian axis (20→24 HO levels), potential quadrature 48→64, position-grid pitch .020→.016 Å, momentum pitch .60→.48 Å⁻¹ and angular GL24×48→GL32×64. It is an empirical combined numerical difference, **not** the continuation band or a rigorous enclosure. See `results/numerical_evidence.json` for integrated distribution differences and all density/measure crossing decisions. A 10⁻⁶ crossing is printed only if its local difference is below 2×10⁻⁷ of the peak and its location is stable within two output steps; all other crossings remain unresolved. Resolved numerical contours outside .48 Å are still extrapolated physical predictions.','',
    '### Mass and temperature sensitivity in the same PdD electronic surface','',
    '| Isotope | T | Old ⟨r²⟩ | New ⟨r²⟩ range (Å²) | Old kinetic | New kinetic range (eV) |',
    '|---|---|---:|---:|---:|---:|']
    for iso in ('H','T'):
        for T in (300,77):
            old=select(iso,T,'old');a=select(iso,T,'quadratic');b=select(iso,T,'quartic')
            lines.append(f"| {iso} | {T} K | {old['mean_r2_A2']:.6f} | {min(a['mean_r2_A2'],b['mean_r2_A2']):.6f}–{max(a['mean_r2_A2'],b['mean_r2_A2']):.6f} | {old['mean_kinetic_eV']:.6f} | {min(a['mean_kinetic_eV'],b['mean_kinetic_eV']):.6f}–{max(a['mean_kinetic_eV'],b['mean_kinetic_eV']):.6f} |")
    lines+=['','These are isolated mass substitutions, not H/T-specific electronic calculations or new coupled ensembles. They use the main solver resolution; the explicit higher-resolution comparison is for D. Temperature changes the canonical ensemble only, not the cell or thermal expansion.','',
    '## Actual joint-host tagged sensitivity','',
    'Coupled048 consumed the same local E(x), all 801 retained frames at 32/128 beads (300/77 K), the actual co-beaded neighbour field, D mass and old occupied graph. The exact finite-P ratio is exp[−β mean_beads(V_new−V_old)] for **one selected site**, not a 32-site refit. Nine contiguous-block ratio errors include paired old/new covariance. Retained frame ESS is 796–800 and block ESS ≈9; participation is not proof of unseen overlap or fully independent blocks.','',
    '| T | Old MSD | Quadratic MSD | Quartic MSD | Quadratic paired change ± SE (Å²) | old→quadratic P(r>.48) |',
    '|---|---:|---:|---:|---:|---:|']
    coupled=json.loads((HERE/'coupled_sensitivity.json').read_text())
    for c in coupled['cases']:
        v=c['variants'];ms=[v[x]['position']['MSD_A2'] for x in ('old','quadratic','quartic')];ext=[v[x]['position']['bead_fraction_gt_0.48_A']['value'] for x in ('old','quadratic')]
        lines.append(f"| {c['temperature_K']} K | {ms[0]['value']:.6f} | {ms[1]['value']:.6f} | {ms[2]['value']:.6f} | {ms[1]['change_vs_old']:+.6f} ± {ms[1]['paired_change_block_se']:.6f} | {100*ext[0]:.3f}%→{100*ext[1]:.3f}% |")
    lines+=['','The warmer ensemble remains broader under both continuations, and the tagged change is small against the warm–cold separation. The saved **true signed [111] displaced-path characteristic** changes at separation .25 Å by +.00267 ± .00080 (300 K) / +.00253 ± .00028 (77 K) for the primary model. At that witness about 27%/33% of its numerator involves shifted paths with a bead beyond .48 Å. This is therefore a continuation-sensitive momentum diagnostic, **not** a validated momentum tail or all-cell replacement. The direct shifted-target numerator and original-polynomial callback controls agree at ~10⁻¹⁵. No Fourier inversion, Gaussian replacement, new i-PI run, or fabricated excitation is supplied.','',
    '## Handoff and boundary','',
    '`tables/report_table.csv` and `results/report_table.json` are the standalone decisive statistics. Position/momentum bands and signed-111 marginals are direct table handoffs. `coupled_position.csv` retains annular measures, sampling SE and paired differences; `coupled_characteristic_111.csv` retains the actual characteristic and overlap. `well.py` is the shared local potential API. [Replay, exact inputs and units](replay.md), `ARCHIVE_FILES.txt`, `SHA256SUMS` and `capsule.json` bind one unpublished scoped capsule. Governing plan: `coherencev1-g10-pdd-refit048-20261124` (sole writer WideWellDelivery).','',
    '**What remains unmeasured:** off-ray angular structure, every exterior continuation, electronic-structure/systematic errors, connected pair/contact response, surfaces/vacancies, damping and well changes during arrival. No 0.80 Å data or NiH additions are implied. This updates equilibrium one-body inputs only; no collision/yield improvement is claimed and no direct g8 excitation comparison is asserted.','']
    (HERE/'README.md').write_text('\n'.join(lines))


def main():
    TABLES.mkdir(exist_ok=True)
    rows=json.loads((OUT/'summary_main.json').read_text())
    extra=OUT/'summary_isotopes.json'
    if extra.exists():rows.extend(json.loads(extra.read_text()))
    potential=potential_evidence();distribution_tables(rows);numerics=numerical_evidence(rows);plots()
    for row in rows:
        z=np.load(OUT/f'{row["isotope"]}_{row["temperature_K"]}_{row["variant"]}_{"main" if row["isotope"]=="D" else "isotopes"}_distributions.npz')
        r=z['position_radius'];lines=z['position_lines'];index=int(np.argmin(abs(r-.48)))
        row['position_line_anisotropy_at_048_max_over_min']=float(lines[index].max()/lines[index].min())
    coupled=json.loads((HERE/'coupled_sensitivity.json').read_text())
    handoff=dict(potential_api='well.energy_force(x,variant)',position_and_momentum_tables='tables/{D,H,T}_{300,77}_{position,momentum}.csv',signed_plane_marginal='tables/D_{300,77}_signed111.csv',full_3D_distributions='results/D_{300,77}_{old,quadratic,quartic}_main_distributions.npz: positive octants; reflect axes',coupled_joint_marginal='coupled_position.csv',coupled_true_characteristic='coupled_characteristic_111.csv',scope='isolated 3D canonical refit plus tagged-site finite-P sensitivity only, not an all-cell replacement',direct_g8_consumer=None)
    dump(OUT/'report_table.json',dict(isolated=rows,numerical_comparison=numerics,coupled=coupled,handoff=handoff))
    keys=['isotope','temperature_K','variant','mean_r2_A2','mean_kinetic_eV','exterior_024','exterior_048','position_line_anisotropy_at_048_max_over_min']
    table(TABLES/'report_table.csv',keys,[[v[k] for v in rows] for k in keys])
    page(rows,numerics,potential)
    print(json.dumps(dict(rows=rows,numerics=numerics),indent=2))


if __name__=='__main__':main()
