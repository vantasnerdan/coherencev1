"""Produce compact tables, continuation bands, actual reader cases and replay.

Commands: exact, coupled, derive, verify. Production stages are resumable from
saved numerical artifacts; verify is one changed-consumer check after settling.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from model import HERE,ROOT,KB,MASS,full,partial_cases,strain_cases,inventory,serializable,kinetic_coefficient

OUT=HERE/'output'


def dump(path,data):
    Path(path).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def table(path,headers,columns):
    np.savetxt(path,np.column_stack(columns),delimiter=',',header=','.join(headers),comments='',fmt='%.12g')


def compact_exact():
    """Retain distributions, not redundant interpolation volumes."""
    exact=OUT/'exact';tables=OUT/'tables';tables.mkdir(parents=True,exist_ok=True)
    summary=json.loads((exact/'exact_summary.json').read_text())
    for row in summary:
        if row['case'].endswith('basis8'):continue
        path=exact/(row['case']+'.npz');z=np.load(path)
        # Old producer volume is consumed once to export directional density;
        # following replays consume the compact curves or recompute the solver.
        data={k:z[k] for k in z.files if k not in ('rho_A_minus3','n_k_A3')}
        for label,gridkey,rhokey,rkey in [('position','x_A','rho_A_minus3','radial_A'),('momentum','k_A_inverse','n_k_A3','momentum_k_A_inverse')]:
            if rhokey not in z:continue
            interp=RegularGridInterpolator((z[gridkey],)*3,z[rhokey],bounds_error=False,fill_value=0.)
            for direction,v in [('100',[1,0,0]),('110',[1,1,0]),('111',[1,1,1])]:
                vec=np.array(v)/np.linalg.norm(v)
                data[f'{label}_direction{direction}']=interp(z[rkey][:,None]*vec)
            # Export all continuation signed marginals without new spectra.
            if label=='position':
                from exact import signed_marginal
                for direction,v in [('100',[1,0,0]),('110',[1,1,0]),('111',[1,1,1])]:
                    target=exact/f'{row["case"]}_normal{direction}.csv'
                    if not target.exists():
                        s,p,line=signed_marginal(z[gridkey],z[rhokey],v)
                        table(target,['s_A','marginal_A_inverse','on_axis_line_density_A_minus3'],[s,p,line])
        z.close();np.savez_compressed(path,**data)
    for material,isotope in [('Pd','H'),('Pd','D'),('Pd','T'),('Ni','H')]:
        for T in (300,77):
            name=f'{material}{isotope}_{T}K';variants=[np.load(exact/f'{name}_{v}.npz') for v in ('joint','energy','force')]
            base=variants[0];coarse=np.load(exact/f'{name}_basis8.npz')
            for label,rkey,dkey,mkey,tkey,unit in [('position','radial_A','radial_density_A_minus3','radial_measure_A_minus1','radial_tail','A'),
                                                 ('momentum','momentum_k_A_inverse','momentum_density_A3','momentum_measure_A','momentum_tail','A_inverse')]:
                density=np.array([z[dkey] for z in variants]);measure=np.array([z[mkey] for z in variants]);tails=np.array([z[tkey] for z in variants])
                ck='radial_density_A_minus3' if label=='position' else 'momentum_density_A3'
                values=[base[rkey],density[0],density.min(axis=0),density.max(axis=0),measure[0],measure.min(axis=0),measure.max(axis=0),tails[0],tails.min(axis=0),tails.max(axis=0),abs(base[dkey]-coarse[ck])]
                headers=[f'coordinate_{unit}','density','fit_continuation_density_min','fit_continuation_density_max','radial_measure','fit_continuation_measure_min','fit_continuation_measure_max','tail_probability','fit_continuation_tail_min','fit_continuation_tail_max','basis8_vs10_density_delta']
                if label=='position':headers+=['beyond_computed_radius_0p24_A'];values+=[base[rkey]>.24]
                else:headers+=['coordinate_support_A_no_direct_k_boundary'];values+=[np.full(len(base[rkey]),.24)]
                for direction in ('100','110','111'):
                    directional=np.array([z[f'{label}_direction{direction}'] for z in variants])
                    headers += ['direction'+direction+'_density','direction'+direction+'_continuation_min','direction'+direction+'_continuation_max']
                    values += [directional[0],directional.min(axis=0),directional.max(axis=0)]
                table(tables/f'{name}_exact_{label}_band.csv',headers,values)
            for direction in ('100','110','111'):
                cuts=[np.genfromtxt(exact/f'{name}_{v}_normal{direction}.csv',delimiter=',',names=True) for v in ('joint','energy','force')]
                p=np.array([c['marginal_A_inverse'] for c in cuts]);lines=np.array([c['on_axis_line_density_A_minus3'] for c in cuts])
                table(tables/f'{name}_g9_normal{direction}_band.csv',
                      ['s_A','marginal_A_inverse','continuation_min_A_inverse','continuation_max_A_inverse','on_axis_line_density_A_minus3','line_min','line_max','beyond_computed_radial_support'],
                      [cuts[0]['s_A'],p[0],p.min(axis=0),p.max(axis=0),lines[0],lines.min(axis=0),lines.max(axis=0),abs(cuts[0]['s_A'])>.24])
    return summary


def spherical_momentum_density(characteristic,k):
    """Angle-averaged n(k), including its k=0 limit; units A^3."""
    s=characteristic.displacement_A
    def kernel(nodes):
        raw=nodes[None,:]**2*np.sinc(k[:,None]*nodes[None,:]/np.pi)/(2*np.pi**2)
        widths=np.diff(nodes)/2
        raw[:,:-1]*=np.r_[widths[0],widths[:-1]+widths[1:]]
        raw[:,-1]*=widths[-1]
        return raw
    transform=kernel(s);value=transform@characteristic.value
    blocks=characteristic.blocks@transform.T;se=blocks.std(axis=0,ddof=1)/np.sqrt(len(blocks))
    shorter=slice(0,1+3*(len(s)-1)//4);coarse=np.r_[np.arange(0,len(s)-1,2),len(s)-1]
    window=np.abs(kernel(s[shorter])@characteristic.value[shorter]-value)
    mesh=np.abs(kernel(s[coarse])@characteristic.value[coarse]-value)
    return dict(value=value,se=se,window=window,mesh=mesh,floor=3*se+window+mesh)


def relative_cutoff(grid,value,floor,overlap):
    peak=int(np.argmax(value));threshold=value[peak]*1e-6
    candidates=np.flatnonzero((np.arange(len(grid))>peak)&(value<=threshold))
    if not len(candidates):return dict(candidate=None,resolved=None,absolute_numerical_floor=None)
    i=int(candidates[0])
    return dict(candidate=float(grid[i]),resolved=float(grid[i]) if overlap and floor[i]<threshold else None,
                absolute_numerical_floor=float(floor[i]),threshold=float(threshold))


def coupled_distributions(rows=None):
    from coupled_momentum import displaced_path,radial_momentum,momentum_component,radial_survival
    from conditioned import radial,quartic
    directory=OUT/'coupled';tables=OUT/'tables'
    if rows is None:rows=json.loads((directory/'trajectory_summary.json').read_text())
    models={m['name']:m for m in [full('Pd'),full('Ni')]+partial_cases()+strain_cases()}
    saved=OUT/'coupled_momentum_summary.json'
    summaries=json.loads(saved.read_text()) if saved.exists() else []
    for row in rows:
        name=row['case'];m=models[row['model']];site=m['site'];mass=MASS[row['isotope']];T=row['temperature_K']
        if sum(r['case']==name for r in summaries)==3:continue
        z=np.load(directory/(name+'_paths.npz'));x=z['beads_A'];field=z['neighbor_field_eV_A']
        basequartic=quartic(x,m['a'],m['b']);variants=[];projections=[]
        for variant in ('joint','energy','force'):
            v=m if variant=='joint' else dict(m,a=full('Ni' if m['name']=='Ni' else 'Pd',variant)['a'],b=full('Ni' if m['name']=='Ni' else 'Pd',variant)['b'])
            logweights=z['log_weights_joint'] if variant=='joint' and 'log_weights_joint' in z else None
            if variant!='joint':
                key='log_weights_'+variant
                logweights=z[key] if key in z else -(quartic(x,v['a'],v['b'])-basequartic).mean(axis=1)/(KB*T)
            rb=radial(x,field,v,mass,T,logweights)
            params=(v['K'][site],v['f'][site],v['a'],v['b'],field)
            characteristic=displaced_path(x[:,:,None,:],mass,T,0,np.linspace(0,2.4,241),site_coefficients=params,
                                          block_size=len(x)//12,frame_batch=64,reference_stiffness_eV_per_A2=max(.5,float(np.trace(v['K'][site])/3)+.8),
                                          log_frame_weights=logweights)
            k=np.linspace(0,42,211);pdf=radial_momentum(characteristic,k);tail=radial_survival(characteristic,k)
            overlap=(characteristic.block_weight_ess>=4 and characteristic.frame_weight_ess>=50)
            n3d=spherical_momentum_density(characteristic,k)
            position_cutoff=relative_cutoff(rb['radius'],rb['density'],3*rb['density_se'],overlap)
            momentum_cutoff=relative_cutoff(k,n3d['value'],n3d['floor'],overlap)
            measure_cutoff=relative_cutoff(k,pdf.value,pdf.detection_floor,overlap)
            if name=='PdD_300K' and variant=='joint':
                finer=displaced_path(x[:,:,None,:],mass,T,0,characteristic.displacement_A,site_coefficients=params,
                                     block_size=len(x)//12,frame_batch=64,angular_order=11,
                                     reference_stiffness_eV_per_A2=max(.5,float(np.trace(v['K'][site])/3)+.8))
                fine_density=spherical_momentum_density(finer,k);i=np.argmin(abs(k-momentum_cutoff['candidate']))
                angular_path=OUT/'angular_sensitivity.json';angular=json.loads(angular_path.read_text()) if angular_path.exists() else {}
                angular['coupled_PdD_300K']=dict(case=name+'_joint',orders=[7,11],
                    maximum_characteristic_difference=float(abs(characteristic.value-finer.value).max()),
                    cutoff_k_A_inverse=float(k[i]),density7_A3=float(n3d['value'][i]),density11_A3=float(fine_density['value'][i]),
                    density_delta_A3=float(abs(n3d['value'][i]-fine_density['value'][i])),density7_floor_A3=float(n3d['floor'][i]),
                    scope='one actual coupled angular control, not an all-case angular error bound')
                dump(angular_path,angular)
            result=dict(case=name,variant=variant,position_normalization=rb['normalization'],
                        conditional_normalization_log_delta=rb['conditional_log_normalizer_order10_vs12'],
                        coordinate_support_A=m['support_A'],coordinate_support_fraction=rb['support_fraction'],
                        frame_weight_ess=characteristic.frame_weight_ess,block_weight_ess=characteristic.block_weight_ess,
                        min_displaced_numerator_block_ess=float(characteristic.numerator_block_ess.min()),
                        importance_overlap_supported=overlap,
                        position_density_1e6=position_cutoff,momentum_density_1e6=momentum_cutoff,
                        momentum_measure_1e6=measure_cutoff,
                        position_density_peak_A=float(rb['radius'][np.argmax(rb['density'])]),
                        position_measure_peak_A=float(rb['radius'][np.argmax(rb['measure'])]),
                        momentum_density_peak_A_inverse=float(k[np.argmax(n3d['value'])]),
                        momentum_measure_peak_A_inverse=float(k[np.argmax(pdf.value)]),
                        k2_k4=characteristic.k2_k4.tolist(),k2_k4_se=characteristic.k2_k4_se.tolist(),
                        mean_kinetic_eV=float(kinetic_coefficient(mass)*characteristic.k2_k4[0]),
                        kinetic_standard_error_eV=float(kinetic_coefficient(mass)*characteristic.k2_k4_se[0]),
                        transform_resolution_A_inverse=pdf.resolution_A_inv,
                        native_kinetic_factor_eV_A2=characteristic.kinetic_factor_eV_A2,
                        native_beta_over_beads_eV_inverse=characteristic.beta_over_beads_eV_inv,
                        positive_signal_last_k_A_inverse=float(k[pdf.positive_signal][-1]) if pdf.positive_signal.any() else None,
                        continuation_scope='all-site quartic path reweight' if 'log_weights_energy' in z else 'source-site quartic only; host continuation held fixed',
                        isotope_semantics=row.get('isotope_semantics','homogeneous isotope replacement'),
                        host_isotope=row.get('host_isotope',row['isotope']),joint_excitation=None)
            summaries.append(result);variants.append((rb,pdf,tail,n3d))
            table(tables/f'{name}_{variant}_characteristic.csv',['displacement_A','characteristic','block_SE','numerator_block_ESS'],[characteristic.displacement_A,characteristic.value,characteristic.se,characteristic.numerator_block_ess])
            projection=displaced_path(x[:,:,None,:],mass,T,0,np.linspace(0,2.4,241),direction=[1,1,1],site_coefficients=params,
                                      block_size=len(x)//12,frame_batch=64,reference_stiffness_eV_per_A2=max(.5,float(np.trace(v['K'][site])/3)+.8),
                                      log_frame_weights=logweights)
            projections.append(momentum_component(projection,np.linspace(-42,42,421)))
            if variant=='joint':
                # g9 receives an actual signed coordinate marginal. Empty bins
                # remain unresolved; not invented zero-density support.
                for direction,vec in [('100',[1,0,0]),('110',[1,1,0]),('111',[1,1,1]),('measured',m['axis'])]:
                    vec=np.array(vec)/np.linalg.norm(vec);s=x@vec;edges=np.linspace(-1,1,201);width=np.diff(edges)
                    hist=[];den=[];bs=len(x)//12
                    fw=np.ones(len(x)) if logweights is None else np.exp(logweights-logweights.max())
                    for b in range(12):
                        sl=slice(b*bs,(b+1)*bs)
                        hist.append(np.histogram(s[sl].ravel(),edges,weights=np.repeat(fw[sl],x.shape[1]))[0]/(bs*x.shape[1])/width)
                        den.append(fw[sl].mean())
                    hist=np.array(hist);den=np.array(den);value=hist.sum(axis=0)/den.sum()
                    linear=value+(hist-den[:,None]*value)/den.mean();se=linear.std(axis=0,ddof=1)/np.sqrt(12)
                    continuations=[value]
                    for alternative in ('energy','force'):
                        av=full('Ni' if m['name']=='Ni' else 'Pd',alternative)
                        key='log_weights_'+alternative
                        lw=z[key] if key in z else -(quartic(x,av['a'],av['b'])-basequartic).mean(axis=1)/(KB*T)
                        fw=np.exp(lw-lw.max())
                        raw=np.histogram(s.ravel(),edges,weights=np.repeat(fw,x.shape[1]))[0]
                        continuations.append(raw/(fw.sum()*x.shape[1])/width)
                    band=np.array(continuations);centers=(edges[:-1]+edges[1:])/2
                    table(tables/f'{name}_g9_sampled_normal{direction}.csv',['s_A','marginal_A_inverse','continuation_min','continuation_max','block_SE','sampled_nonzero','beyond_computed_projected_support'],[centers,value,band.min(axis=0),band.max(axis=0),se,value>0,abs(centers)>m['support_A']])
        rb,pdf,tail,n3d=variants[0];position=np.array([v[0]['density'] for v in variants]);measures=np.array([v[0]['measure'] for v in variants]);survival=np.array([v[0]['tail'] for v in variants]);mom=np.array([v[1].value for v in variants]);mtail=np.array([v[2].value for v in variants]);nvalues=np.array([v[3]['value'] for v in variants])
        overlap=bool(summaries[-3]['importance_overlap_supported'])
        table(tables/f'{name}_coupled_position_band.csv',['r_A','density_A_minus3','continuation_min','continuation_max','block_SE','radial_measure_A_inverse','measure_min','measure_max','tail_probability','tail_min','tail_max','positive_block_signal','beyond_computed_support','importance_overlap_supported'],
              [rb['radius'],position[0],position.min(axis=0),position.max(axis=0),rb['density_se'],measures[0],measures.min(axis=0),measures.max(axis=0),survival[0],survival.min(axis=0),survival.max(axis=0),rb['resolved']&overlap,rb['radius']>m['support_A'],np.full(len(rb['radius']),overlap)])
        table(tables/f'{name}_coupled_momentum_band.csv',['k_A_inverse','radial_measure_A','continuation_min','continuation_max','block_SE','window_delta','mesh_delta','positive_signal','survival_probability','survival_min','survival_max','survival_SE','survival_positive_signal','free_end_action','importance_overlap_supported','density_A3','density_continuation_min','density_continuation_max','density_SE','density_numerical_lower','density_numerical_upper'],
              [pdf.k_A_inv,mom[0],mom.min(axis=0),mom.max(axis=0),pdf.se,pdf.cutoff_delta,pdf.mesh_delta,pdf.positive_signal&overlap,mtail[0],mtail.min(axis=0),mtail.max(axis=0),tail.se,tail.positive_signal&overlap,pdf.free_end_action,np.full(len(pdf.k_A_inv),overlap),n3d['value'],nvalues.min(axis=0),nvalues.max(axis=0),n3d['se'],n3d['value']-n3d['floor'],n3d['value']+n3d['floor']])
        signed=projections[0];projected=np.array([p.value for p in projections])
        table(tables/f'{name}_momentum_projection111.csv',['signed_k_A_inverse','component_density_A','continuation_min','continuation_max','block_SE','window_delta','mesh_delta','positive_signal','coordinate_support_A_no_direct_k_boundary'],[signed.k_A_inv,signed.value,projected.min(axis=0),projected.max(axis=0),signed.se,signed.cutoff_delta,signed.mesh_delta,signed.positive_signal&overlap,np.full(len(signed.k_A_inv),m['support_A'])])
        print(name,'momentum',summaries[-3]['mean_kinetic_eV'],flush=True)
        dump(saved,summaries)
    dump(OUT/'coupled_momentum_summary.json',summaries)
    return summaries


def consumers():
    """Exercise the actual fixed-profile g8 reader and documented g9 mapping."""
    source=ROOT/'g8_why/collision/source';sys.path.insert(0,str(source))
    from bound_forward import BoundSpectatorKernel
    kernel=BoundSpectatorKernel(ROOT,angular_order=32)
    joint=json.loads((OUT/'exact/PdD_300K_joint_g8_joint.json').read_text());states=joint['states']
    # Use actual spectrum/kinetic positive joint nodes, not a marginal with E=0.
    energies=[.25];fold=kernel.fold('PdD',energies,states,removal_eV=0.)
    p=np.array([sum(s['probabilities']) for s in states]);maxexc=max(s['excitation_eV'] for s in states)
    Kmax=max(max(s['kinetic_eV']) for s in states)
    # Match the existing spectator envelope readout with the same joint measure
    # at its denominator and species/depth profile; no profile fit or removal fit.
    cache={}
    def physical(E,ch):
        key=(E,ch)
        if key not in cache:cache[key]=float(kernel.fold('PdD',[E],states,channel=ch,removal_eV=0.)[0])
        return cache[key]
    reported=float(kernel.common.readout('PdD',np.array(energies),physical,current_power=0)[0])
    envelope=json.loads((ROOT/'g8_why/collision/evidence/envelope.json').read_text())
    def find_rows(obj):
        if isinstance(obj,dict):
            if obj.get('material')=='PdD' and obj.get('CM_keV')==.25 and 'reported_upper' in obj:yield obj
            for v in obj.values():yield from find_rows(v)
        elif isinstance(obj,list):
            for v in obj:yield from find_rows(v)
    matching=list(find_rows(envelope));comparison=[dict(nonnegative_common_removal=r['nonnegative_common_removal'],reported_upper=r['reported_upper'],ratio_to_upper={k:reported/float(v) for k,v in r['reported_upper'].items() if k=='0'}) for r in matching]
    g8=dict(actual_reader='g8_why/collision/source/bound_forward.py:BoundSpectatorKernel.fold',material='PdD',temperature_K=300,
            energies_cm_keV=energies,fold_p_per_C=fold.tolist(),reported_current_power0=reported,
            joint_probability=sum(sum(s['probabilities']) for s in states),joint_excitation_support_eV=[0,maxexc],kinetic_node_support_eV=[0,Kmax],
            angular_assumption='isotropized conditional on each exact single-site excitation',removal_eV=0.,
            removal_scope='separate assumed favorable zero, not a binding-energy measurement',spectator_envelope_comparison=comparison)
    normal=np.genfromtxt(OUT/'tables/PdD_300K_g9_normal111_band.csv',delimiter=',',names=True)
    s=normal['s_A'];p=normal['marginal_A_inverse'];normalization=float(np.trapezoid(p,s));Z0=1.0
    # A declared example speed, not a source-model change or beam calibration.
    speed_A_fs=1.;z=Z0+s;delay=s/speed_A_fs
    table(OUT/'tables/g9_mapping_example.csv',['s_A','marginal_A_inverse','Z0_A','z_target_A','speed_A_fs','arrival_shift_fs'],[s,p,np.full(len(s),Z0),z,np.full(len(s),speed_A_fs),delay])
    g9=dict(normalization=normalization,registered_site_A=Z0,normal_in=[float(1/np.sqrt(3))]*3,
            frame='bulk crystal Cartesian x100 y010 z001; n_in=+[111]/sqrt(3); s=displacement dot n_in, positive into chosen metal half-space',
            mapping='z_target=Z0+s, arrival shift=s/v at fixed ballistic speed; reverse sign for the opposite slab normal',
            consumer_scope='table import and frozen-ballistic mapping only; pinned g9 generic lab z has no crystallographic surface or density CSV API',
            transverse_meaning='marginal integrates whole plane perpendicular to n; on_axis_line_density separately retained, not interchangeable',
            surface_scope='not a slab target or target offset; actual slab geometry remains a future input')
    dump(OUT/'consumer_cases.json',dict(g8=g8,g9=g9));return dict(g8=g8,g9=g9)


def pages(exact,coupled,consumer):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    tables=OUT/'tables';fig,axes=plt.subplots(4,2,figsize=(11,14))
    for row,(material,isotope) in enumerate([('Pd','H'),('Pd','D'),('Pd','T'),('Ni','H')]):
        for T,color in [(300,'tab:red'),(77,'tab:blue')]:
            name=f'{material}{isotope}_{T}K'
            for col,label in enumerate(('position','momentum')):
                d=np.genfromtxt(tables/f'{name}_exact_{label}_band.csv',delimiter=',',names=True);x=d[d.dtype.names[0]]
                axes[row,col].plot(x,d['density']/d['density'].max(),color=color,label=f'{T} K')
                axes[row,col].fill_between(x,d['fit_continuation_density_min']/d['density'].max(),d['fit_continuation_density_max']/d['density'].max(),alpha=.25,color=color)
                if col==0:axes[row,col].axvline(.24,color='black',linestyle='--',linewidth=.8);axes[row,col].axvspan(.24,1,color='gray',alpha=.08)
                axes[row,col].set(yscale='log',ylim=(1e-8,1.5),xlim=(0,.9 if col==0 else 36),title=f'{material}{isotope}: {label}',xlabel='r (Å)' if col==0 else '|k| (Å⁻¹)',ylabel='angle-averaged density / peak')
                axes[row,col].legend()
    fig.suptitle('Single-site non-Gaussian tails; shaded bands = energy/force/joint quartic continuations\nPosition dashed line: actual positive computed cuts end at 0.24 Å; no k analogue')
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(HERE/'tails.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for row,material in enumerate(('PdD','NiH')):
        for T,color in ((300,'tab:red'),(77,'tab:blue')):
            name=f'{material}_{T}K'
            for col,kind in enumerate(('position','momentum')):
                isolated=np.genfromtxt(tables/f'{name}_exact_{kind}_band.csv',delimiter=',',names=True)
                d=np.genfromtxt(tables/f'{name}_coupled_{kind}_band.csv',delimiter=',',names=True)
                x=d['r_A'] if col==0 else d['k_A_inverse']
                value=d['density_A_minus3'] if col==0 else d['density_A3']
                lo=d['continuation_min'] if col==0 else d['density_continuation_min']
                hi=d['continuation_max'] if col==0 else d['density_continuation_max']
                floor=3*d['block_SE'] if col==0 else value-d['density_numerical_lower']
                valid=(value>floor)&d['importance_overlap_supported'].astype(bool)
                failures=np.flatnonzero(~valid)
                if len(failures):valid[failures[0]:]=False
                scale=isolated['density'].max();ax=axes[row,col]
                ax.plot(isolated[isolated.dtype.names[0]],isolated['density']/scale,color=color,linestyle='--',label=f'{T} K isolated')
                ax.plot(x,np.where(valid,value/scale,np.nan),color=color,label=f'{T} K coupled')
                ax.fill_between(x,np.where(valid,lo/scale,np.nan),np.where(valid,hi/scale,np.nan),color=color,alpha=.22)
                ax.plot(x,floor/scale,color=color,linestyle=':',linewidth=.7,label=f'{T} K numerical floor')
                if col==0:ax.axvline(.24,color='black',linestyle='-.',linewidth=.8);ax.axvspan(.24,1,color='gray',alpha=.06)
                ax.set(yscale='log',ylim=(1e-8,1.5),xlim=(0,.9 if col==0 else 32),title=f'{material} {kind}',
                       xlabel='r (Å)' if col==0 else '|k| (Å⁻¹)',ylabel='density / isolated peak')
                ax.legend(fontsize=7)
    fig.suptitle('Single-site and coupled anharmonic shapes\nShading: fit-continuation spread; dotted: sampling/window/mesh floor; unresolved lines stop')
    fig.tight_layout(rect=(0,0,1,.93));fig.savefig(HERE/'neighbours.png',dpi=150);plt.close(fig)
    exactmap={r['case']:r for r in exact};trajectory=json.loads((OUT/'coupled/trajectory_summary.json').read_text());cmap={(r['case'],r['variant']):r for r in coupled}
    lines=['# Bulk quantum-site result','', '![Position and momentum continuation bands](tails.png)','', '![Single-site versus coupled non-Gaussian shapes](neighbours.png)','', 'Each position curve marks the 0.24 Å end of the positive computed full-loading cuts. Momentum bands compare the same potential continuations over the entire k range; there is no coordinate-cutoff-to-k boundary. These are finite fitted-well results, not exact DFT tails.','',
           '| Case | T (K) | Single ⟨r²⟩ (Å²) | Coupled ⟨r²⟩ ± block SE (Å²) | Single ⟨K⟩ (eV) | Coupled ⟨K⟩ (eV) | Probability inside 0.24 Å, single |',
           '|---|---:|---:|---:|---:|---:|---:|']
    for tr in trajectory:
        if tr['model'] not in ('Pd','Ni') or '_P' in tr['case']:continue
        er=exactmap[tr.get('exact_case',tr['case'])+'_joint'];cr=cmap[(tr['case'],'joint')]
        lines.append(f'| {tr["case"]} | {tr["temperature_K"]} | {er["position"]["mean_square_radius"]:.6g} | {tr["mean_square_radius_A2"]:.6g} ± {tr["r2_standard_error_A2"]:.2g} | {er["mean_kinetic_eV"]:.6g} | {cr["mean_kinetic_eV"]:.6g} | {er["position"]["support_fraction"]:.5g} |')
    lines += ['', '### Isolated-reference peaks and contours (not coupled crossings)','', '| Isolated case | Position density / measure peak (Å) | Position 10⁻⁶ density / measure (Å) | Momentum density / measure peak (Å⁻¹) | Momentum 10⁻⁶ density / measure (Å⁻¹) |',
              '|---|---:|---:|---:|---:|']
    for er in exact:
        if er.get('variant')!='joint':continue
        p=er['position'];k=er['momentum'];case=er['case'].removesuffix('_joint')
        lines.append(f'| {case} | {p["density_peak_radius"]:.3g} / {p["radial_measure_peak_radius"]:.3g} | {p["density_cutoff_relative_1e6"]:.3g} / {p["radial_measure_cutoff_relative_1e6"]:.3g} | {k["density_peak_radius"]:.3g} / {k["radial_measure_peak_radius"]:.3g} | {k["density_cutoff_relative_1e6"]:.3g} / {k["radial_measure_cutoff_relative_1e6"]:.3g} |')
    lines += ['', '## Axis shape and isotope meaning','', 'Cubic symmetry fixes the full-model mean to zero and its exact second moments to ⟨x²⟩=⟨y²⟩=⟨z²⟩=⟨r²⟩/3; the sampled tensor/ratio below shows finite-chain uncertainty, not a measured breaking of cubic symmetry. The non-Gaussian angular shape still differs along [100]/[110]/[111]. D in PdD and H in NiH are directly sampled homogeneous host references. H/T labelled Dhost are mass tracers, from centroid-scaled paths with normalized importance weights in the D-isotope host; neighbours remain coupled. The sampled homogeneous H-host reference is a different case, not a tracer check.',
              '', '| Case | Type | Sample ⟨x²⟩, ⟨y²⟩, ⟨z²⟩ (Å²) | Sample covariance λmax/λmin (cubic truth 1) | Sample mean along source axis (Å) |',
              '|---|---|---:|---:|---:|']
    models={m['name']:m for m in [full('Pd'),full('Ni')]+partial_cases()+strain_cases()};quantities=[]
    for tr in trajectory:
        if '_P' in tr['case']:continue
        m=models[tr['model']];mean=np.array(tr['mean_A']);second=np.array(tr['second_moment_A2'])
        covariance=second-np.outer(mean,mean);eigenvalues=np.linalg.eigvalsh(covariance)
        ratio=float(eigenvalues[-1]/eigenvalues[0]);axis_mean=float(mean@m['axis'])
        kind=tr.get('isotope_semantics','homogeneous isotope replacement')
        diag=np.diag(second);mom=cmap[(tr['case'],'joint')]
        kparallel=float(m['axis']@m['K'][m['site']]@m['axis'])
        quartic_parallel=(m['a']-m['b'])*np.sum(m['axis']**4)+m['b']
        variance_parallel=float(m['axis']@covariance@m['axis'])
        lines.append(f'| {tr["case"]} | {kind} | {diag[0]:.5g}, {diag[1]:.5g}, {diag[2]:.5g} | {ratio:.4g} | {axis_mean:.5g} |')
        quantities.append(dict(case=tr['case'],case_type=kind,temperature_K=tr['temperature_K'],isotope=tr['isotope'],host_isotope=tr.get('host_isotope',tr['isotope']),
            occupied_sites=len(m['K']),original_probe_index=int(m['original_indices'][m['site']]),
            axis_x=m['axis'][0],axis_y=m['axis'][1],axis_z=m['axis'][2],
            force_parallel_eV_A=float(m['f'][m['site']]@m['axis']),
            K_parallel_eV_A2=float(m['axis']@m['K'][m['site']]@m['axis']),sumJ_source_eV_A2=m['sumJ_source'],
            coupling_over_local_stiffness=m['sumJ_source']/kparallel,
            quartic_over_quadratic_at_axis_variance=quartic_parallel*variance_parallel/kparallel,
            harmonic_quantum_over_kBT=np.sqrt(2*kinetic_coefficient(MASS[tr['isotope']])*kparallel)/(KB*tr['temperature_K']) if kparallel>0 else None,
            vacancy_fraction=1-len(m['K'])/32,
            mean_x_A=mean[0],mean_y_A=mean[1],mean_z_A=mean[2],mean_measured_axis_A=axis_mean,
            x2_A2=diag[0],y2_A2=diag[1],z2_A2=diag[2],r2_A2=float(np.trace(second)),sample_covariance_anisotropy=ratio,
            r2_block_SE_A2=tr['r2_standard_error_A2'],
            moment_scope='finite-path weighted sample; cubic model has exact zero mean and isotropic second tensor',
            mean_k2_A_inverse2=mom['k2_k4'][0],mean_k4_A_inverse4=mom['k2_k4'][1],
            mean_kinetic_eV=mom['mean_kinetic_eV'],kinetic_block_SE_eV=mom['kinetic_standard_error_eV'],
            frame_weight_ess=mom['frame_weight_ess'],block_weight_ess=mom['block_weight_ess'],
            importance_overlap_supported=mom['importance_overlap_supported'],
            coordinate_support_A=m['support_A'],coordinate_support_fraction=mom['coordinate_support_fraction']))
    lines += ['', '## Coupled tail resolution and reweighting','',
              'These are numerical detection limits only; bead, PES and overlap sensitivities remain separate. A missing 10⁻⁶ crossing is printed as unresolved, never extended with a Gaussian. Frame/block participation is not an autocorrelation-adjusted independent-sample count.',
              '', '| Case | Weight ESS frames / blocks | ⟨K⟩ ± block SE (eV) | Position-density 10⁻⁶ crossing (Å) | Momentum-density 10⁻⁶ crossing (Å⁻¹) |',
              '|---|---:|---:|---:|---:|']
    for tr in trajectory:
        if '_P' in tr['case']:continue
        cr=cmap[(tr['case'],'joint')]
        show=lambda value:'unresolved' if value is None else f'{value:.4g}'
        lines.append(f'| {tr["case"]} | {cr["frame_weight_ess"]:.4g} / {cr["block_weight_ess"]:.4g} | {cr["mean_kinetic_eV"]:.5g} ± {cr["kinetic_standard_error_eV"]:.2g} | {show(cr["position_density_1e6"]["resolved"])} | {show(cr["momentum_density_1e6"]["resolved"])} |')
    lines += ['', '## Installed SSCHA: isolated-site anharmonic stiffness','',
              'The free-energy curvature below is not a discrete excitation level. The stationary Gaussian is supplied by analytic self-consistency, then the installed SSCHA force gradient and full Hessian are evaluated in the pinned order. Inversion symmetry makes the odd Hessian correction vanish for this isolated even well.',
              '', '| Case | Bare local Kii (eV/Å²) | SSCHA Hessian trace/3 (eV/Å²) | Gradient residual / error norm (eV/Å²) | Force-rounding gradient change (eV/Å²) |',
              '|---|---:|---:|---:|---:|']
    for ss in json.loads((OUT/'coupled/sscha_reference.json').read_text()):
        material='Ni' if ss['case'].startswith('Ni') else 'Pd';bare=full(material)['K'][0,0,0]
        stiffness=float(np.trace(np.array(ss['free_energy_hessian_relative_eV_A2']))/3)
        lines.append(f'| {ss["case"]} | {bare:.6g} | {stiffness:.6g} | {ss["SSCHA_gradient_stiffness_residual_eV_A2"]:.4g} / {ss["SSCHA_gradient_error_norm_eV_A2"]:.4g} | {ss["force_rounding_1e_minus4_eV_A_gradient_change_eV_A2"]:.3g} |')
    lines += ['', '## Separate partial-loading and strain Gaussian contrasts','',
              'These are stationary Gaussian variational **model contrasts**, not sampled anharmonic tails and not new measured local stiffnesses. The measured source K/F0 and its original signed axis remain separate inputs. Quartics and unmeasured transverse/bond closures are transferred explicitly.',
              '', '| Case | Actual occupied sites | Source original index | Signed source axis | Local K input: measured partial / Γ closure (eV/Å²) | F0 along axis (eV/Å) | Gaussian mean along axis (Å) | Gaussian ⟨r²⟩ (Å²) | Bare minimum eigenvalue (eV/Å²) |',
              '|---|---:|---:|---|---:|---:|---:|---:|---:|']
    contrasts=json.loads((OUT/'coupled/model_contrasts.json').read_text());contrast_rows=[]
    gaussian=json.loads((OUT/'coupled/gaussian_reference.json').read_text())
    reference=next(r for r in gaussian if r['case']=='Pd' and r['temperature_K']==300 and r['mass_amu']==MASS['D'])
    contrasts=[reference]+contrasts
    for contrast in contrasts:
        m=models[contrast['case']];axis=m['axis'];mean=np.array(contrast['local_centroid_A']);C=np.array(contrast['local_covariance_A2'])
        kval=float(axis@m['K'][m['site']]@axis);force=float(m['f'][m['site']]@axis);meanaxis=float(mean@axis);r2=float(np.trace(C)+mean@mean)
        signed_axis='(' + ','.join(f'{v:.4g}' for v in axis) + ')'
        lines.append(f'| {m["name"]} | {len(m["K"])} | {int(m["original_indices"][m["site"]])} | {signed_axis} | {kval:.6g} | {force:.6g} | {meanaxis:.6g} | {r2:.6g} | {contrast["min_bare_graph_eigenvalue_eV_A2"]:.6g} |')
        contrast_rows.append(dict(case=m['name'],case_type='Gaussian variational contrast, not anharmonic distribution',temperature_K=300,isotope='D',
            occupied_sites=len(m['K']),original_probe_index=int(m['original_indices'][m['site']]),axis_x=axis[0],axis_y=axis[1],axis_z=axis[2],
            source_K_parallel_eV_A2=kval,source_force_parallel_eV_A=force,sumJ_source_eV_A2=m['sumJ_source'],
            K_input_kind='bulk Kii + calibrated deltaGamma with frozen J' if m['name'].startswith('g3_') else 'measured source-axis curvature',
            Gaussian_mean_axis_A=meanaxis,Gaussian_r2_A2=r2,min_bare_graph_eigenvalue_eV_A2=contrast['min_bare_graph_eigenvalue_eV_A2'],
            local_support_A=m['support_A'],input_source=m['source']))
    with (tables/'model_contrasts.csv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(contrast_rows[0]));writer.writeheader();writer.writerows(contrast_rows)
    with (tables/'all_quantities.csv').open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(quantities[0]));writer.writeheader();writer.writerows(quantities)
    g8=consumer['g8'];envelope_row=next(r for r in g8['spectator_envelope_comparison'] if not r['nonnegative_common_removal'])
    lines += ['', f'For the actual 300 K PdD single-site joint measure, excitation support is 0–{g8["joint_excitation_support_eV"][1]:.4g} eV and kinetic-node support 0–{g8["kinetic_node_support_eV"][1]:.4g} eV: at 0.25 keV the unchanged g8 fixed-profile readout is {g8["reported_current_power0"]:.5g}, or {envelope_row["ratio_to_upper"]["0"]:.5g} of its zero-removal spectator-envelope upper value over the declared 0–100 eV excitation support (current power 0).']
    lines += ['', '## Meaning of the tables','',
              '`output/tables/*exact_position_band.csv`: angle-averaged spatial density ρ̄(r) (Å⁻³), its per-radius measure 4πr²ρ̄(r) (Å⁻¹), both peaks and 10⁻⁶ relative-to-own-peak cutoffs in `exact_summary.json`, and integrated exterior probability. The directional columns are on-axis densities, not plane marginals. The momentum analogue uses k=p/ℏ (Å⁻¹), density n̄(k) (Å³), per-radius measure 4πk²n̄(k) (Å), and K=ℏ²k²/(2m) with the actual isotope mass. Isolated isotope references change the nucleus mass in the unchanged PdD-source Born–Oppenheimer potential and registered geometry; coupled host/tracer meanings are labelled separately. No isotope-specific DFT is asserted.',
              '', '`*coupled_position_band.csv` is a conditional quartic bead marginal, not a Gaussian fit. Its spring/quadratic factor is integrated against the quartic and normalized at two Hermite orders; block uncertainty and unresolved rare bins remain visible. `*coupled_momentum_band.csv` is an unwindowed displaced-path transform with signed estimates (including negative noisy bins), block SE, mesh/window differences and free-end-action diagnostics. Importance-overlap diagnostics are separate: the reported support flag warns below 50 participating frames or four participating time blocks, not a physical enclosure. Continuation spread is not sampling error. No histogram zero or unresolved Fourier bin is called a resolved 10⁻⁶ tail. The finite basis comparison is recorded at the physical cutoff, not extrapolated into a Coulomb nuclear tail.',
              '', 'The coupled momentum table has **no joint excitation assignment**. Only `output/exact/PdD_*_joint_g8_joint.json` carries the true single-site canonical joint excitation/kinetic quadrature. It is isotropized conditional on each excitation for the existing g8 angular reader; removal is a separate explicit zero assumption. The actual fixed-profile reader case is `output/consumer_cases.json`; its declared excitation support is far below the existing 100 eV spectator support, and the recorded ratio locates this calculation within that envelope without changing density, profile or S factor.',
              '', 'For g9 use `*_g9_normal111_band.csv`: signed s=Δr·n_in in Å, n_in=+[111]/√3, positive into the declared metal half-space; P(s) has units Å⁻¹. [100]/[110] context and on-axis line density are separate columns/files. A plane marginal is not a trace on the incoming ray. `g9_mapping_example.csv` imports the actual table and applies z_target=Z0+s, Δt=s/v for a declared frozen ballistic example. The pinned g9 model has generic lab z only: this does not assign the separately commissioned slab orientation/offset to a bulk material surface.',
              '', '## Inputs, partial loading, strain','',
              '`output/input_inventory.json` retains exact occupied original indices, axes, forces, folded graphs, source hashes, bare eigenvalues and closures. Adopted force Kii is local; scalar J is each source-shell mean multiplied by the actual occupied graph, once. Kii−ΣJ is diagnostic only, never subtracted from Kii again. Finite folded shell multiplicities are 12,3,12,3,1, not infinite FCC counts. Bare negative modes are not clamped; the bounded positive cubic quartic supplies stabilization.',
              '', 'Pd has measured ±0.03 Å [100] and +0.12/+0.24 Å [100]/[111]. No [110] scan: K4,110=(K4,100+3K4,111)/4 is a cubic prediction. NiH uses the updated source column/Kii and one +0.24 Å [100] quartic; [111] is explicitly the older g2 value. The displayed full-loading continuation band spans joint, energy-only and force-only quartic fits, not all possible continuations; missing larger signed cuts remain input uncertainty.',
              '', 'Partial models distinguish actual 0.875 spread clean/near, 0.75 bunched clean/near and 0.875 bunched near probes; there is no invented missing clean case or relaxation. The signed original probe vector and nonzero full F0 are retained. Only its measured parallel curvature is changed; transverse and unprobed K use full-loading closure, quartics are transferred from x1. Source-shell means are symmetrically transferred over the actual occupied set; a source shell with no occupied neighbour has no measured mean, and its other-site bonds explicitly use the x1 mean. Agreement at the measured source is not a measured global Hessian. Their 0.03 Å small-probe support is marked separately. The separate Gaussian contrasts do not supply anharmonic tail distributions.',
              '', 'Strain contrasts use only calibrated changes of the g3_skew collective Γ tensors: local K=bulk force-column Kii+[KΓ(strain)−KΓ(reference)] under frozen unstrained J, with transferred quartics. The k12 collective reference is not substituted for the k3 local Kii. No extrapolated EMT map enters. Missing strained local large displacements are an input uncertainty. Partial/strain rows are 300 K D Gaussian model contrasts, not additional chemical-profile production cases.',
              '', '## Numerical reading','',
              'Angular comparisons are in `output/angular_sensitivity.json`: the higher-rule NiH77 isolated check and one actual coupled PdD300 seven/eleven-order comparison. At the latter numerical cutoff the angular change lies below the sampling/window/mesh floor; this is not an all-case angular error bound. Inside the 0.24 Å radial bracket is not the same as measured coverage at every direction.',
              'The low temperature is 77 K: βℏω is larger, so the full cases use 128 primitive beads versus 32 at 300 K; PdD half-bead comparisons retain the same physical consumer. This controls the imaginary-time step at the local vibrational scale, not a guarantee for an unresolved high-k bin. Exact reference uses 10 parity states per axis versus 8, projected coordinate powers and Fourier-Hermite support; reported contours are finite-grid estimates (0.02 Å coordinate and 0.7 Å⁻¹ momentum interpolation pitch). Installed SSCHA provides local free-energy Hessian/gradient; analytic coupled Gaussian stationarity is separately labelled and never replaces the non-Gaussian distributions. Force-rounding sensitivity is a numerical provider check, not a DFT-force error bar. Hessian follows force evaluation and gradient, because the installed Hessian routine converts ensemble units.',
              '', '### Three largest omissions','', '1. Electron–phonon damping.', '2. Equilibrium distributions, not real-time dynamics.', '3. No change of well during an ion’s arrival.',
              '', '### Input/model uncertainties (separate)','', 'Missing [110] and larger signed full-loading cuts, old-settings Ni [111], transferred partial/strained quartics, assumed transverse/bond closures, fixed-metal environment and registered unrelaxed partial sites; source mesh/finite-cell/fit uncertainties are not numerical solver error bars. Temperature changes the nuclear ensemble only: supplied geometry/loading remain fixed, not a thermal-expansion prediction. Nuclei are quantum distinguishable localized-site degrees of freedom; inter-well exchange/hopping is not included in this fixed-well Hamiltonian. Future surface/vacancy inputs are not part of this bulk deliverable.',
              '', 'Replay and the exact archive whitelist are in [replay.md](replay.md).']
    verification=OUT/'verification.json'
    if verification.exists():
        evidence=json.loads(verification.read_text())
        low=next(r for r in evidence['PdD_bead_sensitivity'] if r['temperature_K']==77)
        tracer=next(r for r in evidence['coupled_consumer_estimators'] if r['case']=='PdH_Dhost_77K')
        lines += ['', '### Observed numerical sensitivity','',
                  f'At 77 K the PdD 64→128 bead kinetic change is {1e3*low["mean_kinetic_difference_eV"]:.3g} meV, versus combined block SE {1e3*low["combined_kinetic_block_SE_eV"]:.3g} meV; its far per-radius tail is unresolved at both bead counts.',
                  f'The 77 K H tracer momentum/centroid-virial estimates differ by {1e3*tracer["kinetic_difference_eV"]:.3g} meV, versus a conservative sum of block SEs {1e3*tracer["conservative_sum_kinetic_block_SE_eV"]:.3g} meV. This is live sampling/reweighting uncertainty, not claimed continuum agreement. The tabulated crossings describe numerical detection only.']
    (HERE/'README.md').write_text('\n'.join(lines)+'\n')


def derive():
    OUT.mkdir(exist_ok=True);dump(OUT/'input_inventory.json',inventory())
    from coupled import gaussian_lattice,observable_stats
    trajectory_path=OUT/'coupled/trajectory_summary.json';trajectory=json.loads(trajectory_path.read_text())
    for row in trajectory:
        with np.load(OUT/f'coupled/{row["case"]}_paths.npz') as z:
            r2=np.mean(np.sum(z['beads_A']**2,axis=-1),axis=1)
            logweights=z['log_weights_joint'] if 'log_weights_joint' in z else None
            row['r2_standard_error_A2']=float(observable_stats(r2,logweights)[1])
    dump(trajectory_path,trajectory)
    dump(OUT/'coupled/model_contrasts.json',[gaussian_lattice(m,MASS['D'],300) for m in partial_cases()+strain_cases()])
    exact=compact_exact();coupled=coupled_distributions()
    consumer=json.loads((OUT/'consumer_cases.json').read_text()) if (OUT/'consumer_cases.json').exists() else consumers()
    pages(exact,coupled,consumer)


def verify():
    """Final changed-consumer evidence, without replaying production science."""
    inventory_rows=json.loads((OUT/'input_inventory.json').read_text());checks={}
    for m in [full('Pd'),full('Ni')]+partial_cases():
        source=m['site'];checks[m['name']+'_source_sumJ_error']=float(abs(m['J'][source].sum()-m['sumJ_source']))
        if not np.allclose(m['J'],m['J'].T):raise AssertionError('Nonsymmetric graph')
        if checks[m['name']+'_source_sumJ_error']>2e-6:raise AssertionError('Occupied-shell multiplicity mismatch')
    exact=json.loads((OUT/'exact/exact_summary.json').read_text());rows=[]
    for row in exact:
        if row['case'].endswith('basis8'):continue
        base=row['case'].rsplit('_',1)[0];label=row['case'].rsplit('_',1)[1]
        if label!='joint':continue
        for kind,key in [('position','radial_density_A_minus3'),('momentum','momentum_density_A3')]:
            data=np.genfromtxt(OUT/f'tables/{base}_exact_{kind}_band.csv',delimiter=',',names=True)
            ratio=data['basis8_vs10_density_delta']/np.maximum(data['density'],np.finfo(float).tiny)
            cutoff=row[kind]['density_cutoff_relative_1e6'];idx=np.argmin(abs(data[data.dtype.names[0]]-cutoff)) if cutoff else -1
            rows.append(dict(case=base,quantity=kind,cutoff_relative_1e6=cutoff,basis8_vs10_relative_density_delta_at_cutoff=float(ratio[idx]),normalization=row[kind]['normalization']))
    trajectory=json.loads((OUT/'coupled/trajectory_summary.json').read_text())
    moment_rows=json.loads((OUT/'coupled_momentum_summary.json').read_text())
    bead_comparison=[]
    for T,P in ((300,16),(77,64)):
        name=f'PdD_{T}K';fine=next(r for r in trajectory if r['case']==name);coarse=next(r for r in trajectory if r['case']==name+f'_P{P}')
        fm=next(r for r in moment_rows if r['case']==name and r['variant']=='joint')
        cm=next(r for r in moment_rows if r['case']==name+f'_P{P}' and r['variant']=='joint')
        ft=np.genfromtxt(OUT/f'tables/{name}_coupled_momentum_band.csv',delimiter=',',names=True)
        ct=np.genfromtxt(OUT/f'tables/{name}_P{P}_coupled_momentum_band.csv',delimiter=',',names=True)
        reference=next(r for r in exact if r['case']==name+'_joint')
        cutoff=reference['momentum']['radial_measure_cutoff_relative_1e6'];i=np.argmin(abs(ft['k_A_inverse']-cutoff))
        bead_comparison.append(dict(temperature_K=T,beads_coarse=P,beads_fine=2*P,
            r2_difference_A2=fine['mean_square_radius_A2']-coarse['mean_square_radius_A2'],
            conservative_sum_r2_block_SE_A2=fine['r2_standard_error_A2']+coarse['r2_standard_error_A2'],
            mean_kinetic_difference_eV=fm['mean_kinetic_eV']-cm['mean_kinetic_eV'],
            combined_kinetic_block_SE_eV=float(np.hypot(fm['kinetic_standard_error_eV'],cm['kinetic_standard_error_eV'])),
            tail_k_A_inverse=float(ft['k_A_inverse'][i]),tail_measure_fine_A=float(ft['radial_measure_A'][i]),
            tail_measure_coarse_A=float(ct['radial_measure_A'][i]),
            fine_tail_block_SE_A=float(ft['block_SE'][i]),coarse_tail_block_SE_A=float(ct['block_SE'][i]),
            tail_numerical_signal_fine=bool(ft['positive_signal'][i]),tail_numerical_signal_coarse=bool(ct['positive_signal'][i]),
            scope='same-PES bead comparison; agreement of moments is not a resolved 1e-6 coupled tail'))
    checks['PdD_bead_sensitivity']=bead_comparison
    checks['SSCHA_local_reference']=json.loads((OUT/'coupled/sscha_reference.json').read_text())
    from coupled import tracer_paths
    tracer_checks=[]
    for row in trajectory:
        if row.get('host_isotope')!='D' or row['isotope'] not in ('H','T'):continue
        name=row['case'];source_name=row['source_case']
        with np.load(OUT/f'coupled/{source_name}_paths.npz') as source,np.load(OUT/f'coupled/{name}_paths.npz') as target:
            transformed,field,weights=tracer_paths(source,full('Pd'),row['isotope'],row['temperature_K'])
            path_error=float(np.max(abs(transformed-target['beads_A'])))
            weight_error=float(np.max(abs(weights['joint']-target['log_weights_joint'])))
            source_spring=MASS['D']*np.sum(np.diff(np.concatenate([source['beads_A'],source['beads_A'][:,:1]],axis=1),axis=1)**2,axis=(1,2))
            target_spring=MASS[row['isotope']]*np.sum(np.diff(np.concatenate([target['beads_A'],target['beads_A'][:,:1]],axis=1),axis=1)**2,axis=(1,2))
            spring_error=float(np.max(abs(source_spring-target_spring)))
        moment=next(r for r in moment_rows if r['case']==name and r['variant']=='joint')
        mass_error=abs(moment['native_kinetic_factor_eV_A2']-kinetic_coefficient(MASS[row['isotope']]))
        if max(path_error,weight_error,spring_error)>1e-10 or mass_error>1e-15:raise AssertionError('Tracer mass/action contract')
        tracer_checks.append(dict(case=name,path_transform_error_A=path_error,logweight_error=weight_error,spring_mass_invariant_error_amu_A2=spring_error,native_kinetic_mass_factor_error_eV_A2=float(mass_error),frame_weight_ess=moment['frame_weight_ess'],block_weight_ess=moment['block_weight_ess']))
    checks['tracer_mass_consistency']=tracer_checks
    density_identities=[];estimator_comparisons=[]
    for row in trajectory:
        data=np.genfromtxt(OUT/f'tables/{row["case"]}_coupled_momentum_band.csv',delimiter=',',names=True)
        error=float(np.max(abs(data['radial_measure_A']-4*np.pi*data['k_A_inverse']**2*data['density_A3'])))
        if error>2e-10:raise AssertionError('Radial density versus per-radius measure')
        density_identities.append(dict(case=row['case'],max_radial_measure_identity_error_A=error))
        momentum=next(r for r in moment_rows if r['case']==row['case'] and r['variant']=='joint')
        position=np.genfromtxt(OUT/f'tables/{row["case"]}_coupled_position_band.csv',delimiter=',',names=True)
        estimator_comparisons.append(dict(case=row['case'],
            momentum_radial_normalization=float(np.trapezoid(data['radial_measure_A'],data['k_A_inverse'])),
            conditional_position_normalization=float(np.trapezoid(position['radial_measure_A_inverse'],position['r_A'])),
            conditional_r2_A2=float(np.trapezoid(position['radial_measure_A_inverse']*position['r_A']**2,position['r_A'])),
            direct_bead_r2_A2=row['mean_square_radius_A2'],
            momentum_kinetic_eV=momentum['mean_kinetic_eV'],centroid_virial_kinetic_eV=row['mean_kinetic_per_site_eV'],
            kinetic_difference_eV=momentum['mean_kinetic_eV']-row['mean_kinetic_per_site_eV'],
            conservative_sum_kinetic_block_SE_eV=momentum['kinetic_standard_error_eV']+row['kinetic_standard_error_eV'],
            scope='different estimators on the same closed trajectory; source-site momentum versus homogeneous site-average virial where stated, not independent continuum validation'))
    checks['coupled_radial_measure_identity']=density_identities
    checks['coupled_consumer_estimators']=estimator_comparisons
    checks['angular_sensitivity']=json.loads((OUT/'angular_sensitivity.json').read_text())
    consumer=json.loads((OUT/'consumer_cases.json').read_text())
    if abs(consumer['g8']['joint_probability']-1)>1e-10:raise AssertionError('Joint measure')
    if abs(consumer['g9']['normalization']-1)>.02:raise AssertionError('g9 table normalization')
    if consumer['g8']['joint_excitation_support_eV'][1]>100:raise AssertionError('Spectator support')
    checks['g8_actual_reader']=consumer['g8'];checks['g9_actual_table_import']=consumer['g9'];checks['finite_basis_cutoff_comparison']=rows
    dump(OUT/'verification.json',checks)
    print(json.dumps(dict(source_graphs=len(inventory_rows),basis_cutoff_rows=len(rows),g9_normalization=consumer['g9']['normalization'],g8_joint_probability=consumer['g8']['joint_probability']),indent=2))


def capsule():
    """Package exactly the scientific directory whitelist, without publishing."""
    import tarfile
    names=['model.py','exact.py','coupled.py','conditioned.py','coupled_momentum.py',
           'coupled_momentum_notes.md','deliver.py','README.md','replay.md','tails.png','neighbours.png']
    paths=[HERE/name for name in names]
    paths += sorted(p for p in OUT.rglob('*') if p.is_file() and p.suffix in ('.json','.npz','.csv'))
    whitelist=HERE/'ARCHIVE_FILES.txt';checksum=HERE/'SHA256SUMS'
    paths += [whitelist,checksum]
    relative=sorted(str(p.relative_to(ROOT)) for p in paths)
    whitelist.write_text('\n'.join(relative)+'\n')
    checksum.write_text(''.join(f'{hashlib.sha256((ROOT/p).read_bytes()).hexdigest()}  {p}\n' for p in relative if ROOT/p!=checksum))
    archive=HERE/'g10_bulk.tar.gz'
    with tarfile.open(archive,'w:gz',compresslevel=6) as tar:
        for path in relative:
            source=ROOT/path;info=tar.gettarinfo(str(source),arcname=path)
            info.uid=info.gid=0;info.uname=info.gname='';info.mtime=0
            with source.open('rb') as stream:tar.addfile(info,stream)
    print(json.dumps(dict(archive=str(archive.relative_to(ROOT)),sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),whitelist=str(whitelist.relative_to(ROOT)),files=len(relative))))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('stage',choices=['exact','coupled','derive','verify','capsule']);args=parser.parse_args()
    if args.stage=='exact':
        from exact import run
        run(OUT/'exact')
    elif args.stage=='coupled':
        from coupled import run
        run(OUT/'coupled')
    else:globals()[args.stage]()


if __name__=='__main__':main()
