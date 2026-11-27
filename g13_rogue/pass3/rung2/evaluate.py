"""Entry_4166 additions; all beam dynamics remain in the pinned signed producer.
Run from repository root. This file never writes into rung1b or uploads anything.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,sys,tarfile
import numpy as np
from scipy.constants import Boltzmann as KB
from scipy.optimize import brentq
from scipy.special import gammainc
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'rung1b'))
import signed as s
BASE='650a2ebd57c35f371ecd23c4db975ebd6e3cc55a'
DRIVER='f6e408a107fd24c767eaba4ab67983f7d3439a76d718e480f6b397dcc3c4e891'
FREQUENCIES='0b4661cdaed7c2c84e2af7bbc82507edd950dc2b139673454a725be77631381c'
CORE=s.BEAMS[:3]
T=300.;KT=KB*T

def pinned():
    if s.sha(s.HERE/'signed.py')!=DRIVER or s.sha(s.HERE/'frequencies.json')!=FREQUENCIES:
        raise AssertionError('frozen signed dependency changed')

def read(path):
    with path.open() as f:return list(csv.DictReader(f))

def built(beam,loss):
    """Zero random load operator, not a copied or edited time integrator.

    Frozen auxiliary mark strata still allocate fine arrival positions. They
    carry no force/work here; their diagnostic is explicitly relabelled below.
    Beta is inactive, so one propagation serves both beta labels.
    """
    pinned();records=list(range(16));model=s.load_model(2);model=dict(model,random=np.zeros_like(model['random']))
    original=s.load_model
    try:
        s.load_model=lambda rung:model if rung==2 else original(rung)
        m,z=s.run_arm(2,beam,0.,loss,True,records)
    finally:s.load_model=original
    if np.any(z['x'][:,:,0]) or np.any(z['y'][:,:,0]) or np.any(z['work'][:,0]):
        raise AssertionError('built-only injected random launch')
    names=s.BEAMS
    try:
        s.BEAMS=tuple(name+'__built_only' for name in names)
        s.read_arm(2,beam,0.,loss,True,records,m,z)
    finally:s.BEAMS=names
    stem=f'r2_{names[beam]}__built_only_b0.0_{loss}_combined_records0-15'
    # Preserve, rather than disguise, the frozen producer's unused-stratum
    # diagnostic. It is not a physically launched energy budget in this case.
    path=HERE/(stem+'_energy.csv');rows=read(path)
    for row in rows:
        row['latent_sampling_stratum_self_J']=row.pop('marked_self_J')
        row['latent_sampling_stratum_host_fraction']=row.pop('marked_host_fraction')
        row['physical_random_launch_budget_J']=0.
        row['source_case']='built_only; beta inactive, shared labels 0/.35'
    s.table(path,rows)
    totalx=z['x'].sum(axis=2);totaly=z['y'].sum(axis=2)
    energy=.5*np.sum(totalx*totalx+totaly*totaly,axis=2)
    s.table(HERE/(stem+'_band_energy.csv'),[
        dict(record=record,source_case='built_only',beta='inactive: shared 0/.35',mean_dynamic_band_J=float(energy[:,i].mean()),peak_dynamic_band_J=float(energy[:,i].max()),thermal_300K_band_mean_J=len(m['omega'])*KT,mean_over_thermal_mean=float(energy[:,i].mean()/(len(m['omega'])*KT)),peak_over_thermal_mean=float(energy[:,i].max()/(len(m['omega'])*KT))) for i,record in enumerate(records)])

def thermal():
    """Stationary classical modal equilibrium; not all-band lattice heat.

    Exact Gaussian moments of the same spatial quadratics. Finite-window mean
    variances use the viscous oscillator two-time covariance at the actual DT.
    No extra stochastic forcing is inserted into any beam trajectory.
    """
    pinned();m=s.load_model(2);w=m['omega'];n=s.SAMPLE;lag=np.arange(n)*s.DT;rows=[];projections=[]
    for loss in ('retained','Cu_proxy'):
        gamma=1/50e-6+m['speed']/.2 if loss=='retained' else 4460*np.log(10)/10*32*w/(2*np.pi*1e6)
        wd=np.sqrt(w*w-gamma*gamma/4)
        d=np.exp(-lag[:,None]*gamma/2)*(np.cos(lag[:,None]*wd)+gamma/(2*wd)*np.sin(lag[:,None]*wd))
        weights=np.r_[n,2*(n-np.arange(1,n))]/n**2
        for region in range(3):
            a=m['regional'][region,0];a=(a+a.T)/2
            cov=.5*KT**2*np.sum((d@(a*a))*d,axis=1)
            varmean=float(weights@cov)
            rows.append(dict(loss=loss,region=region,depth_m=(30e-6,50e-6,200e-6)[region],thermal_mean_elastic_J_m3=float(.5*KT*np.trace(a)),thermal_instant_sd_J_m3=float(KT*np.sqrt(np.sum(a*a)/2)),thermal_window_mean_sd_J_m3=float(np.sqrt(varmean)),samples=n,dt_s=s.DT,exposure_s=n*s.DT,convention='stationary equilibrium in the actual retained modes; same regional elastic quadratic, not a D vibration threshold'))
        for pi,bound in enumerate(m['region_bounds']):
            for oi,name in enumerate(('epsilon_xx','u_z')):
                p=m['region_projectors'][pi,oi];variance=KT*np.sum(p*p)
                meanvar=float(KT*(weights@d)@(p*p))
                if meanvar<0:raise AssertionError('negative thermal projected mean variance')
                projections.append(dict(loss=loss,projection=pi,observable=name,side_m=float(bound[0]),z0_m=float(bound[1]),z1_m=float(bound[2]),mean=0.,instant_rms=float(np.sqrt(variance)),window_mean_rms=float(np.sqrt(meanvar)),exposure_s=n*s.DT))
    s.table(HERE/'thermal_300K_regions.csv',rows);s.table(HERE/'thermal_300K_projections.csv',projections)
    s.dump(HERE/'thermal_300K.json',dict(temperature_K=T,frequency_Hz=[float(w.min()/(2*np.pi)),float(w.max()/(2*np.pi))],modes=len(w),mean_total_band_J=len(w)*KT,sd_total_band_J=float(np.sqrt(len(w))*KT),mean_elastic_band_J=.5*len(w)*KT,sd_elastic_band_J=float(np.sqrt(len(w)/2)*KT),occupation='classical mass-normalized q/v Gaussian equilibrium; hbar omega/kBT below 4.1e-8',time_convention='stationary ambient reference, 2000 samples at 50 ns over 100 us; no beam-onset transient or extra bath added to source cases',limits='thermal mean and fluctuations, not an actually sampled thermal maximum/own-mean distribution, all-band heat, or a local D energy threshold'))

def local_thermal():
    """Pointwise reference matches the three same-event center/depth operators."""
    pinned();m=s.load_model(2);w=m['omega'];n=s.SAMPLE;lag=np.arange(n)*s.DT;rows=[]
    weights=np.r_[n,2*(n-np.arange(1,n))]/n**2
    global_rows=[]
    for loss in ('retained','Cu_proxy'):
        gamma=1/50e-6+m['speed']/.2 if loss=='retained' else 4460*np.log(10)/10*32*w/(2*np.pi*1e6)
        wd=np.sqrt(w*w-gamma*gamma/4)
        d=np.exp(-lag[:,None]*gamma/2)*(np.cos(lag[:,None]*wd)+gamma/(2*wd)*np.sin(lag[:,None]*wd))
        total_cov=KT**2*np.sum(np.exp(-lag[:,None]*gamma)*(1+gamma**2/(2*wd**2)*np.sin(lag[:,None]*wd)**2),axis=1)
        global_rows.append(dict(loss=loss,modes=len(w),thermal_mean_total_band_J=len(w)*KT,thermal_instant_sd_total_band_J=float(np.sqrt(len(w))*KT),thermal_window_mean_sd_total_band_J=float(np.sqrt(weights@total_cov)),samples=n,dt_s=s.DT,convention='same global dynamic modal-energy quadratic; stationary FDT reference, not all-band heat'))
        for di,depth in enumerate((30e-6,50e-6,200e-6)):
            pi=9*di+4;strain=m['strain'][pi];a=strain.T@m['cs'][pi]@strain;a=(a+a.T)/2
            cov=.5*KT**2*np.sum((d@(a*a))*d,axis=1)
            rows.append(dict(loss=loss,depth_m=depth,point_index=pi,thermal_mean_elastic_J_m3=float(.5*KT*np.trace(a)),thermal_instant_sd_J_m3=float(KT*np.sqrt(np.sum(a*a)/2)),thermal_window_mean_sd_J_m3=float(np.sqrt(weights@cov)),samples=n,dt_s=s.DT,convention='same fixed center-point strain quadratic; not region average or a local D threshold'))
    s.table(HERE/'thermal_300K_local.csv',rows)
    s.table(HERE/'thermal_300K_global.csv',global_rows)

def spectrum():
    """Actual current extent plus a priced, unexecuted, fixed-budget addition."""
    pinned();c,*_=s.old.accepted();sp=s._space();m=s.load_model(2);metric=np.load(HERE/'mode_metrics.npz')
    selected=metric['vertical_mass']>.5;wavelength=2*np.pi/metric['vertical_k_per_m'][selected]
    edges=[]
    for theta in (0.,np.pi/4):
        roots={}
        for hz in (250000.,475000.,625000.):
            roots[hz]=brentq(lambda lam:sp.plane_branches(c,lam,theta)[0]-hz,.0011,.008,xtol=1e-9)
        edges.append(dict(direction_rad=theta,Q4_uniform_A0_wavelength_at_250kHz_m=roots[250000.],candidate_low_wavelength_m=roots[475000.],candidate_high_wavelength_m=roots[625000.]))
    # Same illustrative 3-D Gaussian spectrum used by accepted pulse_sensitivities:
    # f^2 exp[-(2 pi f hypot(tau,a/v))^2]. No renormalization into the new band.
    a=20e-9;tau=.5e-12;v=4460.;time=np.hypot(tau,a/v)
    frac=lambda lo,hi:float(gammainc(1.5,(2*np.pi*hi*time)**2)-gammainc(1.5,(2*np.pi*lo*time)**2))
    added=.99*frac(475000.,625000.)/(1-float(gammainc(1.5,(2*np.pi*250000*time)**2)))
    area=float(np.prod(c.DIMS[:2]));estimate=area/(4*np.pi)*np.mean([(2*np.pi/r['candidate_high_wavelength_m'])**2-(2*np.pi/r['candidate_low_wavelength_m'])**2 for r in edges])
    # Include the two low-order in-plane branches in this indicative Weyl price.
    velocity=sp.plane_branches(c,.01,0.)[1:]*.01
    estimate+=area*np.pi*(625000.**2-475000.**2)*np.sum(1/velocity**2)
    new=int(np.ceil(estimate));total=len(m['omega'])+new
    s.table(HERE/'actual_band_extent.csv',edges)
    s.dump(HERE/'band_addition_price.json',dict(status='priced proposal only; not propagated or launched',current_frequency_Hz=[float(m['omega'].min()/(2*np.pi)),float(m['omega'].max()/(2*np.pi))],current_modes=len(m['omega']),vertical_mass_above_half_modes=int(selected.sum()),out_of_plane_component_RMS_lateral_wavelength_m=[float(wavelength.min()),float(wavelength.max())],wavelength_scope='actual grained-mode component RMS, not a unique Bloch wavelength or hard spatial cutoff',projection_screen_m=.01373147/9,target_uniform_reference_A0_Hz=[491424.27691095,588157.570569],candidate_disjoint_Hz=[475000.,625000.],fixed_source=dict(spatial_sigma_m=a,duration_sigma_s=tau,bulk_proxy_speed_m_s=v,law='normalized f^2 Gaussian recoil envelope; illustrative, not measured launch',original_5_250kHz_host_fraction=.01,extra_host_fraction=added,extra_nominal_power_W=added*.087123,remaining_energy='outside the old band, never rescaled into the candidate band'),estimated_added_modes=new,estimated_total_modes=total,additional_eigenvector_bytes=int(63375*new*8),modal_step_cost_factor=total/len(m['omega']),dense_regional_readout_cost_factor=(total/len(m['omega']))**2,price_limits='two-direction uniform-reference Weyl estimate, not an actual grained spectrum or runtime benchmark; requires a bracketed new spectral window and a band-edge phase check. Existing Q4 errors near the upper new edge are not a general accuracy guarantee. No job or later transport model is authorized by this estimate.'))

def refine_price():
    """A bounded homogeneous branch check prices, but does not build, an upgrade."""
    c,*_=s.old.accepted();sp=s._space();price=json.loads((HERE/'band_addition_price.json').read_text());rows=[];counts=[]
    for theta in (0.,np.pi/4):
        roots={}
        for branch in range(3):
            for hz in (250000.,475000.,625000.):
                roots[branch,hz]=brentq(lambda lam:sp.uniform_reference(c,lam,theta,8)[0][branch]-hz,.0008,.05,xtol=1e-9)
        counts.append(float(np.prod(c.DIMS[:2])/(4*np.pi)*sum((2*np.pi/roots[b,625000.])**2-(2*np.pi/roots[b,475000.])**2 for b in range(3))))
        lam=roots[0,625000.]
        for degree in (4,5):
            hz=float(sp.plane_branches(c,lam,theta,degree=degree,cells=16)[0])
            rows.append(dict(direction_rad=theta,lateral_degree=degree,wavelength_m=lam,reference_P8_frequency_Hz=625000.,numerical_frequency_Hz=hz,relative_phase_error=hz/625000.-1,scope='homogeneous A0-like upper-band-edge witness, not a general grained-foil bound'))
        rows.append(dict(direction_rad=theta,lateral_degree='continuum P8',wavelength_m=roots[0,250000.],reference_P8_frequency_Hz=250000.,scope='current band-edge uniform A0-like wavelength'))
    price.pop('additional_eigenvector_bytes',None)
    new=int(np.ceil(np.mean(counts)));total=303+new;ndof=3*81*81*5;nnz=(16*36-15)**2*225
    price.update(estimated_added_modes=new,estimated_total_modes=total,priced_global_space=dict(lateral_degree=5,lateral_cells=16,thickness_degree=4,ndof=ndof,estimated_stiffness_entries=nnz,minimum_raw_CSC_bytes=12*nnz+8*(ndof+1),full_modal_vector_bytes=ndof*total*8,projection='Q4 subspace is nested in Q5 at the same cells: existing mass capture cannot decrease in exact arithmetic'),modal_step_cost_factor=total/303,dense_regional_readout_cost_factor=(total/303)**2,time_price='At 625 kHz a 50-ns left-bin phase scale is 0.098 rad. Validate actual averaged-thermal-force transfer; conservative 10-ns nested parent bins would cost 5x steps, or use a verified higher-order/exact thermal convolution. Neither is implemented here.',source_budget_note='Old 5–250-kHz 1% packet remains unchanged as a separate component; candidate Gaussian share comes only from the remaining host reservoir, with no electronic stopping reclassified as extra acoustic launch. Built-only keeps zero random launch in all bands.',resource_limits='Weyl counts and raw storage are estimates, not measured runtime/peak factor memory. A new bracketed global spectral window and controlled band-edge phase are required; a larger direct factorization is not promised affordable. No dense production or later transport model has been run.')
    s.table(HERE/'band_upgrade_reference.csv',rows);s.dump(HERE/'band_addition_price.json',price)

def finish():
    pinned();local_thermal();refine_price();thermal_rows=read(HERE/'thermal_300K_regions.csv');eq={(r['loss'],int(r['region'])):r for r in thermal_rows}
    sets=[];summary=[];contrasts=[];comparisons=[];widths=[];missing=[]
    for beam,name in enumerate(CORE):
        for loss in ('retained','Cu_proxy'):
            for beta,source in [(b,src) for b in (0.,.35) for src in ('random','combined')]+[(0.,'built_only')]:
                arm=f'r2_{name}'+('__built_only' if source=='built_only' else '')+f'_b{beta}_{loss}_'+('combined' if source=='built_only' else source)
                stem=HERE/(arm+'_records0-15')
                paths={kind:Path(str(stem)+'_'+kind+'.csv') for kind in ('records','peaks','tails','energy')}
                if not all(p.exists() for p in paths.values()):missing.append(arm);continue
                rr=read(paths['records']);pk=read(paths['peaks']);ta=read(paths['tails']);en=read(paths['energy'])
                if sorted(int(r['record']) for r in rr)!=list(range(16)):raise AssertionError('record coverage')
                entry=dict(beam=name,loss=loss,beta=('inactive/shared' if source=='built_only' else beta),source=source,arm=arm,records=rr,peaks=pk,tails=ta,energy=en);sets.append(entry)
                for metric in ('peak_over_own_mean','peak_J_m3','duration_s','random_diagonal_J_m3','momentum_diagonal_J_m3','thermal_diagonal_J_m3','all_component_cross_J_m3','signed_diagonal_fraction'):
                    summary.append(dict(beam=name,loss=loss,beta=entry['beta'],source=source,region='leading event',observable=metric,**s.stats([r[metric] for r in pk])))
                measured=[r for r in pk if r['width_x_m']!='']
                width=dict(beam=name,loss=loss,beta=entry['beta'],source=source,records=16,width_records=len(measured),duration_censored_records=sum(r['duration_censored']=='True' for r in pk),width_foil_censored_records=sum(r['width_foil_censored']=='True' for r in measured),width_band_censored_records=sum(r['width_band_censored']=='True' for r in measured),width_unresolved_records=sum(r['width_map_unresolved']=='True' for r in measured),scope='widths only preselected records 0/1; no 16-record width precision claim')
                for axis in ('x','y'):
                    values=[float(r[f'width_{axis}_m']) for r in measured]
                    width.update({f'width_{axis}_mean_m':float(np.mean(values)),f'width_{axis}_min_m':min(values),f'width_{axis}_max_m':max(values)})
                widths.append(width)
                for key in (key for key in rr[0] if key.endswith('_fraction') or key=='B_initiation_rate_s'):
                    summary.append(dict(beam=name,loss=loss,beta=entry['beta'],source=source,observable=key,**s.stats([r[key] for r in rr])))
                for key in (key for key in rr[0] if key.endswith('_parents')):
                    summary.append(dict(beam=name,loss=loss,beta=entry['beta'],source=source,observable=key[:-8]+'_parent_strike_rate_s',**s.stats([float(r[key])/float(r['exposure_s']) for r in rr])))
                for region in range(3):
                    rt=[r for r in ta if int(r['region'])==region];reference=eq[loss,region]
                    for metric in ('mean_J_m3','max_J_m3','peak_over_own_mean','elastic_eV_D','conditional_eV_D'):
                        summary.append(dict(beam=name,loss=loss,beta=entry['beta'],source=source,region=region,observable=metric,**s.stats([r[metric] for r in rt])))
                    for metric in ('mean_J_m3','max_J_m3'):
                        comparisons.append(dict(beam=name,loss=loss,beta=entry['beta'],source=source,region=region,observable=metric+'_over_equilibrium_mean',thermal_mean_J_m3=reference['thermal_mean_elastic_J_m3'],thermal_instant_sd_J_m3=reference['thermal_instant_sd_J_m3'],thermal_window_mean_sd_J_m3=reference['thermal_window_mean_sd_J_m3'],**s.stats([float(r[metric])/float(reference['thermal_mean_elastic_J_m3']) for r in rt])))
    if missing:raise AssertionError(('unfinished core source cases',missing))
    for current in sets:
        if current['source']=='built_only':continue
        controls=[]
        if current['source']=='combined':controls.append(next(r for r in sets if r['beam']==current['beam'] and r['loss']==current['loss'] and r['beta']==current['beta'] and r['source']=='random'))
        controls.append(next(r for r in sets if r['beam']==current['beam'] and r['loss']==current['loss'] and r['source']=='built_only'))
        for control in controls:
            for region in ('leading event',0,1,2):
                aa=current['peaks'] if region=='leading event' else [r for r in current['tails'] if int(r['region'])==region]
                bb=control['peaks'] if region=='leading event' else [r for r in control['tails'] if int(r['region'])==region]
                for metric in (('peak_over_own_mean',) if region=='leading event' else ('mean_J_m3','max_J_m3','peak_over_own_mean')):
                    a={int(r['record']):float(r[metric]) for r in aa};b={int(r['record']):float(r[metric]) for r in bb}
                    contrasts.append(dict(beam=current['beam'],loss=current['loss'],beta=current['beta'],source=current['source'],control=control['source'],region=region,observable=metric,**s.stats([a[i]/b[i] for i in range(16)])))
        if current['source']=='combined':
            path=s.HERE/(current['arm'].replace('r2_','r1_',1)+'_records0-15_peaks.csv');prior=read(path)
            a={int(r['record']):float(r['peak_over_own_mean']) for r in current['peaks']};b={int(r['record']):float(r['peak_over_own_mean']) for r in prior}
            contrasts.append(dict(beam=current['beam'],loss=current['loss'],beta=current['beta'],source='combined',control='same-sign rung1',region='leading event',observable='peak_over_own_mean',**s.stats([a[i]/b[i] for i in range(16)])))
    counts={r['record']:r['all_parents'] for r in sets[0]['records']}
    if any(any(counts[r['record']]!=r['all_parents'] for r in case['records']) for case in sets):raise AssertionError('parent mismatch')
    built_cases=[r for r in sets if r['source']=='built_only']
    if any(any(float(e['random_work_J'])!=0 or float(e['physical_random_launch_budget_J'])!=0 for e in case['energy']) for case in built_cases):raise AssertionError('nonzero built random work')
    s.table(HERE/'summary_sources_r2.csv',summary);s.table(HERE/'contrasts_sources_r2.csv',contrasts);s.table(HERE/'thermal_level_comparisons.csv',comparisons)
    s.table(HERE/'width_censoring_sources_r2.csv',widths)
    s.dump(HERE/'verification_sources_r2.json',dict(delivery_base=BASE,consumed_driver_sha256=DRIVER,primary_beams=list(CORE),physical_source_arms=len(sets),logical_beta_source_arms=36,independent_records=16,raw_core_record_rows=16*len(sets),built_only_beta='inactive; one shared propagation, never pooled as independent beta replicates',same_parent_totals=True,built_random_work_and_budget_zero=True,max_relative_energy_residual=max(abs(float(e['balance_residual_J']))/float(e['incident_J']) for case in sets for e in case['energy']),thermal_reference='same actual modes, spatial quadratics/projectors and 2000x50-ns window; analytic stationary mean/fluctuations, not sampled thermal maxima',source_scenarios='not rigorous physical lower/upper bounds',ripple_status='previously authorized results retained as adjuncts only; no built-only ripple propagation'))
    # Compare already recorded clip marks to completed paired core records; no replay.
    clips_path=HERE/'clips_r2.json';clips=json.loads(clips_path.read_text())
    for clip in clips:
        loss='Cu_proxy' if 'Cu_proxy' in clip['path'] else 'retained'
        case=next(r for r in sets if r['beam']=='hot_bursts' and r['loss']==loss and r['beta']==.35 and r['source']=='combined')
        clip['ensemble_mark_digest_equal']=next(r for r in case['records'] if r['record']=='0')['draw_sha256']==clip['draw_sha256']
    s.dump(clips_path,clips)

def archive(path,names,manifest_name,manifest):
    s.dump(HERE/manifest_name,manifest)
    with path.open('wb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as gz,tarfile.open(fileobj=gz,mode='w|') as tf:
        for name in names+[manifest_name]:
            source=HERE/name;info=tarfile.TarInfo(str(source.relative_to(ROOT)));info.size=source.stat().st_size;info.mode=0o644;info.mtime=0
            with source.open('rb') as payload:tf.addfile(info,payload)
    expected={row['path']:row for row in manifest['files']};seen=False
    with tarfile.open(path,'r|gz') as tf:
        for member in tf:
            if not member.isfile():raise AssertionError('non-regular member')
            with tf.extractfile(member) as payload:
                if member.name==str((HERE/manifest_name).relative_to(ROOT)):
                    if json.load(payload)!=manifest:raise AssertionError('manifest mismatch')
                    seen=True
                else:
                    row=expected.pop(member.name)
                    if member.size!=row['bytes'] or hashlib.file_digest(payload,'sha256').hexdigest()!=row['sha256']:raise AssertionError('payload mismatch')
    if expected or not seen:raise AssertionError('incomplete archive')
    return dict(path=str(path.relative_to(ROOT)),sha256=s.sha(path),bytes=path.stat().st_size,members=len(names)+1,all_member_hashes_read_back=True)

def package():
    pinned();v=json.loads((HERE/'verification_sources_r2.json').read_text())
    if v['physical_source_arms']!=30:raise AssertionError('incomplete primary source matrix')
    excluded={'manifest.json','evidence_manifest.json','capsule.tar.gz','histories.tar.gz','delivery.json'}
    names=sorted(p.name for p in HERE.iterdir() if p.is_file() and p.name not in excluded)
    evidence=[name for name in names if name.endswith('.npz') or (name.startswith('clip_') and name.endswith('.png') and not name.endswith('_preview.png'))]
    main=[name for name in names if name not in evidence]
    def entries(selected):return [dict(path=str((HERE/name).relative_to(ROOT)),bytes=(HERE/name).stat().st_size,sha256=s.sha(HERE/name)) for name in selected]
    dependency=dict(delivery_base=BASE,input_epoch_base=s.BASE,required_shared_files=[dict(path='g13_rogue/pass3/rung1b/signed.py',sha256=DRIVER),dict(path='g13_rogue/pass3/rung1b/frequencies.json',sha256=FREQUENCIES)],meaning='installed accepted driver is consumed unchanged; no full rung1b history archive is needed')
    evidence_manifest=dict(unpublished=True,stage='local readback only',dependency=dependency,files=entries(evidence))
    data=archive(HERE/'histories.tar.gz',evidence,'evidence_manifest.json',evidence_manifest)
    main=sorted(main+['evidence_manifest.json'])
    manifest=dict(unpublished=True,stage='local readback only; remote stage owned by Main',dependency=dependency,required_evidence_asset=data,evidence_manifest_sha256=s.sha(HERE/'evidence_manifest.json'),files=entries(main))
    result=archive(HERE/'capsule.tar.gz',main,'manifest.json',manifest)
    if result['bytes']>=50*1024**2:raise AssertionError('main capsule exceeds 50 MiB')
    receipt=dict(stage='local archive authentication, not remote GET',main=result,evidence=data,dependency=dependency)
    s.dump(HERE/'delivery.json',receipt);print(json.dumps(receipt),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['built','thermal','spectrum','finish','package']);parser.add_argument('--beam',type=int);parser.add_argument('--loss',choices=['retained','Cu_proxy']);a=parser.parse_args()
    if a.action=='built':built(a.beam,a.loss)
    else:globals()[a.action]()
if __name__=='__main__':main()
