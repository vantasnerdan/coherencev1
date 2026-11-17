"""All six ALD families through the shared physical merit/compare consumer.

Kinetically unknown material properties remain explicit scenario endpoints.
Seals use actual transport.Model; cap retention is not a manufactured rate.
"""
from pathlib import Path
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
import csv
import gzip
import json
import os
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from g4_beam.transport import transport as T, bca
from g4_beam.merit.api import Target, reference, uniform_reference, compare, evaluate, ENERGIES, FOUR, CAUSES, write_csv
from g4_beam.merit import api as merit_api
from g4_beam.forward.profile_thermal import transport_stopping
HERE=Path(__file__).resolve().parent
MATERIALS={r['id']:r for r in json.loads((HERE/'inputs/seal-cap.json').read_text())['materials']}
WORKERS=int(os.environ.get('G5_WORKERS','4'))
if not 1<=WORKERS<=4:raise ValueError('G5_WORKERS must be1–4')
OTHER=json.loads((HERE/'inputs/other-designs.json').read_text())
NS=1.53e19


@dataclass(frozen=True)
class Stopping:
    # entries=(effectiveZ,effectiveM,((Z,M,host_atomic_fraction),...))
    entries: tuple=()
    scale: float=1.
    def __call__(self,e,layer):
        _,_,nh,nd,z,m,_=layer
        species=next((s for zz,mm,s in self.entries if abs(z-zz)<1e-7 and abs(m-mm)<1e-6),None)
        if species is None:
            host=(bca.se_ls(e*1000,z)+bca.sn_zbl(e*1000,z,m))
        else:
            host=sum(f*(bca.se_ls(e*1000,zz)+bca.sn_zbl(e*1000,zz,mm)) for zz,mm,f in species)
        # Scale only non-Pd material stopping; the common Pd substrate stays matched.
        factor=self.scale if z!=46 else 1.
        return float((nh*host*factor+nd*(bca.se_ls(e*1000,1)+bca.sn_zbl(e*1000,1,2.0141)))*1e-18)


def specs(material):
    species=[s for s in material['stopping']['species'] if not (s['Z']==1 and s['M_amu']<3)]
    count=sum(s['stoich'] for s in species)
    z=sum(s['Z']*s['stoich'] for s in species)/count
    mass=sum(s['M_amu']*s['stoich'] for s in species)/count
    fractions=tuple((s['Z'],s['M_amu'],s['stoich']/count) for s in species)
    return z,mass,fractions,count


def sputter(material,e):
    row=material['sputter']
    if 'E_lab_keV' not in row:
        row=MATERIALS['TiD2_cap' if material['id'].endswith('_cap') else 'Al2O3_seal']['sputter']
    x=row['E_lab_keV'];y=row['Y_atoms_per_D']
    return (.85*np.interp(2*e,x,y)+.30*np.interp(e,x,y))/1.15


def erosion(material,e,thickness_nm):
    den=material['density']
    nh=den['metal_atoms_m3'] if material['id'].endswith('_cap') else den['host_atoms_m3']
    rate=T.beam_flux(2*e)*sputter(material,e)/nh # host atoms sputtered / host atom density
    return rate,thickness_nm*1e-9/rate


def compact_profile(profile):
    layers,nads,state=profile
    return dict(layers=np.asarray(layers).tolist(),nads_m2=float(nads),
                state={k:np.asarray(v).tolist() for k,v in state.items()})


def profile_tuple(row):return row['layers'],row['nads_m2'],row.get('state',{})


def observed_run(model,times):
    """Observe a slow case's actual RHS without changing equations or solver."""
    original=model.rhs;last=-np.inf;start=time.monotonic()
    def rhs(t,y,beam_on=True):
        nonlocal last
        now=time.monotonic()
        if now-last>=30.:
            last=now
            print('MODEL_PROGRESS','E_cm',model.c['E_cm'],'beam_on',beam_on,'solver_t_s',t,
                  'wall_s',now-start,'front_D_per_metal',float(y[0]),
                  'max_solute_D_per_metal',float(np.max(y[:model.N])),
                  'film_nm',float(y[3*model.N]*1e9),flush=True)
        return original(t,y,beam_on)
    model.rhs=rhs
    try:
        return model.run_te(times)
    finally:
        del model.rhs


def seal_model(job):
    material_id,dnm,permscale,e,front_s,jrel=job
    mat=MATERIALS[material_id];z,m,_,count=specs(mat)
    key='g5_'+material_id
    old=T.OVERLAYERS['Al2O3']
    T.OVERLAYERS[key]=dict(old,n=mat['density']['formula_units_m3'],natoms=count,nmet=2,Z=z,M=m,
                          perm=mat['permeability']['central']*permscale,sp='Al2O3')
    case=dict(material='Pd',loaded=True,E_cm=e,overlayer=key,d_ol=dnm*1e-9,j_rel=jrel)
    if front_s is not None:case['s_front']=front_s
    model=T.Model(case)
    # Use sourced host inventory and sputter proxy, not Al2O3 density for Er2O3.
    model.Ysp_ol=sputter(mat,e)
    rate,life=erosion(mat,e,dnm)
    return model,rate,life


def transport_job(job):
    material_id,dnm,permscale,e,front_s,jrel=job
    model,rate,life=seal_model(job)
    t95=(.95*dnm*1e-9+.02e-9*np.log(20))/rate
    times=np.unique(np.clip([0,600,1800,min(life/2,model.c['t_run']),min(t95,model.c['t_run']),model.c['t_run']],0,model.c['t_run']))
    start=time.monotonic();model.run_te(times)
    snapshots=[]
    for i,t in enumerate(times):
        row=compact_profile(model.profile(i));row.update(time_s=float(t),T_K=model.T)
        snapshots.append(row)
    return dict(material=material_id,thickness_nm=dnm,permeability_scale=permscale,E_cm_keV=e,
                lifetime_nominal_s=life,lifetime_95_s=t95,erosion_rate_m_s=rate,
                snapshots=snapshots,wall_s=time.monotonic()-start)


def _surface_case(job):
    name,e,s,j=job
    case=dict(material='Pd',loaded=True,E_cm=e,s_front=s,j_rel=j if j>0 else 1.)
    if j==0:case['eta_abs']=0.
    model=T.Model(case)
    observed_run(model,np.array([0.,600.,1800.,model.c['t_run']]))
    return name,e,compact_profile(model.profile()),float(model.T)


def baseline_history(job):
    e,times=job
    model=T.Model(dict(material='Pd',loaded=True,E_cm=e))
    model.run_te(np.asarray(sorted(set(times+[model.c['t_run']]))))
    return e,{str(float(t)):compact_profile(model.profile(i)) for i,t in enumerate(model.t_beam)}


def matched_histories(designs,stage='seal'):
    path=HERE/('matched_reference_history.json.gz' if stage=='seal' else 'matched_'+stage+'_reference_history.json.gz')
    if path.exists():return json.loads(gzip.decompress(path.read_bytes()))
    jobs=[]
    for e in FOUR:
        times=sorted({float(s['time_s']) for d in designs if d.get('erosion_snapshots') for s in d['erosion_snapshots'][str(float(e))]})
        jobs.append((float(e),times))
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:result={str(e):rows for e,rows in pool.map(baseline_history,jobs)}
    path.write_bytes(gzip.compress(json.dumps(result,separators=(',',':')).encode()))
    return result


def promoter_material(name):
    definitions={'BaO':((56,137.327,1),(8,15.999,1)),
                 'SrO':((38,87.62,1),(8,15.999,1)),
                 'La2O3':((57,138.905,2),(8,15.999,3))}
    species=[dict(Z=z,M_amu=m,stoich=n) for z,m,n in definitions[name]]
    nmetal=2 if name=='La2O3' else 1
    # Explicit one-Pd-site/metal inventory and0.3nm island-height scenario.
    nf=NS/nmetal/.3e-9;count=sum(s['stoich'] for s in species)
    return dict(id=name+'_seal',stopping={'species':species},
                density={'formula_units_m3':nf,'metal_atoms_m3':nf*nmetal,'host_atoms_m3':nf*count},
                sputter={'proxy':'Al2O3'})


def promoter_job(job):
    name,s,e=job;mat=promoter_material(name);z,m,_,count=specs(mat)
    key='g5_'+name
    T.OVERLAYERS[key]=dict(T.OVERLAYERS['Al2O3'],n=mat['density']['formula_units_m3'],
                          natoms=count,nmet=2 if name=='La2O3' else 1,Z=z,M=m,
                          xm=0.,perm=np.inf,s=s,sticking_is_patch_value=True,Eads=T.PD["Eads"])
    model=T.Model(dict(material='Pd',loaded=True,E_cm=e,overlayer=key,d_ol=.3e-9))
    model.Ysp_ol=sputter(mat,e);rate,life=erosion(mat,e,.3)
    times=np.unique(np.clip([0,life/2,life,model.c['t_run']],0,model.c['t_run']))
    observed_run(model,times)
    return dict(material=name,s=s,E_cm_keV=e,lifetime_s=life,
                snapshots=[dict(compact_profile(model.profile(i)),time_s=float(t)) for i,t in enumerate(times)])


def bzy_material(dpf):
    source=next(r for r in OTHER['rows'] if r['id']=='BZY20_back_and_front')
    nf=source['density']['formula_units_m3']
    species=[dict(s,stoich=3. if s['element']=='O' and dpf>0 else s['stoich']) for s in source['stopping_species']]
    return dict(id='BZY_seal',stopping={'species':species},density={'formula_units_m3':nf,
                'metal_atoms_m3':2*nf,'host_atoms_m3':nf*sum(s['stoich'] for s in species)},sputter={'proxy':'Al2O3'})


def bzy_job(job):
    d,perm,e=job;mat=bzy_material(0.);z,m,_,count=specs(mat);key='g5_BZY'
    T.OVERLAYERS[key]=dict(T.OVERLAYERS['Al2O3'],n=mat['density']['formula_units_m3'],
                          natoms=count,nmet=2,Z=z,M=m,xm=0.,perm=3e11*perm)
    model=T.Model(dict(material='Pd',loaded=True,E_cm=e,overlayer=key,d_ol=d*1e-9))
    model.Ysp_ol=sputter(mat,e);rate,life=erosion(mat,e,d)
    times=np.unique(np.clip([0,600.,min(life/2,model.c['t_run']),model.c['t_run']],0,model.c['t_run']))
    observed_run(model,times)
    return dict(d=d,perm=perm,E_cm_keV=e,lifetime_s=life,
                snapshots=[dict(compact_profile(model.profile(i)),time_s=float(t)) for i,t in enumerate(times)])


def _logged_remaining_job(payload):
    stage,job=payload;start=time.monotonic()
    print('TRANSPORT_START',stage,json.dumps(job),flush=True)
    worker={'surface':_surface_case,'promoter':promoter_job,'BZY':bzy_job}[stage]
    result=worker(job)
    print('TRANSPORT_END',stage,json.dumps(job),'wall_s',time.monotonic()-start,flush=True)
    return result


def remaining_transport(stage,jobs):
    """One restartable cache for actual completed cases, never invented rows."""
    path=HERE/'remaining_transport.json.gz'
    fingerprint=sha256(Path(T.__file__).read_bytes()+json.dumps([MATERIALS,OTHER],sort_keys=True).encode()).hexdigest()
    cache=json.loads(gzip.decompress(path.read_bytes())) if path.exists() else dict(input_fingerprint=fingerprint,stages={})
    if cache['input_fingerprint']!=fingerprint:
        raise ValueError('Remaining transport cache belongs to different physical source/inputs')
    saved=cache['stages'].setdefault(stage,{})
    key=lambda job:json.dumps(job,separators=(',',':'))
    pending=[job for job in jobs if key(job) not in saved]
    errors=[]
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        futures={pool.submit(_logged_remaining_job,(stage,job)):job for job in pending}
        for future in as_completed(futures):
            try:result=future.result()
            except Exception as error:
                errors.append((futures[future],repr(error)));print('TRANSPORT_FAILED',stage,errors[-1],flush=True)
                continue
            saved[key(futures[future])]=result
            temporary=path.with_suffix('.tmp')
            temporary.write_bytes(gzip.compress(json.dumps(cache,separators=(',',':')).encode(),compresslevel=1))
            temporary.replace(path)
    if errors:raise RuntimeError('Failed transport cases remain missing; valid cases checkpointed: '+repr(errors))
    return [saved[key(job)] for job in jobs]


def film_D_endpoint(profile,material,D_per_metal,retention):
    result=dict(profile,layers=[list(row) for row in profile['layers']])
    z,m,_,_=specs(material)
    if result['layers'][0][4]!=46:
        layer=result['layers'][0];den=material['density']
        layer[2]=den['host_atoms_m3'];layer[3]=den['metal_atoms_m3']*D_per_metal*retention
        layer[4]=z;layer[5]=m;layer[6]=den['metal_atoms_m3']
        result['nads_m2']=0. # no guessed oxide/cap adsorbate
    return result


def _crop(layers,removed):
    out=[]
    for row in layers:
        a=np.array(row,float);a[0]=max(a[0],removed)-removed;a[1]-=removed
        if a[1]>a[0]+1e-18:out.append(a.tolist())
    return out


def add_film(base,material,thickness_nm,D_per_metal=0.,retained=1.):
    layers,nads,state=base;d=thickness_nm*1e-9
    if d<=1e-12:return compact_profile(base)
    z,m,_,count=specs(material);den=material['density'];nm=den['metal_atoms_m3']
    nh=den.get('host_atoms_m3',nm)
    film=[0.,d,nh,nm*D_per_metal*retained,z,m,nm]
    shifted=[]
    for row in layers:
        a=np.array(row,float);a[:2]+=d;shifted.append(a.tolist())
    # A clean Pd surface sheet does NOT become an oxide/cap outer adsorbate.
    return dict(layers=[film]+shifted,nads_m2=0.,state=dict(state,film_nm=thickness_nm,
                retention_fraction=retained,adsorbate_status='covered Pd adsorbate removed; cap adsorption unknown'))


def _write_profiles(path,designs):
    with Path(path).open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['target','E_cm_keV','column','area_fraction','zt_nm','zb_nm','D_per_metal','n_D_m3','nads_m2','state','n_host_m3','n_metal_m3','Z_host','M_host_amu'])
        for design in designs:
            for e in FOUR:
                p=design['profiles'][str(float(e))];cov=design.get('coverage',1.)
                columns=[('covered',cov,p)]
                if cov<1:columns.append(('open',1-cov,compact_profile(reference().profiles(float(e)))))
                for column,area,profile in columns:
                    for row in profile['layers']:
                        w.writerow([design['name'],e,column,area,row[0]*1e9,row[1]*1e9,row[3]/row[6] if row[6]>0 else '',row[3],profile['nads_m2'],design.get('state_label','nominal run end'),row[2],row[6],row[4],row[5]])


def design_target(design,profiles=None,open_profiles=None,stopping_scale=1.):
    """One physical column/stopping construction for reports and age integrals."""
    if profiles is None:profiles={float(e):profile_tuple(p) for e,p in design['profiles'].items()}
    entries=tuple((a,b,tuple(tuple(s) for s in c)) for a,b,c in design.get('stopping_entries',()))
    target=Target(design['name'],profiles,stopping_fn=Stopping(entries,stopping_scale))
    coverage=design.get('coverage',1.)
    if coverage<1:
        target.columns=((coverage,profiles),(1-coverage,reference().profiles if open_profiles is None else open_profiles))
    return target


def report(designs,stage):
    rows=[];lifecycle=[]
    cachepath=HERE/'remaining_report_cache.json.gz' if stage=='remaining' else None
    if cachepath:
        source=Path(__file__).read_bytes()+Path(merit_api.__file__).read_bytes()
        source+=Path(sys.modules[merit_api.forward_yields.__module__].__file__).read_bytes()
        source+=json.dumps(merit_api.bestfits(),sort_keys=True).encode()
        fingerprint=sha256(source).hexdigest()
        cache=json.loads(gzip.decompress(cachepath.read_bytes())) if cachepath.exists() else dict(source_fingerprint=fingerprint,cases={})
        if cache['source_fingerprint']!=fingerprint:
            raise ValueError('Remaining report cache belongs to different consumer source/fits')
    for k,design in enumerate(designs):
        start=time.monotonic();offset=len(lifecycle)
        design_fingerprint=sha256(json.dumps(design,sort_keys=True,separators=(',',':')).encode()).hexdigest() if cachepath else None
        saved=cache['cases'].get(design['name']) if cachepath else None
        if saved:
            if saved['design_fingerprint']!=design_fingerprint:raise ValueError('Changed physical report input: '+design['name'])
            rows.extend(saved['rows']);lifecycle.extend(saved['sensitivities'])
            print('REUSED_REPORT',stage,k+1,len(designs),design['name'],flush=True)
            continue
        print('REPORT_START',stage,k+1,len(designs),design['name'],flush=True)
        target=design_target(design);profiles=target.profiles;entries=target.stopping_fn.entries
        result=compare(target)
        for row in result:
            row.update({x:design.get(x) for x in ('family','material','thickness_nm','scenario','coverage','feasibility','surface_status')})
            e=row['E_cm_keV'];life=design.get('lifetimes',{}).get(str(e),{})
            row.update(life)
        rows.extend(result)
        # Actual stopping sensitivity, including neutral reaction in the film.
        for scale in (.5,2.):
            t=Target(design['name']+'_stopping_'+str(scale),profiles,stopping_fn=Stopping(entries,scale))
            if target.columns:t.columns=target.columns
            for cause in CAUSES:
                values=evaluate(t,FOUR,cause,reference_E_cm_keV=None)['physical_yield_per_C']/evaluate(reference(),FOUR,cause,reference_E_cm_keV=None)['physical_yield_per_C']
                for e,v in zip(FOUR,values):
                    lifecycle.append(dict(target=design['name'],family=design['family'],scan='stopping',factor=scale,cause=cause,E_cm_keV=float(e),merit_PdD=float(v)))
        wall=time.monotonic()-start
        if cachepath:
            cache['cases'][design['name']]=dict(design_fingerprint=design_fingerprint,rows=result,sensitivities=lifecycle[offset:],wall_s=wall)
            temporary=cachepath.with_suffix('.tmp')
            temporary.write_bytes(gzip.compress(json.dumps(cache,separators=(',',':')).encode(),compresslevel=1));temporary.replace(cachepath)
        print('reported',stage,k+1,len(designs),design['name'],'wall_s',wall,flush=True)
    write_csv(HERE/(stage+'_compare.csv'),rows)
    write_csv(HERE/(stage+'_sensitivities.csv'),lifecycle)
    _write_profiles(HERE/(stage+'_depth_profiles.csv'),designs)
    return rows


def early():
    cache=HERE/'early_profiles.json.gz'
    if cache.exists():
        designs=json.loads(gzip.decompress(cache.read_bytes()))
    else:
        # Seven thicknesses resolve sub-monolayer islands and 1-2nm stopping.
        thickness=[.1,.25,.5,.75,1.,1.5,2.]
        jobs=[(m,d,p,float(e),None,1.) for m in ('Al2O3_seal','Er2O3_seal') for d in thickness for p in (.01,100.) for e in ENERGIES]
        transport_cache=HERE/'early_transport.json.gz'
        outcomes=json.loads(gzip.decompress(transport_cache.read_bytes())) if transport_cache.exists() else []
        completed={(r['material'],r['thickness_nm'],r['permeability_scale'],r['E_cm_keV']) for r in outcomes}
        jobs=[j for j in jobs if tuple(j[:4]) not in completed]
        if '--pilot' in sys.argv:jobs=jobs[:1]
        with ProcessPoolExecutor(max_workers=1 if '--pilot' in sys.argv else WORKERS) as pool:
            for i,row in enumerate(pool.map(transport_job,jobs)):
                outcomes.append(row)
                if (i+1)%14==0 or '--pilot' in sys.argv:
                    print('transported seal',i+1,len(jobs),'last wall_s',row['wall_s'],flush=True)
                    transport_cache.write_bytes(gzip.compress(json.dumps(outcomes,separators=(',',':')).encode()))
        transport_cache.write_bytes(gzip.compress(json.dumps(outcomes,separators=(',',':')).encode()))
        if '--pilot' in sys.argv:return
        designs=[]
        for m in ('Al2O3_seal','Er2O3_seal'):
            mat=MATERIALS[m];z,mass,fr,_=specs(mat)
            for d in thickness:
                for perm in (.01,100.):
                    selected=[r for r in outcomes if r['material']==m and r['thickness_nm']==d and r['permeability_scale']==perm]
                    profiles={str(r['E_cm_keV']):r['snapshots'][-1] for r in selected}
                    life={str(r['E_cm_keV']):dict(lifetime_nominal_h=r['lifetime_nominal_s']/3600,lifetime_low_h=r['lifetime_nominal_s']/36000,lifetime_high_h=r['lifetime_nominal_s']/360,lifetime_95_h=r['lifetime_95_s']/3600) for r in selected}
                    # Actual transported closed patches plus explicit islands.
                    for coverage in ((.25,.75,1.) if d<=.5 else (1.,)):
                        designs.append(dict(name=f'{m}_{d}nm_perm{perm}_coverage{coverage}',family=1,material=m,thickness_nm=d,scenario=f'permeability_x{perm}',coverage=coverage,
                                            profiles=profiles,lifetimes=life,stopping_entries=((z,mass,fr),),
                                            feasibility='oxide ALD known; closure/permeability on Pd not demonstrated',surface_status='pass4 retained estimates',
                                            erosion_snapshots={str(r['E_cm_keV']):r['snapshots'] for r in selected}))
        # Caps use actual PdD feed-state substrate; no unmeasured outgassing rate.
        for m in ('TiD2_cap','ErD2_cap','ScD2_cap','YD2_cap'):
            mat=MATERIALS[m];z,mass,fr,_=specs(mat)
            for d in (1.,2.,5.,10.,20.):
                for retention in (0.,1.):
                    for erosion_factor in (.1,10.):
                        prof={};life={}
                        for e in ENERGIES:
                            base=reference().profiles(float(e));rate,t=erosion(mat,float(e),d)
                            run=base[2]['beam_time_s'];remaining=max(0.,d-rate*erosion_factor*run*1e9)
                            prof[str(float(e))]=add_film(base,mat,remaining,2.,retention)
                            life[str(float(e))]=dict(lifetime_nominal_h=t/3600,lifetime_low_h=t/36000,lifetime_high_h=t/360,
                                                    remaining_film_nm=remaining,retention_fraction=retention,
                                                    diffusion_proxy_m2_s=3e-7*np.exp(-.60/(T.KB*base[2]['T_K'])),
                                                    diffusion_low_m2_s=.01*3e-7*np.exp(-.60/(T.KB*base[2]['T_K'])),
                                                    diffusion_high_m2_s=100*3e-7*np.exp(-.60/(T.KB*base[2]['T_K'])))
                        designs.append(dict(name=f'{m}_{d}nm_ret{retention}_erosion{erosion_factor}',family=2,material=m,thickness_nm=d,
                                            scenario=f'D_retained_fraction{retention};sputter_x{erosion_factor}',coverage=1.,profiles=prof,lifetimes=life,stopping_entries=((z,mass,fr),),
                                            feasibility='D-loaded cap assumed; ALD growth/loading/vacuum survival not demonstrated for all species',surface_status='Pd substrate pass4; cap surface adsorption unknown'))
        cache.write_bytes(gzip.compress(json.dumps(designs,separators=(',',':')).encode()))
    if '--prepare-early' in sys.argv:
        print('prepared early design inputs',len(designs),flush=True)
        return
    rows=report(designs,'early')
    erosion_rows=[]
    histories=matched_histories(designs)
    # Erosion snapshots: actual seal trajectory; cap D/geometry endpoints at
    # fixed substrate run-end state, labelled quasi-frozen rather than kinetics.
    for design in designs:
        if design['family']==1 and design['coverage']==1.:
            for fraction in (0.,.5,1.):
                mapping={};reference_mapping={};ages={}
                for e,snaps in design['erosion_snapshots'].items():
                    life=design['lifetimes'][e]['lifetime_nominal_h']*3600
                    desired=min(fraction*life,snaps[-1]['time_s'])
                    snapshot=min(snaps,key=lambda s:abs(s['time_s']-desired))
                    mapping[float(e)]=profile_tuple(snapshot);ages[float(e)]=snapshot['time_s']
                    if float(e) in FOUR:
                        reference_mapping[float(e)]=profile_tuple(histories[e][str(float(snapshot['time_s']))])
                target=Target(design['name']+'_erosion_'+str(fraction),mapping,
                              stopping_fn=Stopping(tuple((a,b,tuple(tuple(x) for x in c)) for a,b,c in design['stopping_entries'])))
                reference_at_age=Target('PdD_matching_age',reference_mapping,stopping_fn=transport_stopping)
                for row in compare(target,energies=FOUR):
                    e=row['E_cm_keV'];cause=row['cause']
                    same_age=evaluate(target,[e],cause,reference_E_cm_keV=None)['physical_yield_per_C'][0]/evaluate(reference_at_age,[e],cause,reference_E_cm_keV=None)['physical_yield_per_C'][0]
                    erosion_rows.append(dict(target=design['name'],scan='actual seal trajectory divided by PdD at same age and energy',age_lifetime_fraction=fraction,
                                             reference_age_s=ages[e],cause=cause,E_cm_keV=e,merit_PdD=float(same_age),
                                             merit_fixed_nominal_end_reference=row['merit_PdD']))
        if design['family']==2 and 'ret1.0_erosion0.1' in design['name']:
            mat=MATERIALS[design['material']]
            for remaining in (1.,.5,.05,0.):
                profiles={float(e):profile_tuple(add_film(reference().profiles(float(e)),mat,design['thickness_nm']*remaining,2.,1.)) for e in ENERGIES}
                target=Target(design['name']+'_remaining_'+str(remaining),profiles,stopping_fn=Stopping(tuple((a,b,tuple(tuple(x) for x in c)) for a,b,c in design['stopping_entries'])))
                for row in compare(target,energies=FOUR):
                    erosion_rows.append(dict(target=design['name'],scan='quasi-frozen cap erosion; retained-D endpoint',remaining_fraction=remaining,
                                             cause=row['cause'],E_cm_keV=row['E_cm_keV'],merit_PdD=row['merit_PdD']))
    write_csv(HERE/'early_erosion.csv',erosion_rows)
    # Return one useful four-energy decision table, without hiding causes.
    best=[]
    for family in (1,2):
        for cause in CAUSES:
            for e in FOUR:
                r=max((r for r in rows if r['family']==family and r['cause']==cause and r['E_cm_keV']==e),key=lambda r:r['merit_PdD'])
                best.append({k:r[k] for k in ('target','family','cause','E_cm_keV','merit_PdD','lifetime_nominal_h','lifetime_low_h','lifetime_high_h')})
    write_csv(HERE/'early_best.csv',best)
    print(json.dumps(best,indent=2),flush=True)


def remaining():
    """Designs3-6; same compare path and explicit thin-film endpoints."""
    designs=[]
    # Deuterated oxide: requested atomic fraction D/(D+metal+O).
    early_designs=json.loads(gzip.decompress((HERE/'early_profiles.json.gz').read_bytes()))
    for seal in (d for d in early_designs if d['family']==1 and d['coverage']==1. and d['thickness_nm'] in (.25,.5,1.,2.)):
        mat=MATERIALS[seal['material']]
        for atom_fraction in (.02,.15):
            for retention in (0.,1.):
                dm=2.5*atom_fraction/(1-atom_fraction)
                profiles={e:film_D_endpoint(p,mat,dm,retention) for e,p in seal['profiles'].items()}
                snapshots={e:[film_D_endpoint(p,mat,dm,retention) for p in ps] for e,ps in seal['erosion_snapshots'].items()}
                designs.append(dict(name='D_'+seal['name']+f'_at{atom_fraction}_ret{retention}',family=3,material=seal['material'],thickness_nm=seal['thickness_nm'],
                                    scenario='requested initial D/(D+M+O); retained/lost-D endpoints; actual transported seal and erosion',
                                    coverage=1.,profiles=profiles,lifetimes=seal['lifetimes'],stopping_entries=seal['stopping_entries'],erosion_snapshots=snapshots,
                                    feasibility='D2O tracer ALD known; beam-end D fraction unknown',surface_status='no qualified oxide D surface kinetics'))
    # Promoter chemistry is not guessed from work function or barrier. Actual
    # same Model loss-coefficient endpoints are shared by all three candidates.
    cases=[(f'loss0.5_j{j}',float(e),.5,j) for j in (0.,1.,2.) for e in ENERGIES]
    surface={}
    for name,e,p,temp in remaining_transport('surface',cases):surface.setdefault(name,{})[str(e)]=p
    jobs=[(name,s,float(e)) for name in ('BaO','SrO','La2O3') for s in (1e-7,.5) for e in ENERGIES]
    promoters=remaining_transport('promoter',jobs)
    for material in ('BaO','SrO','La2O3'):
        mat=promoter_material(material);z,mass,fr,_=specs(mat)
        for s in (1e-7,.5):
            selected=[r for r in promoters if r['material']==material and r['s']==s]
            profiles={str(r['E_cm_keV']):r['snapshots'][-1] for r in selected}
            life={str(r['E_cm_keV']):dict(lifetime_nominal_h=r['lifetime_s']/3600,
                  lifetime_low_h=r['lifetime_s']/36000,lifetime_high_h=r['lifetime_s']/360,
                  lifetime_status='one Pd-site/metal inventory; oxide sputter proxy, chemistry reverts as film erodes') for r in selected}
            for cov in (.1,.5):
                designs.append(dict(name=f'{material}_coverage{cov}_sticking{s}',family=4,material=material,thickness_nm=.3,coverage=cov,profiles=profiles,lifetimes=life,
                                    stopping_entries=((z,mass,fr),),erosion_snapshots={str(r['E_cm_keV']):r['snapshots'] for r in selected},
                                    scenario=f'chemistry-only loss coefficient{s}; one-site/metal inventory at assumed0.3nm island height; actual stopping and erosion',
                                    feasibility='partial promoter layer; actual surface kinetics unknown',surface_status='conditional pass4 sticking scan; no electronic screening change'))
    jobs=[(d,p,float(e)) for d in (2.,5.,10.,20.) for p in (.01,100.) for e in ENERGIES]
    bzy_rows=remaining_transport('BZY',jobs)
    # BZY: backside supply endpoint, not a predicted low-T conductivity gain.
    for d in (2.,5.,10.,20.):
        for factor in (0.,1.,2.):
            profiles=surface[f'loss0.5_j{factor}']
            lifetimes={str(float(e)):dict(lifetime_nominal_h=None,lifetime_low_h=None,lifetime_high_h=None,
                       lifetime_status='back surface not irradiated; electrochemical/chemical lifetime unknown') for e in ENERGIES}
            designs.append(dict(name=f'BZY_back_{d}nm_supply{factor}',family=5,material='BZY20_back',thickness_nm=d,coverage=1.,profiles=profiles,lifetimes=lifetimes,
                                scenario='blocked/baseline/double electrochemical supply endpoint; no invented300K conductivity law',
                                feasibility='ALD reported; low-T D conductivity and Pd interface not measured',surface_status='pass4 retained estimates; back layer sees no incident beam sputter'))
        for perm in (.01,100.):
            selected=[r for r in bzy_rows if r['d']==d and r['perm']==perm]
            for dpf in (0.,.2):
                material=bzy_material(dpf);z,mass,fr,_=specs(material)
                profiles={str(r['E_cm_keV']):film_D_endpoint(r['snapshots'][-1],material,dpf/2,1.) for r in selected}
                snapshots={str(r['E_cm_keV']):[film_D_endpoint(p,material,dpf/2,1.) for p in r['snapshots']] for r in selected}
                life={str(r['E_cm_keV']):dict(lifetime_nominal_h=r['lifetime_s']/3600,lifetime_low_h=r['lifetime_s']/36000,lifetime_high_h=r['lifetime_s']/360) for r in selected}
                designs.append(dict(name=f'BZY_front_{d}nm_perm{perm}_DperF{dpf}',family=5,material='BZY20_front',thickness_nm=d,coverage=1.,profiles=profiles,lifetimes=life,
                                    stopping_entries=((z,mass,fr),),erosion_snapshots=snapshots,
                                    scenario='actual proxy-permeance Pd transport and erosion; dry/full ideal hydration endpoint including added oxygen',
                                    feasibility='D mobility and stable hydration at beam T unknown',surface_status='front adsorption unknown; hydration is a frozen bound not a rate'))
    designs.extend(_laminate_designs())
    (HERE/'remaining_profiles.json.gz').write_bytes(gzip.compress(json.dumps(designs,separators=(',',':')).encode()))
    rows=report(designs,'remaining')
    write_csv(HERE/'all_designs.csv',list(csv.DictReader((HERE/'early_compare.csv').open()))+rows)
    remaining_erosion(designs)
    print('completed remaining designs',len(designs),flush=True)


def _laminate_designs():
    designs=[]
    # Actual nanolaminate through150nm; inventory endpoints, not a kinetic law.
    for oxide in ('Al2O3_seal','Er2O3_seal'):
        mat=MATERIALS[oxide];z,mass,fr,_=specs(mat)
        for od in (.2,.5):
            for loaded in (False,True):
                profiles={};life={};initial_profiles={}
                for e in ENERGIES:
                    base=reference().profiles(float(e));stack=[];top=0.;nm=6.8e28/(1+.16*.7)
                    for k in range(int(np.ceil(150/(1+od)))):
                        # interface inventory0 or1e19/m2 spread over0.1nm Pd.
                        nd=nm*.7 if loaded else 0.
                        stack.append([top,top+.1e-9,nm,nd+(1e29 if loaded and k>0 else 0),46.,106.42,nm]);top+=.1e-9
                        stack.append([top,top+.8e-9,nm,nd,46.,106.42,nm]);top+=.8e-9
                        stack.append([top,top+.1e-9,nm,nd+(1e29 if loaded else 0),46.,106.42,nm]);top+=.1e-9
                        den=mat['density'];stack.append([top,top+od*1e-9,den['host_atoms_m3'],0.,z,mass,den['metal_atoms_m3']]);top+=od*1e-9
                    film_layers=len(stack)
                    for row in substrate_interface_layers(base[0],loaded):
                        a=np.array(row,float);a[:2]+=top;stack.append(a.tolist())
                    initial_profiles[str(float(e))]=dict(layers=stack,nads_m2=0.,state=dict(base[2],film_total_m=top,
                                                             interface_count=2*int(np.ceil(150/(1+od))) if loaded else 0))
                    # Erosion through actual alternating inventories. Nominal
                    # lifetime is sum layer inventory/host sputter flux.
                    lt=0.
                    for a in stack[:film_layers]:
                        y=(.85*T.sputter_yield(2*e*1000,'Pd')+.30*T.sputter_yield(e*1000,'Pd'))/1.15 if a[4]==46 else sputter(mat,e)
                        lt+=(a[1]-a[0])*a[2]/(T.beam_flux(2*e)*y)
                    removed=0.;time_left=base[2]['beam_time_s']
                    for a in stack[:film_layers]:
                        y=(.85*T.sputter_yield(2*e*1000,'Pd')+.30*T.sputter_yield(e*1000,'Pd'))/1.15 if a[4]==46 else sputter(mat,e)
                        velocity=T.beam_flux(2*e)*y/a[2];duration=(a[1]-a[0])/velocity
                        if time_left<duration:removed+=velocity*time_left;break
                        removed+=a[1]-a[0];time_left-=duration
                    profiles[str(float(e))]=compact_profile(base) if removed>=top-1e-15 else dict(layers=_crop(stack,removed),nads_m2=0.,state=dict(base[2],interface_D_m2=1e19 if loaded else 0.,eroded_depth_m=removed))
                    life[str(float(e))]=dict(lifetime_nominal_h=lt/3600,lifetime_low_h=lt/36000,lifetime_high_h=lt/360)
                designs.append(dict(name=f'Pd1nm_{oxide}{od}nm_interfaces{loaded}',family=6,material=oxide,thickness_nm=150.,coverage=1.,profiles=profiles,initial_profiles=initial_profiles,lifetimes=life,stopping_entries=((z,mass,fr),),scenario='unloaded/no-interface-D or ideal hydride+one-ML interface inventory',feasibility='1nm Pd/ultrathin oxide continuity unproven; lateral percolation omitted',surface_status='interface capacity and transport conditional'))
    return designs

def remaining_erosion(designs):
    history=matched_histories(designs,stage='remaining')
    rows=[]
    for design in designs:
        stop=Stopping(tuple((a,b,tuple(tuple(x) for x in c)) for a,b,c in design.get('stopping_entries',())))
        if design.get('erosion_snapshots'):
            for fraction in (0.,.5,1.):
                mapping={};refmap={};ages={}
                for e,snaps in design['erosion_snapshots'].items():
                    life=design['lifetimes'][e]['lifetime_nominal_h']*3600
                    chosen=min(snaps,key=lambda s:abs(s['time_s']-min(fraction*life,snaps[-1]['time_s'])))
                    mapping[float(e)]=profile_tuple(chosen);ages[float(e)]=chosen['time_s']
                    if float(e) in FOUR:refmap[float(e)]=profile_tuple(history[e][str(float(chosen['time_s']))])
                target=Target(design['name'],mapping,stopping_fn=stop)
                cov=design.get('coverage',1.)
                if cov<1:target.columns=((cov,mapping),(1-cov,refmap))
                ref=Target('same_age_PdD',refmap,stopping_fn=transport_stopping)
                for cause in CAUSES:
                    ratio=evaluate(target,FOUR,cause,reference_E_cm_keV=None)['physical_yield_per_C']/evaluate(ref,FOUR,cause,reference_E_cm_keV=None)['physical_yield_per_C']
                    for e,v in zip(FOUR,ratio):rows.append(dict(target=design['name'],family=design['family'],cause=cause,E_cm_keV=float(e),age_s=ages[float(e)],age_lifetime_fraction=fraction,merit_PdD=float(v),scan='actual transported film trajectory; D-retention endpoint; same-age reference'))
        elif design['family']==6:
            for removed in (0.,.5,.95,1.):
                mapping={}
                for e,p in design['initial_profiles'].items():
                    if removed==1.:mapping[float(e)]=reference().profiles(float(e))
                    else:mapping[float(e)]=(_crop(p['layers'],removed*p['state']['film_total_m']),0.,p['state'])
                target=Target(design['name'],mapping,stopping_fn=stop)
                for row in compare(target,energies=FOUR):rows.append(dict(target=design['name'],family=6,cause=row['cause'],E_cm_keV=row['E_cm_keV'],removed_film_fraction=removed,merit_PdD=row['merit_PdD'],scan='actual top-layer sequence; quasi-frozen density endpoint; fixed run-end reference'))
        else:
            target=Target(design['name'],{float(e):profile_tuple(p) for e,p in design['profiles'].items()},stopping_fn=stop)
            for row in compare(target,energies=FOUR):rows.append(dict(target=design['name'],family=design['family'],cause=row['cause'],E_cm_keV=row['E_cm_keV'],merit_PdD=row['merit_PdD'],scan='backside not beam-eroded; chemical lifetime unknown'))
    write_csv(HERE/'remaining_erosion.csv',rows)


def substrate_interface_layers(layers,loaded):
    """Count the final oxide/substrate interface, never a free-surface trap."""
    out=[];width=.1e-9
    for row in layers:
        a=list(row)
        if loaded and a[0]<width:
            end=min(a[1],width);left=list(a);left[1]=end;left[3]+=1e29;out.append(left)
            if a[1]>width:
                a[0]=width;out.append(a)
        else:out.append(a)
    return out


def migrate_profile_caches():
    """Lossless density/state reuse after the omitted-film coordinate fix."""
    summary={}
    for filename in ('early_profiles.json.gz','early_transport.json.gz','dense_seal_profiles.json.gz','resolved_cache.json.gz'):
        path=HERE/filename
        if not path.exists():continue
        data=json.loads(gzip.decompress(path.read_bytes()));count=renamed=0;max_width_change=0.
        def visit(value):
            nonlocal count,renamed,max_width_change
            if isinstance(value,dict):
                for old,new in (('true_current_A','assumed_incident_current_A'),('run_charge_C','assumed_run_charge_C')):
                    if old in value:value[new]=value.pop(old);renamed+=1
                if 'layers' in value and 'state' in value:
                    rows=value['layers'];state=value['state'];d=float(state.get('d',0.))
                    if rows and rows[0][4]==46 and 0<rows[0][0]<=1e-12 and abs(rows[0][0]-d)<=1e-24:
                        shift=rows[0][0];density=tuple(tuple(row[2:]) for row in rows)
                        widths=[row[1]-row[0] for row in rows]
                        for row in rows:row[0]-=shift;row[1]-=shift
                        if 'zc' in state:state['zc']=[z-shift for z in state['zc']]
                        state['surface_offset_m']=0.
                        assert density==tuple(tuple(row[2:]) for row in rows)
                        max_width_change=max(max_width_change,max(abs(row[1]-row[0]-width) for row,width in zip(rows,widths)))
                        count+=1
                    return
                for child in value.values():
                    if isinstance(child,(dict,list)):visit(child)
            elif isinstance(value,list):
                for child in value:
                    if isinstance(child,(dict,list)):visit(child)
        visit(data)
        if count or renamed:
            temporary=path.with_suffix('.tmp')
            temporary.write_bytes(gzip.compress(json.dumps(data,separators=(',',':')).encode(),compresslevel=1))
            temporary.replace(path)
        summary[filename]=dict(profiles_rebased=count,assumption_fields_renamed=renamed,
                               max_binary64_layer_width_change_m=max_width_change,densities_nads_age_and_Model_solution_unchanged=True)
    (HERE/'coordinate_migration.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    if '--repair-coordinates' in sys.argv:migrate_profile_caches()
    elif '--remaining-probe' in sys.argv:
        index=sys.argv.index('--remaining-probe')
        stage=sys.argv[index+1];job=json.loads(sys.argv[index+2])
        remaining_transport(stage,[job])
        print('PROBE_COMPLETE',stage,json.dumps(job),flush=True)
    elif '--remaining-report-pilot' in sys.argv:
        candidates=[d for d in _laminate_designs() if any(p['state']['interface_count']>0 for p in d['initial_profiles'].values())]
        chosen=max(candidates,key=lambda d:len(d['initial_profiles']['0.25']['layers']))
        report([chosen],'remaining')
    elif '--early' in sys.argv:early()
    elif '--remaining' in sys.argv:remaining()
    else:raise SystemExit('Choose --early, --remaining or --remaining-probe STAGE JSON_JOB')
