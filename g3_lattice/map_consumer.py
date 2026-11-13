"""Consume ACTUAL registered octahedral maps, not Pd-metal rows.
python map_consumer.py lattice_inputs.json --source mismatch_exports --out map_results
Uses TOTAL law-bulk contrast by default, retaining the3.89/4.084 background.
All dynamics reuse lattice.run_case; no EMT, DFT, or geometry replay.
"""
import argparse
import json
from math import ceil, pi, sqrt
from pathlib import Path
import numpy as np
from lattice import CELL, HBAR, HC, shell_kernel, run_case, regions, png


def native_mode(p,delta):
    shape=delta.shape[2:]
    N=int(np.prod(shape))
    fftshape=shape[:2]+(shape[2]+4,)
    kernel=shell_kernel(p,fftshape)
    A=p['K4_100_eV_A4'];B=(3*p['K4_111_eV_A4']-A)/2
    K=p['Kii_eV_A2']+(3*A+2*B)*p['zero_point_axis_variance_A2']
    scalar=np.eye(N)*K
    coords=np.indices(shape).reshape(3,-1).T
    row=np.arange(N)
    for index in np.argwhere(abs(kernel)>0):
        offset=np.array([v if v<=fftshape[d]//2 else v-fftshape[d] for d,v in enumerate(index)])
        targets=coords-offset
        # Open depth, periodic in-plane: same truncation as lattice.run_case.
        valid=(targets[:,2]>=0)&(targets[:,2]<shape[2])
        targets=targets[valid]
        targets[:,:2]%=shape[:2]
        col=np.ravel_multi_index(targets.T,shape)
        scalar[row[valid],col]-=kernel[tuple(index)]
    H=np.kron(np.eye(3),scalar)
    for a in range(3):
        for b in range(3):
            H[a*N+row,b*N+row]+=delta[a,b].ravel()
    eigen,vectors=np.linalg.eigh(H)
    positive=np.flatnonzero(eigen>0)
    if not len(positive):
        raise ValueError('No positive native-cell mode for the requested pump')
    selected=int(positive[0])
    mode=vectors[:,selected].reshape(3,N)
    covariance=mode@mode.T
    weights,axes=np.linalg.eigh(covariance)
    direction=axes[:,-1]
    weight=float(weights[-1])
    stiffness=float(eigen[selected])
    Omega=sqrt(HC*HC/p['mass_eVc2']*stiffness)/HBAR
    G=K-sum(s['J_sum_eV_A2'] for s in p['shells'])
    Gamma2=1e12/p['reference_amplitude_decay_ps']
    layer_weights=(mode.reshape((3,)+shape)**2).sum(axis=(0,1,2))
    return dict(lowest_native_eigenvalue_eV_A2=float(eigen[0]),
        selected_eigenvalue_eV_A2=stiffness,
        frequency_THz=Omega/(2*pi*1e12),
        pump_ratio=sqrt(stiffness/G),pump_polarization=direction.tolist(),
        pump_mode_polarization_weight=weight,
        weak_modulation_h_threshold=4*Gamma2*stiffness/(Omega*G*weight),
        native_mode_layer_weights=layer_weights.tolist(),
        scope='lowest positive in-plane-Bloch-Gamma native motif, open depth; not a search over all Bloch q')


def correlations(fields,delta,period):
    amplitude=np.sqrt((abs(fields['field'])**2).sum(axis=0))
    softness=-np.linalg.eigvalsh(np.moveaxis(delta,(0,1),(-2,-1)))[...,0]
    corr=[];repeat=[];enrich=[]
    for k in range(amplitude.shape[2]):
        y=amplitude[:,:,k];x=softness[:,:,k]
        corr.append(float(np.corrcoef(x.ravel(),y.ravel())[0,1]) if np.ptp(x)>0 and np.ptp(y)>0 else None)
        scale=float(np.mean((y-y.mean())**2))
        repeat.append([float(np.sqrt(np.mean((y-np.roll(y,period,axis=d))**2)/scale))
                       if np.ptp(y)>0 and scale>0 else None for d in (0,1)])
        selected=x>=np.quantile(x,.75)
        total=float((y*y).sum())
        enrich.append(float((y[selected]**2).sum()/total/selected.mean()) if np.ptp(x)>0 and total>0 else None)
    return dict(softness_vs_response_per_layer=corr,
        motif_translation_relative_RMS_per_layer=repeat,
        softest_quartile_response_energy_enrichment_per_layer=enrich,
        interpretation='amplitude envelope, not just phase S(q); soft quartile is descriptive, not a phase threshold')


def picture(folder,label,fields,delta,direction):
    nx,ny,nz=fields['occupied'].shape
    minimum=np.linalg.eigvalsh(np.moveaxis(delta,(0,1),(-2,-1)))[...,0]
    amplitude=np.sqrt((abs(fields['field'])**2).sum(axis=0))
    rows=[]
    for k in (0,nz-1):
        soft=np.clip(minimum[:,:,k],-1,1)
        projected=np.einsum('a,aij->ij',np.asarray(direction),fields['field'][:,:,:,k]).real
        phase=np.clip(projected/.25,-1,1)
        first=np.stack((255*(1-np.maximum(-soft,0)),255*(1-abs(soft)),255*(1-np.maximum(soft,0))),axis=-1)
        second=np.zeros((nx,ny,3));second[:,:,0]=255*np.clip(amplitude[:,:,k]/.15,0,1)
        second[:,:,1]=second[:,:,0]
        third=np.stack((255*(1-np.maximum(-phase,0)),255*(1-abs(phase)),255*(1-np.maximum(phase,0))),axis=-1)
        rows.append(np.concatenate((first,second,third),axis=1))
    png(folder/(label+'.png'),np.concatenate(rows,axis=0))


def refresh_observables(parameters,folder,source_root,capture):
    """Reuse valid dynamics; execute ONLY an explicitly named missing snapshot."""
    path=folder/'map_response.json'
    result=json.loads(path.read_text())
    for index,summary in enumerate(result['results']):
        c=summary['configuration']
        case=c['label'].split('_')[0]
        with np.load(Path(source_root)/case/'fcc_octahedral_grid.npz') as source:
            period=source['DeltaK_law_bulk_eV_A2'].shape[2]
        with np.load(c['tensor_input']) as source:
            delta=source['delta_K_eV_A2']
        snapshot=folder/(c['label']+'.npz')
        fields=None
        if c['label']==capture:
            print('Executing explicitly requested missing-field case '+capture,flush=True)
            updated,fields=run_case(parameters,c)
            updated['linear_positive_lower_bound_eV_A2']=summary['linear_positive_lower_bound_eV_A2']
            summary=updated
            np.savez_compressed(snapshot,**fields)
        elif snapshot.exists():
            with np.load(snapshot) as stored:
                fields={key:stored[key] for key in stored.files}
        if fields is not None:
            summary['map_response']=correlations(fields,delta,period)
            summary['regions']=regions(fields,c['geometry_a_A'])
            if summary['active_fraction']>0:
                picture(folder,c['label'],fields,delta,c['pump_polarization'])
        elif summary['active_fraction']==0:
            summary['regions']={'definition':'No sites exceed the declared amplitude cutoff',
                'component_sizes':[],'largest_shapes':[],'raw_neighbour_polarization_alignment':None}
        else:
            raise ValueError('Missing active-field snapshot; capture this case explicitly: '+c['label'])
        softness=-np.linalg.eigvalsh(np.moveaxis(delta,(0,1),(-2,-1)))[...,0]
        for layer in range(softness.shape[2]):
            if np.ptp(softness[:,:,layer])==0:
                summary['map_response']['softness_vs_response_per_layer'][layer]=None
                summary['map_response']['softest_quartile_response_energy_enrichment_per_layer'][layer]=None
        result['results'][index]=summary
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('parameters');ap.add_argument('--source',required=True)
    ap.add_argument('--out',default='map_results')
    ap.add_argument('--refresh-observables',action='store_true')
    ap.add_argument('--capture-case')
    args=ap.parse_args()
    parameters=json.loads(Path(args.parameters).read_text())
    p=parameters['materials']['PdD']
    folder=Path(args.out);folder.mkdir(exist_ok=True)
    if args.refresh_observables:
        refresh_observables(parameters,folder,args.source,args.capture_case)
        return
    results=[];native_results={}
    for case in ('pd389','pd4084'):
        with np.load(Path(args.source)/case/'fcc_octahedral_grid.npz') as source:
            if not np.allclose(source['CELL_dimensionless'],CELL,rtol=0,atol=1e-12):
                raise ValueError('Map and consumer FCC bases differ')
            if not np.allclose(source['O_to_cubic_lab'],np.eye(3),rtol=0,atol=1e-12):
                raise ValueError('Map is not in the declared cubic Cartesian frame')
            native=source['DeltaK_law_bulk_eV_A2'].copy()
            a=float(source['a_Pd_A'])
        mode=native_mode(p,native)
        native_results[case]=mode
        n=native.shape[2]
        tile=ceil(sqrt(10000/(n*n*native.shape[4])))
        delta=np.tile(native,(1,1,tile,tile,1))
        shape=delta.shape[2:]
        mean=np.broadcast_to(delta.mean(axis=(2,3),keepdims=True),delta.shape).copy()
        full_path=folder/(case+'_full_tensor.npz')
        mean_path=folder/(case+'_layer_mean_tensor.npz')
        np.savez_compressed(full_path,delta_K_eV_A2=delta)
        np.savez_compressed(mean_path,delta_K_eV_A2=mean)
        for suffix,h,path,ratio in (('off',0.,full_path,mode['pump_ratio']),
                                   ('resonant',.08,full_path,mode['pump_ratio']),
                                   ('layer_mean_resonant',.08,mean_path,mode['pump_ratio']),
                                   ('bulk_clock',.08,full_path,1.)):
            config=dict(label=case+'_'+suffix,material='PdD',L=list(shape),x=1.,steps=2048,dt=.3,
                boundary_z='open',tensor_input=str(path),tensor_kind='delta_collective_curvature',
                geometry_a_A=a,pump_h=h,pump_ratio=ratio,pump_polarization=mode['pump_polarization'],
                map_contrast='TOTAL DeltaK_law_bulk; includes3.89/4.084 background compression',
                map_limits='normal-only extrapolation; EMT metal geometry is not PdD elasticity; shear response missing')
            summary,fields=run_case(parameters,config)
            chosen=mean if 'layer_mean' in suffix else delta
            summary['map_response']=correlations(fields,chosen,n)
            kernel=shell_kernel(p,shape[:2]+(shape[2]+4,))
            A=p['K4_100_eV_A4'];B=(3*p['K4_111_eV_A4']-A)/2
            K=p['Kii_eV_A2']+(3*A+2*B)*p['zero_point_axis_variance_A2']
            site_min=float(np.linalg.eigvalsh(np.moveaxis(chosen,(0,1),(-2,-1)))[...,0].min())+K
            summary['linear_positive_lower_bound_eV_A2']=site_min-float(np.fft.fftn(kernel).real.max())
            summary['regions']=regions(fields,a)
            if summary['active_fraction']>0:
                np.savez_compressed(folder/(config['label']+'.npz'),**fields)
                picture(folder,config['label'],fields,chosen,config['pump_polarization'])
            results.append(summary)
            print(json.dumps({k:v for k,v in summary.items() if k not in ('spectrum','regions')}),flush=True)
    (folder/'map_response.json').write_text(json.dumps(dict(
        contrast='TOTAL law-bulk collective contrast transferred additively to unchanged local+sameJ baseline',
        source_site='registered octahedral centers; not metal rows',
        source_limits='all source sites extrapolate outside +/-1% normal-strain data; independent shear response unknown; pure-metal EMT geometry',
        native_modes=native_results,results=results),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
