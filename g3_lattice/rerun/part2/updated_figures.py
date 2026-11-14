"""Export only plotted slices/curves, then render without any simulation.
--export DIR --input FILE exports PdD; --export-ni DIR --ni-input FILE exports Ni.
--input FILE --ni-input FILE --out DIR renders the three combined figures only.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb


def export(folder,destination):
    results=json.loads((folder/'pd_results.json').read_text())['results']
    rows={row['label']:row for row in results};data={}
    def slice_field(prefix,label,interface=False):
        row=rows[label]
        with np.load(folder/(label+'.npz')) as fields:
            layer=fields['occupied'].shape[-1]-1 if interface else fields['occupied'].shape[-1]//2
            data[prefix+'_occupied']=fields['occupied'][:,:,layer]
            data[prefix+'_static']=fields['static_mean_A'][0,:,:,layer]
            data[prefix+'_amplitude']=fields['driven_amplitude_A'][0,:,:,layer]
        reference=row['loading_response']['composition_reference']
        data[prefix+'_a_A']=np.array(row['configuration'].get('geometry_a_A',reference['a_angstrom']))
        data[prefix+'_layer']=np.array(layer)
    for arrangement in ('random','clustered','ordered'):
        slice_field('arr_'+arrangement+'_off',arrangement+'_x0.875_h0')
        slice_field('arr_'+arrangement+'_on',arrangement+'_x0.875_h0.08')
    for name in ('pd389','pd4084'):
        slice_field('map_'+name+'_off','map_'+name+'_h0',True)
        slice_field('map_'+name+'_on','map_'+name+'_h0.08',True)
        with np.load(folder/(name+'_tensor.npz')) as field:
            matrix=np.moveaxis(field['delta_K_eV_A2'][:,:,:,:,-1],(0,1),(-2,-1))
            data['map_'+name+'_contrast_min']=np.linalg.eigvalsh(matrix)[...,0]
    labels=[('base','random_x0.875_h0.08'),('plus','C2_+0.5'),('minus','C2_-0.5')]
    for prefix,label in labels:
        slice_field('sensitivity_'+prefix,label)
        row=rows[label]
        data[prefix+'_frequency_THz']=np.asarray(row['spectrum']['frequency_THz'])
        data[prefix+'_local_PSD']=np.asarray(row['spectrum']['mean_local_PSD_A2_per_THz'])
        if prefix=='base':data['base_collective_PSD']=np.asarray(row['spectrum']['collective_PSD_A2_per_THz'])
        radial=row['half_pump_correlations']['phase_only']['radial']
        data[prefix+'_radius_nm']=np.asarray(radial['radius_nm'])
        data[prefix+'_phase_RMS_C']=np.asarray(radial['direction_RMS_modulus'])
    control=rows['random_x0.875_h0']['spectrum']
    data['off_frequency_THz']=np.asarray(control['frequency_THz'])
    data['off_local_PSD']=np.asarray(control['mean_local_PSD_A2_per_THz'])
    data['half_clock_THz']=np.array(rows['random_x0.875_h0.08']['stiffness_modulation_THz']/2)
    data['record_resolution_THz']=np.array(rows['random_x0.875_h0.08']['finite_record_resolution_THz'])
    data['metadata_json']=np.array(json.dumps(dict(source='Actual updated_batch.py snapshots and numeric curves',
        plane='fixed THIRD primitive index k is(001); horizontal[-110],vertical[110]',
        scope='Finite-window reference-model response; no new trajectories in rendering',
        phase='Joint DC+Omega+2Omega late-record x-displacement amplitude; pi alternation need not be a defect wall',
        maps='Conditional pure-metal EMT geometry, all sites outside PartA calibration; contrast eigenvalue is not full-operator stability')))
    # Retain only arrays actually plotted; all scientific summaries stay float64.
    for key in list(data):
        if key.endswith('_static') and ('_on_' in key or key.startswith('sensitivity_')):
            del data[key]
        elif key.endswith('_amplitude') and '_off_' in key:
            del data[key]
    frequency=data['base_frequency_THz']
    for prefix in ('base','plus','minus','off'):
        if not np.array_equal(data[prefix+'_frequency_THz'],frequency):
            raise ValueError('Figure curves do not share the declared sampling grid')
        del data[prefix+'_frequency_THz']
    data['frequency_THz']=frequency
    errors={}
    for key,value in list(data.items()):
        if key.endswith(('_static','_amplitude','_contrast_min')):
            reduced=value.astype(np.complex64 if np.iscomplexobj(value) else np.float32)
            errors[key]=float(np.max(abs(value-reduced)))
            data[key]=reduced
    metadata=json.loads(str(data['metadata_json']))
    metadata.update(display_storage='float32/complex64 display layers; ALL summaries and plotted numeric curves remain float64',
        display_cast_max_abs_error=errors)
    data['metadata_json']=np.array(json.dumps(metadata))
    np.savez_compressed(destination,**data)
    print(json.dumps({'exported':str(destination),'bytes':destination.stat().st_size,'arrays':len(data)}))


def extent(data,prefix):
    shape=data[prefix+'_occupied'].shape;step=float(data[prefix+'_a_A'])*.1/np.sqrt(2)
    return (0,shape[1]*step,0,shape[0]*step)


def phase_picture(amplitude,mask,cap):
    hue=(np.angle(amplitude)+np.pi)/(2*np.pi)
    brightness=abs(amplitude)/cap if cap else np.zeros(amplitude.shape)
    image=hsv_to_rgb(np.stack((hue,np.ones(hue.shape),brightness),axis=-1))
    image[~mask]=0
    return image


def axes_labels(axis):
    axis.set_xlabel('[-1 1 0] reference nm');axis.set_ylabel('[1 1 0] reference nm')


def render(source,ni_source,folder):
    folder.mkdir(parents=True,exist_ok=True)
    with np.load(ni_source) as stored:
        ni={name:stored[name] for name in stored.files}
    with np.load(source) as data:
        arrangements=('random','clustered','ordered')
        cap=max(float(np.max(abs(data['arr_'+name+'_on_amplitude']))) for name in arrangements)
        limit=max(float(np.max(abs(data['arr_'+name+'_off_static']))) for name in arrangements)
        fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
        for col,name in enumerate(arrangements):
            off='arr_'+name+'_off';on='arr_'+name+'_on'
            image=np.ma.array(data[off+'_static'],mask=~data[off+'_occupied'])
            artist=axes[0,col].imshow(image,origin='lower',extent=extent(data,off),cmap='coolwarm',vmin=-limit,vmax=limit)
            axes[0,col].set_title(name+' vacancies\nno pump')
            axes[1,col].imshow(phase_picture(data[on+'_amplitude'],data[on+'_occupied'],cap),origin='lower',extent=extent(data,on))
            axes[1,col].set_title(name+' vacancies\nx-channel Omega')
            for row in range(2):axes_labels(axes[row,col])
        fig.colorbar(artist,ax=axes[0].tolist(),label='mean u_[100], angstrom')
        fig.suptitle('x=.875; (001), fixed third primitive index. Bottom:hue=phase,brightness=|A_Omega|\n'
            f'Common brightness cap {cap:.4g} angstrom; pi alternation is not automatically a domain wall')
        fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-np.pi,np.pi),cmap='hsv'),
            ax=axes[1].tolist(),label='phase of x-channel Omega,rad')
        fig.savefig(folder/'arrangements.png',dpi=95);plt.close(fig)
        cap=max(float(np.max(abs(data['map_'+name+'_on_amplitude']))) for name in ('pd389','pd4084'))
        fig,axes=plt.subplots(2,3,figsize=(12,7),constrained_layout=True)
        for row,name in enumerate(('pd389','pd4084')):
            off='map_'+name+'_off';on='map_'+name+'_on'
            artist=axes[row,0].imshow(data['map_'+name+'_contrast_min'],origin='lower',extent=extent(data,off),cmap='coolwarm')
            fig.colorbar(artist,ax=axes[row,0],label='min elastic-contrast eigenvalue,eV/A2')
            artist=axes[row,1].imshow(np.ma.array(data[off+'_static'],mask=~data[off+'_occupied']),origin='lower',extent=extent(data,off),cmap='coolwarm')
            fig.colorbar(artist,ax=axes[row,1],label='mean u_[100],angstrom')
            axes[row,2].imshow(phase_picture(data[on+'_amplitude'],data[on+'_occupied'],cap),origin='lower',extent=extent(data,on))
            for col,title in enumerate(('contrast, NOT full stability','no-pump displacement','separated Omega response')):
                axes[row,col].set_title(name+': '+title);axes_labels(axes[row,col])
        fig.suptitle('Conditional Ni-on-Pd(001) registered-map response,interface layer\n'
            'Pure-metal EMT geometry; all sites extrapolate PartA. Phase brightness cap '+f'{cap:.4g} angstrom')
        fig.savefig(folder/'interface_maps.png',dpi=95);plt.close(fig)
        cap=max(float(np.max(abs(data['sensitivity_'+name+'_amplitude']))) for name in ('base','plus','minus'))
        cap=max(cap,float(np.max(abs(ni['ni_amplitude']))))
        fig=plt.figure(figsize=(13,8),constrained_layout=True);grid=fig.add_gridspec(2,8)
        for index,(name,title) in enumerate((('base','PdD C2=0 assumed'),('plus','PdD C2=+C1/2'),('minus','PdD C2=-C1/2'))):
            axis=fig.add_subplot(grid[0,index*2:index*2+2]);prefix='sensitivity_'+name
            axis.imshow(phase_picture(data[prefix+'_amplitude'],data[prefix+'_occupied'],cap),origin='lower',extent=extent(data,prefix))
            axis.set_title(title);axes_labels(axis)
        axis=fig.add_subplot(grid[0,6:])
        axis.imshow(phase_picture(ni['ni_amplitude'],ni['ni_occupied'],cap),origin='lower',extent=extent(ni,'ni'))
        axis.set_title('NiH893/401,mixed');axes_labels(axis)
        spectrum=fig.add_subplot(grid[1,:4]);correlation=fig.add_subplot(grid[1,4:])
        for name,title in (('base','PdD base'),('plus','PdD +C1/2'),('minus','PdD -C1/2')):
            f=data['frequency_THz'];p=data[name+'_local_PSD'];keep=(p>0)&(f<=4*float(data['half_clock_THz']))
            spectrum.semilogy(f[keep],p[keep],label=title)
            correlation.plot(data[name+'_radius_nm'],data[name+'_phase_RMS_C'],label=title)
        f=data['frequency_THz'];p=data['off_local_PSD'];keep=(p>0)&(f<=4*float(data['half_clock_THz']))
        spectrum.semilogy(f[keep],p[keep],':',color='gray',label='PdD same mask,no pump')
        correlation.plot(ni['ni_radius_nm'],ni['ni_phase_RMS_C'],'--',label='NiH893/401,110592sites')
        for multiplier in (1,2):spectrum.axvline(multiplier*float(data['half_clock_THz']),color='black',alpha=.3)
        spectrum.set_xlabel('frequency,THz');spectrum.set_ylabel('mean probe PSD,A2/THz');spectrum.legend()
        correlation.set_xlabel('reference separation,nm');correlation.set_ylabel('direction-RMS |C_phase|');correlation.legend()
        fig.suptitle('PdD:same32^3 random x=.875 mask/seed; C2 scenarios are NOT calibrated bounds\n'
            'NiH:one matched41.45ps case,mixed quartics. Hue=phase;brightness common cap '+f'{cap:.3g} angstrom\n'
            'PdD joint-harmonic fields;NiH matched legacy filter. PSD resolution '+f'{float(data["record_resolution_THz"]):.4g} THz')
        fig.savefig(folder/'sensitivity_spectra.png',dpi=95);plt.close(fig)
    print(json.dumps({'rendered':[str(folder/name) for name in ('arrangements.png','interface_maps.png','sensitivity_spectra.png')]}))


def export_ni(folder,destination):
    result=json.loads((folder/'phase.json').read_text())
    with np.load(folder/'phase_arrays.npz') as source:
        layer=source['occupied'].shape[-1]//2
        amplitude=source['amplitude_x_A'][:,:,layer]
        reduced=amplitude.astype(np.complex64)
        occupied=source['occupied'][:,:,layer]
    radial=result['matched_legacy_filter']['phase_only']['radial']
    parameters=json.loads((folder/'parameters.json').read_text())
    metadata=dict(source='Actual ni_epoch.py ONE trajectory, legacy-filter display layer',
        quartic_epoch=result['quartic_epoch'],plane='(001),fixed THIRD primitive index',
        display_storage='complex64; full numeric phase/quantum summaries remain float64',
        display_cast_max_abs_error_A=float(np.max(abs(amplitude-reduced))))
    np.savez_compressed(destination,ni_amplitude=reduced,ni_occupied=occupied,
        ni_a_A=np.array(parameters['materials']['NiH']['a_A']),ni_layer=np.array(layer),
        ni_radius_nm=np.asarray(radial['radius_nm']),ni_phase_RMS_C=np.asarray(radial['direction_RMS_modulus']),
        metadata_json=np.array(json.dumps(metadata)))
    print(json.dumps({'exported':str(destination),'bytes':destination.stat().st_size}))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--export',type=Path);ap.add_argument('--input',type=Path)
    ap.add_argument('--export-ni',type=Path);ap.add_argument('--ni-input',type=Path)
    ap.add_argument('--out',type=Path,default=Path('figures'));args=ap.parse_args()
    if args.export_ni:
        if args.ni_input is None:raise ValueError('Ni export needs its actual target file')
        export_ni(args.export_ni,args.ni_input)
    elif args.export:
        if args.input is None:raise ValueError('PdD export needs its actual target file')
        export(args.export,args.input)
    else:
        if args.input is None or args.ni_input is None:
            raise ValueError('Final figure consumer needs both actual display inputs')
        render(args.input,args.ni_input,args.out)


if __name__=='__main__':main()
