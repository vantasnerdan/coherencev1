"""Figures from saved joined fields and final consumers; no producer execution."""
import csv
import numpy as np
import sources as s
import native


def render(model):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    extent=[-s.L*500,s.L*500,-s.L*500,s.L*500];arm=s.NOMINAL
    for beam in s.BEAMS:
        field=np.load(model.field_path(beam,arm['name']));off={k:np.load(model.field_path(beam,arm['name'],control=k)) for k in ('temperature_off','strain_off','frozen')}
        panels=((field['t2_omega_over_kappa'][:,:,0],'evolved Ω/κE, dimensionless'),
                (field['t2_bare_K'][:,:,0],'native bare Γ curvature, eV/Å²'),
                (field['t2_retained_shell_J_over_onsite'][:,:,0],'frozen shell sum / evolved onsite'),
                ((field['t2_K']-off['temperature_off']['t2_K'])[:,:,0],'joined − T-off pole, eV/Å²'),
                ((field['t2_K']-off['strain_off']['t2_K'])[:,:,0],'joined − strain-off pole, eV/Å²'),
                ((field['t2_deltaT']-off['frozen']['t2_deltaT'])[:,:,0],'joined − frozen T, K'))
        fig,axes=plt.subplots(2,3,figsize=(13,8))
        for ax,(values,title) in zip(axes.flat,panels):
            image=ax.imshow(values.T,origin='lower',extent=extent,cmap='RdBu_r' if '−' in title else 'viridis')
            fig.colorbar(image,ax=ax,shrink=.75);ax.set(title=title,xlabel='x (mm)',ylabel='y (mm)')
        fig.suptitle(beam+': same source; static medium and feedback are separate\nnominal / onsite, 2.5 keV CM, 2 s; depth-cell aperture 0–250 nm')
        fig.tight_layout();fig.savefig(s.HERE/f'medium_{beam}.png',dpi=155);plt.close(fig)
    # Actual source quadrature support, not a coarse-cell energy divided by n_D.
    fig,axes=plt.subplots(3,3,figsize=(14,11));p=s.inputs()['lattice'];C=2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4']
    for row,closure in enumerate(s.PHYSICAL_ARMS):
        name=closure['name']+'_onsite';data=np.load(s.HERE/f'native_E2.5_{name}_sp1_m16.npz');times=data['time_s']*1e15
        index=int(np.argmax(np.max(data['map_U'],axis=1)));positions=data['positions_r_z_nm']
        image=axes[row,0].scatter(positions[:,0],positions[:,1],c=data['map_U'][index],s=12,vmin=0,vmax=.031)
        fig.colorbar(image,ax=axes[row,0],shrink=.7,label='Ue (eV / fixed formula)');axes[row,0].set(xlabel='radial quadrature r (nm)',ylabel='depth z (nm)',title=closure['name']+f': Ue at {data["map_time_s"][index]*1e15:.2f} fs')
        zoom=axes[row,0].inset_axes([.52,.53,.44,.41])
        zoom.scatter(positions[:,0],positions[:,1],c=data['map_U'][index],s=14,vmin=0,vmax=.031)
        zoom.set(xlim=(0,1),ylim=(0,25),title='same points: zoom')
        zoom.tick_params(labelsize=6);zoom.title.set_fontsize(7)
        local=native.Material(16).value(data['maximum_U_eV_formula']);axes[row,1].plot(times,local,label='local-peak U mapped to ΔKΓ')
        axes[row,1].plot(times,data['projected_deltaK_eVA2'],label='actual focus-projected ΔKΓ')
        axes[row,1].set(xscale='log',xlabel='same electrical-parent age (fs)',ylabel='eV/Å²',title='same curve, different spatial/exposure weighting');axes[row,1].legend(fontsize=7)
        for allocation in ('onsite','coupling'):
            a=np.load(s.HERE/f'native_E2.5_{closure["name"]}_{allocation}_sp1_m16.npz')
            axes[row,2].plot(a['time_s']*1e15,np.linalg.norm(a['group_velocity_selected_m_s'],axis=-1).max(axis=1),label=allocation)
        axes[row,2].set(xscale='log',xlabel='same electrical-parent age (fs)',ylabel='selected optical max |vg| (m/s)',title='actual G3/native branch; not Pd sound speed');axes[row,2].legend(fontsize=8)
    fig.suptitle('Native-G16, mesh16: D₂ source at 2.5 keV CM; low-U mesh sign unresolved; same joined foil')
    fig.tight_layout();fig.savefig(s.HERE/'native_support.png',dpi=155);plt.close(fig)
    with (s.HERE/'yield_mixture.csv').open() as stream:rows=list(csv.DictReader(stream))
    rows=[r for r in rows if r['branch']=='p' and r['comparator']=='static40'];energies=np.array(s.ENERGIES)
    primary=[r for r in rows if r['mesh']=='16' and r['source_coupling']=='native_on']
    lower=[];upper=[];wide_lower=[];wide_upper=[];required=[];absolute=[]
    for e in energies:
        selected=[r for r in primary if float(r['E_cm_keV'])==e];allrows=[r for r in rows if float(r['E_cm_keV'])==e]
        lower.append(min(float(r['fraction_of_G12_residual']) for r in selected));upper.append(max(float(r['fraction_of_G12_residual']) for r in selected))
        wide_lower.append(min(float(r['fraction_of_G12_residual']) for r in allrows));wide_upper.append(max(float(r['fraction_of_G12_residual']) for r in allrows))
        required.append(float(selected[0]['G12_required_extra_per_D']));absolute.append(max(float(r['built_field_addition_per_D']) for r in selected))
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    axes[0].fill_between(energies,wide_lower,wide_upper,color='lightgray',label='plus labelled source mesh/coupling contrasts')
    axes[0].fill_between(energies,lower,upper,alpha=.6,label='three closures × two allocations × three beams')
    axes[0].axhline(1,color='black',ls=':',label='accepted residual');axes[0].set(xlabel='Ecm (keV)',ylabel='built-field addition / G12 residual');axes[0].legend(fontsize=7)
    axes[1].plot(energies,required,'o-',label='required extra, own raw reference');axes[1].plot(energies,absolute,'s-',label='largest primary built-field addition')
    axes[1].set(xlabel='Ecm (keV)',ylabel='reactions per incident D');axes[1].legend(fontsize=7)
    for ax,values in ((axes[0],wide_lower+wide_upper),(axes[1],required+absolute)):
        if min(values)>0:ax.set_yscale('log')
        else:
            nonzero=np.abs(np.asarray(values));nonzero=nonzero[nonzero>0]
            ax.set_yscale('symlog',linthresh=float(nonzero.min()) if len(nonzero) else 1.)
    for ax in axes:ax.axvline(2.25,color='gray',ls='--',lw=.7)
    fig.suptitle('Yield last: physical velocity covariances, quiet once; named model ranges, not confidence/contact bounds\nFrozen profile preparation changes from 12 h to 2 h between 2 and 2.5 keV')
    fig.tight_layout();fig.savefig(s.HERE/'yield_comparison.png',dpi=155);plt.close(fig)
