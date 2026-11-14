"""Plot saved correlation data only; no dynamics or invented phase classifier."""
import argparse
import colorsys
import json
from pathlib import Path
import numpy as np
from lattice import CELL,png


def binned(curve,width,modulus=False):
    r=np.array(curve['radius_nm']);w=np.array(curve['occupied_pairs'],float)
    c=np.array(curve['direction_RMS_modulus'])**2 if modulus else np.array(curve['C_real'])
    bins=np.rint(r/width).astype(int)
    denominator=np.bincount(bins,weights=w)
    keep=denominator>0
    value=np.bincount(bins,weights=w*c)[keep]/denominator[keep]
    return (np.bincount(bins,weights=w*r)[keep]/denominator[keep],
            np.sqrt(value) if modulus else value)


def correlation_svg(folder,output,labels,window=False,modulus=False):
    # Matplotlib is only a rendering dependency; the dynamics remains NumPy-only.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['svg.fonttype']='path'
    fig,ax=plt.subplots(figsize=(9.4,5.5))
    for label in labels:
        row=json.loads((folder/(label+'.json')).read_text())
        curve=(row['windowed_coherence']['phase_radial'] if window else
               row['coherence']['phase_only']['carrier_removed_radial'])
        r,c=binned(curve,.26389225073881956,modulus)
        ax.plot(r,c,label=label)
    title=('Direction-RMS equal-time phase correlation' if modulus else
           ('Equal-time short-envelope phase correlations' if window else
            'Finite-q-carrier-removed phase correlations; plateau retained'))
    ax.set(title=title,xlabel='r (nm)',ylabel='C phase (declared projection)',ylim=(-.2,1.03))
    ax.grid(alpha=.25)
    ax.legend(fontsize=8)
    fig.text(.5,.015,'Actual occupied pairs; visual bins0.264nm. No plateau subtraction or exponential fit; unbinned data retained locally.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.04,1,1))
    fig.savefig(output)
    fig.savefig(output.with_suffix('.png'),dpi=150)
    plt.close(fig)


def phase_picture(folder,label,a_A,from_planes=False):
    result=json.loads((folder/(label+'.json')).read_text())
    kstar=np.array(result['coherence']['amplitude_weighted']['peak_index'])
    wavelength=result['coherence']['amplitude_weighted']['peak_q']['wavelength_nm']
    sigma=.5*wavelength
    if from_planes:
        with np.load(folder/(label+'_planes.npz')) as data:
            raw=data['quadrature'];b=data['envelope'];amp=data['amplitude'];mask=data['occupied']
    else:
        with np.load(folder/(label+'.npz')) as data:
            a=data['amplitude_x_A'];occupied=data['occupied']
        shape=a.shape;grid=np.indices(shape)
        phase=2*np.pi*sum(kstar[d]*grid[d]/shape[d] for d in range(3))
        signed=np.array([np.where(grid[d]<=shape[d]//2,grid[d],grid[d]-shape[d]) for d in range(3)])
        displacement=np.einsum('ab,b...->a...',CELL,signed)*a_A*.1
        kernel=np.exp(-np.sum(displacement**2,axis=0)/(2*sigma*sigma))
        inverse=kernel[np.ix_(*[(-np.arange(n))%n for n in shape])]
        kernel=(kernel+inverse)/2
        weight=np.fft.fftn(kernel)
        denominator=np.fft.ifftn(np.fft.fftn(occupied)*weight).real
        envelope=np.fft.ifftn(np.fft.fftn(a*np.exp(-1j*phase))*weight)/denominator
        layer=shape[2]//2;theta=result['coherence']['phase_axis_radians']
        raw=(a*np.exp(-1j*theta)).real[:,:,layer]
        b=envelope[:,:,layer];amp=abs(a[:,:,layer]);mask=occupied[:,:,layer]
        # Compact visualization planes only; all physical summary/C data remain float64.
        np.savez_compressed(folder/(label+'_planes.npz'),quadrature=raw.astype(np.float32),
            envelope=b.astype(np.complex64),amplitude=amp.astype(np.float32),occupied=mask)
    raw_scale=float(abs(raw).max())
    normalized=raw/raw_scale if raw_scale else raw
    first=np.stack((1-np.maximum(-normalized,0),1-abs(normalized),1-np.maximum(normalized,0)),axis=-1)*255
    first[~mask]=0
    maximum=float(abs(b).max())
    bright=abs(b)/maximum if maximum else abs(b)
    hue=(np.angle(b)+np.pi)/(2*np.pi)
    second=np.array([colorsys.hsv_to_rgb(float(h),.9,float(v)) for h,v in zip(hue.ravel(),bright.ravel())]).reshape(hue.shape+(3,))*255
    amp_scale=float(amp.max())
    third=np.repeat((amp/amp_scale if amp_scale else amp)[:,:,None],3,axis=-1)*255
    third[~mask]=0
    png(folder/(label+'_phase.png'),np.concatenate((first,second,third),axis=1))
    caption=dict(case=label,plane='Cartesian(001) at fixed third primitive index k; raster rows increase along[110], columns along[-110]',
        primitive_basis_columns_fractional=CELL.tolist(),a_A=a_A,
        panels=['raw pump-axis real quadrature,red/blue; black vacancies','complex carrier-envelope hue with amplitude brightness','raw complex amplitude magnitude'],
        carrier_peak_index=kstar.tolist(),carrier_wavelength_nm=wavelength,
        visualization_Gaussian_sigma_nm=sigma,
        distinction='Carrier-envelope smoothing is for the picture only; correlation data are unfiltered. Repeated pi sign alternation of a coherent finite-q standing pattern is not a defect wall. Envelope texture is shown without assigning domain radii or classifying every phase gradient as a Z2 wall.')
    (folder/(label+'_phase_caption.json')).write_text(json.dumps(caption,indent=2)+'\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('folder');ap.add_argument('labels',nargs='+');ap.add_argument('--window',action='store_true')
    ap.add_argument('--phase',action='store_true');ap.add_argument('--a',type=float,default=3.732)
    ap.add_argument('--phase-data',action='store_true',help='render exported visualization planes, without per-site trajectories')
    ap.add_argument('--modulus',action='store_true')
    ap.add_argument('--name',default='correlation.svg')
    args=ap.parse_args();folder=Path(args.folder)
    if args.phase:
        for label in args.labels:phase_picture(folder,label,args.a,args.phase_data)
    else:correlation_svg(folder,folder/args.name,args.labels,args.window,args.modulus)
