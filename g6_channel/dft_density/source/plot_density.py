"""Render the actual saved h.16 density and its source-bound string markers."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(readout,output):
    geometry=json.loads((readout/'geometry.json').read_bytes())
    with np.load(readout/'h016-transverse.npz',allow_pickle=False) as data:
        density=data['n_perp_e_A3'];mask=data['outside_Pd_augmentation']
    lx,ly=geometry['lengths_A'][:2]
    fig,ax=plt.subplots(figsize=(8,6),layout='constrained')
    image=ax.imshow(density.T,origin='lower',extent=(0,lx,0,ly),interpolation='nearest',cmap='viridis')
    x,y=np.meshgrid(np.arange(len(mask))/len(mask)*lx,np.arange(len(mask))/len(mask)*ly,indexing='ij')
    ax.contour(x,y,mask.astype(float),levels=[.5],colors=['white'],linewidths=.6,alpha=.7)
    for label,marker,face in [('Pd','o','white'),('D','s','#ffae42')]:
        rows=np.asarray(geometry['rows'][label]['xy_A'])
        ax.scatter(rows[:,0],rows[:,1],s=60,marker=marker,facecolors=face,edgecolors='black',label=label+' strings',zorder=3)
    centres=np.asarray(geometry['centres_xy_A'])
    ax.scatter(centres[:,0],centres[:,1],s=45,marker='+',color='black',label='Declared channel centres',zorder=4)
    ax.set(xlim=(0,lx),ylim=(0,ly),xlabel='[1−10] (Å)',ylabel='[001] (Å)',
           title='Pristine PdD ⟨110⟩ · Native valence line average\nh=0.16 Å; 48³ Native grid → 128² display, not a finer DFT grid')
    ax.legend(loc='upper right',fontsize=8)
    fig.colorbar(image,ax=ax,label='Valence density (electrons/Å³)')
    fig.savefig(output,dpi=160,metadata={'Software':'P255 source-bound PdD density readout'})
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--readout',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    render(args.readout,args.output)
