"""Render the same-time and separately scaled-window comparisons from saved data."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent


def main():
    data=json.loads((HERE/'results.json').read_text())['runs']
    selected=[]
    for row in data:
        if 'failure' in row:continue
        c=row['config']
        if c['q']!=15 or c['N']!=2 or c['n']!=30000 or c.get('detuning'):continue
        if abs(c['g_over_quantum']/(.085*np.sqrt(.1))-1)>1e-10:continue
        if c['gamma_over_omega'] not in (0.,100.):continue
        selected.append(row)
    order={'none':0,'symmetric':1,'lower':2,'lower_limit':3}
    selected.sort(key=lambda r:order[r['config']['loss']])
    colors=['#233b6e','#7c5688','#ce6b2c','#168071']
    labels=['No loss','Symmetric, Γ/ω=100','Lower only, Γ/ω=100','Projected Γ→∞ limit']
    common=min(row['observation_limit'] for row in selected)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    for scaled_window,filename in ((False,'same_time.png'),(True,'scaled_windows.png')):
        fig,axes=plt.subplots(2,1,sharex=True,figsize=(8.2,6.3),layout='constrained')
        for row,color,label in zip(selected,colors,labels):
            times=np.asarray(row['up']['times']);prob=np.asarray(row['up']['target_unconditioned'])
            x=times/row['observation_limit'] if scaled_window else times
            mask=times>0 if scaled_window else times<=common
            legend=label+(f"; T={row['observation_limit']:.0f}" if scaled_window else '')
            axes[0].plot(x[mask],np.maximum(prob[mask],1e-30),color=color,label=legend)
            for direction,style in (('up','-'),('down','--')):
                sink=np.clip(row[direction]['sink'],0,1)
                axes[1].plot(x[mask],sink[mask],style,color=color,linewidth=1.3)
        axes[0].set_yscale('log');axes[0].set_ylim(1e-30,1)
        axes[0].set_ylabel('Unconditioned destination probability')
        axes[0].legend(loc='lower right',fontsize=8)
        axes[1].set_ylabel('Sink probability')
        axes[1].set_ylim(-.02,1.03)
        axes[1].text(.02,.08,'Sink: up solid; down dashed. Destination curves coincide.',transform=axes[1].transAxes,fontsize=8)
        if scaled_window:
            axes[1].set_xscale('log');axes[1].set_xlim(1e-8,1)
            axes[1].set_xlabel('ωt / T for each model (T differs between curves)')
            axes[0].set_title('Longer windows are separate observations, not a same-time comparison')
        else:
            axes[1].set_xlim(0,common);axes[1].set_xlabel('ωt, one common window')
            axes[0].set_title(f'15-quanta branch, N=2, n₀=30000; common ωt ≤ {common:.2f}')
        for ax in axes:ax.grid(alpha=.18)
        axes[0].text(.02,.03,'Display floor 10⁻³⁰; raw values retained.',transform=axes[0].transAxes,
            fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.85})
        fig.savefig(HERE/filename,dpi=160)
        plt.close(fig)
    print('Rendered same_time.png and scaled_windows.png from saved trajectories')


if __name__=='__main__':main()
