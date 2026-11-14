"""Phase correlations for the accepted lattice.py; no changed force or loss law.

Primary order parameter: pump-x complex half-pump demodulated displacement.
Raw amplitude-weighted and phase-only correlations use actual occupied pairs.
Connected correlations, finite-box coherent fractions and finite-q spectra are
reported separately. A normalized transient is not a surviving oscillation.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from lattice import CELL, NN, HC, HBAR, shell_kernel, physical_q, run_case, png


def pair_count(mask):
    f=np.fft.fftn(mask.astype(float))
    # This convolution is integer-valued; rounding recovers counts, not a floor.
    return np.rint(np.fft.ifftn(f.conj()*f).real).astype(np.int64)


def channel(field,mask):
    a=np.where(mask,field,0j)
    scale=float(np.max(np.abs(a)))
    counts=pair_count(mask)
    if scale==0:
        return None,counts,None
    a=a/scale
    g0=float(np.mean(np.abs(a[mask])**2))
    f=np.fft.fftn(a)
    numerator=np.fft.ifftn(f.conj()*f)
    c=np.full(a.shape,np.nan+1j*np.nan)
    np.divide(numerator,counts*g0,out=c,where=counts>0)
    spectrum=np.abs(f)**2
    peak=np.unravel_index(np.argmax(spectrum),spectrum.shape)
    denominator=int(mask.sum())*float(np.sum(np.abs(a)**2))
    info=dict(rms_A=scale*np.sqrt(g0),scale_A=scale,
        Gamma_coherent_fraction=float(spectrum[0,0,0]/denominator),
        peak_coherent_fraction=float(spectrum[peak]/denominator),
        peak_index=list(map(int,peak)),peak_spectral_share=float(spectrum[peak]/spectrum.sum()))
    return c,counts,info


def radial(c,counts,a_A):
    shape=counts.shape
    i=np.indices(shape)
    signed=np.array([np.where(i[d]<=shape[d]//2,i[d],i[d]-shape[d]) for d in range(3)])
    cart2=np.rint(2*np.einsum('ab,b...->a...',CELL,signed)).astype(np.int64)
    distance2=np.sum(cart2*cart2,axis=0)
    radius_nm=a_A*.05*np.sqrt(distance2)
    inscribed_nm=float(np.min(np.array(shape)/np.linalg.norm(np.linalg.inv(CELL),axis=1))*a_A*.05)
    valid=(radius_nm<inscribed_nm)&(counts>0)
    shell=distance2[valid]
    weight=counts[valid].astype(float)
    denominator=np.bincount(shell,weights=weight)
    numerator_r=np.bincount(shell,weights=weight*c[valid].real)
    numerator_i=np.bincount(shell,weights=weight*c[valid].imag)
    squared=np.bincount(shell,weights=weight*np.abs(c[valid])**2)
    keep=denominator>0
    radii=a_A*.05*np.sqrt(np.flatnonzero(keep))
    real=numerator_r[keep]/denominator[keep]
    imag=numerator_i[keep]/denominator[keep]
    rms=np.sqrt(squared[keep]/denominator[keep])
    # This is a descriptive initial decay distance, never an exponential fit.
    crossing=np.flatnonzero((radii>0)&(rms<np.exp(-1)))
    initial=None
    if crossing.size:
        j=int(crossing[0]); previous=max(j-1,0)
        if rms[previous]!=rms[j]:
            initial=float(radii[previous]+(radii[j]-radii[previous])*(rms[previous]-np.exp(-1))/(rms[previous]-rms[j]))
    tail=radii>.7*inscribed_nm
    return dict(radius_nm=radii.tolist(),C_real=real.tolist(),C_imag=imag.tolist(),
        direction_RMS_modulus=rms.tolist(),occupied_pairs=denominator[keep].astype(np.int64).tolist(),
        inscribed_box_radius_nm=inscribed_nm,initial_e_fold_distance_nm=initial,
        tail_pair_weighted_real=float(np.average(real[tail],weights=denominator[keep][tail])) if tail.any() else None,
        tail_RMS_modulus=float(np.sqrt(np.average(rms[tail]**2,weights=denominator[keep][tail]))) if tail.any() else None,
        interpretation='Initial e-fold distance is not a fitted correlation length if a plateau, oscillation or directional long range survives.')


def directions(c,counts,a_A,maximum_nm):
    basis=[('NN_'+str(tuple(map(int,v))),v) for v in NN if next(x for x in v if x)!=-1]
    basis += [(name,np.rint(np.linalg.inv(CELL)@axis).astype(int)) for name,axis in
              (('100',(1,0,0)),('010',(0,1,0)),('001',(0,0,1)))]
    answer={}
    for name,v in basis:
        step=float(np.linalg.norm(CELL@v)*a_A*.1)
        n=np.arange(int(maximum_nm/step)+1)
        idx=tuple((v[d]*n)%c.shape[d] for d in range(3))
        values=c[idx]; p=counts[idx]
        answer[name]=dict(radius_nm=(n*step).tolist(),occupied_pairs=p.tolist(),
            C_real=[float(z.real) if count else None for z,count in zip(values,p)],
            C_imag=[float(z.imag) if count else None for z,count in zip(values,p)],
            modulus=[float(abs(z)) if count else None for z,count in zip(values,p)])
    return answer


def characterize(amplitude,occupied,a_A):
    scalar=amplitude[0]
    c,counts,weighted=channel(scalar,occupied)
    if weighted is None:
        return dict(signal_defined=False,rms_A=0.,phase_only=None),{}
    defined=occupied&(np.abs(scalar)>0)
    z=np.zeros(scalar.shape,complex)
    z[defined]=scalar[defined]/np.abs(scalar[defined])
    cp,phase_counts,phase_info=channel(z,defined)
    phase_info['unit_phasor_rms']=phase_info.pop('rms_A')
    phase_info.pop('scale_A')
    weighted['peak_q']=physical_q(weighted['peak_index'],scalar.shape,a_A)
    weighted['radial']=radial(c,counts,a_A)
    weighted['directions']=directions(c,counts,a_A,weighted['radial']['inscribed_box_radius_nm'])
    phase_info['radial']=radial(cp,phase_counts,a_A)
    phase_info['directions']=directions(cp,phase_counts,a_A,weighted['radial']['inscribed_box_radius_nm'])
    grid=np.indices(scalar.shape)
    phase=np.exp(-2j*np.pi*sum(weighted['peak_index'][d]*grid[d]/scalar.shape[d] for d in range(3)))
    weighted['carrier_removed_radial']=radial(c*phase,counts,a_A)
    phase_info['carrier_removed_radial']=radial(cp*phase,phase_counts,a_A)
    scale=weighted['scale_A']; scaled=scalar/scale
    mean=scaled[occupied].mean()
    if np.ptp(scaled[occupied].real)==0 and np.ptp(scaled[occupied].imag)==0:
        mean=scaled[occupied][0]
    centered=np.where(occupied,scaled-mean,0j)
    cc,pc,connected=channel(centered,occupied)
    if connected is not None:
        connected['radial']=radial(cc,pc,a_A)
        connected['variance_fraction_of_total']=float(np.mean(abs(centered[occupied])**2)/np.mean(abs(scaled[occupied])**2))
        connected['integral_correlation']=integral_correlation(cc,pc,weighted['peak_index'],a_A)
        connected.pop('rms_A');connected.pop('scale_A')
    second=np.mean(z[defined]**2)
    weighted_second=np.sum(scaled[occupied]**2)/np.sum(abs(scaled[occupied])**2)
    axis=float(np.angle(weighted_second)/2) if weighted_second!=0 else None
    histogram=None; wall_fraction=None; sign=None
    if axis is not None:
        angle=np.angle(z[defined]*np.exp(-1j*axis))
        edges=np.linspace(-np.pi,np.pi,37)
        count,_=np.histogram(angle,bins=edges)
        weight,_=np.histogram(angle,bins=edges,weights=abs(scaled[defined])**2)
        histogram=dict(edges_radians=edges.tolist(),phase_only_count=count.tolist(),amplitude_squared_weight=weight.tolist())
        projection=(scalar*np.exp(-1j*axis)).real
        sign=np.sign(projection)
        eligible=defined&(projection!=0)
        walls=0; bonds=0
        for v in NN:
            if next(x for x in v if x)<0:
                continue
            neighbor=np.roll(eligible,tuple(v),axis=(0,1,2))
            paired=eligible&neighbor
            opposite=sign*np.roll(sign,tuple(v),axis=(0,1,2))<0
            walls+=int(np.sum(paired&opposite));bonds+=int(paired.sum())
        wall_fraction=walls/bonds if bonds else None
    quantiles=np.quantile(np.abs(scalar[occupied]),[0,.1,.5,.9,1])
    weighted['integral_correlation']=integral_correlation(c,counts,weighted['peak_index'],a_A)
    phase_info['integral_correlation']=integral_correlation(cp,phase_counts,weighted['peak_index'],a_A)
    result=dict(signal_defined=True,component='pump-x',rms_A=weighted['rms_A'],
        amplitude_quantiles_A=quantiles.tolist(),amplitude_weighted=weighted,phase_only=phase_info,
        connected_amplitude=connected,Z2_phase_only_second_moment=float(abs(second)),
        Z2_amplitude_weighted_second_moment=float(abs(weighted_second)),phase_axis_radians=axis,
        phase_histogram=histogram,opposite_projected_phase_NN_bond_fraction=wall_fraction,
        phase_note='All exactly nonzero occupied amplitudes are included, without an amplitude or ratio floor. A phase of decaying seed/noise is not sustained order. Projected-sign walls are Z2 walls only when the phase histogram is bimodal.')
    arrays=dict(amplitude_x_A=scalar,occupied=occupied,C_amplitude=c,C_phase=cp,
                occupied_pairs=counts,phase_pairs=phase_counts,phase_sign=sign)
    return result,arrays


def integral_correlation(c,counts,carrier_index,a_A):
    """Full unfiltered pair-normalized C; preserve every coherent plateau."""
    if np.any(counts==0):
        return dict(defined=False,reason='Some displacement pair counts vanish; no full-box integral is assigned.')
    value=np.fft.fftn(c)[tuple(carrier_index)]*(a_A*.1)**3/4
    volume=float(value.real)
    return dict(defined=True,carrier_index=list(map(int,carrier_index)),
        correlation_volume_nm3=volume,imaginary_roundoff_nm3=float(value.imag),
        xi_integral_nm=float(np.cbrt(volume/(8*np.pi))) if volume>0 else None,
        definition='Vq=v_site*sum_d C(d)*exp(-iq_carrier.d), xi_integral=(Vq/(8pi))^(1/3). This convention equals an isotropic exponential decay length only for that shape; no exponential was fitted, no plateau removed, and this is not a domain radius. A plateau makes the integral grow with box size.')


class PhaseWindows:
    """Equal-time envelope correlations, not the phase of a long-time mean.

    Each local envelope is an ordinary least-squares cosine/sine/DC fit over
    the declared number of carrier cycles. Only moments/correlations survive;
    no trajectory is retained. This callback does not alter the dynamics.
    """
    def __init__(self,frequency_THz,samples,a_A):
        self.frequency=frequency_THz
        self.samples=samples
        self.a_A=a_A
        self.count=0
        self.windows=0
        self.gram=np.zeros((3,3))
        self.sums=None
        self.energy=0.
        self.power=None
        self.phase_power=None
        self.phase_pairs=None
        self.R2=[]
        self.times=[]

    def __call__(self,q_A,occupied,time_ps):
        if self.sums is None:
            self.mask=occupied.copy()
            self.sums=np.zeros((3,)+occupied.shape)
            self.power=np.zeros(occupied.shape)
            self.phase_power=np.zeros(occupied.shape)
            self.phase_pairs=np.zeros(occupied.shape,np.int64)
        phase=2*np.pi*self.frequency*time_ps
        basis=np.array([np.cos(phase),np.sin(phase),1.])
        self.sums+=basis[:,None,None,None]*q_A[0]
        self.gram+=np.outer(basis,basis)
        self.count+=1
        if self.count==self.samples:
            coefficients=np.linalg.solve(self.gram,self.sums.reshape(3,-1)).reshape(self.sums.shape)
            a=coefficients[0]-1j*coefficients[1]
            f=np.fft.fftn(a)
            self.power+=abs(f)**2
            self.energy+=float(np.mean(abs(a[self.mask])**2))
            defined=self.mask&(abs(a)>0)
            z=np.zeros(a.shape,complex)
            z[defined]=a[defined]/abs(a[defined])
            self.phase_power+=abs(np.fft.fftn(z))**2
            self.phase_pairs+=pair_count(defined)
            self.R2.append(float(abs(np.mean(z[defined]**2))) if defined.any() else None)
            self.windows+=1
            self.times.append(time_ps)
            self.sums.fill(0);self.gram.fill(0);self.count=0

    def finish(self):
        if not self.windows or self.energy==0:
            return dict(defined=False,complete_windows=self.windows,partial_samples=self.count)
        pairs=pair_count(self.mask)
        ca=np.full(pairs.shape,np.nan+1j*np.nan)
        cp=np.full(pairs.shape,np.nan+1j*np.nan)
        np.divide(np.fft.ifftn(self.power),pairs*self.energy,out=ca,where=pairs>0)
        np.divide(np.fft.ifftn(self.phase_power),self.phase_pairs,out=cp,where=self.phase_pairs>0)
        peak=np.unravel_index(np.argmax(self.power),self.power.shape)
        denominator=int(self.mask.sum())**2*self.energy
        phase_radial=radial(cp,self.phase_pairs,self.a_A)
        return dict(defined=True,complete_windows=self.windows,partial_samples=self.count,
            samples_per_window=self.samples,window_end_times_ps=self.times,
            frequency_THz=self.frequency,rms_A=float(np.sqrt(self.energy/self.windows)),
            Gamma_coherent_fraction=float(self.power[0,0,0]/denominator),
            peak_coherent_fraction=float(self.power[peak]/denominator),
            peak_q=physical_q(peak,pairs.shape,self.a_A),
            phase_R2_per_window=self.R2,
            amplitude_radial=radial(ca,pairs,self.a_A),
            phase_radial=phase_radial,
            phase_directions=directions(cp,self.phase_pairs,self.a_A,phase_radial['inscribed_box_radius_nm']),
            meaning='Temporal mean of equal-time short-envelope pair correlations. Cosine/sine/DC least-squares envelopes; partial final window excluded, no trajectory saved. Long-time pump-locked-mean C is separate.')


def linear_guide(p,shape,h,ratio,dt):
    A=p['K4_100_eV_A4'];B=(3*p['K4_111_eV_A4']-A)/2
    K=p['Kii_eV_A2']+(3*A+2*B)*p['zero_point_axis_variance_A2']
    E=np.sqrt(HC*HC/p['mass_eVc2']*p['Kii_eV_A2'])
    factor=HC*HC/(p['mass_eVc2']*E*E)
    J=np.fft.fftn(shell_kernel(p,shape)).real
    w=np.sqrt((K-J)*factor)
    wref=np.sqrt((K-sum(s['J_sum_eV_A2'] for s in p['shells']))*factor)
    gamma=2/(p['reference_amplitude_decay_ps']*1e-12)/(E/HBAR)
    argument=np.cosh(gamma*dt/2)*(1-.5*dt*dt*w*w)
    if np.any(abs(argument)>1):
        raise ValueError('Linear guide assumes underdamped stable discrete unpumped modes')
    numerical_frequency=np.arccos(argument)/dt
    detuning=numerical_frequency-wref*ratio
    coupling=wref*wref/(4*w)
    thresholds=np.sqrt((gamma/2)**2+detuning**2)/coupling
    gain=-gamma/2+np.sqrt(np.maximum((h*coupling)**2-detuning**2,0))
    return dict(weak_modulation_guide_only=True,minimum_h=float(thresholds.min()),
        maximum_growth_per_ps=float(gain.max()*E/HBAR*1e-12),
        physical_elapsed_ps_per_step=float(dt*HBAR/E*1e12),
        pump_half_frequency_THz=float(wref*ratio*E/HBAR/(2*np.pi)*1e-12),
        description='First-order Mathieu guide, using accepted dt discrete unpumped frequencies; not an exact nonlinear or thermal onset.')


def known_case():
    shape=(8,8,8);i,j,k=np.indices(shape)
    mask=(i<3)&((2*j+k)%3!=0)
    a=mask*np.exp(2j*np.pi*j/shape[1])
    c,p,info=channel(a,mask)
    expected=np.exp(2j*np.pi*j/shape[1])
    error=float(np.max(abs(c[p>0]-expected[p>0])))
    wrong=np.fft.ifftn(abs(np.fft.fftn(a))**2)/mask.sum()
    wrong_error=float(np.max(abs(wrong[p>0]-expected[p>0])))
    if error>1e-12 or wrong_error<.1 or not np.any(p==0):
        raise ValueError('Occupied-pair, finite-q or undefined-pair control failed')
    return dict(pair_normalized_plane_wave_error=error,wrong_zero_filled_error=wrong_error,
                undefined_displacements=int(np.sum(p==0)),finite_q_coherent_fraction=info['peak_coherent_fraction'])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('parameters');ap.add_argument('cases');ap.add_argument('--case');ap.add_argument('--out',default='coherence_results')
    ap.add_argument('--known-case',action='store_true')
    ap.add_argument('--window-cycles',type=int)
    args=ap.parse_args()
    if args.known_case:
        print(json.dumps(known_case(),indent=2));return
    parameters=json.loads(Path(args.parameters).read_text())
    cases=json.loads(Path(args.cases).read_text())['cases']
    folder=Path(args.out);folder.mkdir(exist_ok=True)
    for c in cases:
        if args.case and c['label']!=args.case:
            continue
        started=perf_counter()
        p=parameters['materials'][c.get('material','NiH')]
        extent=c.get('L',24);shape=(extent,)*3 if isinstance(extent,int) else tuple(extent)
        guide=linear_guide(p,shape,c.get('pump_h',0),c.get('pump_ratio',1),c.get('dt',.15))
        observer=None;epoch_observers=None;epoch_bounds=None
        if args.window_cycles:
            dt_sample=guide['physical_elapsed_ps_per_step']*c.get('stride',2)
            samples=round(args.window_cycles/(guide['pump_half_frequency_THz']*dt_sample))
            steps=c.get('phase_epoch_steps')
            if steps:
                epoch_bounds=np.array(steps)*guide['physical_elapsed_ps_per_step']
                epoch_observers=[PhaseWindows(guide['pump_half_frequency_THz'],samples,p['a_A']) for _ in steps[:-1]]
                def observer(q,occupied,t):
                    for index,(low,high) in enumerate(zip(epoch_bounds[:-1],epoch_bounds[1:])):
                        if low<=t<high:
                            epoch_observers[index](q,occupied,t)
                            break
            else:
                observer=PhaseWindows(guide['pump_half_frequency_THz'],samples,p['a_A'])
        summary,fields=run_case(parameters,c,observe=observer)
        correlation,arrays=characterize(fields['driven_amplitude_A'],fields['occupied'],p['a_A'])
        summary['coherence']=correlation
        summary['linear_guide']=guide
        if epoch_observers is not None:
            summary['windowed_coherence_epochs']=[dict(interval_ps=[float(epoch_bounds[j]),float(epoch_bounds[j+1])],**o.finish()) for j,o in enumerate(epoch_observers)]
        elif observer is not None:
            summary['windowed_coherence']=observer.finish()
        summary['coherence_total_elapsed_s']=perf_counter()-started
        (folder/(c['label']+'.json')).write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
        if observer is None:
            np.savez_compressed(folder/(c['label']+'.npz'),**arrays)
        print(json.dumps({k:v for k,v in summary.items() if k not in ('spectrum','coherence')}),flush=True)
        print(json.dumps(dict(label=c['label'],rms_x_A=correlation['rms_A'],
            phase_R2=correlation.get('Z2_phase_only_second_moment'),
            Gamma_fraction=correlation['amplitude_weighted']['Gamma_coherent_fraction'],
            phase_initial_efold_nm=correlation['phase_only']['radial']['initial_e_fold_distance_nm'],
            phase_tail_RMS=correlation['phase_only']['radial']['tail_RMS_modulus'])),flush=True)


if __name__=='__main__':
    main()
