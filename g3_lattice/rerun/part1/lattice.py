"""One-off FCC displacement dynamics; NumPy + Python standard library only.
CLI: python lattice.py lattice_inputs.json --inspect
     python lattice.py lattice_inputs.json --batch cases.json --out results
Fields use layer-adapted FCC primitive coordinates; third index is (001) depth.
Only zero-point Hartree is added; sampled classical thermal noise is not added
again as a thermal Hartree shift. Strain is an explicit diagonal coarse map.
"""
import argparse
import itertools
import json
from math import cos, pi, sqrt
from pathlib import Path
import struct
import time
import zlib
import numpy as np
from measured_environment import load_measured_environment, normalized_bias
from hashlib import sha256
from vacancy_bias import predict_reference_bias

HBAR=6.582119569e-16
HC=1973.269804
KB=8.617333262e-5
CELL=.5*np.array([[1,-1,1],[1,1,0],[0,0,1]],float)  # Cartesian r/a=CELL@index
NN=np.array([(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1),
             (1,-1,-1),(-1,1,1),(0,1,1),(0,-1,-1),(1,0,-1),(-1,0,1)])


def shell_kernel(p,shape,unit_shell=None):
    kernel=np.zeros(shape)
    for number,shell in enumerate(p['shells']):
        if unit_shell is not None and number!=unit_shell:
            continue
        aliases=shell['periodic_alias_offsets_A']
        value=(1. if unit_shell is not None else shell['J_sum_eV_A2'])/len(aliases)
        for raw in aliases:
            r=np.rint(2*np.array(raw)/p['a_A']).astype(int)
            choices=[(-2,2) if abs(v)==2 else (int(v),) for v in r]
            images=list(itertools.product(*choices))
            for x,y,z in images:
                uvw=((x+y-z)//2,(y-x+z)//2,z)
                kernel[tuple(v%shape[d] for d,v in enumerate(uvw))]+=value/len(images)
    inverse=kernel[np.ix_(*[(-np.arange(n))%n for n in shape])]
    return (kernel+inverse)/2


def physical_q(index,shape,a):
    f=np.array([v if v<=shape[d]//2 else v-shape[d] for d,v in enumerate(index)],float)/shape
    transform=np.linalg.inv(CELL).T*2*pi/a
    candidates=[transform@(f+v) for v in itertools.product((-1,0,1),repeat=3)]
    q=min(candidates,key=lambda v:np.dot(v,v))
    norm=np.linalg.norm(q)
    return dict(layer_FCC_index=list(map(int,index)),q_cartesian_A_inverse=q.tolist(),
                wavelength_nm=float(2*pi/norm/10) if norm else None)


def occupation(shape,x,kind,rng):
    n=int(np.prod(shape))
    count=round(x*n)
    if not 0<count<=n:
        raise ValueError('Loading must give occupied sites')
    if kind=='random':
        score=rng.random(shape)
    elif kind=='clustered':
        noise=rng.normal(size=shape)
        near=np.zeros(shape)
        for r in NN:
            near[tuple(r%shape)]+=1
        lap=12-np.fft.rfftn(near).real
        score=np.fft.irfftn(np.fft.rfftn(noise)*np.exp(-.5*9*lap/12),s=shape,axes=(0,1,2))
    elif kind=='ordered':
        i,j,k=np.indices(shape)
        rank=((i+j+k)%2)*4+((j+k)%2)*2+k%2
        score=-rank+rng.random(rank.shape)*.01
    else:
        raise ValueError('Unknown arrangement')
    chosen=np.argpartition(score.ravel(),n-count)[n-count:]
    occupied=np.zeros(n,bool)
    occupied[chosen]=True
    return occupied.reshape(shape)


def environment_law(p,vacant_nn,strain):
    """Frozen-x1 branch for hosts without an adopted loading law.
    K contains T0 Hartree only; its declared linear strain response is retained.
    The PdD864 branch below instead uses its actual loading law and PartA contrast.
    """
    A=p['K4_100_eV_A4']
    B=(3*p['K4_111_eV_A4']-A)/2
    zp=(3*A+2*B)*p['zero_point_axis_variance_A2']
    K=p['Kii_eV_A2']+zp+p['collective_strain_slope_eV_A2']*strain
    return K,A,B,None  # None selects the actually supplied frozen-x1 J kernel.


def coupling_force(state,Jfft,fftshape,J_fields=None,unit_ffts=None):
    """Reciprocal bond J_ij=(J_i,s+J_j,s)/2 times each unit shell kernel."""
    axes=(-3,-2,-1)
    ft=np.fft.rfftn(state,s=fftshape,axes=axes)
    if J_fields is None:
        return np.fft.irfftn(ft*Jfft,s=fftshape,axes=axes)[...,:state.shape[-1]]
    total=np.zeros_like(state)
    for field,unit in zip(J_fields,unit_ffts):
        neighbour=np.fft.irfftn(ft*unit,s=fftshape,axes=axes)[...,:state.shape[-1]]
        weighted=np.fft.rfftn(field*state,s=fftshape,axes=axes)
        other=np.fft.irfftn(weighted*unit,s=fftshape,axes=axes)[...,:state.shape[-1]]
        total+=.5*(field*neighbour+other)
    return total


def power_spectrum(series,dt_ps):
    y=np.asarray(series)
    y=y-y.mean(axis=0,keepdims=True)
    window=np.hanning(len(y))
    ft=np.fft.rfft(y*window.reshape((-1,)+(1,)*(y.ndim-1)),axis=0)
    power=dt_ps*abs(ft)**2/np.dot(window,window)
    if len(power)>2:
        power[1:-1]*=2
    if power.ndim>1:
        power=power.sum(axis=-1)
    return np.fft.rfftfreq(len(y),dt_ps),power


def spatial(field,occupied):
    n=occupied.sum()
    mean=field[:,occupied].mean(axis=1)
    connected=(field-mean[:,None,None,None])*occupied
    ft=np.fft.fftn(connected,axes=(-3,-2,-1))
    S=(abs(ft)**2).sum(axis=0)/n
    S[0,0,0]=0
    total=float(S.sum())
    index=np.unravel_index(np.argmax(S),S.shape)
    density=np.fft.fftn(occupied-occupied.mean())
    Sn=abs(density)**2/occupied.size
    Sn[0,0,0]=0
    return S,Sn,index,float(S[index]/total) if total>0 else 0.


def harmonic_projection(sumq,sumamp,sumpump,gram):
    """Joint DC,Omega,2Omega fit; forced pump-frequency motion is not a subharmonic."""
    if gram[4,4]<5 or np.linalg.matrix_rank(gram)!=5:
        raise ValueError('Record cannot separate DC, fundamental and pump-frequency tones')
    weights=np.linalg.solve(gram,np.eye(5))
    # Basis:cos(theta),sin(theta),cos(2theta),sin(2theta),DC.
    # Imaginary Fourier sums have the opposite sign to the sine correlations.
    correlations=(sumamp.real,sumamp.imag,sumpump.real,sumpump.imag)
    signs=(1.,-1.,1.,-1.)
    scratch=np.empty_like(sumq)
    def component(row,out):
        np.multiply(sumq,weights[row,4],out=out)
        for column,(value,sign) in enumerate(zip(correlations,signs)):
            np.multiply(value,weights[row,column]*sign,out=scratch)
            np.add(out,scratch,out=out)
    dc=np.empty_like(sumq);amplitude=np.empty(sumq.shape,complex);pump=np.empty(sumq.shape,complex)
    component(4,dc);component(0,amplitude.real);component(1,amplitude.imag)
    component(2,pump.real);component(3,pump.imag)
    amplitude.imag*=-1;pump.imag*=-1
    return dc,amplitude,pump


OBSERVATION_DEFAULTS={'dt':.15,'steps':4096,'stride':2}


def run_case(parameters,c,observe=None):
    started=time.perf_counter()
    from loading_response import bulk_reference,response_fields,elastic_tensor
    material=c.get('material','NiH')
    p=parameters['materials'][material]
    loading_spec=parameters['loading_law']
    use_loading=loading_spec['kind']=='anchored_directional_affine_v2' and material in loading_spec['hosts']
    extrapolate=c.get('loading_extrapolate',False)
    composition=bulk_reference(c.get('x',1.),extrapolate=extrapolate) if use_loading else None
    extent=c.get('L',24)
    shape=(extent,)*3 if isinstance(extent,int) else tuple(extent)
    n=int(np.prod(shape))
    if not 10000<=n<=1000000:
        raise ValueError('Assigned size range is 1e4 through1e6 sites')
    open_z=c.get('boundary_z','periodic')=='open'
    fftshape=shape[:2]+(shape[2]+4,) if open_z else shape
    rng=np.random.default_rng(c.get('seed',13))
    physical_bias=None;bias_metadata=None
    if 'measured_environment_input' in c:
        if c.get('arrangement')!='measured_periodic' or 'x' not in c:
            raise ValueError('Measured bias requires its exact periodic occupancy and explicit loading')
        occ,physical_bias,_,bias_metadata=load_measured_environment(
            c['measured_environment_input'],shape,c['x'],CELL,c.get('bias_lattice_shift',(0,0,0)))
    else:
        occ=occupation(shape,c.get('x',1.),c.get('arrangement','random'),rng)
    Nocc=int(occ.sum())
    i,j,k=np.indices(shape)
    geometry_a=c.get('geometry_a_A',composition['a_angstrom'] if composition is not None else p['a_A'])
    if bias_metadata is not None and abs(geometry_a-bias_metadata['registered_a_A'])>1e-10:
        if c.get('bias_geometry_policy')!='frozen_reference_under_strain':
            raise ValueError('Transfer of a bulk bias to another geometry must be explicit')
    if 'bias_law_input' in c:
        law_path=Path(c['bias_law_input']);law_bytes=law_path.read_bytes()
        exact_occupancy_source=bias_metadata
        physical_bias,bias_metadata=predict_reference_bias(
            occ,CELL,json.loads(law_bytes),c.get('bias_shell2_multiplier',0.))
        bias_metadata.update(law_file=law_path.name,law_SHA256=sha256(law_bytes).hexdigest(),
            exact_occupancy_template=exact_occupancy_source)
    knm=2*pi/(shape[0]*geometry_a/sqrt(2)*.1)
    strain=c.get('strain',0.)+c.get('strain_gradient_per_nm',0.)/knm*np.sin(2*pi*i/shape[0])
    vacant=np.zeros(shape,int)
    # The calibrated vacancy descriptor is bulk-periodic occupation.
    # Open depth clamps exterior displacements; it does not add vacancies.
    for r in NN:
        neighbour=np.roll(occ,tuple(r),axis=(0,1,2))
        vacant+=~neighbour
    loading_metadata=None;elastic_metadata=None
    if use_loading:
        K,A,B,J_fields,loading_metadata=response_fields(p,occ,vacant,loading_spec['input'],
            composition['a_angstrom'],extrapolate=extrapolate)
    else:
        K,A,B,J_fields=environment_law(p,vacant,strain)
    tensor=None
    if 'tensor_input' in c and c['tensor_kind']!='delta_collective_curvature':
        raise ValueError('Actual map meaning must be resolved before replacing stiffness')
    if use_loading:
        tensor,elastic_metadata=elastic_tensor(strain,composition['a_angstrom'],loading_spec['strain_law_input'],
            map_path=c.get('tensor_input'),extrapolate=c.get('elastic_extrapolate',False))
    elif 'tensor_input' in c:
        tensor=np.load(c['tensor_input'])['delta_K_eV_A2']
    if tensor is not None and tensor.shape!=(3,3)+shape:
        raise ValueError('Tensor map must be Cartesian3x3 on the actual FCC site array')
    if 'control_tensor_eV_A2' in c:
        tensor=np.array(c['control_tensor_eV_A2'])[:,:,None,None,None]
    if tensor is not None and not np.allclose(tensor,np.swapaxes(tensor,0,1),rtol=0,atol=1e-10):
        raise ValueError('Non-symmetric stiffness is not a Hamiltonian Hessian')
    Eref=sqrt(HC*HC/p['mass_eVc2']*p['Kii_eV_A2'])
    length=HC/sqrt(p['mass_eVc2']*Eref)
    factor=length*length/Eref
    bias_force=normalized_bias(physical_bias,length,Eref) if physical_bias is not None else None
    Jfft=np.fft.rfftn(shell_kernel(p,fftshape)*factor)
    if np.max(abs(Jfft.imag))>1e-10:
        raise ValueError('Non-reciprocal kernel')
    unit_ffts=None
    if J_fields is not None:
        if not np.isrealobj(J_fields) or not np.isfinite(J_fields).all():
            raise ValueError('Static shell couplings must be finite and real')
        J_fields=np.asarray(J_fields)*factor
        if J_fields.shape!=(len(p['shells']),)+shape:
            raise ValueError('Environment-law shell fields must be(shell,x,y,z)')
        unit_ffts=[np.fft.rfftn(shell_kernel(p,fftshape,s)) for s in range(len(p['shells']))]
    kr=K*factor
    tr=tensor*factor if tensor is not None else None
    ar=A*length**4/Eref
    br=B*length**4/Eref
    baseK=p['Kii_eV_A2']+(3*A+2*B)*p['zero_point_axis_variance_A2']
    GammaK=baseK-sum(s['J_sum_eV_A2'] for s in p['shells'])
    omega_ref=sqrt(GammaK*factor)
    gamma=2/(p['reference_amplitude_decay_ps']*1e-12)/(Eref/HBAR)
    theta=KB*c.get('temperature_K',0.)/Eref
    h=c.get('pump_h',0.)
    drive=2*omega_ref*c.get('pump_ratio',1.)
    dt=c.get('dt',OBSERVATION_DEFAULTS['dt'])
    steps=c.get('steps',OBSERVATION_DEFAULTS['steps'])
    stride=c.get('stride',OBSERVATION_DEFAULTS['stride'])
    if c.get('harmonic_control',False):
        ar=br=0.
        gamma=theta=h=0.
    q=rng.normal(0,.01,size=(3,)+shape)*occ
    v=np.zeros_like(q)
    mode_vector=np.array([1.,0.,0.])
    control_w2=GammaK*factor
    if c.get('harmonic_control',False):
        mat=np.eye(3)*GammaK
        if tensor is not None:
            mat+=tensor[:,:,0,0,0]
        eigen,vectors=np.linalg.eigh(mat*factor)
        control_w2=eigen[0];mode_vector=vectors[:,0]
        q[:]=.01*mode_vector[:,None,None,None]*occ
    decay=np.exp(-gamma*dt)
    noise=sqrt(theta*(1-decay*decay))
    probes=rng.choice(np.flatnonzero(occ),size=min(24,Nocc),replace=False)
    sumq=np.zeros_like(q);sumamp=np.zeros(q.shape,complex);sumpump=np.zeros(q.shape,complex)
    samples=0;collective=[];local=[];control=[]
    demodulation_gram=np.zeros((5,5))
    pump_direction=None
    if 'pump_polarization' in c:
        pump_direction=np.array(c['pump_polarization'],float)
        pump_direction/=np.linalg.norm(pump_direction)
    def force(state,t):
        linear=coupling_force(state,Jfft,fftshape,J_fields,unit_ffts)
        squared=state*state
        r2=squared.sum(axis=0)
        f=linear-(kr+br*r2+(ar-br)*squared)*state
        if tr is not None:
            f-=np.einsum('ab...,b...->a...',tr,state)
        if h:
            if pump_direction is None:
                f[0]-=h*omega_ref**2*cos(drive*t)*state[0]
            else:
                projection=np.einsum('a,a...->...',pump_direction,state)
                f-=h*omega_ref**2*cos(drive*t)*pump_direction[:,None,None,None]*projection
        if bias_force is not None:
            f+=bias_force
        return f*occ
    F=force(q,0.)
    for step in range(steps):
        v+=.5*dt*F;q+=.5*dt*v
        v*=decay
        if noise:
            v+=noise*rng.normal(size=v.shape)*occ
        q+=.5*dt*v
        F=force(q,(step+1)*dt);v+=.5*dt*F
        if not np.isfinite(q).all():
            raise FloatingPointError('Nonfinite dynamics; result not usable')
        if c.get('harmonic_control',False):
            control.append(float(np.dot(mode_vector,q.reshape(3,-1).mean(axis=1))))
        if step>=steps//2 and step%stride==0:
            t=(step+1)*dt
            phasor=np.exp(-.5j*drive*t)
            pump_phasor=phasor*phasor
            sumq+=q;sumamp+=q*phasor;sumpump+=q*pump_phasor;samples+=1
            trig=np.array([phasor.real,-phasor.imag,pump_phasor.real,-pump_phasor.imag,1.])
            demodulation_gram+=np.outer(trig,trig)
            collective.append(q.sum(axis=(1,2,3))/Nocc*length)
            local.append(q.reshape(3,-1)[:,probes].T*length)
        if observe is not None and step>=c.get('observe_start_step',steps//2) and step%stride==0:
            t=(step+1)*dt
            observe(q*length,occ,t*HBAR/Eref*1e12)
    dc,amplitude,pump_amplitude=harmonic_projection(sumq,sumamp,sumpump,demodulation_gram)
    meanq=dc*length;amp=amplitude*length;amp_pump=pump_amplitude*length
    static_rms=float(np.sqrt((meanq[:,occ]**2).sum(axis=0).mean()))
    driven_rms=float(np.sqrt((abs(amp[:,occ])**2).sum(axis=0).mean()))
    pump_rms=float(np.sqrt((abs(amp_pump[:,occ])**2).sum(axis=0).mean()))
    field_kind='driven_amplitude' if h else 'static_mean'
    field=amp if field_kind=='driven_amplitude' else meanq
    S,Sn,peak,share=spatial(field,occ)
    magnitudes=np.sqrt((abs(field)**2).sum(axis=0))
    active=float((magnitudes[occ]>.05).mean())
    dt_ps=dt*stride*HBAR/Eref*1e12
    freq,ps=power_spectrum(np.array(collective),dt_ps)
    _,local_power=power_spectrum(np.array(local),dt_ps)
    local_power=local_power.mean(axis=-1)
    fmode=drive/2*Eref/HBAR/(2*pi*1e12)
    lines={}
    for name,ratio in [('half_pump',.5),('pump',1.),('third_half_pump',1.5),('twice_pump',2.)]:
        target=2*fmode*ratio;band=abs(freq-target)<=2*(freq[1]-freq[0])
        covered=(target-2*(freq[1]-freq[0])>=0 and target+2*(freq[1]-freq[0])<=freq[-1])
        lines[name]={'frequency_THz':target,'band_within_sampled_spectrum':bool(covered),
                     'collective_band_power_A2':float(ps[band].sum()*(freq[1]-freq[0])) if covered else None,
                     'mean_local_band_power_A2':float(local_power[band].sum()*(freq[1]-freq[0])) if covered else None}
    state='disordered'
    if static_rms>.05:
        state='vacancy_driven_displacement' if c.get('x',1.)<1 else 'static_displacement'
    if h and driven_rms>.05:state='driven_oscillation'
    if active>.05 and share>.02:state='patterned_'+state
    if c.get('scalability_only',False):state='short_scalability_record_not_phase_classified'
    corr=None
    if np.std(vacant[occ])>0 and np.std(magnitudes[occ])>0:
        corr=float(np.corrcoef(vacant[occ],magnitudes[occ])[0,1])
    layer_order=[];layer_fraction=[]
    for layer in range(shape[2]):
        mask=occ[:,:,layer]
        layer_order.append(float(np.sqrt(np.mean(magnitudes[:,:,layer][mask]**2))) if mask.any() else None)
        layer_fraction.append(float(np.mean(magnitudes[:,:,layer][mask]>.05)) if mask.any() else None)
    summary=dict(label=c['label'],configuration=c,sites=n,occupied=Nocc,actual_loading=Nocc/n,
        elapsed_s=time.perf_counter()-started,state=state,
        classification='finite-window displacement/oscillation amplitude diagnostic; partial-loading static signal is vacancy-driven displacement, not spontaneous ordering; extra stretch/pump ordering needs its matched baseline',
        static_rms_A=static_rms,driven_demodulated_rms_A=driven_rms,pump_frequency_rms_A=pump_rms,active_fraction=active,
        order_field=field_kind,
        displacement_peak=physical_q(peak,shape,geometry_a),peak_fraction=share,
        occupation_peak=physical_q(np.unravel_index(np.argmax(Sn),Sn.shape),shape,geometry_a),
        vacancy_vs_amplitude_correlation=corr,layer_order_amplitude_A=layer_order,
        layer_active_fraction=layer_fraction,
        stiffness_modulation_THz=2*fmode,incident_E_squared_carrier_THz=fmode,
        finite_record_resolution_THz=float(freq[1]-freq[0]),
        spectral_lines=lines,spectrum={'frequency_THz':freq.tolist(),'collective_PSD_A2_per_THz':ps.tolist(),
                                    'mean_local_PSD_A2_per_THz':local_power.tolist()},
        strain_profile='periodic in-plane sinusoid; input maximum |d epsilon/d in-plane coordinate| in1/nm',
        geometry='r/a=CELL@index, CELL=[[1,-1,1],[1,1,0],[0,0,1]]/2; layer index=k is(001); reference geometry',
        occupation_SHA256=sha256(np.ascontiguousarray(occ)).hexdigest(),
        tensor_semantics='optional Cartesian symmetric per-site delta to effective curvature; no generated physical moire map',
        force_bias_policy=bias_metadata.get('kind','measured_periodic_zero_SCF') if bias_metadata is not None else c.get('bias_policy','centred_model_without_source_bias'),
        force_bias_source=bias_metadata,displacement_origin='declared reference lattice sites, never F0/K recentered',
        demodulation='Joint DC+Omega+2Omega least squares; half-pump field is separated from forced pump-frequency texture and DC, not chosen by dominance over DC',
        occupation_descriptor='periodic bulk-reference vacancy count; open depth only clamps exterior displacements',
        loading_response=dict(composition_reference=composition,fields=loading_metadata,elastic=elastic_metadata),
        material_status=p['status'])
    if c.get('harmonic_control',False):
        phase=np.arccos(1-.5*dt*dt*control_w2)
        expected=.01*np.cos(np.arange(1,steps+1)*phase)
        error=float(np.max(abs(np.array(control)-expected)))
        summary['harmonic_control_max_coordinate_error']=error
        if error>1e-9:
            raise ValueError('Harmonic consumer rejects the source-normalized tensor trajectory')
    return summary,dict(occupied=occ,field=field,static_mean_A=meanq,
                        driven_amplitude_A=amp,pump_frequency_amplitude_A=amp_pump,
                        S_displacement=S,S_occupation=Sn,open_z=np.array(open_z))


def regions(fields,a):
    """Amplitude morphology, NOT a claim of phase-coherent domain size.
    Neighbour polarization agreement cannot segment a finite-q spin texture:
    even a perfectly coherent modulated state can rotate at every neighbour.
    """
    occ=fields['occupied'];shape=occ.shape
    field=fields['field'];norm=np.sqrt((abs(field)**2).sum(axis=0))
    active=occ&(norm>.05);seen=np.zeros_like(active);sizes=[];shapes=[]
    alignment_sum=0.;pairs=0
    if active.any():
        for delta in NN:
            other=np.roll(field,tuple(delta),axis=(1,2,3))
            other_active=np.roll(active,tuple(delta),axis=(0,1,2))
            other_norm=np.roll(norm,tuple(delta),axis=(0,1,2))
            valid=active&other_active
            if fields['open_z'] and delta[2]>0:valid[:,:,:delta[2]]=False
            if fields['open_z'] and delta[2]<0:valid[:,:,delta[2]:]=False
            product=(field.conj()*other).sum(axis=0).real
            alignment_sum+=float((product[valid]/(norm[valid]*other_norm[valid])).sum())
            pairs+=int(valid.sum())
    for flat in np.flatnonzero(active):
        origin=np.unravel_index(flat,shape)
        if seen[origin]:continue
        queue=[origin];coords=[np.zeros(3,int)];seen[origin]=True
        ptr=0;coordinate_map={origin:coords[0]};cycles=set()
        while ptr<len(queue):
            here=queue[ptr];xyz=coords[ptr];ptr+=1
            for delta in NN:
                raw=np.array(here)+delta
                if fields['open_z'] and not 0<=raw[2]<shape[2]:continue
                other=tuple(raw%shape)
                if not active[other]:continue
                unwrapped=xyz+delta
                if seen[other]:
                    if other in coordinate_map:
                        difference=unwrapped-coordinate_map[other]
                        if np.any(difference):cycles.add(tuple(difference//shape))
                    continue
                seen[other]=True;queue.append(other);coords.append(unwrapped)
                coordinate_map[other]=unwrapped
        sizes.append(len(queue))
        cart=np.array(coords,float)@CELL.T*a*.1
        centered=cart-cart.mean(axis=0)
        rank=0
        if cycles:
            winding=np.array(list(cycles))*np.array(shape)
            _,singular,axes=np.linalg.svd(winding@CELL.T,full_matrices=False)
            rank=int(np.linalg.matrix_rank(winding))
            parallel=axes[:rank]
            centered-=centered@parallel.T@parallel
        eigen=np.linalg.eigvalsh(centered.T@centered/len(cart))
        shapes.append(dict(sites=len(queue),periodic_winding_rank=rank,
            morphology=('finite_cluster','periodic_filament','periodic_sheet','percolating_volume')[rank],
            finite_transverse_gyration_nm=np.sqrt(np.maximum(eigen[rank:],0)).tolist(),
            volume_nm3=len(queue)*a**3/4/1000))
    return dict(definition='nearest-FCC connected displacement-amplitude regions above.05A; polarization/phase texture reported separately by S(q), NOT coherent domain radii',
        raw_neighbour_polarization_alignment=alignment_sum/pairs if pairs else None,
        component_sizes=sorted(sizes,reverse=True),largest_shapes=sorted(shapes,key=lambda z:-z['sites'])[:12])


def png(path,rgb):
    rgb=np.repeat(np.repeat(rgb,6,axis=0),6,axis=1).astype(np.uint8)
    h,w,_=rgb.shape
    def chunk(tag,data):
        return struct.pack('>I',len(data))+tag+data+struct.pack('>I',zlib.crc32(tag+data)&0xffffffff)
    raw=b''.join(b'\0'+row.tobytes() for row in rgb)
    path.write_bytes(b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',w,h,8,2,0,0,0))+
                     chunk(b'IDAT',zlib.compress(raw))+chunk(b'IEND',b''))


def save_picture(folder,label,fields):
    nx,ny,nz=fields['occupied'].shape
    occ=fields['occupied'][:,:,nz//2]
    value=fields['field'][0,:,:,nz//2].real
    scaled=np.clip(value/.25,-1,1)
    image=np.full((nx,2*ny+1,3),255.,float)
    image[:,:ny]=occ[:,:,None]*200
    image[:,ny+1:,0]=255*(1-np.maximum(-scaled,0))
    image[:,ny+1:,1]=255*(1-abs(scaled))
    image[:,ny+1:,2]=255*(1-np.maximum(scaled,0))
    image[:,ny+1:]*=occ[:,:,None]
    png(folder/(label+'.png'),image)
    np.savez_compressed(folder/(label+'.npz'),**fields)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('parameters');ap.add_argument('--inspect',action='store_true')
    ap.add_argument('--batch');ap.add_argument('--out',default='results')
    args=ap.parse_args()
    parameters=json.loads(Path(args.parameters).read_text())
    if parameters['loading_law']['kind'] not in ('fixed_x1','anchored_directional_affine_v2'):
        raise ValueError('Unsupported declared loading model')
    if args.inspect:
        answer={}
        for name,p in parameters['materials'].items():
            kernel=shell_kernel(p,(24,24,24)); f=np.fft.fftn(kernel).real
            A=p['K4_100_eV_A4'];B=(3*p['K4_111_eV_A4']-A)/2
            K=p['Kii_eV_A2']+(3*A+2*B)*p['zero_point_axis_variance_A2']
            index=np.unravel_index(np.argmax(f),f.shape)
            answer[name]=dict(kernel_sum_eV_A2=float(kernel.sum()),max_Jq_eV_A2=float(f[index]),
                K_zero_point_eV_A2=K,softest_q=physical_q(index,(24,24,24),p['a_A']),
                full_loading_first_linear_strain=(float(f[index])-K)/p['collective_strain_slope_eV_A2'],
                Gamma_linear_strain=(float(f[0,0,0])-K)/p['collective_strain_slope_eV_A2'])
        print(json.dumps(answer,indent=2));return
    batch=json.loads(Path(args.batch).read_text())
    folder=Path(args.out);folder.mkdir(exist_ok=True)
    summaries=[]
    for case in batch['cases']:
        summary,fields=run_case(parameters,case)
        summaries.append(summary)
        p=parameters['materials'][case.get('material','NiH')]
        summary['regions']=regions(fields,case.get('geometry_a_A',p['a_A']))
        print(json.dumps({k:v for k,v in summary.items() if k not in ('spectrum','regions')}),flush=True)
        if case['label'] in batch.get('representatives',[]):
            save_picture(folder,case['label'],fields)
    (folder/'phase_map.json').write_text(json.dumps(dict(
        model=parameters['model'],loading_law=parameters['loading_law'],strain_map=parameters['strain_map'],
        representation='layer-adapted FCC, periodic in-plane and selectable periodic/open depth; isotropic shell continuation, not unique infinite-lattice force constants',
        PSD_definition='one-sided Hann PSD; collective mean plus mean local probe PSD; resolution is record-window limited',
        results=summaries),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
