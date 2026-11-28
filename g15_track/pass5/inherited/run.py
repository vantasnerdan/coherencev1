"""G15 Pass 4 — one conservative whole-foil multirate field model.
python run.py --prepare /path/to/exact/rung2
python run.py --early
python run.py --all
python run.py --consume     # cold, frozen inputs + saved fields only
"""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,csv,hashlib,json,math,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
import numpy as np
from scipy.linalg import expm
from scipy.special import roots_genlaguerre,roots_legendre,gamma
import sources as src
import fields as fld
import native
from sources import HERE,QE,MD,HB,KB,L,H,RHO,CV,KTH,ALPHA,BULK,ENERGIES,BEAMS,ARMS,NOMINAL
from fields import N,DX,X,Z,XX,YY,ZZ,VOL,CAP,CHANNELS,BANDS


def emit(name,value):
    def cv(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    (HERE/name).write_text(json.dumps(value,indent=2,allow_nan=False,default=cv)+'\n')


def table(name,rows):
    if not rows:raise ValueError('Empty scientific table '+name)
    names=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=names);w.writeheader();w.writerows(rows)


def network(arm):
    A=np.zeros((9,9))
    def send(j,tau,parts):
        A[j,j]=-1/tau
        for k,f in parts:A[k,j]+=f/tau
    send(0,arm['end'],[(3,.3),(4,.25),(6,.4),(7,.05)])
    send(1,arm['recoil'],[(4,.05),(5,.35),(6,.6)])
    send(2,.73*arm['recoil'],[(3,.25),(4,.5),(6,.25)])
    send(3,arm['ne'],[(4,arm['share']),(5,(1-arm['share'])*.3),(6,(1-arm['share'])*.7)])
    send(4,arm['opt'],[(8,.2),(6,.8)])
    for j,t in ((5,2*arm['pd']),(6,arm['pd']),(7,arm['struct']),(8,arm['te'])):send(j,t,[])
    return A


def source_totals(tr,sp,arm):
    result=np.zeros(9);result[0]=(sp+1)*tr['terminal_unresolved_eV']
    dep=np.array(tr['components_eV'])*(sp+1)
    result[1:4]=dep[[1,2,0]]
    key=f"{tr['initial_eV']*(sp+1)/2000:g}_{sp}"
    debit=arm.get('source_debits',{}).get(key,dict(host=0.,electron=0.))
    result[1]-=debit['host'];result[3]-=debit['electron']
    if np.any(result[:4]<0):raise ValueError('Physical coherent source exhausts its same-parent carrier budget')
    return result


def mode_parameters():
    p=src.inputs()['lattice'];K=native.Material(16).reference
    _,_,_,dT=native.thermal_pole(p,K,300.)
    return p,K,float(dT),p['collective_strain_slope_eV_A2']


def spatial_source(beam):
    p=np.zeros_like(VOL);p[:,:,0]=src.beam_cells(np.linspace(-L/2,L/2,N+1),beam);p/=p.sum();return p



def field_path(beam,arm,ec=2.5,control='joined'):
    suffix='' if control=='joined' else '_'+control
    return HERE/f'foil_E{ec:g}_{beam}_{arm}{suffix}.npz'

def accepted_band_parameters():
    with (HERE/'frozen/guided_branches.csv').open() as f:rows=list(csv.DictReader(f))
    out=[]
    for branch in range(3):
        rr=[r for r in rows if int(r['branch'])==branch]
        k=np.array([float(r['k_parallel_m_1']) for r in rr]);f=np.array([float(r['frequency_Hz']) for r in rr])
        order=np.argsort(k);k=k[order];f=f[order]
        out.append((k,f,np.gradient(2*np.pi*f,k)))
    ans=[]
    for band in range(7):
        f=math.sqrt(BANDS[band]*BANDS[band+1])
        for pol,c in enumerate((2440.,3300.,4460.)):
            if band==0 or (band==1 and f<1e7):
                kx,fx,gx=out[pol]
                k=float(np.interp(f,fx,kx));vg=float(np.interp(f,fx,gx));vp=2*np.pi*f/k;law='accepted free-face guided branch'
            else:
                cap=(3.4e12,4.7e12,6.8e12)[pol];k=2*np.pi*cap/c*math.asin(f/cap);vg=c*math.sqrt(1-(f/cap)**2);vp=2*np.pi*f/k;law='accepted bounded Pd sin dispersion'
            if vg<=0:raise ValueError('Nonpositive supported group velocity')
            ans.append(dict(band=band,polarization=pol,frequency_Hz=f,wavevector_m_inverse=k,phase_velocity_m_s=vp,group_velocity_m_s=vg,law=law))
    return ans


def early():
    """Separated receiver and boundary return of an actual accepted band source."""
    pars=accepted_band_parameters();b=next(x for x in pars if x['band']==0 and x['polarization']==0)
    acoustic=np.load(HERE/'frozen/acoustic_source.npz')['energy_J'][0]
    energy=float(np.dot([.85,.15],acoustic[:,0,0]));v=b['group_velocity_m_s'];k=b['wavevector_m_inverse']
    source=-.0025;receiver=.0025;direct=receiver-source;return_path=(L/2-source)+(L/2-receiver)
    width=DX/v;end=(return_path+4*DX)/v;t=np.linspace(0,end,1200)
    def pulse(path):
        age=t-path/v;u=age/width
        return np.where((u>=0)&(u<=2),np.sin(np.pi*u/2)**2,0.)
    # The declared narrowband FV strain-energy metric uses restoring phase
    # stiffness rho*vp²; group velocity controls transport, not this modulus.
    vp=b['phase_velocity_m_s'];eps=math.sqrt(2*energy/(RHO*vp*vp*VOL[N//2,N//2,0]))
    rate=1/50e-6+v/.2;direct_eps=eps*math.exp(-rate*direct/(2*v))*pulse(direct)
    returned=-math.sqrt(.8)*eps*math.exp(-rate*return_path/(2*v))*pulse(return_path)
    p,K,dT,dE=mode_parameters();domega=dE/(6*K)*(direct_eps+returned)
    # A subsequent finite-budget momentum insertion reads the arriving field.
    velocity=vp*direct_eps;impulse=math.sqrt(2*MD*5000*QE);work=impulse*velocity/QE
    rows=[dict(time_s=float(tt),direct_signed_strain=float(a),reflected_signed_strain=float(r),absorbing_boundary_strain=float(a),
        coupled_fractional_optical_frequency=float(o),strain_coupling_off_fractional_frequency=0.,temperature_change_K=0.,
        later_impulse_p_dot_v_eV=float(w),later_impulse_frozen_medium_cross_work_eV=0.) for tt,a,r,o,w in zip(t,direct_eps,returned,domega,work)]
    table('first_travel.csv',rows)
    result=dict(decision='A separated finite-volume receiver has a nonzero causal arrival and a signed boundary return on the same accepted source. This establishes propagation; its optical modulation is extremely small, not a hot-spot interference claim.',
        source_x_m=source,receiver_x_m=receiver,direct_path_m=direct,return_path_m=return_path,
        source_energy_J=energy,group_velocity_m_s=v,phase_velocity_m_s=b['phase_velocity_m_s'],
        direct_first_arrival_s=direct/v,return_first_arrival_s=return_path/v,
        maximum_direct_strain=float(direct_eps.max()),maximum_return_abs_strain=float(abs(returned).max()),
        maximum_fractional_optical_frequency=float(abs(domega).max()),maximum_later_p_dot_v_eV=float(work.max()),
        prefront_maximum=float(abs(direct_eps[t<direct/v]).max(initial=0.)),
        finite_resolution='17-cell characteristic observable, compact causal sin-squared envelope; carrier phase inside cell not resolved',
        control='same parent energy, propagation and damping; reflection removed only in absorbing comparator; strain coupling off removes optical response without removing travel')
    emit('first_decision.json',result);print(json.dumps(result),flush=True)


def carrier_kernel(arm):
    import inspect
    path=HERE/f'paths_{arm["name"]}.npz'
    epoch=hashlib.sha256((inspect.getsource(fld.path_green)+inspect.getsource(network)+src.digest(HERE/'native.py')).encode()).hexdigest()
    if path.exists():
        saved=np.load(path);valid='epoch' in saved and str(saved['epoch'])==epoch
    else:valid=False
    if valid:segments=saved['segments'];ends=saved['heat_endpoints']
    else:
        segments,ends=fld.path_green(arm,HERE)
        np.savez_compressed(path,segments=segments,heat_endpoints=ends,epoch=np.array(epoch),
            segment_columns=np.array(['origin','sample','channel','start_s','end_s','x0','y0','z0','x1','y1','z1','weight','omega','phase','kx','ky','kz','speed','elastic_scatter_rate']))
    # Normalize the spatial/time quadrature to the exact inherited carrier
    # residence. The stochastic paths specify geometry, not a new branching fit.
    A=network(arm);res=-np.linalg.inv(A)
    kernels={};metrics=[]
    for ch in range(9):
        for origin in range(4):
            j=(0,1,2,3)[origin];rows=segments[(segments[:,0]==origin)&(segments[:,2]==ch)]
            if not len(rows):continue
            w=(rows[:,4]-rows[:,3])*rows[:,11];raw=float(w.sum());target=float(res[ch,j])
            if target<=0:continue
            weights=w/raw;offset=(rows[:,5:8]+rows[:,8:11])/2;age=(rows[:,3]+rows[:,4])/2
            kernels[ch,origin]=(offset,age,weights,target)
            metrics.append(dict(arm=arm['name'],field=CHANNELS[ch],origin=origin,exact_residence_s=target,raw_path_residence_s=raw,
                geometry_quadrature_relative_residence_difference=raw/target-1,
                mean_age_s=float(age@weights),rms_travel_m=float(np.sqrt(np.sum(weights*np.sum(offset*offset,axis=1)))),
                maximum_sampled_travel_m=float(np.linalg.norm(offset,axis=1).max()),
                lateral_cell_m=DX,first_depth_cell_m=fld.DZ[0],parent_phase='one shared source mark; optical k/omega retained; incoherent fields have no absolute phase'))
    return kernels,metrics









def thermo_capacity(strain,thermal_coefficient):
    return CAP*(1-3*BULK*ALPHA/CV*strain.reshape(VOL.shape))+thermal_coefficient


def source_tiles(modes):
    xy=modes['source_xy']-L/2;axis=np.unique(xy[:,0]);cell=((axis+L/2)/DX).astype(int)
    edges=np.r_[-.005,(axis[:-1]+axis[1:])/2,.005]
    for j in range(len(axis)-1):
        if cell[j]!=cell[j+1]:edges[j+1]=-L/2+(cell[j]+1)*DX
    ix=np.searchsorted(axis,xy[:,0]);iy=np.searchsorted(axis,xy[:,1])
    bounds=np.column_stack((edges[ix],edges[ix+1],edges[iy],edges[iy+1]))
    return edges,ix,iy,bounds


def source_probabilities(modes,beam):
    edges,ix,iy,_=source_tiles(modes)
    return src.beam_cells(edges,beam)[ix,iy]


def damping_distribution(modes):
    displacement=modes['displacement'].reshape(-1,3,len(modes['omega']))
    density=RHO*VOL.ravel()[:,None]*np.sum(displacement*displacement,axis=1)
    return density/density.sum(axis=0)



def macro(arm,beam,calendar,kernels,tracks,control='joined',ec=2.5):
    import carriers,acoustic
    modes=np.load(HERE/'frozen/modes.npz');omega=modes['omega'];S=modes['strain'][:,:3].sum(axis=1)
    pressure_basis=modes['pressure'];p,K,dT,g=mode_parameters();coefficient=native.constitutive(arm,ec,control=control)
    gamma=1/50e-6+modes['speed']/.2 if arm['loss']=='retained' else 4460*np.log(10)/10*32*omega/(2*np.pi*1e6)
    strain_on=control not in ('strain_off','frozen');temperature_on=control not in ('temperature_off','frozen')
    source=spatial_source(beam);rate=src.current(ec)/QE
    totals=sum(fr*source_totals(tracks[f'{ec:g}_{sp}'],sp,arm) for sp,fr in enumerate((.85,.15)))
    carrier=carriers.Operator(arm,ec,kernels,tracks,coefficient['species_residence_ratio'])
    base=carrier.field(source,rate)
    power=src.current(ec)*2000*ec
    radiator=acoustic.cached_operator(arm,ec,beam,tracks,control)
    ac,_,_=radiator.snapshot(calendar,1.)
    bandrows=radiator.rows
    acoustic_power=rate*radiator.total_per_parent
    shortshape=carrier.heat_rate(source,rate);shortshape/=shortshape.sum()
    lossshape,_=radiator.heat_shape();lossshape/=lossshape.sum()
    heatshape=((power-acoustic_power)*shortshape+acoustic_power*lossshape)/power
    bathbase=radiator.bath_snapshot(None,1.)
    stress_ratio=coefficient['optical_strain_ratio'];thermal_ratio=coefficient['thermal_energy_ratio_per_K'] if temperature_on else 0.
    electron_coupling=3*stress_ratio/g*coefficient['dk_per_energy_density']/VOL
    sw=source_probabilities(modes,beam);pressure=pressure_basis@sw
    source_points=np.column_stack((modes['source_xy']-L/2,np.zeros(len(sw))))
    impulse=math.sqrt(2*MD*2000*ec*QE)*(.85+.15*math.sqrt(2))
    localization=damping_distribution(modes)
    heat=fld.Heat();temp=np.zeros_like(VOL);q=np.zeros(len(omega));v=q.copy();oldfast=np.zeros_like(base)
    history=[];snap={};qhistory=[];cumulative_input=0.;exits=np.zeros(4);maxdef=0.;pressure_work=0.;momentum_work=0.;oldac=0.
    dt=.001;stored_ac=float(ac.sum());source_energy=2000*ec*QE;last_temperature=np.zeros_like(VOL);last_q=q.copy()
    def interaction(fast,strain):
        z=stress_ratio*fast[4] if strain_on else np.zeros_like(VOL)
        ee=fast[3]+fast[8]
        cross=electron_coupling*fast[4]*ee
        return z,thermal_ratio*fast[4],cross
    for j in range(2000):
        a=j*dt;b=a+dt;mult=float(calendar.integral(a,b,beam)/dt);input_J=rate*source_energy*mult*dt
        fast=(base+bathbase)*mult;ac_J=radiator.inventory(calendar,b);strain=S@q
        z,thermal,cross=interaction(fast,strain);oz,ot,oc=interaction(oldfast,strain)
        delta_fast=float((fast.sum(axis=0)-oldfast.sum(axis=0)+(z-oz)*strain.reshape(VOL.shape)+cross-oc+(thermal-ot)*temp).sum())
        ac_exit=radiator.exits(calendar,a,b)
        bath_delta=float(bathbase.sum()*mult-(bathbase.sum()*float(calendar.integral(a-dt,a,beam)/dt) if j else 0.))
        short_amount=(power-acoustic_power)*mult*dt-delta_fast+bath_delta
        acoustic_amount=acoustic_power*mult*dt-(ac_J-oldac)-ac_exit-bath_delta
        available_map=short_amount*shortshape+acoustic_amount*lossshape
        cap=thermo_capacity(strain,thermal);u=cap*temp
        u,temp,ex,defect=heat.step(u,available_map,dt,cap)
        pv=coefficient['pressure_ratio']*(fast[3]+fast[8])
        force_int=fld.interaction_force(S,CAP*temp,pv,z)
        force_p=rate*mult*impulse*pressure
        oldq=q.copy();oldv=v.copy();before=.5*float(np.sum(v*v+(omega*q)**2))
        q,v,work=fld.oscillator(q,v,force_int+force_p,omega,gamma,dt)
        after=.5*float(np.sum(v*v+(omega*q)**2));damping=before+work-after
        dq=q-oldq;de=S@dq;electronic_work=pv*de.reshape(VOL.shape);pressure_work+=float(electronic_work.sum())
        # Exact opposite physical work, including the pair mean-field energy.
        electron=fast[3]+fast[8];debit=electronic_work/(1+electron_coupling*fast[4])
        fraction=np.divide(debit,electron,out=np.zeros_like(debit),where=electron>0)
        fast[3]*=1-fraction;fast[8]*=1-fraction
        if min(float(fast[3].min()),float(fast[8].min()))<0:raise ValueError('Electronic pressure work exhausts its own carrier')
        point_work=rate*mult*impulse*sw*(pressure_basis.T@dq);wp=float(point_work.sum());momentum_work+=wp
        mode_weight=gamma*(oldv*oldv+v*v+(omega*oldq)**2+(omega*q)**2)
        spatial=localization@mode_weight
        spatial=spatial/spatial.sum() if spatial.sum()>0 else heatshape.ravel()
        strain=S@q;z,thermal,cross=interaction(fast,strain);newcap=thermo_capacity(strain,thermal)
        correction=damping*spatial.reshape(VOL.shape)-fld.deposit(source_points,point_work)
        temp+=correction/newcap;u=newcap*temp
        exits[:3]+=ex;exits[3]+=ac_exit;cumulative_input+=input_J
        fast_nonT=fast.sum(axis=0)+z*strain.reshape(VOL.shape)+cross
        balance=float(u.sum()+fast_nonT.sum()+ac_J+after+exits.sum()-cumulative_input);maxdef=max(maxdef,abs(balance))
        if j%10==9 or j==0:
            history.append(dict(arm=arm['name'],closure=arm['closure'],allocation=arm['allocation'],beam=beam,E_cm_keV=ec,control=control,
                time_s=b,temporal_aperture_s=dt,input_J=cumulative_input,heat_enthalpy_J=float(u.sum()),
                **{name+'_J':float(fast[k].sum()) for k,name in enumerate(CHANNELS)},acoustic_seven_J=ac_J,low_elastic_kinetic_J=after,
                interaction_J=float(u.sum()-np.sum(CAP*temp)+fast_nonT.sum()-fast.sum()),exit_mount_J=exits[0],exit_faces_J=exits[1],
                exit_radiation_J=exits[2],exit_acoustic_boundary_J=exits[3],external_momentum_work_J=momentum_work,
                native_electron_pressure_work_J=pressure_work,balance_J=balance,minimum_capacity_J_K=float(newcap.min()),max_temperature_K=300+float(temp.max())))
            if j%100==99 or j==0:qhistory.append(np.r_[b,q,v])
        if j in (9,499,1999):
            label=f't{b:g}';eps=strain.reshape(VOL.shape);snap[label+'_deltaT']=temp.copy();snap[label+'_strain']=eps.copy()
            snap[label+'_displacement']=(modes['displacement']@q).reshape(N,N,3,3)
            snap_fast=carrier.field(source,rate,calendar,b,beam)+radiator.bath_snapshot(calendar,b)
            snap_ac,phase_green,normalizer=radiator.snapshot(calendar,b)
            snap[label+'_fast_J']=snap_fast;snap[label+'_acoustic_J']=snap_ac
            snap[label+'_acoustic_spatial_quadrature_normalizer']=np.array(normalizer)
            if ec==2.5 and j==1999:
                suffix='' if control=='joined' else '_'+control
                np.savez_compressed(HERE/f'acoustic_phase_E{ec:g}_{beam}_{arm["name"]}{suffix}.npz',conditional_first_moment_green_s=phase_green,
                    phase_scope=np.array('normal-incidence signs retained; unresolved oblique-boundary phase stochastic; not energy from squared first moment'))
            local_delta=coefficient['dk_per_energy_density']*(snap_fast[3]+snap_fast[8])/VOL
            slow=g*((eps/3 if strain_on else 0.)+(ALPHA*temp if temperature_on else 0.))
            bare=K+slow+local_delta
            evolved,Qthermal,_,_=native.thermal_pole(p,bare,300+temp if temperature_on else np.full_like(temp,300.))
            onsite=p['Kii_eV_A2']+(2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4'])*p['zero_point_axis_variance_A2']+slow
            if arm['allocation']=='onsite':onsite+=local_delta
            snap[label+'_bare_K']=bare;snap[label+'_K']=evolved;snap[label+'_thermal_covariance_A2']=Qthermal
            snap[label+'_omega']=np.sqrt(evolved*QE/1e-20/(p['mass_eVc2']*QE/299792458.**2))
            snap[label+'_bare_omega']=np.sqrt(bare*QE/1e-20/(p['mass_eVc2']*QE/299792458.**2))
            snap[label+'_onsite_K']=onsite+(evolved-bare)
            snap[label+'_J_over_K']=(onsite-bare)/(onsite+evolved-bare)
            snap[label+'_retained_shell_J_over_onsite']=sum(r['J_sum_eV_A2'] for r in p['shells'])/(onsite+evolved-bare)
            snap[label+'_omega_over_kappa']=snap[label+'_omega']*arm['opt'];snap[label+'_kappa_s_inverse']=np.full_like(VOL,1/arm['opt'])
            snap[label+'_supported_remote_h2_over_threshold']=np.zeros_like(VOL)
        if j==1998:last_temperature=temp.copy();last_q=q.copy()
        oldfast=fast;oldac=ac_J
    snap['last_temperature_derivative_K_s']=(temp-last_temperature)/dt;snap['last_modal_derivative']=np.r_[(q-last_q)/dt,v]
    snap['warm_state_deltaT_K']=temp.ravel().copy()
    snap={key:(value.astype(np.float32) if isinstance(value,np.ndarray) and value.dtype.kind=='f' and value.ndim>=2 else value) for key,value in snap.items()}
    np.savez_compressed(field_path(beam,arm['name'],ec,control),**snap,source_fraction=source,volume_m3=VOL,x_m=X,z_m=Z,
        macro_mode_history=np.array(qhistory),macro_mode_history_columns=np.array(['time_s']+[f'q_{i}' for i in range(len(q))]+[f'p_{i}' for i in range(len(q))]))
    return history,bandrows,dict(arm=arm['name'],beam=beam,E_cm_keV=ec,control=control,max_balance_J=maxdef,relative_balance=maxdef/cumulative_input),snap,dict(base=base,bath=bathbase,carrier=carrier,radiator=radiator,coefficient=coefficient,heat_power=power*heatshape)



def fine_window(arm,beam,calendar,snap,context):
    """Same-foil count innovations and conservative delayed source interfaces."""
    modes=np.load(HERE/'frozen/modes.npz');w=modes['omega'];S=modes['strain'][:,:3].sum(axis=1)
    pressure=modes['pressure'];U=modes['displacement'];p,K,dT,g=mode_parameters()
    coefficient=context['coefficient'];base=context['base']+context['bath'];carrier=context['carrier'];radiator=context['radiator']
    gamma=1/50e-6+modes['speed']/.2 if arm['loss']=='retained' else 4460*np.log(10)/10*32*w/(2*np.pi*1e6)
    prob=source_probabilities(modes,beam);xy=modes['source_xy']-L/2
    ids=((xy[:,0]+L/2)/DX).astype(int)*N+((xy[:,1]+L/2)/DX).astype(int)
    lookup=np.full(N*N,-1);lookup[carrier.cells]=np.arange(carrier.ns);ids=lookup[ids]
    dt=50e-9;burn=4000;samples=2400;steps=burn+samples;prefix=1024;start=2-steps*dt
    if carrier.delay_s>prefix*dt:raise ValueError('Short-source history exceeds shared fine prefix')
    times=start+(np.arange(steps+prefix)-prefix)*dt;mult=calendar.integral(times,times+dt,beam)/dt;rate=src.current(2.5)/QE
    rng=np.random.default_rng(5151);count=rng.poisson(rate*dt*mult[:,None]*prob).astype(np.int32)
    dp=rng.binomial(count,.85).astype(np.int32);molecules=count-dp
    cells=np.zeros((2,len(times),carrier.ns))
    for node,cell in enumerate(ids):cells[0,:,cell]+=dp[:,node];cells[1,:,cell]+=molecules[:,node]
    deposited=carrier.delayed_heat(cells,dt)[prefix:];del cells
    aprefix=int(np.ceil(.004/dt));atimes=start+(np.arange(aprefix+steps)-aprefix)*dt
    earlier=np.random.default_rng(5150);total=earlier.poisson(rate*calendar.integral(atimes,atimes+dt,beam))
    atomic=earlier.binomial(total,.85);rain=np.stack((atomic,total-atomic)).astype(float)
    rain[:,-len(times):]=np.stack((dp.sum(axis=1),molecules.sum(axis=1)))
    if radiator.maximum_delay>aprefix*dt:raise ValueError('Acoustic history exceeds shared physical prefix')
    deposited+=radiator.delayed_heat(rain,dt,start=aprefix);del rain
    source_hash=hashlib.sha256(count.tobytes()+dp.tobytes()).hexdigest();controls=('joined','temperature_off','strain_off','frozen')
    impulse=math.sqrt(2*MD*5000*QE);mean_p=impulse*(.85+.15*np.sqrt(2));mean_force=rate*mean_p*(pressure@prob)
    warm=snap['warm_state_deltaT_K'].reshape(VOL.shape) if 'warm_state_deltaT_K' in snap else np.asarray(snap['t2_deltaT'],float)
    temp=warm-steps*dt*snap['last_temperature_derivative_K_s']
    stress=coefficient['optical_strain_ratio']*base[4];pv=coefficient['pressure_ratio']*(base[3]+base[8])
    fint=fld.interaction_force(S,CAP*temp,pv,stress)
    static_q=np.broadcast_to((fint+mean_force)/(w*w),(4,len(w))).copy()
    static_q[2:]=(fld.interaction_force(S,CAP*temp,pv,np.zeros_like(stress))+mean_force)/(w*w)
    import prehistory
    initial_q,initial_v=prehistory.state(calendar,start,beam,modes,temp,context['heat_power'],mean_force+S.T@(pv-stress).ravel(),gamma)
    off_q,off_v=prehistory.state(calendar,start,beam,modes,temp,context['heat_power'],mean_force+S.T@pv.ravel(),gamma)
    q=np.array([initial_q,initial_q,off_q,off_q]);v=np.array([initial_v,initial_v,off_v,off_v])
    heat=fld.Heat();localization=damping_distribution(modes)
    rows=[];measure={c:[] for c in controls};maps=[];temperature=[];maptimes=[];state_times=[];states=[];state_counts=[];state_temperature=[]
    selected=set((burn+np.sort(np.random.default_rng(741).choice(samples,64,replace=False))).tolist())
    maximum_defect=0.;minimum_budget=np.inf
    for j in range(steps):
        k=j+prefix;t=start+j*dt
        if j in selected:
            state_times.append(t);states.append(np.concatenate((q,v),axis=1))
            state_counts.append(np.stack((dp[k],molecules[k])));state_temperature.append(temp[:,:,0].copy())
        scale=float(count[k].sum()/(rate*dt));fast=base*scale
        z=coefficient['optical_strain_ratio']*fast[4];bt=coefficient['thermal_energy_ratio_per_K']*fast[4]
        cap=thermo_capacity(S@q[0],bt);_,temp,_,defect=heat.step(cap*temp,deposited[j],dt,cap)
        maximum_defect=max(maximum_defect,abs(defect))
        pv=coefficient['pressure_ratio']*(fast[3]+fast[8]);force0=fld.interaction_force(S,CAP*temp,pv,z)
        force=np.broadcast_to(force0,q.shape).copy();force[2:]=fld.interaction_force(S,CAP*temp,pv,np.zeros_like(z))
        mean=rate*mult[k]*dt*prob;kick=pressure@(impulse*(dp[k]+np.sqrt(2)*molecules[k]-mean*(.85+.15*np.sqrt(2))))
        kickwork=v@kick+.5*float(kick@kick);v+=kick
        minimum_budget=min(minimum_budget,float(count[k].sum()*5000*QE-kickwork.max()))
        if minimum_budget<0:raise ValueError('Parent momentum insertion exhausts beam budget')
        applied=rate*mult[k]*mean_p*(pressure@prob);oldq=q.copy();oldv=v.copy()
        before=.5*float(np.sum(v[0]**2+(w*q[0])**2))
        for ci in range(4):q[ci],v[ci],_=fld.oscillator(q[ci],v[ci],force[ci]+applied,w,gamma,dt)
        after=.5*float(np.sum(v[0]**2+(w*q[0])**2));dq=q[0]-oldq[0]
        damping=before+float((force[0]+applied)@dq)-after;source_work=float(applied@dq)+kickwork[0]
        pressure_work=pv*(S@dq).reshape(VOL.shape);electron=fast[3]+fast[8]
        fraction=np.divide(pressure_work,electron,out=np.zeros_like(electron),where=electron>0)
        fast[3]*=1-fraction;fast[8]*=1-fraction
        if min(fast[3].min(),fast[8].min())<0:raise ValueError('Reciprocal fine work depletes electrons')
        distribution=localization@(gamma*(oldv[0]**2+v[0]**2+(w*oldq[0])**2+(w*q[0])**2));distribution/=distribution.sum()
        capnew=thermo_capacity(S@q[0],bt)
        temp+=(damping*distribution.reshape(VOL.shape)-source_work*spatial_source(beam)-pressure_work)/capnew
        if j>=burn and (j-burn)%8==0:
            dynamic=(S@(q[0]-static_q[0])).reshape(VOL.shape)
            maps.append(dynamic[:,:,0].astype(np.float32));temperature.append(temp[:,:,0].astype(np.float32));maptimes.append(t+dt)
            velocity=(U@v[0]).reshape(N,N,3,3)
            rows.append(dict(arm=arm['name'],beam=beam,time_s=t+dt,window_s=samples*dt,dt_s=dt,parent_count=int(count[k].sum()),
                signed_strain_peak=float(dynamic.max()),signed_strain_hole=float(dynamic.min()),maximum_velocity_m_s=float(abs(velocity).max()),
                elastic_kinetic_J=after,maximum_transport_temperature_K=300+float(temp.max()),minimum_remaining_parent_budget_J=minimum_budget,
                maximum_heat_step_defect_J=maximum_defect,source_parent_hash=source_hash,native_pressure_work_J=float(pressure_work.sum()),
                source_spatial_moment='short carriers per cell; acoustic heat conditional spatial mean, same total parents'))
            site=np.unravel_index(np.argmax(abs(dynamic)),VOL.shape)
            for ci,control in enumerate(controls):
                ee=0. if ci>=2 else float((S@q[ci]).reshape(VOL.shape)[site]);tt=0. if ci in (1,3) else float(temp[site])
                pole,_,_,_=native.thermal_pole(p,K+g*(ee/3+ALPHA*tt),300+tt);quiet,_,_,_=native.thermal_pole(p,K,300.)
                vel=float((U@v[ci]).reshape(N,N,3,3)[site+(2,)])
                measure[control].append((np.sqrt(pole/quiet),impulse*vel/QE,float(np.linalg.norm(q[ci]-q[0]))))
    np.savez_compressed(HERE/f'phase_{beam}_{arm["name"]}.npz',time_s=state_times,modal_phase_space=np.array(states),
        controls=np.array(controls),source_counts=np.array(state_counts),source_xy_m=xy,source_tile_bounds_m=source_tiles(modes)[3],
        surface_temperature_K=300+np.array(state_temperature),map_time_s=maptimes,signed_surface_strain=np.array(maps),
        surface_temperature_history_K=300+np.array(temperature),omega=w,static_reference_q=static_q[0],
        physical_source_window_start_s=start,exposure_s=samples*dt,phase_sampling_exposure_s=len(selected)*dt,dt_s=dt,source_seed=5151,source_parent_hash=np.array(source_hash))
    measured=[]
    for control,values in measure.items():
        a=np.array(values);measured.append(dict(arm=arm['name'],beam=beam,control=control,source_parent_hash=source_hash,
            min_frequency_ratio=float(a[:,0].min()),max_frequency_ratio=float(a[:,0].max()),max_abs_later_impulse_cross_work_eV=float(abs(a[:,1]).max()),
            maximum_reentrant_modal_displacement_difference=float(a[:,2].max()),exposure_s=samples*dt,
            interpretation='conditional fine modes on common transported T; macro files contain fully evolved medium-off controls'))
    return rows,measured


def optical_histogram(tr,sp,arm,sample_times=None):
    """Exact row/constituent source boundaries on the physical fast calendar."""
    A=network(arm);Ai=np.linalg.inv(A);z=np.r_[0.,[r['depth_A']*1e-10 for r in tr['timeline']]]
    t=np.r_[0.,[r['time_fs']*1e-15 for r in tr['timeline']]]
    dep=np.diff(np.array([[0.,0.,0.]]+[[r[k] for k in ('electron_eV','host_recoil_eV','D_recoil_eV')] for r in tr['timeline']]),axis=0)
    # The tracked sound budget is removed once, from its host source.
    totals=source_totals(tr,sp,arm)/(sp+1)
    dep[:,1]*=totals[1]/dep[:,1].sum();dep[:,0]*=totals[3]/dep[:,0].sum()
    delay=.74e-10/math.sqrt(2*tr['initial_eV']*QE/MD) if sp else 0.;offsets=(0.,delay) if sp else (0.,)
    end=t[-1]+delay+20*max(arm['opt'],arm['ne'],arm['end'])
    edges=np.unique(np.r_[np.concatenate([t+off for off in offsets]),end,[] if sample_times is None else sample_times])
    y=np.zeros(9);hist=[np.r_[y,0.]];times=[0.];integral=np.zeros(9);input_eV=0.;sample_count=0
    for left,right in zip(edges[:-1],edges[1:]):
        mid=(left+right)/2;force=np.zeros(9)
        for off in offsets:
            k=np.searchsorted(t,mid-off,side='right')-1
            if 0<=k<len(dep):force[[3,1,2]]+=dep[k]/(t[k+1]-t[k])
        count=max(1,int(math.ceil((right-left)/.5e-15)));dt=(right-left)/count
        E=expm(A*dt);F=Ai@(E-np.eye(9))
        for j in range(count):
            new=E@y+F@force;integral+=Ai@(new-y-force*dt);y=new;input_eV+=float(force.sum()*dt);sample_count+=1
            if sample_count%8==0 or j==count-1:
                hist.append(np.r_[y,y[4]+integral[4]/arm['opt']]);times.append(left+(j+1)*dt)
        for off in offsets:
            if abs(right-t[-1]-off)<1e-27:y[0]+=tr['terminal_unresolved_eV'];input_eV+=tr['terminal_unresolved_eV']
    h=np.array(hist);tt=np.array(times);target=np.linalg.solve(-A,source_totals(tr,sp,arm))[4]
    expected=float(source_totals(tr,sp,arm).sum())
    return tt,h,z[-1],dict(source_discretization_relative_error=input_eV/expected-1,
        optical_residence_relative_error=float(integral[4]/target-1),sampled_final_untransferred_eV=float(y.sum()),D2_delay_s=delay)


def sigma(E,ch='p',U=.04):
    E=np.asarray(E,float);out=np.zeros_like(E);ok=E>0;e=E[ok]
    coeff=([5.5576e4,2.1054e2,-3.2638e-2,1.4987e-6,1.8181e-10] if ch=='p' else [5.3701e4,3.3027e2,-1.2706e-1,2.9327e-5,-2.5151e-9])
    S=sum(c*e**i for i,c in enumerate(coeff));out[ok]=S/np.sqrt(e*(e+U))*np.exp(-31.397/np.sqrt(e+U))*1e-31;return out
GX,GW=roots_genlaguerre(8,.5);MU,MW=roots_legendre(8)


def moving_rate(lab,theta,ch,U):
    lab,theta=np.broadcast_arrays(np.asarray(lab),np.asarray(theta));e=lab[...,None,None]/2
    thermal=GX[:,None]*theta[...,None,None]/2000;er=e+thermal-2*np.sqrt(e*thermal)*MU[None,:]
    return np.sum(sigma(er,ch,U)*np.sqrt(er/e)*GW[:,None]/gamma(1.5)*MW[None,:]/2,axis=(-2,-1))


def pump_rate(lab,energy,ch,theta,U):
    lab,energy=np.broadcast_arrays(np.asarray(lab),np.asarray(energy));phi=(MU+1)*np.pi/4
    kin=energy[...,None,None]*np.sin(phi[:,None])**2/1000
    effective=lab[...,None,None]+kin-2*np.sqrt(lab[...,None,None]*kin)*MU[None,:]
    return np.sum(moving_rate(effective,theta,ch,U)*np.sqrt(effective/lab[...,None,None])*MW[:,None]*MW[None,:]/4,axis=(-2,-1))


def pump_delta(lab,energy,ch,theta,U):
    """Same physical quadratures with a centered smooth Maxwell-kernel lookup.

    This retains each oscillator phase/direction, not its mean kinetic energy.
    A connected relative kernel avoids subtracting two near-equal tiny rates.
    """
    from scipy.interpolate import CubicSpline
    lab=np.asarray(lab);energy=np.asarray(energy);quiet=moving_rate(lab,theta,ch,U)
    maximum=float(energy.max(initial=0.))/1000
    if maximum==0:return np.zeros((len(energy),len(lab)))
    width=maximum/lab+2*np.sqrt(maximum/lab)
    nodes=np.linspace(-1.,1.,129)
    sampled=moving_rate(lab[:,None]*(1+width[:,None]*nodes),theta,ch,U)/quiet[:,None]-1
    spline=CubicSpline(nodes,sampled,axis=1)
    phi=(MU+1)*np.pi/4
    kinetic=energy[:,None,None,None]*np.sin(phi[None,None,:,None])**2/1000
    delta=kinetic/lab[None,:,None,None]-2*np.sqrt(kinetic/lab[None,:,None,None])*MU[None,None,None,:]
    x=delta/width[None,:,None,None]
    index=np.minimum(np.maximum(np.searchsorted(nodes,x,side='right')-1,0),len(nodes)-2)
    coordinate=x-nodes[index];coeff=spline.c[:,index,np.arange(len(lab))[None,:,None,None]]
    ratio=((coeff[0]*coordinate+coeff[1])*coordinate+coeff[2])*coordinate+coeff[3]
    root=np.sqrt(1+delta)
    connected=ratio*root+delta/(root+1)
    return quiet[None,:]*np.sum(connected*MW[None,None,:,None]*MW[None,None,None,:]/4,axis=(-2,-1))


def consume(write=True):
    """Last-stage actual-covariance consumer; no producer execution."""
    import yielding,sys
    return yielding.consume(sys.modules[__name__],write)


def figures():
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    extent=[-L/2*1e3,L/2*1e3,-L/2*1e3,L/2*1e3]
    for beam in BEAMS:
        arm=NOMINAL;d=np.load(field_path(beam,arm['name']))
        fig,ax=plt.subplots(3,4,figsize=(14,10))
        for i,t in enumerate((.01,.5,2.)):
            for j,(suffix,label) in enumerate((('deltaT','temperature excess (K)'),('strain','signed div u'),('K','conditional cell-mean pole (eV/Å²)'),('J_over_K','J sum / K'))):
                a=d[f't{t:g}_'+suffix][:,:,0];im=ax[i,j].imshow(a.T,origin='lower',extent=extent,aspect='equal');fig.colorbar(im,ax=ax[i,j],shrink=.7)
                ax[i,j].set_title(f'{t:g} s: {label}',fontsize=9);ax[i,j].set_xlabel('x (mm)');ax[i,j].set_ylabel('y (mm)')
        fig.suptitle(beam+': whole foil, surface 0–250 nm FV average; no atomic peak inferred');fig.tight_layout();fig.savefig(HERE/f'maps_{beam}.png',dpi=145);plt.close(fig)
        fig,ax=plt.subplots(3,4,figsize=(14,10));a=d['t2_fast_J']
        for j,name in enumerate(CHANNELS):
            im=ax.flat[j].imshow((a[j,:,:,0]/VOL[:,:,0]).T,origin='lower',extent=extent);fig.colorbar(im,ax=ax.flat[j],shrink=.7);ax.flat[j].set_title(name+' (J/m³)',fontsize=9)
        for spare in ax.flat[len(CHANNELS):]:spare.set_visible(False)
        fig.suptitle(beam+': every short carrier, same 2 s state; sparse phases not implied by these means');fig.tight_layout();fig.savefig(HERE/f'carriers_{beam}.png',dpi=145);plt.close(fig)
        fig,ax=plt.subplots(2,4,figsize=(14,7));a=d['t2_acoustic_J'].sum(axis=1)
        for band in range(7):
            im=ax.flat[band].imshow((a[band,:,:,0]/VOL[:,:,0]).T,origin='lower',extent=extent);fig.colorbar(im,ax=ax.flat[band],shrink=.7)
            ax.flat[band].set_title(f'{BANDS[band]:.2g}–{BANDS[band+1]:.2g} Hz, J/m³',fontsize=8)
        phase=np.load(HERE/f'phase_{beam}_{arm["name"]}.npz');a=phase['signed_surface_strain'];im=ax.flat[7].imshow(a[np.argmax(np.max(abs(a),axis=(1,2)))].T,origin='lower',extent=extent,cmap='RdBu_r');fig.colorbar(im,ax=ax.flat[7],shrink=.7);ax.flat[7].set_title('resolved signed low-mode wave strain')
        fig.suptitle(beam+': retarded seven-band energy envelopes and resolved standing modes');fig.tight_layout();fig.savefig(HERE/f'acoustic_{beam}.png',dpi=145);plt.close(fig)
    with (HERE/'energy_history.csv').open() as f:rr=list(csv.DictReader(f))
    fig,ax=plt.subplots(1,2,figsize=(12,4))
    for beam in BEAMS:
        r=[x for x in rr if x['arm']==NOMINAL['name'] and x['beam']==beam and float(x['E_cm_keV'])==2.5 and x['control']=='joined'];t=[float(x['time_s']) for x in r]
        ax[0].plot(t,[float(x['heat_enthalpy_J']) for x in r],label=beam+' heat')
        ax[0].plot(t,[float(x['input_J']) for x in r],ls=':',label=beam+' input')
        for key in ('exit_mount_J','exit_faces_J','exit_radiation_J','low_elastic_kinetic_J','acoustic_seven_J','D_optical_J'):
            ax[1].plot(t,[max(float(x[key]),1e-30) for x in r],label=beam+' '+key)
    ax[0].set(xlabel='whole-foil calendar (s)',ylabel='J');ax[1].set(xlabel='whole-foil calendar (s)',ylabel='J',yscale='log');ax[0].legend(fontsize=6);ax[1].legend(fontsize=5,ncol=2);fig.tight_layout();fig.savefig(HERE/'energy_history.png',dpi=160);plt.close(fig)
    with (HERE/'first_travel.csv').open() as f:r=list(csv.DictReader(f))
    fig,ax=plt.subplots(figsize=(9,4));t=np.array([float(x['time_s']) for x in r])*1e6
    for key in ('direct_signed_strain','reflected_signed_strain','absorbing_boundary_strain'):ax.plot(t,[float(x[key]) for x in r],label=key)
    ax.set(xlabel='same-parent age (µs)',ylabel='signed finite-cell strain',title='Separated 5 mm receiver: finite front and signed return, same damping');ax.legend(fontsize=8);fig.tight_layout();fig.savefig(HERE/'first_travel.png',dpi=160);plt.close(fig)
    with (HERE/'optical_phase.csv').open() as f:r=[x for x in csv.DictReader(f) if x['E_cm_keV']=='2.5' and x['arm']==NOMINAL['name'] and x['beam']=='hot_bursts']
    fig,ax=plt.subplots(1,2,figsize=(11,4));t=np.array([float(x['time_since_parent_s']) for x in r])*1e12
    ax[0].plot(t,[float(x['signed_displacement_m'])*1e10 for x in r]);ax[0].set(xlabel='same global parent age (ps)',ylabel='real optical displacement (Å)')
    ax[1].plot(t,[float(x['medium_minus_frozen_phase_rad']) for x in r]);ax[1].set(xlabel='same global parent age (ps)',ylabel='medium − frozen phase (rad)')
    fig.suptitle('G3 dispersive optical receiver, one lifetime downstream; physical sparse source projection');fig.tight_layout();fig.savefig(HERE/'optical_receiver.png',dpi=160);plt.close(fig)
    import views,sys
    views.render(sys.modules[__name__])


def run_case(arm,beam,ec,calendar,tracks,controls=(),resolve=True):
    """One joined entrypoint; persist settled macro, window and phase stages."""
    import inspect,carriers,acoustic
    # Cache policy is operational, not a change to the identical physical ray
    # operator. Preserve settled stages when only eligible kernel reuse changes.
    macro_physics=inspect.getsource(macro).replace('acoustic.cached_operator(','acoustic.Operator(')
    def physics_digest(name):
        if name=='acoustic.py':return hashlib.sha256((HERE/name).read_bytes().split(b'# Operational kernel reuse;')[0]).hexdigest()
        return src.digest(HERE/name)
    field_signature=hashlib.sha256((macro_physics+inspect.getsource(source_probabilities)+''.join(
        physics_digest(name) for name in ('native.py','fields.py','carriers.py','acoustic.py','retarded.py','response.py','sources.py','launch.py'))).encode()).hexdigest()
    cache=HERE/f'.case_E{ec:g}_{beam}_{arm["name"]}.json'
    result=json.loads(cache.read_text()) if cache.exists() else {}
    if result.get('field_signature')!=field_signature:
        result=dict(field_signature=field_signature,controls=[],energy=[],bands=[],checks=[],travel=[],fine=[],diagnostics=[],optical=[],moments=[],shuffled=[])
    kernels,travel=carrier_kernel(arm);result['travel']=travel;context=None;snap=None
    for control in ('joined',)+tuple(controls):
        if control in result['controls']:continue
        energy,bands,check,current,ctx=macro(arm,beam,calendar,kernels,tracks,control,ec)
        result['energy'].extend(energy);result['bands'].extend(bands);result['checks'].append(check);result['controls'].append(control)
        if control=='joined':snap=current;context=ctx
        else:del ctx;__import__('gc').collect()
        emit(cache.name,result)
    fine_signature=hashlib.sha256((field_signature+inspect.getsource(fine_window)+src.digest(HERE/'prehistory.py')).encode()).hexdigest()
    if resolve and ec==2.5 and result.get('fine_signature')!=fine_signature:
        if snap is None:snap=np.load(field_path(beam,arm['name'],ec))
        if context is None:
            coefficient=native.constitutive(arm,ec);carrier=carriers.Operator(arm,ec,kernels,tracks,coefficient['species_residence_ratio'])
            source=spatial_source(beam);rate=src.current(ec)/QE;base=carrier.field(source,rate)
            radiator=acoustic.cached_operator(arm,ec,beam,tracks);bath=radiator.bath_snapshot(None,1.)
            short=carrier.heat_rate(source,rate);short/=short.sum();loss,_=radiator.heat_shape();loss/=loss.sum()
            power=src.current(ec)*2000*ec;acoustic_power=rate*radiator.total_per_parent
            context=dict(base=base,bath=bath,carrier=carrier,radiator=radiator,coefficient=coefficient,heat_power=(power-acoustic_power)*short+acoustic_power*loss)
        result['fine'],result['diagnostics']=fine_window(arm,beam,calendar,snap,context);result['fine_signature']=fine_signature
        emit(cache.name,result)
    phase_signature=hashlib.sha256((fine_signature+src.digest(HERE/'phase.py')+src.digest(HERE/'envelope.py')).encode()).hexdigest()
    phase_path=HERE/f'conditioned_{beam}_{arm["name"]}.npz'
    phase_ready=False
    if ec==2.5 and phase_path.exists():
        with np.load(phase_path) as saved:phase_ready='phase_ledger_schema' in saved and int(saved['phase_ledger_schema'])==3
    if resolve and ec==2.5 and (result.get('phase_signature')!=phase_signature or not phase_ready):
        import phase,envelope,importlib
        importlib.reload(envelope);importlib.reload(phase)
        result['optical'],result['moments'],result['shuffled']=phase.conditioned(arm,beam,calendar,tracks,mode_parameters())
        result['phase_signature']=phase_signature;emit(cache.name,result)
    print(json.dumps(dict(integrated_case=True,arm=arm['name'],beam=beam,E_cm_keV=ec,resolved_windows=bool(resolve and ec==2.5),
        maximum_relative_energy_defect=max(x['relative_balance'] for x in result['checks']),
        maximum_later_launch_medium_work_eV=max((abs(x['medium_contribution_to_same_impulse_work_eV']) for x in result['moments']),default=0.))),flush=True)
    return result


def first_case():
    native.load()
    calendar=fld.Calendar();tracks=json.loads((HERE/'frozen/tracks.json').read_text())
    result=run_case(NOMINAL,'hot_bursts',2.5,calendar,tracks,('frozen',),resolve=False)
    emit('first_integrated.json',dict(checks=result['checks'],source_rows=native.load()['rows'],scope='two-second joined and matched frozen whole-foil fields; fine phase windows are a separate retained stage'))
    return result


def all_fields(closure=None):
    import workflow,sys
    return workflow.all_fields(sys.modules[__name__],closure=closure)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--prepare');parser.add_argument('--refresh-optical',action='store_true');parser.add_argument('--early',action='store_true');parser.add_argument('--all',action='store_true');parser.add_argument('--consume',action='store_true');parser.add_argument('--case',action='store_true');parser.add_argument('--native',action='store_true');parser.add_argument('--batch',choices=[a['name'] for a in src.PHYSICAL_ARMS]);args=parser.parse_args()
    if args.prepare:src.prepare(args.prepare)
    if args.refresh_optical:src.refresh_optical()
    if args.early:early()
    if args.native:native.prepare()
    if args.case:first_case()
    if args.all:all_fields()
    if args.batch:all_fields(args.batch)
    if args.consume:
        mixture,checks=consume();emit('cold_consume.json',dict(success=True,rows=len(mixture),checks=checks,
            source_index_sha256=src.digest(HERE/'source_index.json'),source_member_hashes_verified=all(src.digest(HERE/r['path'])==r['sha256'] for r in json.loads((HERE/'source_index.json').read_text())['source_files'])))
        print(json.dumps(dict(cold_consume=True,mixture_rows=len(mixture),physical_velocity_kernel=True)),flush=True)
if __name__=='__main__':main()
