"""G15 Pass 1: one conservative conditional impact-to-foil model.
Run: python3 g15_track/pass1/run.py --all
Only this directory is written. No inherited producer is executed.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse, csv, hashlib, json, math, sys, time
from pathlib import Path
from functools import lru_cache
sys.dont_write_bytecode=True
import numpy as np
from scipy.linalg import eigh, expm
from scipy.integrate import solve_ivp
from scipy.special import roots_genlaguerre, roots_legendre, gamma
import bca

HERE=Path(__file__).resolve().parent
QE=1.602176634e-19; MD=2.0141*1.66053906660e-27
KB=8.617333262145e-5; HB=6.582119569e-16
RHO=12007.; CP=25.89/.1064; CV=RHO*CP; KTH=71.8
ALPHA=11.77e-6; BULK=(227.1e9+2*176.04e9)/3
L=np.sqrt(np.pi)*.024/2; H=.00025; AREA=1e-4
CHANNELS=['moving_cutoff_D','moving_Pd_recoil','moving_struck_D','nonthermal_electrons','D_optical','Pd_incoherent_low','Pd_incoherent_high','Pd_coherent','thermal_electrons','local_lattice','stopped_structural']
SCENARIOS=[dict(name='fast_acoustic',opt_share=.25,tau_ne=15e-15,tau_opt=.0802e-12,tau_pdlo=.3e-12,tau_pdhi=.15e-12,tau_te=.3e-12,tau_recoil=.05e-12,tau_end=.1e-12,tau_struct=1.3e-9,radius=8e-9,De=.001),
 dict(name='nominal_mixed',opt_share=.5,tau_ne=.1e-12,tau_opt=.5e-12,tau_pdlo=2e-12,tau_pdhi=1e-12,tau_te=.6e-12,tau_recoil=.2e-12,tau_end=.5e-12,tau_struct=100e-9,radius=20e-9,De=.001),
 dict(name='optical_moment_long',opt_share=.92,tau_ne=.5e-12,tau_opt=2.068018471814257e-12,tau_pdlo=10e-12,tau_pdhi=5e-12,tau_te=.3e-12,tau_recoil=1e-12,tau_end=2e-12,tau_struct=1e-6,radius=60e-9,De=.001)]
CONTACTS=[dict(name='weak_mount_faces',edge=100.,face=2.,emissivity=.05),dict(name='nominal_mount',edge=1e4,face=5.,emissivity=.1),dict(name='strong_mount_faces',edge=1e6,face=100.,emissivity=.3)]
BEAMS=['even_steady','hot_steady','hot_bursts']
DATA=json.loads((HERE/'inputs.json').read_text())
PROFILES={float(r['case']['E_cm_keV']):r for r in DATA['profiles']}


def emit(name,obj):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    (HERE/name).write_text(json.dumps(obj,indent=2,allow_nan=False,default=convert)+'\n')


def table(name,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


@lru_cache(maxsize=None)
def track(ec,species):
    """Exact central saved endpoints; same accepted LS/ZBL equations otherwise."""
    saved=next((r for r in DATA['tracks'] if r['E_cm_keV']==ec),None)
    if saved:return saved['tracks'][species]
    energy=2000*ec/(species+1); initial=energy; dep=np.zeros(3);t=0.;rows=[]
    for zt,zb,nh,nd,z,m,_ in PROFILES[ec]['layers']:
        if energy<=100.00001:break
        def rhs(zz,y):
            e=max(y[0],100.)
            split=np.array([nh*bca.se_ls(e,z)+nd*bca.se_ls(e,1),nh*bca.sn_zbl(e,z,m),nd*bca.sn_zbl(e,1,bca.M_D)])*1e-18
            return np.r_[-split.sum(),split,1e15/np.sqrt(2*e*QE/MD)]
        def stop(zz,y):return y[0]-100.
        stop.terminal=True;stop.direction=-1
        sol=solve_ivp(rhs,(zt,zb),np.r_[energy,dep,t],events=stop,rtol=2e-7,atol=1e-8,max_step=2e-9)
        if not sol.success:raise RuntimeError(sol.message)
        energy=float(sol.y[0,-1]);dep=sol.y[1:4,-1];t=float(sol.y[4,-1]);depth=float(sol.t[-1])
        rows.append(dict(depth_A=depth*1e10,time_fs=t,energy_eV=energy,electron_eV=dep[0],host_recoil_eV=dep[1],D_recoil_eV=dep[2]))
    return dict(species='D+' if species==0 else 'one D from D2+',initial_eV=initial,terminal_unresolved_eV=energy,terminal_time_fs=t,terminal_depth_A=depth,components_eV=dep.tolist(),timeline=rows)


def source(tr,species):
    multiplicity=species+1
    t=np.r_[0.,[r['time_fs']*1e-15 for r in tr['timeline']]]
    z=np.r_[0.,[r['depth_A']*1e-10 for r in tr['timeline']]]
    dep=np.diff(np.array([[0.,0.,0.]]+[[r[k] for k in ('electron_eV','host_recoil_eV','D_recoil_eV')] for r in tr['timeline']]),axis=0)*multiplicity
    # The saved central components close to their own projectile energy, not the
    # five-point beam-average ledger. Do not rescale a track to that other mean.
    coh=float(next(r for r in DATA['coherent'] if r['scenario']=='track' and r['species']==('D+' if species==0 else 'D2+'))['acoustic_J'])/QE
    coh*=tr['initial_eV']*multiplicity/5000
    f=coh/dep[:,1].sum()
    u=np.zeros((len(dep),len(CHANNELS)));u[:,3]=dep[:,0];u[:,1]=dep[:,1]*(1-f);u[:,2]=dep[:,2];u[:,7]=dep[:,1]*f
    return t,z,dep,u


def network(s):
    A=np.zeros((len(CHANNELS),len(CHANNELS)))
    def send(j,tau,parts):
        A[j,j]=-1/tau
        for k,f in parts:A[k,j]+=f/tau
    a=s['opt_share']
    send(0,s['tau_end'],[(3,.3),(4,.25),(6,.4),(10,.05)])
    send(1,s['tau_recoil'],[(4,.05),(5,.35),(6,.6)])
    send(2,s['tau_recoil']*.73,[(3,.25),(4,.5),(6,.25)])
    send(3,s['tau_ne'],[(4,a),(5,(1-a)*.3),(6,(1-a)*.7)])
    send(4,s['tau_opt'],[(8,.2),(6,.8)])
    send(5,s['tau_pdlo'],[(9,1.)]);send(6,s['tau_pdhi'],[(9,1.)])
    send(7,50e-6,[(9,1.)]);send(8,s['tau_te'],[(9,1.)])
    send(9,1e-9,[]);send(10,s['tau_struct'],[(9,1.)])
    return A


def spatial_axis(edges,conductivity,left,right):
    widths=np.diff(edges);centers=(edges[1:]+edges[:-1])/2;n=len(widths)
    cap=CV*widths;G=np.zeros((n,n))
    for j in range(n-1):
        g=conductivity/(centers[j+1]-centers[j]);G[j,j]+=g;G[j+1,j+1]+=g;G[j,j+1]-=g;G[j+1,j]-=g
    # Contact and half-cell conduction are in series.
    gl=1/(1/max(left,1e-300)+widths[0]/(2*conductivity))
    gr=1/(1/max(right,1e-300)+widths[-1]/(2*conductivity))
    G[0,0]+=gl;G[-1,-1]+=gr
    eigval,V=eigh(G/np.sqrt(cap[:,None]*cap[None,:]))
    return dict(x=centers,w=widths,cap=cap,gamma=eigval,V=V,left=gl,right=gr)


def geometry(contact):
    rad=4*contact['emissivity']*5.670374419e-8*300**3
    x=spatial_axis(np.linspace(-L/2,L/2,18),KTH,contact['edge'],contact['edge'])
    z=spatial_axis(np.r_[0.,.25e-6,1e-6,5e-6,20e-6,60e-6,125e-6,190e-6,230e-6,H],KTH,contact['face']+rad,contact['face']+rad)
    # Orthonormal modes in sqrt(volume) coordinates, energy density / CV.
    U=np.einsum('ia,jb,kc->ijkabc',x['V'],x['V'],z['V']).reshape(17*17*9,-1)
    volume=np.einsum('i,j,k->ijk',x['w'],x['w'],z['w']).ravel()
    gam=(x['gamma'][:,None,None]+x['gamma'][None,:,None]+z['gamma'][None,None,:]).ravel()
    X,Y,Z=np.meshgrid(x['x'],x['x'],z['x'],indexing='ij')
    ex=np.zeros_like(X);ez=np.zeros_like(X)
    ex[0]+=x['left']*x['w'][:,None]*z['w'][None,:];ex[-1]+=x['right']*x['w'][:,None]*z['w'][None,:]
    ex[:,0]+=x['left']*x['w'][:,None]*z['w'][None,:];ex[:,-1]+=x['right']*x['w'][:,None]*z['w'][None,:]
    ez[:,:,0]=(contact['face']+rad)/(1+(contact['face']+rad)*z['w'][0]/(2*KTH))*x['w'][:,None]*x['w'][None,:]
    ez[:,:,-1]=(contact['face']+rad)/(1+(contact['face']+rad)*z['w'][-1]/(2*KTH))*x['w'][:,None]*x['w'][None,:]
    exits=np.array([ex.ravel(),ez.ravel()*contact['face']/(contact['face']+rad),ez.ravel()*rad/(contact['face']+rad)])
    return dict(U=U,vol=volume,gamma=gam,X=X,Y=Y,Z=Z,exits=exits,x=x,z=z,contact=contact)


def beam_weight(x,y,beam):
    inside=(abs(x)<=.005)&(abs(y)<=.005)
    if beam=='even_steady':return inside.astype(float)
    w=np.exp(-(x*x+y*y)/(2*.003**2))
    for a,b,g in ((.0025,-.0025,5.),(-.0025,.0025,3.)):
        w+=g*np.exp(-((x-a)**2+(y-b)**2)/(2*.0012**2))
    return w*inside


def distribution(g,tr,beam):
    t,z,dep,_=source(tr,0)
    pz=np.zeros(9)
    for mid,e in zip((z[:-1]+z[1:])/2,dep.sum(axis=1)):
        pz[min(np.searchsorted(np.r_[0.,.25e-6,1e-6,5e-6,20e-6,60e-6,125e-6,190e-6,230e-6,H],mid)-1,8)]+=e
    pz[0]+=tr['terminal_unresolved_eV'];pz/=pz.sum()
    xy=beam_weight(g['X'][:,:,0],g['Y'][:,:,0],beam)*np.outer(g['x']['w'],g['x']['w']);xy/=xy.sum()
    p=(xy[:,:,None]*pz).ravel()
    return p,(g['U'].T@(p/np.sqrt(g['vol'])))/CV


def propagator(A,gam):
    """Confluent-safe exact affine/thermal Green, including repeated poles."""
    n=len(A);I=np.eye(n);Ai=np.linalg.inv(A);c=np.zeros(n);c[9]=1e9
    shifted=A[None,:,:]+gam[:,None,None]*I
    R=np.linalg.solve(np.swapaxes(shifted,1,2),np.broadcast_to(c,(len(gam),n))[...,None])[...,0]
    last=[None,None]
    def green(dt):
        if dt==last[0]:return last[1]
        if np.max(abs(np.diag(A)))*dt<.5:
            # Entire-function series avoid small-dt cancellation and cover
            # equal carrier poles without changing any physical rate.
            ad=A*dt;E=I.copy();F=I*dt;term=I.copy()
            for j in range(1,19):
                term=term@ad/j;E+=term;F+=term*dt/(j+1)
            S=np.broadcast_to(c,(len(gam),n)).copy()
            G=S*dt;J=S*dt*dt/2;power=np.ones(len(gam))
            coefficient=dt
            for j in range(1,19):
                power*=-gam*dt;S=S@ad+power[:,None]*c
                coefficient/=j+1;G+=coefficient*S;J+=coefficient*dt/(j+2)*S
        else:
            E=expm(A*dt);F=Ai@(E-I)
            im=-np.expm1(-gam*dt)/gam
            G=R@E-np.exp(-gam*dt)[:,None]*R
            J=(G-im[:,None]*c)@Ai
        value=(E,F,G,J)
        last[:]=[dt,value]
        return value
    return green


def impulse(tr,species,s,g):
    t,z,dep,src=source(tr,species);A=network(s);green=propagator(A,g['gamma'])
    _,b=distribution(g,tr,'even_steady')
    # A centered single spatial cell, with the same vertical source and operator.
    p=np.zeros(g['X'].shape);p[8,8,0]=1
    b=QE*g['U'].T@(p.ravel()/np.sqrt(g['vol']))/CV
    y=np.zeros(len(A));theta=np.zeros(len(b));ex=np.zeros(3);time_now=0.;qsource=0.;rows=[];all_states=[];all_times=[];mode_samples=[]
    exitmode=g['exits']@(g['U']/np.sqrt(g['vol'])[:,None])/QE
    energy_mode=CV*np.sqrt(g['vol'])@g['U']/QE
    def advance(dt,force):
        nonlocal y,theta,ex,time_now
        E,F,G,J=green(dt)
        newy=E@y+F@force
        B=G@y+J@force
        old=theta.copy();theta=theta*np.exp(-g['gamma']*dt)+b*B
        # Exactly c.integral(E dt), from this generator's sole foil port.
        heat_int=float(sum(y)+sum(force)*dt-sum(newy))
        tint=(old+b*heat_int-theta)/g['gamma']
        tiny=g['gamma']*dt<1e-7
        tint[tiny]=(old[tiny]+theta[tiny])*dt/2
        ex+=exitmode@tint
        y=newy;time_now+=dt
    def save():
        stored=float(energy_mode@theta)
        # Integrate boundary exits independently; roundoff is retained as defect.
        row=dict(scenario=s['name'],species='D+' if species==0 else 'D2+',time_s=time_now,incoming_eV=max(0.,tr['initial_eV']*(species+1)-qsource),source_cumulative_eV=qsource,foil_lattice_eV=stored,exit_mount_eV=ex[0],exit_faces_eV=ex[1],exit_radiation_eV=ex[2])
        row.update(zip(CHANNELS,y));row['balance_eV']=row['incoming_eV']+sum(y)+stored+sum(ex)-tr['initial_eV']*(species+1)
        rows.append(row);all_states.append(y.copy());all_times.append(time_now)
        if time_now>=1e-9:mode_samples.append((time_now,theta.copy()))
    save()
    delay=.74e-10/np.sqrt(2*tr['initial_eV']*QE/MD) if species else 0.
    offsets=(0.,delay) if species else (0.,)
    edges=np.unique(np.concatenate([t+offset for offset in offsets]))
    terminal=tr['terminal_unresolved_eV']*(species+1)
    for ta,tb in zip(edges[:-1],edges[1:]):
        force=np.zeros(len(A))
        for offset in offsets:
            mid=(ta+tb)/2-offset
            j=np.searchsorted(t,mid,side='right')-1
            if 0<=j<len(src):force+=src[j]/((species+1)*(t[j+1]-t[j]))
        n=max(1,int(np.ceil((tb-ta)/.25e-15)));dt=(tb-ta)/n
        for _ in range(n):
            advance(dt,force);qsource+=force.sum()*dt;save()
        for offset in offsets:
            if abs(tb-(t[-1]+offset))<1e-28:
                y[0]+=terminal/(species+1);qsource+=terminal/(species+1);save()
    for target in np.geomspace(max(time_now*1.01,1e-15),300.,460):
        advance(target-time_now,np.zeros(len(A)));save()
        if sum(y)+float(energy_mode@theta)<tr['initial_eV']*(species+1)*1e-8:break
    # Integrated residence from the exact same source and stable inverse DAG.
    total=src.sum(axis=0);total[0]+=terminal
    residence=np.linalg.solve(-A,total)
    return dict(rows=rows,t=np.array(all_times),y=np.array(all_states),residence=residence,total=total,A=A,source_t=t,source_z=z,source_u=src,theta=mode_samples,b=b,exitmode=exitmode,energy_mode=energy_mode)


def thermal_fields(g,impulses):
    rows=[];arrays={};U=g['U'];gam=g['gamma'];vol=g['vol'];tempfac=U/np.sqrt(vol)[:,None]
    for beam in BEAMS:
        for s in SCENARIOS:
            for species in (0,1):
                imp=impulses[s['name'],species];tr=track(2.5,species);p,b=distribution(g,tr,beam)
                # All excess energy reaches this same spatial source; residence
                # upstream changes standing carrier stocks, not total DC power.
                power=138e-6/QE*(.85 if species==0 else .15)*5000*QE
                steady=tempfac@(b*power/gam)
                A=network(s);green=propagator(A,gam)
                force=imp['total']*power/5000
                for snap in (.001,.1,1.,10.,300.):
                    _,_,_,J=green(snap)
                    modal=b*(J@force)
                    field=tempfac@modal
                    arrays[f'{beam}_{s["name"]}_D{species+1}_t{snap:g}']=field.reshape(g['X'].shape)
                arrays[f'{beam}_{s["name"]}_D{species+1}_steady']=steady.reshape(g['X'].shape)
                # Source-row arrival widths are <0.5ps; macroscopic beam
                # snapshots use the same integrated parent impulse and operator.
                exits=g['exits']@steady
                grad=np.gradient(steady.reshape(g['X'].shape),g['x']['x'],g['x']['x'],g['z']['x'],edge_order=1)
                rows.append(dict(boundary=g['contact']['name'],beam=beam,scenario=s['name'],species=species+1,power_W=power,steady_max_delta_K=steady.max(),steady_min_delta_K=steady.min(),mean_delta_K=steady@vol/vol.sum(),max_lateral_gradient_K_m=max(abs(grad[0]).max(),abs(grad[1]).max()),max_depth_gradient_K_m=abs(grad[2]).max(),free_max_linear_strain=ALPHA*steady.max(),fully_clamped_compressive_hydrostatic_Pa=-3*BULK*ALPHA*steady.max(),plane_clamped_inplane_Pa=-120e9/(1-.39)*ALPHA*steady.max(),exit_mount_W=exits[0],exit_faces_W=exits[1],exit_radiation_W=exits[2],steady_balance_W=exits.sum()-power))
                for k,name in enumerate(CHANNELS):
                    density=(p/vol)*(imp['residence'][k]*QE*138e-6/QE*(.85 if species==0 else .15))
                    arrays[f'{beam}_{s["name"]}_D{species+1}_{name}_J_m3']=density.reshape(g['X'].shape)
    arrays.update(x_m=g['x']['x'],z_m=g['z']['x'],volume_m3=vol.reshape(g['X'].shape),ambient_K=np.array(300.))
    np.savez_compressed(HERE/f'foil_{g["contact"]["name"]}.npz',**arrays)
    return rows


def nd_column(ec,length):
    return sum(max(0.,min(zb,length)-zt)*nd for zt,zb,nh,nd,*_ in PROFILES[ec]['layers'])


def short_value(imp,time_s,k):
    return float(np.interp(time_s,imp['t'],imp['y'][:,k],left=0.,right=0.))


def local_and_sibling(impulses):
    rows=[];sibling=[]
    for s in SCENARIOS:
        for species in (0,1):
            imp=impulses[s['name'],species];tr=track(2.5,species);length=tr['timeline'][-1]['depth_A']*1e-10
            column=nd_column(2.5,length);r=s['radius'];N=column*np.pi*r*r
            for ti in (1e-15,10e-15,100e-15,1e-12,10e-12,100e-12,1e-9,1e-6):
                for rad in (0.,r,3*r):
                    width=min(r,2e-9+min(1e6*ti,np.sqrt(4*s['De']*ti)))
                    kernel=np.exp(-rad*rad/(2*width*width))/(2*np.pi*width*width*length)
                    for k in (1,2,3,4,5,6,7,8,9,10):
                        rows.append(dict(scenario=s['name'],species=species+1,time_s=ti,radius_m=rad,channel=CHANNELS[k],excess_eV_m3=short_value(imp,ti,k)*kernel,optical_sites_within_radius=N if k==4 else '',source_depth_m=length,radial_kernel_sigma_m=width,condition='retarded local positive-energy packet; electron spreading bounded by1e6m/s; no nonthermal temperature'))
            # A constituent is the first half of D2 parent's linear response.
            if species==1:
                v=np.sqrt(2*tr['initial_eV']*QE/MD);delay=.74e-10/v
                t,z,dep,src=source(tr,species)
                A=network(s)
                for orientation,dt,offset in [('longitudinal',delay,0.),('transverse',0.,.74e-10)]:
                    for frac in (.25,.5,.75,1.):
                        zz=length*frac;arrival=np.interp(zz,z,t)+dt
                        energy=np.zeros(len(A))
                        for j in range(len(src)):
                            start=t[j];end=min(t[j+1],arrival)
                            if end<=start:continue
                            age=arrival-(start+end)/2
                            width=2e-9+min(1e6*age,np.sqrt(4*s['De']*age))
                            dz=(z[j]+z[j+1])/2-zz
                            kernel=np.exp(-(offset*offset+dz*dz)/(2*width*width))/(2*np.pi*width*width)**1.5
                            amount=src[j]/2*(end-start)/(t[j+1]-t[j])
                            energy+=(expm(A*age)@amount)*kernel
                        nd=next(layer[3] for layer in PROFILES[2.5]['layers'] if layer[1]>=zz)
                        sibling.append(dict(scenario=s['name'],orientation=orientation,launch_separation_A=.74,charged_gas_ground_state_proxy_A=1.06,alternate_axial_delay_fs=1.06e-10/v*1e15,velocity_m_s=v,arrival_delay_fs=dt*1e15,depth_m=zz,first_constituent_nonthermal_eV_m3=energy[3],first_constituent_optical_eV_D=energy[4]/nd,first_constituent_structural_eV_m3=energy[10],first_constituent_thermal_lattice_eV_m3=energy[9],second_self_response_included=False,molecule_at_site='dissociated constituent scenario, not intact molecule'))
    table('local_space_time.csv',rows);table('second_constituent_prior.csv',sibling)


def cox_times(T,rng):
    cuts=[0.];states=[rng.random()<.1]
    while cuts[-1]<T:
        cuts.append(cuts[-1]+rng.exponential(2e-6 if states[-1] else 18e-6));states.append(not states[-1])
    return np.array(cuts),np.array(states)


def fluctuations(impulses):
    """Finite real-flux receiver windows and event-triggered packet observations."""
    rng=np.random.default_rng(151127);rows=[];marks=[];arrays={};records=64
    xx=np.linspace(-.005,.005,301);X,Y=np.meshgrid(xx,xx,indexing='ij')
    for s in SCENARIOS:
        r=s['radius'];disk=np.pi*r*r
        for beam in BEAMS:
            norm=beam_weight(X,Y,beam).mean()*AREA
            for receiver,(rx,ry) in [('center',(0.,0.)),('hotspot',(.0025,-.0025)),('edge',(.0045,.0045))]:
                density=float(beam_weight(np.array(rx),np.array(ry),beam))/norm
                lam=138e-6/QE*disk*density
                T=256/(records*lam)
                counts=[];uniform=[];triggered=[];signed=[];correlated=[];shuffled=[];overlaps=0
                expected=lam*T*records
                for rec in range(records):
                    cuts,state=cox_times(T,rng) if beam=='hot_bursts' else (None,None)
                    phase=rng.uniform(0,2*np.pi)
                    maximum=10.7/1.9 if beam=='hot_bursts' else 1.
                    n=rng.poisson(lam*T*maximum);times=np.sort(rng.uniform(0,T,n))
                    intensity=np.ones(n)
                    if beam=='hot_bursts':
                        intensity=(1+.7*np.cos(2*np.pi*1e5*times+phase)+9*state[np.searchsorted(cuts,times,side='right')-1])/1.9
                        times=times[rng.random(n)<intensity/maximum]
                    counts.append(len(times));species=rng.random(len(times))>=.85
                    overlaps+=int(np.count_nonzero(np.diff(times)<8*s['tau_opt']+1e-12))
                    # Uniform observation times include all empty draws. Future
                    # event fields are never used to determine arrival intensity.
                    observations=rng.uniform(0,T,128)
                    for at in observations:
                        E=0.
                        for sp in (0,1):
                            ages=at-times[species==sp];ages=ages[(ages>=0)&(ages<200*s['tau_opt']+2e-12)]
                            E+=sum(short_value(impulses[s['name'],sp],age,4) for age in ages)
                        length=track(2.5,0)['timeline'][-1]['depth_A']*1e-10
                        uniform.append(E/(nd_column(2.5,length)*disk))
                    for it,sp in enumerate(species):
                        imp=impulses[s['name'],int(sp)];length=track(2.5,int(sp))['timeline'][-1]['depth_A']*1e-10
                        N=nd_column(2.5,length)*disk
                        # Random point in the actual mark footprint and event age.
                        age=rng.uniform(0,8*s['tau_opt']+1e-12);rr=r*np.sqrt(rng.random())
                        amp=short_value(imp,age,4)/N
                        triggered.append(amp);ph=rng.uniform(0,2*np.pi)
                        x=np.sqrt(max(0,2*amp*QE/(MD*(.037/HB)**2)))*np.cos(ph)
                        signed.append(x);correlated.append(x)
                        shuffled.append(np.sqrt(max(0,2*amp*QE/(MD*(.037/HB)**2)))*np.cos(rng.uniform(0,2*np.pi)))
                        marks.append(dict(scenario=s['name'],beam=beam,receiver=receiver,record=rec,time_s=times[it],species=int(sp)+1,radius_m=rr,age_s=age,sites_in_actual_profile=N,optical_eV_D=amp,signed_displacement_m=x,phase_rad=ph,prior='exogenous Cox; event response applied after arrival'))
                uniform=np.array(uniform);triggered=np.array(triggered)
                # Campbell mean and variance of marked positive population, not
                # inferred absence of rare events from sparse uniform draws.
                mean=lam*sum(f*impulses[s['name'],sp]['residence'][4]/(nd_column(2.5,track(2.5,sp)['timeline'][-1]['depth_A']*1e-10)*disk) for sp,f in [(0,.85),(1,.15)])
                occupancy=lam*(8*s['tau_opt']+1e-12)
                q=np.quantile(triggered,[0,.5,.95,.99,1]) if len(triggered) else np.zeros(5)
                rows.append(dict(scenario=s['name'],beam=beam,receiver=receiver,rate_per_s=lam,records=records,window_s=T,total_receiver_seconds=T*records,expected_parents=expected,observed_parents=sum(counts),empty_records=sum(c==0 for c in counts),uniform_draws=len(uniform),uniform_empty_draws=int(np.count_nonzero(uniform==0)),uniform_mean_eV_D=uniform.mean(),campbell_mean_eV_D=mean,finite_window_peak_eV_D=max(uniform.max(),q[-1]),event_conditioned_min=q[0],event_conditioned_median=q[1],event_conditioned_q95=q[2],event_conditioned_q99=q[3],event_conditioned_max=q[4],geometric_overlap_probability=-np.expm1(-occupancy),marked_overlap_count=overlaps,signed_min_m=min(signed,default=0),signed_max_m=max(signed,default=0),holes='zero optical excess, not negative occupation; signed troughs exported separately',phase='declared random phase; D2 siblings share parent, shuffled contrast not a measured phase'))
                arrays[f'{s["name"]}_{beam}_{receiver}_uniform_eV_D']=uniform
                arrays[f'{s["name"]}_{beam}_{receiver}_conditional_eV_D']=triggered
                arrays[f'{s["name"]}_{beam}_{receiver}_signed_m']=np.array(signed)
                arrays[f'{s["name"]}_{beam}_{receiver}_common_parent_m']=np.array(correlated)
                arrays[f'{s["name"]}_{beam}_{receiver}_shuffled_m']=np.array(shuffled)
    table('fluctuations.csv',rows);table('packet_marks.csv',marks)
    np.savez_compressed(HERE/'finite_windows.npz',**arrays)
    return rows


def ratios(impulses,fluct):
    rows=[];p=DATA['lattice'];K=p['Kii_eV_A2'];quartic=p['K4_100_eV_A4'];omega=np.sqrt(K*QE/1e-20/MD)
    for s in SCENARIOS:
        kappa=1/s['tau_opt'];threshold=2*kappa/omega
        for sp in (0,1):
            imp=impulses[s['name'],sp];length=track(2.5,sp)['timeline'][-1]['depth_A']*1e-10
            N=nd_column(2.5,length)*np.pi*s['radius']**2
            maxe=max(imp['y'][:,4])/N
            for locked in (0.,.01,1.):
                # q=A cos(Omega t), q^2=A^2[1+cos(2Omega t)]/2;
                # delta K=3 K4 q^2; only the locked phasor contributes at2Omega.
                A2=2*maxe/(K);h=3*quartic*A2/(2*K)*locked
                strain=ALPHA*1. # conditional one-kelvin thermal strain, DC only
                rows.append(dict(scenario=s['name'],species=sp+1,mode='same-Kii harmonic local projection; band37meV stock is separate',quantum_eV=omega*HB,frequency_Hz=omega/(2*np.pi),energy_loss_s_inverse=kappa,amplitude_loss_s_inverse=kappa/2,omega_over_energy_loss=omega/kappa,f_over_energy_loss=omega/(2*np.pi*kappa),max_packet_E_over_kBT=maxe/(KB*300),max_packet_quanta_D=maxe/(omega*HB),nominal_pump_quanta_D_s=138e-6/QE*(.85 if sp==0 else .15)*imp['residence'][4]*kappa/(N*omega*HB),phase_locked_fraction=locked,anharmonic_fractional_stiffness_at_2omega=h,mathieu_resonant_threshold=threshold,h_over_threshold=h/threshold,neighbour_over_well=sum(r['J_sum_eV_A2'] for r in p['shells'])/K,quartic_ZPM_over_K=quartic*p['zero_point_axis_variance_A2']/K,thermal_DC_strain_per_K=strain,thermal_fractional_stiffness_per_K=p['collective_strain_slope_eV_A2']*strain/K,acoustic_resonant_modulation=0.,acoustic_reason='Pd support ends6.8THz below twice optical frequency; DC strain is not pump',critical_strain_coupling_per_unit_strain=threshold/strain,conditional_critical_kelvin_for_DC_shift=K/abs(p['collective_strain_slope_eV_A2']*ALPHA),threshold_scope='conditional all optical energy in same-Kii local Mathieu mode; lock fraction explicit; broadband stock not phase locked'))
    table('project_ratios.csv',rows)


def sigma(E,ch='p',U=.04):
    E=np.asarray(E,float);out=np.zeros_like(E);ok=E>0;e=E[ok]
    coeff=([5.5576e4,2.1054e2,-3.2638e-2,1.4987e-6,1.8181e-10] if ch=='p' else [5.3701e4,3.3027e2,-1.2706e-1,2.9327e-5,-2.5151e-9])
    S=sum(c*e**i for i,c in enumerate(coeff))
    out[ok]=S/np.sqrt(e*(e+U))*np.exp(-31.397/np.sqrt(e+U))*1e-31
    return out


def moving_rate(lab,theta,ch,U=.04):
    """Thermal quiet-background relative flux, vectorized, in square metres."""
    lab=np.asarray(lab,float)
    x,w=roots_genlaguerre(8,.5);mu,mw=roots_legendre(8)
    e=lab[...,None,None]/2
    er=e+x[:,None]*theta/2000-2*np.sqrt(e*x[:,None]*theta/2000)*mu[None,:]
    return np.sum(sigma(er,ch,U)*np.sqrt(er/e)*w[:,None]/gamma(1.5)*mw[None,:]/2,axis=(-2,-1))
 
 
def pump_rate(lab,energy_site,ch,quiet_theta,U):
    """Random oscillator phase/direction, not fixed kinetic energy E/2.
    K=E sin(phi)^2 at each event, hence <K>=E/2. Quiet thermal motion is
    convolved once about the resulting relative velocity, not added as energy.
    """
    lab,energy_site=np.broadcast_arrays(lab,energy_site)
    q,w=roots_legendre(8);phi=(q+1)*np.pi/4
    kinetic=energy_site[...,None,None]*np.sin(phi[:,None])**2/1000
    effective=lab[...,None,None]+kinetic-2*np.sqrt(lab[...,None,None]*kinetic)*q[None,:]
    out=moving_rate(effective,quiet_theta,ch,U)*np.sqrt(effective/lab[...,None,None])
    return np.sum(out*w[:,None]*w[None,:]/4,axis=(-2,-1))
 
 
def optical_spectrum(ec,sp,s):
    """Residence-conservative time histogram of the actual causal carrier ODE."""
    tr=track(ec,sp);t,z,dep,src=source(tr,sp);A=network(s);Ai=np.linalg.inv(A)
    green=propagator(A,np.array([1.]))
    y=np.zeros(len(A));energies=[];durations=[]
    delay=.74e-10/np.sqrt(2*tr['initial_eV']*QE/MD) if sp else 0.
    offsets=(0.,delay) if sp else (0.,)
    edges=np.unique(np.concatenate([t+offset for offset in offsets]))
    def step(dt,force):
        nonlocal y
        E,F,_,_=green(dt);new=E@y+F@force
        integral=(Ai@(new-y-force*dt))[4]
        if integral>0:
            energies.append(integral/dt);durations.append(dt)
        y=new
    for ta,tb in zip(edges[:-1],edges[1:]):
        force=np.zeros(len(A))
        for off in offsets:
            j=np.searchsorted(t,(ta+tb)/2-off,side='right')-1
            if 0<=j<len(src):force+=src[j]/((sp+1)*(t[j+1]-t[j]))
        n=max(1,int(np.ceil((tb-ta)/.25e-15)))
        for _ in range(n):step((tb-ta)/n,force)
        for off in offsets:
            if abs(tb-t[-1]-off)<1e-28:y[0]+=tr['terminal_unresolved_eV']
    previous=edges[-1]
    for at in np.geomspace(previous*1.001,previous+50*max(s['tau_opt'],s['tau_ne'],s['tau_end']),400):
        step(at-previous,np.zeros(len(A)));previous=at
    en=np.array(energies);wt=np.array(durations)
    bins=np.r_[0.,np.geomspace(en.max()*1e-7,en.max()*1.000001,13)]
    index=np.clip(np.searchsorted(bins,en,side='right')-1,0,len(bins)-2)
    durations=np.bincount(index,weights=wt,minlength=len(bins)-1)
    residence=np.bincount(index,weights=wt*en,minlength=len(bins)-1)
    nonzero=durations>0
    means=residence[nonzero]/durations[nonzero]
    total=src.sum(0);total[0]+=tr['terminal_unresolved_eV']*(sp+1)
    exact=np.linalg.solve(-A,total)[4]
    return dict(energy=means,duration=durations[nonzero],exact=exact,integrated=residence.sum(),length=z[-1],sites=nd_column(ec,z[-1])*np.pi*s['radius']**2)


def optical_beam_increment(tr,ec,ch,U,energy,sites,pump_length,quiet_theta):
    """Depth/profile weighted beam encounter; optical reservoir excludes sheet."""
    z=np.r_[0.,[r['depth_A']*1e-10 for r in tr['timeline']]]
    en=np.r_[tr['initial_eV'],[r['energy_eV'] for r in tr['timeline']]]/1000
    labs=[];columns=[]
    for j in range(len(en)-1):
        length=max(0.,min(z[j+1],pump_length)-z[j])
        if not length:continue
        zz=(z[j]+min(z[j+1],pump_length))/2
        nd=next(r[3] for r in PROFILES[ec]['layers'] if r[1]>=zz)
        labs.append((en[j]+en[j+1])/2);columns.append(nd*length)
    if not labs:return np.zeros_like(energy)
    lab=np.array(labs);col=np.array(columns);out=np.zeros(len(energy))
    for j in range(0,len(lab),12):
        ll=lab[j:j+12]
        pumped=pump_rate(ll[None,:],(energy/sites)[:,None],ch,quiet_theta,U)
        quiet=moving_rate(ll,quiet_theta,ch,U)
        out+=(pumped-quiet[None,:])@col[j:j+12]
    return out


def track_yield(tr,ec,ch,theta=0.,U=.04):
    rows=tr['timeline'];z=np.r_[0.,[r['depth_A']*1e-10 for r in rows]];en=np.r_[tr['initial_eV'],[r['energy_eV'] for r in rows]]/1000
    total=float(PROFILES[ec]['nads_m2'])*(moving_rate(en[0],theta,ch,U) if theta else float(sigma(en[0]/2,ch,U)))
    for j in range(len(rows)):
        zz=(z[j]+z[j+1])/2;nd=next(r[3] for r in PROFILES[ec]['layers'] if r[1]>=zz)
        rate=moving_rate((en[j]+en[j+1])/2,theta,ch,U) if theta else float(sigma((en[j]+en[j+1])/4,ch,U))
        total+=nd*(z[j+1]-z[j])*rate
    return total


def sibling_prior(ec,s):
    """First constituent only, sampled before the second source is inserted."""
    tr=track(ec,1);t,z,dep,src=source(tr,1);A=network(s)
    delay=.74e-10/np.sqrt(2*tr['initial_eV']*QE/MD)
    depths=[0.];values=[0.]
    for zz in z[-1]*np.array([.25,.5,.75,1.]):
        arrival=np.interp(zz,z,t)+delay;energy=0.
        for j in range(len(src)):
            end=min(t[j+1],arrival)
            if end<=t[j]:continue
            age=arrival-(t[j]+end)/2
            width=2e-9+min(1e6*age,np.sqrt(4*s['De']*age))
            dz=(z[j]+z[j+1])/2-zz
            kernel=np.exp(-dz*dz/(2*width*width))/(2*np.pi*width*width)**1.5
            amount=src[j]/2*(end-t[j])/(t[j+1]-t[j])
            energy+=(expm(A*age)@amount)[4]*kernel
        nd=next(r[3] for r in PROFILES[ec]['layers'] if r[1]>=zz)
        depths.append(zz);values.append(energy/nd)
    return np.array(depths),np.array(values)


def low_energy_table(spectrum):
    """Use the same finite-energy relative-mode distribution as pair consumers."""
    groups={}
    for r in spectrum:
        key=(float(r['E_cm_keV']),r['scenario'],int(r['species']))
        groups.setdefault(key,[]).append(r)
    low=[]
    for (ec,scenario,sp),records in groups.items():
        packet=np.array([float(r['optical_stored_eV']) for r in records])
        duration=np.array([float(r['time_weight_s']) for r in records])
        sites=float(records[0]['sites'])
        for phase,divisor in [('common_phase',0.),('independent_phase',1.),('antiphase',2.)]:
            nu=max(2.,round(sites)/max(divisor,1.))
            for energy in (.001,.01,.1,1.,10.,100.,500.,1000.,5000.):
                positive=packet>energy;survival=np.zeros(len(packet))
                if divisor:survival[positive]=np.exp((nu-1)*np.log1p(-energy/packet[positive]))
                tail=float(survival@duration/sum(duration))
                loglargest=(nu-1)*np.log1p(-energy/max(packet)) if divisor and max(packet)>energy else None
                for ch in ('p','n'):
                    low.append(dict(E_cm_keV=ec,scenario=scenario,species=sp,branch=ch,joint_phase=phase,relative_energy_eV=energy,time_weighted_pair_tail_probability=tail,largest_bin_log_survival=loglargest,conditional_BH_continuation_m2=float(sigma(energy/1000,ch)),zero_below_evaluator_support_m2=float(sigma(energy/1000,ch)) if energy>=500 else 0.,support='finite-energy time histogram; continuation/zero below0.5keV are conditional scenarios, not measured bounds',no_zero_relative_energy_static_evaluation=True))
    table('low_energy_scenarios.csv',low)


def yields():
    rows=[];pairs=[];spectrum=[];checks=[]
    residual={float(r['E_cm_keV']):r for r in DATA['residual']}
    rr0=residual[2.5];denom=float(rr0['model_corrected_raw_reactions_per_true_C'])/float(rr0['R_model'])
    xx=np.linspace(-.005,.005,301);X,Y=np.meshgrid(xx,xx,indexing='ij')
    for ec in (.25,.5,1.,1.5,2.,2.5):
        rr=residual[ec];I=np.interp(2*ec,np.arange(.5,15.01,.5),[38,68,94,110,118,123,130,134,136,138,140,145,152,156,160,165,171,176,180,182,192,207,217,222,228,238,251,263,275,283])*1e-6
        residual_R=float(rr['R_observed'])-float(rr['R_model'])
        required=residual_R*denom*QE/1.15
        for s in SCENARIOS:
            spectra=[optical_spectrum(ec,sp,s) for sp in (0,1)]
            for sp,ss in enumerate(spectra):
                defect=abs(ss['integrated']/ss['exact']-1)
                checks.append(dict(E_cm_keV=ec,scenario=s['name'],species=sp+1,optical_integral_eV_s=ss['integrated'],exact_residence_eV_s=ss['exact'],relative_defect=defect,phase_mean_K_over_E=.5))
                for k,(e,dt) in enumerate(zip(ss['energy'],ss['duration'])):
                    spectrum.append(dict(E_cm_keV=ec,scenario=s['name'],species=sp+1,bin=k,optical_stored_eV=e,time_weight_s=dt,optical_eV_D=e/ss['sites'],sites=ss['sites'],source_depth_m=ss['length'],exact_residence_eV_s=ss['exact'],total_energy_not_all_kinetic=True))
                for env in ('bulk','vacancy'):
                    cuts=DATA['bulk_pair'] if env=='bulk' else DATA['vacancy_pair']['rows']
                    for cut in cuts:
                        if env=='bulk':
                            dmin=float(cut['D_D_separation_A']);barrier=float(cut['total_approach_free_energy_eV']);d0=2.88782409436586
                        else:
                            dmin=cut['separation_A'];barrier=cut['AB_eV'];d0=DATA['vacancy_pair']['rest_separation_A']
                        # A finite-energy microcanonical complex-mode ensemble:
                        # one relative mode energy is Beta(1,nu-1)*E_packet.
                        # Independent site phases:nu=N. Pair antiphase:nu=N/2,
                        # common displacement:relative energy identically zero.
                        for phase,divisor in [('common_phase',0.),('independent_phase',1.),('antiphase',2.)]:
                            nu=max(2.,round(ss['sites'])/max(divisor,1.))
                            accessible=ss['energy']>barrier
                            tail=np.zeros(len(ss['energy']));speed=np.zeros(len(tail))
                            if divisor:
                                tail[accessible]=np.exp((nu-1)*np.log1p(-barrier/ss['energy'][accessible]))
                                # K_rel at the supported inward crossing, not a
                                # new potential continued to nuclear contact.
                                speed[accessible]=np.sqrt(4*(ss['energy'][accessible]-barrier)/nu*QE/MD)
                            weighted=tail*ss['duration']
                            attempts=.037/(2*np.pi*HB)*float(sum(weighted))
                            crossing_speed=float(speed@weighted/sum(weighted)) if sum(weighted)>0 else 0.
                            candidates=ss['sites']/2;source_fraction=.85 if sp==0 else .15
                            need=required*1.15/(source_fraction*candidates*attempts) if attempts>1e-290 else None
                            for eta in (0.,1e-12,1.):
                                probability=-np.expm1(-eta*attempts)
                                for ch in ('p','n'):
                                    for beam in BEAMS:
                                        pairs.append(dict(E_cm_keV=ec,scenario=s['name'],species=sp+1,beam=beam,branch=ch,environment=env,joint_phase=phase,supplied_min_separation_A=dmin,rest_separation_A=d0,computed_energy_at_min_eV=barrier,expected_inward_shell_attempts=attempts,mean_conditional_relative_crossing_speed_m_s=crossing_speed,contact_efficiency=eta,conditional_reactions_per_occupied_candidate=probability,candidate_count_if_all_pumped_sites_paired=candidates,conditional_yield_per_incident_D=probability*candidates/(sp+1),mixture_contribution_per_incident_D=source_fraction*probability*candidates/1.15,required_candidate_abundance_times_efficiency=need,G12_required_extra_per_incident_D=required,geometry='opposed neighbour axis, one inward shell crossing per cycle above computed approach energy; no path below supplied cut',PAW_overlap=bool(env=='vacancy' or dmin<2.),normalization='finite-energy microcanonical relative-mode ensemble, conditional occupied candidate/branch efficiency; vacancy abundance unspecified'))
            quiet_theta=KB*PROFILES[ec]['case']['T_K']
            for receiver_sp in (0,1):
                tr=track(ec,receiver_sp);t,z,dep,src=source(tr,receiver_sp)
                for ch in ('p','n'):
                    for U,label in ((0.,'bare'),(.04,'static40')):
                        quiet=track_yield(tr,ec,ch,quiet_theta,U)
                        background=0.;stock=0.;fraction=0.
                        for source_sp,ss in enumerate(spectra):
                            rate=I/QE*(.85 if source_sp==0 else .15)*np.pi*s['radius']**2/AREA
                            increments=optical_beam_increment(tr,ec,ch,U,ss['energy'],ss['sites'],ss['length'],quiet_theta)
                            background+=rate*float(increments@ss['duration'])
                            stock+=rate*ss['exact']/ss['sites']
                            fraction+=rate*sum(ss['duration'])
                        # Same-parent prior at every nominal energy; never use
                        # the second constituent's post-injection response.
                        sibling=0.
                        if receiver_sp==1:
                            depths,ee=sibling_prior(ec,s)
                            for j in range(len(z)-1):
                                zz=(z[j]+z[j+1])/2;nd=next(r[3] for r in PROFILES[ec]['layers'] if r[1]>=zz);lab=(tr['initial_eV'] if j==0 else tr['timeline'][j-1]['energy_eV'])/1000
                                eprior=np.interp(zz,depths,ee)
                                sibling+=.5*nd*(z[j+1]-z[j])*(float(pump_rate(lab,eprior,ch,quiet_theta,U))-float(moving_rate(lab,quiet_theta,ch,U)))
                        sec=0.;sub=0.;secondary_energy=0.;count_total=0.
                        for j,row in enumerate(dep):
                            projectile=tr['initial_eV'] if j==0 else tr['timeline'][j-1]['energy_eV']
                            zz=(z[j]+z[j+1])/2;layer=next(r for r in PROFILES[ec]['layers'] if r[1]>=zz);nh,nd=layer[2:4]
                            for u in (.05,.2,.5,.9):
                                T=projectile*u;count=row[2]/(4*T)
                                stopping=(nh*(bca.se_ls(T,46)+bca.sn_zbl(T,46,106.42))+nd*(bca.se_ls(T,1)+bca.sn_zbl(T,1,bca.M_D)))*1e-18
                                R=T/stopping;path=np.linspace(0,R,25);ener=np.maximum(T*(1-path/R),0)/2000
                                sig=sigma(ener,ch,U)
                                sec+=count*np.trapezoid(nd*sig,path)
                                sub+=count*np.trapezoid(nd*sig*(ener<.5),path)
                                secondary_energy+=count*T;count_total+=count
                        for beam in BEAMS:
                            w=beam_weight(X,Y,beam);spatial=float(np.mean(w*w)/np.mean(w)**2)
                            temporal=(1+.7**2/2+18*.1+81*.1)/1.9**2 if beam=='hot_bursts' else 1.
                            factor=spatial*temporal;addition=background*factor+sibling
                            nuclear_weight=(.85 if receiver_sp==0 else .30)/1.15
                            rows.append(dict(E_cm_keV=ec,scenario=s['name'],species=receiver_sp+1,branch=ch,comparator=label,beam=beam,per_incident_D_quiet=quiet,per_incident_D_beam_moving_addition=addition,per_incident_D_other_parent_optical=background*factor,per_incident_D_same_parent_sibling=sibling,mixture_contribution_per_incident_D=nuclear_weight*addition,per_incident_D_struck_D=sec/(receiver_sp+1),struck_D_supported_only_per_D=(sec-sub)/(receiver_sp+1),struck_D_subsupport_share=sub/sec if sec else 0.,struck_D_energy_eV_per_D=secondary_energy/(receiver_sp+1),struck_D_mean_count_per_D=count_total/(receiver_sp+1),arrival_conditioned_intensity_factor=factor,conditional_occupied_fraction=fraction*factor,actual_mean_optical_eV_D=stock*factor,G12_observed_R=float(rr['R_observed']),G12_model_R=float(rr['R_model']),G12_residual_R=residual_R,G12_raw_own_reference_per_true_C=denom,G12_required_extra_per_incident_D=required,added_beam_fraction_of_residual=addition/required,profile_exposure_s=PROFILES[ec]['case']['beam_time_s'],fixed_sheet_once=True,pump_sheet='none: optical reservoir assigned to bulk occupied sites',sibling_scope='first-constituent prior at this same fixed-profile energy',velocity='random oscillator phase and direction, thermal quiet convolved once; mean kinetic E/2',low_energy_rule='BH continuation below0.5keV conditional; supported-only struck contribution retained'))
    table('yield_consumers.csv',rows);table('pair_consumers.csv',pairs);table('optical_spectra.csv',spectrum);low_energy_table(spectrum);table('consumer_spectrum_join.csv',checks)
    proof=dict(phase_mean_kinetic_over_total=.5,maximum_optical_residence_relative_defect=max(r['relative_defect'] for r in checks),cases=len(checks),velocity_phase_nodes=8,velocity_direction_nodes=8,thermal_background_nodes=[8,8],spectrum_bins=13,profile_and_beam_intensity_weighted=True)
    emit('consumer_refinement.json',proof)
    if proof['maximum_optical_residence_relative_defect']>1e-6:raise ValueError('Actual optical spectrum does not join stored inventory')
    print(json.dumps(proof))
    return rows,pairs


def figure(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sel=[r for r in rows if r['scenario']=='nominal_mixed' and r['species']=='D+'];t=np.array([r['time_s'] for r in sel])
    fig,ax=plt.subplots(2,1,figsize=(10,8),sharex=True)
    for name in ['incoming_eV']+CHANNELS+['foil_lattice_eV']:
        val=np.array([r[name] for r in sel]);ax[0].plot(t*1e15,np.maximum(val,1e-15),label=name)
    ax[0].set(xscale='log',yscale='log',ylabel='stored / incoming energy (eV)',ylim=(1e-10,1e4));ax[0].legend(fontsize=6,ncol=3)
    for name in ('source_cumulative_eV','exit_mount_eV','exit_faces_eV','exit_radiation_eV'):
        ax[1].plot(t*1e15,[r[name] for r in sel],label=name)
    ax[1].set(xscale='log',xlabel='time since electrical parent arrival (fs)',ylabel='integrated source or exit (eV)');ax[1].legend(fontsize=8)
    fig.suptitle('One5keV D+ impact: positive carrier stocks and separate cumulative flows')
    fig.tight_layout();fig.savefig(HERE/'single_impact_channels.png',dpi=170);plt.close(fig)


def joined_beam_histories(g,impulses):
    """Actual same-network/Green convolution of a charge-normalized Cox path."""
    dt=50e-9;steps=2000;rng=np.random.default_rng(153377)
    arrays={'time_s':np.arange(steps+1)*dt};rows=[];bands=[]
    clock=dt*(np.arange(steps)+.5)
    phase=rng.uniform(0,2*np.pi);on=rng.random()<.1
    burst=np.empty(steps)
    for j in range(steps):
        burst[j]=1+.7*np.cos(2*np.pi*1e5*clock[j]+phase)+9*on
        if rng.random()<-np.expm1(-dt/(2e-6 if on else 18e-6)):on=not on
    burst/=burst.mean()
    tempfac=g['U']/np.sqrt(g['vol'])[:,None]
    for beam in BEAMS:
        multiplier=burst if beam=='hot_bursts' else np.ones(steps)
        arrays[f'{beam}_intensity']=multiplier
        for s in SCENARIOS:
            A=network(s);EE,II,G,J=propagator(A,g['gamma'])(dt)
            for sp in (0,1):
                imp=impulses[s['name'],sp];tr=track(2.5,sp);p,b=distribution(g,tr,beam)
                rate=138e-6/QE*(.85 if sp==0 else .15)
                force0=imp['total']*QE*rate
                y=imp['residence']*QE*rate;theta=b*sum(force0)/g['gamma']
                history=[y.copy()];snapshots=[(tempfac@theta).reshape(g['X'].shape)]
                for j,mult in enumerate(multiplier):
                    ff=force0*mult
                    theta=theta*np.exp(-g['gamma']*dt)+b*(G@y+J@ff)
                    y=EE@y+II@ff;history.append(y.copy())
                    if (j+1)%400==0:snapshots.append((tempfac@theta).reshape(g['X'].shape))
                label=f'{beam}_{s["name"]}_D{sp+1}'
                arrays[label+'_carrier_J']=np.array(history);arrays[label+'_deltaT_K']=np.array(snapshots)
                stock=rate*imp['residence'][4]*QE
                influx=stock/s['tau_opt']
                rows.append(dict(beam=beam,scenario=s['name'],species=sp+1,parent_rate_s=rate,nominal_optical_inventory_J=stock,optical_influx_W=influx,identity_inventory_minus_tau_influx_J=stock-s['tau_opt']*influx,history_duration_s=steps*dt,history_dt_s=dt,intensity_mean=multiplier.mean(),optical_history_min_J=min(v[4] for v in history),optical_history_max_J=max(v[4] for v in history),source_insertion='same causal integrated parent; subps row smearing negligible at50ns history sampling',initial='standing beam state; exogenous history not conditioned on resulting field'))
                for k,name in enumerate(CHANNELS):
                    lo,hi=(0.,0.)
                    if k==4:lo,hi=.030/HB/(2*np.pi),.065/HB/(2*np.pi)
                    if k==5:lo,hi=0.,2.5e12
                    if k==6:lo,hi=2.5e12,6.8e12
                    if k==7:lo,hi=5e3,2.5e12
                    bands.append(dict(scenario=s['name'],species=sp+1,channel=name,source_eV=imp['total'][k],visited_flow_eV=-A[k,k]*imp['residence'][k],residence_eV_s=imp['residence'][k],beam_stock_J=imp['residence'][k]*QE*rate,band_lo_Hz=lo,band_hi_Hz=hi,group_support='Pd sin dispersion caps3.4/4.7/6.8THz; incoherent energy not signed acoustic phase' if k in (5,6) else 'optical measured band, not homogeneous linewidth' if k==4 else 'named finite carrier'))
    arrays['snapshot_times_s']=np.arange(6)*400*dt
    np.savez_compressed(HERE/'joined_beam_history.npz',**arrays)
    table('steady_impulse_join.csv',rows);table('carrier_bands.csv',bands)


def moving_particle_readout():
    rows=[]
    for ec in (.25,.5,1.,1.5,2.,2.5):
        for sp in (0,1):
            tr=track(ec,sp);t,z,dep,u=source(tr,sp)
            for s in SCENARIOS:
                v0=np.sqrt(2*tr['terminal_unresolved_eV']*QE/MD)
                for age in (0.,.1,1.,3.,10.,30.):
                    elapsed=age*s['tau_end'];E=tr['terminal_unresolved_eV']*np.exp(-age)
                    distance=2*v0*s['tau_end']*(-np.expm1(-age/2))
                    rows.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],particle='incoming D below accepted100eV cutoff',age_s=elapsed,kinetic_eV=E,speed_m_s=v0*np.exp(-age/2),depth_m=z[-1]+distance,maximum_forward_extension_m=2*v0*s['tau_end'],already_transferred_eV=tr['terminal_unresolved_eV']-E,stopped_fraction=-np.expm1(-age),fate='same end-ion network:30%electrons25%optical40%Pd5%structural; each transferred eV then drains'))
                for j,row in enumerate(dep):
                    energy=(tr['initial_eV'] if j==0 else tr['timeline'][j-1]['energy_eV'])
                    for mass,k,label in ((MD,2,'struck D'),(106.42*1.66053906660e-27,1,'Pd recoil')):
                        Tmax=energy if k==2 else energy*4*MD*mass/(MD+mass)**2
                        T=.2*Tmax;count=row[k]/max(T,1e-300)
                        vv=np.sqrt(2*T*QE/mass);tau=s['tau_recoil']*(.73 if k==2 else 1.)
                        rows.append(dict(E_cm_keV=ec,species=sp+1,scenario=s['name'],particle=label,source_row=j,source_depth_m=(z[j]+z[j+1])/2,mean_mark_energy_eV=T,energy_normalized_count=count,source_recoil_eV=row[k],conditional_forward_cosine=np.sqrt(.2),uniform_azimuth=True,exponential_mean_energy_lifetime_s=tau,maximum_straight_extension_m=2*vv*tau,measured_details='CSDA mean only; conditional kick energy/direction/range, not an observed cascade'))
    table('moving_particles.csv',rows)


def main():
    start=time.perf_counter();g=geometry(CONTACTS[1]);impulses={};rows=[];proof=[]
    for s in SCENARIOS:
        for species in (0,1):
            imp=impulse(track(2.5,species),species,s,g);impulses[s['name'],species]=imp;rows+=imp['rows']
            balance=max(abs(r['balance_eV']) for r in imp['rows'])/5000
            last=imp['rows'][-1]
            proof.append(dict(scenario=s['name'],species=species+1,maximum_relative_energy_defect=balance,remaining_fraction=(sum(imp['y'][-1])+last['foil_lattice_eV'])/5000,last_time_s=last['time_s'],minimum_stored_eV=float(imp['y'].min()),minimum_exit_increment_eV=min(np.diff([r[k] for r in imp['rows']]).min() for k in ('exit_mount_eV','exit_faces_eV','exit_radiation_eV')),source_energy_eV=imp['total'].sum(),steady_residence_identity=float(np.max(abs(imp['A']@imp['residence']+imp['total'])))))
    table('single_channels.csv',rows);table('conservation.csv',proof)
    thermal=[]
    for c in CONTACTS:thermal+=thermal_fields(g if c==CONTACTS[1] else geometry(c),impulses)
    table('whole_foil.csv',thermal)
    joined_beam_histories(g,impulses)
    moving_particle_readout()
    local_and_sibling(impulses);fluct=fluctuations(impulses);ratios(impulses,fluct);yy,pp=yields();figure(rows)
    arrays={}
    for (name,sp),imp in impulses.items():
        arrays[f'{name}_D{sp+1}_time_s']=imp['t'];arrays[f'{name}_D{sp+1}_stored_eV']=imp['y'];arrays[f'{name}_D{sp+1}_residence_eV_s']=imp['residence']
        # Full thermal modal histories permit exact field reconstruction:
        # deltaT=U@theta/sqrt(volume). No negative population interpretation.
        arrays[f'{name}_D{sp+1}_thermal_time_s']=np.array([t for t,q in imp['theta']]);arrays[f'{name}_D{sp+1}_thermal_modes']=np.array([q for t,q in imp['theta']])
    arrays['thermal_Vxy']=g['x']['V'];arrays['thermal_Vz']=g['z']['V'];arrays['thermal_volume_m3']=g['vol'];arrays['thermal_decay_s_inverse']=g['gamma']
    np.savez_compressed(HERE/'single_response.npz',**arrays)
    pd=np.array([r['steady_balance_W'] for r in thermal]);defect=max(p['maximum_relative_energy_defect'] for p in proof)
    result=dict(model='finite conditional positive carrier network plus finite-geometry thermal Green; self-review, not independent experimental validation',scenarios=SCENARIOS,boundaries=CONTACTS,channels=CHANNELS,conservation_max_relative=defect,steady_exit_max_defect_W=float(abs(pd).max()),remaining_tolerance_fraction=1e-8,all_cases=len(proof),thermal_cases=len(thermal),yield_cases=len(yy),pair_conditional_cases=len(pp),finite_receiver_seconds_per_case=[r['total_receiver_seconds'] for r in fluct],source_first_step_max_fs=.25,thermal_conductivity_split=dict(electron_W_m_K=67.8,lattice_W_m_K=4.,total_W_m_K=71.8),prompt_rule='uniform causal injection within accepted row interval; no deposit before source row begins',spatial_resolution=dict(whole_foil_xy_m=L/17,whole_foil_z_edges_m=[0,.25e-6,1e-6,5e-6,20e-6,60e-6,125e-6,190e-6,230e-6,H],local_packet_radii_m=[s['radius'] for s in SCENARIOS]),interpretation='carrier energy and optical occupancy not phase-locked pump; named pair contact efficiencies and vacancy abundance remain conditional',runtime_s=time.perf_counter()-start)
    emit('results.json',result)
    if defect>2e-4 or abs(pd).max()>1e-7:raise ValueError('Joined energy or steady exit decision failed')
    if any(p['remaining_fraction']>1e-8 for p in proof):raise ValueError('Single impact has not drained to requested tolerance')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--all',action='store_true')
    mode.add_argument('--yields-only',action='store_true')
    args=parser.parse_args()
    if args.all:main()
    else:
        prior=json.loads((HERE/'results.json').read_text())
        if prior['conservation_max_relative']>2e-4:raise ValueError('No accepted unchanged carrier fields')
        yy,pp=yields()
        prior.update(yield_cases=len(yy),pair_conditional_cases=len(pp),consumer_refinement='consumer_refinement.json')
        emit('results.json',prior)
