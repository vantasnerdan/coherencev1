"""Retarded weak-target Green used at actual later nuclear encounters.
This off-grid D-mode approximation is not a replay of the nonlinear motifs.
Energetic mechanical marks and states outside the bound chart remain explicit.
"""
import hashlib,math
import numpy as np
from mechanics import ForceHistory
import whole
_CACHE={};_REG={};_DIRECT={}
_FLIGHT_GL,_FLIGHT_GW=np.polynomial.legendre.leggauss(12)


def matched_shift(m,ec,sp,s,env,branch,xyz):
    key=(ec,sp,s['name'],env)
    if key not in _REG:
        with np.load(m.HERE/f'peak_{ec:g}_{sp+1}_{s["name"]}_{env}.npz',allow_pickle=False) as d:
            pos=d['position_A'];weight=d['occupied_weight']
        nb=len(pos)//int(np.prod(m.SHAPE));rng=np.random.default_rng(152027+int(ec*100)+sp)
        shifts=rng.uniform(0,2*np.pi/m.OMEGA,size=(len(m.BRANCHES),nb)+m.SHAPE).reshape(len(m.BRANCHES),-1)
        keep=weight>0;_REG[key]=(pos[keep],shifts[:,keep])
    pos,shifts=_REG[key];nearest=np.sum((xyz[:,None,:]-pos[None,:,:])**2,axis=2).argmin(axis=1)
    b=next(i for i,(name,_) in enumerate(m.BRANCHES) if name==branch)
    return shifts[b,nearest]


def components(m,ec,sp,s,xyz,at,env,full_parent,geometry,control,branch):
    xyz=np.asarray(xyz,float);at=np.broadcast_to(at,(len(xyz),)).astype(float)
    if control==1:at=at-matched_shift(m,ec,sp,s,env,branch,xyz)
    key=(ec,sp,s['name'],env,full_parent,geometry,hashlib.sha256(xyz.tobytes()+at.tobytes()).hexdigest())
    if key in _CACHE:return _CACHE[key]
    qg=np.zeros_like(xyz);vg=np.zeros_like(xyz);qp=np.zeros_like(xyz);vp=np.zeros_like(xyz)
    qm=np.zeros_like(xyz);vm=np.zeros_like(xyz);strong=np.zeros(len(xyz),bool);T=np.zeros(len(xyz))
    tr=m.p1.track(ec,sp);delay=.74/math.sqrt(2*tr['initial_eV']/m.M) if sp and geometry=='axial' else 0.
    constituents=2 if sp and full_parent else 1
    direct_key=(ec,sp,s['name'],env)
    if direct_key not in _DIRECT:
        _DIRECT[direct_key]=m.response_map(ec,sp,s,env,0.,'axial',False)['direct_scale']
    direct=_DIRECT[direct_key]
    for constituent in range(constituents):
        offset=np.array([.74,0.,0.]) if constituent and geometry=='transverse' else np.zeros(3)
        points=xyz-offset;time=at-constituent*delay
        u,gradient,t0,gclock=m.electronic_coordinate(ec,sp,s,points)
        inside=(points[:,2]>=0)&(points[:,2]<=tr['timeline'][-1]['depth_A'])&(np.sum(points[:,:2]**2,axis=1)<(s['radius']*1e10)**2)
        F=gradient*inside[:,None];G=-u[:,None]*gclock*inside[:,None]
        q,v=whole.finite_motion(time-t0,F,np.zeros_like(F),m.M,m.OMEGA,
            1/(s['tau_opt']*1e15),1/(s['tau_opt']*1e15),1/(s['tau_ne']*1e15),G)
        qg+=q;vg+=v
        # The radial cage-force coefficient is unmeasured; scalar susceptibility
        # response is kept separately for all six original site-family vectors.
        polar=np.c_[u/.1,np.zeros((len(u),2))]*inside[:,None]
        q,v=whole.finite_motion(time-t0,polar,np.zeros_like(F),m.M,m.OMEGA,
            1/(s['tau_opt']*1e15),1/(s['tau_opt']*1e15),1/(s['tau_ne']*1e15))
        qp[:,0]+=q[:,0];vp[:,0]+=v[:,0]
        for start in range(0,len(xyz),128):
            sl=slice(start,min(start+128,len(xyz)))
            history=ForceHistory(m,ec,sp,points[sl],env,direct)
            q,v=history.motion(time[sl],m.OMEGA,1/(s['tau_opt']*1e15))
            qm[sl]+=q;vm[sl]+=v;strong[sl]|=history.strong;T[sl]=np.maximum(T[sl],history.T)
    result=dict(q_gradient=qg,v_gradient=vg,q_polar=qp[:,0],v_polar=vp[:,0],q_mechanical=qm,v_mechanical=vm,
        energetic_mechanical_mark=strong,candidate_transfer_eV=T)
    _CACHE[key]=result
    return result


def state(m,ec,sp,s,branch,xyz,at,env='bulk',full_parent=False,geometry='axial',site=0,control=0):
    d=components(m,ec,sp,s,xyz,at,env,full_parent,geometry,control,branch)
    chi=dict(m.BRANCHES)[branch]
    q=d['q_mechanical']+chi*d['q_gradient'];v=d['v_mechanical']+chi*d['v_gradient']
    if env=='vacancy':
        native=np.array(m.VREF['positions_A'])[63:69];normal=native-native.mean(0);normal/=np.linalg.norm(normal,axis=1)[:,None]
        q=q+chi*d['q_polar'][:,None]*normal[site];v=v+chi*d['v_polar'][:,None]*normal[site]
    action=np.sum(q*q+m.M/m.K*v*v,axis=1);cap=.36 if env=='bulk' else .24
    valid=(action<=cap*cap)&~d['energetic_mechanical_mark']
    return dict(q=q,v=v,supported=valid,action_A2=action,
        energetic_mechanical_mark=d['energetic_mechanical_mark'],candidate_transfer_eV=d['candidate_transfer_eV'])


def sibling(m,ec,s,branch,ch,U,theta,env,control):
    tr=m.p1.track(ec,1);ts,z,dep,_=m.p1.source(tr,1);z=z*1e10;ts=ts*1e15
    delay=.74/math.sqrt(2*tr['initial_eV']/m.M)
    output={}
    for geometry in ('axial','transverse'):
        value=0.;unsupported=0.;total=0.;maximum=0.
        for j in range(len(dep)):
            zz=(z[j]+z[j+1])/2;lab=(tr['initial_eV'] if j==0 else tr['timeline'][j-1]['energy_eV'])/1000
            xyz=np.array([[0.,0.,zz]]) if geometry=='axial' else np.array([[.74,0.,zz]])
            at=np.interp(zz,z,ts)+(delay if geometry=='axial' else 0.)
            column=.5*m.profile_nd(ec,np.array([zz]))[0]*(z[j+1]-z[j])*1e-10
            site_values=[];site_valid=[]
            for site in range(6 if env=='vacancy' else 1):
                d=state(m,ec,1,s,branch,xyz,at,env,False,geometry,site,control)
                delta=float(m.signed_rate(lab,d['v'],ch,theta,U)[0]-m.p1.moving_rate(lab,theta,ch,U)) if d['supported'][0] else 0.
                site_values.append(delta);site_valid.append(d['supported'][0])
                maximum=max(maximum,float(np.linalg.norm(d['v'][0])))
            value+=column*np.mean(site_values);unsupported+=column*(1-np.mean(site_valid));total+=column
        output[geometry]=value
        output[geometry+'_unsupported_column_fraction']=unsupported/max(total,1e-100)
        output[geometry+'_maximum_candidate_speed_A_fs']=maximum
        output[geometry+'_total_increment_available']=output[geometry+'_unsupported_column_fraction']==0.
    output['meaning']='supported weak first-field increment only; energetic same-ray collisional exchange is unresolved, not zero force'
    return output


def struck(m,ec,sp,s,branch,ch,U,theta,env,control):
    tr=m.p1.track(ec,sp);ts,z,dep,_=m.p1.source(tr,sp);z=z*1e10;ts=ts*1e15
    delta=0.;baseline=0.;supported_baseline=0.;unsupported=0.;minimum=0.;maximum=0.;speed=0.
    gw=_FLIGHT_GW;fraction=(_FLIGHT_GL+1)/2
    groups=np.array_split(np.arange(len(dep)),min(12,len(dep)))
    for group in groups:
        j=group[len(group)//2];zz=(z[j]+z[j+1])/2;lab=tr['initial_eV'] if j==0 else tr['timeline'][j-1]['energy_eV']
        layer=next(r for r in m.p1.PROFILES[ec]['layers'] if r[1]*1e10>=zz);nh,nd=layer[2:4]
        for transfer in (.05,.2,.5,.9):
            T=lab*transfer;count=float(dep[group,2].sum())/(4*T)
            stopping=(nh*(m.p1.bca.se_ls(T,46)+m.p1.bca.sn_zbl(T,46,106.42))+
                nd*(m.p1.bca.se_ls(T,1)+m.p1.bca.sn_zbl(T,1,m.p1.bca.M_D)))*1e-18
            R=T/stopping;distance=R*fraction;energy=T*(1-fraction)
            direction=np.array([math.sqrt(1-transfer),0.,math.sqrt(transfer)])
            xyz=np.array([0.,0.,zz])+distance[:,None]*1e10*direction
            flight=2*R*1e10/np.sqrt(2*T/m.M)*(1-np.sqrt(1-fraction));at=np.interp(zz,z,ts)+flight
            projectile=np.sqrt(2*energy/m.M)[:,None]*direction
            quiet=m.p1.moving_rate(energy/1000,theta,ch,U);weight=count*m.profile_nd(ec,xyz[:,2])*R/2*gw
            baseline+=float(weight@quiet);minimum=min(minimum,float(flight.min()));maximum=max(maximum,float(flight.max()))
            changes=[];validities=[]
            for site in range(6 if env=='vacancy' else 1):
                d=state(m,ec,sp,s,branch,xyz,at,env,True,'axial',site,control)
                valid=d['supported'];rate=quiet.copy()
                relative=projectile[valid]-d['v'][valid];effective=.5*m.M*np.sum(relative*relative,axis=1)/1000
                rate[valid]=m.p1.moving_rate(effective,theta,ch,U)*np.sqrt(effective/(energy[valid]/1000))
                changes.append(float(weight@(rate-quiet)))
                validities.append(d['supported']);speed=max(speed,float(np.linalg.norm(d['v'],axis=1).max()))
            delta+=np.mean(changes);valid=np.mean(validities,axis=0)
            supported_baseline+=float(weight@(quiet*valid));unsupported+=float(weight@(quiet*(1-valid)))
    return dict(delta=delta/(sp+1),recomputed_baseline=baseline/(sp+1),supported_baseline=supported_baseline/(sp+1),
        unsupported_baseline=unsupported/(sp+1),flight_min_fs=minimum,flight_max_fs=maximum,target_speed_A_fs=speed,
        model='same four energy-transfer quadrature marks, direction cos=sqrt(T/E), finite LS/ZBL range and analytic flight clock',
        matched_control='same source seed/clock shifts, nearest occupied coarse-regulator site off-grid; no fitted phase parameter',
        total_signed_increment_available=unsupported==0.,
        scope='supported signed target-field correction; energetic/off-chart target continuation is separate and unresolved')
