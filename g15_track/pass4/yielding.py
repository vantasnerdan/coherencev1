"""Last consumer: positive sparse-population mixtures of actual velocity covariances.
No nuclear contact probability is assigned. The fixed G12 sheet is quiet once.
"""
import csv,json,math
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import erf
import sources as s
import native
import fields as f


def square_weights(name):
    edges=np.linspace(-s.L/2,s.L/2,f.N+1);lo=np.maximum(edges[:-1],-.005);hi=np.minimum(edges[1:],.005)
    if name=='even_steady':return np.outer(np.maximum(hi-lo,0.),np.maximum(hi-lo,0.))/1e-8
    terms=((0.,0.,1.,.003),(.0025,-.0025,5.,.0012),(-.0025,.0025,3.,.0012));out=np.zeros((f.N,f.N))
    for ax,ay,ag,aa in terms:
        for bx,by,bg,bb in terms:
            variance=aa*aa*bb*bb/(aa*aa+bb*bb);sigma=np.sqrt(variance)
            cx=(ax*bb*bb+bx*aa*aa)/(aa*aa+bb*bb);cy=(ay*bb*bb+by*aa*aa)/(aa*aa+bb*bb)
            amplitude=ag*bg*np.exp(-((ax-bx)**2+(ay-by)**2)/(2*(aa*aa+bb*bb)))
            ix=np.where(hi>lo,sigma*np.sqrt(np.pi/2)*(erf((hi-cx)/(np.sqrt(2)*sigma))-erf((lo-cx)/(np.sqrt(2)*sigma))),0.)
            iy=np.where(hi>lo,sigma*np.sqrt(np.pi/2)*(erf((hi-cy)/(np.sqrt(2)*sigma))-erf((lo-cy)/(np.sqrt(2)*sigma))),0.)
            out+=amplitude*np.outer(ix,iy)
    return out/s.beam_moments(name)[0]**2


def clock_second(calendar):
    a=calendar.edges[:-1];b=calendar.edges[1:];w=calendar.omega;p=calendar.phase;T=calendar.seconds
    cosine=(np.sin(w*T+p)-np.sin(p))/w
    cosine2=T/2+(np.sin(2*w*T+2*p)-np.sin(2*p))/(4*w)
    on=calendar.state;burst=float(np.sum((b-a)*on));mixed=float(np.sum((np.sin(w*b+p)-np.sin(w*a+p))*on)/w)
    return (T+1.4*cosine+.49*cosine2+99*burst+12.6*mixed)/(calendar.norm**2*T)


class Kernel:
    """Stable connected difference of the *same* physical Maxwell quadrature."""
    def __init__(self,model,lab,theta,minimum,maximum,branch,U):
        self.lab=np.asarray(lab);self.reference=np.asarray(model.moving_rate(lab,theta,branch,U))
        lower=theta+minimum;upper=theta+maximum
        if lower<=0:raise ValueError('Nonpositive actual target velocity covariance')
        if upper<=lower:upper=lower*(1+1e-6)
        self.nodes=np.geomspace(lower,upper,385)
        ratio=model.moving_rate(self.lab[:,None],self.nodes[None,:],branch,U)/self.reference[:,None]
        self.spline=CubicSpline(self.nodes,ratio,axis=1)
    def difference(self,j,theta,delta):
        start,end=np.broadcast_arrays(np.asarray(theta),np.asarray(theta)+np.asarray(delta));d=end-start
        nodes=self.nodes;c=self.spline.c
        ia=np.minimum(np.maximum(np.searchsorted(nodes,start,side='right')-1,0),len(nodes)-2)
        ib=np.minimum(np.maximum(np.searchsorted(nodes,end,side='right')-1,0),len(nodes)-2)
        def local(i,x,h):return h*(c[2,i,j]+c[1,i,j]*(2*x+h)+c[0,i,j]*(3*x*x+3*x*h+h*h))
        x=start-nodes[ia];out=local(ia,x,d)
        different=ia!=ib
        if np.any(different):
            boundary=np.where(d>0,nodes[np.minimum(ia+1,len(nodes)-1)],nodes[ia])
            first=local(ia,x,boundary-start);second=local(ib,boundary-nodes[ib],end-boundary)
            nearby=abs(ia-ib)==1
            xe=end-nodes[ib]
            coarse=((c[0,ib,j]*xe+c[1,ib,j])*xe+c[2,ib,j])*xe+c[3,ib,j]-(((c[0,ia,j]*x+c[1,ia,j])*x+c[2,ia,j])*x+c[3,ia,j])
            out=np.where(different,np.where(nearby,first+second,coarse),out)
        return out*self.reference[j]


def trajectory(cache,arm,p):
    t=np.r_[0.,cache['time_s']];original=cache['k_m_inverse'];gamma=int(np.argmin(np.linalg.norm(original,axis=1)))
    indices=np.r_[np.delete(np.arange(len(original)),gamma),np.full(6,gamma)]
    k=original[indices].copy();k[-6:]=np.r_[np.eye(3),-np.eye(3)]/arm['radius']
    mode_weight=np.r_[np.full(len(original)-1,1/len(original)),np.full(6,1/(6*len(original)))]
    _,v,D,_=native.dispersion(p,k,arm['allocation'])
    C=2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4'];mass=p['mass_eVc2']*s.QE/299792458.**2
    _,Q0,_,_=native.thermal_pole(p,D,300.);sel=cache['selected_mode_indices'];Ds=cache['base_K'][sel]
    _,Qs,_,_=native.thermal_pole(p,Ds,300.);nearest=np.argmin(abs(D[:,None]-Ds[None,:]),axis=1)
    Q=np.vstack((Q0,Q0[None,:]+cache['variance_q_selected_A2'][:,nearest]-Qs[nearest][None,:]))
    if np.any(Q<=0):raise ValueError('Selected covariance transport projection is not positive')
    delta=np.r_[0.,cache['projected_deltaK_eVA2']];shell=np.asarray(p['shells'][0]['periodic_alias_offsets_A'])*1e-10
    F=np.cos(k@shell.T).mean(axis=1);gradF=-np.einsum('kb,bi->ki',np.sin(k@shell.T),shell)/len(shell)
    sensitivity=np.ones(len(k)) if arm['allocation']=='onsite' else F
    omega=np.sqrt((D[None,:]+delta[:,None]*sensitivity+C*Q)*s.QE/1e-20/mass)
    baseomega=np.sqrt(D*s.QE/1e-20/mass);gradient=2*mass*baseomega[:,None]*v/(s.QE/1e-20)
    gradient=gradient[None,:,:]+(delta[:,None,None]*gradF[None,:,:] if arm['allocation']=='coupling' else 0.)
    group=gradient*s.QE/1e-20/(2*mass*omega[:,:,None])
    shift=np.vstack((np.zeros((1,len(k))),np.cumsum((group[:-1,:,2]+group[1:,:,2])/2*np.diff(t)[:,None],axis=0)))
    wt=np.r_[np.diff(t)[0]/2,(np.diff(t)[:-1]+np.diff(t)[1:])/2,np.diff(t)[-1]/2]
    return t,shift,wt[:,None]*mode_weight[None,:],indices


def population_histogram(delta,weights,bins=192):
    scale=float(abs(delta).max(initial=0.))
    if scale==0:return np.array([0.]),np.array([weights.sum()])
    positive=np.geomspace(max(scale*1e-12,1e-18),scale*(1+1e-12),bins)
    edges=np.r_[-positive[::-1],0.,positive];index=np.minimum(np.maximum(np.searchsorted(edges,delta,side='right')-1,0),len(edges)-2)
    mass=np.bincount(index.ravel(),weights=weights.ravel(),minlength=len(edges)-1)
    first=np.bincount(index.ravel(),weights=(weights*delta).ravel(),minlength=len(edges)-1);ok=mass>0
    return first[ok]/mass[ok],mass[ok]


def consume(model,write=True):
    d=s.inputs();p=d['lattice'];tracks=json.loads((s.HERE/'frozen/tracks.json').read_text());calendar=f.Calendar()
    residual={float(r['E_cm_keV']):r for r in d['residual']};reference=float(residual[2.5]['model_corrected_raw_reactions_per_true_C'])/float(residual[2.5]['R_model'])
    rows=[];checks=[];phase_rows=[];thermal_clock={b:clock_second(calendar) if b=='hot_bursts' else 1. for b in s.BEAMS}
    for arm in s.ARMS:
        for ec in s.ENERGIES:
            profile=next(r for r in d['profiles'] if r['case']['E_cm_keV']==ec);rr=residual[ec]
            required=(float(rr['R_observed'])-float(rr['R_model']))*reference*s.QE/1.15
            layers=np.array(profile['layers'])
            def density(z):return layers[np.minimum(np.searchsorted(layers[:,1],z,side='right'),len(layers)-1),3]
            for receiver in (0,1):
                track=tracks[f'{ec:g}_{receiver}'];z=np.r_[0.,[r['depth_A']*1e-10 for r in track['timeline']]];mid=(z[:-1]+z[1:])/2
                energy=np.r_[track['initial_eV'],[r['energy_eV'] for r in track['timeline']]]/1000;lab=(energy[:-1]+energy[1:])/2
                nd=density(mid);column=nd*np.diff(z);source_data={}
                for mesh in (16,12):
                    for sp in (0,1):
                        cache=np.load(s.HERE/f'native_E{ec:g}_{arm["name"]}_sp{sp}_m{mesh}.npz');t,shift,wt,indices=trajectory(cache,arm,p)
                        length=float(tracks[f'{ec:g}_{sp}']['terminal_depth_A']);source_data[mesh,sp]=(cache,t,shift,wt,length,indices)
                # Geometry/time covariance quadratures do not depend on beam
                # temperature or nuclear branch. Build them once, not 24 times.
                statistics={};all_minimum=0.;all_maximum=0.
                for (mesh,sp),(cache,t,shift,wt,length,indices) in source_data.items():
                    for coupling,key in (('native_on','velocity_theta_excess_eV'),('native_e_off','velocity_theta_no_e_excess_eV')):
                        delta=np.vstack((np.zeros((1,len(indices))),cache[key][:,indices]))
                        minimum=float(delta.min());maximum=float(delta.max());histograms=[]
                        all_minimum=min(all_minimum,minimum);all_maximum=max(all_maximum,maximum)
                        for j,zz in enumerate(mid):
                            pre=zz-shift;image=-zz-shift
                            source_density=np.where((pre>=0)&(pre<=length),density(pre),0.)+np.where((image>=0)&(image<=length),density(image),0.)
                            histograms.append(population_histogram(delta,wt*source_density/nd[j]))
                        residence=float(np.trapezoid(abs(cache['optical_excess_eV']),cache['time_s']))/float(abs(cache['optical_excess_eV']).max())
                        statistics[mesh,sp,coupling]=(histograms,minimum,maximum,float(abs(shift).max()),residence)
                for beam in s.BEAMS:
                    fields=np.load(model.field_path(beam,arm['name'],ec));fraction=fields['source_fraction'][:,:,0].ravel();square=square_weights(beam).ravel();active=fraction>0
                    theta=s.KB*(300+fields['t2_deltaT'][:,:,0].ravel()[active]);area=square[active];quiet_theta=s.KB*profile['case']['T_K']
                    for branch in ('p','n'):
                        for U,label in ((0.,'bare'),(.04,'static40')):
                            quiet=model.moving_rate(lab,quiet_theta,branch,U)
                            quiet_y=float(quiet@column)+profile['nads_m2']*float(model.moving_rate(energy[0],quiet_theta,branch,U))
                            thermal=float(fraction[active]@((model.moving_rate(lab[None,:],theta[:,None],branch,U)-quiet)@column))
                            kernel=Kernel(model,lab,float(theta.mean()),all_minimum+float(theta.min()-theta.mean()),all_maximum+float(theta.max()-theta.mean()),branch,U)
                            for mesh in (16,12):
                                for coupling in ('native_on','native_e_off'):
                                    addition=0.;weighted_occupation=0.;maxshift=0.
                                    for sp in (0,1):
                                        histograms,minimum,maximum,displacement,residence=statistics[mesh,sp,coupling]
                                        maxshift=max(maxshift,displacement);encounter=0.
                                        for j,(values,mass) in enumerate(histograms):
                                            difference=kernel.difference(j,theta[:,None],values[None,:])
                                            encounter+=float(area@(difference@mass))*column[j]
                                        source_rate=s.current(ec)/s.QE*(.85 if sp==0 else .15)
                                        addition+=source_rate*np.pi*arm['radius']**2*thermal_clock[beam]*encounter
                                        weighted_occupation+=source_rate*np.pi*arm['radius']**2*area.sum()*thermal_clock[beam]*residence
                                        if ec==2.5 and beam=='even_steady' and receiver==0 and mesh==16 and coupling=='native_on' and sp==0 and branch=='p' and U==.04:
                                            extreme=maximum;j=0;actual=float(kernel.difference(j,theta.mean(),extreme))
                                            direct=float(model.moving_rate(lab[j],theta.mean()+extreme,branch,U)-model.moving_rate(lab[j],theta.mean(),branch,U))
                                            error=abs(actual-direct)/max(abs(direct),abs(float(quiet[j]))*1e-12)
                                            checks.append(dict(arm=arm['name'],kind='actual_covariance_kernel',maximum_relative_connected_error=error,
                                                source_covariance_eV=extreme,mesh=mesh,retarded_maximum_depth_shift_m=maxshift))
                                            if error>1e-3:raise ValueError('Actual covariance kernel changes the consumer')
                                    nw=(.85 if receiver==0 else .30)/1.15
                                    rows.append(dict(E_cm_keV=ec,arm=arm['name'],beam=beam,mesh=mesh,source_coupling=coupling,receiver_species=receiver+1,branch=branch,comparator=label,
                                        G12_quiet_temperature_K=profile['case']['T_K'],built_beam_temperature_min_K=float(theta.min()/s.KB),built_beam_temperature_max_K=float(theta.max()/s.KB),
                                        quiet_per_incident_D=quiet_y,transported_heat_delta_per_D=thermal,other_parent_optical_velocity_delta_per_D=addition,
                                        mixture_contribution_per_D=nw*(thermal+addition),G12_required_extra_per_D=required,fraction_of_residual=nw*(thermal+addition)/required,
                                        G12_raw_own_reference_per_true_C=reference,G12_observed_R=float(rr['R_observed']),G12_model_R=float(rr['R_model']),
                                        arrival_conditioned_energy_weighted_occupation=weighted_occupation,maximum_retarded_optical_depth_shift_m=maxshift,
                                        fixed_sheet_once=True,ordinary_quiet_subtracted_once=True,sheet_temperature_unchanged=True,
                                        kernel='compound-Poisson sparse mixture of 64 actual Gaussian velocity covariances; signed time histogram; reciprocal G3 transport; all beam-cell temperatures',
                                        source_scope='mesh/coupling velocity contrast on primary built fields and transported geometry; not an independently evolved no-e ray; no artificial contact probability',low_energy_continuation='BH below0.5keV conditional, not measured'))
    groups={}
    for row in rows:
        key=tuple(row[k] for k in ('E_cm_keV','arm','beam','mesh','source_coupling','branch','comparator'));groups.setdefault(key,[]).append(row)
    mixtures=[]
    for key,rr in groups.items():
        value=sum(r['mixture_contribution_per_D'] for r in rr);need=rr[0]['G12_required_extra_per_D']
        mixtures.append(dict(E_cm_keV=key[0],arm=key[1],beam=key[2],mesh=key[3],source_coupling=key[4],branch=key[5],comparator=key[6],
            built_field_addition_per_D=value,fraction_of_G12_residual=value/need,G12_required_extra_per_D=need,
            supported_inflight='unchanged source-owned struck-D support in frozen prior table',native_neighbour_contact='unresolved; no nuclear probability assigned'))
    for arm in s.ARMS:
        for beam in s.BEAMS:
            cached=np.load(s.HERE/f'conditioned_{beam}_{arm["name"]}.npz');history=cached['history'][0]
            for row in history:
                phase_rows.append(dict(E_cm_keV=2.5,arm=arm['name'],beam=beam,time_since_parent_s=row[0],signed_displacement_m=row[1],complex_envelope_real_m=row[1],complex_envelope_imag_m=row[2],
                    source_x_m=float(cached['source_positions_m'][0,0]),source_y_m=float(cached['source_positions_m'][0,1]),source_z_m=float(cached['source_positions_m'][0,2]),
                    receiver_x_m=float(cached['receiver_positions_m'][0,0]),receiver_y_m=float(cached['receiver_positions_m'][0,1]),receiver_z_m=float(cached['receiver_positions_m'][0,2]),
                    global_parent_arrival_s=float(cached['global_arrival_s'][0]),group_centroid_arrival_age_s=float(cached['group_centroid_arrival_age_s'][0]),
                    signed_velocity_m_s=row[3],medium_minus_frozen_phase_rad=row[4],group_speed_m_s=row[5],instantaneous_K_eV_A2=row[6],omega_rad_s=row[7],
                    coherent_energy_eV=row[8],dephased_incoherent_remnant_eV=row[9],quiet_cross_eV=row[10],incoherent_remainder_eV=row[11],
                    electron_to_coherent_work_eV=row[15],incoherent_optical_to_coherent_work_eV=row[16],causal_heat_correction_eV=row[17],conditional_packet_balance_eV=row[18],physical_packet_loss_eV=row[20],
                    energy_scope='packet energy, not per-D energy or receiver density; total optical budget enters Gaussian covariance consumer once',
                    velocity_scope='canonical optical oscillator momentum/mass, not derivative of an open-system first moment',
                    correlation_bandwidth_Hz=float(cached['correlation_bandwidth_Hz']),supported_remote_slow_h2_over_threshold=0.,
                    pump_scope='supported remote acoustic frequencies below 2Omega; DC/slow shifts are not a resonant pump'))
    if write:
        model.table('yield_consumers.csv',rows);model.table('yield_mixture.csv',mixtures);model.table('optical_phase.csv',phase_rows);model.emit('consumer_checks.json',checks)
    return mixtures,checks
