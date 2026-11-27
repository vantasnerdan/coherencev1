"""Pass5 portable CLI. Reuses accepted whole-foil Green and RNG, never writes Pass4.
Actions: launch, receivers, replay, refine, feedback, check, readout, package, upload.
"""
import argparse
import csv
import hashlib
import json
import numpy as np
from scipy.spatial import cKDTree
import source as s

t, r = s.t, s.r
HERE, ROOT = s.HERE, s.ROOT
BEAMS = ('even_steady','hot_steady','hot_bursts')
LAWS = ('Cu_extrapolated','Al_transfer')


def readrows(name):
    with (HERE/name).open() as f:
        return list(csv.DictReader(f))


def angle_lookup(n, data, directions):
    # Four local directional quadrature neighbours; positive energy interpolation.
    distance, index = cKDTree(directions).query(n, k=4)
    weight = 1/np.maximum(distance,1e-12)**2
    weight /= weight.sum(axis=1)[:,None]
    return np.einsum('ni,ni...->n...',weight,data[index])


def source_phase(mean):
    return np.exp(-1j*np.angle(mean))


def make_kernel(xy, band, energy, law, scenario, source, frequency=None, bandwidth=None,
                angular=None, variance=None, geometry=None, center=None):
    base = r.kernel(xy,band,np.ones((2,3)),law,frequency=frequency,bandwidth=bandwidth) if geometry is None else geometry
    scale=energy[:,np.repeat(np.arange(3),5)].T
    with np.errstate(divide='ignore'):
        ker=dict(base,weight=base['weight']*scale[None,None,:,:],
                 logweight=base['logweight']+np.log(scale)[None,None,:,:])
    ker['source_mu'] = np.ones_like(ker['weight'],complex)
    ker['source_var'] = np.zeros_like(ker['weight'])
    if scenario >= 2 or band == 0:
        return ker
    angular = source['angular_J_sr'][scenario,:,band] if angular is None else angular
    variance = source['angular_variance_J_sr'][scenario,:,band] if variance is None else variance
    center = source['center_mean'][scenario,:,band] if center is None else center
    for d,z in enumerate(t.DEPTHS):
        for path in range(15):
            b, echo = divmod(path,5)
            zi = z if echo == 0 else 2*((echo+1)//2)*t.H+(-1 if echo%2 else 1)*z
            delta = np.column_stack([t.L/2-xy,np.full(len(xy),zi)])
            n = delta/np.linalg.norm(delta,axis=1)[:,None]
            old = 3*n[:,2]**2 if b == 2 else 1.5*(1-n[:,2]**2)
            for species in range(2):
                aj = angle_lookup(n,angular[species,:,b],source['directions'])
                vj = angle_lookup(n,variance[species,:,b],source['directions'])
                phase = angle_lookup(n,center[species,:,b],source['directions'])
                new = 4*np.pi*aj/energy[species,b] if energy[species,b]>0 else np.zeros(len(xy))
                ratio = new/np.maximum(old,1e-30)
                ker['weight'][:,d,path,species] *= ratio
                with np.errstate(divide='ignore'):
                    ker['logweight'][:,d,path,species] += np.log(ratio)
                frac = np.clip(np.divide(vj,aj,out=np.zeros_like(vj),where=aj>0),0,1)
                # Receiver uses +omega*t-k*r, conjugate of source transform.
                ker['source_mu'][:,d,path,species] = np.sqrt(1-frac)*source_phase(phase)
                ker['source_var'][:,d,path,species] = frac
    return ker


def source_moments(ker, prob):
    mean,cov = r.moments(ker,prob,None)
    if not np.any(ker['source_var']):return mean,cov
    marked = 2*np.sqrt(ker['source_var'])*ker['source_mu'].real
    branch = np.repeat(np.arange(3),5)
    same = branch[:,None] == branch[None,:]
    for i in range(4):
        for j in range(i+1):
            dt=ker['arrival'][:,i,:,None]-ker['arrival'][:,j,None,:]
            product=np.einsum('nps,nqs,s->npq',
                ker['weight'][:,i]*marked[:,i],
                ker['weight'][:,j]*marked[:,j],t.FRACTION)
            extra=t.RATE*np.einsum('n,npq,npq,pq->',prob,product,
                np.exp(-dt*dt/(4*ker['sigma']**2)),same)*np.sqrt(np.pi)*ker['sigma']
            cov[i,j]+=extra
            if i!=j:cov[j,i]+=extra
    return mean,cov


def marked_measure(ker,prob,depth):
    if not np.any(ker['source_var']):return r.packet_jump_measure(ker,prob,depth)
    measures=[]
    for code in range(8):
        signs=np.array([1. if code&(1<<b) else -1. for b in range(3)])[np.repeat(np.arange(3),5)]
        coefficient=ker['source_mu']+np.sqrt(ker['source_var'])*signs[None,None,:,None]
        local=dict(ker,weight=ker['weight']*abs(coefficient)**2)
        q=r.packet_jump_measure(local,prob,depth)
        if q is not None:measures.append(q)
    if not measures:return None
    if any(not np.array_equal(q['age'],measures[0]['age']) for q in measures):
        raise ValueError('Marked packet age support changed')
    return dict(rate=np.concatenate([q['rate']/8 for q in measures],axis=1),
        level=np.concatenate([q['level'] for q in measures],axis=1),
        hot_age=np.concatenate([q['hot_age']/8 for q in measures],axis=0),
        age=measures[0]['age'],max_jump_ratio=max(q['max_jump_ratio'] for q in measures))


def packet_record(ker, catalog, xy, record, band, scenario):
    """One physical parent, shared signed source mark, all its retained paths.
    Two-point centered source quadrature preserves each directional raw second
    moment. It is a conditional marked-parent closure, not a saved collision.
    """
    time = catalog['time_s']; tau = ker['sigma']; delay = 4*tau
    nt = len(time); offsets = catalog['offset_m']
    ids = catalog['parent_source_index']; times = catalog['parent_time_s']
    species = catalog['parent_species']; pid = catalog['parent_id']
    field = np.zeros((nt,4,9,3),complex); shuffled = np.zeros_like(field)
    sumw = np.zeros((nt,4,9)); sumw2 = np.zeros_like(sumw); counts = np.zeros_like(sumw,int)
    p = t.Portable()
    u = p.uniform(record,f'pass5-source-marks-{band}',0,len(pid)*4).reshape(-1,4)
    signs = np.where(u[:,:3]<.5,-1.,1.)
    coefficients = np.empty((len(pid),4,15),complex)
    for e,(ii,ts,ss) in enumerate(zip(ids,times,species)):
        own = np.zeros_like(field); support = np.zeros((nt,4,9),bool)
        for d,z in enumerate(t.DEPTHS):
            for path in range(15):
                b,echo = divmod(path,5)
                mu = ker['source_mu'][ii,d,path,ss]
                vf = ker['source_var'][ii,d,path,ss]
                coefficient = mu + np.sqrt(vf)*signs[e,b]
                coefficients[e,d,path] = coefficient
                amp = ker['weight'][ii,d,path,ss]
                if amp<=0: continue
                zi = z if echo==0 else 2*((echo+1)//2)*t.H+(-1 if echo%2 else 1)*z
                dist = np.sqrt(np.sum((t.L/2-xy[ii])**2)+zi*zi)
                shifted = np.sqrt(np.sum((t.L/2-xy[ii]+offsets)**2,axis=1)+zi*zi)
                ds = shifted-dist
                arrival = ts+ker['arrival'][ii,d,path]+delay
                arg = (time[:,None]-arrival-ds[None,:]/ker['vg'][d,b])/tau
                env = np.where(abs(arg)<=4,np.exp(-.25*arg*arg),0.)
                phase = 2*np.pi*ker['f']*(time-ts-delay)-ker['phase'][ii,d,path]+echo*np.pi
                wave = np.sqrt(amp)*coefficient*env*np.exp(1j*(phase[:,None]-ker['k'][d,b]*ds[None,:]))
                own[:,d,:,b] += wave
                support[:,d] |= abs(arg)<=4
        field += own
        shuffled += own*np.exp(2j*np.pi*u[e,3])
        parent_power=np.sum(abs(own)**2,axis=-1)
        sumw += parent_power; sumw2 += parent_power**2; counts += support
    participation = np.divide(sumw**2,sumw2,out=np.zeros_like(sumw),where=sumw2>0)
    conditional = bool(catalog['conditional'])
    elastic = .5*np.sum(abs(field)**2,axis=-1)
    shuffle = .5*np.sum(abs(shuffled)**2,axis=-1)
    return dict(time_s=time,depth_m=t.DEPTHS,offset_m=offsets,
                parent_id=pid,parent_species=species,parent_source_index=ids,parent_time_s=times,
                parent_source_marks=signs,source_coefficients=coefficients,
                signed_branch_amplitude=field,elastic_J_m3=elastic,shuffled_J_m3=shuffle,
                effective_count=participation[:,:,4],counts=counts[:,:,4],
                spatial_effective_count=participation,spatial_counts=counts,
                conditional=np.array(conditional),expected_parents=catalog['expected_parents'],
                unconditional_elastic_J_m3=np.zeros_like(elastic) if conditional else elastic)


def selected_bins(m, source):
    """Same eight frozen frequency nodes split once, no range-refinement ladder."""
    path=HERE/'selected_bins.npz'
    if path.exists(): return dict(np.load(path))
    directions=source['directions'];dw=source['angular_weight_sr']
    angular=np.zeros((2,2,2,2,len(directions),3))
    variance=np.zeros_like(angular);freq=np.empty((2,2))
    tracks=s.tracks(m)
    for bi,band in enumerate((4,5)):
        x,w=np.polynomial.legendre.leggauss(8)
        fs=np.exp(np.log(t.BANDS[band])+(x+1)/2*np.log(10.))
        dfs=w/2*np.log(10.)*fs
        for bn in range(2):freq[bi,bn]=np.exp(np.mean(np.log(fs[bn*4:(bn+1)*4])))
        for scenario,scale in enumerate((1.,.25)):
            for sp,tr in enumerate(tracks):
                mark=s.marks(tr,scale)
                for j,(f,df) in enumerate(zip(fs,dfs)):
                    ap,vv,at,k,vg=s.amplitudes(f,directions,tr,mark,source['heat_scale'][scenario,sp])
                    dos=k*k/(2*t.RHO*(2*np.pi)**3*vg)*(2*np.pi*df)
                    angular[scenario,sp,bi,j//4] += dos*(abs(ap+at)**2+vv)
                    variance[scenario,sp,bi,j//4] += dos*vv
    energy=np.einsum('...nb,n->...b',angular,dw)
    np.savez_compressed(path,angular=angular,variance=variance,energy=energy,frequency=freq)
    return dict(angular=angular,variance=variance,energy=energy,frequency=freq)


def receivers(m, ref):
    source=s.load(); energy=source['energy_J']
    xy,prob,area=r.sources(m);nd=r.population();p=t.Portable()
    s.table('source_quadrature.csv',[dict(x_m=x,y_m=y,area_m2=a,even_probability=prob[0,i],hot_probability=prob[1,i]) for i,((x,y),a) in enumerate(zip(xy,area))])
    # Exact accepted low histories/tables are copied, not regenerated.
    old=t.HERE;t.HERE=HERE
    try:r.slow_background(ref)
    finally:t.HERE=old
    allfront=readrows('inherited_front_ancestry.csv')
    front={(x['beam'],int(x['band']),float(x['depth_m'])):float(x['unit_intensity_mean_distinct']) for x in allfront}
    stats=[];peaks=[];distribution={};refined=[];ledger=[]
    bins=selected_bins(m,source)
    for band in range(7):
        for law in LAWS:
            geometry=r.kernel(xy,band,np.ones((2,3)),law)
            kernels=[make_kernel(xy,band,energy[sc,:,band],law,sc,source,geometry=geometry) for sc in range(5)]
            for beam in BEAMS:
                pp=prob[int(beam!='even_steady')]
                base=kernels[0]
                lam=r.count_masks(base,pp)
                core=np.array([sum(lam[j] for j in range(1,16) if j&(1<<d)) for d in range(4)])
                catalogs=[]
                if np.any(core<200):
                    selected=dict(base,weight=base['weight'].copy());selected['weight'][:,core>=200]=0
                    for record in range(4):
                        catalog,_=r.sparse_record(p,selected,xy,pp,beam,record,band,'track',law,nd,0.)
                        catalogs.append(catalog)
                for sc,name in enumerate(s.SHAPES):
                    ker=kernels[sc];mean,cov=source_moments(ker,pp)
                    eff=np.divide(mean**2,np.diag(cov),out=np.zeros(4),where=np.diag(cov)>0)
                    factors=r.retarded_factor(ker,pp,mean,beam)
                    unc=dict(ker,weight=np.exp(ker['logweight']))
                    fullmean,_=r.moments(unc,pp,None)
                    empirical=[[] for _ in range(4)]
                    for record,cat in enumerate(catalogs):
                        if cat is None:continue
                        selected=dict(ker,weight=ker['weight'].copy())
                        selected['weight'][:,core>=200]=0
                        payload=packet_record(selected,cat,xy,record,band,sc)
                        np.savez_compressed(HERE/f'{name}_b{band}_{law}_{beam}_r{record}.npz',**payload)
                        e=payload['unconditional_elastic_J_m3']
                        for d in range(4):
                            if core[d]>=200:continue
                            empirical[d].append(e[:,d,4])
                            raw=payload['elastic_J_m3'][:,d]
                            ti,point=np.unravel_index(raw.argmax(),raw.shape)
                            own=raw[:,point]
                            peaks.append(dict(scenario=name,band=band,loss=law,beam=beam,record=record,
                                depth_m=t.DEPTHS[d],conditional=bool(payload['conditional']),
                                sample_time_s=float(payload['time_s'][ti]),offset_x_m=float(payload['offset_m'][point,0]),offset_y_m=float(payload['offset_m'][point,1]),
                                peak_J_m3=float(own[ti]),own_mean_J_m3=float(own.mean()),
                                tallest_over_mean=float(own[ti]/own.mean()) if own.mean()>0 else 0.,
                                peak_eV_D=float(own[ti]/(t.QE*nd[d])),
                                distinct_core_at_peak=int(payload['spatial_counts'][ti,d,point]),N_eff_at_peak=float(payload['spatial_effective_count'][ti,d,point]),
                                shuffled_same_peak_J_m3=float(payload['shuffled_J_m3'][ti,d,point]),
                                exposure_s=float(payload['time_s'][-1]-payload['time_s'][0]+payload['time_s'][1]-payload['time_s'][0]),
                                physical_catalog_sha256=hashlib.sha256(payload['parent_id'].tobytes()+payload['parent_time_s'].tobytes()).hexdigest()))
                    for d,z in enumerate(t.DEPTHS):
                        if core[d]<200:
                            values=np.concatenate(empirical[d]) if empirical[d] else np.zeros(2048)
                            model='shared signed marked-parent packets; unconditional zeros retained'
                        elif eff[d]<100:
                            measure=marked_measure(ker,prob,d)
                            values,_=r.draw_packet_energy(measure,p,beam,band,d)
                            model='positive compound-Poisson parent-packet marginal'
                        else:
                            normals=np.array([r.draw_normal(p,j,f'pass5-paired-{band}-{law}',4) for j in range(128)])
                            val,vec=np.linalg.eigh(cov);root=vec*np.sqrt(np.maximum(val,0))
                            values=.5*(mean[d]*factors[:,d]+(normals@root.T)[:,d]*np.sqrt(factors[:,d]))
                            if np.any(values<0):raise ValueError('Invalid central energy approximation')
                            model='large-participation packet-energy central projection'
                        distribution[f'{name}_{band}_{law}_{beam}_{d}']=values
                        stats.append(dict(scenario=name,band=band,lo_Hz=t.BANDS[band],hi_Hz=t.BANDS[band+1],loss=law,beam=beam,depth_m=z,
                            actual_front_ancestry=front[beam,band,z],energetic_core_parents=core[d],N_eff=eff[d],
                            mean_elastic_J_m3=.5*mean[d],uncensored_mean_J_m3=.5*fullmean[d],
                            mean_eV_D=.5*fullmean[d]/(t.QE*nd[d]),sample_max_J_m3=float(values.max()),
                            tallest_over_sample_mean=float(values.max()/values.mean()) if values.mean()>0 else 0.,
                            samples=len(values),zeros=int(np.count_nonzero(values==0)),
                            exposure_s=4*64/ker['f'] if core[d]<200 else 128*15*50e-9,
                            distribution=model,tail_underflow=bool(fullmean[d]==0),
                            depth_reach_cutoff_J_m3=1e-20,above_reach_cutoff=bool(.5*fullmean[d]>=1e-20)))
                    for sp in range(2):
                        power=t.RATE*t.FRACTION[sp]*energy[sc,sp,band].sum()
                        gamma=float(t.loss(ker['f'],np.mean(ker['vg']),law));edge=.1*np.mean(ker['vg'])/(t.L/2)
                        ledger.append(dict(scenario=name,band=band,loss=law,beam=beam,species=sp,launch_W=power,
                            inventory_J=power/(gamma+edge),heat_W=power*gamma/(gamma+edge),escape_W=power*edge/(gamma+edge),balance_W=0.))
            # Two bins are needed only where low-edge survival changes the depth decision.
            if band in (4,5):
                bi=band-4
                for sc in range(2):
                    ks=[]
                    for bn in range(2):
                        ee=bins['energy'][sc,:,bi,bn]
                        ks.append(make_kernel(xy,band,ee,law,sc,source,frequency=bins['frequency'][bi,bn],
                            bandwidth=(t.BANDS[band+1]-t.BANDS[band])/2,
                            angular=bins['angular'][sc,:,bi,bn],variance=bins['variance'][sc,:,bi,bn]))
                    for beam in BEAMS:
                        pp=prob[int(beam!='even_steady')]
                        total=np.zeros(4)
                        for ker in ks:
                            um,_=r.moments(dict(ker,weight=np.exp(ker['logweight'])),pp,None)
                            total+=.5*um
                        for d,z in enumerate(t.DEPTHS):
                            refined.append(dict(scenario=s.SHAPES[sc],band=band,loss=law,beam=beam,depth_m=z,
                                uncensored_mean_J_m3=total[d],mean_eV_D=total[d]/(t.QE*nd[d]),
                                cutoff_J_m3=1e-20,above_cutoff=bool(total[d]>=1e-20),
                                frequency_bins=2,source_nodes=8))
            s.table('readout.csv',stats);s.table('same_record_peaks.csv',peaks)
            s.table('frequency_refinement.csv',refined);s.table('transport_ledger.csv',ledger)
            print(json.dumps(dict(receiver_band=band,loss=law,rows=len(stats))),flush=True)
    np.savez_compressed(HERE/'distributions.npz',**distribution)
    s.emit('receiver_model.json',dict(whole_foil='accepted 24-mm equivalent square, 250-um thick, fixed 8x8 grain map',
        energy_metric='elastic cycle-averaged envelope .5 sum |A_branch|^2, total packet energy twice elastic',
        angular='source frequency-integrated directional work, four-neighbour positive interpolation; common band-centre phase',
        source_marks='one independent real centered two-point source projection per parent/branch, variance from compound-Poisson recoil raw moment; shared across paths',
        observation='same physical parent catalog and source marks across all five source comparisons; 4x512 common packet samples per sparse cell',
        support='same first five images/eight-sigma fronts; 12-energy-length core; inherited finite ancestry separate from participation',
        loss='Cu per actual path length, Al per transit time; no transferable measured Pd GHz row in inspected sources',
        reach='mean elastic inventory >=1e-20 J/m3; four sampled depths, with two-bin correction in bands4/5',
        remainder='mode-diagonal conditional envelope, finite modeled paths, coarse source angular grid; unresolved weak-front phases and late cavity/grain conversion remain open'))

def replay(m):
    """Replay frozen physical catalogs without rerunning source or dense draws."""
    source=s.load();xy,prob,_=r.sources(m);nd=r.population()
    rows=readrows('readout.csv')
    peaks=[x for x in readrows('same_record_peaks.csv') if x['scenario'] not in s.SCENARIOS]
    with np.load(HERE/'distributions.npz') as d:distributions={k:d[k] for k in d.files}
    for band in range(7):
        for law in LAWS:
            available=list(HERE.glob(f'track_b{band}_{law}_*_r[0-3].npz'))
            if not available:continue
            geometry=r.kernel(xy,band,np.ones((2,3)),law)
            for sc,name in enumerate(s.SCENARIOS):
                ker=make_kernel(xy,band,source['energy_J'][sc,:,band],law,sc,source,geometry=geometry)
                for beam in BEAMS:
                    cells=[x for x in rows if x['scenario']==name and int(x['band'])==band and x['loss']==law and x['beam']==beam]
                    sparse=np.array([float(next(x for x in cells if float(x['depth_m'])==z)['energetic_core_parents'])<200 for z in t.DEPTHS])
                    selected=dict(ker,weight=ker['weight'].copy())
                    selected['weight'][:,~sparse]=0
                    empirical=[[] for _ in range(4)]
                    for path in sorted(HERE.glob(f'{name}_b{band}_{law}_{beam}_r[0-3].npz')):
                        record=int(path.stem[-1])
                        with np.load(path) as d:catalog={k:d[k] for k in d.files}
                        payload=packet_record(selected,catalog,xy,record,band,sc)
                        np.savez_compressed(path,**payload)
                        for dep,z in enumerate(t.DEPTHS):
                            if not sparse[dep]:continue
                            empirical[dep].append(payload['unconditional_elastic_J_m3'][:,dep,4])
                            raw=payload['elastic_J_m3'][:,dep]
                            ti,point=np.unravel_index(raw.argmax(),raw.shape);own=raw[:,point]
                            peaks.append(dict(scenario=name,band=band,loss=law,beam=beam,record=record,
                                depth_m=z,conditional=bool(payload['conditional']),
                                sample_time_s=float(payload['time_s'][ti]),offset_x_m=float(payload['offset_m'][point,0]),offset_y_m=float(payload['offset_m'][point,1]),
                                peak_J_m3=float(own[ti]),own_mean_J_m3=float(own.mean()),
                                tallest_over_mean=float(own[ti]/own.mean()) if own.mean()>0 else 0.,
                                peak_eV_D=float(own[ti]/(t.QE*nd[dep])),
                                distinct_core_at_peak=int(payload['spatial_counts'][ti,dep,point]),N_eff_at_peak=float(payload['spatial_effective_count'][ti,dep,point]),
                                shuffled_same_peak_J_m3=float(payload['shuffled_J_m3'][ti,dep,point]),
                                exposure_s=float(payload['time_s'][-1]-payload['time_s'][0]+payload['time_s'][1]-payload['time_s'][0]),
                                physical_catalog_sha256=hashlib.sha256(payload['parent_id'].tobytes()+payload['parent_time_s'].tobytes()).hexdigest()))
                    for dep,z in enumerate(t.DEPTHS):
                        if not sparse[dep]:continue
                        values=np.concatenate(empirical[dep]) if empirical[dep] else np.zeros(2048)
                        distributions[f'{name}_{band}_{law}_{beam}_{dep}']=values
                        cell=next(x for x in cells if float(x['depth_m'])==z)
                        cell.update(sample_max_J_m3=float(values.max()),
                            tallest_over_sample_mean=float(values.max()/values.mean()) if values.mean()>0 else 0.,
                            zeros=int(np.count_nonzero(values==0)),samples=len(values))
            print(json.dumps(dict(frozen_catalog_replay_band=band,loss=law)),flush=True)
    s.table('readout.csv',rows);s.table('same_record_peaks.csv',peaks)
    np.savez_compressed(HERE/'distributions.npz',**distributions)
    model=json.loads((HERE/'receiver_model.json').read_text())
    model['phase_convention']='source exp(+iomega*t-ik.r), receiver conjugate carrier; complex-conjugate source phasor'
    s.emit('receiver_model.json',model)



def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['launch','receivers','replay','refine','feedback','check','readout','package','upload'])
    p.add_argument('--reference-dir',default=str(ROOT/'g13_rogue/pass3/rung2'))
    p.add_argument('--previous-evidence-dir',help='Verified accepted Pass4 evidence extraction for unchanged actual-packet comparison')
    a=p.parse_args()
    if a.action in ('package','upload'):
        import delivery
        getattr(delivery,a.action)();return
    m=s.reference(a.reference_dir)
    if a.action=='launch':s.spectrum(m)
    elif a.action=='receivers':receivers(m,a.reference_dir)
    elif a.action=='replay':replay(m)
    elif a.action=='refine':
        import refine
        refine.run(m,a.previous_evidence_dir)
    elif a.action=='feedback':
        import feedback
        feedback.run(m)
    elif a.action=='check':
        import readout
        readout.check(m,a.reference_dir,a.previous_evidence_dir)
    else:
        import readout
        readout.write_pages(m)

if __name__=='__main__':main()
