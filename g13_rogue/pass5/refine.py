"""The accepted two spectral filters, with updated source and shared tail catalogs.
No new frequency/range quadrature: consumes Pass4 frozen spectrum and Pass5 bins.
"""
from pathlib import Path
import hashlib
import numpy as np
import source as s
from run import make_kernel,packet_record,source_moments,readrows

t,r=s.t,s.r
HERE,ROOT=s.HERE,s.ROOT


def counts(kernels,prob,full=False):
    out=np.zeros(4)
    for i,p in enumerate(prob):
        for d in range(4):
            intervals=[]
            for ker in kernels:
                take=np.isfinite(ker['logweight'][i,d]).any(-1) if full else ker['weight'][i,d].sum(-1)>0
                intervals.extend((a,a+8*ker['sigma']) for a in ker['arrival'][i,d,take])
            out[d]+=t.RATE*p*sum(b-a for a,b in r.unions(intervals))
    return out


def combined(kernels,prob,full=False):
    scale=np.zeros(4)
    if full:
        scale=np.max(np.stack([k['logweight'] for k in kernels]),axis=(0,1,3,4))
        kernels=[dict(k,weight=np.exp(k['logweight']-scale[None,:,None,None])) for k in kernels]
    moments=[source_moments(k,prob) for k in kernels]
    mean=sum(x[0] for x in moments);variance=sum(np.diag(x[1]) for x in moments)
    a,b=kernels;sa,sb=a['sigma'],b['sigma'];den=sa*sa+sb*sb
    same=np.repeat(np.arange(3),5)[:,None]==np.repeat(np.arange(3),5)[None,:]
    for d in range(4):
        wa,wb=a['weight'][:,d],b['weight'][:,d]
        ma=2*np.sqrt(a['source_var'][:,d])*a['source_mu'][:,d].real
        mb=2*np.sqrt(b['source_var'][:,d])*b['source_mu'][:,d].real
        product=np.einsum('nps,nqs,s->npq',wa,wb,t.FRACTION)
        product+=np.einsum('nps,nqs,s->npq',wa*ma,wb*mb,t.FRACTION)*same
        dt=a['arrival'][:,d,:,None]+4*sa-b['arrival'][:,d,None,:]-4*sb
        covariance=t.RATE*np.einsum('n,npq,npq->',prob,product,np.exp(-dt*dt/(2*den)))*np.sqrt(2*np.pi)*sa*sb/np.sqrt(den)
        variance[d]+=2*covariance
    neff=np.divide(mean*mean,variance,out=np.zeros(4),where=variance>0)
    with np.errstate(divide='ignore'):
        logarithm=np.log(.5*mean)+scale
    return .5*mean*np.exp(scale),neff,logarithm


def assumed_bin(m,lo,hi):
    scale=np.hypot(.5e-12,20e-9/4460.)
    aa=(2*np.pi*lo*scale)**2;bb=(2*np.pi*hi*scale)**2
    frac=r.gammainc(1.5,bb)-r.gammainc(1.5,aa) if bb<1 else r.gammaincc(1.5,aa)-r.gammaincc(1.5,bb)
    norm=r.gammaincc(1.5,(2*np.pi*2.5e5*scale)**2)
    return np.repeat((m['heatenergy'][:,1]*.99*.1*frac/norm/3)[:,None],3,axis=1)


def run(m,previous=None):
    source=s.load();selected=dict(np.load(HERE/'selected_bins.npz'))
    old=dict(np.load(ROOT/'g13_rogue/pass4/selected_spectrum.npz'))
    xy,prob,_=r.sources(m);nd=r.population();tracks=s.tracks(m)
    readout=readrows('readout.csv');report=[];tail=[];matched=[];phase_data={}
    for band in (4,5):
        bi=band-4;f=old[f'band{band}_frequency_Hz'];terms=old[f'band{band}_terms_J']
        split=np.sqrt(t.BANDS[band]*t.BANDS[band+1]);edges=(t.BANDS[band],split,t.BANDS[band+1])
        frequencies=[];energies=np.zeros((5,2,2,3));centers=np.empty((2,2,2,len(source['directions']),3),complex)
        for bn in range(2):
            take=(f<split) if bn==0 else (f>=split)
            w=np.einsum('fsb,s->f',terms[take,0].sum(-1),t.FRACTION)
            fc=float(f[take]@w/w.sum());frequencies.append(fc)
            energies[:2,:,bn]=selected['energy'][:,:,bi,bn]
            for shape in range(2):energies[shape+2,:,bn]=terms[take,shape].sum(axis=(0,-1))
            energies[4,:,bn]=assumed_bin(m,edges[bn],edges[bn+1])
            # Evaluate the existing angular source amplitude at the frozen carrier;
            # this does not change source work, frequency nodes, or range model.
            for sc,scale in enumerate((1.,.25)):
                for sp,tr in enumerate(tracks):
                    ap,var,at,_,_=s.amplitudes(fc,source['directions'],tr,s.marks(tr,scale),source['heat_scale'][sc,sp])
                    centers[sc,sp,bn]=ap+at
        phase_data[f'band{band}_frequency_Hz']=np.array(frequencies)
        phase_data[f'band{band}_source_mean']=centers
        for law in ('Cu_extrapolated','Al_transfer'):
            geometry=[r.kernel(xy,band,np.ones((2,3)),law,frequency=frequencies[bn],bandwidth=edges[bn+1]-edges[bn]) for bn in range(2)]
            bank=[]
            for sc in range(5):
                kk=[]
                for bn in range(2):
                    kwargs={} if sc>=2 else dict(angular=selected['angular'][sc,:,bi,bn],variance=selected['variance'][sc,:,bi,bn],center=centers[sc,:,bn])
                    kk.append(make_kernel(xy,band,energies[sc,:,bn],law,sc,source,geometry=geometry[bn],**kwargs))
                bank.append(kk)
            for beam in ('even_steady','hot_steady','hot_bursts'):
                pp=prob[int(beam!='even_steady')]
                core=counts(bank[0],pp);front=counts(bank[0],pp,True)
                for sc,name in enumerate(s.SHAPES):
                    coremean,neff,_=combined(bank[sc],pp);fullmean,fullneff,logmean=combined(bank[sc],pp,True)
                    for d,z in enumerate(t.DEPTHS):
                        report.append(dict(scenario=name,band=band,loss=law,beam=beam,depth_m=z,
                            uncensored_mean_J_m3=fullmean[d],mean_eV_D=fullmean[d]/(t.QE*nd[d]),
                            core_mean_J_m3=coremean[d],energetic_core_parents=core[d],actual_front_ancestry=front[d],
                            N_eff=neff[d],uncensored_N_eff=fullneff[d],cutoff_J_m3=1e-20,above_cutoff=bool(fullmean[d]>=1e-20),
                            ln_mean_over_1_J_m3=logmean[d],numeric_underflow=bool(fullmean[d]==0 and np.isfinite(logmean[d])),
                            low_carrier_Hz=frequencies[0],high_carrier_Hz=frequencies[1],
                            frequency_bins=2,source_nodes=8 if sc<2 else 12,
                            count_scope='same-parent union, finite first-five-image/eight-sigma filters; unit-intensity expectation'))
                        row=next(x for x in readout if x['scenario']==name and int(x['band'])==band and x['loss']==law and x['beam']==beam and float(x['depth_m'])==z)
                        row.update(original_bandcenter_uncensored_mean_J_m3=row['uncensored_mean_J_m3'],
                            original_bandcenter_core_parents=row['energetic_core_parents'],
                            uncensored_mean_J_m3=fullmean[d],mean_elastic_J_m3=coremean[d],mean_eV_D=fullmean[d]/(t.QE*nd[d]),
                            energetic_core_parents=core[d],actual_front_ancestry=front[d],N_eff=neff[d],uncensored_N_eff=fullneff[d],
                            ln_uncensored_mean_over_1_J_m3=logmean[d],
                            sample_frequency_scope='band-centre signed records; two-filter mean/count columns',
                            tail_underflow=bool(fullmean[d]==0),above_reach_cutoff=bool(fullmean[d]>=1e-20))
                if band!=5:continue
                # The upper filter is beyond the 12-length core at all deep
                # receivers; retain its uncensored mean above, not ghost packets.
                assert np.all(bank[0][1]['weight'][:,1:]==0)
                for record in range(4):
                    filename=f'tail_b5_{law}_{beam}_r{record}.npz'
                    preserved=HERE/('old_'+filename)
                    if not preserved.exists():
                        if previous is None:raise ValueError('Supply --previous-evidence-dir for accepted tail catalog bytes')
                        original=Path(previous)/filename
                        preserved.write_bytes(original.read_bytes())
                    with np.load(preserved) as d:catalog={k:d[k] for k in d.files}
                    payloads=[]
                    for sc,name in enumerate(s.SHAPES):
                        ker=dict(bank[sc][0],weight=bank[sc][0]['weight'].copy())
                        ker['weight'][:,0]=0
                        payload=packet_record(ker,catalog,xy,record,5,sc);payloads.append(payload)
                        payload['carrier_Hz']=np.array(frequencies[0]);payload['subband_edges_Hz']=np.array(edges[:2])
                        np.savez_compressed(HERE/f'tail_{name}_b5_{law}_{beam}_r{record}.npz',**payload)
                    for dep in range(1,4):
                        baseline=payloads[2]['elastic_J_m3'][:,dep]
                        ti,point=np.unravel_index(baseline.argmax(),baseline.shape)
                        for sc,name in enumerate(s.SHAPES):
                            payload=payloads[sc];y=payload['elastic_J_m3'][:,dep]
                            pi,ppoint=np.unravel_index(y.argmax(),y.shape);own=y[:,ppoint]
                            tail.append(dict(scenario=name,band=5,loss=law,beam=beam,record=record,depth_m=t.DEPTHS[dep],
                                conditional=bool(payload['conditional']),peak_J_m3=float(y[pi,ppoint]),own_mean_J_m3=float(own.mean()),
                                tallest_over_own_mean=float(y[pi,ppoint]/own.mean()) if own.mean()>0 else 0.,peak_eV_D=float(y[pi,ppoint]/(t.QE*nd[dep])),
                                peak_time_s=float(payload['time_s'][pi]),peak_offset_x_m=float(payload['offset_m'][ppoint,0]),peak_offset_y_m=float(payload['offset_m'][ppoint,1]),
                                old_anchor_time_s=float(payload['time_s'][ti]),old_anchor_offset_x_m=float(payload['offset_m'][point,0]),old_anchor_offset_y_m=float(payload['offset_m'][point,1]),
                                actual_core_parents=int(payload['spatial_counts'][pi,dep,ppoint]),N_eff_at_peak=float(payload['spatial_effective_count'][pi,dep,ppoint]),
                                shuffled_same_peak_J_m3=float(payload['shuffled_J_m3'][pi,dep,ppoint]),
                                paired_old_peak_J_m3=float(y[ti,point]),old_peak_J_m3=float(baseline[ti,point]),
                                exposure_s=64/frequencies[0],carrier_Hz=frequencies[0],
                                physical_catalog_sha256=hashlib.sha256(payload['parent_id'].tobytes()+payload['parent_time_s'].tobytes()).hexdigest()))
                    for sc,name in enumerate(s.SHAPES):
                        for dep in range(1,4):
                            matched.append((name,law,beam,record,dep,payloads[sc]['unconditional_elastic_J_m3'][:,dep,4]))
            print('matched filters',band,law,flush=True)
    for row in readout:
        if int(row['band'])!=5 or float(row['depth_m'])==t.DEPTHS[0]:continue
        dep=int(np.flatnonzero(t.DEPTHS==float(row['depth_m']))[0])
        values=np.concatenate([x[-1] for x in matched if x[0]==row['scenario'] and x[1]==row['loss'] and x[2]==row['beam'] and x[4]==dep])
        row.update(sample_max_J_m3=float(values.max()),tallest_over_sample_mean=float(values.max()/values.mean()) if values.mean()>0 else 0.,
            zeros=int(np.count_nonzero(values==0)),samples=len(values),exposure_s=4*64/float(next(x['low_carrier_Hz'] for x in report if int(x['band'])==5)),
            sample_frequency_scope='matched accepted low-filter signed core; upper filter outside deep core')
    s.table('frequency_refinement.csv',report);s.table('tail_comparison.csv',tail);s.table('readout.csv',readout)
    np.savez_compressed(HERE/'refined_source_phase.npz',**phase_data)
