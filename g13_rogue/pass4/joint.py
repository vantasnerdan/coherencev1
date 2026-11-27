"""One Poisson parent measure shared across all seven bands and four depths.
Receiver kernels are reused; no new launch or frequency exponents are computed.
"""
import numpy as np
import transport as t
import receiver as r


def census(kernels,xy,prob,p,law,full_support=False):
    frozen=np.load(t.ROOT/'g13_rogue/pass2/replication/inputs.npz')
    timing=frozen['time_weights'][2,:,0];dt=50e-9
    sample_bins=4704+np.arange(128)*15
    arrays={'time_s':np.arange(128)*15*dt,'depth_m':t.DEPTHS,'band_edges_Hz':t.BANDS}
    report=[];hist={}
    for i in range(len(xy)):
        endpoints=[]
        for b,k in enumerate(kernels):
            sigma=k['sigma']
            for d in range(4):
                supported=np.isfinite(k['logweight'][i,d]).any(axis=-1) if full_support else k['weight'][i,d].sum(axis=-1)>0
                for a,z in r.unions([(age,age+8*sigma) for age in k['arrival'][i,d,supported]]):
                    flag=1<<(4*b+d);endpoints.extend([(a,flag),(z,-flag)])
        state=0;before=0.
        for age,delta in sorted(endpoints):
            if state and age>before:
                lo=int(np.floor(before/dt));hi=int(np.floor(age/dt))
                for binid in range(lo,hi+1):
                    width=max(0.,min(age,(binid+1)*dt)-max(before,binid*dt))
                    if width:
                        key=(state,binid)
                        if key not in hist:hist[key]=np.zeros(2)
                        hist[key]+=prob[:,i]*width
            state+=delta;before=age
    membership=sorted({k[0] for k in hist});index={v:j for j,v in enumerate(membership)}
    even=np.zeros(len(membership));hot=np.zeros(len(membership));burst=np.zeros((128,len(membership)))
    for (mask,agebin),weights in hist.items():
        j=index[mask];even[j]+=t.RATE*weights[0];hot[j]+=t.RATE*weights[1]
        emit=sample_bins-agebin
        if emit.min()<0 or emit.max()>=len(timing):raise ValueError('Joint source history outside accepted finite record')
        burst[:,j]+=t.RATE*weights[1]*timing[emit]
    for beam,lam in (('even_steady',np.tile(even,(128,1))),('hot_steady',np.tile(hot,(128,1))),('hot_bursts',burst)):
        counts=np.zeros((128,7,4),np.int64)
        for rr in range(128):
            for j,mask in enumerate(membership):
                n=p.poisson(float(lam[rr,j]),rr,500000+j)
                if n:
                    for bit in range(28):
                        if mask&(1<<bit):counts[rr,bit//4,bit%4]+=n
        arrays[law+'_'+beam]=counts;arrays[law+'_'+beam+'_mask_intensity']=lam
        unit=even if beam=='even_steady' else hot
        boundfactor=float(frozen['time_weights'][2].max()) if beam=='hot_bursts' else 1.
        for b in range(7):
            for d,z in enumerate(t.DEPTHS):
                values=counts[:,b,d];bit=1<<(4*b+d)
                nominal=sum(unit[j] for j,mask in enumerate(membership) if mask&bit)
                report.append(dict(loss=law,beam=beam,band=b,depth_m=z,mean_distinct=float(values.mean()),min_distinct=int(values.min()),median_distinct=float(np.median(values)),max_distinct=int(values.max()),zero_snapshots=int(np.count_nonzero(values==0)),snapshots=128,time_step_s=15*dt,unit_intensity_mean_distinct=float(nominal),Cox_max_mean_distinct_bound=float(nominal*boundfactor),meaning='uncensored modeled-front ancestry; finite first-five-image lower bound, not detectable phonon count' if full_support else '12-energy-length importance census; not all weak parent tails'))
    arrays[law+'_membership_mask']=np.array(membership,np.int64)
    print('joint',law,'membership subsets',len(membership),flush=True)
    return arrays,report


def run(m,full_support=False):
    xy,prob,area=r.sources(m);p=t.Portable();energy=np.load(t.HERE/'launch.npz')['built'][0]
    arrays={};report=[]
    for law in (('full_support',) if full_support else ('Cu_extrapolated','Al_transfer')):
        kernels=[r.kernel(xy,b,energy[:,b],'Cu_extrapolated' if full_support else law) for b in range(7)]
        data,rows=census(kernels,xy,prob,p,law,full_support)
        arrays.update(data);report.extend(rows)
    stem='band_depth_full_support' if full_support else 'band_depth_joint'
    np.savez_compressed(t.HERE/(stem+'.npz'),**arrays);t.table(stem+'.csv',report)
    t.dump(stem+'_model.json',dict(parent_measure='one finite Poisson source/emission-time measure across28 band/depth cells',burst='retarded accepted50-ns Cox intensity in post-burn observation',causal_support='modeled front at group transit; packet centre four sigma later; finite eight-sigma observation envelope',observation_span_s=128*15*50e-9,modeled_max_prehistory_s=max(float(k['arrival'].max()+8*k['sigma']) for k in kernels),support='all nonzero modeled rays, irrespective of attenuation; lower bound if later cavity paths included' if full_support else '12 conditional energy lengths; importance census only',same_parent='D2+ and all echoes one parent',limitations='not an infinite-history count, detector occupancy, coherence or experimentally observed rate'))

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference-dir',default=str(t.ROOT/'g13_rogue/pass3/rung2'));p.add_argument('--full-support',action='store_true');a=p.parse_args();run(t.reference(a.reference_dir),a.full_support)
