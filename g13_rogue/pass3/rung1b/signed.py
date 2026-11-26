"""Signed whole-foil pressure / diffusing eigenstrain / old random source.
Run at repository root. The rung2 spatial module changes only the global space.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,sys,tarfile
import numpy as np
from numpy.polynomial.legendre import leggauss,legval,legder
from scipy.constants import elementary_charge as QE
from scipy.stats import t as student
from scipy.ndimage import label
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'rung1'))
import rungs as old
DT=5e-8;BURN=4704;SAMPLE=2000;NREC=16
CP=25.89/.1064;KTH=71.8;ALPHA=11.77e-6;RHO=12007.
ATH=KTH/(RHO*CP);NZ=64
MD=3.3435837768e-27
MOMENTUM=np.sqrt(2*MD*5000*QE)*np.array([1.,np.sqrt(2)])
BASE='fe514a174214b01558c228fc9f8982ef3791fb06'
BEAMS=('even_steady','hot_steady','hot_bursts','even_100k','even_tuned','even_between')

def dump(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def table(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dest(rung):return HERE if rung<2 else HERE.parent/'rung2'

_SPATIAL=None
def _space():
    global _SPATIAL
    if _SPATIAL is None:
        sys.path.insert(0,str(HERE.parent/'rung2'))
        import spatial
        _SPATIAL=spatial
    return _SPATIAL

def stats(x):
    x=np.asarray(x,float);n=len(x);mu=float(x.mean());sd=float(x.std(ddof=1)) if n>1 else 0.
    half=float(student.ppf(.975,n-1)*sd/np.sqrt(n)) if n>1 and sd else 0.
    return dict(n_records=n,mean=mu,sd=sd,ci95_lo=mu-half,ci95_hi=mu+half)

def maps(c,m,points):
    if m['rung']==2:return _space().maps(np.asarray(points),m['v'],c.DIMS)
    phi,b=c.basis(np.asarray(points),degree=m['degree']);n=phi.shape[1]
    return np.array([phi@m['v'][i*n:(i+1)*n] for i in range(3)]).transpose(1,0,2),b@m['v']

def beam_sites(c,rung):
    sub=3 if rung==2 else 1
    axis=(c.DIMS[0]-.01)/2+(.5+np.arange(6*sub))*.01/(6*sub)
    xy=np.stack(np.meshgrid(axis,axis,indexing='ij'),axis=-1).reshape(-1,2)
    parent=(np.arange(6*sub)//sub)[:,None]*6+np.arange(6*sub)[None,:]//sub
    group=(xy[:,0]>=c.DIMS[0]/2).astype(int)*2+(xy[:,1]>=c.DIMS[1]/2)
    return xy,parent.ravel(),group

def profile(c,xy,hot):
    p=xy-c.DIMS[:2]/2
    if not hot:return np.ones(len(xy))/len(xy)
    w=np.exp(-np.sum(p*p,axis=1)/(2*.003**2))
    for center,gain in (([.0025,-.0025],5.),([-.0025,.0025],3.)):
        w+=gain*np.exp(-np.sum((p-center)**2,axis=1)/(2*.0012**2))
    return w/w.sum()

def beam_weights(c,m,hot,frozen=None):
    weights=profile(c,m['xy'],hot)
    if m['rung']==2:
        if frozen is None:frozen=old.accepted()[3]
        coarse=frozen['site_weights'][2 if hot else 0]
        subtotal=np.bincount(m['parent'],weights,minlength=36)
        weights*=coarse[m['parent']]/subtotal[m['parent']]
    return weights


def heat_profiles(c,frozen,nz=NZ):
    _,stop=c.profiles_and_stopping();s=stop['PdD',2.5];h=c.DIMS[2];ks=np.arange(nz)*np.pi/h
    cols=('electron_eV','host_recoil_eV','D_recoil_eV');cos=[];central=[]
    for species,tr in enumerate(s['tracks']):
        z=np.array([0.]+[row['depth_A']*1e-10 for row in tr['timeline']]);mid=(z[1:]+z[:-1])/2
        cc=[];ee=[]
        for col in cols:
            inc=np.diff([0.]+[row[col] for row in tr['timeline']]);ee.append(inc.sum()*(1 if species==0 else 2)*QE)
            cc.append(np.cos(ks[:,None]*mid)@inc/inc.sum())
        cc.append(np.cos(ks*tr['terminal_depth_A']));ee.append(tr['terminal_unresolved_eV']*(1 if species==0 else 2)*QE)
        cos.append(cc);central.append(ee)
    energy=np.array(central)
    # Accepted integrated beam partition owns the mean; central tracks own depth shapes.
    target=np.array([s['Se_J_per_true_C'],s['host_recoil_J_per_true_C'],s['direct_D_recoil_J_per_true_C'],s['residual_cutoff_J_per_true_C']])*QE
    energy*=target/np.dot([.85,.15],energy)
    energy[:,1]=frozen['eh']
    energy[:,0]=5000*QE-energy[:,1:].sum(axis=1)
    return np.array(cos),energy

def signed_inputs(c,l,frozen,m):
    xy,parent,groups=beam_sites(c,m['rung']);ns=len(xy);nm=len(m['omega']);h=c.DIMS[2]
    vel,_=maps(c,m,np.column_stack([xy,np.zeros(ns)]));pressure=vel[:,2].T
    # Point/subgrid projection; 20-nm form factor is unity to <1e-8 here.
    random=[];_,st=c.profiles_and_stopping()
    for tr in st['PdD',2.5]['tracks']:
        vv,_=maps(c,m,np.column_stack([xy,np.full(ns,c.deposition_depth(tr))]));a=vv[:,2].T
        random.append(a/np.linalg.norm(a,axis=0))
    # Exact z-polynomial moments of cosine diffusion, sampled at enough Gauss points.
    zz,ww=leggauss(8);zz=(zz+1)*h/2;ww=ww*h/2
    divergence=[]
    for z in zz:
        _,b=maps(c,m,np.column_stack([xy,np.full(ns,z)]));divergence.append(b[:,:3].sum(axis=1))
    divergence=np.array(divergence) # z,site,mode
    # Compress the polynomial depth dependence before every live update.
    poly=np.array([legval(2*zz/h-1,np.eye(7)[j]) for j in range(7)]).T
    coefficients=np.linalg.lstsq(poly,divergence.reshape(8,-1),rcond=None)[0].reshape(7,ns,nm)
    q,w=leggauss(256);z=(q+1)*h/2;w=w*h/2
    pol=np.array([legval(q,np.eye(7)[j]) for j in range(7)])
    jcos=(pol*w)@np.cos(np.arange(NZ)[:,None]*np.pi*z/h).T
    heatmap=ALPHA*(c.C11+2*c.C12)*(.01**2/ns)*coefficients.transpose(1,0,2)
    cos,energy=heat_profiles(c,frozen)
    m.update(xy=xy,parent=parent,source_groups=groups,random=np.array(random),pressure=pressure,heatmap=heatmap,jcos=jcos,heatcos=cos,heatenergy=energy)
    return m

def make_model(rung,c,l,frozen):
    if rung<2:m=old.consumer_model(rung,c,l,frozen)
    else:
        spatial=_space()
        d=np.load(dest(2)/'model.npz',allow_pickle=False);m=dict(rung=2,v=np.ascontiguousarray(d['eigenvectors']),omega=d['omega'],degree=(spatial.PX,spatial.PX,spatial.PZ))
        nm=len(m['omega']);omega=m['omega'];xy,groups=l.sites();h=c.DIMS[2]
        # Common fixed physical boxes and exact-enough composite cell quadrature.
        def box(side,z0,z1,center,quadratic=False):
            axes=[];weights=[]
            for a,(lo,hi) in enumerate(((center[0]-side/2,center[0]+side/2),(center[1]-side/2,center[1]+side/2),(z0,z1))):
                cuts=[lo,hi]
                if a<2:cuts=sorted([lo,hi]+[v for v in np.linspace(0,c.DIMS[a],spatial.CELLS+1) if lo<v<hi])
                q,w=leggauss(6 if a==2 else (spatial.PX+1 if quadratic else (spatial.PX+2)//2))
                axes.append(np.concatenate([(u+v)/2+(v-u)*q/2 for u,v in zip(cuts[:-1],cuts[1:])]))
                weights.append(np.concatenate([(v-u)*w/2 for u,v in zip(cuts[:-1],cuts[1:])]))
            pts=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3);ww=np.einsum('i,j,k->ijk',*weights).ravel()/(side**2*(z1-z0))
            velocity,b=maps(c,m,pts);strain=b/omega
            if not quadratic:return np.array([ww@strain[:,0],ww@velocity[:,2]/omega])
            cc=c.tensors(pts,c.DIMS[0]/8,'grains')
            return np.array([np.einsum('qim,qij,qjn,q->mn',strain,cc,strain,ww,optimize=True),c.RHO*np.einsum('qim,qin,q->mn',velocity,velocity,ww,optimize=True)])
        probes=np.array([box(.0025,0,50e-6,p) for p in xy]);m['eps']=probes[:,0];m['uz']=probes[:,1];m['groups']=groups
        m['region_bounds']=np.array([(side,z0,z1) for side in (.0025,.01,float(c.DIMS[0])) for z0,z1 in ((0,50e-6),(25e-6,75e-6),(175e-6,225e-6))])
        m['region_projectors']=np.array([box(*b,c.DIMS[:2]/2) for b in m['region_bounds']])
        m['regional']=np.array([box(.0025,z0,z1,c.DIMS[:2]/2,True) for z0,z1 in ((0,50e-6),(25e-6,75e-6),(175e-6,225e-6))])
        m['points']=np.array([[*(c.DIMS[:2]/2+[dx,dy]),z] for z in (30e-6,50e-6,200e-6) for dx in (-.0035,0,.0035) for dy in (-.0035,0,.0035)])
        m['velocity'],b=maps(c,m,m['points']);m['strain']=b/omega;m['cs']=c.tensors(m['points'],c.DIMS[0]/8,'grains')
        center=np.array([c.DIMS[0]/2+.0065,c.DIMS[1]/2,.000125]);half=np.array([.0005,.0005,50e-6]);normals=np.r_[np.eye(3),-np.eye(3)]
        vm,b=maps(c,m,center+normals*half);cc=c.tensors(center+normals*half,c.DIMS[0]/8,'grains');_,bc=maps(c,m,center[None]);bc=bc[0]/omega
        receiver=(np.einsum('qij,qjm->qim',cc,b/omega),vm,normals,np.array([4*half[1]*half[2],4*half[0]*half[2],4*half[0]*half[1]]*2),bc.T@c.tensors(center[None],c.DIMS[0]/8,'grains')[0]@bc,center,2*half)
        m['receiver']=receiver
        # Exactly the same integrated polarization/branch proxy as accepted rung1.
        m['speed']=np.load(dest(2)/'mode_metrics.npz',allow_pickle=False)['speed']
    signed_inputs(c,l,frozen,m)
    if rung==2:
        wt=np.array([.85,.15])*frozen['eh'];wt/=wt.sum();ref=l.source_reference(dict(omega=m['omega']),np.concatenate(m['random'],axis=1),np.repeat(wt/len(m['xy']),len(m['xy'])),m['eps'],m['uz'],float(frozen['ecanon'])/50e-6);m['ref']=ref
    return m

def save_model(m):
    fields=('omega','pressure','random','heatmap','jcos','heatcos','heatenergy','xy','parent','source_groups','eps','uz','groups','region_bounds','region_projectors','regional','points','strain','velocity','cs','speed')
    arrays={k:m[k] for k in fields};arrays.update({'ref_'+str(i):np.asarray(v) for i,v in enumerate(m['ref'])});arrays.update({'receiver_'+str(i):np.asarray(v) for i,v in enumerate(m['receiver'])})
    np.savez_compressed(dest(m['rung'])/f'signed_r{m["rung"]}.npz',**arrays)

def load_model(rung):
    d=dict(np.load(dest(rung)/f'signed_r{rung}.npz',allow_pickle=False));d['ref']=tuple(d.pop('ref_'+str(i)) for i in range(5));d['receiver']=tuple(d.pop('receiver_'+str(i)) for i in range(7));d['rung']=rung
    return d

def prepare(rung):
    c,l,r,f=old.accepted();m=make_model(rung,c,l,f);save_model(m)
    if rung==1:
        freq=m['omega']/(2*np.pi);participation=m['pressure']@profile(c,m['xy'],False)
        # Nearest stable mode with appreciable normal overlap (>=10% of band max).
        near=(freq>=90000)&(freq<=110000)
        keep=near&(abs(participation)>=.1*max(abs(participation[near])));j=np.where(keep)[0][np.argmin(abs(freq[keep]-100000))]
        # Literal adjacent stable eigenfrequencies around 100 kHz, not claimed antiresonance.
        lo=freq[freq<100000].max();hi=freq[freq>100000].min();between=(lo+hi)/2
        dump(HERE/'frequencies.json',dict(tuned_Hz=float(freq[j]),tuned_mode=int(j),normal_participation_per_sqrt_kg=float(participation[j]),between_Hz=float(between),between_neighbors_Hz=[float(lo),float(hi)],criterion='nearest to 100kHz in 90–110kHz with >=10% maximum even-spot participation in that local band',same_for_rungs=[0,1,2]))
    dump(dest(rung)/f'inputs_r{rung}.json',dict(base=BASE,rung=rung,producer_sha256=sha(HERE/'signed.py'),numeric_sha256=sha(dest(rung)/f'signed_r{rung}.npz'),accepted_inputs_sha256=old.EXPECTED,n_modes=len(m['omega']),source_cells=len(m['xy']),source_spacing_m=.01/np.sqrt(len(m['xy'])),source_footprint_sigma_m=20e-9,source_temporal_sigma_s=.5e-12,heat_modes=NZ,cp_J_kg_K=CP,conductivity_W_m_K=KTH,alpha_linear_per_K=ALPHA,diffusivity_m2_s=ATH,skin_100k_m=np.sqrt(ATH/(np.pi*1e5)),lateral_diffusion_length_record_m=np.sqrt(ATH*(BURN+SAMPLE)*DT),incident_power_W=.690,mean_incident_force_N=float(138e-6/QE*np.dot([.85,.15],MOMENTUM)),depth_component_energy_J=m['heatenergy'].tolist()))
    print(json.dumps(dict(prepared_rung=rung,modes=len(m['omega']),sites=len(m['xy']))),flush=True)

def timing(f,portable,record,beam):
    if beam<3:return f['time_weights'][2 if beam==2 else 0,:,record]
    freq=json.loads((HERE/'frequencies.json').read_text());hz=(100000,freq['tuned_Hz'],freq['between_Hz'])[beam-3]
    phase=portable.uniform(record,'beam',0,1)[0]*2*np.pi
    return 1+.7*np.cos(2*np.pi*hz*np.arange(BURN+SAMPLE)*DT+phase)

def transition(omega,gamma):
    a=gamma/2;wd=np.sqrt(omega**2-a*a);d=np.exp(-a*DT);s=np.sin(wd*DT)/wd;c=np.cos(wd*DT)
    return d*(c+a*s),d*s,-d*omega**2*s,d*(c-a*s)

def thermal_force(m,state):
    moments=state@m['jcos'].T
    force=np.empty((len(state),4,len(m['omega'])))
    for g in range(4):
        keep=m['source_groups']==g
        force[:,g]=moments[:,keep].reshape(len(state),-1)@m['heatmap'][keep].reshape(-1,len(m['omega']))
    return force

def run_arm(rung,beam,beta,loss,coherent,records):
    c,l,r,f=old.accepted();m=load_model(rung);p=r.Portable();nr=len(records);nm=len(m['omega']);ns=len(m['xy']);w=m['omega'];h=c.DIMS[2];area=.01**2/ns
    gamma=1/50e-6+m['speed']/.2 if loss=='retained' else 4460*np.log(10)/10*32*w/(2*np.pi*1e6)
    aa,ab,ba,bb=transition(w,gamma);weights=profile(c,m['xy'],beam in (1,2))
    # Preserve the exact accepted parent spatial totals even when child quadrature refines.
    parentweights=f['site_weights'][2 if beam in (1,2) else 0]
    childprob=np.array([weights[m['parent']==i]/weights[m['parent']==i].sum() for i in range(36)])
    if rung==2:weights=beam_weights(c,m,beam in (1,2),f)
    count=r.CountRng(p,records);total=r.TotalRng(p,records);marks=r.MarkRng(p,records,BURN)
    times=np.empty((BURN+SAMPLE,nr),np.int64)
    for interval,part in enumerate((slice(0,BURN),slice(BURN,BURN+SAMPLE))):
        for ri,record in enumerate(records):
            n=total.poisson(138e-6/QE*(part.stop-part.start)*DT);wt=timing(f,p,record,beam)[part];times[part,ri]=count.multinomial(n,wt/wt.sum())
    q=np.zeros((nr,3,4,nm));v=np.zeros_like(q);temperature=np.zeros((nr,ns,NZ));prevQ=np.zeros((nr,4,nm))
    lam=ATH*(np.pi*np.arange(NZ)/h)**2;dec=np.exp(-lam*DT);u=np.ones(NZ);u[1:]=-np.expm1(-lam[1:]*DT)/(lam[1:]*DT)
    gain=np.ones(NZ);gain[1:]=2;gain/=RHO*CP*area*h
    ecanon=float(m['ref'][1]);ve=np.broadcast_to(m['ref'][2]**2,(nr,36)).copy();vu=np.broadcast_to(m['ref'][3]**2,(nr,36)).copy();adapt=-np.expm1(-DT/(5*max(1/gamma)))
    names=('all','A_disturbed','C_any','C_strain_positive','C_strain_negative','C_displacement_positive','C_displacement_negative','A_and_C','A_or_C');counts={n:np.zeros(nr,np.int64) for n in names}
    sx=np.empty((SAMPLE,nr,3,nm));sy=np.empty_like(sx);Bs=np.empty((SAMPLE,nr),bool);sig=np.empty((SAMPLE,nr,18));zs=np.empty_like(sig);rv=np.broadcast_to(ecanon*np.einsum('im,mn,in->i',m['region_projectors'].reshape(-1,nm),m['ref'][0][:nm,:nm],m['region_projectors'].reshape(-1,nm)),(nr,18)).copy()
    inc=np.zeros(nr);hostsum=np.zeros(nr);works=np.zeros((nr,3));losses=np.zeros(nr);marked=np.zeros(nr);cross=np.zeros(nr);initial=None;capped=np.zeros(nr,int);max_temp=np.zeros(nr)
    rp=m['region_projectors'].reshape(-1,nm);frac=.01/.865*np.array([.1,1.,10.]);Qold=np.zeros((nr,4,nm))
    # Fixed mean inward force is the acceleration-balanced static reference.
    mean_force=np.zeros((4,nm))
    for g in range(4):
        keep=m['source_groups']==g;mean_force[g]=138e-6/QE*np.dot([.85,.15],MOMENTUM)*(m['pressure'][:,keep]@weights[keep])
    heat_source=m['heatenergy'][:,:,None]*m['heatcos'];mean_pressure=mean_force if coherent else np.zeros_like(mean_force)
    for ti in range(BURN+SAMPLE):
        Qold=thermal_force(m,temperature) if coherent else np.zeros_like(Qold)
        bias=np.stack([np.zeros_like(Qold),np.broadcast_to(mean_pressure,(nr,4,nm)),Qold],axis=1)/w**2
        if ti==0:bias[:]=0
        dyn=q-bias
        vd=v.copy();vd[:,2]-=(Qold-prevQ)/(DT*w**2);prevQ=Qold.copy()
        pre=dyn.sum(axis=(1,2))*w
        ev=pre@m['eps'].T;uv=pre@m['uz'].T
        masks=l.source_masks(ev,uv,np.sqrt(ve),np.sqrt(vu),m['ref'][2],m['ref'][3]);ve=(1-adapt)*ve+adapt*ev*ev;vu=(1-adapt)*vu+adapt*uv*uv
        reg=pre@rp.T;zreg=reg/np.sqrt(np.maximum(rv,np.finfo(float).tiny));rv=(1-adapt)*rv+adapt*reg*reg
        if ti>=BURN:
            j=ti-BURN
            if initial is None:initial=.5*np.sum((q.sum(axis=(1,2))*w)**2+v.sum(axis=(1,2))**2,axis=1)+RHO*CP*area*h*temperature[:,:,0].sum(axis=1)
            sx[j]=dyn.sum(axis=2)*w;sy[j]=vd.sum(axis=2);sig[j]=reg;zs[j]=zreg
            Bs[j]=l.convergence(dyn.sum(axis=1)*w,vd.sum(axis=1),m['receiver'],1.)[0]
        parents=np.empty((nr,36,2,2),np.int64)
        for ri,record in enumerate(records):
            site=count.multinomial(int(times[ti,ri]),parentweights);dp=count.binomial(site,.85);species=np.column_stack([dp,site-dp]);pos=count.binomial(species,.5);parents[ri,:,:,0]=pos;parents[ri,:,:,1]=species-pos
        highp=.01+beta*.5*(np.tanh(ev/m['ref'][2])**2+np.tanh(uv/m['ref'][3])**2)
        low=marks.binomial(parents,.25);high=marks.binomial(parents-low,highp[:,:,None,None]/.75);real=np.stack([low,parents-low-high,high],axis=-1)
        if ns==36:child=real
        else:
            child=np.empty((nr,ns,2,2,3),np.int64)
            ids=np.array([np.where(m['parent']==cell)[0] for cell in range(36)])
            # Vectorized conditional binomials preserve every old parent/mark.
            for ri,record in enumerate(records):
                remaining=real[ri].copy();prob_left=np.ones(36)
                for sub in range(ids.shape[1]-1):
                    prob=childprob[:,sub]/prob_left
                    draw=p.binomial(remaining,prob[:,None,None,None],record,'subsite',ti*16+sub)
                    child[ri,ids[:,sub]]=draw;remaining-=draw;prob_left-=childprob[:,sub]
                child[ri,ids[:,-1]]=remaining
        signed=child[:,:,:,0]-child[:,:,:,1];amp=np.sum(signed*np.sqrt(2*f['eh'][None,None,:,None]*frac),axis=-1)
        dr=np.zeros((nr,4,nm));dm=np.zeros_like(dr)
        byspecies=child.sum(axis=(3,4))
        for g in range(4):
            keep=m['source_groups']==g
            for sp in range(2):dr[:,g]+=amp[:,keep,sp]@m['random'][sp][:,keep].T
            if coherent:dm[:,g]=(byspecies[:,keep]@MOMENTUM)@m['pressure'][:,keep].T
        vt=v.sum(axis=(1,2));dmt=dm.sum(axis=1);wpm=np.sum(vt*dmt,axis=1)+.5*np.sum(dmt*dmt,axis=1);avail=np.sum(parents*f['eh'][None,None,:,None],axis=(1,2,3))-wpm
        drt=dr.sum(axis=1);a=.5*np.sum(drt*drt,axis=1);b=np.sum((vt+dmt)*drt,axis=1);cap=np.ones(nr);over=a+b>avail
        cap[over]=2*avail[over]/(b[over]+np.sqrt(b[over]**2+4*a[over]*avail[over]));dr*=cap[:,None,None];drt=dr.sum(axis=1)
        dtotal=drt+dmt
        # Symmetric work allocation preserves random/pressure interference explicitly.
        wr=np.sum(vt*drt,axis=1)+.5*np.sum(drt*dtotal,axis=1);wp=np.sum(vt*dmt,axis=1)+.5*np.sum(dmt*dtotal,axis=1)
        v[:,0]+=dr;v[:,1]+=dm
        deposit=byspecies@heat_source.reshape(2,-1);deposit=deposit.reshape(nr,ns,4,NZ).sum(axis=2)
        # Remove actual mechanical work once, not marked self energy a second time.
        remove=(wr+wp)[:,None,None]*weights[None,:,None]*np.average(m['heatcos'][:,1],axis=0,weights=[.85,.15])[None,None,:]
        deposit-=remove
        source=deposit*gain
        average=temperature*u+source*np.where(lam>0,(1-u)/np.where(lam>0,lam*DT,1),.5)
        temperature=temperature*dec+source*u
        thermal=thermal_force(m,average) if coherent else np.zeros_like(Qold)
        force=np.zeros_like(q);force[:,2]=thermal
        before=.5*np.sum((q.sum(axis=(1,2))*w)**2+v.sum(axis=(1,2))**2,axis=1)
        oldq=q.copy();eq=force/w**2
        q=aa*(q-eq)+ab*v+eq;v=ba*(oldq-eq)+bb*v
        wt=np.sum(thermal.sum(axis=1)*(q-oldq).sum(axis=(1,2)),axis=1)
        after=.5*np.sum((q.sum(axis=(1,2))*w)**2+v.sum(axis=(1,2))**2,axis=1)
        damping=before+wt-after
        temperature[:,:,0]-=wt[:,None]*weights/(RHO*CP*area*h)
        if ti>=BURN:
            totalparents=parents.sum(axis=(1,2,3));counts['all']+=totalparents
            for name in names[1:]:counts[name]+=np.sum(parents*masks[name][:,:,None,None],axis=(1,2,3))
            inc+=totalparents*5000*QE;hostsum+=np.sum(parents*f['eh'][None,None,:,None],axis=(1,2,3));works+=np.column_stack([wr,wp,wt]);losses+=damping
            marked+=np.sum(real*f['eh'][None,None,:,None,None]*frac*cap[:,None,None,None,None]**2,axis=(1,2,3,4));capped+=over;cross+=np.sum(drt*dmt,axis=1)
            max_temp=np.maximum(max_temp,np.max(np.sum(temperature,axis=2),axis=1))
    final=.5*np.sum((q.sum(axis=(1,2))*w)**2+v.sum(axis=(1,2))**2,axis=1)+RHO*CP*area*h*temperature[:,:,0].sum(axis=1)
    residual=final-initial-inc+losses
    if max(abs(residual)/inc)>2e-7:raise AssertionError(('energy balance',residual/inc))
    return m,dict(x=sx,y=sy,B=Bs,signed=sig,z=zs,counts=counts,incident=inc,host=hostsum,work=works,damping=losses,marked=marked,cross=cross,residual=residual,capped=capped,max_temp=max_temp,gamma=gamma,digests=[v.hexdigest() for v in marks.sha],final_temperature=temperature)

def read_arm(rung,beam,beta,loss,coherent,records,m,z):
    c,l,_,_=old.accepted();w=m['omega'];arm=f'r{rung}_{BEAMS[beam]}_b{beta}_{loss}_{"combined" if coherent else "random"}'
    if rung==2:spatial=_space()
    x=z['x'].sum(axis=2);y=z['y'].sum(axis=2);nr=len(records);root=np.sqrt(z['gamma'])
    energy=.5*np.einsum('trm,pmn,trn->trp',x,m['regional'][:,0],x,optimize=True)
    production=.5*(np.einsum('trm,pmn,trn->trp',x*root,m['regional'][:,0],x*root,optimize=True)+np.einsum('trm,pmn,trn->trp',y*root,m['regional'][:,1],y*root,optimize=True))
    localstrain=np.einsum('pim,trm->trpi',m['strain'],x);local=.5*np.einsum('trpi,pij,trpj->trp',localstrain,m['cs'],localstrain)
    se=np.einsum('pim,trm->trpi',m['strain'],x*root);sv=np.einsum('pim,trm->trpi',m['velocity'],y*root)
    localrate=.5*np.einsum('trpi,pij,trpj->trp',se,m['cs'],se)+.5*RHO*np.sum(sv*sv,axis=-1)
    prof=c.profiles_and_stopping()[0]['PdD',2.5];rows=[];tails=[];signed=[];ledger=[];same=[];peaks=[];full=None
    for ri,record in enumerate(records):
        row=dict(arm=arm,rung=rung,beam=beam,beta=beta,loss=loss,coherent=coherent,record=record,exposure_s=SAMPLE*DT,draw_sha256=z['digests'][ri])
        for n,val in z['counts'].items():row[n+'_parents']=int(val[ri]);row[n+'_fraction']=float(val[ri]/z['counts']['all'][ri])
        row.update({'B_'+k:v for k,v in l.occurrence(z['B'][:,ri:ri+1],DT).items() if k in ('time_fraction','complete_events','boundary_censored_events','initiation_rate_s')});rows.append(row)
        led=dict(arm=arm,record=record,incident_J=z['incident'][ri],host_J=z['host'][ri],random_work_J=z['work'][ri,0],momentum_work_J=z['work'][ri,1],thermal_work_J=z['work'][ri,2],random_pressure_kick_cross_J=z['cross'][ri],marked_self_J=z['marked'][ri],marked_host_fraction=z['marked'][ri]/z['host'][ri],loss_J=z['damping'][ri],balance_residual_J=z['residual'][ri],max_source_column_surface_rise_K=z['max_temp'][ri],capped_bins=z['capped'][ri]);ledger.append(led)
        for region,(depth,z0,z1) in enumerate(((30e-6,0,50e-6),(50e-6,25e-6,75e-6),(200e-6,175e-6,225e-6))):
            u=energy[:,ri,region];mean=u.mean();at=int(u.argmax());nd=c.column_density(prof['layers'],z0,z1)
            tail=dict(arm=arm,record=record,region=region,depth_m=depth,mean_J_m3=float(mean),max_J_m3=float(u[at]),peak_over_own_mean=float(u[at]/mean),elastic_eV_D=float(u[at]/(QE*nd)),conditional_eV_D=float(production[at,ri,region]*4.1e-12/(QE*nd)))
            for quant in (.5,.9,.95,.99,.999):tail['q'+str(quant)+'_over_mean']=float(np.quantile(u,quant)/mean)
            for level in (3,4,5):
                occ=l.occurrence((u>=level*mean)[:,None],DT)
                for key in ('time_fraction','complete_events','boundary_censored_events','initiation_rate_s','max_duration_us'):tail[f'tail{level}_{key}']=occ[key]
            tails.append(tail)
        at,region=np.unravel_index(energy[:,ri].argmax(),energy[:,ri].shape);u=energy[:,ri,region];lo=hi=at
        while lo and u[lo-1]>=u[at]/2:lo-=1
        while hi<SAMPLE-1 and u[hi+1]>=u[at]/2:hi+=1
        comp=z['x'][at,ri];a=m['regional'][region,0]
        diagonal=.5*np.einsum('gm,mn,gn->g',comp,a,comp);crossvalue=float(u[at]-diagonal.sum())
        peak=dict(arm=arm,record=record,region=int(region),sample=int(at),time_s=at*DT,peak_J_m3=float(u[at]),peak_over_own_mean=float(u[at]/u.mean()),duration_s=(hi-lo+1)*DT,duration_censored=bool(lo==0 or hi==SAMPLE-1),random_diagonal_J_m3=float(diagonal[0]),momentum_diagonal_J_m3=float(diagonal[1]),thermal_diagonal_J_m3=float(diagonal[2]),all_component_cross_J_m3=crossvalue,signed_diagonal_fraction=float((diagonal[1]+diagonal[2])/u[at]),B_distinct_origin_window=bool(z['B'][at,ri]),attribution='mixed/unresolved; disjoint origin sectors are not unique individual impacts')
        # Representative fixed records for expensive spatial maps, every arm retained.
        if record in (0,1):
            if full is None:
                if rung==0:
                    full=dict(rung=0,v=old.accepted()[3]['eigenvectors'],degree=(5,5,2))
                else:
                    data=np.load((HERE.parent/'rung1' if rung==1 else dest(2))/'model.npz',allow_pickle=False)
                    vectors=data['eigenvectors']
                    full=dict(rung=rung,v=vectors,degree=tuple(data['degree']) if rung==1 else (spatial.PX,spatial.PX,spatial.PZ))
            snapshot=dict(full,v=(full['v']@(x[at,ri]/w))[:,None])
            axis=np.linspace(0,c.DIMS[0],129);xx,yy=np.meshgrid(axis,axis,indexing='ij');vals=[]
            for start in range(0,xx.size,128):
                pts=np.column_stack([xx.ravel()[start:start+128],yy.ravel()[start:start+128],np.full(min(128,xx.size-start),(30e-6,50e-6,200e-6)[region])]);_,bb=maps(c,snapshot,pts);strain=bb[:,:,0];cc=c.tensors(pts,c.DIMS[0]/8,'grains');vals.extend(.5*np.einsum('qi,qij,qj->q',strain,cc,strain))
            field=np.array(vals).reshape(xx.shape);inside=(abs(xx-c.DIMS[0]/2)<=.00125)&(abs(yy-c.DIMS[1]/2)<=.00125);ij=np.unravel_index(np.where(inside,field,-np.inf).argmax(),field.shape);lab,_=label(field>=field[ij]/2);coords=np.argwhere(lab==lab[ij]);width=(coords.max(axis=0)-coords.min(axis=0))*(axis[1]-axis[0]);lam=(.01373147,.01373147/3,.01373147/9)[rung]
            peak.update(width_x_m=float(width[0]),width_y_m=float(width[1]),width_foil_censored=bool(np.any(coords==0)|np.any(coords==128)),width_band_censored=bool(min(width)<lam/2),width_map_unresolved=bool(min(width)<2*(axis[1]-axis[0])))
        peaks.append(peak)
        for di,(depth,z0,z1) in enumerate(((30e-6,0,50e-6),(50e-6,25e-6,75e-6),(200e-6,175e-6,225e-6))):
            pi=9*di+4;nd=c.layer_value(prof['layers'],depth);ndr=c.column_density(prof['layers'],z0,z1)
            same.append(dict(arm=arm,record=record,sample=int(at),depth_m=depth,local_elastic_eV_D=float(local[at,ri,pi]/(QE*nd)),local_conditional_eV_D=float(localrate[at,ri,pi]*4.1e-12/(QE*nd)),region_elastic_eV_D=float(energy[at,ri,di]/(QE*ndr)),region_conditional_eV_D=float(production[at,ri,di]*4.1e-12/(QE*ndr))))
        for proj,bound in enumerate(m['region_bounds']):
            for oi,name in enumerate(('epsilon_xx','u_z')):
                if np.linalg.norm(m['region_projectors'][proj,oi])<1e-10*np.linalg.norm(m['region_projectors'][:,oi]):continue
                values=z['signed'][:,ri,2*proj+oi];zz=z['z'][:,ri,2*proj+oi]
                for sign,lab in ((1,'crest'),(-1,'hole')):
                    occ=l.occurrence((sign*zz>=3)[:,None],DT);signed.append(dict(arm=arm,record=record,side_m=bound[0],z0_m=bound[1],z1_m=bound[2],observable=name,sign=lab,min=float(values.min()),max=float(values.max()),**{k:occ[k] for k in ('time_fraction','complete_events','boundary_censored_events','initiation_rate_s','max_duration_us')}))
    out=dest(rung);tag=f'{arm}_records{records[0]}-{records[-1]}'
    for name,data in (('records',rows),('tails',tails),('peaks',peaks),('signed',signed),('energy',ledger),('same_event',same)):table(out/f'{tag}_{name}.csv',data)
    # Full sampled fixed-region tails retained, rather than only selected quantiles.
    np.savez_compressed(out/f'{tag}_histories.npz',energy_J_m3=energy,production_W_m3=production,local_J_m3=local,local_production_W_m3=localrate)
    print(json.dumps(dict(arm=arm,records=records,max_relative_energy_residual=float(max(abs(z['residual'])/z['incident'])))),flush=True)

def run(rung,beam,beta,loss,coherent,first,count):
    records=list(range(first,first+count));m,z=run_arm(rung,beam,beta,loss,coherent,records);read_arm(rung,beam,beta,loss,coherent,records,m,z)

def summarize(rung):
    out=dest(rung)
    def read(kind):return [row for path in sorted(out.glob(f'r{rung}_*_records*-*_{kind}.csv')) for row in csv.DictReader(path.open())]
    rr=read('records');tails=read('tails');peaks=read('peaks');energy=read('energy');summary=[]
    for arm in sorted({p['arm'] for p in rr}):
        records=[p for p in rr if p['arm']==arm]
        for metric in [key for key in records[0] if key.endswith('_fraction') and key!='all_fraction']+['B_initiation_rate_s']:
            summary.append(dict(arm=arm,observable=metric,**stats([p[metric] for p in records])))
        for key in (key for key in records[0] if key.endswith('_parents')):
            summary.append(dict(arm=arm,observable=key[:-8]+'_parent_strike_rate_s',**stats([float(p[key])/float(p['exposure_s']) for p in records])))
        for metric in ('peak_over_own_mean','duration_s','signed_diagonal_fraction','all_component_cross_J_m3'):
            summary.append(dict(arm=arm,region='leading event',observable=metric,**stats([p[metric] for p in peaks if p['arm']==arm])))
        for region in range(3):
            for metric in ('mean_J_m3','peak_over_own_mean','elastic_eV_D','conditional_eV_D'):
                summary.append(dict(arm=arm,region=region,observable=metric,**stats([p[metric] for p in tails if p['arm']==arm and int(p['region'])==region])))
    table(out/f'summary_r{rung}.csv',summary)
    table(out/f'energy_r{rung}.csv',energy)
    # Paired source effect and spatial-rung effect use record as the unit.
    contrasts=[]
    for arm in sorted({p['arm'] for p in tails if p['arm'].endswith('_combined')}):
        controls=[arm.replace('_combined','_random')]
        if rung:controls.append('r'+str(rung-1)+arm[2:])
        for control in controls:
            other=tails if control.startswith(f'r{rung}_') else [row for path in sorted(dest(rung-1).glob(f'{control}_records*-*_tails.csv')) for row in csv.DictReader(path.open())]
            for region in range(3):
                for metric in ('mean_J_m3','peak_over_own_mean','elastic_eV_D'):
                    a={int(p['record']):float(p[metric]) for p in tails if p['arm']==arm and int(p['region'])==region};b={int(p['record']):float(p[metric]) for p in other if p['arm']==control and int(p['region'])==region};keys=sorted(set(a)&set(b))
                    if keys:contrasts.append(dict(arm=arm,control=control,region=region,observable=metric,comparison='paired per-record ratio',**stats([a[k]/b[k] for k in keys])))
            otherpeaks=peaks if control.startswith(f'r{rung}_') else [row for path in sorted(dest(rung-1).glob(f'{control}_records*-*_peaks.csv')) for row in csv.DictReader(path.open())]
            for metric in ('peak_over_own_mean','duration_s'):
                a={int(p['record']):float(p[metric]) for p in peaks if p['arm']==arm};b={int(p['record']):float(p[metric]) for p in otherpeaks if p['arm']==control};keys=sorted(set(a)&set(b))
                if keys:
                    result=stats([a[k]/b[k] for k in keys])
                    contrasts.append(dict(arm=arm,control=control,region='leading event',observable=metric,comparison='paired per-record ratio',within_10percent_equivalence_window=bool(result['ci95_lo']>=.9 and result['ci95_hi']<=1.1),**result))
    table(out/f'contrasts_r{rung}.csv',contrasts)
    # The delivered peak diagonals plus the total cross term are exact. Bound
    # unresolved pressure/thermal cross allocation without pretending diagonals add.
    share=[]
    for p in peaks:
        e=float(p['peak_J_m3']);er=float(p['random_diagonal_J_m3']);em=float(p['momentum_diagonal_J_m3']);et=float(p['thermal_diagonal_J_m3'])
        lo=(np.sqrt(max(em,0))-np.sqrt(max(et,0)))**2;hi=(np.sqrt(max(em,0))+np.sqrt(max(et,0)))**2
        share.append(dict(arm=p['arm'],record=p['record'],signed_only_elastic_fraction_lower=lo/e,signed_only_elastic_fraction_upper=hi/e,random_signed_cross_J_m3_lower=e-er-hi,random_signed_cross_J_m3_upper=e-er-lo,signed_symmetric_allocation_fraction_lower=(e-er+lo)/(2*e),signed_symmetric_allocation_fraction_upper=(e-er+hi)/(2*e),interpretation='algebraic Cauchy bounds at the same peak; not a positive exclusive origin probability'))
    table(out/f'same_sign_share_bounds_r{rung}.csv',share)
    ids={arm:sorted(int(p['record']) for p in rr if p['arm']==arm) for arm in sorted({p['arm'] for p in rr})}
    expected={f'r{rung}_{beam}_b{beta}_{loss}_{source}' for beam in BEAMS for beta in (0.,.35) for loss in ('retained','Cu_proxy') for source in ('random','combined')}
    m=load_model(rung);null=[]
    for proj,bound in enumerate(m['region_bounds']):
        for oi,name in enumerate(('epsilon_xx','u_z')):
            relative=float(np.linalg.norm(m['region_projectors'][proj,oi])/np.linalg.norm(m['region_projectors'][:,oi]))
            if relative<1e-10:null.append(dict(side_m=float(bound[0]),z0_m=float(bound[1]),z1_m=float(bound[2]),observable=name,relative_operator_norm=relative,reason='near-null fixed projection: normalized signed-tail statistic undefined, not an omitted source arm'))
    dump(out/f'verification_r{rung}.json',dict(base=BASE,rung=rung,arms=len(ids),records=len(rr),independent_records=sorted({int(p['record']) for p in rr}),record_ids_by_arm=ids,complete_expected_matrix=(set(ids)==expected and all(values==list(range(NREC)) for values in ids.values())),null_signed_projections=null,max_relative_energy_residual=max(abs(float(p['balance_residual_J']))/float(p['incident_J']) for p in energy),width_representative_records=[0,1],uncertainty='pointwise Student t intervals across independent records, not time bins; no simultaneous/model interval',same_total_parent_counts=len({(p['record'],p['all_parents']) for p in rr})==len({p['record'] for p in rr})))

def package(rung):
    out=dest(rung);prefix=str(out.relative_to(ROOT))+'/'
    names=sorted(p.name for p in out.iterdir() if p.is_file() and p.name not in ('manifest.json','capsule.tar.gz','delivery.json'))
    files=[dict(path=prefix+n,bytes=(out/n).stat().st_size,sha256=sha(out/n)) for n in names]
    dependency=None
    if rung==2:
        delivered=json.loads((HERE/'delivery.json').read_text());prior=json.loads((HERE/'manifest.json').read_text())
        shared={str((HERE/name).relative_to(ROOT)) for name in ('signed.py','frequencies.json')}
        required=[row for row in prior['files'] if row['path'] in shared]
        if not delivered['all_member_hashes_read_back'] or {row['path'] for row in required}!=shared:raise AssertionError('unverified or incomplete rung1b dependency')
        if any(sha(ROOT/row['path'])!=row['sha256'] for row in required):raise AssertionError('delivered rung1b shared input changed')
        dependency=dict(capsule={key:delivered[key] for key in ('path','sha256','bytes')},manifest_sha256=sha(HERE/'manifest.json'),required_shared_files=required)
    manifest=dict(base=BASE,unpublished=True,files=files,dependency=dependency);dump(out/'manifest.json',manifest)
    archive=out/'capsule.tar.gz'
    with archive.open('wb') as raw:
        with gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as gz:
            with tarfile.open(fileobj=gz,mode='w|') as tf:
                for name in names+['manifest.json']:
                    path=out/name;info=tarfile.TarInfo(prefix+name);info.size=path.stat().st_size;info.mode=0o644;info.mtime=0
                    with path.open('rb') as payload:tf.addfile(info,payload)
    # Actual bytes read back and every member authenticated, not writer-return alone.
    expected={row['path']:row for row in files};manifest_seen=False
    with tarfile.open(archive,'r|gz') as tf:
        for member in tf:
            if not member.isfile():raise AssertionError('non-regular capsule member')
            if member.name==prefix+'manifest.json':
                with tf.extractfile(member) as payload:stored=json.load(payload)
                if stored!=manifest:raise AssertionError('manifest readback')
                manifest_seen=True;continue
            row=expected.pop(member.name)
            with tf.extractfile(member) as payload:digest=hashlib.file_digest(payload,'sha256').hexdigest()
            if member.size!=row['bytes'] or digest!=row['sha256']:raise AssertionError('capsule readback')
    if expected or not manifest_seen:raise AssertionError('incomplete capsule readback')
    receipt=dict(path=prefix+'capsule.tar.gz',sha256=sha(archive),bytes=archive.stat().st_size,members=len(names)+1,all_member_hashes_read_back=True)
    dump(out/'delivery.json',receipt);print(json.dumps(receipt),flush=True)

def transfer(rung,m=None,c=None):
    """Read the actual signed modal gains and the thermal projection truncation."""
    if m is None:
        c,l,r,f=old.accepted();m=load_model(rung)
    else:f=old.accepted()[3]
    w=m['omega'];out=dest(rung);frequencies=json.loads((HERE/'frequencies.json').read_text())
    cos128,energy128=heat_profiles(c,f,128);hzquad,hzweight=leggauss(512);h=c.DIMS[2]
    pol=np.array([legval(hzquad,np.eye(7)[j]) for j in range(7)])
    extended=dict(m,jcos=(pol*hzweight*h/2)@np.cos(np.arange(128)[:,None]*np.pi*(hzquad+1)/2).T)
    thermal_checks=[]
    transfer_rows=[]
    for hot in (False,True):
        weights=beam_weights(c,m,hot,f)
        pressure=138e-6/QE*np.dot([.85,.15],MOMENTUM)*(m['pressure']@weights)
        source=np.dot([.85,.15],(m['heatenergy'][:,:,None]*m['heatcos']).sum(axis=1))
        for hz in (100000.,frequencies['tuned_Hz'],frequencies['between_Hz']):
            om=2*np.pi*hz;lam=ATH*(np.pi*np.arange(NZ)/c.DIMS[2])**2
            state=np.broadcast_to((138e-6/QE*source*np.r_[1.,np.full(NZ-1,2.)]/(RHO*CP*(.01**2)*c.DIMS[2])/(lam+1j*om))[None,None,:],(1,len(weights),NZ)).copy()
            state*=weights[None,:,None]*len(weights)
            # Complex linearity; same operator, not a separate thermal-force law.
            thermal=thermal_force(m,state.real).sum(axis=1)[0]+1j*thermal_force(m,state.imag).sum(axis=1)[0]
            src128=np.dot([.85,.15],(energy128[:,:,None]*cos128).sum(axis=1))
            lam128=ATH*(np.pi*np.arange(128)/h)**2
            t128=138e-6/QE*src128*np.r_[1.,np.full(127,2.)]/(RHO*CP*.01**2*h)/(lam128+1j*om)
            state128=weights[None,:,None]*len(weights)*t128[None,None,:]
            fine=thermal_force(extended,state128.real).sum(axis=1)[0]+1j*thermal_force(extended,state128.imag).sum(axis=1)[0]
            for loss in ('retained','Cu_proxy'):
                gamma=1/50e-6+m['speed']/.2 if loss=='retained' else 4460*np.log(10)/10*32*w/(2*np.pi*1e6)
                response=(pressure+thermal)/(w*w-om*om+1j*gamma*om)
                j=int(np.argmax(abs(response)))
                selected=int(np.argmin(abs(w/(2*np.pi)-frequencies['tuned_Hz'])))
                transfer_rows.append(dict(rung=rung,spatial='hot' if hot else 'even',frequency_Hz=hz,loss=loss,largest_response_mode=j,mode_frequency_Hz=w[j]/(2*np.pi),normal_participation=pressure[j]/(138e-6/QE*np.dot([.85,.15],MOMENTUM)),thermal_to_pressure_abs=float(abs(thermal[j])/max(abs(pressure[j]),1e-300)),thermal_phase_rad=float(np.angle(thermal[j])),settled_modal_amplitude_sqrtkg_m=float(abs(response[j])),ringup_amplitude_s=2/gamma[j],finite_exposure_s=SAMPLE*DT,selected_near_mode=selected,selected_near_frequency_Hz=w[selected]/(2*np.pi),selected_normal_participation=pressure[selected]/(138e-6/QE*np.dot([.85,.15],MOMENTUM)),selected_thermal_force_abs=float(abs(thermal[selected])),selected_thermal_phase_rad=float(np.angle(thermal[selected])),selected_modal_amplitude_sqrtkg_m=float(abs(response[selected]))))
                denom=w*w-om*om+1j*gamma*om
                error=float(np.linalg.norm((fine-thermal)/denom)/np.linalg.norm(fine/denom))
                thermal_checks.append(dict(spatial='hot' if hot else 'even',frequency_Hz=hz,loss=loss,relative_64_vs_128_thermal_displacement_norm=error,uniform_spot_surface_128_peak_K=float(abs(t128.sum())),semi_infinite_surface_peak_K=.690/.01**2/KTH*np.sqrt(ATH/om),surface_phase_rad=float(np.angle(t128.sum())),source_modulation='unit peak modulation; actual ripple amplitude .7'))
    table(out/f'transfer_r{rung}.csv',transfer_rows)
    table(out/f'thermal_projection_r{rung}.csv',thermal_checks)

def static(rung,m=None,c=None):
    """Balanced full-space DC solve with an explicitly ideal thermal flange."""
    from scipy.linalg import eigh,solve
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    from scipy.interpolate import RegularGridInterpolator
    if m is None:c,l,r,f=old.accepted();m=load_model(rung)
    if rung==2:spatial=_space()
    out=dest(rung);axis=np.linspace(0,c.DIMS[0],129);step=axis[1]-axis[0]
    xx,yy=np.meshgrid(axis,axis,indexing='ij');rows=[];n=127
    one=sparse.diags([-np.ones(n-1),2*np.ones(n),-np.ones(n-1)],[-1,0,1])/step**2
    lap=sparse.kron(one,sparse.eye(n))+sparse.kron(sparse.eye(n),one)
    degree=(5,5,2) if rung==0 else (15,15,4)
    if rung<2:
        dcK=old.stiffness(c,degree);lam0,rigid=eigh(dcK,subset_by_index=(0,5),check_finite=False)
        scale=float(np.trace(dcK)/len(dcK));balanced=dcK+scale*(rigid@rigid.T);ndc=len(dcK)
    else:
        global_solver=spatial.StaticSolver(c.DIMS);ndc=global_solver.k.shape[0]
        degree='full global Q4 x Q4 x P4; six rigid gauges, no frequency cutoff'
    qz,wz=leggauss(8);zpts=(qz+1)*c.DIMS[2]/2;wz*=c.DIMS[2]/2
    qxy,wxy=leggauss(5);edges=np.linspace(0,c.DIMS[0],17)
    xs=np.concatenate([(a+b)/2+(b-a)*qxy/2 for a,b in zip(edges[:-1],edges[1:])]);ws=np.tile(wxy*c.DIMS[0]/32,16)
    xy=np.stack(np.meshgrid(xs,xs,indexing='ij'),axis=-1).reshape(-1,2);ww=np.outer(ws,ws).ravel()
    for hot in (False,True):
        pxy=np.column_stack([xx.ravel(),yy.ravel()]);density=profile(c,pxy,hot).reshape(xx.shape) if hot else np.ones_like(xx)
        spot=(abs(xx-c.DIMS[0]/2)<=.005)&(abs(yy-c.DIMS[1]/2)<=.005);density*=spot;density/=density.sum()*step**2
        power=.690*density;T=np.zeros_like(xx);T[1:-1,1:-1]=spsolve(lap,(power[1:-1,1:-1]/(KTH*c.DIMS[2])).ravel()).reshape(n,n)
        tint=RegularGridInterpolator((axis,axis),T);pint=RegularGridInterpolator((axis,axis),power)
        xybeam=m['xy'];force=138e-6/QE*np.dot([.85,.15],MOMENTUM)
        if rung<2:
            thermalQ=np.zeros(ndc)
            for zz,wzz in zip(zpts,wz):
                pts=np.column_stack([xy,np.full(len(xy),zz)]);temp=tint(xy)+pint(xy)/KTH*(zz*zz/(2*c.DIMS[2])-zz+c.DIMS[2]/3)
                _,b=c.basis(pts,degree=degree)
                thermalQ+=ALPHA*(c.C11+2*c.C12)*np.einsum('q,qm->m',ww*wzz*temp,b[:,:3].sum(axis=1))
            phi,_=c.basis(np.column_stack([xybeam,np.zeros(len(xybeam))]),degree=degree);nb=phi.shape[1]
            pressureQ=np.zeros(ndc);pressureQ[2*nb:]=force*(beam_weights(c,m,hot)@phi)
            rhs=np.column_stack([pressureQ,thermalQ,pressureQ+thermalQ]);rhs-=rigid@(rigid.T@rhs)
            solution=solve(balanced,rhs,assume_a='pos',check_finite=False)
            residual=float(np.linalg.norm(dcK@solution-rhs)/np.linalg.norm(rhs))
            points=np.array([[c.DIMS[0]/2,c.DIMS[1]/2,0.],[.001,.001,0.]])
            ph,_=c.basis(points,degree=degree);uz=ph@solution[2*nb:]
            rigidthermal=float(np.linalg.norm(rigid.T@thermalQ))
        else:
            thermalQ=ALPHA*(c.C11+2*c.C12)*spatial.eigenstrain_load(c.DIMS,tint,lambda pts:pint(pts)/KTH)
            pressureQ=force*spatial.pressure_load(c.DIMS,xybeam,beam_weights(c,m,hot))
            rhs=np.column_stack([pressureQ,thermalQ,pressureQ+thermalQ])
            solution,residual=global_solver.solve(rhs)
            uz=spatial.maps(np.array([[c.DIMS[0]/2,c.DIMS[1]/2,0.],[.001,.001,0.]]),solution,c.DIMS)[0][:,2]
            rigidthermal=float(np.linalg.norm(global_solver.rigid.T@thermalQ))
        bends=uz[0]-uz[1];mass=c.RHO*np.prod(c.DIMS)
        centroid=beam_weights(c,m,hot)@(xybeam-c.DIMS[:2]/2);inertia=mass*(c.DIMS[0]**2+c.DIMS[2]**2)/12
        rows.append(dict(rung=rung,beam='hot_steady' if hot else 'even_steady',thermal_boundary='ideal isothermal lateral flange; insulated faces except incident heat; thin-aspect steady transfer',static_retained_heat_power_W=.690,max_temperature_rise_K=float(T.max()),mean_temperature_rise_K=float(T.mean()),CM_balanced_pressure_bend_m=float(bends[0]),CM_balanced_thermal_bend_m=float(bends[1]),CM_balanced_total_bend_m=float(bends[2]),rigid_translation_acceleration_m_s2=float(force/mass),rigid_rotation_acceleration_x_rad_s2=float(force*centroid[1]/inertia),rigid_rotation_acceleration_y_rad_s2=float(-force*centroid[0]/inertia),static_equilibrium_relative_residual=residual,mechanical_convention='distributed inertial reaction removes six rigid modes; not a measured clamped flange',DC_space=str(degree),thermal_rigid_force_norm=rigidthermal))
    table(out/f'static_r{rung}.csv',rows)

def temporal_discretization(rung):
    """Quantify actual pre-kick sampling, including its bounded DC artifact."""
    c,*_=old.accepted();m=load_model(rung);w=m['omega'];frequencies=json.loads((HERE/'frequencies.json').read_text());rows=[]
    F=138e-6/QE*np.dot([.85,.15],MOMENTUM)*(m['pressure']@beam_weights(c,m,False))
    for loss in ('retained','Cu_proxy'):
        gamma=1/50e-6+m['speed']/.2 if loss=='retained' else 4460*np.log(10)/10*32*w/(2*np.pi*1e6)
        a,b,d,e=transition(w,gamma);step=np.stack([np.stack([a,b*w],axis=-1),np.stack([d/w,e],axis=-1)],axis=1);right=step[:,:,1]*(w*DT)[:,None]
        steady=np.linalg.solve(np.eye(2)-step,right[:,:,None])[:,:,0]
        for hz in (100000.,frequencies['tuned_Hz'],frequencies['between_Hz']):
            om=2*np.pi*hz;sampled=np.linalg.solve(np.exp(1j*om*DT)*np.eye(2)-step,right[:,:,None])[:,0,0];analytic=w*w/(w*w-om*om+1j*gamma*om);ratio=sampled/analytic
            rows.append(dict(loss=loss,drive_Hz=hz,max_unit_force_harmonic_relative_difference=float(np.max(abs(ratio-1))),max_unit_force_phase_difference_rad=float(np.max(abs(np.angle(ratio)))),max_DC_strobe_displacement_relative_error=float(np.max(abs(steady[:,0]-1))),even_mean_pressure_DC_strobe_kinetic_J=float(.5*np.sum((F/w*steady[:,1])**2)),even_mean_pressure_DC_bias_error_elastic_J=float(.5*np.sum((F/w*(steady[:,0]-1))**2)),interpretation='50-ns parent-bin impulse approximation; bounded stationary sampling artifact, not growing material vibration'))
    table(dest(rung)/f'pressure_sampling_r{rung}.csv',rows)

def proof(rung):
    """Settled consumer checks: static balance, phase, work, causality, projection."""
    c,l,r,f=old.accepted();m=load_model(rung);w=m['omega'];out=dest(rung)
    frequencies=json.loads((HERE/'frequencies.json').read_text());transfer(rung,m,c)
    temporal_discretization(rung)
    # A real affine step must preserve a forced static equilibrium and match its
    # analytic harmonic transfer; wrong-sign thermal or all-positive mode fixes fail.
    gamma=1/50e-6+m['speed']/.2;aa,ab,ba,bb=transition(w,gamma)
    F=138e-6/QE*np.dot([.85,.15],MOMENTUM)*(m['pressure']@beam_weights(c,m,False));eq=F/w**2
    equilibrium_error=float(np.max(abs((aa*(eq-eq)+eq)-eq)))
    # First pre-bin field is zero. Its masks cannot depend on new-bin counts.
    zeros=np.zeros((1,36));mask=l.source_masks(zeros,zeros,m['ref'][2][None],m['ref'][3][None],m['ref'][2],m['ref'][3])
    if np.any(mask['A_disturbed']) or np.any(mask['C_any']):raise AssertionError('new source self-qualified')
    # A fixed harmonic oscillator, exact viscous propagator, forced near resonance.
    j=int(np.argmin(abs(w/(2*np.pi)-frequencies['tuned_Hz'])));om=w[j];g=gamma[j];a,b,d,e=transition(np.array([om]),np.array([g]))
    z=np.exp(1j*om*DT);mat=np.array([[a[0],b[0]],[d[0],e[0]]]);forcing=np.array([(1-a[0])/om**2,-d[0]/om**2])
    discrete=np.linalg.solve(z*np.eye(2)-mat,forcing)[0]
    analytic=1/(om**2-om**2+1j*g*om)
    # Left-constant bin forcing differs by the explicitly sampled half-bin phase.
    sampled_relative=float(abs(discrete/analytic-1))
    if sampled_relative>.08:raise AssertionError('near-mode affine transfer')
    static(rung,m,c)
    dump(out/f'behavior_r{rung}.json',dict(first_bin_A_C_empty=True,static_affine_equilibrium_absolute_error=equilibrium_error,known_mode_sampled_transfer_relative_difference=sampled_relative,force_signs_both_present=bool(np.any(F<0)&np.any(F>0)),thermal_zero_mode_insulated_growth_K_s=.690/(RHO*CP*np.prod(c.DIMS)),temperature_steady_state_under_insulation=False,flange_static_is_conditional=True,static_rung2_compliance_includes_modes_below_5kHz=(rung==2),dynamic_full_bulk_relation=True))
    # Full bulk Fourier-through-thickness witness at the shared kh~1 screen.
    branches=[]
    for theta in (0.,np.pi/4):
        values=[]
        for pz in (4,8):
            zz,ww=leggauss(pz+2);h=c.DIMS[2];v=[];d=[]
            for j in range(pz+1):
                coef=np.zeros(j+1);coef[j]=1;scale=np.sqrt((2*j+1)/(RHO*h))
                v.append(legval(zz,coef)*scale);d.append(legval(zz,legder(coef))*scale*2/h)
            v=np.array(v).T;d=np.array(d).T;n=pz+1;k=2*np.pi/(.01373147/9);b=np.zeros((len(zz),6,3*n),complex)
            b[:,0,:n]=1j*k*np.cos(theta)*v;b[:,1,n:2*n]=1j*k*np.sin(theta)*v;b[:,2,2*n:]=d
            b[:,3,n:2*n]=d;b[:,3,2*n:]=1j*k*np.sin(theta)*v;b[:,4,:n]=d;b[:,4,2*n:]=1j*k*np.cos(theta)*v;b[:,5,:n]=1j*k*np.sin(theta)*v;b[:,5,n:2*n]=1j*k*np.cos(theta)*v
            op=np.einsum('qim,ij,qjn,q->mn',b.conj(),c.cubic(),b,ww*h/2,optimize=True)
            values.append(np.sqrt(np.linalg.eigvalsh(op)[:3])/(2*np.pi))
        for branch in range(3):branches.append(dict(direction_rad=theta,branch_order=branch,kh=k*c.DIMS[2],P4_frequency_Hz=values[0][branch],P8_frequency_Hz=values[1][branch],relative_P4_change=values[0][branch]/values[1][branch]-1,domain='uniform cubic Pd reference, full bulk free-face thickness; not a grainwise continuum-convergence bound'))
    table(out/f'thickness_r{rung}.csv',branches)
    if rung==2:
        spatial=_space()
        phase=[];numerical={angle:spatial.plane_branches(c,.01373147/9,angle) for angle in (0.,np.pi/4)}
        for row in branches:
            hz=float(numerical[row['direction_rad']][row['branch_order']])
            phase.append(dict(direction_rad=row['direction_rad'],branch_order=row['branch_order'],wavelength_m=.01373147/9,Q4_P4_frequency_Hz=hz,continuum_P8_thickness_frequency_Hz=row['P8_frequency_Hz'],relative_uniform_reference_phase_velocity_error=hz/row['P8_frequency_Hz']-1,scope='uniform-reference dispersion, not a grained-foil continuum error bound'))
        table(out/'phase_r2.csv',phase)

def clips(rung,losses=('retained','Cu_proxy')):
    """Actual signed modal trajectories, not synthetic moving pictures."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation,PillowWriter
    from PIL import Image
    import scipy
    c,l,r,f=old.accepted();out=dest(rung)
    if rung==2:spatial=_space()
    receipt_path=out/f'clips_r{rung}.json'
    receipts=json.loads(receipt_path.read_text()) if receipt_path.exists() else []
    for loss in losses:
        m,z=run_arm(rung,2,.35,loss,True,[0]);w=m['omega']
        if rung==0:full=dict(rung=0,v=f['eigenvectors'],degree=(5,5,2))
        else:
            data=np.load((HERE.parent/'rung1' if rung==1 else out)/'model.npz',allow_pickle=False)
            vectors=data['eigenvectors']
            full=dict(rung=rung,v=np.ascontiguousarray(vectors) if rung==2 else vectors,degree=tuple(data['degree']) if rung==1 else (spatial.PX,spatial.PX,spatial.PZ))
        axis=np.linspace(0,c.DIMS[0],65);xx,yy=np.meshgrid(axis,axis,indexing='ij')
        points=np.column_stack([xx.ravel(),yy.ravel(),np.full(xx.size,30e-6)])
        velocity,strain=maps(c,full,points);projectors=np.array([velocity[:,2]/w,strain[:,0]/w])
        frames=np.arange(0,SAMPLE,10);components=z['x'][frames,0]
        components=np.concatenate([components,components.sum(axis=1,keepdims=True)],axis=1)
        fields=np.einsum('opm,tcm->otcp',projectors,components,optimize=True).reshape(2,len(frames),4,65,65)
        scales=np.max(abs(fields),axis=(1,3,4));scales=np.maximum(scales,np.finfo(float).tiny)
        fig,axes=plt.subplots(2,4,figsize=(12,5),dpi=80,layout='constrained')
        names=('random','inward momentum','thermal eigenstrain','combined');images=[]
        for oi,observable in enumerate(('u_z [m]','epsilon_xx')):
            row=[]
            for ci,name in enumerate(names):
                image=axes[oi,ci].imshow(fields[oi,0,ci].T,origin='lower',extent=(0,c.DIMS[0]*1e3,0,c.DIMS[1]*1e3),cmap='RdBu_r',vmin=-scales[oi,ci],vmax=scales[oi,ci],interpolation='nearest')
                axes[oi,ci].set_title(name+chr(10)+f'{observable}, sampled clip +/-{scales[oi,ci]:.2e}',fontsize=8);axes[oi,ci].set_xlabel('x [mm]',fontsize=8);axes[oi,ci].set_ylabel('y [mm]',fontsize=8);axes[oi,ci].tick_params(labelsize=7);row.append(image)
            images.append(row)
        title=fig.suptitle('',fontsize=9)
        def animate(index):
            for oi in range(2):
                for ci in range(4):images[oi][ci].set_data(fields[oi,index,ci].T)
            title.set_text(f'rung {rung}, hot spots + bursts, beta .35, {loss}, record 0; t={frames[index]*DT*1e6:.2f} us'+chr(10)+'Dynamic field at 30 um; separate fixed physical scales per panel; zero is white')
            return [image for row in images for image in row]+[title]
        name=f'clip_r{rung}_hot_bursts_b035_{loss}.png'
        FuncAnimation(fig,animate,frames=len(frames),interval=50,blit=False).save(out/name,writer=PillowWriter(fps=20))
        animate(len(frames)//2);fig.savefig(out/name.replace('.png','_preview.png'));plt.close(fig)
        with Image.open(out/name) as image:
            actual=image.n_frames
            image.seek(actual-1);image.load()
        if actual!=len(frames):raise AssertionError('clip frames collapsed or missing')
        batch=out/f'r{rung}_hot_bursts_b0.35_{loss}_combined_records0-15_records.csv'
        match=None
        if batch.exists():
            row=next(row for row in csv.DictReader(batch.open()) if row['record']=='0')
            match=row['draw_sha256']==z['digests'][0]
        receipts=[row for row in receipts if row['path']!=str((out/name).relative_to(ROOT))]
        receipts.append(dict(path=str((out/name).relative_to(ROOT)),sha256=sha(out/name),bytes=(out/name).stat().st_size,frames=actual,physical_frame_step_s=10*DT,physical_first_time_s=0.,physical_last_time_s=float(frames[-1]*DT),playback_fps=20,source='actual run_arm signed dynamic modal coefficients',record=0,parents=int(z['counts']['all'][0]),draw_sha256=z['digests'][0],ensemble_mark_digest_equal=match,sampled_clip_scales=scales.tolist(),decoded_final_frame=True,media_type='image/apng',spatial_map_spacing_m=float(axis[1]-axis[0]),limitations='representative solo record; finite-band field, not atomic impact resolution'))
        digest=hashlib.sha256()
        for part in (z['x'],z['y']):digest.update(memoryview(part).cast('B'))
        receipts[-1].update(dynamic_fields_sha256=digest.hexdigest(),numerical_runtime=dict(numpy=np.__version__,scipy=scipy.__version__))
    dump(out/f'clips_r{rung}.json',receipts)

def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','summarize','transfer','static','proof','clips','package']);p.add_argument('--rung',type=int,default=1);p.add_argument('--beam',type=int,default=0);p.add_argument('--beta',type=float,default=0.);p.add_argument('--loss',default='retained');p.add_argument('--random-only',action='store_true');p.add_argument('--first',type=int,default=0);p.add_argument('--count',type=int,default=NREC);a=p.parse_args()
    if a.action=='run':run(a.rung,a.beam,a.beta,a.loss,not a.random_only,a.first,a.count)
    else:globals()[a.action](a.rung)
if __name__=='__main__':main()
