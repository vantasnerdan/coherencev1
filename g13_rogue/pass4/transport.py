"""Pass-4 whole-foil conditional Green transport. Run from repository root.
Only this additive directory is written. Existing resolved producers are not rebuilt.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys
import zipfile
import numpy as np
from numpy.polynomial.legendre import leggauss, legval, legder
from scipy.linalg import eigh
from scipy.spatial.transform import Rotation

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'g13_rogue/pass2'))
sys.path.insert(0,str(ROOT/'g13_rogue/pass2/replication'))
import coarse as accepted
from replicate import Portable
QE=1.602176634e-19
MD=3.3435837768e-27
RHO=12007.
CP=25.89/.1064
KTH=71.8
ALPHA=11.77e-6
DTH=KTH/(RHO*CP)
BETA=ALPHA*(accepted.C11+2*accepted.C12)
RATE=138e-6/QE
L=float(accepted.DIMS[0]); H=.00025
FRACTION=np.array([.85,.15])
PARENT_P=np.sqrt(2*MD*5000*QE)*np.array([1.,np.sqrt(2)])
BANDS=np.array([2.5e5*10**j for j in range(8)])
DEPTHS=np.array([63e-9,30e-6,50e-6,200e-6])
ROT=Rotation.random(64,random_state=173).as_matrix()
CS=np.array([accepted.cubic(r) for r in ROT])
PAIRS=((0,0),(1,1),(2,2),(1,2),(0,2),(0,1))
CT=np.empty((64,3,3,3,3))
for a,(i,j) in enumerate(PAIRS):
    for b,(k,l) in enumerate(PAIRS):
        for ii,jj in ((i,j),(j,i)):
            for kk,ll in ((k,l),(l,k)):
                CT[:,ii,jj,kk,ll]=CS[:,a,b]
FCAP=np.array([3.4e12,4.7e12,6.8e12])
EXPECTED={'model.npz':'59697212d397a6233f7aee9a4d02bd07ac773c8d3d380d50eaba116dffbcef8d','signed_r2.npz':'55d3d2b85f66d12d77bee585e1dcb385e7e4f307fb2372646f85426ec7849548'}


def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def dump(name,x):
    def scalar(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    (HERE/name).write_text(json.dumps(x,indent=2,allow_nan=False,default=scalar)+'\n')


def table(name,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (HERE/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def reference(path):
    path=Path(path)
    for name,digest in EXPECTED.items():
        if sha(path/name)!=digest:raise ValueError('Wrong exact rung2 '+name)
    m=dict(np.load(path/'signed_r2.npz',allow_pickle=False))
    with np.load(path/'model.npz',allow_pickle=False) as d:
        if not np.array_equal(d['omega'],m['omega']):raise ValueError('Rung2 mode epoch mismatch')
    with zipfile.ZipFile(path/'model.npz') as archive:
        with archive.open('eigenvectors.npy') as header:
            if np.lib.format.read_magic(header)!=(1,0):raise ValueError('Unexpected frozen basis header format')
            shape=list(np.lib.format.read_array_header_1_0(header)[0])
    dump('dependencies.json',dict(base='650a2ebd57c35f371ecd23c4db975ebd6e3cc55a',accepted_dependency_commit='2afd47ffb85df8e373abb30dde577620e450c251',files=[dict(path='g13_rogue/pass3/rung2/'+k,sha256=v) for k,v in EXPECTED.items()],access='Use exact accepted dependency commit 2afd47ffb85df8e373abb30dde577620e450c251 and its rung2 evidence recipe: draft 460035582 evidence asset 798321285, SHA e61ab2af46343f68cf682cab4ace136e6fdc8fbddbf5f37148605d87e5845638. Extract into repository root or pass --reference-dir to a verified read-only extraction. No dependency bytes are duplicated here.',numeric_basis_shape=shape))
    return m


def loss(f,vg,law):
    f=np.asarray(f)
    if law=='retained':return 1/50e-6+np.asarray(vg)/.2
    if law=='Cu_extrapolated':return np.asarray(vg)*np.log(10)/10*32*f/1e6
    if law=='Al_transfer':return f/(34e-12*240e9)
    raise ValueError(law)


def overlap(m):
    """Matched one-parent band-limited impulse, exact spatial energy quadratic."""
    src=int(np.argmin(np.sum((m['xy']-L/2)**2,axis=1)))
    p=PARENT_P[0]; a=p*m['pressure'][:,src]; w=m['omega']
    work=float(.5*a@a)
    t=np.linspace(0,100e-6,2001)
    rows=[]; arrays={'time_s':t}
    sys.path.insert(0,str(ROOT/'g13_rogue/pass3/rung1b'))
    import signed as accepted_signed
    for law in ('retained','Cu_extrapolated'):
        gamma=loss(w/(2*np.pi),m['speed'] if law=='retained' else 4460.,law)
        wd=np.sqrt(w*w-gamma*gamma/4)
        q=a*np.exp(-t[:,None]*gamma/2)*np.sin(t[:,None]*wd)/wd
        v=a*np.exp(-t[:,None]*gamma/2)*(np.cos(t[:,None]*wd)-gamma/(2*wd)*np.sin(t[:,None]*wd))
        x=q*w
        aa,ab,ba,bb=accepted_signed.transition(w,gamma)
        qa=np.zeros(len(w));va=a.copy()
        accepted_q=[qa.copy()];accepted_v=[va.copy()]
        for _ in range(len(t)-1):
            qa,va=aa*qa+ab*va,ba*qa+bb*va
            accepted_q.append(qa.copy());accepted_v.append(va.copy())
        accepted_q=np.array(accepted_q);accepted_v=np.array(accepted_v)
        # An explicit parent-phase shuffle, same amplitudes and quadratic consumer.
        ph=2*np.pi*((np.arange(len(w))*.6180339887498949)%1)
        xsh=a*np.exp(-t[:,None]*gamma/2)*np.sin(t[:,None]*wd+ph)*w/wd
        for depth in (0,1,2):
            j=depth*9+4;z=m['points'][j,2]
            strain=x@m['strain'][j].T
            elast=.5*np.einsum('ti,ij,tj->t',strain,m['cs'][j],strain)
            kin=.5*RHO*np.sum((v@m['velocity'][j].T)**2,axis=1)
            accepted_strain=(accepted_q*w)@m['strain'][j].T
            accepted_energy=.5*np.einsum('ti,ij,tj->t',accepted_strain,m['cs'][j],accepted_strain)+.5*RHO*np.sum((accepted_v@m['velocity'][j].T)**2,axis=1)
            arrays[f'{law}_z{depth}_accepted_transition_total']=accepted_energy
            reconstruction_relative=float(np.max(abs(accepted_energy-elast-kin))/max(float((elast+kin).max()),1e-300))
            rows.append(dict(comparison='matched exact guided Green versus accepted transition',loss=law,depth_m=z,impulse_N_s=p,launch_work_J=work,accepted_mean_total_J_m3=float(accepted_energy.mean()),green_mean_total_J_m3=float((elast+kin).mean()),accepted_sd_total_J_m3=float(accepted_energy.std()),green_sd_total_J_m3=float((elast+kin).std()),max_relative_reconstruction_error=reconstruction_relative,first_sample_phase='q=0,v=normal impulse; shared 50-ns sampling',reflecting_boundaries='same frozen whole-foil six-free-face operator',decision='retain exact guided kernel, same receiver metric/normalization'))
            shuffled=.5*np.einsum('ti,ij,tj->t',xsh@m['strain'][j].T,m['cs'][j],xsh@m['strain'][j].T)
            d=np.linalg.norm(np.r_[m['xy'][src]-L/2,z]);vg=4460.
            # Same total band work; finite Gaussian packet and reflected free-face images.
            tau=1/(2*np.pi*245e3)
            bulk=np.zeros_like(t)
            for n in range(-4,5):
                for sign in (-1,1):
                    dist=np.linalg.norm(np.r_[m['xy'][src]-L/2,2*n*H+sign*z])
                    tt=dist/vg
                    bulk+=work/(4*np.pi*(dist*dist+(20e-9)**2)*vg*np.sqrt(2*np.pi)*tau)*np.exp(-loss(1e5,vg,law)*tt-.5*((t-tt)/tau)**2)
            arrays[f'{law}_z{depth}_resolved_elastic']=elast
            arrays[f'{law}_z{depth}_resolved_total']=elast+kin
            arrays[f'{law}_z{depth}_bulk_total']=bulk
            rows.append(dict(loss=law,depth_m=z,source_x_m=m['xy'][src,0],source_y_m=m['xy'][src,1],impulse_N_s=p,launch_work_J=work,resolved_mean_total_J_m3=float((elast+kin).mean()),bulk_mean_total_J_m3=float(bulk.mean()),bulk_over_resolved_mean=float(bulk.mean()/(elast+kin).mean()),resolved_elastic_sd_J_m3=float(elast.std()),shuffled_elastic_sd_J_m3=float(shuffled.std()),resolved_peak_s=float(t[np.argmax(elast+kin)]),bulk_first_transit_s=d/vg,bulk_thickness_roundtrip_s=2*H/vg,resolved_slowest_group_m_s=float(m['speed'].min()),resolved_fastest_group_m_s=float(m['speed'].max()),resolved_temporal_integral_J_s_m3=float(np.trapezoid(elast+kin,t)),bulk_temporal_integral_J_s_m3=float(np.trapezoid(bulk,t)),decision='retain exact resolved low-band kernel; bulk images do not preserve guided modal receiver quadratics'))
    table('overlap.csv',rows)
    np.savez_compressed(HERE/'overlap.npz',**arrays)
    dump('overlap_decision.json',dict(common_band_Hz=[5000,250000],normalized_forcing='one D+ normal surface impulse; exact accepted source column',matched_comparison='exact guided Green versus unmodified accepted signed.transition; same phase, dt, free-face operator and spatial quadratic',bulk_illustration='finite front/rear bulk-image sum only, not a matched all-face or converged resonance calculation; its disagreement is not a universal transport no-go',decision='retain exact resolved kernel below 250 kHz and free-face guided thickness branches; do not replace by uncontrolled bulk kernel',ratios=[r['bulk_over_resolved_mean'] for r in rows if 'bulk_over_resolved_mean' in r]))
    print(json.dumps(dict(overlap='completed',guided_max_relative_error=max(r['max_relative_reconstruction_error'] for r in rows if 'max_relative_reconstruction_error' in r),illustrative_bulk_ratios=[r['bulk_over_resolved_mean'] for r in rows if 'bulk_over_resolved_mean' in r])),flush=True)


def track_data(m):
    _,st=accepted.profiles_and_stopping();ts=st['PdD',2.5]['tracks'];out=[]
    for s,tr in enumerate(ts):
        z=np.array([0.]+[r['depth_A']*1e-10 for r in tr['timeline']]);time=np.array([0.]+[r['time_fs']*1e-15 for r in tr['timeline']])
        mid=(z[1:]+z[:-1])/2;tm=(time[1:]+time[:-1])/2
        heat=[]
        for col in ('electron_eV','host_recoil_eV','D_recoil_eV'):
            a=np.diff([0.]+[r[col] for r in tr['timeline']]);heat.append(a/a.sum())
        en=np.array([tr['initial_eV']]+[r['energy_eV'] for r in tr['timeline']])
        moment=-np.diff(np.sqrt(np.maximum(en,0)));moment[-1]+=np.sqrt(en[-1]);moment/=moment.sum()
        out.append(dict(z=mid,t=tm,heat=heat,momentum=moment,terminal_z=tr['terminal_depth_A'],terminal_t=tr['terminal_time_fs']*1e-15,energy=m['heatenergy'][s]))
    return out


def christoffel(n,grain=27):
    n=np.asarray(n)
    g=np.einsum('ijkl,...j,...l->...ik',CT[grain],n,n)/RHO
    val,e=np.linalg.eigh(g)
    return np.sqrt(np.maximum(val,0)),e


def dispersion(f,c):
    ratio=np.asarray(f)/FCAP
    k=2*np.pi*FCAP/c*np.arcsin(ratio)
    vg=c*np.sqrt(1-ratio*ratio)
    return k,vg


def plate_matrices(k,theta,grain,degree=14):
    z,ww=leggauss(degree+3);ww*=H/2
    ph=np.array([legval(z,np.eye(degree+1)[j])*np.sqrt((2*j+1)/(RHO*H)) for j in range(degree+1)]).T
    dz=np.array([legval(z,legder(np.eye(degree+1)[j]))*2/H*np.sqrt((2*j+1)/(RHO*H)) for j in range(degree+1)]).T
    b=np.zeros((len(z),6,3*(degree+1)),complex);n=degree+1
    kx=k*np.cos(theta);ky=k*np.sin(theta)
    b[:,0,:n]=1j*kx*ph;b[:,1,n:2*n]=1j*ky*ph;b[:,2,2*n:]=dz
    b[:,3,n:2*n]=dz;b[:,3,2*n:]=1j*ky*ph
    b[:,4,:n]=dz;b[:,4,2*n:]=1j*kx*ph
    b[:,5,:n]=1j*ky*ph;b[:,5,n:2*n]=1j*kx*ph
    kk=np.einsum('zai,ab,zbj,z->ij',b.conj(),CS[grain],b,ww,optimize=True)
    val,v=eigh(kk,subset_by_index=[0,min(17,len(kk)-1)])
    omega=np.sqrt(np.maximum(val,0))
    return omega,v,ph,dz,z,ww,b


def bulk_spectrum(tracks,lo,hi,nfreq=12):
    """The sole bulk DOS/source implementation, also used by selected refinement."""
    shapes=((20e-9,.5e-12),(60e-9,2e-12))
    a0,aw=leggauss(32);u=(a0+1)/2
    mu=np.r_[-1+u*u,1-u*u];mw=np.r_[aw*u,aw*u]
    theta=2*np.pi*(np.arange(12)+.5)/12
    directions=np.array([[np.sqrt(1-u*u)*np.cos(a),np.sqrt(1-u*u)*np.sin(a),u] for u in mu for a in theta])
    dw=np.repeat(mw,12)*2*np.pi/12;c,ev=christoffel(directions)
    logf,lw=leggauss(nfreq);freq=np.exp((logf+1)/2*np.log(hi/lo)+np.log(lo));wf=lw/2*np.log(hi/lo)*freq
    output=np.zeros((nfreq,2,2,3,3))
    for j,(f,df) in enumerate(zip(freq,wf)):
        w=2*np.pi*f;k,vg=dispersion(f,c)
        dos=k*k/(2*RHO*(2*np.pi)**3*vg)*dw[:,None]*(2*np.pi*df)
        for s,tr in enumerate(tracks):
            for ish,(a,tau) in enumerate(shapes):
                ff=np.exp(-.5*(a*a*k*k*(1-directions[:,2,None]**2)+(w*tau)**2))
                phase=np.exp(1j*w*tr['t'][None,None,:]-1j*k[:,:,None]*directions[:,2,None,None]*tr['z'])
                fp=phase@tr['momentum']
                fq=sum((phase@q)*tr['energy'][i] for i,q in enumerate(tr['heat']))
                fq+=tr['energy'][3]*np.exp(1j*w*tr['terminal_t']-1j*k*directions[:,2,None]*tr['terminal_z'])
                ap=PARENT_P[s]*ev[:,2,:]*fp*ff
                longitudinal=np.einsum('ni,nij->nj',directions,ev)
                thermal_decay=DTH*np.minimum(k*k,(2*np.pi/20e-9)**2)
                at=-1j*BETA/(RHO*CP)*k*longitudinal*fq/(thermal_decay-1j*w)*ff
                terms=np.stack([abs(ap)**2,abs(at)**2,2*np.real(ap*at.conj())],axis=-1)
                output[j,ish,s]=np.sum(dos[:,:,None]*terms,axis=0)
    return freq,output


def built_launch(m,nk=120,nang=8):
    tracks=track_data(m);energies=np.zeros((2,2,7,3,3))
    # radius, duration are assumed, not CSDA track observables.
    shapes=((20e-9,.5e-12),(60e-9,2e-12))
    # Real finite-thickness branches in the first and transition band; eight
    # source-grain/direction pairs span the full seeded foil, not separate patches.
    grainbank=(9,12,18,23,40,45,51,54)
    records=[]
    ks=np.geomspace(.1,1.2e5,nk);dk=np.gradient(ks)
    for ia in range(nang):
        theta=2*np.pi*(ia+.5)/nang;grain=grainbank[ia%8]
        for ki,k in enumerate(ks):
            w,v,ph,dz,z,ww,b=plate_matrices(k,theta,grain)
            freq=w/(2*np.pi)
            surf=np.array([legval(-1,np.eye(15)[j])*np.sqrt((2*j+1)/(RHO*H)) for j in range(15)])
            up=surf@v[30:]
            div=(b[:,:3].sum(axis=1)@v).conj()
            modes=np.arange(64)*np.pi/H
            integ=(div.T*ww)@np.cos(np.outer((z+1)*H/2,modes))
            for s,tr in enumerate(tracks):
                hc=m['heatcos'][s]
                hproj=integ*((hc.T*tr['energy']).sum(axis=1)*np.r_[1.,np.full(63,2.)]/H)[None,:]
                th=BETA/(RHO*CP)*np.sum(hproj/(DTH*(k*k+modes*modes)[None,:]-1j*w[:,None]),axis=1)
                for ish,(a,tau) in enumerate(shapes):
                    ff=np.exp(-.5*(a*k)**2-.5*(w*tau)**2)
                    ap=PARENT_P[s]*up.conj()*ff
                    at=th*ff
                    # Free-face eigenstress is projected by divergence, not a net pressure.
                    terms=np.stack([abs(ap)**2,abs(at)**2,2*np.real(ap*at.conj())],axis=-1)
                    de=.5*k*dk[ki]/(2*np.pi*nang)*terms
                    for j,f in enumerate(freq):
                        band=int(np.searchsorted(BANDS,f,side='right')-1)
                        if 0<=band<2 and f<1e7:
                            # lowest branch A0, next two guided, higher thickness families.
                            branch=min(j,2)
                            energies[ish,s,band,branch]+=de[j]
            if ia==0:
                for j,f in enumerate(freq):
                    if 1e4<f<1e7:
                        records.append(dict(k_parallel_m_1=k,branch=j,frequency_Hz=f,kh=k*H,grain=grain,theta_rad=theta,normal_surface_participation_per_sqrt_kg_m2=float(abs(up[j]))))
    table('guided_branches.csv',records)
    # Bulk continuum only above 10 MHz. Coarse all-decade quadrature first.
    for band in range(1,7):
        freq,terms=bulk_spectrum(tracks,max(BANDS[band],1e7),BANDS[band+1])
        energies[:,:,band]+=terms.sum(axis=0)
    rows=[]
    total=energies.sum(axis=-1)
    if np.any(total < -1e-35):raise ValueError('Negative built acoustic energy')
    for ish,(a,tau) in enumerate(shapes):
        for s in range(2):
            for band in range(7):
                for branch in range(3):
                    p,t,cross=energies[ish,s,band,branch]
                    rows.append(dict(scenario='built',shape=ish,species=('D+','D2+')[s],band=band,branch=branch,lo_Hz=BANDS[band],hi_Hz=BANDS[band+1],radius_m=a,assumed_sigma_t_s=tau,momentum_J=p,thermoelastic_J=t,cross_J=cross,total_J=p+t+cross,model='free-face guided thickness spectrum below 10 MHz; dispersive crystal bulk DOS above'))
    np.savez_compressed(HERE/'launch.npz',built_terms=energies,built=total)
    table('launch.csv',rows)
    dump('launch_inputs.json',dict(shapes=shapes,thermal_diffusivity_m2_s=DTH,thermal_high_k='D k^2 capped at wavelength 20 nm; conditional mesoscopic relaxation continuation, not atomic Fourier measurement',caps_Hz=FCAP,cap_model='omega=omega_cap sin(c0 k/omega_cap); approximate interpolation to visually read pure-Pd branches, not a fitted PdD dispersion',grain_samples=grainbank,plate_degree=14,plate_k_nodes=nk,angular_nodes=nang,bulk_frequency_nodes_per_decade=12,bulk_mu_nodes=64,bulk_azimuth_nodes=12,track=[dict(terminal_m=t['terminal_z'],transit_s=t['terminal_t'],energy_J=t['energy']) for t in tracks]))
    print(json.dumps(dict(launch='completed',total_J=total.sum(axis=(2,3)).tolist())),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['overlap','launch','receivers'])
    p.add_argument('--reference-dir',default=str(ROOT/'g13_rogue/pass3/rung2'))
    p.add_argument('--beam',choices=['even_steady','hot_steady','hot_bursts'])
    a=p.parse_args();m=reference(a.reference_dir)
    if a.action=='overlap':overlap(m)
    elif a.action=='launch':built_launch(m)
    elif a.action=='receivers':
        import receiver
        receiver.run(m,a.reference_dir,a.beam)

if __name__=='__main__':main()
