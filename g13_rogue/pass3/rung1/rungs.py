"""Unpublished G13 whole-foil support and in-band refinement; run at repo root."""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import io
import json
import sys
import tarfile
import time
import numpy as np
from numpy.polynomial.legendre import leggauss, legval, legder
from scipy.linalg import eigh
from scipy.optimize import brentq
from scipy.special import spherical_jn, gammainc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = 'b0a67c08661ec765c3916f1dad02c47f07cf3bdd'
PREFIX = 'g13_rogue/pass3/rung1/'
DT = 5e-8
BURN = 4704
SAMPLE = 2000
RECORDS = [0, 1]
INPUT = ROOT/'g13_rogue/pass2/replication/inputs.npz'
EXPECTED = '4ca83f716d56cf76d6c416571fdd01f3d0bc7e70b9b39a2eb4adf485e6f7a9b1'


def dump(name, obj):
    (HERE/name).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')


def table(name, rows):
    if rows:
        with (HERE/name).open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def accepted():
    sys.path.insert(0, str(ROOT/'g13_rogue/pass2'))
    sys.path.insert(0, str(ROOT/'g13_rogue/pass2/replication'))
    import coarse
    import live
    import replicate
    if sha(INPUT) != EXPECTED:
        raise ValueError('Frozen accepted numeric inputs changed')
    ledger=json.loads((INPUT.parent/'source_ledger.json').read_text())
    dependencies=dict(ledger['accepted_sources'])
    manifest=json.loads((INPUT.parent/'manifest.json').read_text())
    for f in manifest['files']:
        if f['path'].endswith(('/replicate.py','/portable.c')):
            dependencies[f['path']]=f['sha256']
    for path,expected in dependencies.items():
        if sha(ROOT/path)!=expected:
            raise ValueError('Accepted dependency changed: '+path)
    return coarse, live, replicate, dict(np.load(INPUT, allow_pickle=False))


def capture(degree, wavelength, length):
    """Real sin/cos plane waves, physical-volume L2; sampled phase/direction min."""
    theta = np.linspace(0, np.pi/2, 37)
    phase = np.array([0., np.pi/4, np.pi/2])
    q = np.pi*length/wavelength
    n = np.arange(degree+1)
    ax = np.sqrt(2*n+1)*1j**n*spherical_jn(n[None, :], q*np.cos(theta[:, None]))
    ay = np.sqrt(2*n+1)*1j**n*spherical_jn(n[None, :], q*np.sin(theta[:, None]))
    coeff = ax[:, :, None]*ay[:, None, :]
    real = np.real(np.exp(1j*phase[:, None, None, None])*coeff[None])
    norm = .5*(1+np.cos(2*phase[:, None])*np.sinc(2*q*np.cos(theta)/np.pi)*np.sinc(2*q*np.sin(theta)/np.pi))
    scores = np.sum(real*real, axis=(2, 3))/norm
    j = np.unravel_index(np.argmin(scores), scores.shape)
    return float(scores[j]), float(theta[j[1]]), float(phase[j[0]])


def support(length):
    lam0 = brentq(lambda lam:capture(5, lam, length)[0]-.95, .010, .020, xtol=1e-13)
    degree = next(p for p in range(6, 40) if capture(p, lam0/3, length)[0] >= .95)
    rows = []
    for p in (5, degree):
        for lam in sorted(set([.00255, .003519192520091718, .005, .007, .010, .012, .015, .020, lam0, lam0/3, lam0/9])):
            value, angle, phase = capture(p, lam, length)
            rows.append(dict(degree_x=p, degree_y=p, wavelength_m=lam, min_mass_capture=value,
                             worst_sampled_direction_rad=angle, worst_sampled_phase_rad=phase,
                             passes_95percent=value >= .95-1e-12))
    table('projection_support.csv', rows)
    return lam0, degree


def one_axis(p, length, cells):
    x, w = leggauss(p+1)
    result = np.empty((cells, 2, 2, p+1, p+1))
    for cell in range(cells):
        z = (cell+(x+1)/2)*length/cells
        xi = 2*z/length-1
        v, d = [], []
        for n in range(p+1):
            a = np.zeros(n+1); a[n] = 1
            norm = np.sqrt((2*n+1)/length)
            v.append(norm*legval(xi, a))
            d.append(norm*legval(xi, legder(a))*2/length)
        vv = [np.array(v).T, np.array(d).T]
        for a in range(2):
            for b in range(2):
                result[cell, a, b] = vv[a].T@((w*length/(2*cells))[:, None]*vv[b])
    return result


def stiffness(c, degree):
    """Exact per-grain polynomial integration, tensor contractions, natural free BC."""
    from scipy.spatial.transform import Rotation
    p, z = degree[0], degree[2]
    ax = one_axis(p, c.DIMS[0], 8)
    ay = one_axis(p, c.DIMS[1], 8)
    az = one_axis(z, c.DIMS[2], 1)[0]
    tensors = np.array([c.cubic(r) for r in Rotation.random(64, random_state=173).as_matrix()]).reshape(8, 8, 6, 6)
    smap = np.array([[0, 5, 4], [5, 1, 3], [4, 3, 2]])
    n = (p+1)**2*(z+1)
    k = np.empty((3*n, 3*n))
    for i in range(3):
        for j in range(i, 3):
            block = np.zeros((n, n))
            for a in range(3):
                for b in range(3):
                    xx = ax[:, int(a == 0), int(b == 0)]
                    yy = ay[:, int(a == 1), int(b == 1)]
                    cc = tensors[:, :, smap[i, a], smap[j, b]]
                    xy = np.einsum('ab,aij,bkl->ikjl', cc, xx, yy, optimize=True).reshape((p+1)**2, (p+1)**2)
                    block += np.kron(xy, az[int(a == 2), int(b == 2)])/c.RHO
            k[i*n:(i+1)*n, j*n:(j+1)*n] = block
            if i != j:
                k[j*n:(j+1)*n, i*n:(i+1)*n] = block.T
    return (k+k.T)/2


def embed_indices(p, z):
    return np.arange(3*(p+1)**2*(z+1)).reshape(3, p+1, p+1, z+1)[:, :6, :6, :3].ravel()


def gradient_metrics(c, v, degree):
    nx, ny, nz = [p+1 for p in degree]
    vv = v.reshape(3, nx, ny, nz, -1)
    gx = one_axis(nx-1, c.DIMS[0], 1)[0, 1, 1]
    gy = one_axis(ny-1, c.DIMS[1], 1)[0, 1, 1]
    k2 = np.einsum('caijm,ab,cbijm->m', vv, gx, vv)+np.einsum('caijm,ib,cabjm->m', vv, gy, vv)
    mass = np.sum(vv*vv, axis=(1, 2, 3))
    return np.sqrt(k2), mass


def prepare():
    start = time.monotonic()
    c, l, r, frozen = accepted()
    lam0, p = support(c.DIMS[0])
    degree = (p, p, 4)
    k = stiffness(c, degree)
    indices = embed_indices(p, 4)
    oldk = k[np.ix_(indices, indices)]
    v0, w0 = frozen['eigenvectors'], frozen['omega']
    nested_residual = np.linalg.norm(oldk@v0-v0*w0**2)/np.linalg.norm(v0*w0**2)
    if nested_residual > 1e-8:
        raise AssertionError('Tensor assembly does not recover accepted Ritz restriction')
    val, v = eigh(k, subset_by_value=((2*np.pi*5000)**2, (2*np.pi*250000)**2), driver='evr', check_finite=False)
    omega = np.sqrt(val)
    residual = np.linalg.norm(k@v-v*val)/np.linalg.norm(v*val)
    kk0, mass0 = gradient_metrics(c, v0, (5, 5, 2))
    kk, mass = gradient_metrics(c, v, degree)
    overlap = v[indices].T@v0
    # One representative mode at the screened long-wave scale, not a continuum bound.
    choice = int(np.argmin(abs(2*np.pi/kk0-lam0)))
    match = int(np.argmax(overlap[:, choice]**2))
    rows = []
    for rung, om, kval, masses in ((0, w0, kk0, mass0), (1, omega, kk, mass)):
        for m in range(len(om)):
            rows.append(dict(rung=rung, mode=m, frequency_Hz=om[m]/(2*np.pi), full_vector_k_per_m=kval[m],
                             full_vector_gradient_lambda_m=2*np.pi/kval[m], mass_x=masses[0,m], mass_y=masses[1,m], mass_z=masses[2,m], kh=kval[m]*c.DIMS[2]))
    table('mode_support.csv', rows)
    np.savez_compressed(HERE/'model.npz', omega=omega, eigenvectors=v, degree=degree)
    summary = dict(base=BASE, frozen_inputs_sha256=EXPECTED, criterion='minimum sampled real-plane-wave physical-volume L2 projection mass >= .95; 37 directions, three phases',
        lambda0_screen_m=lam0, lambda1_target_m=lam0/3, lambda2_target_m=lam0/9, rung1_degree=degree,
        lambda1_capture=capture(p,lam0/3,c.DIMS[0])[0], nominal_pass2_fmax_Hz=250000, actual_old_max_f_Hz=max(w0)/(2*np.pi),
        old_full_gradient_min_lambda_m=min(2*np.pi/kk0), nested_operator_residual=nested_residual, enriched_operator_residual=residual,
        old_modes=len(w0), enriched_modes=len(omega), ndof=len(k), elapsed_s=time.monotonic()-start,
        representative_dispersion=dict(old_mode=choice, enriched_mode=match, old_frequency_Hz=w0[choice]/(2*np.pi), enriched_frequency_Hz=omega[match]/(2*np.pi),
            squared_mass_overlap=overlap[match,choice]**2, relative_frequency_change=omega[match]/w0[choice]-1,
            old_full_gradient_lambda_m=2*np.pi/kk0[choice]))
    dump('support.json', summary)
    print(json.dumps(summary), flush=True)


def linear_box(c, v, omega, degree, side, z0, z1, center):
    av, dv = [], []
    for axis, (lo, hi) in enumerate(((center[0]-side/2,center[0]+side/2),
                                    (center[1]-side/2,center[1]+side/2),(z0,z1))):
        x, w = leggauss(degree[axis]+1)
        pts = (lo+hi)/2+(hi-lo)*x/2
        vals, ders = [], []
        for n in range(degree[axis]+1):
            a=np.zeros(n+1); a[n]=1
            norm=np.sqrt((2*n+1)/c.DIMS[axis])
            vals.append(float(w@legval(2*pts/c.DIMS[axis]-1,a)/2)*norm)
            ders.append(float(w@legval(2*pts/c.DIMS[axis]-1,legder(a))/2)*norm*2/c.DIMS[axis])
        av.append(np.array(vals)); dv.append(np.array(ders))
    value=np.einsum('i,j,k->ijk',*av).ravel()/np.sqrt(c.RHO)
    deriv=np.einsum('i,j,k->ijk',dv[0],av[1],av[2]).ravel()/np.sqrt(c.RHO)
    n=len(value)
    return deriv@v[:n]/omega, value@v[2*n:]/omega


def modal_speeds(c, v, omega, degree):
    k, mass=gradient_metrics(c,v,degree)
    n=v.shape[0]//3
    # Exact vertical-only diagnostic, not a cutoff for tiny vertical leakage.
    vz=v[2*n:]
    kz,_=gradient_metrics(c,np.vstack([np.zeros_like(vz),np.zeros_like(vz),vz]),degree)
    kz/=np.sqrt(np.maximum(mass[2],np.finfo(float).tiny))
    bulk=(c.C11+2*c.C12)/3
    gv=(c.C11-c.C12+3*c.C44)/5
    gr=5*(c.C11-c.C12)*c.C44/(4*c.C44+3*(c.C11-c.C12))
    shear=(gv+gr)/2
    young=9*bulk*shear/(3*bulk+shear)
    nu=(3*bulk-2*shear)/(2*(3*bulk+shear))
    alpha=np.sqrt(young*c.DIMS[2]**2/(12*(1-nu*nu)*c.RHO))
    axes=[one_axis(p,length,1)[0] for p,length in zip(degree,c.DIMS)]
    divergence=np.zeros(len(omega))
    for i in range(3):
        for j in range(3):
            factors=[axes[a][int(i==a),int(j==a)] for a in range(3)]
            g=np.kron(np.kron(factors[0],factors[1]),factors[2])
            divergence+=np.sum(v[i*n:(i+1)*n]*(g@v[j*n:(j+1)*n]),axis=0)
    bf=bulk/c.RHO*divergence/omega**2
    speed=np.where(bf>.5,np.sqrt((bulk+4*shear/3)/c.RHO),np.sqrt(shear/c.RHO))
    flex=(mass[2]>.5)&(kz*c.DIMS[2]<.5)
    speed[flex]=2*alpha*kz[flex]
    return speed


def regional_energy(c,v,omega,degree,z0,z1):
    axes=[]; weights=[]
    for axis,(lo,hi) in enumerate(((c.DIMS[0]/2-.00125,c.DIMS[0]/2+.00125),
                                  (c.DIMS[1]/2-.00125,c.DIMS[1]/2+.00125),(z0,z1))):
        q,w=leggauss(degree[axis]+1)
        cuts=sorted([lo,hi]+([g for g in np.arange(c.DIMS[axis]/8,c.DIMS[axis],c.DIMS[axis]/8) if lo<g<hi] if axis<2 else []))
        axes.append(np.concatenate([(a+b)/2+(b-a)*q/2 for a,b in zip(cuts[:-1],cuts[1:])]))
        weights.append(np.concatenate([(b-a)*w/2 for a,b in zip(cuts[:-1],cuts[1:])]))
    pts=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3)
    ww=np.einsum('i,j,k->ijk',*weights).ravel()/(.0025**2*(z1-z0))
    a=np.zeros((len(omega),len(omega))); kinetic=np.zeros_like(a)
    for start in range(0,len(pts),128):
        p=pts[start:start+128]; w=ww[start:start+128]
        phi,b=c.basis(p); n=phi.shape[1]
        bm=b@v/omega
        cc=c.tensors(p,c.DIMS[0]/8,'grains')
        z=np.einsum('qji,qjm->qim',np.linalg.cholesky(cc),bm)*np.sqrt(w[:,None,None])
        a+=z.reshape(-1,len(omega)).T@z.reshape(-1,len(omega))
        for i in range(3):
            velocity=(phi@v[i*n:(i+1)*n])*np.sqrt(c.RHO*w[:,None])
            kinetic+=velocity.T@velocity
    return a,kinetic


def consumer_model(rung,c,l,frozen):
    if rung==0:
        degree=(5,5,2); v=frozen['eigenvectors']; omega=frozen['omega']
    else:
        d=np.load(HERE/'model.npz',allow_pickle=False)
        degree=tuple(int(i) for i in d['degree']); v=d['eigenvectors']; omega=d['omega']
    epoch=hashlib.sha256(v.astype('<f8').tobytes()+omega.astype('<f8').tobytes()).hexdigest()
    cache=HERE/f'consumer_r{rung}.npz'
    fields=('source','eps','uz','speed','strain','velocity','cs','points','groups','region_bounds','region_projectors','regional')
    if cache.exists():
        with np.load(cache,allow_pickle=False) as d:
            if str(d['epoch'])!=epoch:
                raise ValueError('Frozen consumer belongs to a different numeric mode epoch')
            result={key:d[key] for key in fields}
            result.update(ref=tuple(d['ref_'+str(i)] for i in range(5)),
                          receiver=tuple(d['receiver_'+str(i)] for i in range(7)))
        return dict(result,rung=rung,degree=degree,v=v,omega=omega)
    original=c.basis
    c.basis=lambda points:original(points,degree=degree)
    try:
        model=dict(v=v,omega=omega,grain=c.DIMS[0]/8,kind='grains')
        xy,groups=l.sites()
        projectors=[linear_box(c,v,omega,degree,.0025,0,50e-6,p) for p in xy]
        eps,uz=np.array([p[0] for p in projectors]),np.array([p[1] for p in projectors])
        _,stop=c.profiles_and_stopping()
        tracks=stop['PdD',2.5]['tracks']
        # Finite 20-nm Gaussian spatial traction; Gauss-Hermite integrates its
        # projection. The 0.5-ps temporal Gaussian is an impulse in this band.
        hx,hw=np.polynomial.hermite.hermgauss(3); hw/=np.sqrt(np.pi)
        source=[]
        for track in tracks:
            columns=[]
            depth=c.deposition_depth(track)
            for point in xy:
                pts=np.array([[*(point+np.sqrt(2)*20e-9*np.array([a,b])),depth] for a in hx for b in hx])
                phi,_=c.basis(pts)
                force=v[2*(v.shape[0]//3):].T@(np.outer(hw,hw).ravel()@phi)
                force*=np.exp(-.5*(omega*.5e-12)**2)
                columns.append(force/np.linalg.norm(force))
            source.append(np.array(columns).T)
        weights=np.array([.85,.15])*frozen['eh']; weights/=weights.sum()
        ref=l.source_reference(model,np.concatenate(source,axis=1),np.repeat(weights/36,36),eps,uz,float(frozen['ecanon'])/50e-6)
        receiver=l.off_beam_receiver(model)
        # Common fixed physical point bank; local energy, not square of mean strain.
        points=np.array([[*(c.DIMS[:2]/2+[dx,dy]),z] for z in (30e-6,50e-6,200e-6) for dx in (-.0035,0,.0035) for dy in (-.0035,0,.0035)])
        phi,b=c.basis(points); n=phi.shape[1]
        strain=b@v/omega
        velocity=np.array([phi@v[i*n:(i+1)*n] for i in range(3)]).transpose(1,0,2)
        cs=c.tensors(points,model['grain'],'grains')
        region_bounds=np.array([(side,z0,z1) for side in (.0025,.01,float(c.DIMS[0]))
                                for z0,z1 in ((0.,50e-6),(25e-6,75e-6),(175e-6,225e-6))])
        region_projectors=np.array([linear_box(c,v,omega,degree,side,z0,z1,c.DIMS[:2]/2)
                                    for side,z0,z1 in region_bounds])
        speed=frozen['group_velocity'] if rung==0 else modal_speeds(c,v,omega,degree)
        result=dict(rung=rung,degree=degree,v=v,omega=omega,source=source,eps=eps,uz=uz,ref=ref,receiver=receiver,
                    speed=speed,strain=strain,velocity=velocity,cs=cs,points=points,groups=groups,
                    region_bounds=region_bounds,region_projectors=region_projectors,
                    regional=np.array([regional_energy(c,v,omega,degree,z0,z1) for z0,z1 in ((0.,50e-6),(25e-6,75e-6),(175e-6,225e-6))]))
        arrays={key:np.asarray(result[key]) for key in fields}
        arrays.update({'ref_'+str(i):np.asarray(a) for i,a in enumerate(ref)})
        arrays.update({'receiver_'+str(i):np.asarray(a) for i,a in enumerate(receiver)})
        np.savez_compressed(cache,epoch=np.array(epoch),**arrays)
        return result
    finally:
        c.basis=original


def freeze():
    c,l,_,frozen=accepted()
    result=[]
    for rung in (0,1):
        m=consumer_model(rung,c,l,frozen)
        error=max(float(np.max(abs(np.linalg.norm(s,axis=0)-1))) for s in m['source'])
        if error>1e-12:
            raise AssertionError('Whole physical-band source normalization failed')
        result.append(dict(rung=rung,modes=len(m['omega']),source_unit_mass_norm_error=error,
                           consumer_path=PREFIX+f'consumer_r{rung}.npz',sha256=sha(HERE/f'consumer_r{rung}.npz')))
    dump('numeric_consumers.json',dict(consumers=result,
        provenance='Reconstructed after propagation from the same frozen modes/routines/runtime; no new wave propagation or cross-runtime equality claim'))
    print(json.dumps(result),flush=True)


def simulate(m,frozen,portable,r,l,spatial,timing,beta,loss):
    omega=m['omega']; nm=len(omega); nr=len(RECORDS)
    ecanon=m['ref'][1]
    # Same finite physical burn/exposure for every arm. No stationarity claim.
    if loss=='retained':
        gamma=1/50e-6+m['speed']/.2
    else:
        # Measured Cu bulk rate transferred as a mode-frequency proxy, not an
        # added grain-scattering law. Mount loss is independently declared zero.
        gamma=4460*np.log(10)/10*32*(omega/(2*np.pi*1e6))
    co,si=np.cos(omega*DT),np.sin(omega*DT)
    damp=np.exp(-gamma*DT/2)
    count,total=r.CountRng(portable,RECORDS),r.TotalRng(portable,RECORDS)
    marks=r.MarkRng(portable,RECORDS,BURN)
    times=np.empty((BURN+SAMPLE,nr),dtype=np.int64)
    for interval,part in enumerate((slice(0,BURN),slice(BURN,BURN+SAMPLE))):
        for rep,record in enumerate(RECORDS):
            n=total.poisson(138e-6/l.QE*(part.stop-part.start)*DT)
            weights=frozen['time_weights'][timing,part,record]
            times[part,rep]=count.multinomial(n,weights/weights.sum())
    x=np.zeros((nr,4,nm)); y=np.zeros_like(x)
    ve=np.broadcast_to(m['ref'][2]**2,(nr,36)).copy()
    vu=np.broadcast_to(m['ref'][3]**2,(nr,36)).copy()
    alpha=-np.expm1(-DT/(5*max(1/gamma)))
    names=('all','A_disturbed','C_any','C_strain_positive','C_strain_negative','C_displacement_positive','C_displacement_negative','A_and_C','A_or_C')
    counts={name:np.zeros(nr,dtype=np.int64) for name in names}
    sx=np.empty((SAMPLE,nr,4,nm)); sy=np.empty_like(sx)
    bmask=np.empty((SAMPLE,nr),dtype=bool)
    masks_saved=np.empty((SAMPLE,nr,3),dtype=float)
    work=np.zeros(nr); incoming=np.zeros(nr); self_energy=np.zeros(nr); damping=np.zeros(nr)
    classwork=np.zeros((nr,3)); capped=np.zeros(nr,dtype=np.int64)
    rp=m['region_projectors'].reshape(-1,nm)
    rv=ecanon*np.einsum('im,mn,in->i',rp,m['ref'][0][:nm,:nm],rp)
    rv=np.broadcast_to(rv,(nr,len(rp))).copy()
    signed_region=np.empty((SAMPLE,nr,len(rp)))
    signed_z=np.empty_like(signed_region)
    initial=None
    fractions=.01/.865*np.array([.1,1.,10.])
    for t in range(BURN+SAMPLE):
        pre=x.sum(axis=1)
        ev=np.sqrt(ecanon)*(pre@m['eps'].T); uv=np.sqrt(ecanon)*(pre@m['uz'].T)
        masks=l.source_masks(ev,uv,np.sqrt(ve),np.sqrt(vu),m['ref'][2],m['ref'][3])
        ve=(1-alpha)*ve+alpha*ev*ev; vu=(1-alpha)*vu+alpha*uv*uv
        region_values=np.sqrt(ecanon)*(pre@rp.T)
        region_z=region_values/np.sqrt(np.maximum(rv,np.finfo(float).tiny))
        rv=(1-alpha)*rv+alpha*region_values**2
        if t>=BURN:
            j=t-BURN
            if initial is None:
                initial=ecanon/2*np.sum(pre**2+y.sum(axis=1)**2,axis=1)
            sx[j]=x; sy[j]=y
            bmask[j]=l.convergence(x,y,m['receiver'],ecanon)[0]
            signed_region[j]=region_values; signed_z[j]=region_z
        parents=np.empty((nr,36,2,2),dtype=np.int64)
        for rep in range(nr):
            site=count.multinomial(int(times[t,rep]),frozen['site_weights'][spatial])
            dp=count.binomial(site,.85)
            species=np.column_stack([dp,site-dp])
            pos=count.binomial(species,.5)
            parents[rep,:,:,0]=pos; parents[rep,:,:,1]=species-pos
        highp=.01+beta*.5*(np.tanh(ev/m['ref'][2])**2+np.tanh(uv/m['ref'][3])**2)
        low=marks.binomial(parents,.25)
        high=marks.binomial(parents-low,highp[:,:,None,None]/.75)
        realized=np.stack([low,parents-low-high,high],axis=-1)
        signed=realized[:,:,:,0]-realized[:,:,:,1]
        amp=np.sum(signed*np.sqrt(2*frozen['eh'][None,None,:,None]*fractions),axis=-1)
        delta=np.zeros_like(y)
        for g in range(4):
            keep=m['groups']==g
            for species in range(2):
                delta[:,g]+=amp[:,keep,species]@m['source'][species][:,keep].T
        oldp=np.sqrt(ecanon)*y.sum(axis=1); dp=delta.sum(axis=1)
        aa=.5*np.sum(dp*dp,axis=1); bb=np.sum(oldp*dp,axis=1)
        host=np.sum(parents*frozen['eh'][None,None,:,None],axis=(1,2,3))
        cap=np.ones(nr); over=aa+bb>host
        cap[over]=2*host[over]/(bb[over]+np.sqrt(bb[over]**2+4*aa[over]*host[over]))
        actual=aa*cap**2+bb*cap
        if np.any(actual-host>1e-18):
            raise AssertionError('Collective incoming recoil budget exceeded')
        y+=delta*cap[:,None,None]/np.sqrt(ecanon)
        if t>=BURN:
            total_parent=parents.sum(axis=(1,2,3))
            counts['all']+=total_parent
            for name in names[1:]:
                counts[name]+=np.sum(parents*masks[name][:,:,None,None],axis=(1,2,3))
            masks_saved[j]=np.column_stack([np.sum(parents*masks[name][:,:,None,None],axis=(1,2,3))/total_parent for name in ('A_disturbed','C_any','A_and_C')])
            work+=actual; incoming+=host
            self_energy+=np.sum(realized*frozen['eh'][None,None,:,None,None]*fractions*cap[:,None,None,None,None]**2,axis=(1,2,3,4))
            capped+=over
            # Overlapping classes carry signed pre-field work and half of each
            # new-kick cross term; their contributions must not be added.
            for ci,name in enumerate(('A_disturbed','C_any','A_and_C')):
                dc=np.zeros((nr,nm))
                for species in range(2):
                    dc+=(amp[:,:,species]*masks[name])@m['source'][species].T
                classwork[:,ci]+=cap*np.sum(oldp*dc,axis=1)+.5*cap**2*np.sum(dp*dc,axis=1)
        before=ecanon/2*(x.sum(axis=1)**2+y.sum(axis=1)**2)
        x,y=x*co+y*si,y*co-x*si
        x*=damp; y*=damp
        if t>=BURN:
            damping+=np.sum(before*(1-damp*damp),axis=1)
    final=ecanon/2*np.sum(x.sum(axis=1)**2+y.sum(axis=1)**2,axis=1)
    residual=final-initial-work+damping
    if max(abs(residual)/incoming)>1e-9:
        raise AssertionError('Source/propagation energy ledger failed')
    return dict(x=sx,y=sy,B=bmask,counts=counts,masks=masks_saved,gamma=gamma,
                work=work,incoming=incoming,marked=self_energy,loss=damping,residual=residual,capped=capped,
                draw_sha256=[h.hexdigest() for h in marks.sha],classwork=classwork,
                signed_region=signed_region,signed_z=signed_z)


def peak_width(c,m,coefficient,depth):
    # Actual reconstructed full-vector elastic energy on a common 129x129 map.
    from scipy.ndimage import label
    axis=np.linspace(0,c.DIMS[0],129)
    xx,yy=np.meshgrid(axis,axis,indexing='ij')
    pts=np.column_stack([xx.ravel(),yy.ravel(),np.full(xx.size,depth)])
    values=np.empty(len(pts)); original=c.basis
    try:
        c.basis=lambda p:original(p,degree=m['degree'])
        for start in range(0,len(pts),128):
            p=pts[start:start+128]; _,b=c.basis(p)
            e=b@coefficient
            cc=c.tensors(p,c.DIMS[0]/8,'grains')
            values[start:start+len(p)]=.5*np.einsum('qi,qij,qj->q',e,cc,e)
    finally:
        c.basis=original
    z=values.reshape(xx.shape); peak=np.unravel_index(z.argmax(),z.shape)
    labels,_=label(z>=z[peak]/2)
    component=labels==labels[peak]
    ij=np.argwhere(component); step=axis[1]-axis[0]
    widths=(ij.max(axis=0)-ij.min(axis=0))*step
    censored=bool(np.any(ij==0)|np.any(ij==128))
    return dict(map_peak_J_m3=float(z[peak]),width_x_m=float(widths[0]),width_y_m=float(widths[1]),
                map_spacing_m=float(step),foil_censored=censored,
                width_unresolved=bool(np.any(widths<2*step))),z


def compare():
    start=time.monotonic(); c,l,r,frozen=accepted(); portable=r.Portable()
    records=[]; tails=[]; peaks=[]; signed=[]; hist={}; ledger=[]
    for rung in (0,1):
        m=consumer_model(rung,c,l,frozen)
        for spatial,timing in ((0,0),(2,0),(2,2)):
            for beta in (0.,.35):
                for loss in ('retained','Cu_proxy'):
                    arm=f'r{rung}_s{spatial}_t{timing}_b{beta}_{loss}'
                    result=simulate(m,frozen,portable,r,l,spatial,timing,beta,loss)
                    ecanon=m['ref'][1]
                    x=result['x'].sum(axis=2); y=result['y'].sum(axis=2)
                    strain=np.einsum('pim,trm->trpi',m['strain'],x)*np.sqrt(ecanon)
                    u=.5*np.einsum('trpi,pij,trpj->trp',strain,m['cs'],strain)
                    region_u=.5*ecanon*np.einsum('trm,pmn,trn->trp',x,m['regional'][:,0],x,optimize=True)
                    # Positive nonlocal quadratic loss-density proxy. Its full-
                    # domain integral is the exact modal energy-decay rate.
                    root=np.sqrt(result['gamma'])
                    el=np.einsum('pim,trm->trpi',m['strain'],x*root)*np.sqrt(ecanon)
                    vl=np.einsum('pim,trm->trpi',m['velocity'],y*root)*np.sqrt(ecanon)
                    production=.5*np.einsum('trpi,pij,trpj->trp',el,m['cs'],el)+.5*c.RHO*np.sum(vl*vl,axis=-1)
                    region_production=.5*ecanon*(np.einsum('trm,pmn,trn->trp',x*root,m['regional'][:,0],x*root,optimize=True)
                        +np.einsum('trm,pmn,trn->trp',y*root,m['regional'][:,1],y*root,optimize=True))
                    hist[arm+'_energy_J_m3']=u
                    hist[arm+'_loss_proxy_W_m3']=production
                    hist[arm+'_fixed_region_energy_J_m3']=region_u
                    hist[arm+'_fixed_region_loss_proxy_W_m3']=region_production
                    profiles,_=c.profiles_and_stopping()
                    profile=profiles['PdD',2.5]
                    for rep,record in enumerate(RECORDS):
                        row=dict(arm=arm,rung=rung,record=record,spatial=spatial,timing=timing,beta=beta,loss=loss,exposure_s=SAMPLE*DT)
                        for name,n in result['counts'].items():
                            row[name+'_parents']=int(n[rep])
                            if name!='all': row[name+'_fraction']=float(n[rep]/result['counts']['all'][rep])
                        b=l.occurrence(result['B'][:,rep:rep+1],DT)
                        row.update({('B_'+key):b[key] for key in ('time_fraction','complete_events','boundary_censored_events','initiation_rate_s')})
                        row['draw_sha256']=result['draw_sha256'][rep]
                        records.append(row)
                        ledger.append(dict(arm=arm,record=record,incoming_host_J=float(result['incoming'][rep]),source_work_J=float(result['work'][rep]),
                            marked_self_J=float(result['marked'][rep]),marked_physical_fraction=float(result['marked'][rep]/result['incoming'][rep]),
                            damping_J=float(result['loss'][rep]),balance_residual_J=float(result['residual'][rep]),capped_bins=int(result['capped'][rep])))
                        for ri,(side,z0,z1) in enumerate(m['region_bounds']):
                            for oi,observable in enumerate(('epsilon_xx','u_z')):
                                values=result['signed_region'][:,rep,2*ri+oi]
                                zz=result['signed_z'][:,rep,2*ri+oi]
                                for sign,label in ((1,'positive_peak'),(-1,'negative_hole')):
                                    occ=l.occurrence((sign*zz[:,None]>=3),DT)
                                    signed.append(dict(arm=arm,record=record,side_m=side,z0_m=z0,z1_m=z1,observable=observable,sign=label,
                                        threshold_rms=3,physical_min=float(values.min()),physical_max=float(values.max()),
                                        **{key:occ[key] for key in ('time_fraction','complete_events','boundary_censored_events','initiation_rate_s','max_duration_us')}))
                        for pi,point in enumerate(m['points']):
                            energy=u[:,rep,pi]; mean=float(energy.mean()); peak=int(energy.argmax())
                            depth=float(point[2]); nd=c.layer_value(profile['layers'],depth)
                            tail=dict(arm=arm,record=record,point=pi,depth_m=depth,n_D_m3=nd,mean_J_m3=mean,max_over_mean=float(energy[peak]/mean),
                                max_elastic_eV_per_D=float(energy[peak]/(l.QE*nd)),
                                conditional_converted_eV_per_D=float(production[peak,rep,pi]*4.1e-12/(l.QE*nd)))
                            for q in (.5,.9,.95,.99,.999):
                                tail['quantile_'+str(q)+'_over_mean']=float(np.quantile(energy,q)/mean)
                            for h in (3.,4.,5.):
                                occurrence=l.occurrence((energy[:,None]>=h*mean),DT)
                                tail['time_fraction_ge_'+str(h)]=occurrence['time_fraction']
                                tail['complete_events_ge_'+str(h)]=occurrence['complete_events']
                            tails.append(tail)
                        for ri,(z0,z1) in enumerate(((0.,50e-6),(25e-6,75e-6),(175e-6,225e-6))):
                            energy=region_u[:,rep,ri]; peak=int(energy.argmax()); mean=float(energy.mean())
                            nd=c.column_density(profile['layers'],z0,z1)
                            tail=dict(arm=arm,record=record,point=27+ri,depth_m=(30e-6,50e-6,200e-6)[ri],n_D_m3=nd,
                                mean_J_m3=mean,max_over_mean=float(energy[peak]/mean),max_elastic_eV_per_D=float(energy[peak]/(l.QE*nd)),
                                conditional_converted_eV_per_D=float(region_production[peak,rep,ri]*4.1e-12/(l.QE*nd)))
                            for q in (.5,.9,.95,.99,.999):
                                tail['quantile_'+str(q)+'_over_mean']=float(np.quantile(energy,q)/mean)
                            for h in (3.,4.,5.):
                                occurrence=l.occurrence((energy[:,None]>=h*mean),DT)
                                tail['time_fraction_ge_'+str(h)]=occurrence['time_fraction']
                                tail['complete_events_ge_'+str(h)]=occurrence['complete_events']
                            tails.append(tail)
                        # Peak over the same three central 2.5-mm-square,
                        # 50-um-depth regions; never an unseen spacetime max.
                        t,ri=np.unravel_index(region_u[:,rep].argmax(),region_u[:,rep].shape)
                        pi=9*ri+4
                        point=m['points'][pi]
                        coefficient=m['v']@(np.sqrt(ecanon)*x[t,rep]/m['omega'])
                        width,z=peak_width(c,m,coefficient,float(point[2]))
                        hist[arm+f'_record{record}_peak_map_J_m3']=z
                        energy=region_u[:,rep,ri]; threshold=energy[t]/2
                        lo=hi=t
                        while lo>0 and energy[lo-1]>=threshold: lo-=1
                        while hi<SAMPLE-1 and energy[hi+1]>=threshold: hi+=1
                        separate=.5*ecanon*np.einsum('gm,mn,gn->',result['x'][t,rep],m['regional'][ri,0],result['x'][t,rep])
                        peaks.append(dict(arm=arm,record=record,fixed_region=int(ri),fixed_region_peak_J_m3=float(energy[t]),
                            positive_denominator_J_m3=float(energy.mean()),fixed_region_peak_over_mean=float(energy[t]/energy.mean()),
                            duration_half_height_s=float((hi-lo+1)*DT),duration_censored=bool(lo==0 or hi==SAMPLE-1),
                            cross_sector_elastic_J_m3=float(energy[t]-separate),
                            A_signed_source_work_J=float(result['classwork'][rep,0]),C_signed_source_work_J=float(result['classwork'][rep,1]),
                            A_parent_fraction_in_peak_bin=float(result['masks'][t,rep,0]),C_parent_fraction_in_peak_bin=float(result['masks'][t,rep,1]),
                            B_remote_window=bool(result['B'][t,rep]),attribution='mixed evidence; unique parent/echo cause not resolved',**width))
                    print(json.dumps(dict(arm=arm,modes=len(m['omega']),elapsed_s=time.monotonic()-start)),flush=True)
    table('records.csv',records); table('tails.csv',tails); table('peaks.csv',peaks); table('signed.csv',signed); table('energy_ledger.csv',ledger)
    np.savez_compressed(HERE/'histories.npz',**hist)
    dump('comparison.json',dict(independent_records=RECORDS,exposure_s=SAMPLE*DT,burn_s=BURN*DT,dt_s=DT,arms=24,
        quiet_old_band_physical_fraction=.01,new_disjoint_frequency_band_fraction=0,
        max_relative_energy_residual=max(abs(p['balance_residual_J'])/p['incoming_host_J'] for p in ledger),elapsed_s=time.monotonic()-start,
        limits=['fixed central physical-region tail, not whole-domain spacetime extremum','two records do not resolve population rates or uncertainty',
                'Cu proxy deliberately includes no old intrinsic or mount loss','same source quadrature; no impact mechanics resolution']))


def pulse_and_loss():
    # Explicit illustrative bulk Gaussian impulse, never a measured Lamb launch law.
    rows = []
    for a in (1e-8, 1e-7, 1e-6):
        for tau in (1e-13, 5e-13, 1e-11):
            t = np.hypot(tau, a/4460.)
            def frac(lo, hi):
                return float(gammainc(1.5,(2*np.pi*hi*t)**2)-gammainc(1.5,(2*np.pi*lo*t)**2))
            f0 = frac(5000,250000)
            high_normalizer=1-float(gammainc(1.5,(2*np.pi*250000*t)**2))
            rows.append(dict(lateral_gaussian_sigma_m=a, pulse_sigma_s=tau, proxy_bulk_speed_m_s=4460.,
                low_band_fraction_of_whole_gaussian=f0, total_recoil_fraction_needed_for_one_percent=.01/f0,
                disjoint_250kHz_to_1MHz_fraction=frac(250000,1e6), disjoint_1_to_10MHz_fraction=frac(1e6,1e7),
                quiet_250kHz_to_1MHz_share_if_high_reservoir_10percent=.99*.1*frac(250000,1e6)/high_normalizer,
                quiet_250kHz_to_1MHz_share_if_high_reservoir_maximal=.99*frac(250000,1e6)/high_normalizer,
                whole_source_quiet_share_if_high_reservoir_10percent=.01+.99*.1,
                status='assumed bulk illustrative spectrum; separate calibrated low component required'))
    table('pulse_sensitivities.csv', rows)
    rows = []
    for f in (5000,10000,50000,100000,250000,1e6,5e6,10e6):
        db = 32*f/1e6
        aa = db*np.log(10)/20
        rows.append(dict(frequency_Hz=f, alpha_dB_per_m=db, alpha_amplitude_Np_per_m=aa, energy_length_m=1/(2*aa),
            bulk_proxy_energy_decay_per_s=2*aa*4460, bulk_proxy_lifetime_s=1/(2*aa*4460),
            domain='Cu longitudinal bulk proxy; not measured Pd/PdD or Lamb; sub-MHz extrapolation'))
    table('loss_law.csv', rows)


def readout():
    from scipy.ndimage import label
    from scipy.stats import t as student
    s=json.loads((HERE/'support.json').read_text())
    pulse_and_loss()
    def rows(path):
        with path.open() as f:
            return list(csv.DictReader(f))
    record=rows(HERE/'records.csv'); tails=rows(HERE/'tails.csv'); peaks=rows(HERE/'peaks.csv')
    history=np.load(HERE/'histories.npz',allow_pickle=False)
    # Consumer-only width correction: use a peak IN the fixed central region,
    # not an unrelated remote maximum elsewhere in the same frame.
    for p in peaks:
        z=history[p['arm']+f"_record{p['record']}_peak_map_J_m3"]
        axis=np.linspace(0,np.sqrt(np.pi)*.012,129)
        inside=abs(axis-axis[-1]/2)<=.00125
        valid=inside[:,None]&inside[None,:]
        at=np.unravel_index(np.where(valid,z,-np.inf).argmax(),z.shape)
        labels,_=label(z>=z[at]/2)
        ij=np.argwhere(labels==labels[at])
        widths=(ij.max(axis=0)-ij.min(axis=0))*(axis[1]-axis[0])
        p['map_peak_J_m3']=float(z[at]); p['width_x_m']=float(widths[0]); p['width_y_m']=float(widths[1])
        p['foil_censored']=bool(np.any(ij==0)|np.any(ij==128))
        p['width_unresolved']=bool(np.any(widths<2*(axis[1]-axis[0])))
        lam=s['lambda0_screen_m']/(3 if p['arm'].startswith('r1_') else 1)
        p['band_width_censored']=bool(np.min(widths)<lam/2)
        p['attribution']='unattributed; source-class and signed cross-sector evidence are nonexclusive'
    c,l,_,frozen=accepted()
    profiles,_=c.profiles_and_stopping()
    profile=profiles['PdD',2.5]
    event_energy=[]
    for p in peaks:
        record_id=int(p['record']); rep=RECORDS.index(record_id)
        uu=history[p['arm']+'_fixed_region_energy_J_m3']
        at,ri=np.unravel_index(uu[:,rep].argmax(),uu[:,rep].shape)
        p['peak_sample_index']=int(at)
        p['peak_time_after_burn_s']=float(at*DT)
        for depth_index,(depth,z0,z1) in enumerate(((30e-6,0.,50e-6),(50e-6,25e-6,75e-6),(200e-6,175e-6,225e-6))):
            nd=c.layer_value(profile['layers'],depth)
            nd_region=c.column_density(profile['layers'],z0,z1)
            local=history[p['arm']+'_energy_J_m3'][at,rep,9*depth_index+4]
            rate=history[p['arm']+'_loss_proxy_W_m3'][at,rep,9*depth_index+4]
            region_rate=history[p['arm']+'_fixed_region_loss_proxy_W_m3'][at,rep,depth_index]
            event_energy.append(dict(arm=p['arm'],record=record_id,peak_sample_index=int(at),depth_m=depth,n_D_local_m3=nd,
                local_elastic_J_m3=float(local),local_elastic_eV_per_D=float(local/(l.QE*nd)),
                local_conditional_converted_eV_per_D=float(rate*4.1e-12/(l.QE*nd)),region_n_D_m3=nd_region,
                region_elastic_eV_per_D=float(uu[at,rep,depth_index]/(l.QE*nd_region)),
                region_conditional_converted_eV_per_D=float(region_rate*4.1e-12/(l.QE*nd_region))))
    table('leading_event_depth_energy.csv',event_energy)
    table('peaks.csv',peaks)
    def stats(values):
        x=np.asarray(values,dtype=float); n=len(x); mean=float(x.mean())
        sd=float(x.std(ddof=1)) if n>1 else None
        half=float(student.ppf(.975,n-1)*sd/np.sqrt(n)) if n>1 and sd>0 else None
        return dict(n_independent_records=n,mean=mean,sample_sd=sd,
                    ci95_lo=mean-half if half is not None else '',ci95_hi=mean+half if half is not None else '')
    summary=[]
    for arm in sorted(set(p['arm'] for p in record)):
        rr=[p for p in record if p['arm']==arm]
        for name in ('A_disturbed','C_any','C_strain_positive','C_strain_negative','C_displacement_positive','C_displacement_negative','A_and_C','A_or_C'):
            summary.append(dict(arm=arm,observable=name,denominator='physical charge-parent impacts',quantity='fraction',**stats([p[name+'_fraction'] for p in rr])))
            summary.append(dict(arm=arm,observable=name,denominator='physical record seconds',quantity='impact_rate_per_s',
                                **stats([float(p[name+'_parents'])/(SAMPLE*DT) for p in rr])))
        for name in ('B_time_fraction','B_initiation_rate_s'):
            summary.append(dict(arm=arm,observable=name,denominator='off-impact receiver seconds',quantity=name,**stats([p[name] for p in rr])))
    table('class_summary.csv',summary)
    trend=[]
    for arm in sorted(set(p['arm'] for p in tails if p['arm'].startswith('r1_'))):
        old='r0_'+arm[3:]
        for point in (27,28,29):
            for metric in ('mean_J_m3','max_over_mean','max_elastic_eV_per_D','conditional_converted_eV_per_D'):
                aa={int(p['record']):float(p[metric]) for p in tails if p['arm']==arm and int(p['point'])==point}
                bb={int(p['record']):float(p[metric]) for p in tails if p['arm']==old and int(p['point'])==point}
                ratio=[aa[r]/bb[r] for r in RECORDS]
                trend.append(dict(new_arm=arm,old_arm=old,point=point,quantity=metric,comparison='paired rung1/rung0 ratio',
                                  **stats(ratio)))
        for metric in ('fixed_region_peak_over_mean','width_x_m','width_y_m','duration_half_height_s'):
            aa={int(p['record']):float(p[metric]) for p in peaks if p['arm']==arm}
            bb={int(p['record']):float(p[metric]) for p in peaks if p['arm']==old}
            trend.append(dict(new_arm=arm,old_arm=old,point='region maximum',quantity=metric,comparison='paired rung1/rung0 ratio',
                              **stats([aa[r]/bb[r] if bb[r] else np.nan for r in RECORDS])))
    table('trend.csv',trend)
    # Omit physically unsupported whole-foil average translation, rather than
    # divide roundoff by roundoff and count artificial standardized holes.
    c,l,r,frozen=accepted()
    signed=rows(HERE/'signed.csv'); omitted=[]; support=[]
    for rung in (0,1):
        if rung==0:
            degree=(5,5,2); v=frozen['eigenvectors']; omega=frozen['omega']
        else:
            d=np.load(HERE/'model.npz',allow_pickle=False); degree=tuple(d['degree']); v=d['eigenvectors']; omega=d['omega']
        bounds=[(side,z0,z1) for side in (.0025,.01,float(c.DIMS[0])) for z0,z1 in ((0.,50e-6),(25e-6,75e-6),(175e-6,225e-6))]
        proj=np.array([linear_box(c,v,omega,degree,*b,c.DIMS[:2]/2) for b in bounds])
        norms=np.linalg.norm(proj,axis=2)
        for i,b in enumerate(bounds):
            for oi,name in enumerate(('epsilon_xx','u_z')):
                relative=norms[i,oi]/norms[:,oi].max()
                support.append(dict(rung=rung,side_m=b[0],z0_m=b[1],z1_m=b[2],observable=name,relative_projector_norm=relative,
                                    supported=bool(relative>1e-10)))
                if relative<=1e-10:
                    omitted.append((rung,*b,name))
    kept=[]
    for p in signed:
        rung=int(p['arm'][1])
        key=(rung,float(p['side_m']),float(p['z0_m']),float(p['z1_m']),p['observable'])
        if key not in omitted:
            kept.append(p)
    table('signed.csv',kept); table('observable_support.csv',support)
    # Compare changed consumer against already portable accepted observations;
    # do not replay the historical producer for an acceptance ritual.
    old_records=rows(ROOT/'g13_rogue/pass2/replication/records.csv')
    control=[]
    for p in record:
        if not p['arm'].startswith('r0_') or p['loss']!='retained':
            continue
        q=next(q for q in old_records if q['record']==p['record'] and q['spatial']==p['spatial']
               and q['timing']==p['timing'] and float(q['beta'])==float(p['beta']))
        control.append(dict(arm=p['arm'],record=p['record'],parent_difference=int(p['all_parents'])-int(q['parents']),
            A_parent_difference=int(p['A_disturbed_parents'])-int(q['A_disturbed_parents']),
            C_parent_difference=int(p['C_any_parents'])-int(q['C_any_parents']),
            B_fraction_difference=float(p['B_time_fraction'])-float(q['B_time_fraction']),
            addressed_mark_digest_equal=p['draw_sha256']==q['draw_sha256']))
    table('accepted_consumer_comparison.csv',control)
    matched={}
    for p in record:
        key=(p['record'],p['spatial'],p['timing'])
        matched.setdefault(key,set()).add(int(p['all_parents']))
    if any(len(v)!=1 for v in matched.values()):
        raise AssertionError('Physical record charge is not matched across rungs/loss/feedback')
    quiet={}
    for p in record:
        if float(p['beta'])==0:
            quiet.setdefault((p['record'],p['spatial'],p['timing']),set()).add(p['draw_sha256'])
    if any(len(v)!=1 for v in quiet.values()):
        raise AssertionError('Quiet fixed-state marks diverged across numerical spaces')
    result=dict(physical_parent_totals_match=True,quiet_addressed_mark_digests_match=True,
        accepted_control=control,signed_unsupported_rows_omitted=len(signed)-len(kept),
        uncertainty='two independent records; pointwise t intervals are descriptive and not model uncertainty',
        trend_rows=len(trend),class_summary_rows=len(summary))
    dump('consumer_verification.json',result)
    print(json.dumps(result),flush=True)


def package():
    names = sorted(p.name for p in HERE.iterdir() if p.is_file() and p.name not in ('manifest.json','capsule.tar.gz','delivery.json'))
    files = [dict(path=PREFIX+n, bytes=(HERE/n).stat().st_size, sha256=sha(HERE/n)) for n in names]
    dump('manifest.json', dict(base=BASE, unpublished=True, files=files))
    out = io.BytesIO()
    with gzip.GzipFile(fileobj=out, mode='wb', mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode='w') as tar:
            for n in names+['manifest.json']:
                data=(HERE/n).read_bytes(); item=tarfile.TarInfo(PREFIX+n)
                item.size=len(data); item.mode=0o644; item.mtime=0
                tar.addfile(item, io.BytesIO(data))
    data=out.getvalue(); (HERE/'capsule.tar.gz').write_bytes(data)
    receipt=dict(path=PREFIX+'capsule.tar.gz', sha256=hashlib.sha256(data).hexdigest(), bytes=len(data), members=len(names)+1)
    dump('delivery.json',receipt)
    print(json.dumps(receipt), flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare','freeze','compare','readout','package'))
    a=p.parse_args()
    globals()[a.action]()


if __name__ == '__main__':
    main()
