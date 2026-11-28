"""Frozen source adapters. No predecessor producer is invoked or written."""
from pathlib import Path
import hashlib, importlib.util, itertools, json, math, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import erf
from functools import lru_cache
HERE=Path(__file__).resolve().parent
QE=1.602176634e-19; MD=2.0141*1.66053906660e-27; HB=6.582119569e-16; KB=8.617333262145e-5
L=np.sqrt(np.pi)*.012; H=.00025; RHO=12007.; CV=RHO*25.89/.1064; KTH=71.8; ALPHA=11.77e-6
BULK=(227.1e9+2*176.04e9)/3
ENERGIES=(.25,.5,1.,1.5,2.,2.5); BEAMS=('even_steady','hot_steady','hot_bursts')
ARMS=[dict(name='fast_Cu',share=.25,ne=15e-15,opt=.0802e-12,pd=.15e-12,te=.3e-12,recoil=.05e-12,end=.1e-12,struct=1.3e-9,radius=8e-9,lambda_e=10e-9,lambda_pd=2e-9,dephase=.08e-12,loss='Cu'),
 dict(name='nominal_retained',share=.5,ne=.1e-12,opt=.5e-12,pd=1e-12,te=.6e-12,recoil=.2e-12,end=.5e-12,struct=100e-9,radius=20e-9,lambda_e=20e-9,lambda_pd=10e-9,dephase=.5e-12,loss='retained'),
 dict(name='long_retained',share=.92,ne=.5e-12,opt=2.068018471814257e-12,pd=5e-12,te=.3e-12,recoil=1e-12,end=2e-12,struct=1e-6,radius=60e-9,lambda_e=30e-9,lambda_pd=30e-9,dephase=2e-12,loss='retained')]
PHYSICAL_ARMS=ARMS
ARMS=[dict(a,name=a['name']+'_'+allocation,closure=a['name'],allocation=allocation,mesh=16)
      for a in PHYSICAL_ARMS for allocation in ('onsite','coupling')]
NOMINAL=next(a for a in ARMS if a['closure']=='nominal_retained' and a['allocation']=='onsite')

def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def load_file(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def inputs():return json.loads((HERE/'frozen/pass1_inputs.json').read_text())

def beam(x,y,name):
    mask=(abs(x)<=.005)&(abs(y)<=.005)
    if name=='even_steady':return mask.astype(float)
    w=np.exp(-(x*x+y*y)/(2*.003**2))
    for a,b,g in ((.0025,-.0025,5.),(-.0025,.0025,3.)):
        w+=g*np.exp(-((x-a)**2+(y-b)**2)/(2*.0012**2))
    return w*mask

def beam_cells(edges,name):
    """Exact Gaussian/cropped-square integral, including partial boundary cells."""
    lo=np.maximum(np.asarray(edges[:-1]),-.005);hi=np.minimum(np.asarray(edges[1:]),.005)
    width=np.maximum(hi-lo,0.)
    if name=='even_steady':return np.outer(width,width)/1e-4
    result=np.zeros((len(lo),len(lo)))
    for a,b,g,sigma in ((0.,0.,1.,.003),(.0025,-.0025,5.,.0012),(-.0025,.0025,3.,.0012)):
        ix=np.where(hi>lo,sigma*np.sqrt(np.pi/2)*(erf((hi-a)/(np.sqrt(2)*sigma))-erf((lo-a)/(np.sqrt(2)*sigma))),0.)
        iy=np.where(hi>lo,sigma*np.sqrt(np.pi/2)*(erf((hi-b)/(np.sqrt(2)*sigma))-erf((lo-b)/(np.sqrt(2)*sigma))),0.)
        result+=g*np.outer(ix,iy)
    return result/result.sum()

@lru_cache(maxsize=3)
def beam_moments(name):
    nodes,weights=np.polynomial.legendre.leggauss(64);nodes*=.005;weights*=.005
    x,y=np.meshgrid(nodes,nodes,indexing='ij');w=beam(x,y,name);area=np.outer(weights,weights)
    norm=float(np.sum(w*area));second=float(np.sum(w*w*area))
    return norm,1e-4*second/(norm*norm)

def current(ec):
    return np.interp(2*ec,np.arange(.5,15.01,.5),[38,68,94,110,118,123,130,134,136,138,140,145,152,156,160,165,171,176,180,182,192,207,217,222,228,238,251,263,275,283])*1e-6

def tracks(d):
    bca=load_file(HERE/'frozen/bca.py','frozen_bca');out={}
    for ec in ENERGIES:
        profile=next(r for r in d['profiles'] if r['case']['E_cm_keV']==ec)
        for sp in (0,1):
            saved=next((r for r in d['tracks'] if r['E_cm_keV']==ec),None)
            if saved:tr=saved['tracks'][sp]
            else:
                energy=2000*ec/(sp+1);dep=np.zeros(3);t=0.;rows=[]
                for zt,zb,nh,nd,z,m,_ in profile['layers']:
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
                tr=dict(initial_eV=2000*ec/(sp+1),terminal_unresolved_eV=energy,terminal_time_fs=t,terminal_depth_A=depth,timeline=rows,components_eV=dep.tolist())
            out[f'{ec:g}_{sp}']=tr
    return out

def alias_bonds(p):
    """Exactly lattice.py:shell_kernel sign-image weights, in Cartesian metres."""
    rr=[];jj=[]
    for shell in p['shells']:
        aliases=shell['periodic_alias_offsets_A'];base=shell['J_sum_eV_A2']/len(aliases)
        for raw in aliases:
            r=np.rint(2*np.asarray(raw)/p['a_A']).astype(int)
            images=list(itertools.product(*[(-2,2) if abs(v)==2 else (int(v),) for v in r]))
            for image in images:
                rr.append(np.asarray(image)*p['a_A']/2*1e-10);jj.append(base/len(images))
    return np.array(rr),np.array(jj)

def optical(p,k):
    r,j=alias_bonds(p);A=p['K4_100_eV_A4'];B=(3*p['K4_111_eV_A4']-A)/2
    K=p['Kii_eV_A2']+(3*A+2*B)*p['zero_point_axis_variance_A2']
    phase=np.asarray(k)@r.T;D=K-np.cos(phase)@j
    if np.any(D<=0):raise ValueError('Unsupported unstable optical projection')
    mass=p['mass_eVc2']*QE/(299792458.**2)
    w=np.sqrt(D*QE/1e-20/mass)
    grad=np.einsum('...b,b,bi->...i',np.sin(phase),j,r)*QE/1e-20/(2*mass*w[...,None])
    return w,grad,D

def refresh_optical():
    """Rebuild only the small exact G3 projection; never repeat the elastic solve."""
    p=inputs()['lattice']
    import native
    k=native.primitive_basis(p);w,vg,D=optical(p,k)
    native.soft_basis()
    np.savez_compressed(HERE/'frozen/optical_branches.npz',k=k,omega=w,vg=vg,D=D,mass_kg=p['mass_eVc2']*QE/(299792458.**2))
    meta=json.loads((HERE/'source_index.json').read_text());meta['optical_mass_owner']='G3 mass_eVc2, not upstream atomic stopping mass'
    meta['optical_mass_kg']=p['mass_eVc2']*QE/(299792458.**2);meta['source_files']=[]
    for path in sorted((HERE/'frozen').iterdir()):
        if path.is_file():meta['source_files'].append(dict(path='frozen/'+path.name,bytes=path.stat().st_size,sha256=digest(path)))
    (HERE/'source_index.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(dict(optical_projection_refreshed=True,mass_kg=meta['optical_mass_kg'],maximum_group_m_s=float(np.linalg.norm(vg,axis=1).max()))),flush=True)

def prepare(reference):
    root=HERE.parents[1];reference=Path(reference)
    expected={'model.npz':'59697212d397a6233f7aee9a4d02bd07ac773c8d3d380d50eaba116dffbcef8d','signed_r2.npz':'55d3d2b85f66d12d77bee585e1dcb385e7e4f307fb2372646f85426ec7849548'}
    for name,sha in expected.items():
        if digest(reference/name)!=sha:raise ValueError('Wrong accepted modal epoch: '+name)
    d=inputs();tr=tracks(d)
    (HERE/'frozen/tracks.json').write_text(json.dumps(tr,indent=2)+'\n')
    m=np.load(reference/'model.npz');s=np.load(reference/'signed_r2.npz')
    if not np.array_equal(m['omega'],s['omega']):raise ValueError('Modal epochs differ')
    spatial=load_file(root/'g13_rogue/pass3/rung2/spatial.py','frozen_modal_projection')
    edges=np.linspace(0,L,18);x=(edges[:-1]+edges[1:])/2;z=np.array([.125e-6,62.625e-6,187.5e-6])
    points=np.array(list(itertools.product(x,x,z)))
    u,e=spatial.maps(points,m['eigenvectors'],np.array([L,L,H]))
    # Keep the exact 303-dimensional modal metric; project observable rows only.
    np.savez_compressed(HERE/'frozen/modes.npz',omega=m['omega'],displacement=u,strain=e,points=points,
        source_xy=s['xy'],pressure=s['pressure'],heatmap=s['heatmap'],jcos=s['jcos'],heatcos=s['heatcos'],heatenergy=s['heatenergy'],speed=s['speed'])
    rng=np.random.default_rng(4401);k=rng.uniform(-np.pi,np.pi,(64,3))/(d['lattice']['a_A']*1e-10)
    w,vg,D=optical(d['lattice'],k)
    np.savez_compressed(HERE/'frozen/optical_branches.npz',k=k,omega=w,vg=vg,D=D)
    meta=dict(base='fd76b0cc3a26ac5e77f9da0577d46c9c583b7eb1',original_modes=expected,
        projection='exact Q4 x Q4 x P4 interpolation and physical derivatives at 17 x 17 x 3 FV centers; all 303 native mass-normalized coefficients retained',
        mode_count=len(m['omega']),point_count=len(points),original_eigenvector_shape=list(m['eigenvectors'].shape),
        optical_alias_bonds=len(alias_bonds(d['lattice'])[1]),source_files=[])
    for path in sorted((HERE/'frozen').iterdir()):
        if path.is_file():meta['source_files'].append(dict(path='frozen/'+path.name,bytes=path.stat().st_size,sha256=digest(path)))
    (HERE/'source_index.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps(dict(prepared=True,modes=len(m['omega']),points=len(points),max_sampled_optical_group_m_s=float(np.linalg.norm(vg,axis=1).max()))),flush=True)


def refresh_source_points(reference):
    """Additional exact observable rows; no repeated modal solve/projection."""
    root=HERE.parents[1];reference=Path(reference)
    data=np.load(reference/'model.npz');signed=np.load(reference/'signed_r2.npz')
    spatial=load_file(root/'g13_rogue/pass3/rung2/spatial.py','frozen_source_projection')
    points=np.column_stack((signed['xy'],np.zeros(len(signed['xy']))))
    displacement,strain=spatial.maps(points,data['eigenvectors'],np.array([L,L,H]))
    np.savez_compressed(HERE/'frozen/source_modes.npz',points=points,displacement=displacement,strain=strain,
        original_sha256=np.array('59697212d397a6233f7aee9a4d02bd07ac773c8d3d380d50eaba116dffbcef8d'))
    refresh_optical()
