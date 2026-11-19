"""Resolved axial-continuum / stochastic-ZBL pilot. Units: eV, nm, amu.
Run serially; see README for the finite-model interpretation and limitations.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import resource
import time

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import brentq
from scipy.special import kv

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent

def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

bca = load('channel_bca', 'g4_beam/transport/bca.py')
fwd = load('channel_forward', 'g4_beam/forward/model.py')
MD, MT, E2, A0 = 2.0141, 3.01605, 1.439964, .0529177
COEF = np.array([.18175, .50986, .28022, .028171])
BETA = np.array([3.19980, .94229, .40290, .20162])
# Atomic exponential density rho=exp(-r/b)/(8*pi*b^3); b is a model input.
# Slater effective-valence radial-scale estimates, NOT a DFT charge density.
MATS = {
    'PdD': dict(a=.4084, Z=46, M=106.42, x=1, b=.045, u=.008, phonon=.037),
    'TiD2': dict(a=.444, Z=22, M=47.867, x=2, b=.065, u=.008, phonon=.050),
    'LiD': dict(a=.4072, Z=3, M=6.94, x=1, b=.080, u=.012, phonon=.070),
}
FCC = np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
OCT = (FCC + [.5,0,0]) % 1
TET = np.array(np.meshgrid(*[[.25,.75]]*3)).reshape(3,-1).T
AXES = {'100': [[0,1,0],[0,0,1],[1,0,0]],
        '110': [[1,-1,0],[0,0,1],[1,1,0]],
        '111': [[1,-1,0],[1,1,-2],[1,1,1]]}
ENERGIES = [2000., 13000., 50000., 200000.]
DIVERGENCES = [0., .5, 2.]


def thermal_u(mass, phonon):
    return np.sqrt(4.18016e-5/(2*mass*phonon)/np.tanh(phonon/(2*8.617333262e-5*300)))


def geometry(mat, axis):
    a = mat['a']
    vectors = np.array(AXES[axis], float)
    lengths = np.linalg.norm(vectors, axis=1)*a
    unit = vectors/np.linalg.norm(vectors, axis=1)[:,None]
    translations = np.array(np.meshgrid(*[np.arange(-3,4)]*3)).reshape(3,-1).T
    result = []
    for basis in [FCC, TET if mat['x']==2 else OCT]:
        points = (translations[:,None,:]+basis).reshape(-1,3)*a
        coords = points @ unit.T
        ok = np.all((coords >= -1e-9) & (coords < lengths-1e-9), axis=1)
        coords = coords[ok]
        xy, counts = np.unique(np.round(coords[:,:2],10), axis=0, return_counts=True)
        result.append((xy, counts/lengths[2]))
    return lengths, result


def fields(mat, axis, mass, grid, density_scale=1., vibration_scale=1.):
    lengths, rows = geometry(mat, axis)
    kx = 2*np.pi*np.fft.fftfreq(grid, d=lengths[0]/grid)
    ky = 2*np.pi*np.fft.fftfreq(grid, d=lengths[1]/grid)
    KX, KY = np.meshgrid(kx,ky,indexing='ij')
    k2 = KX*KX+KY*KY
    potentials, electron, nuclear = [], [], []
    ns = [4/mat['a']**3, mat['x']*4/mat['a']**3]
    fuel_phonon = mat['phonon']*np.sqrt(MD/mass)
    us = [mat['u']*vibration_scale, thermal_u(mass,fuel_phonon)*vibration_scale]
    X,Y = np.meshgrid(np.arange(grid)/grid*lengths[0],np.arange(grid)/grid*lengths[1],indexing='ij')
    for j, ((positions, linear), Z, u, b) in enumerate(zip(
            rows, [mat['Z'],1], us, [mat['b']*density_scale,A0/2*density_scale])):
        nk = np.zeros((grid,grid), complex)
        for r, line in zip(positions,linear):
            nk += line/(lengths[0]*lengths[1])*np.exp(-1j*(KX*r[0]+KY*r[1]))
        dw = np.exp(-.5*u*u*k2)
        az = .8854*A0/(1+Z**.23)
        uk = 4*np.pi*Z*E2*nk*dw*sum(c/(k2+(beta/az)**2) for c,beta in zip(COEF,BETA))
        potentials.append(uk)
        electron.append(np.fft.ifft2(nk/ns[j]*dw/(1+(b*b)*k2)**2).real*grid**2)
        # Positive real-space Gaussian avoids Fourier truncation ringing of
        # narrow nuclear rows. Periodic nearest images suffice here (u << cell).
        nuclei = np.zeros((grid,grid))
        for r,line in zip(positions,linear):
            dx=(X-r[0]+lengths[0]/2)%lengths[0]-lengths[0]/2
            dy=(Y-r[1]+lengths[1]/2)%lengths[1]-lengths[1]/2
            nuclei += line/(2*np.pi*u*u)*np.exp(-(dx*dx+dy*dy)/(2*u*u))/ns[j]
        nuclear.append(nuclei/nuclei.mean())
    uk = sum(potentials)
    U = np.fft.ifft2(uk).real*grid**2
    fx = np.fft.ifft2(-1j*KX*uk).real*grid**2
    fy = np.fft.ifft2(-1j*KY*uk).real*grid**2
    if min(x.min() for x in electron) < 0:
        raise ValueError('Underresolved projected electronic density')
    stacked = np.array([U,fx,fy,*electron,*nuclear])
    # Distance to nearest host row: geometrical continuum validity boundary.
    distance = np.full((grid,grid),np.inf)
    for r in rows[0][0]:
        dx = (X-r[0]+lengths[0]/2)%lengths[0]-lengths[0]/2
        dy = (Y-r[1]+lengths[1]/2)%lengths[1]-lengths[1]/2
        distance = np.minimum(distance,np.hypot(dx,dy))
    d = 1/np.mean(rows[0][1])
    Z = mat['Z']; az = .8854*A0/(1+Z**.23)
    def curvature(r):
        q=BETA/az
        return 2*Z*E2/d*np.sum(COEF*q*q*(kv(0,q*r)+kv(2,q*r))/2)
    egrid=np.geomspace(10.,200000.,256)
    rc=np.array([brentq(lambda r:curvature(r)-8*e/d**2,1e-6,10.) for e in egrid])
    rc=np.sqrt(rc*rc+(2*us[0])**2)
    barrier=[]; open_fraction=[]
    h=2*max(lengths[:2])/grid
    for r in rc:
        valid=distance>r
        shell=np.abs(distance-r)<h
        barrier.append(float(U[shell].min()) if valid.any() and shell.any() else float(U.min()))
        open_fraction.append(float(valid.mean()))
    return dict(array=stacked, lengths=lengths, n=np.array(ns), u=us, d=d,
                egrid=egrid, barrier=np.array(barrier), rc=rc, open_fraction=open_fraction,
                distance=distance, Umin=float(U.min()), Umax=float(U.max()), grid=grid)


def sample(field, xy, indexes=None):
    coords=(xy/field['lengths'][:2]*field['grid']).T
    arrays=field['array'] if indexes is None else field['array'][indexes]
    return np.array([map_coordinates(a,coords,order=1,mode='grid-wrap',prefilter=False) for a in arrays])


def cross_sections(E, mass):
    ecm=E/1000*mass/(MD+mass)
    if mass==MD:
        p=fwd.reaction_sigma(ecm,0.,'p','proton')
        n=fwd.reaction_sigma(ecm,0.,'n','proton')
        return np.stack([p,n],axis=-1)*1e-13
    e=ecm
    S=(6.927e4+e*(7.454e8+e*(2.050e6+e*5.2002e4)))/(1+e*(6.38e1+e*(-.995+e*(6.981e-5+e*1.728e-4))))
    sig=S/e*np.exp(-34.3827/np.sqrt(e))
    return np.stack([sig,np.zeros_like(sig)],axis=-1)*1e-13


class Scatter:
    def __init__(self, mat, mass):
        self.mat=mat; self.mass=mass
        self.n=np.array([4/mat['a']**3,mat['x']*4/mat['a']**3])
        self.N=self.n.sum(); self.L=self.N**(-1/3)
        self.pmax=1/np.sqrt(np.pi*self.N*self.L)
        self.eg=np.geomspace(10.,200000.,192)
        self.bg=(np.arange(192)+.5)/192
        self.tables=[]
        for Z,M in [(mat['Z'],mat['M']),(1,mass)]:
            a=.8854*A0/(1+Z**.23)
            E=np.repeat(self.eg,192)
            eps=a*E*M/(MD+M)/(Z*E2)
            B=np.tile(self.pmax*np.sqrt(self.bg)/a,192)
            self.tables.append(bca.cm_angle(eps,B).reshape(192,192))

    def angle(self, E, j, r):
        coords=np.array([(np.log(E)-np.log(10.))/np.log(200000./10.)*191,r*192-.5])
        return map_coordinates(self.tables[j],coords,order=1,mode='nearest',prefilter=False)

    def electronic(self, E):
        return np.array([bca.se_ls(E,self.mat['Z']),bca.se_ls(E,1)])


def rotate(v, psi, phi):
    ref=np.zeros_like(v); ref[:,2]=1
    polar=np.abs(v[:,2])>.9
    ref[polar]=[1,0,0]
    e1=np.cross(v,ref); e1/=np.linalg.norm(e1,axis=1)[:,None]
    e2=np.cross(v,e1)
    return v*np.cos(psi)[:,None]+(e1*np.cos(phi)[:,None]+e2*np.sin(phi)[:,None])*np.sin(psi)[:,None]


def run_batch(mat, axis, mass, scatter, nions, seed, grid=96, step=.025,
              density_scale=1., vibration_scale=1., switch_angle=.20, only=None):
    rng=np.random.default_rng(seed)
    cases=[(e,d) for e in ENERGIES for d in DIVERGENCES] if only is None else only
    labels=np.repeat(np.arange(len(cases)),nions); n=len(labels)
    E0=np.array([cases[i][0] for i in labels]); E=E0.copy()
    divergence=np.deg2rad([cases[i][1] for i in labels])
    v=np.zeros((n,3)); v[:,:2]=rng.normal(size=(n,2))*divergence[:,None]/np.sqrt(2)
    v[:,2]=1; v/=np.linalg.norm(v,axis=1)[:,None]
    random=axis=='random'
    field=None if random else fields(mat,axis,mass,grid,density_scale,vibration_scale)
    xy=rng.random((n,2))*(1 if random else field['lengths'][:2])
    z=np.zeros(n); path=np.zeros(n); y=np.zeros((n,2)); encounter=np.zeros(n)
    eloss=np.zeros(n); nloss=np.zeros((n,2)); back=np.zeros(n,bool)
    captured=np.zeros(n,bool); exited=np.zeros(n,bool); exit_length=np.full(n,np.nan)
    hybrid=np.full(n,random); active=np.ones(n,bool)
    optical_depth=rng.exponential(size=(n,2))
    initial_ep=np.zeros(n); initial_bar=np.zeros(n)
    if not random:
        vals=sample(field,xy)
        initial_ep=E*(v[:,0]**2+v[:,1]**2)+vals[0]
        initial_bar=np.interp(E,field['egrid'],field['barrier'])
        captured=initial_ep<initial_bar
    iterations=0
    start=time.monotonic()
    while active.any():
        ids=np.flatnonzero(active); ee=E[ids]; vv=v[ids].copy(); xx=xy[ids].copy()
        crystal=~hybrid[ids]
        vals=np.ones((7,len(ids))); vals[:3]=0
        if crystal.any():
            vals[:,crystal]=sample(field,xx[crystal])
        density=vals[3:5]
        nucdensity=vals[5:7]
        se=scatter.electronic(ee)*scatter.n[:,None]*density
        rate=scatter.n[:,None]/scatter.N/scatter.L*nucdensity
        # Resolve deterministic transverse work, projected densities, electronic loss,
        # and nuclear collision hazard. No trajectory is held at a step budget.
        ds=np.minimum(np.where(crystal,step,2*scatter.L),.025*ee/np.maximum(se.sum(0),1e-20))
        # Event clocks give exact homogeneous Poisson free paths. Crystal steps
        # additionally resolve spatially varying optical depth.
        ds=np.minimum(ds,np.min(optical_depth[ids].T/np.maximum(rate,1e-100),axis=0))
        if crystal.any():
            force=np.hypot(vals[1],vals[2])
            ds=np.minimum(ds,np.sqrt(.008*ee/np.maximum(force,1e-20)))
            speedxy=np.linalg.norm(vv[:,:2],axis=1)
            ds[crystal]=np.minimum(ds[crystal],.3*min(field['u'])/np.maximum(speedxy[crystal],1e-5))
        ds=np.minimum(ds,2*(ee-10)/np.maximum(se.sum(0),1e-20))
        # Midpoint energy reaction integral on the SAME position/direction/path.
        dE=se.sum(0)*ds
        emid=np.maximum(ee-dE/2,10.)
        y[ids]+=cross_sections(emid,mass)*(scatter.n[1]*nucdensity[1]*ds)[:,None]
        encounter[ids]+=nucdensity[1]*ds
        eloss[ids]+=dE
        ee-=dE
        # Conservative transverse continuum force: d(v_hat_perp)/ds=F/(2 E).
        # It changes direction, not kinetic energy: total longitudinal work averages out.
        dv=np.zeros_like(vv); dv[:,:2]=vals[1:3].T*ds[:,None]/(2*np.maximum(ee,10)[:,None])
        dv-=vv*np.sum(dv*vv,axis=1)[:,None]
        vv+=dv; vv/=np.linalg.norm(vv,axis=1)[:,None]
        xx+=vv[:,:2]*ds[:,None]
        z[ids]+=vv[:,2]*ds
        path[ids]+=ds
        # Independent local thermal-nucleus / fuel-nucleus kicks and their exact recoil.
        for j,M in enumerate([mat['M'],mass]):
            optical_depth[ids,j]-=rate[j]*ds
            hit=(optical_depth[ids,j]<=1e-10)&(ee>10)
            if hit.any():
                h=np.flatnonzero(hit)
                theta=scatter.angle(ee[h],j,rng.random(len(h)))
                recoil=4*MD*M/(MD+M)**2*ee[h]*np.sin(theta/2)**2
                psi=np.arctan2(np.sin(theta),np.cos(theta)+MD/M)
                vv[h]=rotate(vv[h],psi,2*np.pi*rng.random(len(h)))
                ee[h]-=recoil; nloss[ids[h],j]+=recoil
                optical_depth[ids[h],j]=rng.exponential(size=len(h))
        if crystal.any():
            h=np.flatnonzero(crystal)
            end=sample(field,xx[h],indexes=[0])[0]
            ep=ee[h]*(vv[h,0]**2+vv[h,1]**2)+end
            bar=np.interp(np.maximum(ee[h],10),field['egrid'],field['barrier'])
            lost=(ep>=bar)&captured[ids[h]]&~exited[ids[h]]
            exit_length[ids[h[lost]]]=path[ids[h[lost]]]
            exited[ids[h[lost]]]=True
            # Physical continuum exit: no closed allowed transverse region, or
            # oblique flight outside the paraxial approximation AND above barrier.
            invalid=(bar<=field['Umin']+1e-8)|((np.hypot(vv[h,0],vv[h,1])>switch_angle)&(ep>=bar))
            hybrid[ids[h[invalid]]]=True
        E[ids]=ee; v[ids]=vv; xy[ids]=xx
        reflected=z[ids]<0
        back[ids[reflected]]=True
        active[ids[(ee<=10.000001)|reflected]]=False
        iterations+=1
        if iterations>2000000:
            raise RuntimeError('Unfinished trajectories: no censored output is permitted')
    records=[]
    for i,(e,d) in enumerate(cases):
        sel=labels==i; cap=sel&captured; lost=cap&exited
        Y=y[sel].mean(0); Q=np.array([4.033e6,3.269e6]) if mass==MD else np.array([17.589e6,0.])
        total=y[sel].sum(1); ret=y[sel]@Q/e
        row=dict(material=next(k for k,vv in MATS.items() if vv['Z']==mat['Z']), isotope='D' if mass==MD else 'T',
                 axis=axis,E_lab_keV=e/1000,divergence_rms_deg=d,nions=nions,
                 yield_primary=Y[0],yield_secondary=Y[1],yield_total=Y.sum(),
                 yield_se=float(total.std(ddof=1)/np.sqrt(nions)),energy_return=float(ret.mean()),
                 energy_return_se=float(ret.std(ddof=1)/np.sqrt(nions)),
                 fuel_flux=float(encounter[sel].sum()/path[sel].sum()),
                 electronic_eV_nm=float(eloss[sel].sum()/path[sel].sum()),
                 nuclear_host_eV_nm=float(nloss[sel,0].sum()/path[sel].sum()),
                 nuclear_fuel_eV_nm=float(nloss[sel,1].sum()/path[sel].sum()),
                 stopping_eV_nm=float((eloss[sel]+nloss[sel].sum(1)).sum()/path[sel].sum()),
                 path_nm=float(path[sel].mean()),projected_endpoint_nm=float(z[sel].mean()),
                 implanted_projected_nm=float(z[sel&~back].mean()) if (sel&~back).any() else None,
                 reflected_fraction=float(back[sel].mean()),captured_fraction=float(captured[sel].mean()),
                 dechannel_fraction_of_captured=float(exited[cap].mean()) if cap.any() else None,
                 dechannel_length_observed_mean_nm=float(exit_length[lost].mean()) if lost.any() else None,
                 captured_stopped_before_exit=int((cap&~exited).sum()),
                 survival_restricted_mean_nm=float(np.where(exited[cap],exit_length[cap],path[cap]).mean()) if cap.any() else None,
                 hybrid_fraction=float(hybrid[sel].mean()))
        if not random:
            psi=np.sqrt(2*mat['Z']*E2/(e*field['d']))
            row.update(lindhard_psi_deg=float(np.rad2deg(psi)),
                       potential_angle_deg=float(np.rad2deg(np.sqrt(max(0,np.interp(e,field['egrid'],field['barrier'])-field['Umin'])/e))),
                       initial_open_area_fraction=float(np.interp(e,field['egrid'],field['open_fraction'])),
                       host_row_period_nm=float(field['d']),fuel_u_nm=float(field['u'][1]))
        records.append(row)
    return records,dict(seconds=time.monotonic()-start,iterations=iterations,
                        peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)


def references(nions, seed):
    results=[]
    for name,x in [('PdD',.7),('TiD2',1.7)]:
        original=MATS[name]
        forward_mat='PdD' if name=='PdD' else 'TiD'
        density=fwd.host_density(forward_mat)*1e-21
        mat=dict(original,x=x,a=(4/density)**(1/3))
        # Historical transport fills, explicitly not the production stoichiometry.
        sc=Scatter(mat,MD)
        rows,usage=run_batch(mat,'random',MD,sc,nions,seed,step=.025,only=[(2000.,0.),(13000.,0.)])
        for row in rows:
            e=row['E_lab_keV']*1000
            comps=[(mat['Z'],mat['M'],1/(1+x),25),(1,MD,x/(1+x),0)]
            ref=bca.run(e,comps,sc.N,nions=nions,Ecut=10,zmax=10000,dz=.25,seed=seed+19)
            target=fwd.Target('PdD' if name=='PdD' else 'TiD',E_cut=.01,E_max=200,N=4096,reaction='proton')
            p=target.Y_const(e/1000,0.,'p'); nn=target.Y_const(e/1000,0.,'n')
            results.append(dict(material=name,fill=x,E_lab_keV=e/1000,
                                our_random=row,bca_implanted_projected_nm=ref['Rp'],bca_reflection=ref['RN'],
                                forward_uniform_DDp_probability=float(p),forward_uniform_DDn_probability=float(nn),
                                forward_path_nm=float(target.path_range(e/1000)),
                                dynamic_to_forward_yield=float(row['yield_total']/(p+nn)),
                                projected_to_bca=float(row['implanted_projected_nm']/ref['Rp']),usage=usage))
    return results


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['pilot','production','refine','references'],required=True)
    parser.add_argument('--nions',type=int,default=48)
    parser.add_argument('--grid',type=int,default=96)
    parser.add_argument('--step',type=float,default=.025)
    parser.add_argument('--output',default=None)
    parser.add_argument('--material',choices=list(MATS))
    parser.add_argument('--axis',choices=['100','110','111'])
    parser.add_argument('--isotope',choices=['D','T'])
    parser.add_argument('--energy',type=float,help='Selected lab keV')
    parser.add_argument('--divergence',type=float)
    parser.add_argument('--density-scale',type=float,default=1.)
    parser.add_argument('--vibration-scale',type=float,default=1.)
    parser.add_argument('--switch-angle',type=float,default=.20)
    args=parser.parse_args(); start=time.monotonic()
    rows=[]; usage=[]
    if args.stage=='references':
        payload=dict(references=references(args.nions,8741))
    else:
        for mi,(name,mat) in enumerate(MATS.items()):
            if args.stage=='pilot' and name!='TiD2': continue
            if args.material and name!=args.material: continue
            for mass in ([MD] if args.stage=='pilot' else [MD,MT]):
                if args.isotope and (mass==MD)!=(args.isotope=='D'): continue
                scatter=Scatter(mat,mass)
                for ai,axis in enumerate(['random','100','110','111']):
                    if args.axis and axis not in ['random',args.axis]: continue
                    cases=[(13000.,0.),(200000.,0.)] if args.stage=='pilot' else None
                    if args.energy is not None or args.divergence is not None:
                        cases=[(e,d) for e in ENERGIES for d in DIVERGENCES
                               if (args.energy is None or e==args.energy*1000)
                               and (args.divergence is None or d==args.divergence)]
                    rr,uu=run_batch(mat,axis,mass,scatter,args.nions,913+mi*1000+ai*100+(mass==MT)*17,
                                    grid=args.grid,step=args.step,only=cases,
                                    density_scale=args.density_scale,vibration_scale=args.vibration_scale,
                                    switch_angle=args.switch_angle)
                    rows+=rr; usage.append(dict(material=name,isotope='D' if mass==MD else 'T',axis=axis,**uu))
        refs={(r['material'],r['isotope'],r['E_lab_keV'],r['divergence_rms_deg']):r for r in rows if r['axis']=='random'}
        for r in rows:
            ref=refs[(r['material'],r['isotope'],r['E_lab_keV'],r['divergence_rms_deg'])]
            for quantity in ['yield_total','energy_return','fuel_flux','stopping_eV_nm','electronic_eV_nm']:
                r[quantity+'_relative_random']=r[quantity]/ref[quantity] if ref[quantity] else None
            # Independent ensembles; this is a delta-method standard error, not a bound.
            r['yield_ratio_se']=0. if r is ref else r['yield_total_relative_random']*np.hypot(r['yield_se']/r['yield_total'],ref['yield_se']/ref['yield_total'])
        payload=dict(rows=rows,usage=usage)
    payload.update(stage=args.stage,nions=args.nions,grid=args.grid,step_nm=args.step,
                   total_seconds=time.monotonic()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   pinned_base='47f090229bc6f02c6273ae5e6ce6e6d66c6b5c26')
    payload['settings']={key:value for key,value in vars(args).items() if key!='output'}
    output=Path(args.output) if args.output else OUT/(args.stage+'.json')
    output.write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    if rows:
        with output.with_suffix('.csv').open('w',newline='') as stream:
            keys=list(dict.fromkeys(k for r in rows for k in r))
            writer=csv.DictWriter(stream,fieldnames=keys); writer.writeheader(); writer.writerows(rows)
    print(json.dumps(dict(output=output.name,rows=len(rows),seconds=payload['total_seconds'],peak_rss_MiB=payload['peak_rss_MiB'])))


if __name__=='__main__':
    main()
