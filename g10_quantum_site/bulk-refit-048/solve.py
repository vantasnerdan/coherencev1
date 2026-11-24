"""Matched 3D Cartesian HO Galerkin calculation with nonseparable quadrature.

No Gaussian density approximation: all saved distributions are canonical sums
of squared position/Fourier eigenfunctions. Parity/permutation sectors complete.
"""
import argparse
import itertools
import json
from pathlib import Path
import time
import numpy as np
from scipy.constants import physical_constants,hbar,electron_volt,atomic_mass
from scipy.linalg import eigh
from scipy.special import roots_hermite
from scipy.interpolate import RegularGridInterpolator
from well import HERE,K,A,B,VARIANTS,energy_force

KB=physical_constants['Boltzmann constant in eV/K'][0]
MASS={s:physical_constants[n+' mass'][0]/atomic_mass for s,n in [('H','proton'),('D','deuteron'),('T','triton')]}
PARITIES=((0,0,0),(1,0,0),(1,1,0),(1,1,1))
MULT=(1,3,3,1)
OUT=HERE/'results'


def kinetic(mass):
    return hbar*hbar/(2*mass*atomic_mass*electron_volt*1e-20)


def hermite(x,n,length,momentum=False):
    q=np.asarray(x)*(length if momentum else 1/length)
    out=np.empty((len(q),n),complex if momentum else float)
    out[:,0]=np.pi**(-.25)*np.exp(-q*q/2)
    if n>1:out[:,1]=np.sqrt(2)*q*out[:,0]
    for j in range(2,n):out[:,j]=np.sqrt(2/j)*q*out[:,j-1]-np.sqrt((j-1)/j)*out[:,j-2]
    out*=np.sqrt(length) if momentum else 1/np.sqrt(length)
    if momentum:out*=(-1j)**np.arange(n)
    return out


def length_scale(mass):
    c=2*kinetic(mass)
    roots=np.roots([1.,0.,-c*K,-(3*A+2*B)*c*c/2])
    quantum=max(v.real for v in roots if abs(v.imag)<1e-9 and v.real>0)
    return np.sqrt(c/quantum)


def solve(isotope,variant,npar,order):
    mass=MASS[isotope];length=length_scale(mass);n=2*npar
    q,weights=roots_hermite(order)
    positive=q>0;q=q[positive];weights=weights[positive]
    # Full even products reduce exactly to 8 times positive-octant rule.
    maps=hermite(q,n,1.)*np.sqrt(2*weights*np.exp(q*q))[:,None]
    xx=np.stack(np.meshgrid(q*length,q*length,q*length,indexing='ij'),axis=-1)
    V=energy_force(xx,variant)[0].ravel()
    ladder=np.diag(np.sqrt(np.arange(1,n+2)),1)
    Q=(ladder+ladder.T)/np.sqrt(2)
    T=(np.diag(2*np.arange(n+2)+1)-Q@Q)*kinetic(mass)/length**2
    I=np.eye(npar);blocks=[];start=time.monotonic()
    for par,mult in zip(PARITIES,MULT):
        indices=[np.arange(p,n,2) for p in par]
        bs=[maps[:,i] for i in indices]
        design=np.einsum('ai,bj,ck->abcijk',*bs,optimize=True).reshape(len(V),-1)
        H=design.T@(V[:,None]*design)
        ts=[T[np.ix_(i,i)] for i in indices]
        H+=np.kron(np.kron(ts[0],I),I)+np.kron(np.kron(I,ts[1]),I)+np.kron(np.kron(I,I),ts[2])
        e,v=eigh(H,subset_by_value=(-np.inf,1.05),check_finite=False,overwrite_a=True)
        # A matrix residual is a solve diagnostic, not continuum error.
        blocks.append(dict(par=par,mult=mult,energies=e,vectors=v))
    print('spectrum',isotope,variant,npar,order,round(time.monotonic()-start,2),flush=True)
    return dict(isotope=isotope,variant=variant,npar=npar,order=order,length=length,mass=mass,blocks=blocks)


def save_spectrum(s,path):
    arrays={k:s[k] for k in ('isotope','variant','npar','order','length','mass')}
    for i,b in enumerate(s['blocks']):
        arrays[f'e{i}']=b['energies'];arrays[f'v{i}']=b['vectors']
    np.savez_compressed(path,**arrays)


def load_spectrum(path):
    z=np.load(path);s={k:z[k].item() for k in ('isotope','variant','npar','order','length','mass')}
    s['blocks']=[dict(par=p,mult=m,energies=z[f'e{i}'],vectors=z[f'v{i}']) for i,(p,m) in enumerate(zip(PARITIES,MULT))]
    return s


def thermal(s,T):
    e0=s['blocks'][0]['energies'][0]
    Z=sum(b['mult']*np.exp(-(b['energies']-e0)/(KB*T)).sum() for b in s['blocks'])
    for b in s['blocks']:
        w=np.exp(-(b['energies']-e0)/(KB*T))/Z
        active=w>2e-13
        yield b['par'],b['mult'],b['energies'][active]-e0,b['vectors'][:,active],w[active]


def values(coeff,maps):
    return np.einsum('ai,bj,ck,ijks->abcs',*maps,coeff,optimize=True)


def density(s,T,grid,momentum=False):
    maps=hermite(grid,2*s['npar'],s['length'],momentum)
    rho=np.zeros((len(grid),)*3);omitted=1.
    for par,mult,en,v,w in thermal(s,T):
        representative=np.zeros_like(rho)
        for start in range(0,len(w),8):
            c=v[:,start:start+8].reshape((s['npar'],)*3+(-1,))
            val=values(c,[maps[:,p::2] for p in par])
            representative+=np.einsum('ijks,s->ijk',abs(val)**2,w[start:start+8])
        rho+=representative
        if mult==3:rho+=representative.transpose(1,2,0)+representative.transpose(2,0,1)
        omitted-=mult*w.sum()
    return rho,max(0.,float(omitted))


def point_density(s,T,points,momentum=False):
    shape=points.shape[:-1];points=points.reshape(-1,3)
    maps=[hermite(points[:,i],2*s['npar'],s['length'],momentum) for i in range(3)]
    rho=np.zeros(len(points))
    for par,mult,en,v,w in thermal(s,T):
        for perm in sorted(set(itertools.permutations(range(3))) if mult==3 else {(0,1,2)}):
            # Only three distinct permutations of a parity representative.
            permpar=tuple(par[i] for i in perm)
            if mult==3 and perm!=min(p for p in itertools.permutations(range(3)) if tuple(par[i] for i in p)==permpar):continue
            coeff=v.reshape((s['npar'],)*3+(-1,)).transpose(tuple(perm)+(3,))
            for start in range(0,len(w),16):
                val=np.einsum('ai,aj,ak,ijks->as',*[maps[i][:,p::2] for i,p in enumerate(permpar)],coeff[...,start:start+16],optimize=True)
                rho+=np.einsum('as,s->a',abs(val)**2,w[start:start+16])
    return rho.reshape(shape)


def radial_distribution(s,T,r,grid,rho,momentum=False,angular=24):
    from scipy.ndimage import map_coordinates,spline_filter
    mu,muw=np.polynomial.legendre.leggauss(angular);phi=np.arange(2*angular)*np.pi/angular
    u=np.stack(np.broadcast_arrays(np.sqrt(1-mu[:,None]**2)*np.cos(phi),np.sqrt(1-mu[:,None]**2)*np.sin(phi),mu[:,None]+phi*0),axis=-1).reshape(-1,3)
    w=np.repeat(muw,2*angular)/(4*angular)
    # Tensor eigenfunction contraction is much cheaper than pointwise 3D
    # contraction. Even cubic-spline interpolation is checked at finer pitch.
    coefficients=spline_filter(rho,order=3,mode='mirror')
    points=abs(r[:,None,None]*u[None,:,:])
    vals=map_coordinates(coefficients,(points/(grid[1]-grid[0])).reshape(-1,3).T,order=3,mode='mirror',prefilter=False).reshape(len(r),-1)
    d=vals@w;measure=4*np.pi*r*r*d
    tail=np.zeros_like(r);tail[:-1]=np.cumsum(((measure[:-1]+measure[1:])/2*np.diff(r))[::-1])[::-1]
    directions=np.array([[1.,0,0],[1.,1,0],[1.,1,1]]);directions/=np.linalg.norm(directions,axis=1)[:,None]
    lines=point_density(s,T,r[:,None,None]*directions[None,:,:],momentum)
    return dict(radius=r,density=d,measure=measure,tail=tail,lines=lines)


def moments(s,T):
    n=2*s['npar'];ladder=np.diag(np.sqrt(np.arange(1,n+2)),1);Q=(ladder+ladder.T)/np.sqrt(2)
    Q2=Q@Q;P2=np.diag(2*np.arange(n+2)+1)-Q2
    r2=k2=0.
    for par,mult,en,v,w in thermal(s,T):
        c=v.reshape((s['npar'],)*3+(-1,))
        for axis,p in enumerate(par):
            idx=np.arange(p,n,2)
            for mat,ispos in [(Q2,True),(P2,False)]:
                transformed=np.moveaxis(np.tensordot(mat[np.ix_(idx,idx)],c,axes=(1,axis)),0,axis)
                val=np.einsum('ijks,ijks,s->',c,transformed,w)*mult
                if ispos:r2+=val*s['length']**2
                else:k2+=val/s['length']**2
    return float(r2),float(k2*kinetic(s['mass']))


def summary(s,T,rd,md):
    r2,ke=moments(s,T)
    def stat(d):
        r,y,m=d['radius'],d['density'],d['measure']
        def cutoff(a):
            hit=np.flatnonzero((np.arange(len(a))>=np.argmax(a))&(a<=a.max()*1e-6))
            return None if not len(hit) else float(r[hit[0]])
        return dict(normalization=float(np.trapezoid(m,r)),density_peak=float(r[np.argmax(y)]),measure_peak=float(r[np.argmax(m)]),candidate_density_1e6=cutoff(y),candidate_measure_1e6=cutoff(m))
    return dict(isotope=s['isotope'],temperature_K=T,variant=s['variant'],npar=s['npar'],quadrature_order=s['order'],length_A=s['length'],ground_energy_eV=float(s['blocks'][0]['energies'][0]),mean_r2_A2=r2,mean_kinetic_eV=ke,position=stat(rd),momentum=stat(md),exterior_024=float(np.interp(.24,rd['radius'],rd['tail'])),exterior_048=float(np.interp(.48,rd['radius'],rd['tail'])),thermal_retained_states=sum(mult*len(w) for par,mult,en,v,w in thermal(s,T)))


def analyze(s,tag,full=True):
    rows=[]
    for T in (300,77):
        r=np.linspace(0,1.6,321);k=np.linspace(0,48,241)
        size=101 if tag=='check' else 81
        angular=32 if tag=='check' else 24
        x=np.linspace(0,1.6,size);kg=np.linspace(0,48,size)
        rho,omit=density(s,T,x);nk,komit=density(s,T,kg,True)
        rd=radial_distribution(s,T,r,x,rho,angular=angular)
        md=radial_distribution(s,T,k,kg,nk,True,angular)
        row=summary(s,T,rd,md);rows.append(row)
        row.update(position_grid_pitch_A=float(x[1]),momentum_grid_pitch_A_inverse=float(kg[1]),angular_GL_order=angular,discarded_thermal_weight=max(omit,komit),finite_basis_energy_ceiling_thermal_weight_bound=float(8*s['npar']**3*np.exp(-(1.05-row['ground_energy_eV'])/(KB*T))))
        arrays={f'position_{key}':val for key,val in rd.items()};arrays.update({f'momentum_{key}':val for key,val in md.items()})
        if full:
            arrays.update(x_positive_A=x,rho_positive_octant_A_minus3=rho,k_positive_A_inverse=kg,n_positive_octant_A3=nk)
            # Signed [111] plane marginal: integrate both transverse axes.
            n=np.ones(3)/np.sqrt(3);a=np.array([1.,-1,0])/np.sqrt(2);b=np.cross(n,a)
            signed=np.linspace(-1.2,1.2,241);trans=np.linspace(-1.4,1.4,141)
            aa,bb=np.meshgrid(trans,trans,indexing='ij');interp=RegularGridInterpolator((x,)*3,rho,bounds_error=False,fill_value=0.)
            marginal=[]
            for z in signed:
                vals=interp(abs(z*n+aa[...,None]*a+bb[...,None]*b))
                marginal.append(np.trapezoid(np.trapezoid(vals,trans,axis=1),trans))
            arrays.update(s111_A=signed,marginal111_A_inverse=marginal)
            row['marginal111_normalization']=float(np.trapezoid(marginal,signed))
        np.savez_compressed(OUT/f'{s["isotope"]}_{T}_{s["variant"]}_{tag}_distributions.npz',**arrays)
        print('observables',s['isotope'],T,s['variant'],row['mean_r2_A2'],row['mean_kinetic_eV'],flush=True)
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('--isotopes',nargs='+',default=['D']);p.add_argument('--variants',nargs='+',default=list(VARIANTS));p.add_argument('--npar',type=int,default=10);p.add_argument('--order',type=int,default=40);p.add_argument('--tag',default='main');p.add_argument('--compact',action='store_true');args=p.parse_args()
    OUT.mkdir(exist_ok=True);rows=[]
    for iso in args.isotopes:
        for var in args.variants:
            path=OUT/f'{iso}_{var}_{args.tag}_spectrum.npz'
            s=solve(iso,var,args.npar,args.order);save_spectrum(s,path)
            rows.extend(analyze(s,args.tag,not args.compact))
            (OUT/f'summary_{args.tag}.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


if __name__=='__main__':main()
