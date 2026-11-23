"""Parity-exact projected cubic-well Hamiltonian and Fourier-Hermite densities.

Finite-basis calculation, not DFT exact. Coordinate powers are projected from
an extended basis; Fourier phases are (-i)^n, physical k in A^-1.
"""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.special import eval_hermite, gammaln, roots_hermite
from scipy.interpolate import RegularGridInterpolator
from model import KB, MASS, kinetic_coefficient, full

PARITIES=((0,0,0),(1,0,0),(1,1,0),(1,1,1))
DEGENERACY=(1,3,3,1)


def hermite(x,n,length,momentum=False):
    q=np.asarray(x)*(length if momentum else 1/length)
    vals=np.empty((len(q),n),complex if momentum else float)
    vals[:,0]=np.pi**(-.25)*np.exp(-q*q/2)
    if n>1: vals[:,1]=np.sqrt(2)*q*vals[:,0]
    for j in range(2,n):
        vals[:,j]=np.sqrt(2/j)*q*vals[:,j-1]-np.sqrt((j-1)/j)*vals[:,j-2]
    vals*=np.sqrt(length) if momentum else 1/np.sqrt(length)
    if momentum:vals*=(-1j)**np.arange(n)
    return vals


def solve(k,a,b,mass,npar):
    c=2*kinetic_coefficient(mass)
    roots=np.roots([1.,0.,-c*k,-(3*a+2*b)*c*c/2])
    quantum=max(v.real for v in roots if abs(v.imag)<1e-9 and v.real>0)
    length=np.sqrt(c/quantum); n=2*npar; nf=n+4
    ladder=np.diag(np.sqrt(np.arange(1,nf)),1)
    Q=length/np.sqrt(2)*(ladder+ladder.T);Q2=Q@Q;Q4=Q2@Q2
    h=np.diag(quantum*(np.arange(nf)+.5))+.5*(k-quantum**2/c)*Q2+a*Q4/4
    I=np.eye(npar);blocks=[]
    for par,multiplicity in zip(PARITIES,DEGENERACY):
        indices=[np.arange(p,n,2) for p in par]
        hs=[h[np.ix_(i,i)] for i in indices];qs=[Q2[np.ix_(i,i)] for i in indices]
        kron=lambda A,B,C:np.kron(np.kron(A,B),C)
        H=kron(hs[0],I,I)+kron(I,hs[1],I)+kron(I,I,hs[2])
        H+=b/2*(kron(qs[0],qs[1],I)+kron(qs[0],I,qs[2])+kron(I,qs[1],qs[2]))
        energies,vectors=eigh(H,overwrite_a=True,check_finite=False)
        blocks.append((par,multiplicity,energies,vectors))
    return dict(length=length,n=n,npar=npar,blocks=blocks,mass=mass,k=k,a=a,b=b)


def tensor_values(coeff,maps):
    return np.einsum('ai,bj,ck,ijks->abcs',*maps,coeff,optimize=True)


def thermal_blocks(spectrum,T):
    E0=spectrum['blocks'][0][2][0]
    Z=sum(mult*np.exp(-(en-E0)/(KB*T)).sum() for _,mult,en,_ in spectrum['blocks'])
    for par,mult,en,v in spectrum['blocks']:
        w=np.exp(-(en-E0)/(KB*T))/Z
        retained=w>1e-13
        yield par,mult,en[retained]-E0,v[:,retained],w[retained]


def density_grid(spectrum,T,grid,momentum=False):
    maps=hermite(grid,spectrum['n'],spectrum['length'],momentum)
    rho=np.zeros((len(grid),)*3)
    omitted=1.
    for par,mult,en,v,w in thermal_blocks(spectrum,T):
        axis_maps=[maps[:,p::2] for p in par]
        representative=np.zeros_like(rho)
        for start in range(0,len(w),8):
            coeff=v[:,start:start+8].reshape((spectrum['npar'],)*3+(-1,))
            vals=tensor_values(coeff,axis_maps)
            representative+=np.einsum('ijks,s->ijk',abs(vals)**2,w[start:start+8])
        if mult==3:
            rho+=representative+representative.transpose(1,2,0)+representative.transpose(2,0,1)
        else:rho+=representative
        omitted-=mult*w.sum()
    return rho,max(0.,float(omitted))


def distributions(grid,rho,support,radii=None,angular_order=16):
    if radii is None:radii=np.linspace(0,float(grid[-1]),181)
    mu,muw=np.polynomial.legendre.leggauss(angular_order);nphi=2*angular_order;phi=np.arange(nphi)*2*np.pi/nphi
    directions=np.stack(np.broadcast_arrays(np.sqrt(1-mu[:,None]**2)*np.cos(phi),
                           np.sqrt(1-mu[:,None]**2)*np.sin(phi),mu[:,None]+phi*0),axis=-1).reshape(-1,3)
    aw=np.repeat(muw,nphi)/(2*nphi)
    interp=RegularGridInterpolator((grid,)*3,rho,bounds_error=False,fill_value=0.)
    sampled=interp((radii[:,None,None]*directions[None,:,:]).reshape(-1,3)).reshape(len(radii),-1)
    density=sampled@aw; measure=4*np.pi*radii*radii*density
    tail=np.zeros(len(radii));tail[:-1]=np.cumsum(((measure[:-1]+measure[1:])/2*np.diff(radii))[::-1])[::-1]
    cutoff=lambda y:float(radii[np.flatnonzero((radii>=radii[np.argmax(y)])&(y<=y.max()*1e-6))[0]]) if np.any((radii>=radii[np.argmax(y)])&(y<=y.max()*1e-6)) else None
    norm=float(np.trapezoid(measure,radii));r2=float(np.trapezoid(measure*radii*radii,radii))
    support_fraction=float(np.interp(support,radii,norm-tail)) if support is not None else None
    return dict(radius=radii,density=density,measure=measure,tail=tail),dict(
        normalization=norm,mean_square_radius=r2,density_peak_radius=float(radii[np.argmax(density)]),
        radial_measure_peak_radius=float(radii[np.argmax(measure)]),
        density_cutoff_relative_1e6=cutoff(density),radial_measure_cutoff_relative_1e6=cutoff(measure),
        support_fraction=support_fraction,exterior_fraction=None if support_fraction is None else norm-support_fraction,
        grid_face_density_over_peak=float(rho[-1].max()/rho.max()))


def signed_marginal(grid,rho,direction):
    n=np.asarray(direction,float);n/=np.linalg.norm(n)
    t=np.cross(n,[0,1,0]);t/=np.linalg.norm(t);u=np.cross(n,t)
    s=np.linspace(-.9,.9,181);trans=np.linspace(-.9,.9,101)
    aa,bb=np.meshgrid(trans,trans,indexing='ij')
    interp=RegularGridInterpolator((grid,)*3,rho,bounds_error=False,fill_value=0.)
    marginal=[];line=[]
    for z in s:
        points=z*n+aa[...,None]*t+bb[...,None]*u
        vals=interp(points.reshape(-1,3)).reshape(aa.shape)
        marginal.append(np.trapezoid(np.trapezoid(vals,trans,axis=1),trans))
        line.append(float(interp((z*n)[None,:])[0]))
    return s,np.array(marginal),np.array(line)


def joint_states(spectrum,T,nbins=40):
    order=spectrum['n']+6;q,weights=roots_hermite(order)
    maps=hermite(q/spectrum['length'],spectrum['n'],spectrum['length'],True)
    quadrature=(weights*np.exp(q*q)/spectrum['length'])
    measure=quadrature[:,None,None]*quadrature[None,:,None]*quadrature[None,None,:]
    kinetic=kinetic_coefficient(spectrum['mass'])/spectrum['length']**2*(q[:,None,None]**2+q[None,:,None]**2+q[None,None,:]**2)
    edges=np.linspace(0,float(kinetic.max())+1e-12,nbins+1)
    binid=np.clip(np.digitize(kinetic.ravel(),edges)-1,0,nbins-1)
    states=[]
    for par,mult,en,v,w in thermal_blocks(spectrum,T):
        for start in range(0,len(w),8):
            values=tensor_values(v[:,start:start+8].reshape((spectrum['npar'],)*3+(-1,)),[maps[:,p::2] for p in par])
            for j in range(values.shape[-1]):
                prob=(abs(values[...,j])**2*measure).ravel()
                p=np.bincount(binid,weights=prob,minlength=nbins)
                pk=np.bincount(binid,weights=prob*kinetic.ravel(),minlength=nbins)
                active=p>1e-16
                states.append(dict(excitation_eV=float(en[start+j]),kinetic_eV=(pk[active]/p[active]).tolist(),
                                   probabilities=(mult*w[start+j]*p[active]).tolist()))
    total=sum(sum(row['probabilities']) for row in states)
    for row in states:row['probabilities']=(np.array(row['probabilities'])/total).tolist()
    return states,dict(quadrature_order=order,probability_before_roundoff_normalization=total,
                       kinetic_bin_width_eV=float(edges[1]-edges[0]),isotropization='spherical radial folding, angular correlations not supplied')


def run(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True);summary=[]
    for material,isotope in [('Pd','H'),('Pd','D'),('Pd','T'),('Ni','H')]:
        for variant in ('joint','energy','force'):
            m=full(material,variant);spec=solve(m['K'][0,0,0],m['a'],m['b'],MASS[isotope],10)
            for T in (300,77):
                name=f'{material}{isotope}_{T}K_{variant}'
                x=np.linspace(-1.,1.,101);k=np.linspace(-42.,42.,121)
                rho,missing=density_grid(spec,T,x);nk,_=density_grid(spec,T,k,True)
                radial,rs=distributions(x,rho,.24);momentum,ms=distributions(k,nk,None,np.linspace(0,42,211))
                if material=='Ni' and T==77 and variant=='joint':
                    angular={}
                    for kind,grid,density,curve,statistics,support in (
                            ('position',x,rho,radial,rs,.24),('momentum',k,nk,momentum,ms,None)):
                        finer,fs=distributions(grid,density,support,curve['radius'],angular_order=32)
                        i=np.argmin(abs(curve['radius']-statistics['density_cutoff_relative_1e6']))
                        angular[kind]=dict(case=name,angular_rule='GL16 x phi32 versus GL32 x phi64 on unchanged density grid',
                            relative_density_change_at_1e6=float(abs(finer['density'][i]-curve['density'][i])/finer['density'][i]),
                            cutoff_coarse=statistics['density_cutoff_relative_1e6'],cutoff_fine=fs['density_cutoff_relative_1e6'],
                            normalization_change=fs['normalization']-statistics['normalization'])
                    (output.parent/'angular_sensitivity.json').write_text(json.dumps(angular,indent=2)+'\n')
                np.savez_compressed(output/(name+'.npz'),x_A=x,rho_A_minus3=rho,k_A_inverse=k,n_k_A3=nk,
                                    radial_A=radial['radius'],radial_density_A_minus3=radial['density'],radial_measure_A_minus1=radial['measure'],radial_tail=radial['tail'],
                                    momentum_k_A_inverse=momentum['radius'],momentum_density_A3=momentum['density'],momentum_measure_A=momentum['measure'],momentum_tail=momentum['tail'])
                row=dict(case=name,material=material,isotope=isotope,temperature_K=T,variant=variant,npar=10,
                         length_A=spec['length'],thermal_omitted_probability=missing,position=rs,momentum=ms,
                         mean_kinetic_eV=kinetic_coefficient(MASS[isotope])*ms['mean_square_radius'])
                if variant=='joint':
                    for direction,vec in [('100',[1,0,0]),('110',[1,1,0]),('111',[1,1,1])]:
                        s,p,line=signed_marginal(x,rho,vec)
                        np.savetxt(output/f'{name}_normal{direction}.csv',np.column_stack([s,p,line]),delimiter=',',header='s_A,marginal_A_inverse,on_axis_line_density_A_minus3',comments='')
                    if material=='Pd' and isotope=='D':
                        states,jmeta=joint_states(spec,T)
                        (output/f'{name}_g8_joint.json').write_text(json.dumps(dict(states=states,metadata=jmeta),separators=(',',':'))+'\n')
                summary.append(row)
                print(name,rs['mean_square_radius'],ms['density_cutoff_relative_1e6'],flush=True)
            if variant=='joint':
                coarse=solve(m['K'][0,0,0],m['a'],m['b'],MASS[isotope],8)
                for T in (300,77):
                    nk,_=density_grid(coarse,T,k,True);rho,_=density_grid(coarse,T,x)
                    rr,rs=distributions(x,rho,.24);mm,ms=distributions(k,nk,None,np.linspace(0,42,211))
                    np.savez_compressed(output/f'{material}{isotope}_{T}K_basis8.npz',radial_density_A_minus3=rr['density'],momentum_density_A3=mm['density'])
                    summary.append(dict(case=f'{material}{isotope}_{T}K_basis8',position=rs,momentum=ms))
    (output/'exact_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    return summary
