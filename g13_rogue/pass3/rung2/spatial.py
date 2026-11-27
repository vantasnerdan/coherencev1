"""Globally coupled Q4 x Q4 x P4 variational foil; no local wave patches."""
from pathlib import Path
import json
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, splu
from numpy.polynomial.legendre import leggauss, legval, legder

HERE=Path(__file__).resolve().parent
CELLS=16
PX=4
PZ=4
_NODES=np.linspace(-1,1,PX+1)
_COEFF=np.array([np.polynomial.polynomial.polyfromroots(np.delete(_NODES,i))/np.prod(node-np.delete(_NODES,i)) for i,node in enumerate(_NODES)])
_DER=np.polynomial.polynomial.polyder(_COEFF,axis=1)

def axis(length, cells=CELLS):
    """Exact consistent mass and derivative bilinear forms, split by old grain."""
    x,w=leggauss(PX+1)
    vals=np.polynomial.polynomial.polyval(x,_COEFF.T).T
    ders=np.polynomial.polynomial.polyval(x,_DER.T).T*2*cells/length
    out=[]
    for grain in range(8):
        forms=[]
        for a in range(2):
            row=[]
            for b in range(2):
                local=(ders if a else vals).T@((w*length/(2*cells))[:,None]*(ders if b else vals))
                rr=[];cc=[];vv=[]
                for e in range(grain*cells//8,(grain+1)*cells//8):
                    ids=np.arange(PX*e,PX*e+PX+1)
                    rr.extend(np.repeat(ids,PX+1));cc.extend(np.tile(ids,PX+1));vv.extend(local.ravel())
                row.append(sparse.coo_matrix((vv,(rr,cc)),shape=(PX*cells+1,)*2).tocsr())
            forms.append(row)
        out.append(forms)
    return out

def zforms(h):
    z,w=leggauss(PZ+2);v=[];d=[]
    for j in range(PZ+1):
        c=np.zeros(j+1);c[j]=1
        v.append(legval(z,c)*np.sqrt((2*j+1)/h))
        d.append(legval(z,legder(c))*np.sqrt((2*j+1)/h)*2/h)
    v=np.array(v).T;d=np.array(d).T
    return [[sparse.csr_matrix(a.T@((w*h/2)[:,None]*b)) for b in (v,d)] for a in (v,d)]

def maps(points,v,dims):
    """Full modal displacement and engineering strain at physical points."""
    nxy=PX*CELLS+1;n=nxy*nxy*(PZ+1);nm=v.shape[1]
    out=np.zeros((len(points),3,nm));strain=np.zeros((len(points),6,nm))
    local_fields={}
    for k,p in enumerate(points):
        ind=[];val=[];der=[]
        for a in range(2):
            e=min(CELLS-1,max(0,int(p[a]*CELLS/dims[a])));s=2*(p[a]*CELLS/dims[a]-e)-1
            ind.append(np.arange(PX*e,PX*e+PX+1));val.append(np.polynomial.polynomial.polyval(s,_COEFF.T))
            der.append(np.polynomial.polynomial.polyval(s,_DER.T)*2*CELLS/dims[a])
        zz=[];zd=[]
        for j in range(PZ+1):
            c=np.zeros(j+1);c[j]=1;scale=np.sqrt((2*j+1)/dims[2]);s=2*p[2]/dims[2]-1
            zz.append(legval(s,c)*scale);zd.append(legval(s,legder(c))*scale*2/dims[2])
        ids=((ind[0][:,None,None]*nxy+ind[1][None,:,None])*(PZ+1)+np.arange(PZ+1)[None,None,:]).ravel()
        f=np.einsum('i,j,k->ijk',val[0],val[1],zz).ravel()
        df=[np.einsum('i,j,k->ijk',der[0],val[1],zz).ravel(),np.einsum('i,j,k->ijk',val[0],der[1],zz).ravel(),np.einsum('i,j,k->ijk',val[0],val[1],zd).ravel()]
        key=(int(ind[0][0]//PX),int(ind[1][0]//PX))
        block=local_fields.get(key)
        if block is None:
            block=v[np.arange(3)[:,None]*n+ids[None,:]]
            local_fields[key]=block
        grad=np.einsum('ai,cim->cam',np.array(df),block)
        out[k]=np.einsum('i,cim->cm',f,block)
        strain[k,:3]=np.array([grad[j,j] for j in range(3)])
        strain[k,3:]=np.array([grad[1,2]+grad[2,1],grad[0,2]+grad[2,0],grad[0,1]+grad[1,0]])
    return out,strain

def projection_capture(length,wavelength):
    # Consistent L2 projection onto the actual conforming lateral space.
    forms=axis(length);mass=sum(g[0][0] for g in forms).tocsc();solve=splu(mass).solve
    x,w=leggauss(12);n=PX*CELLS+1;theta=np.linspace(0,np.pi/2,37)
    scores=[]
    for angle in theta:
        factors=[]
        for wave in (2*np.pi/wavelength*np.cos(angle),2*np.pi/wavelength*np.sin(angle)):
            rhs=np.zeros(n,dtype=complex)
            for e in range(CELLS):
                xx=(e+(x+1)/2)*length/CELLS;phi=np.polynomial.polynomial.polyval(x,_COEFF.T)
                rhs[PX*e:PX*e+PX+1]+=phi@(w*length/(2*CELLS)*np.exp(1j*wave*xx))
            factors.append(solve(rhs.real)+1j*solve(rhs.imag))
        a,b=factors;z=np.outer(a,b)
        for phase in (0,np.pi/4,np.pi/2):
            real=np.real(np.exp(1j*phase)*z)
            norm=.5*length**2*(1+np.real(np.exp(2j*phase+1j*2*np.pi/wavelength*length*(np.cos(angle)+np.sin(angle)))*np.sinc(2*length/wavelength*np.cos(angle))*np.sinc(2*length/wavelength*np.sin(angle))))
            scores.append(np.sum(real*(mass@real@mass))/norm)
    return float(min(scores))

def branch_comparison(c,mass,v,omega):
    """Compare unit-mass modes, not unnormalised nodal interpolants."""
    nxy=PX*CELLS+1;n=nxy*nxy*(PZ+1)
    old=np.load(Path(__file__).resolve().parent.parent/'rung1'/'model.npz',allow_pickle=False)
    nodes=np.linspace(-1,1,nxy);pol=np.polynomial.legendre.legvander(nodes,15)*np.sqrt((2*np.arange(16)+1)/c.DIMS[0])
    embedded=np.einsum('ia,jb,cabzm->cijzm',pol,pol,old['eigenvectors'].reshape(3,16,16,5,-1),optimize=True).reshape(3*n,-1)/np.sqrt(c.RHO)
    weighted=mass@embedded;norm=np.einsum('im,im->m',embedded,weighted)
    overlap=(v.T@weighted)/np.sqrt(norm)[None,:];matches=np.argmax(abs(overlap),axis=0)
    np.savez_compressed(HERE/'branch_comparison.npz',old_frequency_Hz=old['omega']/(2*np.pi),new_frequency_Hz=omega[matches]/(2*np.pi),mass_overlap_squared=overlap[matches,np.arange(len(matches))]**2,old_interpolant_mass_norm_squared=norm)

def prepare(c):
    from scipy.spatial.transform import Rotation
    x=axis(c.DIMS[0]);y=axis(c.DIMS[1]);z=zforms(c.DIMS[2]);nxy=PX*CELLS+1;n=nxy*nxy*(PZ+1)
    tensors=np.array([c.cubic(r) for r in Rotation.random(64,random_state=173).as_matrix()]).reshape(8,8,6,6)
    smap=np.array([[0,5,4],[5,1,3],[4,3,2]])
    blocks=[]
    for i in range(3):
        row=[]
        for j in range(3):
            block=sparse.csr_matrix((n,n))
            for a in range(3):
                for b in range(3):
                    xy=sparse.csr_matrix((nxy*nxy,)*2)
                    for g in range(8):
                        for h in range(8):
                            xy+=tensors[g,h,smap[i,a],smap[j,b]]*sparse.kron(x[g][int(a==0)][int(b==0)],y[h][int(a==1)][int(b==1)],format='csr')
                    block+=sparse.kron(xy,z[int(a==2)][int(b==2)],format='csr')
            row.append(block)
        blocks.append(row)
    k=sparse.bmat(blocks,format='csc');k=(k+k.T)*.5
    mx=sum(g[0][0] for g in x);my=sum(g[0][0] for g in y)
    scalar=sparse.kron(sparse.kron(mx,my),z[0][0],format='csc')*c.RHO
    mass=sparse.block_diag([scalar]*3,format='csc')
    print(json.dumps(dict(stage='assembled',ndof=k.shape[0],nnz=k.nnz)),flush=True)
    # One bounded shift-invert extraction of the entire old physical band.
    sparse.save_npz(HERE/'stiffness.npz',k);sparse.save_npz(HERE/'mass.npz',mass)
    val,v=eigsh(k,k=400,M=mass,sigma=(2*np.pi*1000)**2,which='LM',tol=2e-8,ncv=820)
    return finish(c,k,mass,val,v,(x,y,z,scalar,mx,my))

def finish(c,k,mass,val,v,geometry):
    """Freeze the solved mode epoch before deriving consumers or diagnostics."""
    x,y,z,scalar,mx,my=geometry
    nxy=PX*CELLS+1;n=nxy*nxy*(PZ+1)
    order=np.argsort(val);val=val[order];v=v[:,order]
    freq=np.sqrt(np.maximum(val,0))/(2*np.pi);keep=(freq>=5000)&(freq<=250000)
    np.savez_compressed(HERE/'extraction.npz',frequency_Hz=freq)
    print(json.dumps(dict(stage='extracted',highest_frequency_Hz=float(freq[-1]),old_band_modes=int(keep.sum()))),flush=True)
    if freq[-1]<=250000:raise ValueError('Extraction did not bracket full physical band')
    val=val[keep];v=v[:,keep];omega=np.sqrt(val)
    np.savez_compressed(HERE/'model.npz',omega=omega,eigenvectors=v)
    residual=np.linalg.norm(k@v-(mass@v)*val)/np.linalg.norm((mass@v)*val)
    merror=np.max(abs(v.T@(mass@v)-np.eye(len(val))))
    target=.01373147/9;capture=projection_capture(c.DIMS[0],target)
    if capture<.95:raise AssertionError('Actual global conforming space fails shared screen')
    branch_comparison(c,mass,v,omega)
    # Same branch/group-speed loss proxy as rung1, using consistent global forms.
    xx=[[sum(g[a][b] for g in x) for b in range(2)] for a in range(2)]
    yy=[[sum(g[a][b] for g in y) for b in range(2)] for a in range(2)]
    grad=sparse.kron(sparse.kron(xx[1][1],my),z[0][0])+sparse.kron(sparse.kron(mx,yy[1][1]),z[0][0])
    vz=v[2*n:];massz=np.sum(vz*(scalar@vz),axis=0)
    kz=np.sqrt(c.RHO*np.sum(vz*(grad@vz),axis=0)/np.maximum(massz,1e-300))
    divergence=np.zeros(len(omega))
    for i in range(3):
        for j in range(3):
            op=sparse.kron(sparse.kron(xx[int(i==0)][int(j==0)],yy[int(i==1)][int(j==1)]),z[int(i==2)][int(j==2)])
            divergence+=np.sum(v[i*n:(i+1)*n]*(op@v[j*n:(j+1)*n]),axis=0)
    bulk=(c.C11+2*c.C12)/3
    speed=speed_proxy(c,omega,massz,kz,bulk*divergence/omega**2)
    np.savez_compressed(HERE/'mode_metrics.npz',speed=speed,vertical_mass=massz,vertical_k_per_m=kz,bulk_energy_fraction=bulk*divergence/omega**2)
    result=dict(space='global conforming Q4 x Q4 x P4',lateral_cells=CELLS,lateral_degree=PX,lateral_node_spacing_m=float(c.DIMS[0]/(PX*CELLS)),thickness_degree=PZ,ndof=k.shape[0],nnz=k.nnz,modes=len(val),operator_residual=float(residual),mass_error=float(merror),projection_target_m=target,min_sampled_mass_capture=capture,extraction_bracket_frequency_Hz=float(freq[-1]),free_faces='natural variational',grain_cells_per_axis=8,disjoint_band_added=False)
    (HERE/'support.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

def lateral_basis(points,length):
    """Nodal values and physical derivatives, including shared interface nodes."""
    n=PX*CELLS+1;phi=np.zeros((len(points),n));der=np.zeros_like(phi)
    for i,p in enumerate(points):
        e=min(CELLS-1,max(0,int(p*CELLS/length)));s=2*(p*CELLS/length-e)-1
        phi[i,PX*e:PX*e+PX+1]=np.polynomial.polynomial.polyval(s,_COEFF.T)
        der[i,PX*e:PX*e+PX+1]=np.polynomial.polynomial.polyval(s,_DER.T)*2*CELLS/length
    return phi,der

def eigenstrain_load(dims,temperature,power_over_conductivity):
    """Geometric integral of div(test) times the stated steady T(x,y,z)."""
    q,w=leggauss(PX//2+1);axis=[];weights=[]
    # Resolve every bilinear cell of the common 129-square temperature field.
    for length in dims[:2]:
        edges=np.linspace(0,length,129)
        axis.append(np.concatenate([(a+b)/2+(b-a)*q/2 for a,b in zip(edges[:-1],edges[1:])]))
        weights.append(np.tile(w*length/256,128))
    xx,yy=np.meshgrid(*axis,indexing='ij');xy=np.column_stack([xx.ravel(),yy.ravel()])
    x,dx=lateral_basis(axis[0],dims[0]);y,dy=lateral_basis(axis[1],dims[1])
    lateral=[]
    for field in (temperature,power_over_conductivity):
        f=field(xy).reshape(xx.shape)*weights[0][:,None]*weights[1][None,:]
        lateral.append((dx.T@f@y,x.T@f@dy,x.T@f@y))
    zz,ww=leggauss(5);h=dims[2];z=(zz+1)*h/2;ww*=h/2
    value=[];derivative=[]
    for j in range(PZ+1):
        coef=np.zeros(j+1);coef[j]=1;scale=np.sqrt((2*j+1)/h)
        value.append(legval(zz,coef)*scale);derivative.append(legval(zz,legder(coef))*scale*2/h)
    value=np.array(value);derivative=np.array(derivative);shape=z*z/(2*h)-z+h/3
    zvalue=[value@ww,value@(ww*shape)];zder=[derivative@ww,derivative@(ww*shape)]
    out=[]
    for comp in range(3):
        block=sum(lateral[j][comp][:,:,None]*(zder[j] if comp==2 else zvalue[j]) for j in range(2))
        out.append(block.ravel())
    return np.concatenate(out)

def pressure_load(dims,points,weights):
    x,_=lateral_basis(points[:,0],dims[0]);y,_=lateral_basis(points[:,1],dims[1])
    xy=np.einsum('qi,q,qj->ij',x,weights,y)
    front=np.array([(-1)**j*np.sqrt((2*j+1)/dims[2]) for j in range(PZ+1)])
    n=(PX*CELLS+1)**2*(PZ+1)
    return np.r_[np.zeros(2*n),(xy[:,:,None]*front).ravel()]

class StaticSolver:
    """Six gauge conditions remove only rigid freedom; physical force is balanced."""
    def __init__(self,dims):
        self.k=sparse.load_npz(HERE/'stiffness.npz').tocsc()
        self.mass=sparse.load_npz(HERE/'mass.npz').tocsr()
        nxy=PX*CELLS+1;h=dims[2];n=nxy*nxy*(PZ+1)
        rigid=np.zeros((3,nxy,nxy,PZ+1,6));root=np.sqrt(h)
        for comp in range(3):rigid[comp,:,:,0,comp]=root
        x,y=np.meshgrid(np.linspace(-dims[0]/2,dims[0]/2,nxy),np.linspace(-dims[1]/2,dims[1]/2,nxy),indexing='ij')
        rigid[1,:,:,1,3]=-h**1.5/np.sqrt(12);rigid[2,:,:,0,3]=y*root
        rigid[0,:,:,1,4]=h**1.5/np.sqrt(12);rigid[2,:,:,0,4]=-x*root
        rigid[0,:,:,0,5]=-y*root;rigid[1,:,:,0,5]=x*root
        rigid=rigid.reshape(3*n,6)
        chol=np.linalg.cholesky(rigid.T@(self.mass@rigid))
        self.rigid=np.linalg.solve(chol,rigid.T).T;self.mr=self.mass@self.rigid
        fixed=[0,n,2*n,n+(nxy-1)*nxy*(PZ+1),2*n+(nxy-1)*nxy*(PZ+1),2*n+(nxy-1)*(PZ+1)]
        self.free=np.setdiff1d(np.arange(3*n),fixed)
        self.factor=splu(self.k[self.free][:,self.free].tocsc())

    def solve(self,force):
        balanced=force-self.mr@(self.rigid.T@force)
        u=np.zeros_like(force);u[self.free]=self.factor.solve(balanced[self.free])
        u-=self.rigid@(self.rigid.T@(self.mass@u))
        residual=np.linalg.norm(self.k@u-balanced)/np.linalg.norm(balanced)
        return u,float(residual)

def speed_proxy(c,omega,massz,kz,bulk_fraction):
    """Accepted modal speed convention: kz is the lateral gradient of u_z."""
    bulk=(c.C11+2*c.C12)/3
    mu=((c.C11-c.C12+3*c.C44)/5+5*(c.C11-c.C12)*c.C44/(4*c.C44+3*(c.C11-c.C12)))/2
    young=9*bulk*mu/(3*bulk+mu);nu=(3*bulk-2*mu)/(2*(3*bulk+mu))
    speed=np.where(bulk_fraction>.5,np.sqrt((bulk+4*mu/3)/c.RHO),np.sqrt(mu/c.RHO))
    flex=(massz>.5)&(kz*c.DIMS[2]<.5)
    speed[flex]=2*np.sqrt(young*c.DIMS[2]**2/(12*(1-nu*nu)*c.RHO))*kz[flex]
    return speed

def uniform_reference(c,wavelength,theta,pz=4):
    """Full-vector free-thickness reference, independent of the lateral mesh."""
    from scipy.linalg import eigh
    qz,wz=leggauss(pz+2);h=c.DIMS[2];phi=[];der=[]
    for j in range(pz+1):
        coef=np.zeros(j+1);coef[j]=1;scale=np.sqrt((2*j+1)/h)
        phi.append(legval(qz,coef)*scale);der.append(legval(qz,legder(coef))*scale*2/h)
    phi=np.array(phi).T;der=np.array(der).T;nz=pz+1;kx,ky=2*np.pi/wavelength*np.array([np.cos(theta),np.sin(theta)])
    b=np.zeros((len(qz),6,3*nz),complex)
    b[:,0,:nz]=1j*kx*phi;b[:,1,nz:2*nz]=1j*ky*phi;b[:,2,2*nz:]=der
    b[:,3,nz:2*nz]=der;b[:,3,2*nz:]=1j*ky*phi;b[:,4,:nz]=der;b[:,4,2*nz:]=1j*kx*phi;b[:,5,:nz]=1j*ky*phi;b[:,5,nz:2*nz]=1j*kx*phi
    exact=np.einsum('qim,ij,qjn,q->mn',b.conj(),c.cubic(),b,wz*h/2,optimize=True)
    values,reference=eigh((exact+exact.conj().T)/2,c.RHO*np.eye(3*nz),subset_by_index=(0,2),check_finite=False)
    return np.sqrt(values)/(2*np.pi),reference

def plane_branches(c,wavelength,theta,degree=None,cells=None,details=False):
    """Uniform-reference Bloch dispersion for choosing the lateral hp space."""
    from scipy.linalg import eigh
    cells=CELLS if cells is None else cells
    degree=PX if degree is None else degree
    smap=np.array([[0,5,4],[5,1,3],[4,3,2]]);q,w=leggauss(degree+2)
    nodes=np.linspace(-1,1,degree+1);coeff=[]
    for i,node in enumerate(nodes):
        roots=np.delete(nodes,i);coeff.append(np.polynomial.polynomial.polyfromroots(roots)/np.prod(node-roots))
    coeff=np.array(coeff);val=np.polynomial.polynomial.polyval(q,coeff.T).T
    derivative=np.polynomial.polynomial.polyder(coeff,axis=1)
    forms=[]
    for length,k in zip(c.DIMS[:2]/cells,2*np.pi/wavelength*np.array([np.cos(theta),np.sin(theta)])):
        der=np.polynomial.polynomial.polyval(q,derivative.T).T*2/length
        transform=np.zeros((degree+1,degree),complex);transform[:degree]=np.eye(degree);transform[-1,0]=np.exp(1j*k*length)
        forms.append([[transform.conj().T@(a.T@((w*length/2)[:,None]*b))@transform for b in (val,der)] for a in (val,der)])
    x,y=forms;z=[[a.toarray() for a in row] for row in zforms(c.DIMS[2])]
    n=degree*degree*(PZ+1);tensor=c.cubic();blocks=[]
    for i in range(3):
        row=[]
        for j in range(3):
            block=np.zeros((n,n),complex)
            for a in range(3):
                for b in range(3):
                    block+=tensor[smap[i,a],smap[j,b]]*np.kron(np.kron(x[int(a==0)][int(b==0)],y[int(a==1)][int(b==1)]),z[int(a==2)][int(b==2)])
            row.append(block)
        blocks.append(row)
    stiffness=np.block(blocks);mass=c.RHO*np.kron(np.eye(3),np.kron(np.kron(x[0][0],y[0][0]),z[0][0]))
    values,vectors=eigh((stiffness+stiffness.conj().T)/2,mass,check_finite=False)
    # Match physical k, not the lowest folded Bloch eigenvalues. Higher-p cells
    # can represent several wavelengths per cell, so frequency order aliases k.
    _,reference=uniform_reference(c,wavelength,theta,PZ);nz=PZ+1
    kx,ky=2*np.pi/wavelength*np.array([np.cos(theta),np.sin(theta)])
    xn=np.arange(degree)*c.DIMS[0]/(cells*degree);yn=np.arange(degree)*c.DIMS[1]/(cells*degree)
    phase=np.exp(1j*(kx*xn[:,None]+ky*yn[None,:]))
    target=np.einsum('ij,czb->cijzb',phase,reference.reshape(3,nz,3)).reshape(3*n,3)
    overlap=abs(vectors.conj().T@(mass@target))**2
    norm=np.real(np.sum(target.conj()*(mass@target),axis=0))
    match=np.argmax(overlap,axis=0);scores=overlap[match,np.arange(3)]/norm
    frequency=np.sqrt(values[match])/(2*np.pi)
    return (frequency,scores) if details else frequency
