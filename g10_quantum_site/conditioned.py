"""Rao-Blackwell bead-coordinate marginal, not a Gaussian physical tail.

Condition one selected bead on its two path neighbours and all other nuclei.
The Gaussian spring/quadratic factor is integrated against the actual quartic;
normalization uses Gauss-Hermite, with a second-order comparison recorded.
"""
import numpy as np
from scipy.special import roots_hermite,logsumexp
from model import KB,kinetic_coefficient


def quartic(x,a,b):
    xx=x*x
    return (a-b)/4*np.sum(xx*xx,axis=-1)+b/4*np.sum(xx,axis=-1)**2


def conditional_data(x,field,model,mass,T,order):
    F,P,_=x.shape;site=model['site'];beta=1/(KB*T);t=kinetic_coefficient(mass)
    beads=np.unique(np.linspace(0,P-1,min(P,8),dtype=int))
    precision=beta/P*model['K'][site]+P/(beta*t)*np.eye(3)
    covariance=np.linalg.inv(precision);chol=np.linalg.cholesky(covariance)
    linear=beta/P*(model['f'][site]+field[:,beads])+P/(2*beta*t)*(np.roll(x,1,axis=1)[:,beads]+np.roll(x,-1,axis=1)[:,beads])
    centers=np.einsum('ij,fpj->fpi',covariance,linear)
    q,w=roots_hermite(order)
    xyz=np.array(np.meshgrid(q,q,q,indexing='ij')).reshape(3,-1).T
    offsets=np.sqrt(2)*xyz@chol.T
    weights=(w[:,None,None]*w[None,:,None]*w[None,None,:]).ravel()/np.pi**1.5
    lognorm=np.empty(centers.shape[:-1])
    flat=centers.reshape(-1,3)
    for start in range(0,len(flat),128):
        vals=quartic(flat[start:start+128,None,:]+offsets[None,:,:],model['a'],model['b'])
        lognorm.reshape(-1)[start:start+128]=logsumexp(np.log(weights)-beta/P*vals,axis=1)
    return centers,precision,lognorm


def radial(x,field,model,mass,T,logweights=None):
    F,P,_=x.shape;beta=1/(KB*T);dt=beta/P
    centers,precision,lognorm=conditional_data(x,field,model,mass,T,10)
    _,_,comparison=conditional_data(x,field,model,mass,T,12)
    lognorm_delta=float(np.max(abs(lognorm-comparison)))
    radii=np.linspace(0,1,181);mu,muw=np.polynomial.legendre.leggauss(10);phi=np.arange(20)*2*np.pi/20
    dirs=np.stack(np.broadcast_arrays(np.sqrt(1-mu[:,None]**2)*np.cos(phi),np.sqrt(1-mu[:,None]**2)*np.sin(phi),mu[:,None]+phi*0),axis=-1).reshape(-1,3)
    angular=np.repeat(muw,20)/40
    points=(radii[:,None,None]*dirs[None,:,:]).reshape(-1,3)
    exponent=-.5*np.einsum('pa,ab,pb->p',points,precision,points)-dt*quartic(points,model['a'],model['b'])
    prefactor=np.linalg.det(precision)**.5/(2*np.pi)**1.5
    weights=np.ones(F) if logweights is None else np.exp(logweights-np.max(logweights))
    nb=12;bs=F//nb;curves=[];den=[]
    for block in range(nb):
        sl=slice(block*bs,(block+1)*bs);center=centers[sl].reshape(-1,3)
        cc=-.5*np.einsum('fa,ab,fb->f',center,precision,center)-lognorm[sl].ravel()
        ww=np.repeat(weights[sl],centers.shape[1]);curve=np.empty(len(points))
        for start in range(0,len(points),256):
            p=points[start:start+256]
            values=np.exp(exponent[start:start+256,None]+p@precision@center.T+cc)
            curve[start:start+256]=prefactor*(values@ww)/(bs*centers.shape[1])
        curves.append(curve.reshape(len(radii),-1)@angular);den.append(weights[sl].mean())
    curves=np.array(curves);den=np.array(den);value=curves.sum(axis=0)/den.sum()
    blocks=value+(curves-value*den[:,None])/den.mean()
    se=blocks.std(axis=0,ddof=1)/np.sqrt(nb)
    measure=4*np.pi*radii*radii*value;measure_se=4*np.pi*radii*radii*se
    survival=np.zeros(len(radii));survival[:-1]=np.cumsum(((measure[:-1]+measure[1:])/2*np.diff(radii))[::-1])[::-1]
    return dict(radius=radii,density=value,density_se=se,measure=measure,measure_se=measure_se,tail=survival,
                resolved=value>3*se,normalization=float(np.trapezoid(measure,radii)),
                conditional_log_normalizer_order10_vs12=lognorm_delta,
                support_fraction=float(1-np.interp(model['support_A'],radii,survival)),
                frame_weight_ess=float(weights.sum()**2/(weights@weights)),
                scope='finite-P conditional quartic coordinate marginal; quadrature-normalized, block SE; rare-bin precision does not repair unsampled path-centre overlap')
