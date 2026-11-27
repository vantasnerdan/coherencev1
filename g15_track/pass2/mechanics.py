"""Weak, bound screened-force histories and separate energetic scattering marks."""
import numpy as np
from scipy.integrate import cumulative_trapezoid


class ForceHistory:
    """Cumulative impulse on the accepted piecewise-linear projectile clock.
    This is a weak straight-ray force model, not a second recoil energy source.
    Equal-mass ZBL scattering bounds its amplitude; energetic marks stay outside
    the fitted bound chart and retain their own original moving-D channel.
    """
    def __init__(self,m,ec,sp,positions,env,direct_scale):
        self.m=m;self.pos=np.asarray(positions);self.n=len(positions)
        tr=m.p1.track(ec,sp)
        self.z=np.r_[0.,[r['depth_A'] for r in tr['timeline']]]
        self.t=np.r_[0.,[r['time_fs'] for r in tr['timeline']]]
        self.v0=np.sqrt(2*tr['initial_eV']/m.M);self.vend=np.sqrt(2*tr['terminal_unresolved_eV']/m.M)
        en=np.interp(self.pos[:,2],self.z,np.r_[tr['initial_eV'],[r['energy_eV'] for r in tr['timeline']]])
        self.b=np.maximum(np.linalg.norm(self.pos[:,:2],axis=1),.003)
        unit=np.zeros((self.n,3));unit[:,:2]=self.pos[:,:2]/self.b[:,None]
        self.radial=unit
        self.eta=np.linspace(-np.arcsinh(24.),np.arcsinh(24.),1025)
        s=self.b[:,None]*np.sinh(self.eta)[None,:]
        projectile=self.pos[:,2,None]+s
        tm=np.interp(projectile.ravel(),self.z,self.t).reshape(s.shape)
        tm=np.where(projectile<0,projectile/self.v0,tm)
        tm=np.where(projectile>self.z[-1],self.t[-1]+(projectile-self.z[-1])/self.vend,tm)
        r=np.sqrt(self.b[:,None]**2+s*s);a=.8854*.529177/2
        force=14.39964*(m.p1.bca.zbl_phi(r/a)/r**2-m.p1.bca.zbl_dphi(r/a)/(a*r))
        force*=((projectile>=0)&(projectile<=self.z[-1]))
        fr=force*self.b[:,None]/r;fz=-force*s/r
        ir=cumulative_trapezoid(fr,tm,axis=1,initial=0.)
        iz=cumulative_trapezoid(fz,tm,axis=1,initial=0.)
        _,T,raw,theta=m.force_signature(self.pos,en,m.OMEGA)
        cap=.5*m.K*(.36 if env=='bulk' else .24)**2+.25*m.A*(.36 if env=='bulk' else .24)**4
        norm=np.minimum(1.,np.sqrt(2*m.M*T)/np.maximum(np.hypot(ir[:,-1],iz[:,-1]),1e-100))*np.sqrt(direct_scale)
        norm[T>cap]=0.
        self.ir=ir*norm[:,None];self.iz=iz*norm[:,None]
        self.T=T;self.raw=raw;self.theta=theta;self.strong=T>cap
        self.clock=np.interp(self.pos[:,2],self.z,self.t)
        self.transit_fs=self.b/np.sqrt(2*en/m.M)
        self.tm=tm;self.fr=fr*norm[:,None];self.fz=fz*norm[:,None]

    def indices(self,at):
        at=np.asarray(at,float)
        zp=np.interp(at.ravel(),self.t,self.z).reshape(at.shape)
        zp=np.where(at<0,at*self.v0,zp)
        zp=np.where(at>self.t[-1],self.z[-1]+(at-self.t[-1])*self.vend,zp)
        eta=np.arcsinh((zp-self.pos[:,2])/self.b)
        index=np.clip((eta-self.eta[0])/(self.eta[1]-self.eta[0]),0,len(self.eta)-1)
        lo=np.minimum(index.astype(int),len(self.eta)-2)
        return lo,index-lo

    def lookup(self,a,indices):
        lo,f=indices;site=np.arange(self.n)
        return a[site,lo]*(1-f)+a[site,lo+1]*f

    def cumulative(self,at):
        indices=self.indices(at)
        return self.lookup(self.ir,indices)[...,None]*self.radial+self.lookup(self.iz,indices)[...,None]*np.array([0.,0.,1.])

    def interval(self,ta,tb,shift):
        return self.cumulative(tb-shift)-self.cumulative(ta-shift)

    def motion(self,at,omega,kappa):
        """Finite-time damped Green, retaining the field before closest approach."""
        alpha=kappa/2;wd=np.sqrt(omega*omega-alpha*alpha)
        age=np.asarray(at)-self.clock;clock=self.tm-self.clock[:,None]
        cosine=np.exp(alpha*clock)*np.cos(wd*clock)
        sine=np.exp(alpha*clock)*np.sin(wd*clock)
        indices=self.indices(at)
        cc=[];ss=[]
        for force in (self.fr,self.fz):
            cc.append(self.lookup(cumulative_trapezoid(force*cosine,self.tm,axis=1,initial=0.),indices))
            ss.append(self.lookup(cumulative_trapezoid(force*sine,self.tm,axis=1,initial=0.),indices))
        C=cc[0][...,None]*self.radial+cc[1][...,None]*np.array([0.,0.,1.])
        S=ss[0][...,None]*self.radial+ss[1][...,None]*np.array([0.,0.,1.])
        co=np.cos(wd*age)[...,None];si=np.sin(wd*age)[...,None];decay=np.exp(-alpha*age)[...,None]
        q=decay*(si*C-co*S)/(self.m.M*wd)
        v=decay*((co-alpha/wd*si)*C+(si+alpha/wd*co)*S)/self.m.M
        return q,v
