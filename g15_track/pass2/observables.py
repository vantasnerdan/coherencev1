"""Consumer sufficient measures integrated on the actual unequal source steps.
Gaussian phase-symmetric ordinary-mode population is a declared approximation;
its energy, unsupported centroid tails and directional pair flux stay separate.
"""
import numpy as np
from scipy.special import ndtr,log_ndtr,gammaincc


def normal_interval(a,b):
    out=np.empty_like(a)
    positive=a>=0;negative=b<=0;middle=~(positive|negative)
    x=log_ndtr(-a[positive]);y=log_ndtr(-b[positive])
    out[positive]=np.exp(x)*(-np.expm1(y-x))
    x=log_ndtr(b[negative]);y=log_ndtr(a[negative])
    out[negative]=np.exp(x)*(-np.expm1(y-x))
    out[middle]=ndtr(b[middle])-ndtr(a[middle])
    return out


def positive_normal(mean,sigma):
    out=np.maximum(mean,0.)
    good=sigma>0
    a=mean[good]/sigma[good]
    out[good]=sigma[good]*np.exp(-.5*a*a)/np.sqrt(2*np.pi)+mean[good]*ndtr(a)
    return out


class ConsumerMeasures:
    def __init__(self,m,ec,sp,N,weights,pair_weights,shells,env):
        self.m=m;self.N=N;self.w=weights;self.pw=pair_weights;self.shells=shells
        self.population=np.linspace(0,1,33)**2*(2000*ec/N)
        cap=.36 if env=='bulk' else .24
        vmax=np.sqrt(m.K/m.M)*cap*1.000001
        self.labs=np.array([2*ec,ec]);self.vb=np.sqrt(2*self.labs*1000/m.M)
        margin=2*vmax/self.vb+(vmax/self.vb)**2
        self.scale=self.labs*margin
        self.moments=np.zeros((len(m.BRANCHES)*3,2,33,9))
        self.pair=np.zeros((len(m.BRANCHES)*3,len(shells),7))
        self.population_tail=np.zeros((len(m.BRANCHES)*3,2))
        self.Ecap=.5*m.K*cap*cap;self.width=.05

    def add(self,dt,velocity,pool,rvec,vrel,pair_active):
        m=self.m;L=len(velocity);ep=np.maximum(pool,0.)/self.N
        pi=np.clip(np.searchsorted(self.population,ep,side='right')-1,0,31)
        fraction=(ep-self.population[pi])/(self.population[pi+1]-self.population[pi])
        if fraction.max()>1+1e-8:raise ValueError('Population exceeds the incoming parent budget')
        lane=np.arange(L)
        velocity2=np.sum(velocity*velocity,axis=2)
        for incoming,vb in enumerate(self.vb):
            # Expand about the exact quiet energy: no subtraction of vb².
            x=m.M*(.5*velocity2-vb*velocity[:,:,2])/(1000*self.scale[incoming])
            if abs(x).max()>1+1e-6:raise ValueError('A bound velocity escaped the consumer interval')
            power=np.ones_like(x)
            for order in range(9):
                value=dt*np.sum(self.w[None,:]*power,axis=1)
                for pp,wp in ((pi,1-fraction),(pi+1,fraction)):
                    np.add.at(self.moments[:,incoming,:,order],(lane,pp),value*wp)
                power*=x
        # Relative phase-space covariance of two independent ordinary D modes.
        sigmaq2=2*ep/(3*m.K);sigmav2=2*ep/(3*m.M)
        tail=np.zeros(L);energy_tail=np.zeros(L);nonzero=ep>0
        tail[nonzero]=gammaincc(3,3*self.Ecap/ep[nonzero])
        energy_tail[nonzero]=gammaincc(4,3*self.Ecap/ep[nonzero])
        self.population_tail[:,0]+=dt*self.N*tail
        self.population_tail[:,1]+=dt*pool*energy_tail
        r=np.linalg.norm(rvec,axis=2);speed2=np.sum(vrel*vrel,axis=2)
        parallel=np.sum(rvec*vrel,axis=2)/r
        for k,shell in enumerate(self.shells):
            deterministic=(abs(r-shell)<=self.width/2)
            pband=deterministic.astype(float)
            mean=parallel.copy();angular_variance=np.zeros_like(r)
            good=np.broadcast_to(nonzero[:,None],r.shape)
            if np.any(good):
                sig=np.sqrt(np.broadcast_to(sigmaq2[:,None],r.shape)[good]);rr=r[good]
                lo=shell-self.width/2;hi=shell+self.width/2
                a=(lo-rr)/sig;b=(hi-rr)/sig;c=(lo+rr)/sig;d=(hi+rr)/sig
                probability=normal_interval(a,b)+normal_interval(c,d)
                probability+=sig/(rr*np.sqrt(2*np.pi))*(np.exp(-d*d/2)-np.exp(-b*b/2)-np.exp(-c*c/2)+np.exp(-a*a/2))
                if probability.min() < -1e-10 or probability.max()>1+1e-10:raise ValueError('Noncentral radial probability failed')
                pband[good]=np.clip(probability,0,1)
                concentration=shell*rr/sig**2
                langevin=1/np.tanh(concentration)-1/concentration
                small=concentration<.001
                langevin[small]=concentration[small]/3-concentration[small]**3/45
                T=langevin/concentration
                mean[good]=parallel[good]*langevin
                angular_variance[good]=np.maximum(0,speed2[good]*T+parallel[good]**2*(1-3*T-langevin**2))
            sigma=np.sqrt(np.broadcast_to(sigmav2[:,None],r.shape))
            lower=positive_normal(mean,sigma)
            upper=np.minimum(.5*(mean+np.sqrt(mean*mean+angular_variance+sigma*sigma)),lower+np.sqrt(angular_variance))
            central=positive_normal(mean,np.sqrt(sigma*sigma+angular_variance))
            weights=self.pw[None,:]*pair_active*dt/self.width
            for column,flux in enumerate((central,lower,upper)):
                self.pair[:,k,column]+=np.sum(weights*pband*flux,axis=1)
            self.pair[:,k,3]+=np.sum(weights*deterministic*np.maximum(parallel,0.),axis=1)
            # A union/Cauchy bound, not a claimed rate of an unsampled rare event.
            unresolved=np.sqrt((speed2+sigma*sigma)*np.minimum(1.,2*tail[:,None]))
            self.pair[:,k,4]+=np.sum(weights*unresolved,axis=1)
            self.pair[:,k,5]+=np.sum(weights*self.width*pband*(.25*m.M*(speed2+3*sigma*sigma)),axis=1)
            self.pair[:,k,6]+=np.sum(weights*self.width*pband,axis=1)

    def arrays(self):
        return dict(beam_joint_energy_moments_D_fs=self.moments,beam_relative_energy_scale_keV=self.scale,
            incoming_energy_keV=self.labs,
            population_energy_grid_eV_D=self.population,pair_band_measures=self.pair,
            pair_band_columns=np.array(['Gaussian_angular_moment_flux','angular_Jensen_lower','angular_Cauchy_upper',
                'coherent_same_band_flux','unsupported_mode_tail_flux_upper','relative_KE_band_eV_fs','band_residence_pair_fs']),
            pair_band_width_A=self.width,population_unsupported_D_fs=self.population_tail[:,0],
            population_unsupported_energy_eV_fs=self.population_tail[:,1])

def inward_crossing(previous,current,shells):
    """Entry roots of |r0+s*(r1-r0)|=R, open-left / closed-right intervals.
    A secant relative velocity is the derivative of this same interpolant.
    A tangency has zero inward flux and is not an entry; exit roots are omitted.
    """
    dr=current-previous
    a=np.sum(dr*dr,axis=-1)[...,None]
    b=2*np.sum(previous*dr,axis=-1)[...,None]
    c=np.sum(previous*previous,axis=-1)[...,None]-np.asarray(shells)**2
    discriminant=b*b-4*a*c
    tolerance=32*np.finfo(float).eps*(b*b+4*a*abs(c))
    denominator=-b+np.sqrt(np.maximum(discriminant,0.))
    root=np.full_like(c,2.)
    np.divide(2*c,denominator,out=root,where=denominator>0)
    entry=(discriminant>tolerance)&(b<0)&(root>1e-10)&(root<=1+1e-10)
    return entry,np.clip(root,0,1)
