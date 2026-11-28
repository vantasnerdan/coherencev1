"""Full-Cox-history initialization of retained warm thermoelastic modes.

This is a resolvent identity for the same FV heat and native modal operators,
not an independent foil or source. It avoids erasing millisecond Cu memory when
conditioning a short phase window. Fine evolution subsequently uses the actual
conservative delayed source ports and common Poisson count innovations.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
import fields as f
import sources as s


def filtered_clock(calendar,time,poles,name):
    poles=np.asarray(poles,complex);decay=np.exp(poles*time)
    constant=np.expm1(poles*time)/poles
    if name!='hot_bursts':return constant
    omega=calendar.omega;phase=calendar.phase
    ripple=.35*(np.exp(1j*phase)*(np.exp(1j*omega*time)-decay)/(1j*omega-poles)+
                np.exp(-1j*phase)*(np.exp(-1j*omega*time)-decay)/(-1j*omega-poles))
    a=calendar.edges[:-1];b=np.minimum(calendar.edges[1:],time);mask=calendar.state&(a<time)
    a=a[mask];b=b[mask];burst=np.zeros(len(poles),complex)
    for j in range(0,len(poles),12):
        alpha=poles[j:j+12,None]
        burst[j:j+12]=np.sum(np.exp(alpha*(time-b))*np.expm1(alpha*(b-a))/alpha,axis=1)
    return (constant+ripple+9*burst)/calendar.norm


def state(calendar,time,beam,modes,temperature,heat_power,direct_force,gamma):
    """q,p from the same complete conditional history, at the actual warm T.

    For T'=-A T+h C and alpha=-gamma/2+i*wd,
      integral exp(alpha(t-s))T(s) ds = (A+alpha)^-1(h Z-T(t)).
    The boundary heat flux is retained in A. Short-carrier timing is resolved by
    the subsequent window; only its prehistory uses the accumulated heat source.
    """
    omega=modes['omega'];wd=np.sqrt(omega*omega-gamma*gamma/4);poles=-gamma/2+1j*wd
    clock=filtered_clock(calendar,time,poles,beam);heat=f.Heat();cap=f.CAP.ravel()
    A=sparse.diags(1/cap)@heat.G;h=np.asarray(heat_power).ravel()/cap;T=np.asarray(temperature).ravel()
    S=modes['strain'][:,:3].sum(axis=1);beta=3*s.BULK*s.ALPHA*f.VOL.ravel()
    convolution=np.zeros(len(omega),complex)
    for j,pole in enumerate(poles):
        response=splu((A+sparse.eye(len(cap))*pole).tocsc()).solve(h*clock[j]-T)
        convolution[j]=direct_force[j]*clock[j]+np.dot(beta*S[:,j],response)
    q=convolution.imag/wd;velocity=convolution.real-gamma*q/2
    return q,velocity
