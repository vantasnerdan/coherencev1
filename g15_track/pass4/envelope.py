"""Physical Duffing packet, reciprocal parameter work and existing carrier ports.

The coherent packet replaces a finite part of the same optical budget. Signed
corrections, not new power, propagate through the existing positive generator.
"""
import numpy as np
from scipy.linalg import expm
import sources as s


def step(q,p,K,A,M,dt):
    """Discrete-gradient Hamiltonian step: kinetic plus K q²/2+A q⁴/4."""
    omega=np.sqrt(K/M)
    if omega*abs(dt)>=np.pi:raise ValueError('Optical carrier under-resolved')
    # Phase-fitted discrete gradient: exact harmonic rotation and the same
    # exact nonlinear energy identity, without a secular Cayley phase error.
    dt=2*np.tan(omega*dt/2)/omega
    old=q;new=q+dt*p/M
    for _ in range(12):
        force=.5*K*(new+old)+.25*A*(new+old)*(new*new+old*old)
        derivative=.5*K+.25*A*(3*new*new+2*new*old+old*old)
        correction=(new-old-dt*p/M+dt*dt*force/(2*M))/(1+dt*dt*derivative/(2*M))
        new-=correction
        if abs(correction)<2e-16:break
    else:raise ValueError('Coherent Hamiltonian step did not converge')
    return new,2*M*(new-old)/dt-p


def trace(arm,cache,birth,ages,Db,sensitivity,Qoffset,slow,q0,v0,sites,mass,quartic,covariance_coefficient,network,frozen=False):
    dt=2e-15;t=np.arange(int(np.ceil(ages[-1]/dt))+1)*dt;clock=birth+t
    M=mass*1e10/s.QE;q=q0*1e10;p=M*v0*1e-5
    selected=int(np.argmin(abs(cache['base_K'][cache['selected_mode_indices']]-Db)))
    Q=np.interp(clock,cache['time_s'],cache['variance_q_selected_A2'][:,selected])+Qoffset
    C=covariance_coefficient
    electronic=sensitivity*np.interp(clock,cache['time_s'],cache['projected_deltaK_eVA2']) if not frozen else np.zeros(len(t));population=C*Q
    hot_reference=np.interp(clock,cache['time_s'],cache['hot_e_eV'])
    thermal_reference=np.interp(clock,cache['time_s'],cache['thermal_e_excess_eV'])
    optical_reference=np.interp(clock,cache['time_s'],cache['optical_excess_eV'])
    K=Db+slow+electronic+population
    if np.min(K)<=0:raise ValueError('Nonpositive propagated coherent curvature')
    def energy(q,p,k):return sites*(p*p/(2*M)+.5*k*q*q+.25*quartic*q**4)
    initial=energy(q,p,K[0]);reference=initial;correction=np.zeros(9);heat=0.;we=0.;wo=0.;loss=0.;balance=0.;minimum=np.inf
    propagation=expm(network*dt);attenuation=np.exp(-dt/(2*arm['opt']));rows=[];phase=0.;old_angle=np.arctan2(p/np.sqrt(M*K[0]),q)
    for j in range(len(t)):
        E=energy(q,p,K[j]);w=np.exp(-2*t[j]/arm['dephase']);omega=np.sqrt(K[j]/M)*(1+3*quartic*E/(4*sites*K[j]**2))*1e15
        angle=np.arctan2(p/np.sqrt(M*K[j]),q);increment=(old_angle-angle+np.pi)%(2*np.pi)-np.pi
        if j:phase+=increment
        old_angle=angle
        defect=E-reference+correction.sum()+heat;balance=max(balance,abs(defect))
        hot=float(hot_reference[j])+correction[3]
        thermal=float(thermal_reference[j])+correction[8]
        optical=float(optical_reference[j])-reference+correction[4]
        minimum=min(minimum,hot,thermal,optical)
        if minimum < -1e-10:raise ValueError('Coherent work exhausts a physical source carrier')
        rows.append((q*1e-10,p/M*1e5,E*w,E*(1-w),we,wo,heat,defect,phase,omega,E,loss,minimum,optical))
        if j==len(t)-1:break
        before=correction.sum();correction=propagation@correction;heat+=before-correction.sum()
        for half in (0,1):
            if half:
                q,p=step(q,p,(K[j]+K[j+1])/2,quartic,M,dt*1e15)
            ee=.25*sites*q*q*(electronic[j+1]-electronic[j]);oo=.25*sites*q*q*(population[j+1]-population[j])
            denominator=hot+thermal
            fraction=hot/denominator if denominator>0 else .5
            correction[3]-=ee*fraction;correction[8]-=ee*(1-fraction);correction[4]-=oo;we+=ee;wo+=oo
            # Linear strain/temperature are held at the actual pre-impact state
            # during this picosecond conditional window, so no hidden T work.
        before=energy(q,p,K[j+1]);q*=attenuation;p*=attenuation
        debit=before-energy(q,p,K[j+1]);base_debit=reference*(-np.expm1(-dt/arm['opt']));reference-=base_debit
        correction[6]+=.8*(debit-base_debit);correction[8]+=.2*(debit-base_debit);loss+=debit
    array=np.array(rows)
    sampled=np.column_stack([np.interp(ages,t,array[:,j]) for j in range(array.shape[1])])
    return sampled,dict(initial_coherent_eV=initial,maximum_reciprocal_balance_eV=balance,minimum_source_carrier_eV=minimum,
        electron_work_eV=we,incoherent_optical_work_eV=wo,physical_damping_eV=loss,
        relative_roundoff=balance/max(initial,1e-300),source_birth_age_s=birth,internal_step_s=dt)
