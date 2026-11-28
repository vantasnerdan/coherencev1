"""Retarded D-covariance loading of the *same* acoustic and optical fields.

The equilibrium real stiffness is a counterterm, not a second elastic baseline.
The source is the G3 Hamiltonian derivative; no attenuation coefficient is fitted.
All frequencies below are cyclic Hz at the interface and angular fs^-1 internally.
"""
import numpy as np
import native
import sources as s


def covariance(arm,frequency,temperature=300.,linear_strain=0.,electronic_u=0.,frozen=False):
    p=s.inputs()['lattice'];k=native.primitive_basis(p)
    _,_,D,_=native.dispersion(p,k,arm['allocation'],electronic_u,frozen=frozen)
    g=p['collective_strain_slope_eV_A2'];C=2*p['K4_100_eV_A4']+3*p['K4_111_eV_A4']
    D=D+g*(linear_strain+s.ALPHA*(temperature-300.))
    _,Q,_,_=native.thermal_pole(p,D,temperature)
    theta=s.KB*temperature;mass=p['mass_eVc2']*s.QE/299792458.**2*1e10/s.QE
    V=theta/mass;dQ=-Q/(D+2*C*Q);kappa=1/(arm['opt']*1e15)
    A=np.zeros((len(k),3,3));A[:,0,0]=-kappa;A[:,0,1]=2
    A[:,1,0]=-(D+2*C*Q)/mass;A[:,1,1]=-kappa;A[:,1,2]=1
    A[:,2,1]=-2*(D+C*Q)/mass;A[:,2,2]=-kappa
    b=np.stack((kappa*dQ,-Q/mass,np.zeros_like(Q)),axis=-1)
    w=2*np.pi*np.asarray(frequency)*1e-15
    operator=-1j*w[...,None,None,None]*np.eye(3)-A
    chi=np.linalg.solve(operator,b[...,None])[...,0]
    derivative=np.linalg.solve(operator,(1j*chi)[...,None])[...,0]*1e-15
    factor=s.QE*g*g/6
    bulk=factor*(chi[...,0].mean(axis=-1)-dQ.mean())
    d_bulk=factor*derivative[...,0].mean(axis=-1)
    # Exact quadratic relative Gaussian free energy above instantaneous
    # equilibrium. This is stored bound response, not another positive source.
    x=chi[...,0]-dQ
    availability=(.125*(theta/(Q*Q)+C)*abs(x)**2+
        .125*theta/(V*V)*abs(chi[...,2])**2+.25*theta/(Q*V)*abs(chi[...,1])**2).mean(axis=-1)
    availability*=s.QE*g*g/3
    return dict(bulk_per_D_J=bulk,bulk_derivative_per_D_J_s=d_bulk,
        availability_per_D_J=availability,
        Q_per_linear_strain_A2=chi[...,0].mean(axis=-1)*g/3,
        Q_nonequilibrium_per_linear_strain_A2=(chi[...,0].mean(axis=-1)-dQ.mean())*g/3,
        theta_per_linear_strain_eV=mass*chi[...,2].mean(axis=-1)*g/3,
        frozen_reference=frozen)


def coefficients(arm,frequency,number_density,phase_velocity,group_velocity,hydro_ratio,
                 temperature=300.,linear_strain=0.,frozen=False,hydro_logarithmic_k_slope=None):
    """First-order retained-mode projection of the passive susceptibility.

    hydro_ratio is <div u>^2 / total bare mode energy density, in Pa^-1.
    Group delay and energy transport are distinct in an absorptive medium.
    The latter follows the projected energy flux and remains a finite front.
    """
    r=covariance(arm,frequency,temperature,linear_strain,frozen=frozen)
    B=number_density*r['bulk_per_D_J'];Bd=number_density*r['bulk_derivative_per_D_J_s']
    psi=number_density*r['availability_per_D_J'];w=2*np.pi*frequency
    fraction=float(B.real*hydro_ratio)
    if 1+fraction<=0:raise ValueError('Unstable coupled acoustic projection')
    evolved=w*np.sqrt(1+fraction)
    energy_factor=1+2*hydro_ratio*(psi-.25*B.real)/(1+fraction)
    if energy_factor<=0:raise ValueError('Nonpositive coupled acoustic energy metric')
    # Fixed retained host k and shape. A q-resolved native PdD dispersion was
    # not supplied; derivatives of its weak shape correction are not invented.
    slope=2*(1-group_velocity/phase_velocity) if hydro_logarithmic_k_slope is None else hydro_logarithmic_k_slope
    numerator=w*group_velocity*(1+fraction)+.5*w*phase_velocity*fraction*slope
    energy_velocity=numerator/(evolved*energy_factor)
    group=numerator/(evolved-.5*w*w*hydro_ratio*Bd.real)
    optical_rate=-evolved*B.imag*hydro_ratio/((1+fraction)*energy_factor)
    if optical_rate < -1e-9 or energy_velocity<=0:raise ValueError('Nonpassive reciprocal optical response')
    return dict(omega=evolved,phase_velocity=phase_velocity*np.sqrt(1+fraction),
        group_velocity=float(group),energy_velocity=float(energy_velocity),
        optical_energy_loss_rate=float(optical_rate),energy_factor=float(energy_factor),
        D_availability_fraction=float(2*hydro_ratio*psi/((1+fraction)*energy_factor)),
        hydro_per_energy=float(hydro_ratio/((1+fraction)*energy_factor)),
        Q_response=r['Q_per_linear_strain_A2'],theta_response=r['theta_per_linear_strain_eV'])


def relative_free_energy(Q,Cov,V,Qeq,Veq,theta,quartic):
    """Positive, cancellation-safe Gaussian availability, eV per coordinate."""
    x=(Q-Qeq)/Qeq;y=(V-Veq)/Veq;z2=Cov*Cov/(Q*V)
    if np.any(x<=-1) or np.any(y<=-1) or np.any(z2>=1):raise ValueError('Nonpositive covariance')
    def remainder(a):
        small=abs(a)<1e-4
        series=a*a*(.5+a*(-1/3+a*(.25+a*(-.2+a/6))))
        return np.where(small,series,a-np.log1p(a))
    return .5*theta*(remainder(x)+remainder(y)-np.log1p(-z2))+.25*quartic*(Q-Qeq)**2
