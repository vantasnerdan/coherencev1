"""LSS Eq(4.1) projected-range transport with ZBL binary scattering.

Lindhard/Scharff/Schiott1963 pp24-25: infinite homogeneous-medium first
projected moment. This is NOT a half-space implanted-only conditional mean.
Continuous LS electronic stopping; discrete nuclear energy/angle transfers.
Units: eV, nm, amu, atoms/nm^3. No Monte Carlo or range calibration.
"""
import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

E2_EV_NM = 1.43996454784255
A0_NM = .05291772105638411
MASS_D_AMU = 2.0141


def zbl_screen(x):
    return (.18175*np.exp(-3.1998*x)+.50986*np.exp(-.94229*x)
            +.28022*np.exp(-.40290*x)+.02817*np.exp(-.20162*x))


def cm_angles(epsilon, impact, radial_order=64, screen=zbl_screen):
    nodes, weights = leggauss(radial_order)
    t = (nodes+1)/2
    u = 1-t*t
    roots = []
    for b in impact:
        equation = lambda x: 1-(b/x)**2-screen(x)/(epsilon*x)
        roots.append(brentq(equation, b, b+1/epsilon+2,
                             xtol=1e-13, rtol=1e-13))
    closest = np.asarray(roots)
    denominator = (1-(impact[:, None]*u/closest[:, None])**2
                   -screen(closest[:, None]/u)*u/(epsilon*closest[:, None]))
    if not np.isfinite(denominator).all() or np.min(denominator) <= 0:
        raise ValueError('Radial scattering quadrature crossed its turning point')
    integral = np.sum((weights/2)*(2*t)/np.sqrt(denominator), axis=1)
    angle = np.pi-2*impact/closest*integral
    if np.min(angle) < -1e-7 or np.max(angle) > np.pi+1e-7:
        raise ValueError('Scattering angle outside physical range beyond roundoff')
    return np.clip(angle, 0, np.pi)  # Only the checked <=1e-7 roundoff interval.


def ls_cross_section(E_eV, Z, scale=1.0):
    # RustBCA3.0.0 material.rs314-319 / consts.rs35:1.212 eV A^2 sqrt(E/amu).
    return scale*.01212*Z/(1+Z**(2/3))**1.5*np.sqrt(E_eV/MASS_D_AMU)


def projected_range(components, energy_max=13000., cutoff_eV=5.,
                    energy_nodes=256, impact_order=64, radial_order=64,
                    impact_max=60.):
    """components=(Z,mass_amu,number_density_nm^-3,LS_scale)."""
    energies = np.geomspace(cutoff_eV, energy_max, energy_nodes)
    projected = np.zeros(energy_nodes)
    path = np.zeros(energy_nodes)
    nodes, weights = leggauss(impact_order)
    lo, hi = np.log(1e-5), np.log(impact_max)
    impact = np.exp(lo+(nodes+1)*(hi-lo)/2)
    log_weights = weights*(hi-lo)/2
    for i in range(1, energy_nodes):
        E = energies[i]
        delta = E-energies[i-1]
        se = sum(n*ls_cross_section(E, Z, scale) for Z, mass, n, scale in components)
        diagonal = np.array([se/delta, se/delta])
        rhs = np.array([1+se/delta*projected[i-1], 1+se/delta*path[i-1]])
        for Z, mass, density, scale in components:
            a = .88534*A0_NM/(1+Z**.23)
            epsilon = a*E*mass/(MASS_D_AMU+mass)/(Z*E2_EV_NM)
            angle = cm_angles(epsilon, impact, radial_order)
            fraction = 4*MASS_D_AMU*mass/(MASS_D_AMU+mass)**2*np.sin(angle/2)**2
            remaining = E*(1-fraction)
            cosine = np.cos(np.arctan2(np.sin(angle), np.cos(angle)+MASS_D_AMU/mass))
            sigma = density*2*np.pi*a*a*impact*impact*log_weights
            current = remaining >= energies[i-1]
            alpha = np.zeros_like(remaining)
            alpha[current] = (remaining[current]-energies[i-1])/delta
            for k, (values, angular) in enumerate([(projected, cosine),
                                                  (path, np.ones_like(cosine))]):
                known = np.interp(remaining, energies[:i], values[:i], left=0.)
                known[current] = (1-alpha[current])*values[i-1]
                diagonal[k] += np.sum(sigma*(1-angular*alpha))
                rhs[k] += np.sum(sigma*angular*known)
        values = rhs/diagonal
        if not np.isfinite(values).all() or np.min(values) < 0:
            raise ValueError('Projected-range Volterra solution lost positivity')
        projected[i], path[i] = values
    return dict(energy_eV=energies, projected_nm=projected, path_nm=path)
