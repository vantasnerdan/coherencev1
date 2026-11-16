"""Czerski arXiv:2409.02112v3 Eqs 7-9/21, SI throughout.
The mobile/trapped pair density is N0*N, NOT N**2/2.
Track exposure is V*tau*Gamma_p/Gamma; it is not a free additive yield.
A narrow-resonance model is conditional on E_R and the hot-track screening.
"""
import numpy as np

QE = 1.602176634e-19
KB = 1.380649e-23
HBAR = 1.054571817e-34
MU_KG = 1.6717915e-27
BARRIERS = {'Pd': .150, 'Ti': .680, 'Zr': .413}


def reactivity_m3_s(T_K, Ue_hot_keV, E_R_eV=1.0, radius_fm=7.0):
    """Eq 21 prefactor, excluding Arrhenius and Gamma_p/Gamma (m^3/s)."""
    if min(T_K, Ue_hot_keV, E_R_eV, radius_fm) <= 0:
        raise ValueError('Thermal parameters must be positive')
    kt = KB*T_K
    penetration = np.sqrt(986.0/Ue_hot_keV)
    return (8*np.pi**1.5*HBAR**3*np.sqrt(E_R_eV*QE)
            /(MU_KG**2*radius_fm*1e-15*kt**1.5)
            *penetration*np.exp(-penetration-E_R_eV*QE/kt))


def probability(n_D_m3, *, Ue_hot_keV, E_A_eV, exposure_m3_s,
                T_K=1080., E_R_eV=1., radius_fm=7.):
    """Thermal reactions per incident D = N0*N <sigma v> V*tau*branch.

    n_D is actual total D density, and exposure_m3_s = V*tau*Gamma_p/Gamma.
    Source nominal V*tau = 1870 nm^3 * 100 ps = 1.87e-34 m^3 s.
    E_R=1 eV is a declared fit convention, not measured by these yield data;
    its amplitude degeneracy with exposure is retained, not called a lifetime.
    """
    nd = np.asarray(n_D_m3, float)
    if np.any(~np.isfinite(nd)) or np.any(nd < 0) or E_A_eV < 0 or exposure_m3_s < 0:
        raise ValueError('Nonnegative finite density/barrier/exposure required')
    return (nd**2*np.exp(-E_A_eV*QE/(KB*T_K))
            *reactivity_m3_s(T_K,Ue_hot_keV,E_R_eV,radius_fm)*exposure_m3_s)
