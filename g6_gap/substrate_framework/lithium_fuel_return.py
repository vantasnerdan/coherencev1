"""Finite-domain D-on-6LiT arithmetic with explicit stopping inputs.

Bosch/Hale, Nucl.Fusion32(1992)620-621, Eqs8-9/TableIV.
Electronic baseline: Geant4v11.3.0 G4BraggModel::ElectronicStoppingPower,
ICRU49/Ziegler table, rows H and Li. Nuclear loss is the unchanged ZBL
formula imported from coherencev1@47f0902 g4_beam/forward/model.py.
These are model imports, not measured compound LiT stopping or a threshold.
"""
from __future__ import annotations

import numpy as np

__all__ = ['dt_eq7_sigma_mb', 'lit_icru_electronic', 'lit_inherited_sqrt_electronic', 'lit_zbl_nuclear',
           'stopping_per_formula', 'thick_target_return', 'break_even_barns']

MASS_D = 2.01410  # Keep the inherited forward-model projectile mass convention.
MASS_T = 3.0160492779  # NIST isotope table:3.0160492779(24)u, atomic mass.
MASS_LI6 = 6.0151228874  # NIST isotope table:6.0151228874(16)u, atomic mass.
DT_Q_KEV = 17590.0
DT_BG = 34.3827


def _positive_energy(energy):
    energy = np.asarray(energy, dtype=float)
    if not np.isfinite(energy).all() or np.any(energy <= 0):
        raise ValueError('Finite positive deuteron lab energies required')
    return energy


def dt_eq7_sigma_mb(energy_lab_keV, screening_keV=0.):
    """Requested Eq7 form using the source D-T S-factor, in millibarns.

    TableIV fit domain is0.5<=E_cm<=550keV. Evaluations outside that domain
    are mathematical extrapolation; the caller must label them as such.
    This formula never shifts the S-factor or treats it as a measurement.
    """
    lab = _positive_energy(energy_lab_keV)
    if not np.isfinite(screening_keV) or screening_keV < 0:
        raise ValueError('Nonnegative finite screening energy required')
    e = lab*MASS_T/(MASS_D+MASS_T)
    numerator = 6.927e4+e*(7.454e8+e*(2.050e6+e*5.2002e4))
    denominator = 1+e*(63.8+e*(-.995+e*(6.981e-5+e*1.728e-4)))
    s = numerator/denominator
    return s/np.sqrt(e*(e+screening_keV))*np.exp(-DT_BG/np.sqrt(e+screening_keV))

def lit_inherited_sqrt_electronic(energy_lab_keV):
    """Unchanged forward Se_atom convention extended by the tabulated Li A1.

    The low-velocity coefficient law was applied to the entire energy range
    in the inherited consumer. Its continuation above10keV/u is explicitly
    a baseline extrapolation, not the full source parametrization.
    """
    return (1.411+1.254)*np.sqrt(_positive_energy(energy_lab_keV)*1.007276/MASS_D)



def lit_icru_electronic(energy_lab_keV):
    """Bragg H+Li electronic cross section, eV/(1e15formula/cm2).

    Equal-velocity proton energy divided by its mass gives E_D/M_D in
    keV/u. Low-energy metallic sqrt extrapolation is deliberately retained
    as the comparator, not imposed on the first-principles stopping.
    """
    t = _positive_energy(energy_lab_keV)/MASS_D
    x = np.maximum(t, 10.)
    factor = np.sqrt(np.minimum(t/10., 1.))
    total = np.zeros_like(t)
    for slow_coefficient, high_coefficient, inverse_energy, energy_term in (
            (1.440, 242.6, 12000., .1159),
            (1.600, 725.6, 3013., .04578)):
        slow = slow_coefficient*x**.45
        high = high_coefficient/x*np.log1p(inverse_energy/x+energy_term*x)
        total += slow*high/(slow+high)*factor
    return total


def _zbl(energy, target_z, target_mass):
    a = 1+target_z**.23
    eps = 32.53*target_mass*energy/(target_z*(MASS_D+target_mass)*a)
    low = np.log1p(1.1383*eps)/(2*(eps+.01321*eps**.21226+.19593*np.sqrt(eps)))
    high = np.log(eps)/(2*eps)
    reduced = np.where(eps <= 30, low, high)
    return 8.462*target_z*MASS_D*reduced/((MASS_D+target_mass)*a)


def lit_zbl_nuclear(energy_lab_keV):
    """D-on-Li6 plus D-on-T ZBL loss per LiT formula, unchanged law."""
    energy = _positive_energy(energy_lab_keV)
    return _zbl(energy, 3, MASS_LI6)+_zbl(energy, 1, MASS_T)


def stopping_per_formula(electronic_eV_A, volume_per_formula_A3):
    """Convert eV/angstrom to eV/(1e15formula/cm2), one T per formula."""
    if not np.isfinite(volume_per_formula_A3) or volume_per_formula_A3 <= 0:
        raise ValueError('Positive formula-unit volume required')
    electronic = np.asarray(electronic_eV_A, dtype=float)
    if not np.isfinite(electronic).all() or np.any(electronic < 0):
        raise ValueError('Passive finite electronic stopping required')
    return electronic*volume_per_formula_A3/10


def thick_target_return(energy_lab_keV, total_stopping_per_formula, screening_keV=0.):
    """Return G on the supplied finite domain, integrating from its first energy.

    Retains the Target E_cut convention; no integral below the supplied
    floor, hidden stopping join, interpolation, or low-energy extrapolation
    is supplied here. Callers own the measured/computed bridge and report
    compound/Bragg, electronic-model and fit-domain eligibility separately.
    """
    energy = _positive_energy(energy_lab_keV)
    stopping = np.asarray(total_stopping_per_formula, dtype=float)
    if (energy.ndim != 1 or energy.size < 2 or np.any(np.diff(energy) <= 0)
            or stopping.shape != energy.shape or not np.isfinite(stopping).all()
            or np.any(stopping <= 0)):
        raise ValueError('Increasing energy grid and matched positive stopping required')
    # mb/(eV*1e-15cm2) times1000eV/keV = 1e-9 per keV.
    integrand = dt_eq7_sigma_mb(energy, screening_keV)/stopping
    probability = np.empty_like(energy)
    probability[0] = 0.
    np.cumsum(np.diff(energy)*(integrand[1:]+integrand[:-1])/2.,
              out=probability[1:])
    probability *= 1e-9
    return DT_Q_KEV/energy*probability


def break_even_barns(total_stopping_per_formula):
    """Cross section for local G=1, one T nucleus per LiT formula."""
    return np.asarray(total_stopping_per_formula, dtype=float)*1e6/DT_Q_KEV
