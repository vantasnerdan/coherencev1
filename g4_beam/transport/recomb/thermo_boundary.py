"""Convert an occupied-pair event rate to transport.py's D-atom boundary.

n_pairs=n_sites*z*reactive_bond_fraction/2 counts UNORDERED bonds. z=6 counts
all nearest neighbours on a triangular Pd111 site lattice; z=1 counts a chosen
disjoint pairing. The 1x2 cell's one pair/area is the LATTER, not all bonds.
Use a separate sum of channels when their barriers/occupancies differ.
"""
import math
from ase import units
from thermo_adsorption import P_STANDARD_PA
from thermo_harmonic import D_MASS_U, KB, positive


def d2_impingement_per_Pa(temperature_K):
    """Molecular collisions/(m2*s*Pa), isotope-correct SI mass."""
    positive(temperature_K, 'temperature_K')
    return 1 / math.sqrt(2 * math.pi * 2 * D_MASS_U * units._amu * units._k * temperature_K)


def pair_density(site_density_m2, coordination, reactive_bond_fraction):
    positive(site_density_m2, 'site_density_m2')
    positive(coordination, 'coordination')
    if not 0 <= reactive_bond_fraction <= 1:
        raise ValueError('reactive_bond_fraction outside [0,1]')
    return site_density_m2 * coordination * reactive_bond_fraction / 2


def desorption_flux(rate_per_occupied_pair_s, *, site_density_m2, coordination,
                    reactive_bond_fraction, occupied_pair_probability):
    """Explicit pair statistics; theta**2 only for independent sites."""
    if not math.isfinite(rate_per_occupied_pair_s) or rate_per_occupied_pair_s < 0:
        raise ValueError('Finite nonnegative pair event rate required')
    if not 0 <= occupied_pair_probability <= 1:
        raise ValueError('occupied_pair_probability outside [0,1]')
    density = pair_density(site_density_m2, coordination, reactive_bond_fraction)
    events = density * occupied_pair_probability * rate_per_occupied_pair_s
    return dict(unordered_pair_density_m2=density, D2_events_m2_s=events,
                D_atom_outflux_m2_s=2 * events)


def effective_sticking(temperature_K, theta, Eads_eV_per_D, desorption_D2_events_m2_s,
                       *, boundary_gamma_per_Pa, boundary_kT_eV):
    """Detailed balance in ORIGINAL boundary convention, no clipping of s.

    Given the thermodynamic (pre-beam-stripping) theta, p_eq=1Pa*exp(2mu/kT).
    R_ads=R_des*p_gas/p_eq; J_D=2*(R_des-R_ads)
          =2*s_eff*Gamma_boundary*(p_eq-p_gas).
    Pass actual model.gamma1(T) for byte-pinned transport normalization; its
    rounded D2 mass/constants differ slightly from isotope-exact impingement.
    Pass actual model.kT too: the pinned model rounds kB independently of ASE.
    This is a local effective coefficient, not a coverage-independent s0.
    """
    positive(temperature_K, 'temperature_K')
    positive(boundary_gamma_per_Pa, 'boundary_gamma_per_Pa')
    if not 0 < theta < 1:
        raise ValueError('Interior thermodynamic coverage required')
    if not math.isfinite(desorption_D2_events_m2_s) or desorption_D2_events_m2_s < 0:
        raise ValueError('Finite nonnegative molecular event flux required')
    positive(boundary_kT_eV, 'boundary_kT_eV')
    kt = boundary_kT_eV
    mu = kt * (math.log(theta) - math.log1p(-theta)) - Eads_eV_per_D
    peq = P_STANDARD_PA * math.exp(2 * mu / kt)
    s = desorption_D2_events_m2_s / (boundary_gamma_per_Pa * peq)
    return dict(mu_relative_half_gas1pa_eV=mu, equilibrium_pressure_Pa=peq,
                s_effective=float(s), within_probability_interval=bool(0 <= s <= 1),
                boundary_gamma_per_Pa=boundary_gamma_per_Pa)
