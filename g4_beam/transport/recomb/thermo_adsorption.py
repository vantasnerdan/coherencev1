"""Gas1Pa reference and explicitly coarse-grained pair-to-Langmuir free energy.

Assumption: split the specified pair adsorption free energy equally between two
independent equivalent sites in equilibrium with ideal D2 at 1 Pa, with explicit
site degeneracy; use harmonic adsorbates and rigid-rotor ideal gas, kappa=1 TST.

This is an adopted coarse pair/reservoir map, NOT a numerical differentiation
of the 1ML mean energy, nor a disjoint two-state (empty/full) pair isotherm.
ASE's classical rotational convention excludes nuclear-spin ortho/para entropy.
"""
import math
from ase import units
from thermo_harmonic import D_MASS_U, KB, harmonic_parts, positive

P_STANDARD_PA = 1.0


def gas_external_parts(temperature_K, bond_A, *, symmetry_number,
                       electronic_degeneracy):
    """D2 translational+rigid-rotor G-E_el at 1Pa, excluding ALL vibration.

    Uses D mass2.014101778u, I=(m_D/2)*r**2, homonuclear symmetry normally2,
    singlet electronic degeneracy normally1; both choices explicit, not fitted.
    """
    for value, name in [(temperature_K, 'temperature_K'), (bond_A, 'bond_A'),
                        (symmetry_number, 'symmetry_number'),
                        (electronic_degeneracy, 'electronic_degeneracy')]:
        positive(value, name)
    kt_j = units._k * temperature_K
    kt = KB * temperature_K
    mass_kg = 2 * D_MASS_U * units._amu
    inertia = D_MASS_U * units._amu / 2 * (bond_A * 1e-10)**2
    qtrans = (2 * math.pi * mass_kg * kt_j / units._hplanck**2)**1.5 * kt_j / P_STANDARD_PA
    qrot = 8 * math.pi**2 * inertia * kt_j / (symmetry_number * units._hplanck**2)
    # H_trans+rot=(5/2+1)kT cancels the +5/2,+1 entropy constants.
    return dict(translation_G_eV=-kt * math.log(qtrans),
                rotation_G_eV=-kt * math.log(qrot),
                electronic_G_eV=-kt * math.log(electronic_degeneracy),
                G_external_eV=-kt * math.log(qtrans * qrot * electronic_degeneracy),
                mass_D2_kg=mass_kg, inertia_kg_m2=inertia,
                standard_pressure_Pa=P_STANDARD_PA)


def gas1pa(temperature_K, bond_A, vibration_eV, *, symmetry_number,
           electronic_degeneracy):
    """G_D2(T,1Pa)-E_D2,el in eV/molecule, with one supplied stable stretch."""
    positive(vibration_eV, 'D2 stretch h*nu in eV')
    external = gas_external_parts(temperature_K, bond_A,
        symmetry_number=symmetry_number, electronic_degeneracy=electronic_degeneracy)
    vib = harmonic_parts(temperature_K, [vibration_eV])
    return dict(**external, vibration=vib,
                correction_G_eV=external['G_external_eV'] + vib['F_eV'])


def adsorption_from_pair(temperature_K, electronic_adsorption_pair_eV,
                         adsorbed_stable_eV, clean_stable_eV, gas,
                         *, occupied_site_degeneracy, empty_site_degeneracy):
    """Convert E(slab+2D)-E(clean)-E(D2) to adopted finite-T Eads>0 binding.

    Adsorbed/clean spectra must use matched Pd mobile-space backgrounds; adsorbed
    adds six D Cartesian modes. TST initial modes can be reused for adsorbed.
    Configurational mixing is NOT included in harmonic F. Site degeneracy enters
    once as kT ln(g_occ/g_empty), separate from TST path/basin degeneracies.
    """
    if not math.isfinite(electronic_adsorption_pair_eV):
        raise ValueError('Finite pair adsorption electronic energy required')
    positive(occupied_site_degeneracy, 'occupied_site_degeneracy')
    positive(empty_site_degeneracy, 'empty_site_degeneracy')
    if len(adsorbed_stable_eV) != len(clean_stable_eV) + 6:
        raise ValueError('Two adsorbed D require six additional Cartesian modes')
    if gas['standard_pressure_Pa'] != P_STANDARD_PA:
        raise ValueError('Gas reference must be 1Pa')
    ads = harmonic_parts(temperature_K, adsorbed_stable_eV)
    clean = harmonic_parts(temperature_K, clean_stable_eV)
    dg = electronic_adsorption_pair_eV + ads['F_eV'] - clean['F_eV'] - gas['correction_G_eV']
    degeneracy = KB * temperature_K * math.log(occupied_site_degeneracy / empty_site_degeneracy)
    return dict(adsorption_pair_G_eV=dg, pair_to_site_G_eV=dg / 2,
                degeneracy_binding_eV_per_D=degeneracy,
                Eads_eV_per_D=-dg / 2 + degeneracy,
                standard_pressure_Pa=P_STANDARD_PA)


def langmuir_coverage(temperature_K, mu_relative_half_gas1pa_eV, Eads_eV_per_D):
    """No beam stripping: theta/(1-theta)=exp[(mu+Eads)/kT]."""
    positive(temperature_K, 'temperature_K')
    x = (mu_relative_half_gas1pa_eV + Eads_eV_per_D) / (KB * temperature_K)
    if x >= 0:
        return 1 / (1 + math.exp(-x))
    e = math.exp(x)
    return e / (1 + e)
