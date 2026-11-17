"""Coverage-dependent sheet chemical potential and finite-site detailed balance."""
import math
from sheet_thermo import SheetThermo
from thermo_boundary import d2_impingement_per_Pa, effective_sticking
from thermo_harmonic import KB as THERMO_KB


def supplied_sheet(thermodynamics, temperature_K, capture_s0):
    sheet = thermodynamics(temperature_K)
    if not isinstance(sheet, SheetThermo) or sheet.temperature_K != temperature_K:
        raise ValueError('Source-bound sheet thermodynamics at this temperature required')
    if not isinstance(capture_s0, (int, float)) or not 0 < capture_s0 <= 1:
        raise ValueError('Explicit intrinsic capture probability in (0,1] required')
    return sheet


def coverage_from_mu(mu_eV, sheet, kT_eV):
    return sheet.theta_from_mu(mu_eV, kT_eV)


def local_boundary(model, theta, *, sheet, capture_s0):
    """Map the differential, not mean, binding to the original atom boundary."""
    if not 0 < theta < 1 or sheet.temperature_K != model.T:
        raise ValueError('Interior sheet coverage at the actual foil temperature required')
    energy = sheet.differential_binding(theta)
    gamma = d2_impingement_per_Pa(model.T)
    rdes = gamma * capture_s0 * math.exp(-2 * energy / model.kT) * theta**2
    result = effective_sticking(model.T, theta, energy, rdes,
        boundary_gamma_per_Pa=model.G1, boundary_kT_eV=model.kT)
    rads = gamma * capture_s0 * model.p_gas * (1 - theta)**2
    result.update(theta_ML=theta, differential_binding_eV_per_D=energy,
                  mean_binding_eV_per_D=sheet.mean_binding(theta),
                  capture_s0=capture_s0, D2_desorption_events_m2_s=rdes,
                  D2_adsorption_events_m2_s=rads,
                  D_atom_net_outflux_m2_s=2 * (rdes - rads),
                  availability_normalized_s=gamma / model.G1 * capture_s0 * (1 - theta)**2,
                  source_kT_eV=model.kT, accepted_thermo_kT_eV=THERMO_KB * model.T,
                  kB_relative_convention_difference=model.kT / (THERMO_KB * model.T) - 1)
    return result


def gas_equilibrated(model, *, sheet, capture_s0):
    if model.c['loaded']:
        raise ValueError('Loaded foil is bulk-fed, not gas equilibrium')
    mu = .5 * model.kT * math.log(model.p_gas)
    theta = coverage_from_mu(mu, sheet, model.kT)
    result = local_boundary(model, theta, sheet=sheet, capture_s0=capture_s0)
    result.update(interpretation='gas-equilibrated no-bulk regular sheet',
                  beam_stimulated_desorption=False,
                  theta_solver='Derivative of excess plus ideal configurational free energy',
                  coverage_target_ML=1e-4,
                  balance_relative_residual=abs(result['D_atom_net_outflux_m2_s']) /
                    (2 * (result['D2_desorption_events_m2_s'] + result['D2_adsorption_events_m2_s'])),
                  theta_roundtrip_ML=abs(theta - coverage_from_mu(
                    result['mu_relative_half_gas1pa_eV'], sheet, model.kT)))
    return result
