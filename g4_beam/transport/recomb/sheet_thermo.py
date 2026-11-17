"""Anchored regular-solution sheet free energy; all energies in eV per site/D.

Mean binding B(theta)=B1+w(1-theta) defines g_ex=-theta*B(theta).
The differential binding is B1+w-2w*theta, NOT the mean binding.
Ideal configurational entropy is added once, outside the DFT harmonic data.
"""
from dataclasses import dataclass
import math
from scipy.optimize import brentq
from scipy.special import expit


@dataclass(frozen=True)
class SheetThermo:
    temperature_K: float
    full_mean_binding_eV: float
    interaction_eV: float
    source_id: str

    def __post_init__(self):
        if (not math.isfinite(self.temperature_K) or self.temperature_K <= 0
                or not math.isfinite(self.full_mean_binding_eV)
                or not math.isfinite(self.interaction_eV) or not self.source_id):
            raise ValueError('Finite source-bound sheet thermodynamics required')

    def mean_binding(self, theta):
        return self.full_mean_binding_eV + self.interaction_eV * (1 - theta)

    def differential_binding(self, theta):
        return self.full_mean_binding_eV + self.interaction_eV * (1 - 2 * theta)

    def excess_free_energy(self, theta):
        return -theta * self.mean_binding(theta)

    def chemical_potential(self, theta, kT_eV):
        if not 0 < theta < 1 or kT_eV <= 0:
            raise ValueError('Interior coverage and positive kT required')
        return (kT_eV * (math.log(theta) - math.log1p(-theta))
                - self.differential_binding(theta))

    def theta_from_mu(self, mu_eV, kT_eV):
        if not math.isfinite(mu_eV) or not math.isfinite(kT_eV) or kT_eV <= 0:
            raise ValueError('Finite chemical potential and positive kT required')
        if 2 * self.interaction_eV + 4 * kT_eV <= 0:
            raise ValueError('Nonconvex sheet requires a phase-coexistence model')
        a = self.full_mean_binding_eV + self.interaction_eV
        twice_w = 2 * self.interaction_eV
        lower = (mu_eV + a - max(0., twice_w)) / kT_eV - 1
        upper = (mu_eV + a - min(0., twice_w)) / kT_eV + 1
        def residual(log_odds):
            return kT_eV * log_odds + twice_w * expit(log_odds) - a - mu_eV
        root = brentq(residual, lower, upper, xtol=1e-12)
        theta = float(expit(root))
        if not 0 < theta < 1:
            raise ArithmeticError('Sheet coverage is beyond representable interior')
        return theta

    def dtheta_dmu(self, theta, kT_eV):
        if not 0 < theta < 1:
            raise ValueError('Interior coverage required for susceptibility')
        return 1 / (2 * self.interaction_eV + kT_eV / (theta * (1 - theta)))

    def free_energy(self, theta, kT_eV):
        if not 0 <= theta <= 1:
            raise ValueError('Coverage outside physical site interval')
        mixing = 0. if theta in (0, 1) else (
            theta * math.log(theta) + (1 - theta) * math.log1p(-theta))
        return self.excess_free_energy(theta) + kT_eV * mixing


def fit_anchored(temperature_K, full_mean_binding_eV, lower_means, *,
                 source_id, residual_target_eV=.01):
    """Fit w to the supplied (.25,.5) means, with the qualified1ML point exact."""
    if {float(theta) for theta, _ in lower_means} != {.25, .5} or len(lower_means) != 2:
        raise ValueError('Both declared lower-coverage free energies required')
    differences = [(1 - theta, binding - full_mean_binding_eV)
                   for theta, binding in lower_means]
    interaction = math.fsum(x * y for x, y in differences) / math.fsum(
        x * x for x, _ in differences)
    sheet = SheetThermo(temperature_K, full_mean_binding_eV, interaction, source_id)
    residuals = [dict(coverage_ML=theta, supplied_binding_eV=binding,
                     fitted_binding_eV=sheet.mean_binding(theta),
                     residual_eV=binding - sheet.mean_binding(theta))
                 for theta, binding in lower_means]
    maximum = max(abs(row['residual_eV']) for row in residuals)
    if maximum > residual_target_eV:
        raise ValueError(f'Regular-sheet fit residual {maximum:g}eV exceeds {residual_target_eV:g}; model decision required')
    return sheet, dict(interaction_eV=interaction, residuals=residuals,
                       max_residual_eV=maximum, residual_target_eV=residual_target_eV)
