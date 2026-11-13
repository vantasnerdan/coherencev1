"""Evaluate accepted laws, not new symbolic searches or dynamics.

PySR threshold: Fleet entry_194. SINDy features/coefficients: entry_210.
The latter is a post-slaving near-threshold reduction, not a global closure.
"""
import math


SINDY_LINEAR = 0.9702868696288892
SINDY_ALPHA_CUBIC = 1.3428214004150867
SINDY_BETA_CUBIC = -0.9638929752453904


def features(m):
    rho, r = float(m.rho), float(m.r)
    theta = float(m.temperature / m.r)
    w = math.tanh(1 / (2 * theta))
    inverse_threshold_squared = m.f * w - m.eta
    g = float(m.g)
    stiffness = float(m.stiffness)
    A = 4 * g * g * w / r - stiffness
    D = 2 * m.kappa + (4 * g * g * w / r) * m.gamma / (r * r + m.gamma * m.gamma)
    B = 32 * g**4 / r**3 * (w - (1 - w * w) / (2 * theta))
    alpha, beta = A / D, B / D
    margin = 4 * g * g * w / (r * stiffness) - 1
    outside = []
    for name, value, lower, upper in (
        ("r", r, 0.5, 1.5), ("theta", theta, 0.1, 0.6),
        ("f", m.f, 0.4, 1), ("z_eff", m.z_eff, 0.5, 1),
        ("kappa", m.kappa, 0.2, 0.8), ("gamma", m.gamma, 0.1, 0.8),
        ("absolute_threshold_margin", abs(margin), 0.002, 0.015),
    ):
        if not lower <= value <= upper:
            outside.append(name)
    return dict(
        rho=rho, r=r, theta=theta, w=w,
        inverse_rho_c_squared=inverse_threshold_squared,
        rho_c=1 / math.sqrt(inverse_threshold_squared) if inverse_threshold_squared > 0 else None,
        ordered_by_exact_threshold=inverse_threshold_squared * rho * rho > 1,
        alpha=alpha, beta=beta,
        sindy_linear_coefficient=SINDY_LINEAR * alpha,
        sindy_cubic_coefficient=SINDY_ALPHA_CUBIC * alpha + SINDY_BETA_CUBIC * beta,
        outside_declared_sindy_parameter_ranges=outside,
        sindy_scope="Features/coefficients only; no new time trajectory or surrogate power. Empty range failures would still not prove fast slaving.",
    )
