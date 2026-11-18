"""Static bare-minus-screened *mutual* unit-charge potential (eV).

Proposed Coherencev1/#242 readout, not an accepted material prediction.
Yan et al., arXiv:1104.1273v1, Eq. (4) fixes the Fourier measure to
sum_q w_q / Omega; Eq. (8) gives v_G = 4*pi*C/|q+G|**2.  With the
Coulomb-symmetrized inverse E = v**(-1/2) epsilon_inverse v**(1/2),
W = epsilon_inverse v = sqrt(v) E sqrt(v).  Thus the mutual kernel is
v-W = sqrt(v) (I-E) sqrt(v), with no self-energy factor of one half.
Evaluating at R0+r*n and R0 gives exp(+i*Q_G.R0) exp(-i*Q_G'.R0);
the spherical average over n replaces exp(i*Q_G.r*n) by j0(|Q_G|r).
Yan Eq. (11) and its local-field discussion require the full complex E;
PAW augmentation (Eq. (15)-(18)) belongs to the supplied response, not to
an additional correction here.

Subtract the diagonal Gaussian-damped Thomas-Fermi reference in each q
contraction, sum normalized q weights, then add reference_eV ONCE.
All reciprocal vectors must be nonzero: no Gamma value is inferred or
silently zeroed.  Callers own matrix provenance, inversion quality, q
normalization across calls, and finite-basis/q-grid convergence.  C is
injected explicitly in eV Angstrom; no material or constants are fitted.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import quad
from scipy.special import erf, erfcx

__all__ = ["q_residual_eV", "reference_eV"]


def _real_array(value: object, name: str) -> np.ndarray:
    if np.iscomplexobj(value):
        raise ValueError(f"{name} must be real")
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite")
    return array


def _positive(value: object, name: str) -> float:
    array = _real_array(value, name)
    if array.ndim != 0 or array <= 0:
        raise ValueError(f"{name} must be a positive scalar")
    return float(array)


def _radii(value: object) -> np.ndarray:
    radii = _real_array(value, "radii_A")
    if radii.ndim != 1 or radii.size == 0 or np.any(radii < 0):
        raise ValueError("radii_A must be a nonempty nonnegative radius vector")
    return radii


def q_residual_eV(
    Q_Ainv: object,
    inverse_epsilon_sym: object,
    source_A: object,
    radii_A: object,
    volume_A3: float,
    weight: float,
    kappa_Ainv: float,
    width_Ainv: float,
    coulomb_eV_A: float,
) -> np.ndarray:
    """Return one weighted q residual, in the input radius order.

    Q_Ainv is (N,3), source_A is (3,), and inverse_epsilon_sym is (N,N)
    in exactly that Q ordering.  The supplied inverse is Coulomb-symmetrized,
    NOT the unsymmetrized Yan Eq. (11) matrix or a macroscopic head.  A
    finite-broadening complex matrix is allowed; take the real mutual
    potential only after contraction.  Positive weight <= 1 is one member
    of a normalized BZ quadrature; volume_A3 is the primitive cell volume.

    u = exp(-i Q.R0)/|Q|; result = 4*pi*C*w/Omega times
    Re[(j0(|Q|r)*u)^dagger ((I-E)-diag(f)) u], where
    f = kappa**2/(Q**2+kappa**2) exp[-(Q/width)**2].
    Only one matrix-vector product is needed; no dense identity, reference
    matrix, or radius-by-N workspace is constructed.  Inputs are not mutated.
    """
    Q = _real_array(Q_Ainv, "Q_Ainv")
    source = _real_array(source_A, "source_A")
    radii = _radii(radii_A)
    volume = _positive(volume_A3, "volume_A3")
    q_weight = _positive(weight, "weight")
    kappa = _positive(kappa_Ainv, "kappa_Ainv")
    width = _positive(width_Ainv, "width_Ainv")
    coulomb = _positive(coulomb_eV_A, "coulomb_eV_A")
    if Q.ndim != 2 or Q.shape[1] != 3 or Q.shape[0] == 0:
        raise ValueError("Q_Ainv must have shape (N,3) with N > 0")
    if source.shape != (3,):
        raise ValueError("source_A must have shape (3,)")
    if q_weight > 1:
        raise ValueError("weight must be <= 1 in a normalized q quadrature")
    inverse = np.asarray(inverse_epsilon_sym)
    if inverse.shape != (len(Q), len(Q)) or not np.all(np.isfinite(inverse)):
        raise ValueError("inverse_epsilon_sym must be a finite (N,N) matrix")
    norms = np.hypot(np.hypot(Q[:, 0], Q[:, 1]), Q[:, 2])
    if np.any(norms == 0):
        raise ValueError("all Q_Ainv vectors must be nonzero; Gamma needs quadrature")

    u = np.exp(-1j * (Q @ source)) / norms
    response = u - inverse @ u
    with np.errstate(over="ignore"):
        reference = np.exp(-np.square(norms / width)) / (
            1.0 + np.square(norms / kappa)
        )
    response -= reference * u
    radial_density = (u.conj() * response).real
    prefactor = 4.0 * math.pi * coulomb * q_weight / volume
    result = np.empty(radii.size)
    for index, radius in enumerate(radii):
        result[index] = prefactor * np.dot(
            np.sinc(norms * radius / math.pi), radial_density
        )
    return result


def reference_eV(
    radii_A: object,
    kappa_Ainv: float,
    width_Ainv: float,
    coulomb_eV_A: float,
) -> np.ndarray:
    """Return the infinite-continuum reference, to add once after q summation.

    The Fourier measure gives (2*C/pi) integral_0^infinity f(Q)j0(Qr)dQ.
    At contact this is C*kappa*erfcx(kappa/width).  For r>0 the equivalent
    positive heat-kernel integral avoids oscillatory Q quadrature: introduce
    kappa**2/(Q**2+kappa**2) = integral_0^infinity exp(-s)
    exp(-s*Q**2/kappa**2) ds, and integrate the Gaussian j0 transform.
    The result is C/r integral exp(-s) erf(kappa*r/(2*sqrt(s+a**2))) ds,
    a=kappa/width.  Evaluating erf(z)/z keeps the r->0 limit well conditioned.
    Numerical integration errors are surfaced, not replaced by a fallback.
    """
    radii = _radii(radii_A)
    kappa = _positive(kappa_Ainv, "kappa_Ainv")
    width = _positive(width_Ainv, "width_Ainv")
    coulomb = _positive(coulomb_eV_A, "coulomb_eV_A")
    a = kappa / width
    result = np.empty(radii.size)
    for index, radius in enumerate(radii):
        if radius == 0:
            result[index] = coulomb * kappa * erfcx(a)
            continue

        def integrand(s: float) -> float:
            denominator = math.hypot(math.sqrt(s), a)
            z = (0.5 * kappa * radius) / denominator
            if z < 1e-4:
                z2 = z * z
                erf_over_z = (2.0 / math.sqrt(math.pi)) * (
                    1.0 - z2 / 3.0 + z2 * z2 / 10.0
                )
            else:
                erf_over_z = erf(z) / z
            return math.exp(-s) * erf_over_z / (2.0 * denominator)

        integral, error = quad(
            integrand, 0.0, np.inf, epsabs=2e-12, epsrel=2e-11, limit=200
        )
        if not math.isfinite(integral) or error > max(2e-11, abs(integral) * 2e-10):
            raise ArithmeticError("continuum reference quadrature did not converge")
        result[index] = coulomb * kappa * integral
    return result
