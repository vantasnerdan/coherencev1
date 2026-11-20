"""Finite-q observables of a Coulomb-symmetrized microscopic inverse.

Q=q+G is Cartesian inverse length, not fractional q. Full matrix inversion
must precede these contractions. Callers own BZ completeness/convergence.
Campillo et al., cond-mat/9808113v1, Eqs. 3.11--3.15; Yan et al.,
arXiv:1104.1273v1, Eqs. 4, 8, 11. No fitted homogeneous response is used.
"""
from __future__ import annotations

import numpy as np

__all__ = ["mutual_screening_q", "powder_stopping_q", "metallic_screening_length"]


def _vectors(Q):
    Q = np.asarray(Q, dtype=float)
    if Q.ndim != 2 or Q.shape[1] != 3 or not np.isfinite(Q).all():
        raise ValueError("Q must be finite Cartesian (N,3)")
    norms = np.linalg.norm(Q, axis=1)
    if not len(Q) or np.any(norms == 0):
        raise ValueError("Nonzero Q required; Gamma needs an explicit quadrature")
    return Q, norms


def _prefactor(volume, weight, coulomb):
    if not (np.isfinite([volume, weight, coulomb]).all()
            and volume > 0 and 0 < weight <= 1 and coulomb > 0):
        raise ValueError("Positive volume/C and normalized individual weight required")
    return 4 * np.pi * coulomb * weight / volume


def mutual_screening_q(Q, inverse_sym, source, radii, *, volume, weight, coulomb):
    """One q's spherical bare-minus-screened mutual potential contribution.

    Lengths in Angstrom, volume in Angstrom^3, C in eV Angstrom; output eV.
    No self-energy half. This is a finite-basis contribution, NOT a total or
    a continuum tail correction. Sum a complete normalized q rule externally.
    """
    Q, norms = _vectors(Q)
    inverse = np.asarray(inverse_sym)
    source = np.asarray(source, dtype=float)
    radii = np.asarray(radii, dtype=float)
    if (inverse.shape != (len(Q), len(Q)) or not np.isfinite(inverse).all()
            or source.shape != (3,) or not np.isfinite(source).all()
            or radii.ndim != 1 or not np.isfinite(radii).all()
            or np.any(radii < 0)):
        raise ValueError("Malformed matrix, source or radii")
    u = np.exp(-1j * (Q @ source)) / norms
    density = (u.conj() * (u - inverse @ u)).real
    prefactor = _prefactor(volume, weight, coulomb)
    return np.array([prefactor * np.dot(np.sinc(norms * r / np.pi), density)
                     for r in radii])


def powder_stopping_q(Q, energy, loss_wG, hbar_v, *, volume, weight, coulomb):
    """One q's random stopping, averaged over projectile directions.

    loss_wG=-Im(inverse_sym_wGG.diagonal); LF enter via the full inversion.
    energy is hbar*omega in eV, increasing from zero. hbar_v is eV Angstrom.
    Returns eV/Angstrom. For each G evaluate integral_0^1 mu L(Q*v*mu)dmu
    exactly for a piecewise-linear energy interpolant; no angular sampling,
    clipping, extrapolation, or forced low-speed metallic coefficient.

    The all-Q Fourier prefactor is 4pi, equivalent to the positive-energy
    8pi convention. Its homogeneous limit is 2/(pi*v^2) integral dQ/Q
    integral_0^(Q*v) omega*L d(omega) in atomic units.
    """
    _, norms = _vectors(Q)
    energy = np.asarray(energy, dtype=float)
    loss = np.asarray(loss_wG, dtype=float)
    speeds = np.atleast_1d(np.asarray(hbar_v, dtype=float))
    if (energy.ndim != 1 or len(energy) < 2 or energy[0] != 0
            or not np.isfinite(energy).all() or np.any(np.diff(energy) <= 0)
            or loss.shape != (len(energy), len(norms))
            or not np.isfinite(loss).all() or speeds.ndim != 1
            or not np.isfinite(speeds).all() or np.any(speeds <= 0)):
        raise ValueError("Malformed frequency/loss/speed grid")
    if np.any(norms.max() * speeds > energy[-1]):
        raise ValueError("Doppler interval exceeds saved frequency range")
    # Integral of energy*piecewise-linear(loss), retaining the sign as saved.
    slopes = np.diff(loss, axis=0) / np.diff(energy)[:, None]
    x = energy[:-1, None]
    dx = np.diff(energy)[:, None]
    increments = (x * loss[:-1] * dx
                  + (loss[:-1] + x * slopes) * dx**2 / 2
                  + slopes * dx**3 / 3)
    cumulative = np.concatenate((np.zeros((1, len(norms))),
                                 np.cumsum(increments, axis=0)))
    prefactor = _prefactor(volume, weight, coulomb)
    result = np.empty(len(speeds))
    g = np.arange(len(norms))
    for index, speed in enumerate(speeds):
        upper = norms * speed
        interval = np.searchsorted(energy, upper, side="right") - 1
        interval = np.minimum(interval, len(energy) - 2)
        left = energy[interval]
        width = upper - left
        y = loss[interval, g]
        slope = slopes[interval, g]
        integral = (cumulative[interval, g] + left * y * width
                    + (y + left * slope) * width**2 / 2
                    + slope * width**3 / 3)
        result[index] = prefactor * np.sum(integral / (upper**2 * norms))
    return result


def metallic_screening_length(q_Ainv, epsilon_macro, *, coulomb):
    """Small-q head-intercept length; never an insulator or short-range fit.

    Preserve the pinned screen_length.fit_length convention:
    q²(epsilon_M-1)=kappa²+regular*q², with independent two-smallest-q
    and pointwise window spread. The injected C/lambda is explicitly a
    Yukawa-model low-Q proxy, not a microscopic contact potential.
    """
    q = np.asarray(q_Ainv, dtype=float)
    epsilon = np.asarray(epsilon_macro, dtype=complex)
    if (q.ndim != 1 or len(q) not in (2, 3) or epsilon.shape != q.shape
            or not np.isfinite(q).all() or not np.isfinite(epsilon).all()
            or not np.isfinite(coulomb) or coulomb <= 0):
        raise ValueError("Two or three finite q/head samples and positive C required")
    order = np.argsort(q)
    q, epsilon = q[order], epsilon[order]
    if np.any(q <= 0) or np.any(np.diff(q) <= 0):
        raise ValueError("Distinct positive physical q required")
    if np.max(abs(epsilon.imag)) >= 1e-7 * np.max(abs(epsilon.real)):
        raise ValueError("Static imaginary response is not negligible")
    if np.any(epsilon.real <= 1):
        raise ValueError("Metallic Yukawa length is undefined")
    x = q*q
    y = x*(epsilon.real-1)

    def intercept(n):
        design = np.column_stack([np.ones(n), x[:n]/x[n-1]])
        return float(np.linalg.lstsq(design, y[:n], rcond=None)[0][0])

    kappa2 = intercept(len(q))
    if kappa2 <= 0:
        raise ValueError("Nonpositive metallic small-q intercept")
    length = float(1/np.sqrt(kappa2))
    variants = [length, *list(1/np.sqrt(y))]
    if len(q) == 3:
        short = intercept(2)
        if short <= 0:
            raise ValueError("Nonpositive two-smallest-q intercept")
        variants.append(float(1/np.sqrt(short)))
    spread = max(abs(v/length-1) for v in variants)
    return dict(lambda_A=length, kappa_squared_Ainv2=kappa2,
                q_Ainv=q.tolist(), q_squared_epsilon_minus_one_Ainv2=y.tolist(),
                fit="q^2(epsilon_M-1)=kappa^2+regular*q^2",
                q_window_lambda_span_A=[min(variants), max(variants)],
                q_window_relative_spread=float(spread),
                q_window_pass_10percent=bool(spread <= .1),
                Ue_static_Yukawa_proxy_eV=float(coulomb/length),
                model="Valence-PAW static RPA low-q length; C/lambda is a Yukawa proxy, not near-origin screening")
