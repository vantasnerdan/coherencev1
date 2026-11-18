"""Huke Eq. (10) speed closure and an intraband band-RMS velocity.

These are proposal kernels, not a shell-resolved dynamical DFT response.
Huke et al., arXiv:0805.4538v1, rendered p. 10, Eq. (10), gives
U_dyn**2 / U_ad**2, so the positive square root is required here.
"""

from __future__ import annotations

import numpy as np


# For t = 1/u <= 1/2, F(u)^2 = t^2 sum_n t^(2n)/((2n+1)(2n+3)).
# All terms are positive. The omitted tail after 24 terms is bounded by
# t^48 / (49*51*(1-t^2)) in the polynomial, < 2e-18 at t = 1/2.
_HIGH_SPEED_COEFFICIENTS = tuple(
    1.0 / ((2 * n + 1) * (2 * n + 3)) for n in range(24)
)


def speed_factor(u):
    """Return U_dyn/U_ad for finite real speed ratios u >= 0.

    Scalars return a float; array-like inputs preserve their shape. The
    continuous values are F(0)=1 and F(1)=1/sqrt(2). For u -> infinity,
    F(u) ~ 1/(sqrt(3)*u), evaluated without subtracting nearly equal terms
    or squaring a potentially enormous u. No cancellation-error clipping
    is used, including when the result is subnormal in binary64.
    """
    if np.iscomplexobj(u):
        raise ValueError("u must be finite, real and nonnegative")
    ratios = np.asarray(u, dtype=float)
    if not np.all(np.isfinite(ratios)) or np.any(ratios < 0):
        raise ValueError("u must be finite, real and nonnegative")

    result = np.empty_like(ratios)
    result[ratios == 0] = 1.0
    result[ratios == 1] = np.sqrt(0.5)

    low = (ratios > 0) & (ratios < 1)
    x = ratios[low]
    result[low] = np.sqrt(
        0.5 + 0.5 * (1 - x) * (1 + x) * (np.arctanh(x) / x)
    )

    middle = (ratios > 1) & (ratios < 2)
    x = ratios[middle]
    result[middle] = np.sqrt(
        0.5 - 0.5 * ((x - 1) * (x + 1) / x) * np.arctanh(1 / x)
    )

    high = ratios >= 2
    reciprocal = 1 / ratios[high]
    polynomial = np.polynomial.polynomial.polyval(
        reciprocal * reciprocal, _HIGH_SPEED_COEFFICIENTS
    )
    result[high] = reciprocal * np.sqrt(polynomial)
    return result.item() if result.ndim == 0 else result


def fermi_rms_velocity_au(
    plasma_squared_au, dos_per_Hartree_cell, volume_Bohr3
):
    """Return the finite-width intraband RMS speed in atomic units.

    ``plasma_squared_au`` is the actual 3x3 intraband omega_p^2 tensor
    (Hartree**2), NOT its elementwise square root, an interband optical
    plasma fit, or an eV-valued tensor. ``dos_per_Hartree_cell`` is
    D_T = sum_(n,k,s) w_k [-df/dE] with normalized k weights and spin
    counted once, using the SAME bands and occupation derivative as the
    tensor. For FD occupations -df/dE = f(1-f)/width in Hartree units.
    Volume is the corresponding cell volume in Bohr**3.

    With M_ij = sum w_k [-df/dE] v_i v_j, the intraband relation is
    omega_p,ij^2 = (4*pi/Omega)*M_ij. Hence
    v_rms^2 = trace(M)/D_T = Omega*trace(omega_p^2)/(4*pi*D_T).
    The trace, not its one-direction average, gives the total RMS speed.
    For an isotropic parabolic gas at T=0, D=3N/k_F^2 and each diagonal
    omega_p^2=4*pi*N/Omega, giving v_rms=k_F=v_F (electron mass=1).

    This is an isotropic closure velocity for multiband metals, not an
    absolute DFT Fermi energy or a density of 18 free electrons. Localized
    bound shells can contribute static screening without intraband weight;
    applying one Huke factor to their screening is an explicit model
    limitation, not a shell-resolved or full frequency-dependent result.
    """
    tensor = np.asarray(plasma_squared_au, dtype=complex)
    if tensor.shape != (3, 3) or not np.all(np.isfinite(tensor)):
        raise ValueError("plasma_squared_au must be a finite 3x3 tensor")
    if np.iscomplexobj(dos_per_Hartree_cell) or np.iscomplexobj(volume_Bohr3):
        raise ValueError("DOS and volume must be finite positive real scalars")
    dos = np.asarray(dos_per_Hartree_cell, dtype=float)
    volume = np.asarray(volume_Bohr3, dtype=float)
    if (
        dos.ndim != 0 or volume.ndim != 0
        or not np.isfinite(dos) or not np.isfinite(volume)
        or dos <= 0 or volume <= 0
    ):
        raise ValueError("DOS and volume must be finite positive real scalars")
    tolerance = 64 * np.finfo(float).eps * np.max(np.abs(tensor))
    if np.max(np.abs(tensor - tensor.conj().T)) > tolerance:
        raise ValueError("plasma_squared_au must be Hermitian")
    if np.min(np.linalg.eigvalsh(tensor)) < -tolerance:
        raise ValueError("plasma_squared_au must be positive semidefinite")
    trace = np.trace(tensor).real
    if trace < 0:
        raise ValueError("plasma_squared_au must have nonnegative trace")
    return float(np.sqrt(volume) * np.sqrt(trace / (4 * np.pi)) / np.sqrt(dos))
