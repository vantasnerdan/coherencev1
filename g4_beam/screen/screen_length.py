"""Small-q cubic RPA length and explicitly Yukawa-model derived quantities."""
import itertools

import numpy as np
from ase.units import Bohr, Hartree

from screen_io import require

COULOMB_EV_A = Hartree*Bohr


def fit_length(q_Ainv, epsilon_macro):
    q = np.asarray(q_Ainv, float)
    epsilon = np.asarray(epsilon_macro, complex)
    require(q.ndim == 1 and len(q) in [2, 3] and epsilon.shape == q.shape,
            'Length requires two or three scalar small-q samples')
    order = np.argsort(q)
    q, epsilon = q[order], epsilon[order]
    require(np.isfinite(epsilon).all() and np.all(q > 0) and np.all(np.diff(q) > 0),
            'Nonfinite, zero, or duplicate physical q')
    require(np.max(np.abs(epsilon.imag)) < 1e-7*np.max(np.abs(epsilon.real)),
            'Static imaginary response is not negligible')
    require(np.all(epsilon.real > 1), 'Metallic Yukawa length is undefined')
    x = q*q
    y = x*(epsilon.real-1)

    def intercept(n):
        design = np.column_stack([np.ones(n), x[:n]/x[n-1]])
        return float(np.linalg.lstsq(design, y[:n], rcond=None)[0][0])

    kappa2 = intercept(len(q))
    require(kappa2 > 0, 'Nonpositive metallic small-q intercept')
    length = 1/np.sqrt(kappa2)
    variants = [float(length), *list(1/np.sqrt(y))]
    if len(q) == 3:
        short = intercept(2)
        require(short > 0, 'Two-smallest-q intercept is nonpositive')
        variants.append(float(1/np.sqrt(short)))
    spread = max(abs(v/length-1) for v in variants)
    return dict(lambda_A=float(length), kappa_squared_Ainv2=kappa2,
                q_Ainv=q.tolist(), q_squared_epsilon_minus_one_Ainv2=y.tolist(),
                fit='q^2(epsilon_M-1)=kappa^2+regular*q^2',
                q_window_lambda_span_A=[min(variants), max(variants)],
                q_window_relative_spread=float(spread),
                q_window_pass_10percent=bool(spread <= .1),
                Ue_static_Yukawa_proxy_eV=float(COULOMB_EV_A/length),
                model='Valence-PAW static RPA low-q length; C/lambda is a Yukawa proxy, not a full-RPA near-origin shift')


def nearest_periodic_distance(positions_A, cell_A):
    positions = np.asarray(positions_A, float)
    cell = np.asarray(cell_A, float)
    require(len(positions) > 0 and np.linalg.det(cell) > 0,
            'No sites or invalid periodic cell')
    closest = float('inf')
    for shift in itertools.product([-1, 0, 1], repeat=3):
        translation = np.asarray(shift)@cell
        for i, origin in enumerate(positions):
            for j, target in enumerate(positions):
                if i == j and shift == (0, 0, 0):
                    continue
                distance = float(np.linalg.norm(target+translation-origin))
                require(distance > 1e-8, 'Coincident registered sites')
                closest = min(closest, distance)
    return closest


def ratio_and_force(length_A, DD_A, metal_A):
    require(length_A > 0 and DD_A > 0 and metal_A > 0, 'Nonpositive length or spacing')
    x = DD_A/length_A
    return dict(lambda_over_DD=length_A/DD_A, lambda_over_metal=length_A/metal_A,
                DD_spacing_A=DD_A, metal_spacing_A=metal_A,
                Yukawa_radial_force_factor=float((1+x)*np.exp(-x)),
                force_convention='(1+d/lambda)exp(-d/lambda); not a longitudinal Hessian or a second screening of DFT forces')
