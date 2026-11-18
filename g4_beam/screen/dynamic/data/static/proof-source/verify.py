#!/usr/bin/env python3
"""Bounded, pure-CPU evidence for the proposed Coherencev1/#242 kernel.

Required raw Yan text/PDF are SHA-bound; no GPAW, Native, network, or material
results are consumed.  The full complex finite-matrix oracle integrates the
observation direction, rather than repeating the kernel's j0 contraction.
Continuum samples use a normalized spherical reciprocal-volume quadrature;
Gaussian real-space Coulomb and independent radial TF integrals fix its units,
weights and prefactor.  All preregistered radii participate in every curve.

Run with PYTHONPATH=src and --yan-text PATH --yan-pdf PATH.  Optional
--mutation half-prefactor / zero-phase deliberately changes the tested q
kernel; the SAME independent baseline assertions must exit nonzero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Callable

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.special import erf

from substrate_framework import microscopic_screening as microscopic
from substrate_framework.verification import CheckLedger

HERE = Path(__file__).resolve().parent
YAN_BINDINGS = {
    "yan-paw-response.txt": (
        76624, "bd246c3bf3517278c39996e26bf3b3227f57ce65a5f62c9807961da4b76221fb"
    ),
    "yan-paw-response.pdf": (
        443238, "44df983b77c96798d9527af1f84163acd1015f780898f1f255c0b4b2d590a6f8"
    ),
}
QKernel = Callable[..., np.ndarray]


def binding(path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def bind_sources(checks: CheckLedger, text: Path, pdf: Path) -> dict[str, object]:
    inventory = json.loads((HERE / "source_inventory.json").read_text())
    entries = {Path(row["path"]).name: row for row in inventory["sources"]}
    bindings = {}
    for name, supplied in (("yan-paw-response.txt", text), ("yan-paw-response.pdf", pdf)):
        actual = binding(supplied)
        size, digest = YAN_BINDINGS[name]
        expected = {"bytes": size, "sha256": digest}
        checks.check(f"raw primary binding: {name}", actual == expected)
        checks.check(
            f"inventory agrees with primary: {name}",
            {key: entries[name][key] for key in expected} == expected,
        )
        bindings[name] = actual
    for path in (
        Path(microscopic.__file__), Path(__file__), HERE / "preregistration.json",
        HERE / "source_inventory.json",
    ):
        bindings[path.name] = binding(path)
    return bindings


def angular_oracle(
    Q: np.ndarray, inverse: np.ndarray, source: np.ndarray, radii: np.ndarray,
    volume: float, weight: float, kappa: float, width: float, coulomb: float,
) -> np.ndarray:
    """Direct two-charge Fourier potential, integrated over the second charge."""
    mu, mu_weight = leggauss(24)
    phi = np.arange(48) * (2.0 * math.pi / 48)
    directions = np.stack(
        [
            (np.sqrt(1 - mu[:, None] ** 2) * np.cos(phi)).ravel(),
            (np.sqrt(1 - mu[:, None] ** 2) * np.sin(phi)).ravel(),
            np.broadcast_to(mu[:, None], (24, 48)).ravel(),
        ], axis=1,
    )
    angular_weights = np.repeat(mu_weight / (2.0 * 48), 48)
    norms = np.sqrt(np.sum(Q * Q, axis=1))
    values = np.zeros(len(radii), dtype=complex)
    for g in range(len(Q)):
        for h in range(len(Q)):
            coefficient = (1.0 if g == h else 0.0) - inverse[g, h]
            if g == h:
                coefficient -= kappa**2 / (norms[g] ** 2 + kappa**2) * math.exp(
                    -(norms[g] / width) ** 2
                )
            source_phase = np.exp(1j * np.dot(Q[g] - Q[h], source))
            for index, radius in enumerate(radii):
                observation = np.exp(1j * radius * (directions @ Q[g]))
                values[index] += (
                    coefficient * source_phase * np.dot(angular_weights, observation)
                    / (norms[g] * norms[h])
                )
    return 4.0 * math.pi * coulomb * weight / volume * values.real


def check_local_fields(checks: CheckLedger, q_kernel: QKernel, radii: np.ndarray) -> dict:
    Q = np.array([[0.5, 0.2, -0.1], [1.4, -0.4, 0.3], [-0.6, 1.1, 0.8]])
    response = np.array(
        [[0.42, 0.04 + 0.09j, -0.03 + 0.02j],
         [0.04 - 0.09j, 0.25, 0.02 - 0.05j],
         [-0.03 - 0.02j, 0.02 + 0.05j, 0.18]]
    )
    inverse = np.eye(3) - response
    source = np.array([0.31, -0.27, 0.46])
    parameters = (17.0, 0.23, 1.7, 2.6, 17.3)
    expected = angular_oracle(Q, inverse, source, radii, *parameters)
    actual = q_kernel(Q, inverse, source, radii, *parameters)
    error = float(np.max(np.abs(actual - expected)))
    checks.check("full complex finite-radius two-charge angular oracle", error < 2e-12)
    gaps = {}
    mutations = {
        "no source phase": (inverse, np.zeros(3)),
        "reversed source phase": (inverse, -source),
        "discarded imaginary local fields": (inverse.real, source),
        "diagonal-only local fields": (np.diag(np.diag(inverse)), source),
    }
    for name, (matrix, site) in mutations.items():
        changed = q_kernel(Q, matrix, site, radii, *parameters)
        gaps[name] = float(np.max(np.abs(changed - expected)))
        checks.check(f"angular oracle exposes {name}", gaps[name] > 1e-3)
    checks.check(
        "opposite-q conjugation preserves the real readout",
        np.max(np.abs(q_kernel(-Q, inverse.conj(), source, radii, *parameters) - actual))
        < 2e-12,
    )
    permutation = np.array([2, 0, 1])
    checks.check(
        "joint reciprocal-basis permutation preserves the readout",
        np.max(np.abs(q_kernel(
            Q[permutation], inverse[np.ix_(permutation, permutation)], source,
            radii, *parameters,
        ) - actual)) < 2e-12,
    )
    for name, factor in (("self-energy half", 0.5), ("missing 4pi", 1 / (4 * math.pi))):
        gaps[name] = float(np.max(np.abs(factor * actual - expected)))
        checks.check(f"independent angular oracle exposes {name}", gaps[name] > 0.01)
    return {"angular_max_error_eV": error, "mutation_gaps_eV": gaps}


def radial_nodes(order: int, cutoff: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Uniform reciprocal-ball measure, not a hand-inserted radial prefactor."""
    nodes, gauss_weights = leggauss(order)
    Q = 0.5 * cutoff * (nodes + 1.0)
    dQ = 0.5 * cutoff * gauss_weights
    ball_volume = (4.0 * math.pi / 3.0) * cutoff**3
    volume = (2.0 * math.pi) ** 3 / ball_volume
    weights = 4.0 * math.pi * Q**2 * dQ / ball_volume
    return Q, weights, volume


def sampled_curve(
    q_kernel: QKernel, radii: np.ndarray, response: Callable[[float], float],
    kappa: float, width: float, coulomb: float, order: int = 192,
    cutoff: float = 24.0, mutation: str = "none",
) -> np.ndarray:
    nodes, weights, volume = radial_nodes(order, cutoff)
    if mutation == "equal weights":
        weights = np.full(order, 1.0 / order)
    if mutation == "twice volume":
        volume *= 2.0
    if not np.isclose(np.sum(weights), 1.0, rtol=0, atol=3e-14):
        raise AssertionError("independent ball quadrature weights are not normalized")
    residual = np.zeros(len(radii))
    for Q, weight in zip(nodes, weights):
        residual += q_kernel(
            [[Q, 0.0, 0.0]], [[1.0 - response(float(Q))]], [0.37, -0.29, 0.41],
            radii, volume, float(weight), kappa, width, coulomb,
        )
    reference = microscopic.reference_eV(radii, kappa, width, coulomb)
    if mutation == "repeated addback":
        reference *= order
    return residual + reference


def gaussian_oracle(radii: np.ndarray, width: float, coulomb: float) -> np.ndarray:
    result = np.full(len(radii), coulomb * width / math.sqrt(math.pi))
    positive = radii > 0
    result[positive] = coulomb * erf(width * radii[positive] / 2.0) / radii[positive]
    return result


def radial_tf_oracle(radii: np.ndarray, kappa: float, width: float, coulomb: float):
    """Adaptive oscillatory Fourier route, independent of the heat-kernel API."""
    values, errors = [], []
    for radius in radii:
        integral, error = quad(
            lambda Q: (kappa**2 / (Q**2 + kappa**2))
            * math.exp(-(Q / width) ** 2)
            * (math.sin(Q * radius) / (Q * radius) if Q * radius else 1.0),
            0.0, np.inf, epsabs=2e-12, epsrel=2e-12, limit=200,
        )
        values.append(2.0 * coulomb / math.pi * integral)
        errors.append(2.0 * coulomb / math.pi * error)
    return np.array(values), max(errors)


def check_continuum(checks: CheckLedger, q_kernel: QKernel, radii: np.ndarray) -> dict:
    coulomb, physical_width = 17.3, 2.3
    gaussian = lambda Q: math.exp(-(Q / physical_width) ** 2)
    expected = gaussian_oracle(radii, physical_width, coulomb)
    errors, reference_errors = {}, {}
    for kappa in (1.0, 2.0, 3.0):
        for width in (1.5, 3.0):
            name = f"kappa={kappa},width={width}"
            fine = sampled_curve(q_kernel, radii, gaussian, kappa, width, coulomb)
            coarse = sampled_curve(
                q_kernel, radii, gaussian, kappa, width, coulomb, order=96
            )
            shorter = sampled_curve(
                q_kernel, radii, gaussian, kappa, width, coulomb, cutoff=18.0
            )
            error = float(np.max(np.abs(fine - expected)))
            errors[name] = error
            checks.check(f"Gaussian Coulomb analytic curve: {name}", error < 2e-8)
            checks.check(
                f"separate radial refinement and cutoff: {name}",
                np.max(np.abs(coarse - fine)) < 2e-7
                and np.max(np.abs(shorter - fine)) < 2e-8,
            )
            direct, quad_error = radial_tf_oracle(radii, kappa, width, coulomb)
            heat = microscopic.reference_eV(radii, kappa, width, coulomb)
            ref_error = float(np.max(np.abs(heat - direct)))
            reference_errors[name] = {"max_eV": ref_error, "quad_error_eV": quad_error}
            checks.check(
                f"independent damped-TF Fourier integral: {name}",
                ref_error < 2e-9 and quad_error < 2e-9,
            )
    # A different physical response, not an identity reference-minus-itself.
    kappa, width = 1.7, 2.6
    physical_tf = lambda Q: kappa**2 / (Q**2 + kappa**2) * math.exp(-(Q / width) ** 2)
    target_tf, _ = radial_tf_oracle(radii, kappa, width, coulomb)
    tf = sampled_curve(q_kernel, radii, physical_tf, 3.0, 1.5, coulomb)
    checks.check("nonidentical-reference TF continuum reconstruction", np.max(abs(tf - target_tf)) < 2e-8)
    gaps = {}
    for mutation in ("equal weights", "twice volume", "repeated addback"):
        changed = sampled_curve(
            q_kernel, radii, gaussian, 2.0, 3.0, coulomb, mutation=mutation
        )
        gaps[mutation] = float(np.max(np.abs(changed - expected)))
        checks.check(f"analytic continuum exposes {mutation}", gaps[mutation] > 0.1)
    # Removing the Gaussian damping approaches the Yukawa/TF closed form;
    # |j0|<=1 bounds the error by the omitted contact spectral weight.
    tf_limit = np.full(len(radii), coulomb * kappa)
    positive = radii > 0
    tf_limit[positive] = -coulomb * np.expm1(-kappa * radii[positive]) / radii[positive]
    for broad_width in (40.0, 160.0):
        damped = microscopic.reference_eV(radii, kappa, broad_width, coulomb)
        contact_bound = 2.0 * coulomb * kappa**2 / (math.sqrt(math.pi) * broad_width)
        checks.check(
            f"analytic undamped TF limit and uniform error bound: {broad_width}",
            np.max(abs(damped - tf_limit)) <= contact_bound * (1.0 + 1e-8),
        )
    return {"gaussian_errors_eV": errors, "reference_errors": reference_errors,
            "continuum_mutation_gaps_eV": gaps}


def check_edges(checks: CheckLedger, radii: np.ndarray) -> None:
    args = [[[0.5, 0.0, 0.0]], [[0.6]], [0.2, 0.1, -0.3], radii,
            17.0, 0.25, 1.7, 2.6, 17.3]
    invalid = {
        "Gamma": (0, [[0.0, 0.0, 0.0]]),
        "matrix dimension": (1, np.eye(2)),
        "nonfinite matrix": (1, [[np.nan]]),
        "source dimension": (2, [0.0, 0.0]),
        "negative radius": (3, [-0.1]),
        "scalar radius": (3, 0.0),
        "empty radii": (3, []),
        "nonfinite radius": (3, [np.inf]),
        "zero volume": (4, 0.0),
        "zero weight": (5, 0.0),
        "unnormalized individual weight": (5, 1.1),
        "zero kappa": (6, 0.0),
        "negative width": (7, -1.0),
        "complex units": (8, 1.0 + 1j),
    }
    for name, (index, value) in invalid.items():
        candidate = list(args)
        candidate[index] = value
        try:
            microscopic.q_residual_eV(*candidate)
        except ValueError:
            checks.check(f"explicit invalid-input refusal: {name}", True)
        else:
            checks.check(f"explicit invalid-input refusal: {name}", False)
    unordered = np.array([1.0, 0.0, 1e-9, 0.1, 0.0])
    actual = microscopic.reference_eV(unordered, 1.7, 2.6, 17.3)
    direct, _ = radial_tf_oracle(unordered, 1.7, 2.6, 17.3)
    checks.check("contact limit, duplicate and unsorted radii", np.max(abs(actual - direct)) < 2e-9)
    before = [np.array(value, copy=True) for value in args[:4]]
    microscopic.q_residual_eV(*args)
    checks.check("q kernel does not mutate its inputs", all(
        np.array_equal(np.asarray(value), old) for value, old in zip(args[:4], before)
    ))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--yan-text", type=Path, required=True)
    parser.add_argument("--yan-pdf", type=Path, required=True)
    parser.add_argument("--mutation", choices=("none", "half-prefactor", "zero-phase"), default="none")
    args = parser.parse_args()
    checks = CheckLedger("Coherencev1-screen242/proposed-kernel")
    bindings = bind_sources(checks, args.yan_text, args.yan_pdf)
    plan = json.loads((HERE / "preregistration.json").read_text())
    radii = np.array(plan["static_readout_proposal_to_check"]["radii_A"], dtype=float)

    def q_kernel(*values):
        if args.mutation == "zero-phase":
            values = list(values)
            values[2] = np.zeros(3)
        result = microscopic.q_residual_eV(*values)
        return 0.5 * result if args.mutation == "half-prefactor" else result

    evidence = {"bindings": bindings, "radii_A": radii.tolist(), "mutation": args.mutation}
    evidence["local_fields"] = check_local_fields(checks, q_kernel, radii)
    evidence["continuum"] = check_continuum(checks, q_kernel, radii)
    check_edges(checks, radii)
    print(json.dumps(evidence, sort_keys=True, allow_nan=False))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())
