"""Source-bound symbolic/numerical evidence for the two issue-243 APIs.

Run explicitly with the original rendered Huke page, PDF, and installed
chi0_drude.py supplied as arguments. No GPAW import, material data or Native
execution is involved. Hashes bind the inspected source, not an automated
claim that OCR or a hash proves its interpretation. Ordinary high precision
is numerical evidence, not a rigorous floating-point enclosure.

Analytic argument: expand atanh(t)=sum t^(2n+1)/(2n+1) for |t|<1.
In the u>1 branch t=1/u, cancellation of the constant leaves
F^2=sum_(n>=1) t^(2n)/((2n-1)(2n+1)), a strictly positive series.
For t<=1/2 its factored 24-term truncation has remainder bounded by
 t^48/(49*51*(1-t^2)) < 2e-18.
The two one-sided limits at u=1 are finite despite a divergent derivative.

Independence: the numerical oracle evaluates the dimensionful printed
Eq. (10) at adaptive precision, not the production polynomial. The gas
normalization is derived from phase-space/angular integration and dN/dE,
not from restating the RMS API. Synthetic finite-width multiband moments
then check its meaning beyond that zero-temperature gas limit.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from g4_beam.screen.dynamic._upstream.electronic_speed import (
    fermi_rms_velocity_au,
    speed_factor,
)


SOURCE_BINDINGS = {
    "huke_pdf": (653081, "7a428728bc7d86c12bc8aa9bb659c57bb74d9070ec87d91ba0dbb0aea74cd33b"),
    "huke_page": (252672, "948f313fa1a3ea9728d490df93a5536ffc79bed824240bb9f2fd9041f418d50c"),
    "drude_source": (7715, "5ff1fa7ef2facf1dc4624a2bd54fcae72a52a33a9149807edc89207385a32eb2"),
}


def symbolic_evidence():
    u, v, vf, t = sp.symbols("u v vf t", positive=True)
    half = sp.Rational(1, 2)
    low = half + (1-u**2)/(4*u) * sp.log((1+u)/(1-u))
    high = half + (1-u**2)/(4*u) * sp.log((u+1)/(u-1))
    dimensional_low = half + (vf**2-v**2)/(4*vf*v) * sp.log((v+vf)/(vf-v))
    dimensional_high = half + (vf**2-v**2)/(4*vf*v) * sp.log((v+vf)/(v-vf))
    assert sp.simplify(dimensional_low.subs(v, u*vf) - low) == 0
    assert sp.simplify(dimensional_high.subs(v, u*vf) - high) == 0
    limits = [sp.limit(low, u, 0, dir="+"),
              sp.limit(low, u, 1, dir="-"),
              sp.limit(high, u, 1, dir="+")]
    assert limits == [1, half, half]
    reciprocal = sp.simplify(high.subs(u, 1/t))
    series = sp.series(reciprocal, t, 0, 10).removeO().expand()
    expected = sum(t**(2*n)/((2*n-1)*(2*n+1)) for n in range(1, 5))
    assert sp.simplify(series - expected) == 0
    assert sp.limit(reciprocal/t**2, t, 0, dir="+") == sp.Rational(1, 3)
    n = sp.symbols("n", integer=True, positive=True)
    assert sp.simplify((1/(2*n-1)-1/(2*n+1))/2 - 1/((2*n-1)*(2*n+1))) == 0

    # Spin degeneracy g is explicit; normalized k integration is per cell.
    k, kf, omega, g = sp.symbols("k kf Omega g", positive=True)
    theta, phi = sp.symbols("theta phi", real=True)
    measure = g*omega/(2*sp.pi)**3
    number = measure*4*sp.pi*sp.integrate(k**2, (k, 0, kf))
    ef = kf**2/2
    dos = sp.diff(number, kf)/sp.diff(ef, kf)
    angular_zz = sp.integrate(sp.cos(theta)**2*sp.sin(theta),
                              (theta, 0, sp.pi), (phi, 0, 2*sp.pi))
    assert angular_zz == 4*sp.pi/3
    # delta(k^2/2-EF) contributes the radial Jacobian 1/kF.
    plasma_zz = (4*sp.pi/omega)*measure*angular_zz*kf**4/kf
    assert sp.simplify(plasma_zz - 4*sp.pi*number/omega) == 0
    rms_squared = sp.simplify(omega*3*plasma_zz/(4*sp.pi*dos))
    assert sp.simplify(rms_squared - kf**2) == 0
    return {
        "squared_endpoint_limits": [str(x) for x in limits],
        "high_speed_squared_series": str(series),
        "positive_root_high_speed_limit": "u*F(u) -> 1/sqrt(3)",
        "free_gas_N": str(number), "free_gas_D": str(dos),
        "free_gas_plasma_diagonal": str(sp.simplify(plasma_zz)),
        "free_gas_rms_squared": str(rms_squared),
    }


def printed_equation_oracle(ratio):
    if ratio == 0:
        return mp.mpf(1)
    if ratio == 1:
        return mp.sqrt(mp.mpf("0.5"))
    # The printed ratio/log then cancellation loses up to three log10(u) digits.
    digits = max(100, int(4*abs(np.log10(ratio))) + 100)
    with mp.workdps(digits):
        vf = mp.mpf(7)/5
        v = mp.mpf(float(ratio))*vf
        squared = mp.mpf("0.5") + (vf**2-v**2)/(4*vf*v) * mp.log(abs((v+vf)/(v-vf)))
        assert squared > 0
        return +mp.sqrt(squared)


def speed_evidence():
    points = np.unique(np.concatenate((
        np.logspace(-300, 300, 121),
        [0, np.nextafter(0., 1.), 0.125, 0.5,
         np.nextafter(1., 0.), 1., np.nextafter(1., 2.),
         np.nextafter(2., 0.), 2., np.nextafter(2., 3.),
         1e8, 1e12, 1e16, np.finfo(float).max],
    )))
    actual = speed_factor(points)
    expected = np.array([float(printed_equation_oracle(x)) for x in points])
    assert np.all(np.isfinite(actual)) and np.all(actual > 0)
    np.testing.assert_allclose(actual, expected, rtol=3e-14, atol=0)
    np.testing.assert_array_equal(speed_factor(points.reshape(1, -1)), actual.reshape(1, -1))
    assert isinstance(speed_factor(0.5), float)
    assert np.all(np.diff(actual) <= 2e-15)
    relative_error = np.max(np.abs((actual-expected)/expected))

    probe = 0.5
    correct = float(printed_equation_oracle(probe))
    logarithmic_term = (1-probe**2)/(4*probe)*np.log((1+probe)/(1-probe))
    mutants = {
        "missing_positive_root": 0.5 + logarithmic_term,
        "wrong_sign": np.sqrt(0.5 - logarithmic_term),
        "half_log_prefactor": np.sqrt(0.5 + logarithmic_term/2),
    }
    for value in mutants.values():
        assert not np.isclose(value, correct, rtol=1e-6, atol=0)
    huge = 1e16
    naive = np.sqrt(0.5 + (1-huge**2)/(4*huge)*np.log(abs((1+huge)/(1-huge))))
    assert not np.isclose(naive, float(printed_equation_oracle(huge)), rtol=1e-6, atol=0)
    return {"samples": len(points), "rtol": 3e-14, "atol": 0,
            "max_relative_error": float(relative_error),
            "exposed_mutations": list(mutants) + ["direct_high_speed_cancellation"]}


def velocity_evidence():
    gas_errors = []
    for kf, omega, spin in [(0.4, 5., 1), (1.7, 23., 2), (3.2, 11., 4)]:
        number = spin*omega*kf**3/(6*np.pi**2)
        dos = spin*omega*kf/(2*np.pi**2)
        plasma = np.eye(3)*(4*np.pi*number/omega)
        actual = fermi_rms_velocity_au(plasma, dos, omega)
        np.testing.assert_allclose(actual, kf, rtol=2e-15, atol=0)
        gas_errors.append(abs(actual/kf-1))

    # Two synthetic bands, three normalized k points and spin degeneracy 2.
    energies = np.array([-0.7, -0.1, 0.25, 0.6, 1.1, 1.8])
    velocities = np.array([[1., 2., 0.], [-1., 0., 3.], [0., 1., 2.],
                           [2., -1., 1.], [0.5, 2., -1.], [3., 1., 0.5]])
    width, chemical_potential, omega = 0.3, 0.2, 17.
    weights = np.tile(np.array([0.2, 0.3, 0.5]), 2)*2
    f = 1/(1+np.exp((energies-chemical_potential)/width))
    derivative_weights = weights*f*(1-f)/width
    dos = derivative_weights.sum()
    moment = np.einsum("n,ni,nj->ij", derivative_weights, velocities, velocities)
    plasma = 4*np.pi/omega*moment
    actual = fermi_rms_velocity_au(plasma, dos, omega)
    with mp.workdps(80):
        numerator, denominator = mp.mpf(0), mp.mpf(0)
        for energy, velocity, weight in zip(energies, velocities, weights):
            occupation = 1/(1+mp.exp((mp.mpf(float(energy))-mp.mpf(chemical_potential))/mp.mpf(width)))
            dw = mp.mpf(float(weight))*occupation*(1-occupation)/mp.mpf(width)
            denominator += dw
            numerator += dw*sum(mp.mpf(float(component))**2 for component in velocity)
        expected = float(mp.sqrt(numerator/denominator))
    np.testing.assert_allclose(actual, expected, rtol=3e-15, atol=0)
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    np.testing.assert_allclose(fermi_rms_velocity_au(rotation@plasma@rotation.T, dos, omega),
                               expected, rtol=3e-15, atol=0)
    shifted_f = 1/(1+np.exp(((energies+7)-(chemical_potential+7))/width))
    np.testing.assert_allclose(shifted_f, f, rtol=3e-15, atol=0)
    mutants = {"missing_4pi": actual*np.sqrt(4*np.pi),
               "directional_mean_not_trace": actual/np.sqrt(3),
               "missing_volume": actual/np.sqrt(omega),
               "spin_mismatched_DOS": actual*np.sqrt(2),
               "plasma_frequency_not_squared": fermi_rms_velocity_au(np.diag(np.sqrt(np.diag(plasma))), dos, omega)}
    for value in mutants.values():
        assert not np.isclose(value, expected, rtol=1e-6, atol=0)
    return {"free_gas_max_relative_error": max(gas_errors),
            "synthetic_multiband_relative_error": abs(actual/expected-1),
            "exposed_mutations": list(mutants),
            "limitation": "Finite-width multiband isotropic closure; not absolute DFT EF or 18 free electrons."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--huke-pdf", type=Path, required=True)
    parser.add_argument("--huke-page", type=Path, required=True)
    parser.add_argument("--drude-source", type=Path, required=True)
    args = parser.parse_args()
    bindings = {}
    for name, expected in SOURCE_BINDINGS.items():
        raw = getattr(args, name).read_bytes()
        binding = (len(raw), hashlib.sha256(raw).hexdigest())
        if binding != expected:
            raise ValueError(f"Source binding mismatch: {name}")
        bindings[name] = {"bytes": binding[0], "sha256": binding[1]}
    result = {"source_bindings": bindings, "symbolic": symbolic_evidence(),
              "speed_numeric": speed_evidence(), "velocity_numeric": velocity_evidence(),
              "scope": "Algebraic/sampled numerical kernel evidence only; no material U, Native, full dynamical response or physical adequacy claim."}
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
