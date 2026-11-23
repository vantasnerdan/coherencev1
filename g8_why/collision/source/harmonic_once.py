"""One fixed harmonic spectral kernel, not a well-refinement campaign.

Benchmark hbar*omega=27.1 meV is the pinned P2 Gamma optic, not a claim about
an anharmonic local-D well. At 300 K the exact oscillator level/momentum joint
weights retain zero-point motion and level energy; no on-shell kinetic fuel.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import time

import numpy as np
import sympy as sp
from scipy.constants import Boltzmann, electron_volt
from scipy.integrate import quad
from scipy.special import gamma, roots_genlaguerre
from bound_forward import BoundSpectatorKernel
from yield_kernel import FixedProfileShellKernel


def states(omega_eV=.0271, temperature_K=300., order=32, maximum_level=32):
    x, weight = roots_genlaguerre(order, .5)
    weight /= gamma(1.5)
    q = math.exp(-omega_eV/(Boltzmann*temperature_K/electron_volt))
    coefficients = [np.ones_like(x)]
    for n in range(maximum_level):
        term = (2*x-2*n)*coefficients[n]
        if n >= 1:
            term += 3*coefficients[n-1]
        if n >= 2:
            term += (2*n+2-2*x)*coefficients[n-2]
        if n >= 3:
            term += n*coefficients[n-3]
        coefficients.append(term/(n+1))
    result = []
    for n, coefficient in enumerate(coefficients):
        if np.min(coefficient) < 0:
            raise ArithmeticError('Mehler coefficient lost positivity; do not clip it')
        probability = (1-q)**3*q**n*coefficient*weight
        result.append(dict(excitation_eV=n*omega_eV, kinetic_eV=x*omega_eV/2,
                           probabilities=probability))
    mass = sum(s['probabilities'].sum() for s in result)
    kinetic = sum(np.dot(s['probabilities'], s['kinetic_eV']) for s in result)
    expected = 3*omega_eV/4*(1+q)/(1-q)
    if abs(mass-1) > 1e-8 or abs(kinetic/expected-1) > 1e-7:
        raise ArithmeticError('Actual oscillator spectral mass/kinetic moment did not close')
    return result, dict(probability_sum=float(mass), mean_kinetic_eV=float(kinetic),
        exact_mean_kinetic_eV=expected, retained_maximum_level=maximum_level,
        radial_order=order, omega_eV=omega_eV, temperature_K=temperature_K,
        interpretation='Quantum harmonic 3D canonical spectral mixture; no fitted target heating')


def operator_checks(bound, free):
    energies = [.25,.5,1.,2.5,6.5]
    authoritative = np.array([free.common.yield_C('PdD', e) for e in energies])
    stationary = np.array([bound.shell('PdD', e, 0.) for e in energies])
    residual = float(np.max(np.abs(stationary/authoritative-1)))
    L, K = .5, .1
    a, b = (L-K)/2, math.sqrt(L*K)
    upper = min(1., a/b)
    def direct(mu):
        er = a-b*mu
        return float(free.common.DD.sigma(er, .040, 'p')*math.sqrt(2*er/L)/2)
    independent, _ = quad(direct, -1., upper, epsabs=0, epsrel=1e-9, limit=200)
    transformed = (bound.local_samples([L], K*1000)[0] @ bound.angular_weights).item()
    angular = abs(transformed/independent-1)
    p, t, mass, mu, excitation, removal = sp.symbols('p t m mu excitation removal', positive=True)
    beam, target = p*p/(2*mass), t*t/(2*mass)
    available = (beam-target)/2-p*t*mu/(2*mass)+excitation-removal
    total_minus_cm = beam+excitation-removal-(p*p+t*t+2*p*t*mu)/(4*mass)
    identity = sp.simplify(available-total_minus_cm)
    if residual > 2e-13 or angular > 1e-5 or identity != 0:
        raise ArithmeticError('Actual stationary, angular or energy-identity check failed')
    return dict(stationary_consumer_max_relative=residual, selected_threshold_angular_relative=angular,
        symbolic_energy_residual=str(identity),
        scope='Kinematics and the declared spectator on-shell matrix prescription, not an exact nuclear theorem')


def run(output):
    resource.setrlimit(resource.RLIMIT_CPU, (180,180))
    started = time.process_time()
    bound = BoundSpectatorKernel(angular_order=32)
    free = FixedProfileShellKernel(angular_order=32)
    checks = operator_checks(bound, free)
    spectral, state_metadata = states()
    energies = np.array([.25,.5,1.,2.5,6.5])
    baseline = np.array([bound.common.yield_C('PdD', e) for e in energies])
    prediction = bound.fold('PdD', energies, spectral, removal_eV=0.)
    kinetic = spectral[0]['kinetic_eV']
    probabilities = sum(s['probabilities'] for s in spectral)
    onshell = free.mixture('PdD', energies, kinetic, probabilities)
    rows = [dict(E_cm_keV=float(e), Y40=float(y0), bound_spectator_Y_over_Y40=float(yb/y0),
                 free_onshell_diagnostic_Y_over_Y40=float(yf/y0))
            for e,y0,yb,yf in zip(energies,baseline,prediction,onshell)]
    result = dict(schema='g8-harmonic-once-v1', harmonic=state_metadata, rows=rows,
        operator_checks=checks, removal_energy_eV=0.,
        approximation='Energy-conserving frozen-spectator impulse model; favorable zero removal, infinite spectator mass, unchanged +40eV on-shell DD amplitude; excludes subthreshold virtual nuclear amplitudes',
        source_frequency='P2 README, pinned2318593: Gamma collective optic27.1meV, used only as a defined harmonic benchmark, not a local-well measurement',
        model_uncertainty='Actual chemical removal energy and exact off-shell nuclear matrix are not measured by this benchmark',
        process_cpu_s=time.process_time()-started,
        process_total_CPU_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output)
