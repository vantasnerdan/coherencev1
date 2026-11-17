"""Bounded mathematical-fixture proof of the derivative sheet and frozen-v closure.

No lower-coverage material energies are supplied or inferred by this script.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
from pathlib import Path
import resource
import signal
import time
import traceback
import numpy as np
import sympy as sp
from scipy.optimize import brentq, minimize_scalar
from sheet_thermo import SheetThermo, fit_anchored
from coverage_vacancies import read_run_end
from coverage_inputs import load_source, operating_points, make_model, binding
from coverage_boundary import gas_equilibrated, local_boundary
from coverage_consume import consume_cases


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def algebra_and_sensitivity(source, rows):
    theta, a, w, kt = sp.symbols('theta a w kt', positive=True)
    free = -a*theta+w*theta**2+kt*(theta*sp.log(theta)+(1-theta)*sp.log(1-theta))
    mu = sp.diff(free, theta)
    require(sp.simplify(mu-(-a+2*w*theta+kt*(sp.log(theta)-sp.log(1-theta)))) == 0,
            'Chemical potential is not the free-energy derivative')
    differential = a-2*w*theta
    ratio = (theta/(1-theta))**2*sp.exp(-2*differential/kt)
    odds = sp.symbols('odds', positive=True)
    require(sp.simplify((ratio/sp.exp(2*mu/kt)).subs(theta, odds/(1+odds))) == 1,
            'Finite-site detailed balance differs from the derivative chemical potential')
    gasrow = next(row for row in rows if not row['loaded'])
    model = make_model(source, gasrow, s_effective=1.)
    sheet = SheetThermo(model.T, .12, .08, 'declared coupled-sheet mathematical fixture')
    gas = gas_equilibrated(model, sheet=sheet, capture_s0=.3)
    reservoir = .5*model.kT*np.log(model.p_gas)
    oracle = minimize_scalar(lambda value: sheet.free_energy(value, model.kT)-reservoir*value,
        bounds=(1e-9, 1-1e-9), method='bounded', options={'xatol': 1e-13})
    require(abs(oracle.x-gas['theta_ML']) < 1e-7, 'Independent grand-potential minimum disagrees')
    wrong = brentq(lambda value: model.kT*np.log(value/(1-value))
                   -sheet.mean_binding(value)-reservoir, 1e-9, 1-1e-9)
    require(abs(wrong-gas['theta_ML']) > .05, 'Mean-versus-differential binding sensitivity absent')
    differences = []
    for value in (.23, .71):
        step = 1e-6
        fd = (sheet.free_energy(value+step, model.kT)-sheet.free_energy(value-step, model.kT))/(2*step)
        differences.append(abs(fd-sheet.chemical_potential(value, model.kT)))
    require(max(differences) < 1e-9, 'Actual free-energy derivative check failed')
    local = local_boundary(model, .37, sheet=sheet, capture_s0=.3)
    actual = make_model(source, gasrow, s_effective=local['s_effective'])
    flux = float(actual.front_flux(local['mu_relative_half_gas1pa_eV'], 0.)[0])
    relative = abs(flux-local['D_atom_net_outflux_m2_s'])/max(abs(flux), 1.)
    require(relative < 1e-12, 'Original Model atom boundary differs from finite-site flux')
    # A separately evaluated symbolic excess energy supplies the synthetic means.
    means = [(value, float((-(-.2*theta+.08*theta**2)/theta).subs(theta, value)))
             for value in (.25, .5)]
    fitted, fit = fit_anchored(model.T, .12, means, source_id='symbolic fixture means')
    require(abs(fitted.interaction_eV-.08) < 1e-14, 'Anchored mean-energy fit failed')
    return dict(symbolic_derivative=str(mu), detailed_balance_identity=True,
        gas_theta=gas['theta_ML'], independent_minimum_theta=oracle.x,
        wrong_mean_binding_theta=wrong, derivative_max_error_eV=max(differences),
        actual_original_boundary_relative_error=relative, fit=fit)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--transport-source', type=Path, required=True)
    parser.add_argument('--operating-points', type=Path, required=True)
    parser.add_argument('--profiles', type=Path, required=True)
    parser.add_argument('--profiles-sha256', required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpu', type=int, required=True)
    args = parser.parse_args()
    if os.sched_getaffinity(0) != {args.cpu}:
        raise ValueError('Declared singleton proof CPU required')
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.process_time()
    def stop(*unused):
        raise TimeoutError('Frozen-sheet proof60CPU-s/75wall-s cap reached')
    signal.signal(signal.SIGXCPU, stop)
    signal.signal(signal.SIGALRM, stop)
    resource.setrlimit(resource.RLIMIT_CPU, (60, 61))
    signal.alarm(75)
    report = dict(schema='coherence-derivative-sheet-proof-v1', Native_launched=False,
                  scope='Mathematical sheet parameters and actual saved-v/source Model, not PdD binding prediction',
                  failures=[])
    try:
        source = load_source(args.transport_source)
        rows = operating_points(args.operating_points)
        loaded = next(row for row in rows if row['loaded'])
        vacancies = read_run_end(args.profiles, args.profiles_sha256, condition='PdD',
            energy_cm_keV=loaded['E_cm_keV'], temperature_K=loaded['T_K'], source_revision=args.revision)
        report['algebra'] = algebra_and_sensitivity(source, rows)
        fine = dict(h0=.125e-9, g1=np.sqrt(1.05), hmid=1.5e-9, g2=np.sqrt(1.13))
        model = make_model(source, loaded, s_effective=1., grid=fine)
        mapped = vacancies.on_grid(model.e)
        original_inventory = float(vacancies.fraction_per_Pd @ np.diff(vacancies.edges_m))
        mapped_inventory = float(mapped @ model.h)
        error = abs(mapped_inventory-original_inventory)/original_inventory
        require(error < 1e-12, 'Frozen-v remap changed vacancy inventory')
        try:
            vacancies.on_grid(model.e*.5)
        except ValueError:
            wrong_foil_refused = True
        else:
            raise AssertionError('Wrong-depth vacancy profile admitted')
        report['vacancies'] = dict(provenance=vacancies.provenance,
            original_integral_v_m=original_inventory, remapped_integral_v_m=mapped_inventory,
            relative_error=error, wrong_foil_refused=wrong_foil_refused)
        factory = lambda T: SheetThermo(T, .12, .08, 'mathematical fixture only')
        result = consume_cases(source=source, rows=rows, thermodynamics=factory,
            capture_probability=lambda T: .3,
            vacancies=vacancies, provenance={'mathematical_fixture': True})
        report['fixture'] = result
        bulk = result['bulk_fed_frozen_v']
        require(bulk['status'] == 'converged_frozen_v_steady', 'Frozen-v fixture failed residual/grid convergence')
        require(bulk['refined']['max_near_surface_x_plus_xt'] < .95,
                'Recorded1percent vacancy field unexpectedly exceeds capacity')
        require(not bulk['refined']['evidence']['vacancy_equation_enforced'],
                'Frozen profile mislabeled as thermal-vacancy steady')
        report['status'] = 'passed'
    except BaseException as exc:
        report['status'] = 'failed'
        report['failures'].append(dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc()))
    finally:
        signal.alarm(0)
        report['process_cpu_s'] = time.process_time()
        report['post_import_cpu_s'] = time.process_time()-started
        report['code'] = [binding(Path(__file__).with_name(name)) for name in
            ('sheet_proof.py', 'sheet_thermo.py', 'coverage_vacancies.py', 'coverage_boundary.py',
             'coverage_solve.py', 'coverage_steady.py', 'coverage_consume.py', 'coverage_inputs.py')]
        (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        print('SHEET_PROOF', report['status'], args.output, report['process_cpu_s'], flush=True)
        if report['failures']:
            print(json.dumps(report['failures'], indent=2), flush=True)
    if report['status'] != 'passed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
