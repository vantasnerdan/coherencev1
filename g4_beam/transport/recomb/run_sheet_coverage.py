"""Source-bound physical sheet fit and no-beam gas/bulk-fed coverage readout."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
from pathlib import Path
import resource
import signal
import time
import traceback
from coverage_inputs import binding, load_source, operating_points
from coverage_consume import consume_cases
from coverage_vacancies import read_run_end
from sheet_material import build_material
from surface_io import save


def stop(*unused):
    raise TimeoutError('Physical sheet readout60CPU-s/75wall-s cap reached')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ledger', type=Path, required=True)
    parser.add_argument('--ledger-sha256', required=True)
    parser.add_argument('--transport-source', type=Path, required=True)
    parser.add_argument('--operating-points', type=Path, required=True)
    parser.add_argument('--profiles', type=Path, required=True)
    parser.add_argument('--profiles-sha256', required=True)
    parser.add_argument('--profile-revision', required=True)
    parser.add_argument('--capture-probability', type=float, required=True)
    parser.add_argument('--elbow-readout', type=Path)
    parser.add_argument('--elbow-readout-sha256')
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if os.sched_getaffinity(0) != {args.cpu}:
        raise ValueError('Declared singleton readout CPU required')
    if not 0 < args.capture_probability <= 1:
        raise ValueError('Explicit intrinsic capture probability in(0,1] required')
    if bool(args.elbow_readout) != bool(args.elbow_readout_sha256):
        raise ValueError('The optional electronic ridge needs both readout and digest')
    args.output.mkdir(parents=True, exist_ok=False)
    resource.setrlimit(resource.RLIMIT_CPU, (60, 61))
    signal.signal(signal.SIGXCPU, stop)
    signal.signal(signal.SIGALRM, stop)
    signal.alarm(75)
    report = dict(schema='coherence-physical-sheet-coverage-v1', Native_launched=False,
        capture_probability_computed=False, assumed_capture_probability=args.capture_probability,
        failures=[])
    try:
        source = load_source(args.transport_source)
        rows = operating_points(args.operating_points)
        capture = lambda temperature: args.capture_probability
        if args.elbow_readout:
            from clean_capture import capture_from_elbow
            reference = binding(args.elbow_readout)
            if reference['sha256'] != args.elbow_readout_sha256:
                raise ValueError('Electronic route readout identity changed')
            capture, report['restricted_elbow_rate'] = capture_from_elbow(reference,
                source=source, rows=rows, assumed_capture=args.capture_probability)
        models, material = build_material(args.ledger, args.ledger_sha256)
        report['material'] = material
        loaded = next(row for row in rows if row['loaded'])
        vacancies = read_run_end(args.profiles, args.profiles_sha256, condition='PdD',
            energy_cm_keV=loaded['E_cm_keV'], temperature_K=loaded['T_K'],
            source_revision=args.profile_revision)
        result = consume_cases(source=source, rows=rows,
            thermodynamics=lambda temperature: models[temperature],
            capture_probability=capture,
            vacancies=vacancies, provenance=material)
        report['coverage'] = result
        if result['bulk_fed_frozen_v']['status'] != 'converged_frozen_v_steady':
            raise ArithmeticError('Physical bulk-fed coverage failed grid/residual qualification')
        comparison = result['overflow_term_on_comparison']
        if comparison is not None and comparison['status'] != 'converged_frozen_v_steady':
            raise ArithmeticError('Required overflow-on comparison remains unconverged')
        report['inputs'] = [binding(path) for path in (args.ledger, args.transport_source,
            args.operating_points, args.profiles)]
        if args.elbow_readout:
            report['inputs'].append(binding(args.elbow_readout))
        report['status'] = 'qualified_with_declared_local_model_and_assumed_capture'
        if args.elbow_readout:
            report.update(status='numerically_converged_restricted_route_sensitivity',
                transport_ledger_replacement_authorized=False,
                transport_source_defaults_changed=False,
                peak_k_slab_convergence_qualified=False)
    except BaseException as exc:
        report['status'] = 'failed'
        report['failures'].append(dict(type=type(exc).__name__, message=str(exc),
                                      traceback=traceback.format_exc()))
    finally:
        signal.alarm(0)
        report['process_cpu_s'] = time.process_time()
        report['code'] = [binding(Path(__file__).with_name(name)) for name in
            ('run_sheet_coverage.py', 'sheet_material.py', 'sheet_thermo.py',
             'coverage_inputs.py', 'coverage_vacancies.py', 'coverage_boundary.py',
             'coverage_steady.py', 'coverage_solve.py', 'coverage_consume.py',
             'thermo_harmonic.py', 'thermo_boundary.py')]
        if args.elbow_readout:
            report['code'].extend(binding(Path(__file__).with_name(name))
                for name in ('clean_capture.py', 'elbow_path_results.py'))
        save(args.output/'report.json', report)
        print('PHYSICAL_SHEET_READOUT', report['status'], args.output,
              report['process_cpu_s'], flush=True)
        for failure in report['failures']:
            print(failure['traceback'], flush=True)
    if report['status'] == 'failed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
