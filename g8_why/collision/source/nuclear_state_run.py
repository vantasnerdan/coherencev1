"""Source-bound finite-Hermite nuclear-state calculation and saved-output replay."""
import argparse
import json
import os
from pathlib import Path
import resource
import time

import numpy as np

import nuclear_momentum
import quantum_dvr
from well_field import file_binding
from well_model import CalibratedWell


def run(model_path, n, output):
    started = time.process_time()
    output.mkdir(parents=True, exist_ok=False)
    well = CalibratedWell(model_path)
    source_paths = [Path(__file__), Path(quantum_dvr.__file__), Path(nuclear_momentum.__file__),
                    Path(__file__).with_name('well_model.py'), Path(__file__).with_name('well_field.py'),
                    Path('<PRIVATE_PATH>')]
    inputs = [file_binding(p) for p in source_paths] + [well.binding, well.field.binding]
    (output / 'input-bindings.json').write_text(json.dumps(inputs, indent=2) + '\n')
    print(f'NUCLEAR_STATE_START n={n} cpu={time.process_time():.6f}', flush=True)
    spectrum = quantum_dvr.solve(well, n)
    print(f'NUCLEAR_EIGENSPECTRUM_COMPLETE n={n} cpu={time.process_time():.6f}', flush=True)
    weights, log_partition = quantum_dvr.canonical_weights(spectrum.energies, 300.0)
    quantum_dvr.save_spectrum(output / 'nuclear-spectrum.npz', spectrum, weights, log_partition, 300.0)
    operator_check = quantum_dvr.consume_spectrum(output / 'nuclear-spectrum.npz')
    momentum_check = nuclear_momentum.momentum_and_position(spectrum, weights, output)
    readout = nuclear_momentum.consume_momentum(output, well)
    invariant_ok = (operator_check['max_eigenresidual_eV'] <= 1e-8
                    and operator_check['max_vector_norm_error'] <= 1e-10
                    and abs(operator_check['weight_sum'] - 1) <= 1e-12
                    and momentum_check['max_state_normalization_error'] <= 1e-10
                    and momentum_check['kinetic_first_identity_error_eV'] <= 1e-10
                    and momentum_check['kinetic_second_identity_error_eV2'] <= 1e-9
                    and abs(readout['joint_normalization'] - 1) <= 1e-12)
    localization_ok = readout['localization']['probability_r_gt_0.75A'] < 1e-3
    report = dict(schema='g8-nuclear-state-result-v1', n_per_axis=n, state_count=n**3,
                  method='Complete finite-Hermite3D DVR; exact projected p^2; n-point potential collocation',
                  basis_length_A=0.30, deuteron_mass_kg=quantum_dvr.DEUTERON_MASS_KG,
                  temperature_K=300, log_partition=log_partition,
                  ground_energy_eV=float(spectrum.energies[0]), maximum_energy_eV=float(spectrum.energies[-1]),
                  minimum_potential_grid_eV=float(spectrum.potential.min()),
                  input_bindings=inputs, operator_consumer=operator_check,
                  momentum_quadrature_identities=momentum_check, saved_output_consumer=readout,
                  scientific_invariants_pass=invariant_ok, localization_pilot_pass=localization_ok,
                  localization_pilot_criterion='Finite-state quadrature probability r>0.75A <1e-3; a localization screen, not a continuum theorem.',
                  process_cpu_seconds=time.process_time() - started,
                  process_resource_usage=dict(user_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime,
                                              system_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_stime,
                                              max_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
                  affinity=sorted(os.sched_getaffinity(0)),
                  findings=[dict(id='finite-spectrum', disposition='established_as_computed' if invariant_ok else 'invariant_failure',
                                 statement='All n^3 states retained in the canonical partition and joint table; no thermal level truncation.'),
                            dict(id='localization', disposition='localized_screen_pass' if localization_ok else 'localization_failure',
                                 statement='Oversampled finite-state spatial probability inspected before expansion.'),
                            dict(id='physical-tail', disposition='unresolved',
                                 statement='Physical point-Coulomb cusp/continuum moment error is uncontrolled. Finite-state first/second-moment envelopes are numerical, not rigorous physical bounds. Physical high-p tail remains [0,1].'),
                            dict(id='ensemble-scope', disposition='limitation',
                                 statement='Clamped one-body canonical model, not exact correlated PdD; DVR potential collocation is not exact Galerkin integration. No hard wavefunction box or Gaussian physical-tail assumption.')])
    report['output_bindings'] = [file_binding(p) for p in sorted(output.iterdir()) if p.is_file()]
    (output / 'result.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(f'NUCLEAR_STATE_COMPLETE n={n} invariants={invariant_ok} localization={localization_ok} cpu={time.process_time():.6f}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', required=True, type=Path)
    parser.add_argument('--n', required=True, type=int, choices=(14, 18))
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    run(args.model, args.n, args.output)
