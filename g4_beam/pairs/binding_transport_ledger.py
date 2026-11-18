"""Per-block electronic Eb equivalents in the existing finite-T trap convention.

This is not a multi-occupancy grand partition or a new transport trajectory.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import time
import numpy as np
from g4_beam.pairs.surface_io import manifest, save

PROFILE_SHA = 'e2e91077595de602b815ae86f45c1a675d6a553153272d0e1295b72e2a5a9943'
TRANSPORT_SHA = '383aad3068c5d5b62fc5f031ccfe8c35c1c6ab426f41daeca6c18aeba3d4a47e'


def consume(binding_path, profiles_path, transport_path, selection_label):
    if (manifest(profiles_path)['sha256'] != PROFILE_SHA or
            manifest(transport_path)['sha256'] != TRANSPORT_SHA):
        raise ValueError('The actual run-end state and retained trap convention are required')
    ledger = json.loads(Path(binding_path).read_text())
    spec = importlib.util.spec_from_file_location('pair_ledger_transport', transport_path)
    transport = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(transport)
    record = json.loads(Path(profiles_path).read_text())['records'][6]
    temperature = record['case']['T_K']
    iso = transport.Isotherm(transport.PD, temperature)
    mobile_x = np.asarray(record['state']['x'])
    raw_v = np.asarray(record['state']['v'])
    # Exactly the original RHS projection of tiny negative numerical vacancies.
    vacancy = np.maximum(raw_v, 0)
    mu = iso.mu(mobile_x)
    mu0 = iso.mu0()
    capacity = transport.PD['mtrap']
    baseline_Eb = transport.PD['Eb']
    baseline_theta = transport.langmuir_x(mu, mu0-baseline_Eb, 1.0, iso.kT)
    rows = []
    for block in ledger['blocks']:
        Eb = block['binding_eV_per_added_D']
        theta = transport.langmuir_x(mu, mu0-Eb, 1.0, iso.kT)
        # Expose the sign/reference/unit map on the actual unsaturated front cell.
        t, t0 = float(theta[0]), float(baseline_theta[0])
        if not 0 < t < 1 or not 0 < t0 < 1:
            raise ValueError('Front-cell detailed-balance check requires an unsaturated state')
        ratio = (t/(1-t))/(t0/(1-t0))
        expected = np.exp((Eb-baseline_Eb)/iso.kT)
        error = abs(ratio/expected-1)
        if error > 1e-10:
            raise ValueError('Binding-to-source-trap detailed balance failed')
        rows.append(dict(block=block, source_Eb_equivalent_eV=Eb,
                         theta_eq_single_block=theta.tolist(),
                         trapped_D_per_metal_single_block=(capacity*vacancy*theta).tolist(),
                         front_odds_ratio_to_default=ratio,
                         detailed_balance_relative_error=error))
    return dict(schema='coherence-pairs-transport-binding-ledger-v1',
                selection_label=selection_label,
                inputs=[manifest(binding_path), manifest(profiles_path), manifest(transport_path)],
                temperature_K=temperature, source_mu0_eV=float(mu0),
                source_delta_mu_eV=(mu-mu0).tolist(), z_m=record['state']['zc'],
                actual_runend_v=raw_v.tolist(), source_projected_v=vacancy.tolist(),
                source_mtrap=capacity, tested_occupancies=ledger['computed_occupancies'],
                physical_capacity_proven=False, current_transport_Eb_eV=baseline_Eb,
                current_theta_eq=baseline_theta.tolist(),
                current_trapped_D_equilibrium_diagnostic=(capacity*vacancy*baseline_theta).tolist(),
                rows=rows, replace_current_transport=False,
                interpretation=[
                    'Positive B transfers D favorably from the matched bulk interstitial into the vacancy',
                    'For each block separately:Eb=B_el in the original theta=logistic((mu-mu0+Eb)/kT)',
                    'This electronic-only mapping assumes omitted differential ionic ZPE/thermal free energy is zero',
                    'mtrap6 is retained phenomenological source capacity, not inferred from maximum computed occupancy',
                    'D2-to-D4 and D4-to-D6 remain blockmeans; no D3/D5 energy or site degeneracy is invented',
                    'Single-block curves are alternative source-Eb equivalents, not a simultaneous multi-occupancy partition',
                    'Actual7200s vacancy profile stays frozen; equilibrium diagnostics do not replace run-end xt or evolve v_eq',
                    'The default0.23eV is not changed by this diagnostic; D1/D2 size qualification and physical mapping remain separate',
                ])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binding', required=True)
    parser.add_argument('--profiles', required=True)
    parser.add_argument('--transport', required=True)
    parser.add_argument('--selection-label', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    started = time.process_time()
    result = consume(args.binding, args.profiles, args.transport, args.selection_label)
    result['process_cpu_s'] = time.process_time()-started
    save(args.output, result)
    print(json.dumps(dict(selection_label=result['selection_label'],
          source_mtrap=result['source_mtrap'], source_mu0_eV=result['source_mu0_eV'],
          front=[dict(block=r['block'], theta=r['theta_eq_single_block'][0],
                      odds_ratio=r['front_odds_ratio_to_default'],
                      check=r['detailed_balance_relative_error']) for r in result['rows']],
          process_cpu_s=result['process_cpu_s'])))
