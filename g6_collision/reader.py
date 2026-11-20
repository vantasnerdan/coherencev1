#!/usr/bin/env python3
"""NumPy/stdlib-only collision capsule CLI: python reader.py CAPSULE_DIRECTORY."""
import argparse
import json
import math
import sys
import numpy as np
from capsule_io import Capsule, require
from controls import controls
from green import BOHR_A
from replay import close, compare_rows, load_case, pair


def replay(root):
    capsule = Capsule(root)
    bare_name = 'receipts/collision-native-nuclear-kernel-002.json'
    bare = capsule.json(bare_name)
    require(bare['passed'] and bare['native_alpha1_Bohr_inverse2'] == [10., 10.], 'Wrong qualified nuclear operator')
    alpha, beta = bare['native_alpha1_Bohr_inverse2']
    width = math.sqrt(alpha * beta / (alpha + beta)) / BOHR_A
    close(width, bare['pair_width_inverse_A'], 'Native compensation width')
    with np.load(capsule.root / 'receipts/collision-native-nuclear-kernel-002.npz', allow_pickle=False) as data:
        cell = data['cell_A'].copy()
    cases = {phase + '-' + role: load_case(capsule, phase, role) for phase in ('base', 'gamma') for role in ('occ', 'vac')}
    pairs = {}; differences = {}; export_differences = {}
    for phase, export in (('base', 'baseline'), ('gamma', 'gamma')):
        occupied = cases[phase + '-occ']; vacant = cases[phase + '-vac']
        close(cell, occupied['host']['cell_A'], 'Qualified operator cell')
        rows = pair(capsule, occupied, vacant, width)
        reference_name = 'receipts/collision-' + phase + '-periodic-pair-readout-001.json'
        reference = capsule.json(reference_name)
        curve = capsule.json('curves/' + export + '_BO_pair.json')
        require(curve['schema'] == 'coherencev1-target-conditioned-pair-v1' and curve['observable'] == 'target_conditioned_BO_total_interaction', 'Wrong curve observable')
        for role, case in (('occupied', occupied), ('vacant', vacant)):
            for kind, name in case['bindings'].items():
                capsule.bind(reference['source_bindings'][role][kind], name)
                capsule.bind(curve['source_bindings'][role][kind], name)
        capsule.bind(curve['source_bindings']['bare'], bare_name)
        for key, name in (('h', 'collision-h016-static-comparison-001.json'), ('cell', 'collision-image-comparison-001.json'), ('actual_displacement', 'collision-displacement-pair-admission-001.json')):
            capsule.bind(curve['qualification'][key], 'receipts/' + name)
        differences[phase] = compare_rows(rows, reference['rows'], phase)
        ordered = sorted(rows, key=lambda row: row['R_A'])
        export_differences[phase] = dict(R_A=close([r['R_A'] for r in ordered], curve['R_A'], 'Export support'),
            V_pair_total_eV=close([r['V_pair_periodic_point_eV'] for r in ordered], curve['V_pair_total_eV'], 'Exact qualified BO export'))
        pairs[phase] = rows
    geometry = {}
    for role in ('occ', 'vac'):
        base = cases['base-' + role]['ground']; gamma = cases['gamma-' + role]['ground']
        require(base['symbols'] == gamma['symbols'], 'Displacement species changed')
        expected = np.zeros_like(base['positions_A'], dtype=float)
        for i, symbol in enumerate(base['symbols']):
            if symbol == 'H':
                expected[i] = [.0768, 0, 0]
        geometry[role] = close(np.array(gamma['positions_A']) - base['positions_A'], expected, 'All-D shift with fixed Pd')
        close(np.array(gamma['path']['target_A']) - base['path']['target_A'], [.0768, 0, 0], 'Reanchored target')
    contrasts = []
    for base, gamma in zip(pairs['base'], pairs['gamma']):
        contrasts.append(dict(R_A=base['R_A'], base_candidate_lowering_eV=base['candidate_lowering_vs_isolated_C_over_R_eV'],
            gamma_candidate_lowering_eV=gamma['candidate_lowering_vs_isolated_C_over_R_eV'],
            delta_gamma_minus_base_lowering_eV=gamma['candidate_lowering_vs_isolated_C_over_R_eV']-base['candidate_lowering_vs_isolated_C_over_R_eV'],
            delta_gamma_minus_base_encounter_excitation_eV=gamma['encounter_excitation_change_eV']-base['encounter_excitation_change_eV'],
            base_vacancy_mean_loss_eV_A=base['vacancy_mean_loss_eV_A'], gamma_vacancy_mean_loss_eV_A=gamma['vacancy_mean_loss_eV_A']))
    displacement = capsule.json('receipts/collision-displacement-pair-admission-001.json')
    require(displacement['passed'], 'Displacement gate failed')
    contrast_diff = compare_rows(contrasts, displacement['contrasts'], 'All-seven displacement contrast')
    gates = controls(capsule, cases, pairs['base'], width)
    pinned = capsule.json('receipts/collision-qualified-consumer-proof-001.json')
    pinned_differences = {}
    for phase, export in (('base', 'baseline'), ('gamma', 'gamma')):
        capsule.bind(pinned['phases'][export]['source'], 'curves/' + export + '_BO_pair.json')
        pinned_differences[phase] = close(pairs[phase][-1]['candidate_lowering_vs_isolated_C_over_R_eV'], pinned['phases'][export]['minimum_radius_lowering_eV'], 'Pinned consumer Coulomb constant difference')
    return dict(schema='pdd-numpy-saved-output-replay-v1', passed=True, pairs=pairs, contrasts=contrasts, controls=gates,
                replay_differences=dict(pair_columns=differences, exact_BO_exports=export_differences, contrasts=contrast_diff,
                    trajectory={key: case['differences'] for key, case in cases.items()}, geometry_A=geometry,
                    pinned_consumer_different_C_eV=pinned_differences),
                full_trajectory_rows={key: len(case['td']['rows']) for key, case in cases.items()},
                force_work_energy_residuals_eV={key: case['max_force_work_energy_residual_eV'] for key, case in cases.items()},
                qualifications=dict(encounter_excitation='production only; NOT dt-checked', vacancy_drag='finite-path production; halfdt pending, may be second asset',
                    h_cell='representative R0.1 point admits both supplied BO curves; not exact infinite-cell limit',
                    pinned_consumer='2645b19 uses a different Coulomb constant; difference <1e-9eV is documented, not tuned'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('capsule')
    args = parser.parse_args()
    try:
        print(json.dumps(replay(args.capsule), indent=2, sort_keys=True, allow_nan=False))
    except (ValueError, KeyError, OSError) as error:
        print('REFUSED: ' + str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
