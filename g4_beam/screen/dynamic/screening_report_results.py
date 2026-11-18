"""Combine actual static profiles with separately labelled fit benchmarks."""
import argparse
import csv
import json
from pathlib import Path

from g4_beam.screen.dynamic.largeq_bank import binding, require
from g4_beam.screen.dynamic.static_comparison import compare_static


def benchmark_rows(reference):
    for row in reference['pass4_uniform_bulk']:
        yield dict(kind='pass4_uniform_bulk', cause='uniform', reading=row['measurement'],
                   loading_eV=row['loading_PdD_minus_Pd_eV'],
                   fitted_Pd_eV=row['materials']['Pd']['Ue_common_proton_eV'],
                   fitted_PdD_eV=row['materials']['PdD']['Ue_common_proton_eV'])
    for row in reference['transported_refits']:
        yield dict(kind='pass5_cause_parameter', cause=row['cause'],
                   reading='refit power' + str(row['inverse_current_power']),
                   loading_eV=row['bulk_loading_PdD_minus_Pd_eV'],
                   fitted_Pd_eV=row['energy_parameters_eV']['Ue_Pd'],
                   fitted_PdD_eV=row['energy_parameters_eV']['Ue_PdD'],
                   selected_for_Target=row['selected_by_actual_api_min_chi2'])
    for row in reference['Eq7_curve_diagnostics']:
        yield dict(kind='direct_Eq7_curve_diagnostic', cause=row['cause'],
                   reading='additional power' + str(row['additional_inverse_current_power']),
                   loading_eV=row['paper_curve_loading_PdD_minus_Pd_eV'],
                   fitted_Pd_eV=row['paper_curve_Eq7_effective_eV']['Pd'],
                   fitted_PdD_eV=row['paper_curve_Eq7_effective_eV']['PdD'])


def build_report(static, reference, speed):
    comparison = compare_static(static)
    loading = comparison['loading18']['N3_eV'][0]
    rows = []
    for benchmark in benchmark_rows(reference):
        require(benchmark['loading_eV'] > 0, 'A signed loading benchmark requires its positive source step')
        rows.append(dict(**benchmark, computed_static_loading_eV=loading,
                         signed_fraction_of_fitted_step=loading / benchmark['loading_eV'],
                         computed_Pd_U18_eV=comparison['U18']['Pd']['U18_N3_eV'][0],
                         computed_PdD_U18_eV=comparison['U18']['PdD']['U18_N3_eV'][0]))
    proxy_definitions = {}
    for material in ('Pd', 'PdD'):
        old = speed['models'][material]['U0_eV']
        full10 = comparison['matched_largeq_U10'][material]['U10_largeq_N3_eV'][0]
        full18 = comparison['U18'][material]['U18_N3_eV'][0]
        proxy_definitions[material] = dict(historical_lowq_U10_proxy_eV=old,
            largeq_U10_PW500_eV=full10, largeq_U18_PW1200_eV=full18,
            full10_minus_old_proxy_eV=full10 - old,
            full18_minus_old_proxy_NOT_core_eV=full18 - old,
            interpretation='Historical low-q C/lambda proxy is not a matched microscopic shell baseline; '
                           'use the separately reported PW50018-minus10 activation estimate.')
    return dict(schema='screening-static-speed-fit-comparison-v1', static=comparison,
        definition_comparison=proxy_definitions, benchmark_rows=rows,
        step2_four_laboratory_energies=speed['four_laboratory_energies'],
        step2_velocity_models=speed['models'], step2_velocity_controls=speed['velocity_controls'],
        dynamic_scope='Only the established ten-valence low-q Huke closure is used for step2. '
                      'No automatic scalar speed factor is applied to the localized4s/4p response.',
        benchmark_scope='Pass4 uniform bulk,pass5 transported-cause parameters and directEq7 '
                        'curve-fit diagnostics are different phenomenological quantities. '
                        'All are shown, none is added to computed U or used as a consumer override.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--static', required=True, type=Path)
    parser.add_argument('--fit-reference', required=True, type=Path)
    parser.add_argument('--speed-result', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--csv', required=True, type=Path)
    args = parser.parse_args()
    result = build_report(json.loads(args.static.read_text()),
                          json.loads(args.fit_reference.read_text()),
                          json.loads(args.speed_result.read_text()))
    result['source_bindings'] = [binding(path) for path in
        (args.static, args.fit_reference, args.speed_result, Path(__file__))]
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    columns = ['kind', 'cause', 'reading', 'selected_for_Target', 'loading_eV',
               'computed_static_loading_eV', 'signed_fraction_of_fitted_step',
               'fitted_Pd_eV', 'fitted_PdD_eV', 'computed_Pd_U18_eV', 'computed_PdD_U18_eV']
    with args.csv.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(result['benchmark_rows'])
    print(json.dumps(dict(output=binding(args.output), csv=binding(args.csv),
        loading18_contact_eV=result['static']['loading18']['N3_eV'][0],
        sampled_loading_numerics=result['static']['loading18']['empirical_numerics']['all_sampled_changes_within_tolerance'],
        step3=result['static']['step3_trigger'])))


if __name__ == '__main__':
    main()
