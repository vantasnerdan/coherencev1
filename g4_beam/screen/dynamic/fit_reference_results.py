"""Extract all three reading refits and separate their Eq7 diagnostics."""
import argparse
import csv
import json
import math
from pathlib import Path

from g4_beam.screen.dynamic.largeq_bank import binding, require


CAUSES = ('screeningrise', 'surface', 'neutrals')
MATERIALS = ('Pd', 'PdD', 'Ti', 'TiD')


def _one_value(rows, name):
    values = {row[name] for row in rows}
    require(len(values) == 1, 'A reported global fit varies across energy rows: ' + name)
    value = values.pop()
    require(isinstance(value, (float, int)) and math.isfinite(value),
            'Nonfinite or missing reported fit: ' + name)
    return value


def _bulk_uniform_rows(path):
    with path.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    measurements = sorted({row['measurement'] for row in rows})
    require(len(measurements) == 6, 'Expected all six pass4 charge-normalization hypotheses')
    result = []
    for measurement in measurements:
        group = [row for row in rows if row['measurement'] == measurement]
        require(len(group) == 4 and {row['material'] for row in group} == set(MATERIALS),
                'Incomplete uniform-target material set')
        values = {row['material']: dict(
            Ue_common_proton_eV=float(row['Ue_common_proton_keV']) * 1000,
            common_sigma_eV=float(row['common_sigma_keV']) * 1000,
            Ue_Eq7_only_eV=float(row['Ue_effective_eq7_keV']) * 1000) for row in group}
        result.append(dict(measurement=measurement, materials=values,
            loading_PdD_minus_Pd_eV=values['PdD']['Ue_common_proton_eV'] -
                                   values['Pd']['Ue_common_proton_eV']))
    return result


def extract_reference(fits_path, beam_path, speed_path, bulk_path):
    fits = json.loads(fits_path.read_text())
    beam = json.loads(beam_path.read_text())
    speed = json.loads(speed_path.read_text())
    bulk = _bulk_uniform_rows(bulk_path)
    actual = binding(fits_path)
    expected = beam['source_bindings']['transported_cause_fits']
    require(actual['bytes'] == expected['bytes'] and actual['sha256'] == expected['sha256'],
            'Cause refits are not the actual source bound by the saved Target comparison')
    keys = [(row['cause'], row['current_power']) for row in fits]
    require(len(keys) == len(set(keys)) and set(keys) ==
            {(cause, power) for cause in CAUSES for power in range(3)},
            'Expected the complete unique cause-by-three-reading fit set')
    best = {cause: min((row for row in fits if row['cause'] == cause),
                       key=lambda row: row['chi2']) for cause in CAUSES}
    parameter_rows = []
    for row in fits:
        params = row['params']
        converted = {key: value * 1000 for key, value in params.items()
                     if key.startswith('Ue_') or key.startswith('Ues_')}
        parameter_rows.append(dict(cause=row['cause'], inverse_current_power=row['current_power'],
            selected_by_actual_api_min_chi2=(row is best[row['cause']]),
            chi2=row['chi2'], dof=row['dof'], source_parameters=params,
            energy_parameters_eV=converted,
            bulk_loading_PdD_minus_Pd_eV=converted['Ue_PdD'] - converted['Ue_Pd'],
            interpretation='Phenomenological transported-cause parameters; not microscopic screening.'))
    diagnostics = []
    for cause in CAUSES:
        rows = [row for row in beam['rows'] if row['cause'] == cause]
        require(rows and all(row['cause_fit_current_power'] == best[cause]['current_power'] and
                abs(row['cause_fit_chi2'] - best[cause]['chi2']) < 1e-10 for row in rows),
                'Saved Target table did not use the declared best transported fit')
        for power in range(3):
            values = {material: _one_value(rows, f'matched_{material}_Ue_effective_{power}_keV') * 1000
                      for material in MATERIALS}
            diagnostics.append(dict(cause=cause, additional_inverse_current_power=power,
                native_cause_fit_current_power=best[cause]['current_power'],
                target_curve_Eq7_effective_eV=_one_value(rows, f'Ue_effective_{power}_keV') * 1000,
                paper_curve_Eq7_effective_eV=values,
                paper_curve_loading_PdD_minus_Pd_eV=values['PdD'] - values['Pd'],
                interpretation='Literal2.5<=Ecm<6.5keV Eq7 fits of plotted/weighted curves. '
                    'Additional inverse-current powers do not refit the cause parameters. '
                    'Never add these effective values to microscopic U.'))
    return dict(schema='all-reading-phenomenological-reference-v1',
        pass4_uniform_bulk=bulk,
        pass4_bulk_source_commit='853c37738caca76ba70cefd8f4df7350bc204937',
        pass4_bulk_scope='Uniform thick-target2.5–6.5keV fits to28 experimental points; '
                         'fixed common reaction term; all six charge hypotheses.',
        pass4_exact_loading_range_eV=[min(row['loading_PdD_minus_Pd_eV'] for row in bulk),
                                     max(row['loading_PdD_minus_Pd_eV'] for row in bulk)],
        transported_refits=parameter_rows, Eq7_curve_diagnostics=diagnostics,
        step2_four_laboratory_energies=speed['four_laboratory_energies'],
        requested_loading_benchmark_eV=speed['four_laboratory_energies'][0]['fitted_loading_step_eV'],
        requested_benchmark_scope='Rounded pass4 uniform-bulk step; not the transported-refit or direct Eq7 diagnostic.',
        static_or_dynamic_microscopic_U_is_charge_reading_independent=True,
        no_refitting_or_new_forward_run=True,
        source_bindings=[binding(path) for path in (fits_path, beam_path, speed_path, bulk_path, Path(__file__))])


def _comparison_rows(result):
    for row in result['pass4_uniform_bulk']:
        yield dict(kind='pass4_uniform_bulk', cause='uniform', reading=row['measurement'],
            U_Pd_eV=row['materials']['Pd']['Ue_common_proton_eV'],
            U_PdD_eV=row['materials']['PdD']['Ue_common_proton_eV'],
            loading_PdD_minus_Pd_eV=row['loading_PdD_minus_Pd_eV'])
    for row in result['transported_refits']:
        yield dict(kind='pass5_cause_parameter', cause=row['cause'],
            reading='refit inverse-current power' + str(row['inverse_current_power']),
            selected_for_Target=row['selected_by_actual_api_min_chi2'],
            U_Pd_eV=row['energy_parameters_eV']['Ue_Pd'],
            U_PdD_eV=row['energy_parameters_eV']['Ue_PdD'],
            loading_PdD_minus_Pd_eV=row['bulk_loading_PdD_minus_Pd_eV'])
    for row in result['Eq7_curve_diagnostics']:
        yield dict(kind='direct_Eq7_curve_diagnostic', cause=row['cause'],
            reading='additional inverse-current power' + str(row['additional_inverse_current_power']),
            U_Pd_eV=row['paper_curve_Eq7_effective_eV']['Pd'],
            U_PdD_eV=row['paper_curve_Eq7_effective_eV']['PdD'],
            loading_PdD_minus_Pd_eV=row['paper_curve_loading_PdD_minus_Pd_eV'],
            target_curve_effective_eV=row['target_curve_Eq7_effective_eV'])


def write_reference_csv(path, result):
    columns = ['kind', 'cause', 'reading', 'selected_for_Target', 'U_Pd_eV',
               'U_PdD_eV', 'loading_PdD_minus_Pd_eV', 'target_curve_effective_eV']
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(_comparison_rows(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cause-fits', required=True, type=Path)
    parser.add_argument('--beam-result', required=True, type=Path)
    parser.add_argument('--speed-result', required=True, type=Path)
    parser.add_argument('--bulk-screening', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--csv', required=True, type=Path)
    args = parser.parse_args()
    result = extract_reference(args.cause_fits, args.beam_result, args.speed_result, args.bulk_screening)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    write_reference_csv(args.csv, result)
    print(json.dumps(dict(output=binding(args.output), csv=binding(args.csv),
        selected_loading_eV={row['cause']: row['bulk_loading_PdD_minus_Pd_eV']
            for row in result['transported_refits'] if row['selected_by_actual_api_min_chi2']},
        refit_rows=len(result['transported_refits']), diagnostics=len(result['Eq7_curve_diagnostics']),
        pass4_bulk_loading_range_eV=result['pass4_exact_loading_range_eV'])))


if __name__ == '__main__':
    main()
