"""Export six saved BO samples after source-bound point admission, not cell convergence.

No Native calculations or preparation/rate inference are performed. The output
must be a new directory; B is copied byte-for-byte from the admitted reuse.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from readout_defect import RADII, binding, close, load_run, BOHR_A, COULOMB_EV_A
from reuse_baseline import SOURCE_SHA
from capsule_io import require


SCHEMA = 'coherencev1-target-conditioned-pair-v1'
OBSERVABLE = 'target_conditioned_BO_total_interaction'


def read_json(path):
    data = json.loads(Path(path).read_bytes())
    # Reject nonfinite values anywhere, including otherwise unused metadata.
    json.dumps(data, allow_nan=False)
    return data


def digest(path):
    data = Path(path).read_bytes()
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def check_binding(value, label):
    require(isinstance(value, dict), label + ': missing source binding')
    require(isinstance(value.get('path'), str) and value['path'] and
            Path(value['path']).name == value['path'] and value['path'] not in ('.', '..'),
            label + ': unsafe source name')
    require(type(value.get('bytes')) is int and value['bytes'] > 0,
            label + ': invalid byte count')
    sha = value.get('sha256')
    require(isinstance(sha, str) and len(sha) == 64 and
            all(c in '0123456789abcdef' for c in sha), label + ': invalid SHA256')


def number(value, label):
    require(type(value) in (int, float) and math.isfinite(value), label + ': nonfinite/non-numeric')
    return value


def check_rows(actual, expected, label):
    require(isinstance(actual, list) and len(actual) == len(RADII), label + ': require six rows')
    require([row['R_A'] for row in actual] == RADII, label + ': wrong ordered support')
    for row, reference in zip(actual, expected):
        for key, value in reference.items():
            close(number(row[key], label + ' ' + key), value, label + ' ' + key)


def check_pair(data, label):
    require(data['schema'] == SCHEMA and data['observable'] == OBSERVABLE,
            label + ': wrong mapped BO observable')
    require('prepared_state' not in data, label + ': prepared-state export is forbidden')
    require(data['R_A'] == RADII and len(data['V_pair_total_eV']) == len(RADII),
            label + ': require exact six-point support')
    for radius, value in zip(data['R_A'], data['V_pair_total_eV']):
        number(radius, label + ' radius')
        number(value, label + ' total interaction')
    reference = data['reference']
    for key in ('source', 'configuration', 'energy_zero', 'periodic_to_isolated'):
        require(isinstance(reference.get(key), str) and reference[key].strip(),
                label + ': missing explicit reference ' + key)
    require(reference['bare'] == 'isolated_unit_charge_coulomb', label + ': wrong bare reference')


def native_case(directory, case):
    """Replay accepted Green mapping and bind the actual saved Native records."""
    directory = Path(directory)
    result = load_run(directory, case, RADII)
    run, quartet = result['run'], result['quartet']
    require(run['completed'] is True and run['schema'] == 'pd-defect-native-run-v1' and
            run['stage'] == 'production' and quartet['schema'] == 'pd-defect-neutral-quartet-v1',
            case + ': wrong completed Native contract')
    require(quartet['specification']['case'] == case, case + ': mismatched specification')
    check_rows(sorted(run['rows'], key=lambda row: row['R_A']), result['rows'], case + ' run rows')
    relaxation = run['relaxation']
    close(quartet['relaxation_raw_fmax_eV_A'], relaxation['raw_fmax_eV_A'],
          case + ': actual geometry/quartet fmax')
    require(relaxation['passed'] is True and
            0 <= number(relaxation['raw_fmax_eV_A'], case + ' actual fmax') <= .1,
            case + ': unqualified actual geometry')
    geometry = directory / ('geometry-selection.json' if 'source_frame' in relaxation
                            else 'relax/relaxation.json')
    require(read_json(geometry) == relaxation, case + ': changed saved geometry record')
    state_bindings = {}
    for role in ('occupied', 'vacant'):
        host_name = role + '-host'
        host_path = directory / 'bo' / (host_name + '.json')
        require(read_json(host_path) == quartet['states'][host_name],
                case + ': changed actual ' + host_name)
        state_bindings[host_name] = binding(host_path)
        for state in quartet['states'][role + '-probes']:
            stem = role + '-probe-R' + format(state['R_A'], '.3f').replace('.', 'p')
            path = directory / 'bo' / (stem + '.json')
            # Radius/projectile annotations are added after save_state writes JSON.
            saved = {key: value for key, value in state.items() if key not in ('R_A', 'projectile_A')}
            require(read_json(path) == saved, case + ': changed actual probe ' + stem)
            state_bindings[stem] = binding(path)
    for key in ('input', 'frozen_source'):
        check_binding(run[key], case + ' ' + key)
    require(binding(directory / 'run-result.json') == result['run_binding'] and
            binding(directory / 'bo/quartet.json') == result['quartet_binding'],
            case + ': source changed during replay')
    result['geometry_binding'] = binding(geometry)
    result['state_bindings'] = state_bindings
    return result


def qualification(readout, cases, baseline, baseline_path):
    require(readout['schema'] == 'pd-defect-screening-readout-v2' and
            readout['point_qualified'] is True and readout['cell_qualified'] is False,
            'Require a point-qualified v2 readout with explicitly missing cell-size qualification')
    require('qualified' not in readout, 'Obsolete undifferentiated qualification is forbidden')
    require(isinstance(readout['qualification_scope'], str) and readout['qualification_scope'].strip(),
            'Require explicit point-only qualification scope')
    require(set(readout['sources']) == {'A', 'C', 'D'} and
            set(readout['rows']) == {'A', 'B', 'C', 'D'}, 'Wrong readout case inventory')
    for case, result in cases.items():
        expected = dict(run=result['run_binding'], quartet=result['quartet_binding'])
        require(readout['sources'][case] == expected, case + ': mismatched readout Native epoch')
        check_rows(readout['rows'][case], result['rows'], case + ' readout rows')
    require(readout['B_source'] == binding(baseline_path), 'Mismatched exact B reuse epoch')
    b_rows = [dict(R_A=r, V_pair_total_eV=v, lowering_eV=COULOMB_EV_A/r-v)
              for r, v in zip(RADII, baseline['V_pair_total_eV'])]
    check_rows(readout['rows']['B'], b_rows, 'B readout rows')
    controls = readout['C_controls']
    require(set(controls) == {'h', 'cell'}, 'Require the C grid check and explicit cell disposition')
    control = controls['h']
    require(control['kind'] == 'h' and control['passed'] is True and control['R_A'] == .1,
            'Unqualified C h control')
    require(control['base'] == cases['C']['run_binding'], 'Mismatched C h base epoch')
    check_binding(control['control'], 'C h control')
    require(control['control']['sha256'] != control['base']['sha256'], 'C h control reuses its base')
    base = number(control['base_lowering_eV'], 'C h control base lowering')
    other = number(control['control_lowering_eV'], 'C h control lowering')
    close(base, cases['C']['rows'][0]['lowering_eV'], 'C h control actual base lowering')
    difference, tolerance = abs(other-base), max(5., .1*abs(other))
    close(number(control['absolute_difference_eV'], 'C h control difference'), difference,
          'C h control difference')
    close(number(control['tolerance_eV'], 'C h control tolerance'), tolerance, 'C h control tolerance')
    require(difference <= tolerance, 'Failed actual C h comparison')
    cell = controls['cell']
    require(cell['converged'] is False and 'passed' not in cell,
            'Cell convergence is waived, not passed')
    require(isinstance(cell['scope'], str) and cell['scope'].strip(), 'Missing cell disposition scope')
    if cell['status'] == 'omitted':
        require('source' not in cell, 'Omitted cell result cannot claim a sensitivity source')
    else:
        require(cell['status'] == 'unrelaxed_sensitivity' and cell['R_A'] == .1,
                'Only omitted or unrelaxed-sensitivity cell dispositions are authorized')
        check_binding(cell['source'], 'C unrelaxed cell sensitivity')
        require(cell['source']['path'] == 'cell-sensitivity.json', 'Wrong sensitivity source record')
        base = number(cell['base_lowering_eV'], 'C sensitivity base lowering')
        other = number(cell['sensitivity_lowering_eV'], 'C sensitivity lowering')
        close(base, cases['C']['rows'][0]['lowering_eV'], 'C sensitivity actual base lowering')
        difference, tolerance = abs(other-base), max(5., .1*abs(other))
        close(number(cell['absolute_difference_eV'], 'C sensitivity difference'), difference,
              'C sensitivity difference')
        close(number(cell['original_tolerance_eV'], 'C sensitivity original tolerance'), tolerance,
              'C sensitivity original tolerance')
        require(type(cell['within_original_tolerance']) is bool and
                cell['within_original_tolerance'] == (difference <= tolerance),
                'Changed sensitivity comparison metadata')
        # Within-tolerance metadata is not a cell-convergence admission gate.
    at_point = {case: rows[0]['lowering_eV'] for case, rows in readout['rows'].items()}
    triggers = [dict(defect=x, baseline=y) for x in ('C', 'D') for y in ('A', 'B')
                if at_point[x] > 0 and (at_point[y] <= 0 or at_point[x] > 3*at_point[y])]
    require(readout['sample_three_times_conditions'] == triggers and
            type(readout['E_required']) is bool and readout['E_required'] == bool(triggers),
            'Readout E trigger does not match the actual six-point curves')


def material_pair(case, result, readout, readout_digest):
    run, quartet = result['run'], result['quartet']
    relaxation, spec = run['relaxation'], quartet['specification']
    host = quartet['states']['occupied-host']
    cell = np.asarray(host['cell_A'], dtype=float)
    volume = float(abs(np.linalg.det(cell)))
    require(math.isfinite(volume) and volume > 0, case + ': invalid actual volume')
    geometry_provenance = {
        key: relaxation[key] for key in (
            'raw_fmax_eV_A', 'force_threshold_eV_A', 'authority', 'target_index',
            'symbols', 'positions_A', 'cell_A', 'initial', 'vacancy', 'frames',
            'source_frame', 'source_last_completed_frame', 'source_last_completed_fmax_eV_A',
            'free_energy_drift_after_selected_eV', 'scope') if key in relaxation
    }
    rows = result['rows']
    curve = dict(
        schema=SCHEMA, observable=OBSERVABLE,
        R_A=[row['R_A'] for row in rows],
        V_pair_total_eV=[row['V_pair_total_eV'] for row in rows],
        reference=dict(
            bare='isolated_unit_charge_coulomb',
            source=('P255/issue247 actual completed Native ' + case +
                    ' neutral-host quartet; SHA256 ' + result['quartet_binding']['sha256']),
            configuration=(
                'Actual frozen standard-H PAW geometry followed by all four PBE/Pd10/Hae BO references; '
                'actual cell, volume, target ray, grid, k sampling and force provenance below. '
                'Neutral occupied/vacant hosts, target-only neutral removal, charge+1 probes with no added electron.'),
            energy_zero=(
                'E_probe,occupied - E_probe,vacant - E_host,occupied + E_host,vacant on the same frozen '
                'background, using E_internal = E_free - minus_TS at Fermi-Dirac width 0.1 eV. '
                'No fitted offset or far-point subtraction; no standard-H PAW energy enters the Hae quartet.'),
            periodic_to_isolated=(
                'Native Hae alpha1=10 Bohr^-2 nuclear compensation replaced by the same-zero-G periodic '
                'point-minus-Gaussian Green block, width sqrt(5)/Bohr, independently replayed from the actual cell/ray. '
                'Representative C h0.20->0.16 comparison passes at R0.1 under max(5 eV,10% of control lowering). '
                'Converged cell-size qualification is explicitly missing and waived. Any optional 27-site '
                'unrelaxed sensitivity is metadata only, never a passed convergence check or an exact isolated limit. '
                'No bare-only image subtraction.')),
        source_bindings=dict(
            run=result['run_binding'], quartet=result['quartet_binding'],
            geometry=result['geometry_binding'], states=result['state_bindings'],
            native_input=run['input'], frozen_executable_source=run['frozen_source'],
            point_qualified_readout=readout_digest,
            cell_sensitivity=readout['C_controls']['cell'].get('source')),
        configuration=dict(case=case, specification=spec, path=quartet['path'],
                           volume_A3=volume, geometry_H_setup='standard_PAW',
                           BO_setups=host['setups'], actual_kpoints=host['actual_kpoints'],
                           realspace_shape=host['realspace_shape'], nbands=host['nbands']),
        geometry_provenance=geometry_provenance,
        mapping_provenance=dict(
            energy_combination=quartet['energy_zero'], rows=rows,
            Hae_alpha1_Bohr_minus2=10, gaussian_width_inverse_A=math.sqrt(5.) / BOHR_A,
            zero_G='same zero-G periodic convention for point and Gaussian'),
        qualification=dict(point_qualified=True, cell_qualified=False, C_controls=readout['C_controls'],
                           sources=readout['sources'], E_required=readout['E_required'],
                           sample_three_times_conditions=readout['sample_three_times_conditions'],
                           scope=readout['qualification_scope']),
        support=('Only the six supplied 0.1-1.6 A samples. Static BO on the common inherited expanded lattice; '
                 'no U(0), prepared state, rate, thermal/ZPE/dynamics, saturation, global-ground-state '
                 'or infinite-cell conclusion. Missing converged cell-size qualification remains an omission. '
                 'Actual vacancy counts and target displacements are geometry descriptors, not stable VH_n assignments. '
                 'Force thresholds and energy drift are provenance, not a certified force-to-screening error bound.'))
    check_pair(curve, case)
    return curve


def export_curves(readout_path, baseline_path, a_path, c_path, d_path, e_path=None):
    """Return filename->JSON curves from saved Native inputs, without writes."""
    readout_path, baseline_path = Path(readout_path), Path(baseline_path)
    readout_digest = digest(readout_path)
    baseline_digest = digest(baseline_path)
    readout, baseline = read_json(readout_path), read_json(baseline_path)
    check_pair(baseline, 'B')
    reuse = baseline['g7_reuse']
    require(reuse['original_sha256'] == SOURCE_SHA and type(reuse['native_runs']) is int and
            reuse['native_runs'] == 0 and reuse['selected_original_indices'] == [0, 1, 2, 4, 5, 6] and
            reuse['excluded_original_radii_A'] == [.4], 'B must retain exact admitted zero-Native reuse')
    cases = {case: native_case(path, case) for case, path in (('A', a_path), ('C', c_path), ('D', d_path))}
    qualification(readout, cases, baseline, baseline_path)
    if readout['E_required']:
        require(e_path is not None, 'Triggered E requires actual complete E Native output')
        cases['E'] = native_case(e_path, 'E')
    else:
        require(e_path is None, 'E is not triggered; do not introduce an unrequested E state')
    curves = {case + '_BO_pair.json': material_pair(case, result, readout, readout_digest)
              for case, result in cases.items()}
    curves['B_BO_pair.json'] = baseline
    require(digest(readout_path) == readout_digest and digest(baseline_path) == baseline_digest,
            'Readout or B source changed during export')
    return curves


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('readout', 'baseline', 'a', 'c', 'd', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--e', type=Path)
    args = parser.parse_args()
    require(not args.output.exists(), 'Output directory must not already exist (no stale curves)')
    baseline_raw = args.baseline.read_bytes()
    curves = export_curves(args.readout, args.baseline, args.a, args.c, args.d, args.e)
    require(args.baseline.read_bytes() == baseline_raw and
            json.loads(baseline_raw) == curves['B_BO_pair.json'],
            'B changed before exact byte copy')
    # Validate everything before publishing a finite, complete directory.
    with tempfile.TemporaryDirectory(prefix='.defect-pairs-', dir=args.output.parent) as temporary:
        staging = Path(temporary) / 'curves'
        staging.mkdir()
        for name, curve in sorted(curves.items()):
            data = baseline_raw if name == 'B_BO_pair.json' else (
                json.dumps(curve, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
            (staging / name).write_bytes(data)
        require(not args.output.exists(), 'Output directory appeared during export')
        staging.rename(args.output)
    print(json.dumps(dict(output=str(args.output), files=sorted(curves),
                         point_qualified=True, cell_qualified=False), indent=2))


if __name__ == '__main__':
    main()
