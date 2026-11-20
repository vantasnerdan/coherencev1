"""Replay saved neutral/bare LiD and gas BO quartets; qualify both F controls.

CLI writes a new JSON file only. This is a finite-R method comparison, not a gas
U(0), rate calculation, or proof of an exact isolated electronic limit.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from hydride_replay import (RADII, binding, close, compare_json, load_run,
                            require)

CHARGE_SCOPE = (
    'q0: neutral H/D probe adds one electron; q1: bare unit-charge probe adds none. '
    'Occupied and target-removed hosts are neutral, differing by one nucleus and one electron. '
    'Matched same-cell probe charges cancel a gauge/background constant proportional to q^2 '
    'within each quartet, not between unlike charges or cells. Residual finite-cell electronic '
    'and multipole effects are checked at R=0.1 A for BOTH F charges.')
SCOPE = (
    'Static target-conditioned BO total interactions at the six supplied radii only. '
    'F uses freshly force-relaxed finite FD cells at the original PBE PW500/k4 LiD volume; '
    'this is not FD volume/stress reoptimization. G is a D atom in a periodic 10 A box; '
    'isotope label D does not change the nonspin BO electronic energy. '
    'No fitted offset, U0, drag, rate, prepared state, experimental gas screening value, '
    'thermal/ZPE correction, or exact isolated electronic limit is inferred.')


def _complete(result, label):
    for charge in ('q0', 'q1'):
        require([row['R_A'] for row in result['rows'][charge]] == RADII,
                label + ': full six-R curves required for both charges')


def _controls(base, fine, large):
    sb = base['run']['specification']
    sf = fine['run']['specification']
    sl = large['run']['specification']
    require(sb['case'] == sf['case'] == sl['case'] == 'F', 'Controls must belong to F')
    require(sb['h_A'] == sl['h_A'] == .2 and sf['h_A'] == .16 and
            sb['repeat_primitive'] == sf['repeat_primitive'] == 2 and sl['repeat_primitive'] == 3,
            'Require actual h0.20/0.16 and repeat2/repeat3 F controls')
    for control in (fine, large):
        spec = control['run']['specification']
        for field in ('method', 'primitive_symbols', 'primitive_cell_A', 'primitive_positions_A',
                      'volume_source', 'force_tolerance_eV_A', 'ionic_steps', 'kmesh'):
            compare_json(spec[field], sb[field], 'Matched F control ' + field)
        close(control['run']['geometry']['path']['direction'],
              base['run']['geometry']['path']['direction'], 'Matched F approach direction')
    close(fine['run']['geometry']['cell_A'], base['run']['geometry']['cell_A'], 'Same h-control cell')
    close(large['run']['geometry']['cell_A'], 1.5 * np.asarray(base['run']['geometry']['cell_A']),
          'Repeat3 versus repeat2 control cell')
    bg = np.asarray(base['quartet']['states']['occupied-host']['realspace_shape'])
    fg = np.asarray(fine['quartet']['states']['occupied-host']['realspace_shape'])
    require(bool(np.all(fg >= bg) and np.any(fg > bg)), 'Nominal h control did not refine actual grid')
    controls = {'h': {}, 'cell': {}}
    for kind, other in (('h', fine), ('cell', large)):
        for charge in ('q0', 'q1'):
            a, b = base['rows'][charge][0], other['rows'][charge][0]
            require(a['R_A'] == b['R_A'] == .1, 'Missing actual R0.1 control')
            delta = abs(b['lowering_eV'] - a['lowering_eV'])
            tolerance = max(5., .1 * abs(b['lowering_eV']))
            require(delta <= tolerance, 'F ' + charge + ' ' + kind + ' control failed: |delta| exceeds max(5 eV,10% control lowering)')
            controls[kind][charge] = dict(
                kind=kind, charge=int(charge[1]), R_A=.1, passed=True,
                base=base['run_binding'], control=other['run_binding'],
                base_geometry=base['geometry_binding'], control_geometry=other['geometry_binding'],
                base_lowering_eV=a['lowering_eV'], control_lowering_eV=b['lowering_eV'],
                absolute_difference_eV=delta, tolerance_eV=tolerance,
                base_raw_fmax_eV_A=base['run']['geometry']['raw_fmax_eV_A'],
                control_raw_fmax_eV_A=other['run']['geometry']['raw_fmax_eV_A'],
                scope='Actual independently force-qualified finite FD geometry; no inherited cell waiver.')
    return controls


def _ratio(numerator, denominator, numerator_label, denominator_label):
    return dict(numerator=numerator_label, denominator=denominator_label,
                numerator_lowering_eV=numerator, denominator_lowering_eV=denominator,
                delta_eV=numerator - denominator,
                signed_factor=numerator / denominator if denominator != 0 else None,
                positive_lowering_factor=numerator / denominator if numerator > 0 and denominator > 0 else None,
                factor_scope='Positive-lowering factor defined only when both lowerings are positive; signed factor undefined for a zero denominator.')


def provenance(result):
    run, q = result['run'], result['quartet']
    states = q['states']
    geometry = run['geometry']
    charges = {}
    for label in ('q0', 'q1'):
        charge = int(label[1])
        charges[label] = dict(charge=charge, probe_added_electrons=1 - charge,
                             occupied_probe_electrons=[s['electrons'] for s in sorted(states[label]['occupied-probes'], key=lambda s: s['R_A'])],
                             vacant_probe_electrons=[s['electrons'] for s in sorted(states[label]['vacant-probes'], key=lambda s: s['R_A'])])
    return dict(case=q['case'], specification=run['specification'], geometry=geometry,
                volume_A3=float(np.linalg.det(np.asarray(geometry['cell_A']))),
                volume_scope='Inherited original PW volume; fresh FD ionic relaxation only' if q['case'] == 'F' else 'Finite periodic 10 A gas box, not an exact isolated limit',
                original_volume_provenance=result['original_volume_provenance'],
                occupied_host_electrons=states['occupied-host']['electrons'],
                vacant_host_electrons=states['vacant-host']['electrons'],
                vacant_host_backend=states['vacant-host']['backend'], charges=charges,
                actual_kpoints=states['occupied-host']['actual_kpoints'],
                realspace_shape=states['occupied-host']['realspace_shape'],
                nbands=states['occupied-host']['nbands'], BO_setups=states['occupied-host']['setups'])


def summarize(f, g, fine, large, g2=None):
    """Recompute all quartet values and require both h/cell tests for BOTH charges."""
    cases = {name: load_run(path, case) for name, path, case in
             (('F', f, 'F'), ('G', g, 'G'), ('F_h', fine, 'F'), ('F_cell', large, 'F'))}
    if g2 is not None:
        cases['G2'] = load_run(g2, 'G2')
    for name in ('F', 'G', 'G2'):
        if name in cases:
            _complete(cases[name], name)
    require(cases['G']['run']['specification']['h_A'] == .2, 'Gas must use base h0.20')
    if 'G2' in cases:
        require(cases['G2']['run']['specification']['h_A'] == .2, 'G2 must use base h0.20')
    for name, result in cases.items():
        require(result['source_manifest'] == cases['F']['source_manifest'],
                name + ': calculations must use the same complete frozen executable source epoch')
        compare_json(result['run']['specification']['method'], cases['F']['run']['specification']['method'],
                     name + ': same BO method')
    controls = _controls(cases['F'], cases['F_h'], cases['F_cell'])
    rows = {name + '_' + charge: result['rows'][charge]
            for name, result in cases.items() if name in ('F', 'G', 'G2') for charge in ('q0', 'q1')}
    ratios = {}
    for charge in ('q0', 'q1'):
        value = rows['F_' + charge][0]['lowering_eV']
        ratios['F_' + charge + '_over_PdD40p5'] = _ratio(value, 40.5, 'F_' + charge, 'stated PdD 40.5 eV (original bare probe)')
        ratios['F_' + charge + '_over_G_' + charge] = _ratio(value, rows['G_' + charge][0]['lowering_eV'], 'F_' + charge, 'G_' + charge)
    result = dict(schema='hydride-screening-readout-v1', qualified=True, rows=rows,
                  F_controls=controls, ratios_R01=ratios,
                  sources={name: value['sources'] for name, value in cases.items()},
                  provenance={name: provenance(value) for name, value in cases.items()},
                  charge_gauge_scope=CHARGE_SCOPE,
                  PdD_reference=dict(stated_lowering_eV=40.5, original_probe_charge=1,
                                     scope='Stated comparison value, not a new same-code PdD calculation; q0 comparison changes probe charge.'),
                  scope=SCOPE)
    json.dumps(result, allow_nan=False)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('f', 'g', 'fine', 'large', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--g2', type=Path)
    args = parser.parse_args()
    require(not args.output.exists(), 'Output readout must be new')
    result = summarize(args.f, args.g, args.fine, args.large, args.g2)
    with args.output.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + '\n')
    print(json.dumps(dict(output=args.output.name, qualified=True, source=binding(args.output)), allow_nan=False))


if __name__ == '__main__':
    main()
