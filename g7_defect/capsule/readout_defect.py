"""Replay actual neutral quartets and admit only the registered C h/cell comparisons."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

SUPPORT = Path(__file__).resolve().parent
if not (SUPPORT/'green.py').is_file():
    SUPPORT = SUPPORT.parent/'delivery'
sys.path.insert(0, str(SUPPORT))
from green import BOHR_A, COULOMB_EV_A, PeriodicNuclearGreen
from capsule_io import require

RADII = [.1, .2, .3, .5, 1., 1.6]


def binding(path):
    data = path.read_bytes()
    return dict(path=path.name, bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def close(actual, expected, label, tolerance=1e-9):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(), label + ': nonfinite/shape')
    require(bool(np.all(np.abs(a-b) <= tolerance)), label + ': changed values')


def energy(state):
    require(state['scf_converged'] and state['xc'] == 'PBE' and state['setups'] == {'Pd': '10', 'H': 'ae'},
            'Every quartet energy must be converged PBE/Pd10/Hae')
    expected = sum(10 if symbol == 'Pd' else 1 for symbol in state['symbols']) - state['charge']
    close(state['electrons'], expected, 'Actual electron count')
    e = state['energy']
    close(e['internal_energy_eV'], e['free_energy_eV']-e['minus_TS_eV'], 'Internal energy definition')
    return e['internal_energy_eV']


def load_run(directory, case, radii):
    run = json.loads((directory / 'run-result.json').read_bytes())
    require(run['completed'] and run['case'] == case, 'Wrong/incomplete Native case')
    require(binding(directory / 'bo/quartet.json') == run['quartet'], 'Changed Native quartet bytes')
    q = json.loads((directory / 'bo/quartet.json').read_bytes())
    relaxation = run['relaxation']
    require(relaxation['passed'] and relaxation['raw_fmax_eV_A'] <= .1, 'Unqualified geometry under Axiom2518')
    require(q['case'] == case and q['specification']['geometry_H_setup'] == 'standard_PAW', 'Wrong geometry route')
    close(q['path']['target_A'], relaxation['positions_A'][relaxation['target_index']], 'Actual relaxed target')
    close(q['path']['direction'], relaxation['path']['direction'], 'Actual relaxed ray')
    close(np.linalg.norm(q['path']['direction']), 1., 'Unit ray')
    states = q['states']; occupied, vacant = states['occupied-host'], states['vacant-host']
    eo, ev = energy(occupied), energy(vacant)
    target = relaxation['target_index']
    require(occupied['symbols'][target] == 'H', 'Wrong target species')
    require(vacant['symbols'] == occupied['symbols'][:target]+occupied['symbols'][target+1:], 'Target-only removal')
    close(vacant['positions_A'], np.delete(occupied['positions_A'], target, axis=0), 'Frozen vacant host')
    close(occupied['electrons']-vacant['electrons'], 1., 'Neutral target removal')
    require(occupied['charge'] == vacant['charge'] == 0, 'Hosts must be neutral')
    close(occupied['positions_A'], relaxation['positions_A'], 'Frozen relaxed occupied host')
    close(occupied['cell_A'], relaxation['cell_A'], 'Relaxed host cell')
    close(vacant['cell_A'], occupied['cell_A'], 'Frozen vacant-host cell')
    green = PeriodicNuclearGreen(occupied['cell_A']); width = math.sqrt(5.) / BOHR_A
    replayed = []
    require(len(states['occupied-probes']) == len(states['vacant-probes']) == len(q['rows']), 'Unpaired probe support')
    for o, v, saved in zip(states['occupied-probes'], states['vacant-probes'], q['rows']):
        radius = o['R_A']
        require(radius == v['R_A'] == saved['R_A'], 'Mismatched actual radii')
        point = np.asarray(q['path']['target_A']) - radius*np.asarray(q['path']['direction'])
        for probe, host in ((o, occupied), (v, vacant)):
            require(probe['charge'] == 1 and probe['symbols'] == host['symbols']+['H'], 'Wrong charged probe inventory')
            close(probe['positions_A'][:-1], host['positions_A'], 'Frozen probe background')
            close(probe['positions_A'][-1], point, 'Head-on probe position')
            close(probe['cell_A'], occupied['cell_A'], 'Same quartet cell')
            close(probe['electrons'], host['electrons'], 'Probe must not add an electron')
        raw = math.fsum([energy(o), -energy(v), -eo, ev])
        correction = green.point_minus_gaussian_eV(-radius*np.asarray(q['path']['direction']), width)
        total = raw + correction
        row = dict(R_A=radius, BO_quartet_soft_eV=raw, point_minus_gaussian_eV=correction,
                   V_pair_total_eV=total, lowering_eV=COULOMB_EV_A/radius-total)
        for key, value in row.items():
            close(saved[key], value, 'Saved quartet '+key)
        replayed.append(row)
    replayed.sort(key=lambda row: row['R_A'])
    require([row['R_A'] for row in replayed] == radii, 'Wrong requested support')
    return dict(run=run, quartet=q, rows=replayed, run_binding=binding(directory / 'run-result.json'),
                quartet_binding=binding(directory / 'bo/quartet.json'))


def compare_h(base, control):
    qb, qc = base['quartet'], control['quartet']
    sb, sc = qb['specification'], qc['specification']
    require(sb['case'] == sc['case'] == 'C', 'Grid controls must belong to C')
    close(qb['path']['direction'], qc['path']['direction'], 'Matched grid-control ray')
    close(qb['states']['occupied-host']['actual_kpoints'], qc['states']['occupied-host']['actual_kpoints'], 'Matched Gamma sampling')
    require(sb['h_A'] == .2 and sc['h_A'] == .16 and sb['repeat_primitive'] == sc['repeat_primitive'] == 2, 'Wrong two h values')
    close(base['run']['relaxation']['positions_A'], control['run']['relaxation']['positions_A'], 'Fixed h-control geometry')
    close(qb['states']['occupied-host']['cell_A'], qc['states']['occupied-host']['cell_A'], 'Fixed h-control cell')
    gb, gc = np.array(qb['states']['occupied-host']['realspace_shape']), np.array(qc['states']['occupied-host']['realspace_shape'])
    require(bool(np.all(gc >= gb) and np.any(gc > gb)), 'Nominal h change did not refine the actual grid')
    a, b = base['rows'][0]['lowering_eV'], control['rows'][0]['lowering_eV']
    tolerance = max(5., .1*abs(b)); difference = abs(b-a)
    return dict(kind='h', base=base['run_binding'], control=control['run_binding'], R_A=.1,
                base_lowering_eV=a, control_lowering_eV=b, absolute_difference_eV=difference,
                tolerance_eV=tolerance, passed=difference <= tolerance)


def summarize(a, c, d, fine, baseline, cell_sensitivity=None):
    cases = {case: load_run(path, case, RADII) for case, path in (('A', a), ('C', c), ('D', d))}
    h = compare_h(cases['C'], load_run(fine, 'C', [.1]))
    cell = dict(status='omitted', converged=False, authority='Axiom revised cell-scope decision, 2026-11-20',
                scope='Converged two-cell qualification was stopped and waived; cell size remains an omission.')
    if cell_sensitivity is not None:
        from cell_sensitivity import replay
        sensitivity = replay(cell_sensitivity)
        require(sensitivity['base_force_frame'] == cases['C']['run']['relaxation']['source_frame'], 'Sensitivity belongs to another local C geometry')
        close(sensitivity['direction'], cases['C']['quartet']['path']['direction'], 'Same sensitivity ray')
        close(sensitivity['cell_A'], 1.5*np.asarray(cases['C']['quartet']['states']['occupied-host']['cell_A']), '27-site sensitivity cell')
        base_value = cases['C']['rows'][0]['lowering_eV']; other = sensitivity['row']['lowering_eV']
        cell = dict(status='unrelaxed_sensitivity', converged=False, R_A=.1,
                    source=binding(Path(cell_sensitivity)/'cell-sensitivity.json'),
                    base_lowering_eV=base_value, sensitivity_lowering_eV=other,
                    absolute_difference_eV=abs(other-base_value), original_tolerance_eV=max(5., .1*abs(other)),
                    within_original_tolerance=abs(other-base_value)<=max(5., .1*abs(other)),
                    scope='Unrelaxed far field and embedded 8-site local geometry; NOT a converged cell-size check.')
    b = json.loads(baseline.read_bytes())
    require(b['schema'] == 'coherencev1-target-conditioned-pair-v1' and b['observable'] == 'target_conditioned_BO_total_interaction', 'Wrong reused B observable')
    require(b['R_A'] == RADII and len(b['V_pair_total_eV']) == len(RADII), 'B must be exact six-point reuse')
    b_rows = [dict(R_A=r, V_pair_total_eV=v, lowering_eV=COULOMB_EV_A/r-v) for r, v in zip(b['R_A'], b['V_pair_total_eV'])]
    table = {case: value['rows'] for case, value in cases.items()}; table['B'] = b_rows
    at_point = {case: rows[0]['lowering_eV'] for case, rows in table.items()}
    contrasts = {}
    for name, numerator, denominator in (('defect_low', 'C', 'A'), ('defect_high', 'D', 'B'), ('loading_perfect', 'B', 'A'), ('loading_vacancy', 'D', 'C')):
        x, y = at_point[numerator], at_point[denominator]
        contrasts[name] = dict(numerator=numerator, denominator=denominator, delta_eV=x-y,
                               positive_lowering_factor=x/y if x > 0 and y > 0 else None)
    point_qualified = h['passed']
    triggers = [dict(defect=x, baseline=y) for x in ('C', 'D') for y in ('A', 'B')
                if at_point[x] > 0 and (at_point[y] <= 0 or at_point[x] > 3*at_point[y])]
    return dict(schema='pd-defect-screening-readout-v2', point_qualified=point_qualified, cell_qualified=False,
                rows=table, contrasts_R01=contrasts, C_controls=dict(h=h, cell=cell), sample_three_times_conditions=triggers,
                E_required=bool(triggers) if point_qualified else None, B_source=binding(baseline),
                sources={case: dict(run=value['run_binding'], quartet=value['quartet_binding']) for case, value in cases.items()},
                qualification_scope='Force-qualified actual A/C/D and representative C grid check; NO converged cell-size qualification, per Axiom revised scope.',
                scope='Static BO point estimates at the common inherited volume; no U0, rate, thermal ensemble, capacity or infinite-cell conclusion.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('a', 'c', 'd', 'fine', 'baseline', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--cell-sensitivity', type=Path)
    args = parser.parse_args()
    result = summarize(args.a, args.c, args.d, args.fine, args.baseline, args.cell_sensitivity)
    with args.output.open('x') as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(args.output), point_qualified=result['point_qualified'],
                         cell_qualified=False, E_required=result['E_required']), indent=2))
