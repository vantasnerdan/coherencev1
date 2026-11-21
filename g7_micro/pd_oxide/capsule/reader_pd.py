"""Stdlib/NumPy replay of source-bound compact Pd/PdO saved Native outputs.

No Native execution. The energy, state, epoch, Green and public-Z2 operator
implementations are the accepted common reader helpers, not a second model.
"""
import argparse
import itertools
from datetime import datetime, timezone
import json
import math
import re
from pathlib import Path

import numpy as np
import reader
from reader import (RADII, PeriodicNuclearGreen, binding, bound_file, close,
                    compare_json, number, read_json, require)

VALENCE = {'Pd': 10, 'O': 6, 'H': 1, 'He': 2}
CELL = np.diag([6.13639979955674, 6.13639979955674, 22.039528899485433])
FIXED = list(range(10, 15))
ROLES = tuple(v + suffix for suffix in ('', '-h', '-k') for v in ('plain', 'oxygen'))
EPOCH_SOURCES = ('native.py', 'pd_run.py', 'pd_relax.py', 'pd_ray.py', 'geometry_pd.py',
                 'run_zr.py', 'relax_zr.py', 'quartet_zr.py', 'contact_zr.py',
                 'contact_setup.py', 'collision_bare_reference.py')
CRITICAL_SOURCES = (
    'campaigns/P255-pdd-collision/micro_calibration/native.py',
    'campaigns/P255-pdd-collision/micro_calibration/quartet_zr.py',
    'campaigns/P255-pdd-collision/micro_calibration/contact_zr.py',
    'campaigns/P255-pdd-collision/hydride_screening/contact_setup.py',
    'src/substrate_framework/collision_bare_reference.py',
    'campaigns/P255-pdd-collision/collision_ground.py',
    'campaigns/P255-pdd-collision/delivery/green.py',
)
SCOPE = ('Actual finite-cell target-conditioned internal BO screening for compact Pd and a coherent '
         'added 4Pd+4O sheet, not oxygen insertion at fixed Pd inventory. D is unconstrained and its '
         'actual migrated/relaxed site is used. Pd10/O6 frozen-core PAW; q+1 Hae probe adds no electron. '
         'Paired h/k checks apply only at R0.1/contact. Slab-thickness, cell/vacuum and full-radius '
         'h/k controls omitted. Geometry-convergence error is unbounded (geometry_error_bound_eV=null), '
         'including for force-qualified frames; residual forces do not bound pair-lowering error. '
         'No isolated limit, site-retention, bulk-PdO, Ueq, action/Eq7 or experimental phase claim.')



def _recorded_stop(spec):
    stop = datetime.fromisoformat(spec['absolute_stop_utc'].replace('Z', '+00:00'))
    require(spec['authority_entry'] == 'entry_2841' and stop.tzinfo is not None and
            stop <= datetime(2026, 11, 21, 3, tzinfo=timezone.utc),
            'Changed commissioned Pd authority or extended recorded hard stop')

def _spec(spec, run):
    require(spec['material'] == 'pd' and spec['phase'] == 'bo' and spec['size'] == 'compact' and
            spec['variant'] in ('plain', 'oxygen') and spec['case'] == run['case'], 'Wrong commissioned compact Pd case')
    require(spec['h_A'] in (.2, .16) and spec['kmesh'] in ([1, 1, 1], [2, 2, 1]) and
            spec['radii_A'] in (RADII, [.4, .2, .1], [.1], [.1, .2, .4], [.8, 1.2, 1.6]), 'Wrong commissioned Pd numerical control')
    if spec['radii_A'] not in (RADII, [.1]):
        require(run.get('schema') == 'micro-pd-bo-branch-v1' or
                run.get('container_schema') == 'micro-pd-parallel-case-v1',
                'Extra radius subsets require explicit actual branch/composite context')
    require(spec['force_tolerance_eV_A'] == .05 and spec['initial_geometry'] is None and
            spec['selected_geometry'] is not None and spec['probe_ray'] is not None,
            'Missing selected Pd geometry/common ray or changed force policy')
    _recorded_stop(spec)
    require(spec.get('ray_policy_entry') == 'entry_2971' and
            spec.get('geometry_policy_entry') in ('entry_2904', 'entry_3001'), 'Missing commissioned Pd ray/geometry authority')
    require('charge' not in spec or type(spec['charge']) is int and spec['charge'] == 1,
            'Only q+1/no-added-electron Pd probe commissioned')


def _frame(g, variant):
    symbols = ['Pd'] * 15 + (['Pd', 'Pd', 'O', 'O'] * 2 if variant == 'oxygen' else []) + ['H']
    require(g['symbols'] == symbols and g['target_index'] == len(symbols) - 1,
            'Wrong ordered substrate/oxide/D inventory or target')
    positions = np.asarray(g['positions_A'], float)
    cell = np.asarray(g['cell_A'], float)
    masses = np.asarray(g['masses_u'], float)
    require(positions.shape == (len(symbols), 3) and np.isfinite(positions).all(), 'Invalid Pd positions')
    close(cell, CELL, 'Commissioned identical compact cell')
    PeriodicNuclearGreen(cell)
    close(masses, [{'Pd': 106.42, 'O': 15.999, 'H': 2.01410177812}[s] for s in symbols], 'Pd/O/D isotope masses')
    return positions, cell, masses


def _forces(record, shape, *, geometric_witness=False, labelled=False):
    free = [i for i in range(shape[0]) if i not in FIXED]
    require(record['fixed_indices'] == FIXED and record['free_indices'] == free and shape[0]-1 in free,
            'Only bottom five substrate Pd may be fixed; D must be free')
    raw = np.asarray(record['raw_forces_eV_A'], float)
    constrained = np.asarray(record['constrained_forces_eV_A'], float)
    require(raw.shape == constrained.shape == shape and np.isfinite(raw).all() and np.isfinite(constrained).all(),
            'Missing finite raw/constrained force arrays')
    expected = raw.copy()
    expected[FIXED] = 0.
    close(constrained, expected, 'Fixed-bottom-only force operation')
    close(record['forces_eV_A'], raw, 'Reported raw PAW forces')
    norms = np.linalg.norm(raw, axis=1)
    for key, value in dict(raw_fmax_eV_A=norms.max(), free_fmax_eV_A=norms[free].max(),
                           fixed_fmax_eV_A=norms[FIXED].max(), max_force_eV_A=norms.max()).items():
        close(record[key], float(value), 'Actual force diagnostic ' + key)
    free_max, raw_max = float(norms[free].max()), float(norms.max())
    require(geometric_witness or labelled or free_max <= .30, 'Pd actual free force exceeds conditional geometry ceiling')
    result = dict(force_qualified=bool(free_max <= .05) if not geometric_witness else False,
                  relaxation_status=('geometric-witness-not-final' if geometric_witness else
                                     'force-qualified' if free_max <= .05 else 'partly-relaxed'),
                  free_fmax_eV_A=free_max, raw_fmax_eV_A=raw_max, geometry_error_bound_eV=None)
    for key in ('force_qualified', 'relaxation_status', 'geometry_error_bound_eV'):
        if key in record:
            require(type(record[key]) is type(result[key]) and record[key] == result[key],
                    'Mislabeled Pd force record ' + key)
    return result

def _qualification(g):
    """Derive status from arrays, never from a producer's eligibility boolean."""
    positions, _, _ = _frame(g, g['specification']['variant'])
    labelled = isinstance(g.get('labelled_admission'), dict)
    result = _forces(g, positions.shape, labelled=labelled)
    if not result['force_qualified']:
        require(all(key in g for key in ('force_qualified', 'relaxation_status', 'geometry_error_bound_eV')),
                'Partly relaxed Pd geometry lacks explicit qualification labels')
        require(isinstance(g.get('partial_admission'), dict) or labelled, 'Missing explicit Pd geometry admission')
    return result

def _partial_admission(root, g, raw, positions, cell):
    from reader_pd_chain import partial
    partial(root, g, raw, positions, cell)



def _neighbors(positions, cell, target, cutoff):
    delta = positions - positions[target]
    sigma = float(np.linalg.svd(cell, compute_uv=False)[-1])
    shell = math.ceil((cutoff + float(np.linalg.norm(delta, axis=1).max())) / sigma)
    result = {}
    for shift in itertools.product(range(-shell, shell+1), repeat=3):
        distances = np.linalg.norm(delta + np.asarray(shift) @ cell, axis=1)
        for index in np.flatnonzero(distances < cutoff):
            if index != target:
                result[int(index), shift] = float(distances[index])
    return result


def _diagnostics(g, initial, positions, cell):
    target, symbols = g['target_index'], g['symbols']
    displacement = positions - np.asarray(initial['positions_A'], float)
    mic = displacement - np.rint(displacement @ np.linalg.inv(cell)) @ cell
    d = g['displacements']
    close(d['vectors_A'], mic, 'Actual initial-to-relaxed MIC displacement')
    lengths = np.linalg.norm(mic, axis=1)
    close(d['distances_A'], lengths, 'Actual displacement distances')
    close(d['maximum_A'], float(lengths.max()), 'Actual maximum displacement')
    close(d['target_A'], float(lengths[target]), 'Actual D migration distance')
    close(d['unwrapped_target_vector_A'], displacement[target], 'Actual unwrapped D displacement')
    c = g['coordination']
    require(c['target_index'] == target and c['target_symbol'] == 'H' and
            c['physical_site'] == 'actual relaxed site; migration permitted, no retained-site qualification',
            'False retained-site label or wrong target diagnostics')
    close(c['actual_target_A'], positions[target], 'Actual relaxed D site')
    actual = _neighbors(positions, cell, target, 3.5)
    saved = {}
    for row in c['neighbors_within_3_5_A']:
        index, image = row['index'], tuple(row['image'])
        require(type(index) is int and 0 <= index < len(symbols) and len(image) == 3 and
                all(type(i) is int for i in image) and row['symbol'] == symbols[index] and
                (index, image) not in saved, 'Invalid/duplicate Pd coordination neighbor')
        saved[index, image] = row['distance_A']
    require(set(saved) == set(actual), 'Incomplete actual Pd/O target coordination')
    for key, value in actual.items():
        close(saved[key], value, 'Actual periodic target-neighbor distance')
    for symbol in ('Pd', 'O'):
        expected = sorted(i for (i, _), distance in actual.items() if symbols[i] == symbol and distance < 3.)
        require(c['final_' + symbol + '_indices_within_3_A'] == expected, 'Wrong final ' + symbol + ' coordination')
        candidates = [i for i, s in enumerate(symbols) if s == symbol]
        row = c['nearest_' + symbol]
        if not candidates:
            require(row is None, 'Unexpected nearest ' + symbol)
            continue
        require(type(row['index']) is int and row['index'] in candidates and row['symbol'] == symbol and
                len(row['image']) == 3 and all(type(i) is int for i in row['image']), 'Invalid nearest ' + symbol)
        vectors = positions[candidates] - positions[target]
        vectors -= np.rint(vectors @ np.linalg.inv(cell)) @ cell
        closest = float(np.linalg.norm(vectors, axis=1).min())
        witness = positions[row['index']] - positions[target] + np.asarray(row['image']) @ cell
        close(row['distance_A'], closest, 'Actual nearest ' + symbol + ' distance')
        close(float(np.linalg.norm(witness)), closest, 'Actual nearest ' + symbol + ' image')
    before = _neighbors(np.asarray(initial['positions_A'], float), cell, target, 3.)
    require(c['initial_Pd_indices_within_3_A'] == sorted(i for i, _ in before if symbols[i] == 'Pd'),
            'Changed initial target Pd coordination')


def _clearance(saved, g, direction):
    positions, cell = np.asarray(g['positions_A'], float), np.asarray(g['cell_A'], float)
    target = g['target_index']
    minimum = reader._segment_clearance(positions, cell, target, np.asarray(direction), radius_min=0.)
    close(saved['min_clearance_A'], minimum, 'All-image continuous non-target ray clearance')
    require(saved['R_interval_A'] == [0., 1.6], 'Changed continuous clearance interval')
    index, image = saved['spectator_index'], saved['spectator_image']
    require(type(index) is int and 0 <= index < len(positions) and index != target and
            saved['spectator_symbol'] == g['symbols'][index] and len(image) == 3 and
            all(type(i) is int for i in image), 'Invalid clearance witness')
    vector = positions[index] - positions[target] + np.asarray(image) @ cell
    radius = float(np.clip(-vector @ direction, 0., 1.6))
    close(saved['at_R_A'], radius, 'Exact continuous-segment projection radius')
    close(float(np.linalg.norm(vector + radius * np.asarray(direction))), minimum, 'Clearance witness distance')
    spectators = np.delete(positions, target, axis=0) - positions[target]
    spectators -= np.rint(spectators @ np.linalg.inv(cell)) @ cell
    close(saved['search_radius_A'], 1.6 + float(np.linalg.norm(spectators, axis=1).min()) + 1e-9,
          'Proved finite periodic-image search sphere')
    return minimum


def _candidates(original):
    candidates = [('vacuum-side-normal', [0., 0., -1.]), ('source-original-in-plane', original)]
    for theta in range(15, 181, 15):
        for phi in range(0, 360, 15):
            if theta == 180 and phi:
                continue
            t, p = np.deg2rad([theta, phi])
            n = [float(np.sin(t)*np.cos(p)), float(np.sin(t)*np.sin(p)), float(-np.cos(t))]
            if not any(np.allclose(n, row[1], rtol=0, atol=1e-12) for row in candidates):
                candidates.append((f'polar-{theta:03d}-azimuth-{phi:03d}', n))
    return candidates


def _physical_sources(epoch):
    files = epoch['files']
    result = {}
    for path in CRITICAL_SOURCES:
        matches = [item for item in files if Path(item['path']).name == Path(path).name]
        require(len(matches) == 1 and matches[0]['path'] == path,
                'Missing/ambiguous physical source ' + path)
        result[path] = {key: matches[0][key] for key in ('bytes', 'sha256')}
    return result


def _policy_files(root, policy):
    files = policy['evidence_files']
    require(isinstance(files, list) and files, 'Missing immutable ray evidence closure')
    names = set()
    for item in files:
        require(item['path'] not in names, 'Duplicate ray evidence')
        names.add(item['path'])
        bound_file(root, item, 'Ray evidence')
    producer = bound_file(root, policy['producer'], 'Original ray producer')
    require(producer.name == 'pd_ray.py' and policy['producer'] in files, 'Unbound ray producer')


def _witness(root, witness):
    require(witness['schema'] == 'micro-pd-coated-boundary-witness-v1' and
            witness['not_final'] is True and witness['force_qualified'] is False and
            witness['bo_eligible'] is False, 'Coated boundary is geometric evidence only')
    paths = {key: bound_file(root, witness[key], 'Coated witness ' + key)
             for key in ('boundary', 'gpw', 'initial_frame', 'original_input', 'original_source_epoch')}
    boundary, initial = read_json(paths['boundary']), read_json(paths['initial_frame'])
    spec = read_json(paths['original_input'])
    compare_json(boundary['specification'], spec, 'Witness actual consumed specification')
    _recorded_stop(spec)
    require(spec['material'] == 'pd' and spec['phase'] == 'geometry' and spec['variant'] == 'oxygen' and
            spec['size'] == 'compact' and spec['h_A'] == .2 and spec['kmesh'] == [1, 1, 1] and
            spec['initial_geometry'] is None and spec['force_tolerance_eV_A'] == .05,
            'Wrong source-built coated geometric witness')
    require(initial['schema'] == 'micro-pd-initial-frame-v1', 'Wrong witness source initial frame')
    frame = dict(boundary, target_index=initial['target_index'])
    positions, cell, masses = _frame(frame, 'oxygen')
    ip, ic, im = _frame(initial, 'oxygen')
    close(ip[FIXED], positions[FIXED], 'Witness unmoved bottom substrate')
    close(ic, cell, 'Witness initial cell')
    close(im, masses, 'Witness initial masses')
    require(initial['fixed_indices'] == FIXED and initial['free_indices'] == boundary['free_indices'] and
            list(np.flatnonzero(np.abs(ip[:15, 2] - ip[:15, 2].min()) < 1e-10)) == FIXED,
            'Wrong witness fixed/free source frame')
    design = initial['initial_design']
    require(design['variant'] == 'oxygen' and design['pair'] == 'pd-oxide' and
            design['size'] == 'compact' and design['target_index'] == initial['target_index'] and
            design['matching']['metal_Pd'] == 15 and design['matching']['total_Pd'] == 19 and
            design['matching']['O'] == 4 and design['matching']['D'] == 1,
            'Wrong witness source intervention')
    close(design['cell_control']['cell_A'], cell, 'Witness source cell')
    close(ip[initial['target_index']], design['nominal_site']['initial_target_A'], 'Witness source target')
    epoch = read_json(paths['original_source_epoch'])
    require(epoch['schema'] == 'micro-frozen-repository-source-v1' and epoch['files'],
            'Missing witness original source closure')
    require(all(witness['original_source_epoch'][k] == spec['frozen_source_manifest'][k]
                for k in ('bytes', 'sha256')), 'Witness source epoch differs from consumed input')
    names = set()
    for item in epoch['files']:
        require(item['path'] not in names, 'Duplicate witness original source')
        names.add(item['path'])
        bound_file(paths['original_source_epoch'].parent / 'source', item, 'Witness original source')
    require(all(sum(Path(name).name == name0 for name in names) == 1 for name0 in EPOCH_SOURCES),
            'Incomplete witness producer closure')
    _physical_sources(epoch)
    for item in witness['source_files']:
        bound_file(root, item, 'Witness source file')
    expected_sources = []
    for item in epoch['files']:
        source_path = paths['original_source_epoch'].parent / 'source' / item['path']
        expected_sources.append(dict(path=source_path.relative_to(root).as_posix(),
                                     **{k: item[k] for k in ('bytes', 'sha256')}))
    compare_json(witness['source_files'], expected_sources, 'Complete ordered witness original source closure')
    for item in design['structural_sources']:
        name = item.get('relative_name', item.get('path'))
        matches = [row for row in epoch['files'] if row['path'] == name]
        require(len(matches) == 1 and all(matches[0][k] == item[k] for k in ('bytes', 'sha256')),
                'Witness structural builder differs from original epoch')
    reader._no_extra_electrons(boundary)
    require(boundary['schema'] == 'micro-native-state-v1' and boundary['backend'] == 'gpaw-native' and
            boundary['scf_converged'] is True and boundary['geometry_PAW'] is True and
            boundary['charge'] == 0 and boundary['xc'] == 'PBE' and boundary['spinpol'] is False and
            type(boundary['accepted_step']) is int and boundary['accepted_step'] >= 0 and
            number(boundary['process_cpu_s'], 'Witness boundary CPU') > 0, 'Not an accepted actual Native boundary')
    close(boundary['electrons'], sum(VALENCE[s] for s in boundary['symbols']), 'Witness actual Ne')
    close(boundary['occupation_width_eV'], .1, 'Witness actual smearing')
    lengths = 1 / np.linalg.norm(np.linalg.inv(cell), axis=0)
    close(boundary['realspace_shape'], np.maximum(4, (lengths / spec['h_A'] / 4 + .5).astype(int) * 4),
          'Witness actual grid')
    require(boundary['actual_kpoints'] == 1 and boundary['nbands'] == 116 and
            type(boundary['scf_iterations']) is int and 0 < boundary['scf_iterations'] <= 150 and
            spec['eigensolver'] == {'name': 'davidson', 'niter': 2} and boundary['eigensolver'] == 'Davidson',
            'Wrong witness SCF method')
    require(set(boundary['actual_setups']) == {'Pd', 'O', 'H'}, 'Missing witness actual PAW setups')
    for symbol, setup in boundary['actual_setups'].items():
        close(setup['Nv'], VALENCE[symbol], 'Witness PAW valence')
        require(type(setup['ni']) is int and setup['ni'] > 0 and setup['Nc'] >= 0,
                'Wrong witness PAW projectors')
    energy = boundary['energy']
    require(set(energy) == {'internal_energy_eV', 'free_energy_eV', 'minus_TS_eV'}, 'Wrong witness energy fields')
    for key, value in energy.items():
        number(value, 'Witness ' + key)
    close(energy['internal_energy_eV'], energy['free_energy_eV'] - energy['minus_TS_eV'], 'Witness internal energy')
    require(energy['minus_TS_eV'] <= 1e-9 and
            all(witness['gpw'][k] == boundary['gpw'][k] for k in ('bytes', 'sha256')),
            'Wrong witness energy/checkpoint binding')
    _forces(boundary, positions.shape, geometric_witness=True)
    return frame


def _precommit(path, policy, inventory, archive_root):
    root = path.parent
    require(policy['schema'] == 'micro-pd-precommitted-shared-ray-v1' and
            policy['authority_entry'] == 'entry_2971', 'Wrong precommitted ray authority')
    require(policy['label'] == 'precommitted ray; final coated shared-open admission pending',
            'Precommit falsely claims final coated admission')
    _policy_files(root, policy)
    plain_path = bound_file(root, policy['plain_geometry'], 'Qualified bare precommit frame')
    plain = read_json(plain_path)
    _relaxed_geometry(plain_path.parent, plain, plain, plain_path.parent, inventory, archive_root)
    require(_qualification(plain)['force_qualified'], 'Precommit bare frame must be force-qualified')
    coated = _witness(root, policy['coated_witness'])
    frames = {'plain': plain, 'oxygen': coated}
    close(plain['cell_A'], coated['cell_A'], 'Current precommit identical cells')
    definition = policy['candidate_definition']
    require(definition['polar_axis'] == '-z' and definition['polar_degrees'] == list(range(15, 181, 15)) and
            definition['azimuth_degrees'] == list(range(0, 360, 15)) and definition['poles_deduplicated'] is True and
            definition['extra'] == ['vacuum-side normal', 'source original in-plane direction'] and
            definition['R_interval_A'] == [0., 1.6], 'Changed declared precommit candidate set')
    require(definition['objective'] ==
            'maximize minimum clearance across qualified bare and current coated boundary',
            'Changed current-frame precommit objective')
    expected = _candidates(plain['initial_design']['direction'])
    candidates = policy['candidates']
    require(len(candidates) == len(expected), 'Missing/added precommit candidates')
    scores = []
    for row, (name, direction) in zip(candidates, expected):
        require(row['name'] == name and set(row['clearance_by_variant']) == {'plain', 'oxygen'},
                'Wrong deterministic precommit candidate ordering')
        close(row['direction_n_hat'], direction, 'Precommit candidate direction')
        score = min(_clearance(row['clearance_by_variant'][v], frames[v], direction) for v in frames)
        close(row['minimum_paired_clearance_A'], score, 'Current-frame common clearance objective')
        scores.append(score)
    preferred = scores[0] >= 1.
    selected = 0 if preferred else max(range(len(scores)), key=scores.__getitem__)
    chosen, selection = candidates[selected], policy['selection']
    require(selection['normal_preferred'] is preferred and selection['candidate_name'] == chosen['name'] and
            selection['preferred_clearance_A'] == 1. and selection['collision_epsilon_A'] == 1e-6 and
            scores[selected] > 1e-6, 'Precommit normal-first/max-common current-frame ray collides')
    require(selection['rule'] ==
            'normal first if paired clearance >=1 A; otherwise maximum paired clearance over declared finite set; first exact tie wins'
            and selection['claim'] == 'current-frame finite-set selection, not maximum-common-final',
            'False final-frame optimum claim')
    require(policy['direction_n_hat'] == chosen['direction_n_hat'], 'Changed exact precommitted direction')
    compare_json(policy['clearance_by_variant'], chosen['clearance_by_variant'], 'Precommitted clearance')
    close(policy['minimum_paired_clearance_A'], scores[selected], 'Precommitted current-frame objective')
    return plain, coated


def _ray(root, run, original, inventory, archive_root):
    g, spec = run['geometry'], run['specification']
    path = bound_file(root, g['probe_ray'], 'Selected immutable ray policy')
    require(all(g['probe_ray'][k] == spec['probe_ray'][k] for k in ('bytes', 'sha256')), 'Changed input ray policy')
    policy = read_json(path)
    precommit_path = bound_file(root, g['probe_ray_precommit'], 'Original immutable precommit')
    precommit = read_json(precommit_path)
    plain, coated = _precommit(precommit_path, precommit, inventory, archive_root)
    executed = read_json(bound_file(root, run['source_epoch'], 'Executed BO source epoch'))
    witness_epoch = read_json(bound_file(precommit_path.parent,
        precommit['coated_witness']['original_source_epoch'], 'Witness executed geometry epoch'))
    compare_json(_physical_sources(executed), _physical_sources(witness_epoch),
                 'Witness/BO identical physical source bindings')
    members = g['probe_ray_members']
    require(set(members) == ({'plain'} if spec['variant'] == 'plain' else {'plain', 'oxygen'}),
            'Bare case must not claim a final coated member')
    for key in ('bytes', 'sha256'):
        require(members['plain'][key] == precommit['plain_geometry'][key], 'Changed qualified bare ray member')
    compare_json(read_json(bound_file(root, members['plain'], 'Archived bare ray frame')), plain, 'Bare member bytes')
    for key in ('boundary', 'gpw', 'initial_frame', 'original_input', 'original_source_epoch'):
        require(all(g['coated_witness'][key][k] == precommit['coated_witness'][key][k] for k in ('bytes', 'sha256')),
                'Changed archived coated witness ' + key)
    compare_json(_witness(root, g['coated_witness']), coated, 'Actual archived coated geometric witness')
    if spec['variant'] == 'plain':
        require(g['probe_ray'] == g['probe_ray_precommit'] and
                g['direction_status'] == 'precommitted ray; final coated shared-open admission pending',
                'Bare ray must retain pending final-coated status')
        compare_json(plain, original, 'Precommitted selected bare frame')
        clearance = precommit['clearance_by_variant']['plain']
    else:
        require(policy['schema'] == 'micro-pd-precommitted-shared-admission-v1' and
                policy['authority_entry'] == 'entry_2971' and policy['label'] == 'precommitted shared-open' and
                g['direction_status'] == 'precommitted shared-open', 'Wrong final coated shared-open admission')
        _policy_files(path.parent, policy)
        admitted_precommit = bound_file(path.parent, policy['precommit'], 'Admission original precommit')
        require(all(policy['precommit'][k] == g['probe_ray_precommit'][k] for k in ('bytes', 'sha256')),
                'Admission changed precommit identity')
        compare_json(read_json(admitted_precommit), precommit, 'Admission immutable precommit')
        final_path = bound_file(path.parent, policy['coated_geometry'], 'Actual final admitted coated frame')
        compare_json(read_json(final_path), original, 'Admission actual selected coated geometry')
        require(all(members['oxygen'][k] == policy['coated_geometry'][k] for k in ('bytes', 'sha256')),
                'Changed final coated member')
        compare_json(read_json(bound_file(root, members['oxygen'], 'Final coated member')), original,
                     'Final coated member original')
        require(policy['direction_n_hat'] == precommit['direction_n_hat'], 'Changed exact admitted direction')
        clearance = policy['clearance']
        minimum = _clearance(clearance, original, precommit['direction_n_hat'])
        require(minimum > 1e-6, f'Precommitted final coated ray blocked: actual clearance {minimum:.17g} A')
    require(g['direction'] == precommit['direction_n_hat'] and g['direction_n_hat'] == g['direction'],
            'Changed exact BO precommitted direction')
    compare_json(g['ray_clearance'], clearance, 'BO actual own continuous clearance')


def _relaxed_geometry(root, g, original, state_root, inventory, archive_root):
    """Common geometry-only/BO admission; provisional rays never gate relaxation."""
    spec = original['specification']
    require(g['schema'] == 'micro-site-geometry-v1' and g['passed'] is True, 'Unqualified Pd geometry')
    positions, cell, masses = _frame(g, spec['variant'])
    qualification = _qualification(g)
    if spec.get('checkpoint_resume') is None:
        require(g['source_geometry'] is None and g.get('checkpoint_resume') is None,
                'Uncommissioned imported Pd geometry')
    elif 'partial_admission' not in original:
        from reader_pd_chain import resumed_boundaries
        resumed_boundaries(state_root, original)
    _recorded_stop(spec)
    require(spec['material'] == 'pd' and spec['phase'] == 'geometry' and spec['variant'] in ('plain','oxygen') and
            spec['size'] == 'compact' and spec['h_A'] == .2 and spec['kmesh'] == [1,1,1] and
            spec['force_tolerance_eV_A'] == .05 and spec['initial_geometry'] is None,
            'Wrong selected neutral Gamma h0.20 geometry family')
    initial = read_json(bound_file(root, g['initial_frame'], 'Bound initial Pd frame'))
    require(all(g['initial_frame'][k] == original['initial_frame'][k] for k in ('bytes','sha256')) and
            initial['schema'] == 'micro-pd-initial-frame-v1', 'Changed initial-frame source binding')
    ip, ic, im = _frame(initial, spec['variant'])
    design = g['initial_design']
    compare_json(initial['initial_design'], design, 'Bound source initial design')
    require(design['pair'] == 'pd-oxide' and design['size'] == 'compact' and design['variant'] == spec['variant'] and
            design['target_index'] == g['target_index'] and design['matching']['metal_Pd'] == 15 and
            design['matching']['total_Pd'] == (19 if spec['variant'] == 'oxygen' else 15) and
            design['matching']['O'] == (4 if spec['variant'] == 'oxygen' else 0) and design['matching']['D'] == 1,
            'Wrong source Pd/oxide intervention design')
    close(design['cell_control']['cell_A'], cell, 'Source fixed Pd cell')
    close(ip[g['target_index']], design['nominal_site']['initial_target_A'], 'Source intended initial D position only')
    require(initial['fixed_indices'] == FIXED and initial['free_indices'] == g['free_indices'] and
            list(np.flatnonzero(np.abs(ip[:15,2]-ip[:15,2].min()) < 1e-10)) == FIXED,
            'Initial fixed indices are not exactly bottom substrate layer')
    close(positions[FIXED], ip[FIXED], 'Unmoved fixed substrate positions')
    close(im, masses, 'Initial/relaxed isotope frame')
    close(ic, cell, 'Initial/relaxed cell')
    bound_file(root, g['gpw'], 'Final geometry PAW GPW')
    raw = read_json(bound_file(root, g['state_file'], 'Final Pd PAW state'))
    reader._state(state_root, raw, 'geometry-PAW', spec, g, True, inventory, archive_root, valence=VALENCE)
    require(raw['charge'] == 0 and raw['symbols'] == g['symbols'] and raw['gpw'] == original['gpw'] and
            raw['log'] == original['log'], 'Wrong selected PAW charge/inventory/checkpoint/log')
    close(raw['positions_A'], positions, 'Actual final PAW coordinates')
    close(raw['masses_u'], masses, 'Actual final PAW masses')
    _forces(raw, positions.shape, labelled='labelled_admission' in original)
    for key in ('raw_forces_eV_A','constrained_forces_eV_A','raw_fmax_eV_A','free_fmax_eV_A','fixed_fmax_eV_A'):
        compare_json(g[key], raw[key], 'Geometry/actual PAW force ' + key)
    _diagnostics(g, initial, positions, cell)
    if 'partial_admission' in original:
        _partial_admission(state_root, original, raw, positions, cell)
    if 'labelled_admission' in original:
        from reader_pd_label import validate
        validate(state_root, original, raw, positions, cell)
    return positions, cell, masses


def _geometry(root, run, inventory, archive_root):
    g, spec = run['geometry'], run['specification']
    selected = run['provenance']['selected_geometry']
    selected_path = bound_file(root, selected, 'Selected original Pd geometry')
    original = read_json(selected_path)
    source = spec['selected_geometry']
    require(all(selected[k] == source['binding'][k] for k in ('bytes','sha256')), 'Selected input geometry bytes changed')
    require(original['specification']['variant'] == spec['variant'], 'Wrong selected Pd variant')
    for key in ('schema','specification','symbols','positions_A','cell_A','masses_u','target_index','passed','initial_design',
                'source_geometry','raw_forces_eV_A','constrained_forces_eV_A','forces_eV_A','max_force_eV_A',
                'raw_fmax_eV_A','free_fmax_eV_A','fixed_fmax_eV_A','fixed_indices','free_indices',
                'displacements','coordination','force_qualification','accepted_steps'):
        compare_json(g[key], original[key], 'Selected Pd geometry ' + key)
    for key in ('force_qualified', 'relaxation_status', 'geometry_error_bound_eV'):
        if key in original:
            compare_json(g[key], original[key], 'Selected Pd geometry policy ' + key)
    for key in ('checkpoint_resume', 'checkpoint_import', 'accepted_boundaries', 'native_log', 'native_log_status'):
        require((key in g) == (key in original), 'Changed selected resume evidence ' + key)
        if key in original:
            compare_json(g[key], original[key], 'Selected archived resume evidence ' + key)
    require(('partial_admission' in g) == ('partial_admission' in original), 'Changed selected partial admission')
    if 'partial_admission' in original:
        compare_json(g['partial_admission'], original['partial_admission'], 'Selected Pd partial admission')
    require(('labelled_admission' in g) == ('labelled_admission' in original), 'Changed selected labelled admission')
    expected_policy = 'entry_3001' if 'labelled_admission' in original else 'entry_2904'
    require(spec.get('geometry_policy_entry') == expected_policy,
            'BO policy differs from actual nested geometry admission')
    if 'labelled_admission' in original:
        compare_json(g['labelled_admission'], original['labelled_admission'], 'Selected labelled original admission')
        require(spec.get('geometry_policy_entry') == 'entry_3001' and
                g.get('labelled_admission_root') == Path(selected['path']).parent.as_posix(),
                'Missing authority-labelled BO policy/root')
    if not _qualification(g)['force_qualified']:
        if 'labelled_admission' not in original:
            require(spec.get('geometry_policy_entry') == 'entry_2904' and
                    g.get('partial_admission_root') == Path(selected['path']).parent.as_posix(),
                    'Missing conditional BO policy or changed selected admission root')
    for key in ('gpw', 'state_file', 'log', 'sanitized_log', 'initial_frame'):
        require(all(g[key][k] == original[key][k] for k in ('bytes', 'sha256')),
                'Changed selected geometry companion ' + key)
        bound_file(root, g[key], 'Selected geometry companion ' + key)
    require(all(g['gpw'][k] == source['gpw_binding'][k] for k in ('bytes','sha256')), 'Selected PAW checkpoint changed')
    positions, cell, masses = _relaxed_geometry(root, g, original, selected_path.parent, inventory, archive_root)
    _ray(root, run, original, inventory, archive_root)
    return positions, cell, masses, np.asarray(g['direction'], float)


def replay_geometry(directory, *, log_inventory=None, archive_root=None):
    root = Path(directory).resolve()
    envelopes = [root / name for name in ('geometry-run.json', 'geometry-partly-relaxed-run.json', 'geometry-labelled-run.json')
                 if (root / name).is_file()]
    require(len(envelopes) == 1, 'Require exactly one original or adopted Pd geometry-run envelope')
    path = envelopes[0]
    adopted = path.name == 'geometry-partly-relaxed-run.json'
    labelled = path.name == 'geometry-labelled-run.json'
    original_binding = binding(path)
    run = read_json(path)
    expected_schema = ('micro-site-labelled-geometry-run-v1' if labelled else
                       'micro-site-adopted-geometry-run-v1' if adopted else 'micro-site-geometry-run-v1')
    require(run['schema'] == expected_schema and run['case'] == run['specification']['case'], 'Wrong Pd geometry-run envelope')
    g = read_json(bound_file(root, run['geometry'], 'Original completed Pd geometry'))
    compare_json(g['specification'], run['specification'], 'Consumed geometry-run specification')
    require(('partial_admission' in g) is adopted, 'Geometry envelope/admission kind mismatch')
    require(('labelled_admission' in g) is labelled, 'Geometry envelope/labelled admission mismatch')
    if adopted:
        require(run['geometry']['path'] == 'geometry-partly-relaxed.json',
                'Adopted run must bind its immutable partly-relaxed geometry companion')
    reader._no_extra_electrons(run)
    reader._no_extra_electrons(g)
    epoch, _ = reader._epoch(root, dict(run, geometry=g), required_sources=EPOCH_SOURCES)
    _relaxed_geometry(root, g, g, root, log_inventory, archive_root)
    require(binding(path) == original_binding, 'Geometry envelope changed during replay')
    if labelled:
        from reader_pd_label import validate
        raw = read_json(bound_file(root, g['state_file'], 'Labelled actual state'))
        accounting = validate(root, g, raw, np.asarray(g['positions_A']), np.asarray(g['cell_A']))
        compare_json(run['cpu_accounting'], accounting, 'Separate original/failed actual unit counters')
    elif adopted:
        evidence = read_json(bound_file(root, g['partial_admission']['stopping_evidence'], 'Adopted geometry CPU evidence'))
        counter = run['Native_unit_CPUUsageNSec']
        require(type(counter) is str and re.fullmatch(r'[0-9]+', counter) is not None and
                str(int(counter)) == counter and int(counter) == evidence['CPUUsageNSec'],
                'Adopted geometry lacks its lossless actual terminal unit CPU counter')
        accounting = dict(scope='terminal-manager-unit', CPUUsageNSec=counter,
                          derived_cpu_s=int(counter) / 10**9)
    else:
        cpu = number(run['process_cpu_s'], 'Geometry driver process CPU')
        require(cpu > 0, 'Missing original geometry driver process CPU')
        accounting = dict(scope='geometry-driver-process', process_cpu_s=cpu)
    qualification = _qualification(g)
    return dict(schema='micro-pd-geometry-replay-v1', material='Pd', case=run['case'],
                specification=run['specification'], geometry=dict(g, **qualification), source_epoch=run['source_epoch'],
                source_manifest=epoch, case_binding=original_binding, qualified=False,
                cpu_accounting=accounting, envelope_kind='labelled-original-after-failed-import' if labelled else 'adopted-cap' if adopted else 'original',
                geometry_binding=run['geometry'], **qualification,
                geometry_force_qualified=qualification['force_qualified'], paired_hk_qualified=False,
                full_scientific_qualification=False,
                missing_controls=list(ROLES), BO_available=False, common_ray_selected=False,
                recorded_deadline_utc=run['specification']['absolute_stop_utc'],
                deadline_scope='Configured/consumed stop only; runtime authorization is not inferred',
                scope='Conditionally eligible geometry only, not a BO/contact result. ' + SCOPE)


def replay(directory, *, log_inventory=None, archive_root=None):
    if not (Path(directory) / 'case.json').exists():
        return replay_geometry(directory, log_inventory=log_inventory, archive_root=archive_root)
    if read_json(Path(directory) / 'case.json')['schema'] == 'micro-pd-parallel-case-v1':
        from reader_pd_parallel import replay as replay_parallel
        result = replay_parallel(directory, log_inventory=log_inventory, archive_root=archive_root)
    elif read_json(Path(directory) / 'case.json')['schema'] == 'micro-pd-terminal-recovered-case-v1':
        from reader_pd_terminal import replay as replay_terminal
        result = replay_terminal(directory, log_inventory=log_inventory, archive_root=archive_root)
    else:
        result = reader._replay(directory, log_inventory=log_inventory, archive_root=archive_root,
                                valence=VALENCE, required_sources=EPOCH_SOURCES,
                                validate_spec=_spec, geometry_validator=_geometry)
    qualification = _qualification(result['geometry'])
    _physical_sources(result['source_manifest'])
    result['ray_policy'] = _ray_status(result['geometry'])
    result['geometry'] = dict(result['geometry'], **qualification)
    result.update(schema='micro-pd-replay-v1', material='Pd', scope=SCOPE, **qualification,
                  geometry_force_qualified=qualification['force_qualified'], paired_hk_qualified=False,
                  qualified=False, full_scientific_qualification=False,
                  recorded_deadline_utc=result['specification']['absolute_stop_utc'],
                  deadline_scope='Configured/consumed stop only; runtime authorization is not inferred by this reader')
    return result


def contrast_change(reference, control):
    reference, control = number(reference, 'Base contrast'), number(control, 'Control contrast')
    change = control - reference
    tolerance = max(5., .1 * abs(reference))
    return dict(base_contrast_eV=reference, control_contrast_eV=control, control_minus_base_eV=change,
                absolute_change_eV=abs(change), tolerance_eV=tolerance, resolved=abs(change) <= tolerance)

def _ray_status(g):
    witness = g['coated_witness']
    require(witness['not_final'] is True and witness['force_qualified'] is False and
            witness['bo_eligible'] is False, 'Coated witness falsely admitted')
    require(g['direction_status'] == ('precommitted shared-open' if
            g['specification']['variant'] == 'oxygen' else
            'precommitted ray; final coated shared-open admission pending'),
            'Wrong variant-specific precommit/admission status')
    return dict(authority_entry='entry_2971', precommit=g['probe_ray_precommit'],
                consumed_policy=g['probe_ray'], direction_n_hat=g['direction'],
                clearance=g['ray_clearance'],
                status=g['direction_status'],
                final_coated_admitted=g['direction_status'] == 'precommitted shared-open',
                coated_witness=dict(schema=witness['schema'], not_final=True,
                    force_qualified=False, bo_eligible=False, boundary=witness['boundary']),
                final_frame_optimum_claim=False)



def paired_comparisons(cases, *, allow_partial=False):
    require((cases or allow_partial) and set(cases) <= set(ROLES), 'Unknown/empty Pd role map')
    missing = [role for role in ROLES if role not in cases]
    require(allow_partial or not missing, 'Missing paired base/h/k case files; use --allow-partial only for missing files')
    first = next(iter(cases.values()), None)
    method_keys = ('eigensolver','dataset_bindings','library_bindings','force_tolerance_eV_A','ionic_steps')
    geometry_by_role = {}
    variant_frames = {}
    for role, run in cases.items():
        variant, suffix = role.split('-',1) if '-' in role else (role,'')
        spec = run['specification']
        _spec(spec, run)
        qualification = _qualification(run['geometry'])
        for key, value in qualification.items():
            if key in run:
                compare_json(run[key], value, 'Replayed Pd qualification ' + key)
        geometry_by_role[role] = qualification
        if not qualification['force_qualified']:
            expected = 'entry_3001' if 'labelled_admission' in run['geometry'] else 'entry_2904'
            require(spec.get('geometry_policy_entry') == expected, 'BO policy differs from actual admission')
        require(spec['variant'] == variant and spec['h_A'] == (.16 if suffix == 'h' else .2) and
                spec['kmesh'] == ([2,2,1] if suffix == 'k' else [1,1,1]) and
                spec['radii_A'] == ([.1] if suffix else run.get('completed_radii_A', RADII)), 'Wrong actual Pd named role ' + role)
        compare_json({k:spec[k] for k in method_keys}, {k:first['specification'][k] for k in method_keys}, 'Common Pd physical method')
        compare_json(_physical_sources(run['source_manifest']), _physical_sources(first['source_manifest']),
                     'Identical critical physical source bindings')
        for key in ('bytes', 'sha256'):
            require(run['geometry']['probe_ray_precommit'][key] == first['geometry']['probe_ray_precommit'][key],
                    'Paired cases changed immutable precommit identity')
        close(run['geometry']['cell_A'], first['geometry']['cell_A'], 'Paired identical compact cell')
        require(run['geometry']['direction'] == first['geometry']['direction'], 'Paired exact precommitted ray')
        _ray_status(run['geometry'])
        anchor = variant_frames.setdefault(variant, run)
        for key in ('symbols', 'positions_A', 'cell_A', 'masses_u', 'target_index', 'direction',
                    'fixed_indices', 'free_indices'):
            compare_json(run['geometry'][key], anchor['geometry'][key], 'Variant frozen frame ' + key)
        for key in ('bytes', 'sha256'):
            require(spec['selected_geometry']['binding'][key] == anchor['specification']['selected_geometry']['binding'][key]
                    and spec['probe_ray'][key] == anchor['specification']['probe_ray'][key],
                    'Variant controls changed selected geometry/ray policy identity')
        for state_role, records in run['states'].items():
            state = records[0] if isinstance(records, list) else records
            reference = first['states'][state_role]
            reference = reference[0] if isinstance(reference, list) else reference
            for symbol in state['actual_setups'].keys() & reference['actual_setups'].keys():
                compare_json(state['actual_setups'][symbol], reference['actual_setups'][symbol],
                             'Paired actual setup ' + state_role + ':' + symbol)
        compare_json({k: {f: v[f] for f in ('bytes', 'sha256')} for k, v in run['operator_sources'].items()},
                     {k: {f: v[f] for f in ('bytes', 'sha256')} for k, v in first['operator_sources'].items()},
                     'Paired actual operator sources')
        if suffix and variant in cases:
            base = cases[variant]
            for key in ('symbols','positions_A','cell_A','masses_u','target_index','direction','fixed_indices','free_indices'):
                compare_json(run['geometry'][key], base['geometry'][key], 'Frozen-frame ' + role + ':' + key)
            for key in ('bytes','sha256'):
                require(spec['selected_geometry']['binding'][key] == base['specification']['selected_geometry']['binding'][key],
                        'Control changed selected original geometry')
            require(all(spec['probe_ray'][key] == base['specification']['probe_ray'][key]
                        for key in ('bytes', 'sha256')), 'Variant controls changed original ray policy')
            if suffix == 'h':
                coarse, fine = base['states']['occupied-host']['realspace_shape'], run['states']['occupied-host']['realspace_shape']
                require(all(f >= c for f,c in zip(fine,coarse)) and any(f > c for f,c in zip(fine,coarse)), 'Actual h grid not refined')
    pairs = {}
    for suffix in ('','-h','-k'):
        if not all(v+suffix in cases for v in ('plain','oxygen')):
            continue
        plain, oxygen = (cases[v+suffix] for v in ('plain','oxygen'))
        point = {}
        for v, run in (('plain',plain),('oxygen',oxygen)):
            matches = [row for row in run['rows'] if row['R_A'] == .1]
            require(len(matches) == 1, 'Missing actual Pd R0.1 row')
            point[v] = matches[0]['lowering_eV']
        pairs[suffix or 'base'] = dict(plain_R01_lowering_eV=point['plain'], oxygen_R01_lowering_eV=point['oxygen'],
            oxygen_minus_plain_R01_lowering_eV=point['oxygen']-point['plain'],
            plain_signedW0_eV=plain['signedW0_eV'], oxygen_signedW0_eV=oxygen['signedW0_eV'],
            oxygen_minus_plain_signedW0_eV=oxygen['signedW0_eV']-plain['signedW0_eV'],
            geometries={v: geometry_by_role[v+suffix] for v in ('plain', 'oxygen')},
            geometry_error_bound_eV=None,
            geometry_error_status='Unbounded; actual forces do not bound contrast error')
    changes = {}
    for key in ('h','k'):
        if 'base' in pairs and '-'+key in pairs:
            changes[key] = {name: contrast_change(pairs['base']['oxygen_minus_plain_'+name+'_eV'],
                                                 pairs['-'+key]['oxygen_minus_plain_'+name+'_eV'])
                            for name in ('R01_lowering','signedW0')}
    resolved = not missing and len(changes) == 2 and all(x['resolved'] for group in changes.values() for x in group.values())
    geometry_force_qualified = all(v in cases for v in ('plain', 'oxygen')) and all(
        value['force_qualified'] for value in geometry_by_role.values())
    return dict(schema='micro-pd-paired-comparisons-v1', material='Pd', pairs=pairs, changes=changes,
                qualified=bool(resolved and geometry_force_qualified), paired_hk_qualified=bool(resolved),
                geometry_force_qualified=geometry_force_qualified, geometry_by_role=geometry_by_role,
                ray_policy_entry='entry_2971',
                ray_by_role={role: _ray_status(run['geometry']) for role, run in cases.items()},
                shared_precommit=({key: first['geometry']['probe_ray_precommit'][key]
                                   for key in ('bytes', 'sha256')} if first else None),
                final_coated_admitted=any(run['specification']['variant'] == 'oxygen' for run in cases.values()),
                final_frame_optimum_claim=False,
                geometry_error_bound_eV=None, full_scientific_qualification=False,
                resolved=bool(resolved), complete=not missing, missing_controls=missing,
                recorded_deadlines_utc={role: run['specification']['absolute_stop_utc'] for role, run in cases.items()},
                deadline_scope='Configured/consumed stops only; no runtime authorization inferred',
                qualification_scope='Separate paired R0.1/contact h/k criterion and geometry force stationarity; '
                                    'h/k checks do not qualify geometry or bound geometry-convergence error',
                tolerance_definition='max(5 eV, 0.1*abs(base coated-minus-bare contrast)) per observable',
                omitted_controls=['geometry-convergence error bound','slab-thickness','cell/vacuum','full-radius h/k'],
                paired_cell_qualified=False, full_curve_higher_k_qualified=False, isolation_claim=False, scope=SCOPE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', nargs='?')
    parser.add_argument('--case', action='append', metavar='ROLE=DIRECTORY')
    parser.add_argument('--allow-partial', action='store_true', help='Allow missing case files only; never waive geometry admission')
    args = parser.parse_args()
    if args.case:
        require(args.directory is None, 'Use one case directory or a paired role map')
        roots = {}
        for item in args.case:
            role, path = item.split('=',1)
            require(role not in roots, 'Duplicate Pd role')
            roots[role] = path
        cases = {role:replay(path) for role,path in roots.items()}
        result = dict(cases=cases, comparisons=paired_comparisons(cases, allow_partial=args.allow_partial))
    else:
        require(args.directory is not None and not args.allow_partial, 'One original case root required')
        result = replay(args.directory)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
