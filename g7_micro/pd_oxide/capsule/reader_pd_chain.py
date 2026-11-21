"""Independent saved-evidence replay of Pd best-boundary and one-resume chains.

No calculator imports: nonselected boundaries are evidence, never BO admissions.
"""
import math
from pathlib import Path
import re

import numpy as np

import reader
from reader import bound_file, close, compare_json, number, read_json, require

SELECTION_AUTHORITY = dict(bytes=4677, sha256='283e4bb2f45ffb64a7026aa7a243acf9047aad809a421621f8e76ab920eef664')


def _pd():
    import reader_pd
    return reader_pd


def _same_binding(left, right, label):
    require(all(left[k] == right[k] for k in ('bytes', 'sha256')), label)


def _files(root, items, label):
    require(isinstance(items, list) and items, 'Missing ' + label)
    names = set()
    for item in items:
        require(item['path'] not in names, 'Duplicate ' + label)
        names.add(item['path'])
        bound_file(root, item, label)


def epoch(path, spec):
    pd = _pd()
    value = read_json(path)
    require(value['schema'] == 'micro-frozen-repository-source-v1', 'Wrong original source epoch')
    _same_binding(reader.binding(path), spec['frozen_source_manifest'],
                  'Original input/source epoch mismatch')
    _files(path.parent / 'source', value['files'], 'original frozen source')
    require(all(sum(Path(row['path']).name == name for row in value['files']) == 1
                for name in pd.EPOCH_SOURCES), 'Incomplete original source closure')
    pd._physical_sources(value)
    return value


def cap(root, evidence, spec, *, command_key='manager_command'):
    """Replay actual terminal accounting, never interpret an arbitrary SCF failure as a cap."""
    fields = ('manager_accounting', 'native_log', 'original_input', 'original_source_epoch',
              'launch_receipt', 'manager_bounds', command_key)
    paths = {key: bound_file(root, evidence[key], 'Capped original ' + key) for key in fields}
    compare_json(read_json(paths['original_input']), spec, 'Original capped input')
    require(spec['phase'] == 'geometry' and spec['material'] == 'pd' and spec['size'] == 'compact' and
            spec['variant'] in ('plain', 'oxygen') and spec['initial_geometry'] is None and
            spec['selected_geometry'] is None and spec['force_tolerance_eV_A'] == .05 and
            spec['h_A'] == .2 and spec['kmesh'] == [1, 1, 1], 'Wrong capped geometry family')
    _pd()._recorded_stop(spec)
    limit = number(spec['cpu_limit_s'], 'Actual geometry CPU cap')
    require(0 < limit <= (2400 if spec['variant'] == 'plain' else 6000), 'Changed geometry cap')
    epoch(paths['original_source_epoch'], spec)
    def manager(path):
        result = {}
        for line in path.read_text().splitlines():
            key, sep, value = line.partition('=')
            require(sep and key not in result, 'Malformed/duplicate manager field')
            result[key] = value
        return result
    accounting = manager(paths['manager_accounting'])
    require(accounting.get('Result') in ('timeout', 'exit-code') and
            re.fullmatch(r'[0-9]+', accounting.get('CPUUsageNSec', '')) is not None and
            re.fullmatch(r'[0-9a-fA-F]{32}', accounting.get('InvocationID', '')) is not None and
            bool(accounting.get('MemoryPeak')), 'Missing actual closed manager receipt')
    launched = read_json(paths['launch_receipt'])
    require(launched['case'] == spec['case'] and launched['phase'] == 'geometry' and
            launched['unit'] == 'coherence-' + launched['name'] + '.service', 'Wrong capped launch')
    _same_binding(launched['input'], evidence['original_input'], 'Capped launch input mismatch')
    command = read_json(paths[command_key])
    require(isinstance(command, list) and command.count('systemctl') == 1, 'Wrong manager query')
    start = command.index('systemctl')
    require(command[start:start+4] == ['systemctl', '--user', 'show', launched['unit']] and
            'InvocationID' in command[start+4:], 'Manager query names another invocation')
    observed = manager(paths['manager_bounds'])
    require(observed.get('InvocationID') == accounting['InvocationID'] and
            observed.get('CPUAffinity') == str(launched['cpu']) and observed.get('Restart') == 'no',
            'Cap belongs to another manager invocation')
    text = paths['native_log'].read_text()
    require(text.startswith('MICRO_PD_START case=' + spec['case'] + ' phase=geometry\n') and
            'MICRO_PD_COMPLETE' not in text, 'Not the closed capped geometry log')
    cpu = int(accounting['CPUUsageNSec'])
    if accounting['Result'] == 'exit-code':
        require(text.rstrip().endswith('ValueError: Registered whole-process CPU cap reached; no automatic retry') and
                'Traceback (most recent call last):' in text and cpu >= limit * 10**9,
                'Failure is not the actual process CPU cap')
    actual = dict(result=accounting['Result'], CPUUsageNSec=cpu, cpu_limit_s=spec['cpu_limit_s'],
                  InvocationID=accounting['InvocationID'], MemoryPeak=accounting['MemoryPeak'], unit=launched['unit'])
    if evidence.get('schema') == 'micro-pd-cap-evidence-v1':
        require(type(evidence['CPUUsageNSec']) is int,
                'Terminal unit CPU clock must remain a lossless integer')
        for key, value in actual.items():
            compare_json(evidence[key], value, 'Actual cap ' + key)
    return paths, text, actual


def boundary(record, spec, initial):
    """Check every candidate's actual arrays and physics without requiring its GPW copy."""
    pd = _pd()
    compare_json(record['specification'], spec, 'Accepted boundary frozen specification')
    frame = dict(record, target_index=initial['target_index'])
    positions, cell, masses = pd._frame(frame, spec['variant'])
    ip, ic, im = pd._frame(initial, spec['variant'])
    close(cell, ic, 'Accepted original cell')
    close(masses, im, 'Accepted original isotope inventory')
    close(positions[pd.FIXED], ip[pd.FIXED], 'Accepted fixed bottom positions')
    require(initial['fixed_indices'] == pd.FIXED and initial['free_indices'] == record['free_indices'],
            'Accepted bottom-only constraints changed')
    reader._no_extra_electrons(record)
    require(record['schema'] == 'micro-native-state-v1' and record['backend'] == 'gpaw-native' and
            record['scf_converged'] is True and record['geometry_PAW'] is True and
            type(record['charge']) is int and record['charge'] == 0 and record['xc'] == 'PBE' and
            record['spinpol'] is False, 'Not an accepted neutral Native PAW boundary')
    ne = sum(pd.VALENCE[s] for s in record['symbols'])
    close(record['electrons'], ne, 'Accepted actual Ne')
    close(record['occupation_width_eV'], .1, 'Accepted smearing')
    lengths = 1 / np.linalg.norm(np.linalg.inv(cell), axis=0)
    close(record['realspace_shape'], np.maximum(4, (lengths / spec['h_A'] / 4 + .5).astype(int) * 4),
          'Accepted actual grid')
    require(record['actual_kpoints'] == math.prod(spec['kmesh']) and record['nbands'] == math.ceil(ne / 2) + 8 and
            type(record['scf_iterations']) is int and 0 < record['scf_iterations'] <= 150 and
            spec['eigensolver'] == {'name': 'davidson', 'niter': 2} and record['eigensolver'] == 'Davidson',
            'Accepted SCF method changed')
    require(set(record['actual_setups']) == set(record['symbols']), 'Missing accepted PAW setup metadata')
    for symbol, setup in record['actual_setups'].items():
        close(setup['Nv'], pd.VALENCE[symbol], 'Accepted PAW valence')
        require(type(setup['ni']) is int and setup['ni'] > 0 and setup['Nc'] >= 0,
                'Accepted PAW projector family changed')
    energy = record['energy']
    require(set(energy) == {'free_energy_eV', 'internal_energy_eV', 'minus_TS_eV'}, 'Accepted energy fields')
    for key, value in energy.items():
        number(value, 'Accepted ' + key)
    close(energy['internal_energy_eV'], energy['free_energy_eV'] - energy['minus_TS_eV'], 'Accepted energy identity')
    require(energy['minus_TS_eV'] <= 1e-9 and type(record['accepted_step']) is int and
            record['accepted_step'] >= 0 and number(record['process_cpu_s'], 'Accepted process CPU') > 0,
            'Invalid accepted energy/step/CPU')
    return pd._forces(record, positions.shape, geometric_witness=True)['free_fmax_eV_A']


def candidates(root, descriptors, text, spec, initial, *, first=0):
    _files(root, descriptors, 'accepted boundary inventory')
    records = [read_json(bound_file(root, row, 'Accepted boundary')) for row in descriptors]
    marks = re.findall(r'^MICRO_PD_GEOMETRY_STEP accepted_step=(\d+) free_fmax_eV_A=(\S+) '
                       r'raw_fmax_eV_A=(\S+) process_cpu_s=(\S+)\s*$', text, re.MULTILINE)
    indices = list(range(first, first + len(records)))
    require([row['accepted_step'] for row in records] == indices and
            [int(row[0]) for row in marks] == indices, 'Closed log/inventory has missing, extra or reordered boundaries')
    require(type(spec['ionic_steps']) is int and 0 < spec['ionic_steps'] <= 60 and
            records[-1]['accepted_step'] <= (spec['ionic_steps'] if first == 0 else first - 1 + spec['ionic_steps']),
            'Accepted inventory exceeds the bounded optimizer steps')
    forces = []
    for record, mark in zip(records, marks):
        forces.append(boundary(record, spec, initial))
        compare_json(record['actual_setups'], records[0]['actual_setups'], 'Unchanged accepted PAW setups')
        for actual, key in zip(mark[1:], ('free_fmax_eV_A', 'raw_fmax_eV_A', 'process_cpu_s')):
            close(float(actual), record[key], 'Accepted logged ' + key)
    return records, forces


def resume(root, geometry):
    """Return the immutable original chain after independently checking the one charged resume."""
    pd = _pd()
    spec = geometry['specification']
    item = geometry.get('checkpoint_resume')
    require(isinstance(item, dict) and isinstance(spec.get('checkpoint_resume'), dict), 'Unbound checkpoint resume')
    _same_binding(item, spec['checkpoint_resume'], 'Resume manifest differs from consumed fresh input')
    path = bound_file(root, item, 'Archived checkpoint resume manifest')
    manifest = read_json(path)
    base = path.parent
    require(manifest['schema'] == 'micro-pd-checkpoint-resume-v1' and manifest['generation'] == 1 and
            type(manifest['generation']) is int and manifest['authority_entry'] == 'entry_2980', 'Renewed/wrong resume authority')
    keys = ('original_input', 'original_initial_frame', 'original_source_epoch', 'original_boundary',
            'original_gpw', 'manager_accounting', 'native_log', 'launch_receipt', 'manager_bounds',
            'manager_query_argv', 'authority')
    paths = {key: bound_file(base, manifest[key], 'Resume original ' + key) for key in keys}
    require(re.search(r'\bentry_2980\b', paths['authority'].read_text()) is not None, 'Wrong original resume authority')
    _same_binding(manifest['authority'], SELECTION_AUTHORITY, 'Resume requires exact captured entry_2980 original')
    old = read_json(paths['original_input'])
    require(old.get('checkpoint_resume') is None and old['variant'] == 'oxygen' and
            old['ionic_steps'] == 60 and 5880 <= old['cpu_limit_s'] <= 6000 and spec['variant'] == 'oxygen' and
            type(spec['ionic_steps']) is int and 0 < spec['ionic_steps'] <= 60 and
            spec['initial_geometry'] is None and spec['selected_geometry'] is None,
            'Resume is not one bounded coated recovery from original source-built input')
    mutable = {'case', 'cpu_limit_s', 'ionic_steps', 'frozen_source_manifest', 'checkpoint_resume'}
    compare_json({k: v for k, v in spec.items() if k not in mutable},
                 {k: v for k, v in old.items() if k not in mutable}, 'Resume changed commissioned physics or bounds')
    pd._recorded_stop(spec)
    require(0 < number(spec['cpu_limit_s'], 'Charged resume CPU cap') <= 3600,
            'Resume exceeds one charged sixty-minute allowance')
    initial = read_json(paths['original_initial_frame'])
    require(initial['schema'] == 'micro-pd-initial-frame-v1', 'Wrong original initial frame')
    _same_binding(geometry['initial_frame'], manifest['original_initial_frame'], 'Resume reset migration reference')
    _same_binding(geometry['source_geometry'], manifest['original_boundary'], 'Resume source is not original last boundary')
    _, text, receipt = cap(base, manifest, old, command_key='manager_query_argv')
    compare_json(manifest['cap_receipt'], receipt, 'Resume actual original cap receipt')
    require(type(manifest['cap_receipt']['CPUUsageNSec']) is int and
            read_json(paths['launch_receipt'])['core_min_cap'] == 100,
            'Resume requires the actual original 100 core-minute cap and integer CPU clock')
    records, forces = candidates(base, manifest['accepted_boundaries'], text, old, initial)
    require(len(records) >= 3 and all(value > .30 for value in forces), 'Eligible original frame forbids purchasing resume')
    last = records[-1]
    require(type(manifest['original_last_accepted_step']) is int and
            manifest['original_last_accepted_step'] == last['accepted_step'], 'Resume original global step mismatch')
    compare_json(read_json(paths['original_boundary']), last, 'Resume must start from last completed boundary')
    _same_binding(manifest['original_boundary'], manifest['accepted_boundaries'][-1], 'Wrong last source-bound boundary')
    _same_binding(manifest['original_gpw'], last['gpw'], 'Original checkpoint binding changed')
    gate = manifest['gate']
    require(gate['energy_field'] == 'free_energy_eV', 'Wrong actual ASE force-consistent optimizer energy convention')
    energies = [row['energy'][gate['energy_field']] for row in records]
    improving_energy = all(right < left for left, right in zip(energies, energies[1:]))
    recent, earlier = min(forces[-2:]), min(forces[:-2])
    expected = dict(energy_field='free_energy_eV', energy_series_eV=energies,
                    free_fmax_series_eV_A=forces, energy_stepwise_decreasing=improving_energy,
                    last_two_min_fmax=recent, earlier_min_fmax=earlier, force_improving=recent < earlier)
    for key, value in expected.items():
        compare_json(gate[key], value, 'Actual original-cap gate ' + key)
    require(improving_energy and recent < earlier, 'Original actual-cap recovery gate failed')
    _files(base, manifest['source_files'], 'resume original full source closure')
    original_epoch = read_json(paths['original_source_epoch'])
    expected_sources = [dict(path=(paths['original_source_epoch'].parent / 'source' / row['path']).relative_to(base).as_posix(),
                             bytes=row['bytes'], sha256=row['sha256']) for row in original_epoch['files']]
    compare_json(manifest['source_files'], expected_sources, 'Complete original resume source closure')
    _files(base, manifest['evidence_files'], 'resume immutable evidence closure')
    energy_items = [row for row in manifest['evidence_files']
                    if Path(row['path']).parent.name == 'energy-convention']
    energy_sources = {Path(row['path']).name: bound_file(base, row, 'Original optimizer convention')
                      for row in energy_items}
    require(set(energy_sources) == {'bfgslinesearch.py', 'optimize.py'} and len(energy_items) == 2,
            'Missing exact qualified ASE energy dispatch evidence')
    optimize = energy_sources['optimize.py'].read_text()
    require('force_consistent=True' in optimize and 'force_consistent=force_consistent' in optimize,
            'Archived ASE does not establish force-consistent free-energy convention')
    current_epoch = read_json(bound_file(root, dict(path='provenance/source-manifest.json',
                               **{k: spec['frozen_source_manifest'][k] for k in ('bytes', 'sha256')}), 'Resume current epoch'))
    compare_json(pd._physical_sources(current_epoch), pd._physical_sources(original_epoch), 'Resume changed physical source')
    return manifest, records, forces, initial


def partial(root, geometry, raw, positions, cell):
    pd = _pd()
    admission = geometry['partial_admission']
    require(admission['schema'] == 'micro-pd-partial-admission-v2' and
            admission['authority_entry'] == 'entry_2904' and admission['selection_authority_entry'] == 'entry_2980' and
            admission['maximum_free_fmax_eV_A'] == .30, 'Wrong best-accepted partial authority')
    keys = ('source_boundary', 'selected_boundary', 'stopping_evidence', 'authority', 'selection_authority', 'adoption_source')
    paths = {key: bound_file(root, admission[key], 'Partial ' + key) for key in keys}
    for key, entry in (('authority', 'entry_2904'), ('selection_authority', 'entry_2980')):
        require(re.search(r'\b' + entry + r'\b', paths[key].read_text()) is not None, 'Wrong bound selection authority')
    _same_binding(admission['selection_authority'], SELECTION_AUTHORITY, 'Selection requires exact captured entry_2980 original')
    require(paths['adoption_source'].name == 'pd_adopt.py', 'Wrong bound adopter')
    evidence = read_json(paths['stopping_evidence'])
    require(evidence['schema'] == 'micro-pd-cap-evidence-v1', 'Wrong cap evidence schema')
    spec = geometry['specification']
    _, text, _ = cap(root, evidence, spec)
    initial = read_json(bound_file(root, geometry['initial_frame'], 'Partial original initial frame'))
    if spec.get('checkpoint_resume') is not None:
        records, forces = resumed_boundaries(root, geometry)
        compare_json(admission['accepted_boundaries'], geometry['accepted_boundaries'],
                     'Partial/resume complete newly accepted inventory')
        _same_binding(evidence['native_log'], geometry['native_log'], 'Partial/resume closed Native log')
    else:
        records, forces = candidates(root, admission['accepted_boundaries'], text, spec, initial)
    eligible = [(value, index) for index, value in enumerate(forces) if value <= .30]
    require(eligible, 'No eligible accepted geometry exists')
    _, index = min(eligible)
    selected = records[index]
    require(type(admission['selected_accepted_step']) is int and
            admission['selected_accepted_step'] == selected['accepted_step'] == geometry['accepted_steps'],
            'Partial selection is not lowest eligible actual free force (first exact tie)')
    compare_json(admission['selected_boundary'], admission['accepted_boundaries'][index], 'Selected inventory descriptor')
    _same_binding(admission['source_boundary'], admission['selected_boundary'], 'Changed selected source boundary')
    compare_json(read_json(paths['selected_boundary']), selected, 'Selected original accepted boundary')
    added = {'log', 'sanitized_log', 'log_privacy_transforms', 'force_qualified', 'relaxation_status', 'geometry_error_bound_eV'}
    compare_json({key: value for key, value in raw.items() if key not in added}, selected, 'Unaltered selected accepted PAW')
    pd._forces(selected, positions.shape)
    site = admission['site']
    require(isinstance(site['description'], str) and site['description'].strip() and
            site['criterion'] == 'nonempty Pd/O coordination within 3.5 A and no nuclear coincidence', 'Missing partial site criterion')
    for key, value in geometry['coordination'].items():
        compare_json(site[key], value, 'Partial actual coordination ' + key)
    actual = pd._neighbors(positions, cell, geometry['target_index'], 3.5)
    require(actual and min(actual.values()) > 1e-6 and any(geometry['symbols'][i] in ('Pd', 'O') for i, _ in actual),
            'Partial target lacks definable metal/oxide site')


def resumed_boundaries(root, geometry):
    """Replay import evidence separately from newly accepted SCF boundaries."""
    manifest, old_records, _, initial = resume(root, geometry)
    spec = geometry['specification']
    proof = read_json(bound_file(root, geometry['checkpoint_import'], 'Original cached checkpoint import'))
    last = old_records[-1]
    require(proof['schema'] == 'micro-pd-checkpoint-import-v1' and
            proof['cached_import'] is True and proof['new_scf'] is False and
            type(proof['original_accepted_step']) is int and
            proof['original_accepted_step'] == last['accepted_step'] and
            proof['optimizer_history'] == 'reset; no saved Hessian' and
            proof['reuse'] == 'public GPAW.read; identical atoms/electronic parameters; original cached evidence only',
            'Imported checkpoint falsely reported as new Native work or retained optimizer history')
    _same_binding(proof['checkpoint_resume'], geometry['checkpoint_resume'], 'Changed imported manifest')
    _same_binding(proof['source_boundary'], manifest['original_boundary'], 'Changed imported boundary')
    compare_json(proof['energy'], last['energy'], 'Cached import original energy')
    compare_json(proof['raw_forces_eV_A'], last['raw_forces_eV_A'], 'Cached import original forces')
    text = bound_file(root, geometry['native_log'], 'Immutable resumed accepted-boundary log').read_text()
    require(text.startswith('MICRO_PD_START case=' + spec['case'] + ' phase=geometry\n') and
            re.findall(r'^MICRO_PD_CHECKPOINT_IMPORT accepted_step=(\d+)\s*$', text, re.MULTILINE) ==
            [str(last['accepted_step'])], 'Missing/duplicate original checkpoint import mark')
    rows, forces = candidates(root, geometry['accepted_boundaries'], text, spec, initial,
                              first=last['accepted_step'] + 1)
    require(rows[-1]['accepted_step'] <= last['accepted_step'] + spec['ionic_steps'],
            'Resume exceeded fresh bounded optimizer step allowance')
    for row in rows:
        compare_json(row['actual_setups'], last['actual_setups'], 'Resumed PAW setup continuity')
    selected = [row for row in rows if row['accepted_step'] == geometry['accepted_steps']]
    require(len(selected) == 1, 'Selected geometry is not a newly accepted resume boundary')
    row = selected[0]
    for key in ('symbols', 'positions_A', 'cell_A', 'masses_u', 'raw_forces_eV_A', 'constrained_forces_eV_A'):
        compare_json(geometry[key], row[key], 'Resumed selected actual ' + key)
    _same_binding(geometry['gpw'], row['gpw'], 'Resumed selected checkpoint changed')
    if 'partial_admission' not in geometry:
        require(geometry['native_log_status'] == 'accepted-boundary prefix; unit log not yet closed',
                'Full resume must not misrepresent its in-driver log prefix as a closed unit')
        require(row is rows[-1] and forces[-1] <= .05, 'False full resumed geometry')
    return rows, forces
