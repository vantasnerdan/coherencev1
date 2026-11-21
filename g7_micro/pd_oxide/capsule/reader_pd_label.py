"""Independent replay of authority-labelled original frames after a failed import."""
from pathlib import Path
import re

from reader import binding, bound_file, compare_json, read_json, require
import reader_pd_chain as chain

AUTHORITY = dict(bytes=7354, sha256='14b58a2ebb1d0bf6bcec5ba311fa4077675cc4ac45274ddd7859c5bf9d5ce1a5')
SCHEMA = 'micro-pd-labelled-frame-admission-v1'
FRAME_AUTHORITY = dict(bytes=4511, sha256='08c8f29684910ca7bb212d072949aa398ff3de9cedbd55ee07765fdd9806ed6a')


def manager(path):
    result = {}
    for line in path.read_text().splitlines():
        key, sep, value = line.partition('=')
        require(sep and key not in result, 'Malformed/duplicate terminal manager field')
        result[key] = value
    return result


def validate(root, geometry, raw, positions, cell):
    import reader_pd as pd
    root = Path(root)
    admission = geometry['labelled_admission']
    require(admission['schema'] == SCHEMA and admission['authority_entry'] == 'entry_2994' and
            admission['geometry_error_bound_eV'] is None and admission['force_target_eV_A'] == .05 and
            admission['new_resumed_steps'] == 0 and admission['checkpoint_import_qualified'] is False,
            'Wrong labelled-frame authority or invented resumed work')
    authority = bound_file(root, admission['authority'], 'Exact labelled authority')
    chain._same_binding(binding(authority), AUTHORITY, 'Labelled admission requires exact original entry_2994')
    require('entry_2994' in authority.read_text(), 'Wrong captured labelled authority')
    producer = bound_file(root, admission['adoption_source'], 'Labelled producer')
    confirmation = bound_file(root, admission['frame_authority'], 'Specific original-frame authorization')
    chain._same_binding(binding(confirmation), FRAME_AUTHORITY, 'Require exact original entry_3001 confirmation')
    require(admission['frame_authority_entry'] == 'entry_3001' and 'entry_3001' in confirmation.read_text(),
            'Missing specific original-frame authority')
    require(producer.name == 'pd_label.py', 'Wrong labelled producer')
    evidence = read_json(bound_file(root, admission['stopping_evidence'], 'Original terminal cap'))
    spec = geometry['specification']
    require(spec.get('checkpoint_resume') is None and spec['variant'] == 'oxygen', 'Labelled frame must retain original Native spec')
    _, text, original_counter = chain.cap(root, evidence, spec)
    initial = read_json(bound_file(root, geometry['initial_frame'], 'Original source initial frame'))
    records, forces = chain.candidates(root, admission['accepted_boundaries'], text, spec, initial)
    if admission['selection_kind'] == 'lowest-original-free-force':
        index = min(range(len(forces)), key=forces.__getitem__)
        require(records[index]['accepted_step'] == 5, 'Authority base is original BEST step 5')
    else:
        require(admission['selection_kind'] == 'second-original-frame-sensitivity' and
                admission['selected_accepted_step'] in (3, 6), 'Unauthorized second original frame')
        index = next(i for i, row in enumerate(records) if row['accepted_step'] == admission['selected_accepted_step'])
    selected = records[index]
    compare_json(admission['selected_boundary'], admission['accepted_boundaries'][index], 'Actual selected original boundary')
    compare_json(admission['source_boundary'], admission['selected_boundary'], 'Original source boundary')
    compare_json(read_json(bound_file(root, admission['source_boundary'], 'Selected original')), selected, 'Selected original bytes')
    require(admission['selected_accepted_step'] == geometry['accepted_steps'] == selected['accepted_step'], 'Changed original selected step')
    added = {'log', 'sanitized_log', 'log_privacy_transforms', 'force_qualified', 'relaxation_status', 'geometry_error_bound_eV'}
    compare_json({k: v for k, v in raw.items() if k not in added}, selected, 'Unchanged original Native state/spec/clock')
    pd._forces(selected, positions.shape, labelled=True)
    chain._same_binding(geometry['gpw'], selected['gpw'], 'Changed selected original GPW')
    failed = admission['failed_attempt']
    chain._files(root, failed['files'], 'Closed failed attempt closure')
    base = (root / failed['root']).resolve()
    require(base.is_relative_to(root.resolve()), 'Escaping failed-attempt root')
    fresh = read_json(base / 'provenance/input.json')
    epoch = chain.epoch(base / 'provenance/source-manifest.json', fresh)
    compare_json(pd._physical_sources(epoch), pd._physical_sources(read_json(bound_file(root, evidence['original_source_epoch'], 'Original epoch'))), 'Failed attempt physical source continuity')
    resume_path = base / 'checkpoint-resume/manifest.json'
    resume_manifest = read_json(resume_path)
    synthetic = dict(specification=fresh, checkpoint_resume=dict(path='checkpoint-resume/manifest.json', **{k: binding(resume_path)[k] for k in ('bytes', 'sha256')}),
                     initial_frame=resume_manifest['original_initial_frame'], source_geometry=resume_manifest['original_boundary'])
    manifest, old_records, _, _ = chain.resume(base, synthetic)
    compare_json(old_records, records, 'Failed attempt references same original accepted inventory')
    unit = base.parent
    log = (unit / 'native.log').read_text()
    require(log.startswith('MICRO_PD_START case=' + fresh['case'] + ' phase=geometry\n') and
            log.rstrip().endswith('ValueError: Public checkpoint import changed guarded observers/electronic parameters') and
            'MICRO_PD_CHECKPOINT_IMPORT' not in log and 'MICRO_PD_GEOMETRY_STEP' not in log and
            'MICRO_PD_COMPLETE' not in log and not (base / 'checkpoint-import.json').exists() and
            not list(base.glob('geometry-step-*-boundary.json')), 'Not the actual zero-new-step startup failure')
    accounting = manager(bound_file(root, failed['manager_accounting'], 'Failed unit actual accounting'))
    launch = read_json(unit / 'launch-receipt.json')
    observed = manager(unit / '05-manager-bounds.stdout')
    command = read_json(unit / '05-manager-bounds.argv.json')
    require(accounting['Result'] == 'exit-code' and re.fullmatch(r'[0-9]+', accounting['CPUUsageNSec']) and
            re.fullmatch(r'[0-9a-fA-F]{32}', accounting['InvocationID']) and accounting['MemoryPeak'] and
            observed['InvocationID'] == accounting['InvocationID'] and observed['Restart'] == 'no' and
            observed['CPUAffinity'] == str(launch['cpu']), 'Unbound failed-unit terminal accounting')
    start = command.index('systemctl')
    require(command[start:start+4] == ['systemctl', '--user', 'show', launch['unit']] and
            'InvocationID' in command[start+4:] and launch['case'] == fresh['case'] and launch['phase'] == 'geometry', 'Wrong failed-unit invocation')
    chain._same_binding(launch['input'], binding(base / 'provenance/input.json'), 'Failed launch consumed input')
    issued = read_json(bound_file(root, failed['issued'], 'Original actual issuance flag'))
    attempt = read_json(bound_file(root, failed['attempt'], 'Original actual attempt flag'))
    require(issued == dict(generation=1, manifest=fresh['checkpoint_resume']) and attempt['generation'] == 1 and
            attempt['manifest'] == issued['manifest'] and isinstance(attempt['calculation'], str) and
            Path(attempt['calculation']).is_absolute(), 'Not the single actual consumed resume issuance')
    compare_json(geometry['coordination'], {k: admission['site'][k] for k in geometry['coordination']}, 'Labelled actual D site')
    require(admission['site']['description'].strip() and admission['site']['criterion'] ==
            'nonempty Pd/O coordination within 3.5 A and no nuclear coincidence', 'Missing labelled D-site description')
    neighbors = pd._neighbors(positions, cell, geometry['target_index'], 3.5)
    require(neighbors and min(neighbors.values()) > 1e-6 and any(geometry['symbols'][i] in ('Pd', 'O') for i, _ in neighbors), 'Undefined labelled actual D site')
    return dict(scope='separate terminal manager units; not adopter process CPU',
                original_CPUUsageNSec=str(original_counter['CPUUsageNSec']),
                failed_resume_CPUUsageNSec=str(int(accounting['CPUUsageNSec'])),
                new_resumed_steps=0, checkpoint_import_qualified=False)
