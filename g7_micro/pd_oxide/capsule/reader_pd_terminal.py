"""NumPy-only replay of a terminal OOM plus one genuinely fresh H-R0 execution.

The recovered container is a derivation, never a successful original Native case.
Its two descriptors bind complete ``original/`` and ``recovery/`` unit trees;
state origins and all raw state bindings stay relative to each calculation root.
Only the H-R0 state belongs to recovery. Original process CPU is unavailable:
whole-unit counters and the successful recovery process clock remain separate.
"""
from pathlib import Path
import math
import re

import numpy as np
import reader
from reader import binding, bound_file, close, compare_json, number, read_json, require
from reader_pd_label import manager
from reader_pd_parallel import HOSTS, METHOD, _contact

SCHEMA = 'micro-pd-terminal-recovered-case-v1'
PROCESS_SCOPE = ('successful recovery process only; original killed process CPU unavailable; '
                 'original whole-unit CPU retained separately')
UNIT_SCOPE = 'Exact original whole unit; not the collector or an adoption clock'
UNIT_FILES = {'input.json', 'source-manifest.json', 'launch-config.json', 'launch-receipt.json',
              'geometry-admission.json', '05-manager-bounds.stdout',
              '05-manager-bounds.argv.json', 'native.log'}
STATE_ROLES = {*HOSTS, 'occupied-probes', 'vacant-probes', 'occupied-probe-He', 'vacant-probe-H-R0'}


def _component(root, descriptor, identity):
    require(descriptor['id'] == identity and descriptor['root'] == identity + '/calculation',
            'Wrong terminal component root/identity')
    base = (root / descriptor['root']).resolve()
    tree = root / identity
    require(base.is_relative_to(root) and tree.resolve().is_relative_to(root), 'Unsafe terminal root')
    files = descriptor['files']
    require(isinstance(files, list) and files and len({item['path'] for item in files}) == len(files),
            'Missing/duplicate terminal archive binding')
    for item in files:
        path = bound_file(root, item, 'Complete immutable terminal component')
        require(path.is_relative_to(tree.resolve()), 'Component binding escapes its original tree')
    require({item['path'] for item in files} ==
            {path.relative_to(root).as_posix() for path in tree.rglob('*') if path.is_file()},
            'Incomplete terminal whole-unit closure')
    input_path = bound_file(root, descriptor['input'], 'Actual terminal input')
    epoch_path = bound_file(root, descriptor['source_epoch'], 'Actual terminal source epoch')
    require(input_path == base / 'provenance/input.json' and
            epoch_path == base / 'provenance/source-manifest.json', 'Changed terminal provenance location')
    spec = read_json(input_path)
    execution = dict(case=spec['case'], specification=spec, geometry=descriptor['geometry'],
                     provenance=descriptor['provenance'],
                     source_epoch=dict(descriptor['source_epoch'], path='provenance/source-manifest.json'))
    return base, execution


def _unit(root, descriptor, execution, identity):
    receipt = descriptor['unit_receipt']
    require(receipt['counter_scope'] == UNIT_SCOPE, 'Wrong original whole-unit counter scope')
    fields = manager(bound_file(root, receipt['accounting'], 'Actual terminal unit accounting'))
    for key in ('Result', 'CPUUsageNSec', 'InvocationID', 'MemoryPeak'):
        compare_json(receipt[key], fields.get(key), 'Actual terminal whole-unit ' + key)
    counter = fields['CPUUsageNSec']
    expected = 'oom-kill' if identity == 'original' else 'success'
    require(fields['Result'] == expected and isinstance(counter, str) and
            re.fullmatch(r'[0-9]+', counter) and str(int(counter)) == counter and int(counter) > 0 and
            re.fullmatch(r'[0-9a-fA-F]{32}', fields['InvocationID']) and fields.get('MemoryPeak'),
            'Missing/noncanonical actual terminal result/counter/invocation')
    require(set(receipt['files']) == UNIT_FILES, 'Incomplete actual terminal unit receipt')
    files = {name: bound_file(root, item, 'Actual terminal launch ' + name)
             for name, item in receipt['files'].items()}
    for name, path in files.items():
        require(path == root / identity / name, 'Terminal launch file moved outside original unit')
    launch = read_json(files['launch-receipt.json'])
    observed = manager(files['05-manager-bounds.stdout'])
    command = read_json(files['05-manager-bounds.argv.json'])
    start = command.index('systemctl')
    require(command[start:start+4] == ['systemctl', '--user', 'show', launch['unit']] and
            'InvocationID' in command[start+4:] and launch['unit'] == receipt['unit'] ==
            'coherence-' + launch['name'] + '.service' and launch['phase'] == 'bo' and
            launch['case'] == execution['case'] and observed['InvocationID'] == fields['InvocationID'] and
            observed['CPUAffinity'] == str(launch['cpu']) and observed['Restart'] == 'no',
            'Terminal accounting belongs to another launch/invocation')
    compare_json(read_json(files['input.json']), execution['specification'], 'Actual terminal launch input')
    for key in ('bytes', 'sha256'):
        require(launch['input'][key] == binding(files['input.json'])[key] == descriptor['input'][key] and
                launch['frozen_source'][key] == binding(files['source-manifest.json'])[key] ==
                descriptor['source_epoch'][key], 'Terminal launch source/input binding changed')
    text = files['native.log'].read_text()
    marker = ' case=' + execution['case'] + ' phase=bo'
    require(text.count('MICRO_PD_START' + marker) == 1 and
            text.count('MICRO_PD_COMPLETE') == (identity == 'recovery') and
            (identity != 'recovery' or text.count('MICRO_PD_COMPLETE' + marker) == 1),
            'False original success or incomplete actual recovery Native execution')
    events = re.findall(r'^MICRO_SCF_(START|COMPLETE) state=(\S+)$', text, re.MULTILINE)
    stems = ([*HOSTS, *('occupied-probes-' + format(i, '02d') for i in range(6)),
              *('vacant-probes-' + format(i, '02d') for i in range(6)), 'occupied-probe-He']
             if identity == 'original' else ['vacant-probe-H-R0'])
    expected_events = [(event, stem) for stem in stems for event in ('START', 'COMPLETE')]
    if identity == 'original':
        expected_events.append(('START', 'vacant-probe-H-R0'))
    compare_json(events, expected_events, 'Actual complete SCF sequence and terminal interruption')
    return dict(unit=receipt['unit'], Result=expected, scope=UNIT_SCOPE,
                CPUUsageNSec=counter, InvocationID=fields['InvocationID'], MemoryPeak=fields['MemoryPeak'])


def replay(directory, *, log_inventory=None, archive_root=None):
    import reader_pd as pd
    from reader_pd_chain import epoch as validate_epoch

    root = Path(directory).resolve()
    original_binding = binding(root / 'case.json')
    case = read_json(root / 'case.json')
    require(case['schema'] == SCHEMA, 'Wrong recovered terminal schema')
    reader._no_extra_electrons(case)
    view = case['derived_view']
    require(view['schema'] == 'micro-site-case-v1' and view['case'] == case['case'] and
            isinstance(case['scope'], str) and case['scope'], 'Missing recovered derivation identity/scope')
    descriptors = dict(original=case['original_terminal'], recovery=case['recovery_component'])
    roots, executions, epochs, sources, units = {}, {}, {}, {}, {}
    for identity, descriptor in descriptors.items():
        base, execution = _component(root, descriptor, identity)
        roots[identity], executions[identity] = base, execution
        spec = execution['specification']
        if identity == 'recovery':
            path = bound_file(root, descriptor['branch'], 'Actual closed H-R0 branch')
            require(path == base / 'branch.json', 'Recovery branch/root mismatch')
            branch = read_json(path)
            require(branch['schema'] == 'micro-pd-bo-branch-v1' and
                    branch['bo_branch_role'] == spec.get('bo_branch_role') == 'contact-h0' and
                    spec.get('delivery_authority_entry') == 'entry_2994' and
                    branch['execution_radii_A'] == [] and set(branch['states']) == {'vacant-probe-H-R0'} and
                    'contact' not in branch and not branch.get('host_components'),
                    'Recovery must execute exactly fresh H-R0 without adopted hosts/contact')
            for key in ('case', 'specification', 'geometry', 'source_epoch', 'provenance'):
                compare_json(branch[key], execution[key], 'Actual recovery branch ' + key)
            execution = executions[identity] = branch
            require(0 < number(spec['cpu_limit_s'], 'Authorized fresh H-R0 process cap') <= 1200,
                    'Fresh H-R0 cap exceeds authorized twenty core-minutes')
        else:
            require(spec.get('bo_branch_role') is None and spec['radii_A'] == reader.RADII and
                    not (base / 'case.json').exists() and not (base / 'branch.json').exists() and
                    not (base / 'vacant-probe-H-R0.json').exists(),
                    'Original killed monolithic run was rewritten as a successful case/branch/state')
        pd._spec(spec, execution)
        epochs[identity], sources[identity] = reader._epoch(base, execution,
            required_sources=(*pd.EPOCH_SOURCES, *(('pd_bo_worker.py',) if identity == 'recovery' else ())))
        compare_json(validate_epoch(base / 'provenance/source-manifest.json', spec), epochs[identity],
                     'Actual input-bound terminal source epoch')
        pd._geometry(base, execution, log_inventory, archive_root)
        units[identity] = _unit(root, descriptor, execution, identity)
    original, recovery = executions['original'], executions['recovery']
    require(units['original']['InvocationID'] != units['recovery']['InvocationID'] and
            units['original']['unit'] != units['recovery']['unit'], 'Recovery reused original unit/invocation')
    for key in (key for key in METHOD if key != 'delivery_authority_entry'):
        compare_json(original['specification'].get(key), recovery['specification'].get(key),
                     'Unchanged original/recovery physical method ' + key)
    compare_json(pd._physical_sources(epochs['original']), pd._physical_sources(epochs['recovery']),
                 'Unchanged original/recovery physical source bytes')
    compare_json(original['geometry'], recovery['geometry'], 'Identical actual original frame/site/ray')
    completed = case['completed_radii_A']
    require(completed == reader.RADII and case['missing_radii_A'] == [], 'Recovered base must retain all six original radii')
    compare_json(view['specification'], dict(original['specification'], case=case['case'], radii_A=completed),
                 'Derived specification is not an executed recovery specification')
    for key in ('geometry', 'source_epoch', 'provenance'):
        compare_json(view[key], original[key], 'Derived original ' + key + ' anchor')
    config = read_json(bound_file(root, case['collector_provenance']['input'], 'Actual recovery collection input'))
    require(config['schema'] == 'micro-pd-terminal-recovery-input-v1' and config['case'] == case['case'] and
            set(config['unit_accounting']) == {'original', 'recovery'}, 'Wrong actual recovery collection input')
    for identity in descriptors:
        for key in ('bytes', 'sha256'):
            require(config['unit_accounting'][identity][key] == descriptors[identity]['unit_receipt']['accounting'][key],
                    'Collector consumed different terminal accounting')
    for field, item in (('original_input', descriptors['original']['input']),
                        ('recovery_branch', descriptors['recovery']['branch'])):
        require(all(config[field][key] == item[key] for key in ('bytes', 'sha256')),
                'Collector consumed different original input/recovery branch')
    collector_epoch = validate_epoch(bound_file(root, case['collector_source_epoch'], 'Actual recovery collector source'),
        dict(original['specification'], frozen_source_manifest=config['frozen_source_manifest']))
    require(sum(Path(item['path']).name == 'pd_terminal_recover.py' for item in collector_epoch['files']) == 1,
            'Missing actual terminal collector source')
    compare_json(pd._physical_sources(epochs['original']), pd._physical_sources(collector_epoch),
                 'Unchanged collector physical source bytes')

    g, states = view['geometry'], view['states']
    positions, cell, masses = pd._frame(g, original['specification']['variant'])
    direction, target = np.asarray(g['direction'], float), g['target_index']
    require(set(states) == set(case['state_origins']) == STATE_ROLES, 'Missing exact recovered state origins')
    energies, seen, setups, grid = {}, set(), {}, None
    for role, value in states.items():
        records = value if isinstance(value, list) else [value]
        origins = case['state_origins'][role]
        require(len(records) == len(origins) == (6 if role.endswith('-probes') else 1), 'Incomplete recovered state inventory')
        values = []
        for index, (state, origin) in enumerate(zip(records, origins)):
            identity = 'recovery' if role == 'vacant-probe-H-R0' else 'original'
            require(origin['component'] == identity and origin['state_file'] == state['state_file'],
                    'Changed original state identity/root')
            stem = role + '-' + format(index, '02d') if role.endswith('-probes') else role
            require(state['state_file']['path'] == stem + '.json' and (identity, stem) not in seen,
                    'Reordered/duplicate original state file')
            seen.add((identity, stem))
            if identity == 'recovery':
                compare_json(state, recovery['states'][role], 'Unchanged genuinely executed H-R0 state')
            energy = reader._state(roots[identity], state, role, executions[identity]['specification'], g,
                                   False, log_inventory, archive_root, valence=pd.VALENCE)
            ss, pp, mm = list(g['symbols']), positions.copy(), masses.copy()
            if role in ('vacant-host', 'vacant-probes'):
                del ss[target]
                pp, mm = np.delete(pp, target, axis=0), np.delete(mm, target)
            if role.endswith('-probes'):
                radius = completed[index]
                projectile = positions[target] - radius * direction
                close(state['R_A'], radius, 'Original descending radius')
                close(state['projectile_A'], projectile, 'Original unmodified projectile ray')
                ss.append('H')
                pp, mm = np.vstack((pp, projectile)), np.append(mm, masses[target])
            elif role == 'occupied-probe-He':
                ss[target], mm[target] = 'He', 2 * mm[target]
            require(state['sequence_index'] == (index if role.endswith('-probes') else 0) and
                    state['warm_sequence'] == role and
                    state['wavefunction_policy'] == 'one live calculator across same-inventory descending-R states',
                    'Changed actual sequence order/calculator policy')
            if role in ('occupied-probe-He', 'vacant-probe-H-R0'):
                close(state['R_A'], 0., 'Exact contact radius')
                close(state['projectile_A'], positions[target], 'Exact contact location')
            require(state['symbols'] == ss and state['charge'] == (0 if role in HOSTS else 1),
                    'Wrong recovered nuclear/charge role')
            close(state['positions_A'], pp, 'Original frozen nuclear coordinates')
            close(state['masses_u'], mm, 'Original frozen isotope inventory')
            if grid is None:
                grid = state['realspace_shape']
            require(state['realspace_shape'] == grid, 'Unequal recovered quartet grids')
            for symbol, setup in state['actual_setups'].items():
                if symbol in setups:
                    compare_json(setup, setups[symbol], 'Same actual loaded setup ' + symbol)
                else:
                    setups[symbol] = setup
            require(number(state['process_cpu_s'], 'Actual cumulative sequence CPU') > 0,
                    'Missing actual state process clock')
            values.append(energy)
        if role.endswith('-probes'):
            require(all(records[i]['process_cpu_s'] < records[i+1]['process_cpu_s'] for i in range(5)),
                    'Original cumulative sequence clocks rewritten')
        energies[role] = values if role.endswith('-probes') else values[0]
    for item in recovery['logs']:
        reader._log(roots['recovery'], item, log_inventory, archive_root)
    green, rows = reader.PeriodicNuclearGreen(cell), []
    for index, radius in enumerate(completed):
        q = math.fsum([energies['occupied-probes'][index], -energies['vacant-probes'][index],
                       -energies['occupied-host'], energies['vacant-host']])
        correction = green.point_minus_gaussian_eV(radius * direction, math.sqrt(5.) / reader.BOHR_A)
        total = q + correction
        rows.append(dict(R_A=radius, BO_quartet_soft_eV=q, point_minus_gaussian_eV=correction,
                         V_pair_total_eV=total, lowering_eV=reader.COULOMB_EV_A/radius-total))
    compare_json(view['rows'], rows, 'Independent recovered quartet/Green replay')
    cpu = number(recovery['process_cpu_s'], 'Actual successful H-R0 process CPU')
    require(cpu > 0 and int(units['recovery']['CPUUsageNSec']) / 1e9 + 1 >= cpu,
            'Recovery process exceeds actual whole-unit CPU')
    close(view['process_cpu_s'], cpu, 'Only genuinely observed successful recovery process CPU')
    require(case['process_cpu_s_scope'] == PROCESS_SCOPE, 'False original/recovery process clock scope')
    require(binding(root / 'case.json') == original_binding, 'Recovered case changed during replay')
    return dict(schema='micro-site-replay-v1', case=case['case'], specification=view['specification'],
                geometry=g, rows=rows, states=states, source_epoch=original['source_epoch'],
                source_manifest=epochs['original'], case_binding=original_binding,
                raw_quartet_energies_eV=energies, completed_radii_A=completed, missing_radii_A=[],
                container_schema=SCHEMA, original_terminal=descriptors['original'],
                recovery_component=descriptors['recovery'], state_origins=case['state_origins'],
                recovery_scope=case['scope'], cpu_accounting=dict(scope=PROCESS_SCOPE,
                    process_cpu_s=cpu, original_process_cpu_s=None, units=units,
                    components=dict(original=None, recovery=cpu)),
                **_contact(roots['original'], view, sources['original'], states, energies, cell,
                           original['specification']))
