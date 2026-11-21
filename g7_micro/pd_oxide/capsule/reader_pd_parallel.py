"""Replay parallel originals at their actual roots/specs; derived views are not executions."""
from pathlib import Path
import math
import re
import numpy as np
import reader
from reader import binding, bound_file, close, compare_json, number, read_json, require

HOSTS = ('occupied-host', 'vacant-host')
METHOD = ('material', 'variant', 'size', 'h_A', 'kmesh', 'eigensolver', 'force_tolerance_eV_A',
          'ionic_steps', 'selected_geometry', 'probe_ray', 'geometry_policy_entry', 'ray_policy_entry',
          'delivery_authority_entry', 'dataset_bindings', 'library_bindings')


def _unit(root, descriptor, branch):
    from reader_pd_label import manager
    receipt = descriptor['unit_receipt']
    require(receipt['counter_scope'] == 'Exact original whole unit; not the collector or an adoption clock',
            'Wrong terminal unit counter scope')
    fields = manager(bound_file(root, receipt['accounting'], 'Actual successful unit accounting'))
    for key in ('Result', 'CPUUsageNSec', 'InvocationID', 'MemoryPeak'):
        compare_json(receipt[key], fields.get(key), 'Original whole-unit ' + key)
    counter = fields['CPUUsageNSec']
    require(fields['Result'] == 'success' and isinstance(counter, str) and
            re.fullmatch(r'[0-9]+', counter) and str(int(counter)) == counter and
            re.fullmatch(r'[0-9a-fA-F]{32}', fields['InvocationID']), 'Noncanonical/missing successful unit receipt')
    files = {key: bound_file(root, item, 'Closed actual unit ' + key) for key, item in receipt['files'].items()}
    require(set(files) == {'input.json', 'source-manifest.json', 'launch-config.json', 'launch-receipt.json',
                          'geometry-admission.json', '05-manager-bounds.stdout',
                          '05-manager-bounds.argv.json', 'native.log'}, 'Incomplete actual unit closure')
    launch = read_json(files['launch-receipt.json'])
    observed = manager(files['05-manager-bounds.stdout'])
    command = read_json(files['05-manager-bounds.argv.json'])
    start = command.index('systemctl')
    require(command[start:start+4] == ['systemctl', '--user', 'show', launch['unit']] and
            'InvocationID' in command[start+4:] and launch['unit'] == receipt['unit'] ==
            'coherence-' + launch['name'] + '.service' and launch['phase'] == 'bo' and
            launch['case'] == branch['case'] and observed['InvocationID'] == receipt['InvocationID'] and
            observed['CPUAffinity'] == str(launch['cpu']) and observed['Restart'] == 'no',
            'Unit receipt belongs to another invocation')
    compare_json(read_json(files['input.json']), branch['specification'], 'Actual launched branch input')
    for key in ('bytes', 'sha256'):
        require(launch['input'][key] == binding(files['input.json'])[key] and
                launch['frozen_source'][key] == binding(files['source-manifest.json'])[key] ==
                branch['source_epoch'][key], 'Actual launch input/source bindings changed')
    marker = ' case=' + branch['case'] + ' phase=bo'
    text = files['native.log'].read_text()
    require(text.count('MICRO_PD_START' + marker) == text.count('MICRO_PD_COMPLETE' + marker) == 1,
            'Actual Native branch did not close exactly once')
    require(int(counter) / 1e9 + 1 >= branch['process_cpu_s'], 'Actual process exceeds whole-unit CPU')
    return dict(unit=receipt['unit'], scope=receipt['counter_scope'], CPUUsageNSec=counter,
                InvocationID=receipt['InvocationID'])


def _contact(root, branch, sources, states, energies, cell, spec):
    contact = branch['contact']
    convention = contact['absolute_convention']
    for key, expected in dict(alpha1=10., alpha2=300., H_E_reference_Hartree=0., He_E_reference_Hartree=0.,
                              representation='same-alpha all-electron Z2', arbitrary_atomic_offset_eV=0.,
                              energy='internal=free-minus-negativeTS').items():
        require(convention[key] == expected, 'Wrong contact absolute convention ' + key)
    operators = contact['operator_sources']
    require(set(operators) == {'ae', 'setup', 'basis_data', 'hamiltonian', 'density'}, 'Missing contact operator closure')
    opaths = {key: bound_file(root, item, 'Actual contact operator') for key, item in operators.items()}
    mapping = next(path for name, path in sources.items() if Path(name).name == 'collision_bare_reference.py')
    proof = reader.verify_mapping(cell, convention['H_M_Hartree'], convention['He_M_Hartree'],
                                  convention['ASE_Hartree_eV'], opaths['ae'], mapping)
    units = reader.verify_ase_units(bound_file(root, convention['ASE_units_source'], 'Actual ASE units'),
                                    convention['ASE_Hartree_eV'], convention['ASE_CODATA_version'])
    for actual in (operators['ae'], operators['setup'], convention['ASE_units_source']):
        require(any(all(actual[k] == row[k] for k in ('bytes', 'sha256')) for row in spec['library_bindings']),
                'Actual contact operator/unit source not consumed by contact input')
    close(convention['B_native_eV'], proof['B_native_eV'], 'Contact source/native B')
    identity_record = contact['operator_identity']
    require(identity_record['schema'] == 'hydride-contact-operator-identity-v1' and identity_record['source'] == operators,
            'Changed constructed operator sources')
    coincident = identity_record['coincident_native']
    require(coincident['performed'] is False and coincident['reason'], 'Unqualified coincident Native claim')
    guards = {key: bound_file(root, item, 'Unmodified public guard') for key, item in coincident['guard_sources'].items()}
    radial = read_json(bound_file(root, identity_record['radial_file'], 'Constructed radial operator'))
    for role, records in states.items():
        for state in records if isinstance(records, list) else [records]:
            for symbol in ('H', 'He'):
                if symbol in state['actual_setups']:
                    for key in ('Z', 'Nv', 'Nc', 'ni', 'E', 'alpha1', 'alpha2', 'M', 'Delta0'):
                        close(state['actual_setups'][symbol][key], radial[symbol][key], 'Actual common contact operator ' + role + '/' + symbol + '/' + key)
    identity = reader.verify_operator_identity(radial, opaths, guards, convention, 1, states['occupied-probe-He']['electrons'])
    close(identity['constructed_operator_difference_eV'], proof['B_native_eV'], 'Common-density Bnative identity')
    q = math.fsum([energies['occupied-probe-He'], -energies['vacant-probe-H-R0'],
                   -energies['occupied-host'], energies['vacant-host']])
    signed = math.fsum([q, proof['unit_constant_residual_eV'], proof['mapping']['regular_images_eV']])
    close(contact['Q_He_eV'], q, 'Independent contact quartet')
    close(contact['signed_W0_eV'], signed, 'Independent regular signed contact')
    return dict(signedW0_eV=signed, Q_He_eV=q, absolute_convention=convention,
                operator_sources=operators, mapping_proof=proof, unit_source_proof=units,
                functional_identity_replay=identity)


def replay(directory, *, log_inventory=None, archive_root=None):
    import reader_pd as pd
    root = Path(directory).resolve()
    original = binding(root / 'case.json')
    case = read_json(root / 'case.json')
    require(case['schema'] == 'micro-pd-parallel-case-v1', 'Wrong parallel container schema')
    reader._no_extra_electrons(case)
    components, epochs, source_maps, roots = {}, {}, {}, {}
    unit_counters = {}
    energies_by_origin = {}
    all_setups = {}
    grid = None
    for descriptor in case['components']:
        role = descriptor['id']
        require(role == descriptor['bo_branch_role'] and role not in components, 'Duplicate/wrong component identity')
        base = (root / descriptor['root']).resolve()
        require(base.is_relative_to(root) and descriptor['root'] == 'components/' + role + '/calculation', 'Unsafe component root')
        files = descriptor['files']
        require(len({item['path'] for item in files}) == len(files), 'Duplicate archived component member')
        for item in files:
            bound_file(root, item, 'Immutable component original')
        require({item['path'] for item in files} ==
                {path.relative_to(root).as_posix() for path in base.rglob('*') if path.is_file()},
                'Missing/unbound component closure member')
        branch_path = bound_file(root, descriptor['branch'], 'Original closed branch')
        require(branch_path == base / 'branch.json', 'Branch/root mismatch')
        branch = read_json(branch_path)
        spec = branch['specification']
        require(branch['schema'] == 'micro-pd-bo-branch-v1' and branch['bo_branch_role'] == role == spec['bo_branch_role'] and
                branch['case'] == spec['case'] and spec['delivery_authority_entry'] == 'entry_2994', 'Wrong actual branch identity/authority')
        pd._spec(spec, branch)
        compare_json(read_json(bound_file(root, descriptor['input'], 'Actual branch input')), spec, 'Actual input/specification')
        require(bound_file(root, descriptor['input'], 'Branch input location') == base / 'provenance/input.json' and
                bound_file(root, descriptor['source_epoch'], 'Branch epoch location') == base / 'provenance/source-manifest.json', 'Wrong component provenance root')
        epoch, sources = reader._epoch(base, branch, required_sources=(*pd.EPOCH_SOURCES, 'pd_bo_worker.py'))
        compare_json(descriptor['source_epoch'], dict(branch['source_epoch'], path=descriptor['root'] + '/' + branch['source_epoch']['path']), 'Component epoch binding')
        p, cell, masses, direction = pd._geometry(base, branch, log_inventory, archive_root)
        if components:
            first = next(iter(components.values()))
            for key in METHOD:
                compare_json(spec.get(key), first['specification'].get(key), 'Common physical branch method ' + key)
            compare_json(pd._physical_sources(epoch), pd._physical_sources(next(iter(epochs.values()))), 'Unchanged Native physical sources')
            compare_json(branch['geometry'], first['geometry'], 'Identical selected frame/site/ray')
        expected = ({role} if role in HOSTS else {'occupied-probe-He', 'vacant-probe-H-R0'} if role == 'contact' else {role.rsplit('-', 1)[0]})
        require(role in (*HOSTS, 'occupied-probes-near', 'vacant-probes-near', 'occupied-probes-far', 'vacant-probes-far', 'contact') and
                set(branch['states']) == expected, 'Wrong complete branch state roles')
        execution = [] if role in (*HOSTS, 'contact') else spec['radii_A']
        require(branch['execution_radii_A'] == execution and (not execution or execution ==
                ([.1] if spec['radii_A'] == [.1] else [.1, .2, .4] if role.endswith('-near') else [.8, 1.2, 1.6])), 'Changed branch execution order')
        target, symbols = branch['geometry']['target_index'], branch['geometry']['symbols']
        for state_role, value in branch['states'].items():
            records = value if isinstance(value, list) else [value]
            require(len(records) == (len(execution) if execution else 1), 'Incomplete actual branch')
            for ordinal, state in enumerate(records):
                energy = reader._state(base, state, state_role, spec, branch['geometry'], False, log_inventory, archive_root, valence=pd.VALENCE)
                ss, pp, mm = list(symbols), p.copy(), masses.copy()
                if state_role in ('vacant-host', 'vacant-probes'):
                    del ss[target]
                    pp, mm = np.delete(pp, target, axis=0), np.delete(mm, target)
                charge = 0 if state_role in HOSTS else 1
                if execution:
                    radius = execution[ordinal]
                    projectile = p[target] - radius * direction
                    require(state['R_A'] == radius and state['sequence_index'] == ordinal and state['warm_sequence'] == state_role,
                            'Changed actual warm-sequence order')
                    close(state['projectile_A'], projectile, 'Actual projectile ray')
                    ss.append('H')
                    pp, mm = np.vstack((pp, projectile)), np.append(mm, masses[target])
                elif state_role == 'occupied-probe-He':
                    ss[target], mm[target] = 'He', 2 * mm[target]
                if role == 'contact':
                    close(state['R_A'], 0., 'Exact contact R')
                    close(state['projectile_A'], p[target], 'Exact contact location')
                require(state['symbols'] == ss and state['charge'] == charge, 'Wrong actual nuclear/charge role')
                close(state['positions_A'], pp, 'Frozen branch nuclear coordinates')
                close(state['masses_u'], mm, 'Frozen branch isotope inventory')
                if grid is None:
                    grid = state['realspace_shape']
                require(state['realspace_shape'] == grid, 'Unequal physical quartet grids')
                for symbol, setup in state['actual_setups'].items():
                    if symbol in all_setups:
                        compare_json(setup, all_setups[symbol], 'Common actual loaded setup ' + symbol)
                    else:
                        all_setups[symbol] = setup
                key = (role, state['state_file']['path'])
                require(key not in energies_by_origin, 'Duplicate original state')
                energies_by_origin[key] = energy
        for item in branch['logs']:
            reader._log(base, item, log_inventory, archive_root)
        require(number(branch['process_cpu_s'], 'Actual component process CPU') > 0, 'Missing component CPU')
        components[role], epochs[role], source_maps[role], roots[role] = branch, epoch, sources, base
        unit_counters[role] = _unit(root, descriptor, branch)
    require(set(HOSTS) <= components.keys() and 'contact' in components and
            {'occupied-probes-near', 'vacant-probes-near'} <= components.keys(), 'Missing near/host/contact branch')
    completed = case['completed_radii_A']
    require(completed in ([.1], [.4, .2, .1], reader.RADII) and
            case['missing_radii_A'] == [r for r in reader.RADII if r not in completed], 'False completed/missing R inventory')
    require(('occupied-probes-far' in components) == ('vacant-probes-far' in components) == (completed == reader.RADII), 'Incomplete/unreported far branches')
    view = case['derived_view']
    host = components['occupied-host']
    require(view['schema'] == 'micro-site-case-v1' and view['case'] == case['case'], 'Wrong derived container identity')
    compare_json(view['specification'], dict(host['specification'], case=case['case'], radii_A=completed), 'Derived specification is not an executed monolithic spec')
    compare_json(view['geometry'], host['geometry'], 'Derived actual geometry')
    compare_json(view['source_epoch'], host['source_epoch'], 'Derived epoch anchor')
    compare_json(view['provenance'], host['provenance'], 'Derived provenance anchor')
    from reader_pd_chain import epoch as validate_epoch
    collector_path = bound_file(root, case['collector_source_epoch'], 'Actual collector source epoch')
    config_path = bound_file(root, case['collector_provenance']['input'], 'Actual collector input')
    config = read_json(config_path)
    require(config['schema'] == 'micro-pd-bo-collection-input-v1' and config['case'] == case['case'] and
            set(config['branches']) == set(components) == set(config['unit_accounting']),
            'Changed actual collector input/branch inventory')
    collector_epoch = validate_epoch(collector_path, dict(host['specification'],
                                     frozen_source_manifest=config['frozen_source_manifest']))
    require(sum(Path(item['path']).name == 'pd_bo_collect.py' for item in collector_epoch['files']) == 1,
            'Missing actual collector program source')
    for descriptor in case['components']:
        role = descriptor['id']
        for key in ('bytes', 'sha256'):
            require(config['branches'][role][key] == descriptor['branch'][key] and
                    config['unit_accounting'][role][key] == descriptor['unit_receipt']['accounting'][key],
                    'Collector consumed different original branch/accounting')
    states, energies = view['states'], {}
    expected_roles = {*HOSTS, 'occupied-probes', 'vacant-probes', 'occupied-probe-He', 'vacant-probe-H-R0'}
    require(set(states) == set(case['state_origins']) == expected_roles, 'Missing derived exact state roles')
    seen = set()
    for role, value in states.items():
        records = value if isinstance(value, list) else [value]
        origins = case['state_origins'][role]
        require(len(records) == len(origins) == (len(completed) if role.endswith('-probes') else 1), 'Incomplete canonical R/state origin map')
        values = []
        for ordinal, (state, origin) in enumerate(zip(records, origins)):
            component = origin['component']
            require(component in components and role in components[component]['states'], 'Wrong state origin role')
            actual = components[component]['states'][role]
            actual = actual if isinstance(actual, list) else [actual]
            matches = [item for item in actual if item['state_file'] == origin['state_file']]
            require(len(matches) == 1, 'Ambiguous original state origin')
            compare_json(state, matches[0], 'Derived view preserves actual original state verbatim')
            key = (component, origin['state_file']['path'])
            require(key not in seen, 'Repeated original state')
            seen.add(key)
            if role.endswith('-probes'):
                close(state['R_A'], completed[ordinal], 'Canonical mapping by actual R value')
            values.append(energies_by_origin[key])
        energies[role] = values if role.endswith('-probes') else values[0]
    require(seen == set(energies_by_origin), 'Unconsumed/unreported original states')
    contact_branch = components['contact']
    require(case['contact_component'] == 'contact', 'Wrong contact origin')
    compare_json(view['contact'], contact_branch['contact'], 'Unchanged original contact record')
    require({item['bo_branch_role'] for item in contact_branch['host_components']} == set(HOSTS), 'Missing contact actual host bindings')
    for item in contact_branch['host_components']:
        role = item['bo_branch_role']
        compare_json(read_json(bound_file(roots['contact'], item['branch'], 'Consumed original contact host')), components[role], 'Contact reused exact fixed host')
        for key in ('bytes', 'sha256'):
            require(item['original_branch'][key] == binding(roots[role] / 'branch.json')[key], 'Contact host original identity changed')
    g = view['geometry']
    green, direction = reader.PeriodicNuclearGreen(g['cell_A']), np.asarray(g['direction'])
    rows = []
    for i, radius in enumerate(completed):
        q = math.fsum([energies['occupied-probes'][i], -energies['vacant-probes'][i], -energies['occupied-host'], energies['vacant-host']])
        correction = green.point_minus_gaussian_eV(radius * direction, math.sqrt(5.) / reader.BOHR_A)
        total = q + correction
        rows.append(dict(R_A=radius, BO_quartet_soft_eV=q, point_minus_gaussian_eV=correction,
                         V_pair_total_eV=total, lowering_eV=reader.COULOMB_EV_A/radius-total))
    compare_json(view['rows'], rows, 'Independent original-branch quartet/Green replay')
    cpu = math.fsum(branch['process_cpu_s'] for branch in components.values())
    close(view['process_cpu_s'], cpu, 'Separate component process CPU sum')
    require(case['process_cpu_s_scope'] == 'sum of actual component processes; not a unit clock', 'False composite clock scope')
    require(binding(root / 'case.json') == original, 'Composite changed during replay')
    return dict(schema='micro-site-replay-v1', case=case['case'], specification=view['specification'], geometry=g,
                rows=rows, states=states, source_epoch=host['source_epoch'], source_manifest=epochs['occupied-host'],
                case_binding=original, raw_quartet_energies_eV=energies, completed_radii_A=completed,
                missing_radii_A=case['missing_radii_A'], components=case['components'], state_origins=case['state_origins'],
                container_schema=case['schema'], cpu_accounting=dict(scope=case['process_cpu_s_scope'], process_cpu_s=cpu,
                    units=unit_counters,
                    components={role: branch['process_cpu_s'] for role, branch in components.items()}),
                **_contact(roots['contact'], contact_branch, source_maps['contact'], states, energies, np.asarray(g['cell_A']), contact_branch['specification']))
