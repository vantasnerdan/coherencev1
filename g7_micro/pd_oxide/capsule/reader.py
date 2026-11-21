"""Replay actual Zr micro-calibration saved states using NumPy and stdlib only.

Frozen producer/source files are data: they are hashed and parsed, never imported.
No synthetic/analytic state, inferred missing control, or isolated-limit admission.
"""
import argparse
from collections import Counter
import gzip
import xml.etree.ElementTree as ET
import itertools
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True

import numpy as np

HERE = Path(__file__).resolve().parent
SUPPORT = HERE / 'support'
if not SUPPORT.is_dir():
    SUPPORT = HERE.parents[1] / 'hydride_screening'
sys.path.insert(0, str(SUPPORT))
from hydride_replay import (binding, bound_file, check_binding, close, compare_json,
                            number, read_json, require)
from green import BOHR_A, COULOMB_EV_A, PeriodicNuclearGreen
from verify_contact_mapping import verify_ase_units, verify_mapping, verify_operator_identity
from log_privacy import allowed_log_transform, sanitize_log

RADII = [1.6, 1.2, .8, .4, .2, .1]
VALENCE = {'Zr': 12, 'O': 6, 'H': 1, 'He': 2}
DATASET_VALENCE = VALENCE | {'Pd': 10}
EPOCH_SOURCES = ('native.py', 'run_zr.py', 'geometry_zr.py', 'contact_setup.py',
                 'collision_bare_reference.py')
ROW_KEYS = {'R_A', 'BO_quartet_soft_eV', 'point_minus_gaussian_eV',
            'V_pair_total_eV', 'lowering_eV'}
SCOPE = ('Actual frozen finite-periodic-cell target-conditioned internal BO energies; '
         'q+1 probe adds no electron. No exact isolated, h/k/cell limit, Ueq, '
         'action/Eq7, site stability, or experimental phase assignment follows.')


def _no_extra_electrons(value):
    if isinstance(value, dict):
        require(not ({'extraNe', 'extra_electrons', 'added_electrons'} & set(value)),
                'Uncommissioned extra-electron input')
        for child in value.values():
            _no_extra_electrons(child)
    elif isinstance(value, list):
        for child in value:
            _no_extra_electrons(child)


def _log(root, item, inventory, archive_root):
    if inventory is None:
        return bound_file(root, item, 'Final Native log')
    check_binding(item, 'Original Native log')
    name = (Path(root).resolve().relative_to(Path(archive_root).resolve()) / item['path']).as_posix()
    require(name in inventory, 'Missing original/exported log accounting')
    entry = inventory[name]
    require(set(entry) == {'ORIGINAL', 'EXPORTED', 'transformation'} and
            entry['ORIGINAL'] == {k: item[k] for k in ('bytes', 'sha256')},
            'Changed original Native log binding')
    path = bound_file(root, dict(path=item['path'], **entry['EXPORTED']), 'Exported log')
    changes = entry['transformation']
    if changes:
        require(allowed_log_transform(changes) and sanitize_log(path.read_bytes())[1] == [] and
                entry['ORIGINAL']['bytes'] > entry['EXPORTED']['bytes'],
                'Unapproved Native log transformation')
    else:
        require(changes == [] and entry['ORIGINAL'] == entry['EXPORTED'], 'Undeclared log change')
    return path


def _coordination(geometry, positions, cell, target):
    delta = positions - positions[target]
    sigma = float(np.linalg.svd(cell, compute_uv=False)[-1])
    shell = math.ceil((3.5 + float(np.linalg.norm(delta, axis=1).max())) / sigma)
    actual = {}
    oxygen_distances = []
    for shift in itertools.product(range(-shell, shell + 1), repeat=3):
        vectors = delta + np.asarray(shift) @ cell
        distances = np.linalg.norm(vectors, axis=1)
        for index in np.flatnonzero(distances < 3.5):
            if index == target and shift == (0, 0, 0):
                continue
            actual[index, shift] = float(distances[index])
        oxygen = geometry['initial_design']['nominal_site']['oxygen_index']
        if oxygen is not None:
            oxygen_distances.append(float(distances[oxygen]))
    coordination = geometry['coordination']
    saved = {}
    for row in coordination['neighbors_within_3_5_A']:
        key = row['index'], tuple(row['image'])
        require(key not in saved and row['symbol'] == geometry['symbols'][row['index']],
                'Invalid/duplicate reported coordination neighbor')
        saved[key] = row['distance_A']
    require(set(saved) == set(actual), 'Incomplete actual target coordination')
    for key, distance in actual.items():
        close(saved[key], distance, 'Actual periodic neighbor distance')
    final_zr = sorted(index for (index, _), distance in actual.items()
                      if geometry['symbols'][index] == 'Zr' and distance < 3.)
    require(coordination['final_Zr_indices_within_3_A'] == final_zr, 'Wrong final Zr coordination')
    displacements = geometry['displacements']
    lengths = np.linalg.norm(np.asarray(displacements['vectors_A'], dtype=float), axis=1)
    close(displacements['distances_A'], lengths, 'Reported displacement vector lengths')
    close(displacements['target_A'], float(lengths[target]), 'Reported target displacement')
    close(displacements['maximum_A'], float(lengths.max()), 'Reported maximum displacement')
    retained = coordination['initial_Zr_indices_within_3_A'] == final_zr and lengths[target] < .5
    require(coordination['retained_by_reported_geometric_criterion'] is bool(retained),
            'Wrong reported retention/migration criterion')
    expected_site = ('retained under reported geometric criterion' if retained else
                     'migrated/recoordinated target; no intended-site label')
    require(coordination['physical_site'] == expected_site, 'False intended physical site label')
    if oxygen_distances:
        close(coordination['target_O_distance_A'], min(oxygen_distances), 'Actual target-O MIC distance')

def _epoch(root, run, *, required_sources=EPOCH_SOURCES):
    item = run['source_epoch']
    path = bound_file(root, item, 'Frozen source epoch')
    require(item['path'] == 'provenance/source-manifest.json', 'Nonportable source epoch')
    epoch = read_json(path)
    require(set(epoch) == {'schema', 'files'} and
            epoch['schema'] == 'micro-frozen-repository-source-v1' and
            isinstance(epoch['files'], list) and epoch['files'], 'Empty/wrong source epoch')
    files = {}
    for entry in epoch['files']:
        require(entry['path'] not in files, 'Duplicate source epoch member')
        files[entry['path']] = bound_file(path.parent / 'source', entry, 'Frozen source')
    for basename in required_sources:
        require(sum(Path(name).name == basename for name in files) == 1,
                'Missing/ambiguous source ' + basename)
    require(read_json(root / 'provenance/input.json') == run['specification'],
            'Consumed input differs from case specification')
    for field in ('dataset_bindings', 'library_bindings'):
        archived = run['provenance'][field]
        require(isinstance(archived, list) and len(archived) == len(run['specification'][field]),
                'Missing actual consumed ' + field)
        for item, original in zip(archived, run['specification'][field]):
            require(set(item) == {'original', 'local'} and item['original'] == original,
                    'Changed original consumed binding')
            copied = bound_file(root, item['local'], 'Consumed original ' + field)
            require(all(item['local'][key] == original[key] for key in ('bytes', 'sha256')),
                    'Consumed source snapshot differs from original')
            if field == 'dataset_bindings':
                atom = ET.fromstring(gzip.decompress(copied.read_bytes())).find('atom')
                require(atom is not None and atom.attrib == original['atom'], 'Changed original setup atom metadata')
                symbol = atom.attrib['symbol']
                require(symbol in DATASET_VALENCE and float(atom.attrib['valence']) == DATASET_VALENCE[symbol],
                        'Wrong source setup valence')
    for item in run['geometry']['initial_design']['structural_sources']:
        name = item.get('relative_name', item.get('path'))
        require(name in files and all(binding(files[name])[key] == item[key] for key in ('bytes', 'sha256')),
                'Actual structural primary missing/changed in frozen epoch')
    return epoch, files


def _segment_clearance(positions, cell, target, direction, *, radius_min=.1, radius_max=1.6):
    """Enumerate the finite image ball using the cell's smallest singular value."""
    spectators = np.delete(positions, target, axis=0) - positions[target]
    fractional = spectators @ np.linalg.inv(cell)
    spectators -= np.rint(fractional) @ cell
    upper = float(np.min(np.linalg.norm(spectators + radius_min * direction, axis=1)))
    bound = radius_max + upper + 1e-9
    sigma = float(np.linalg.svd(cell, compute_uv=False)[-1])
    require(sigma > 0, 'Degenerate cell')
    shell = math.ceil((bound + float(np.linalg.norm(spectators, axis=1).max())) / sigma)
    minimum = math.inf
    for shift in itertools.product(range(-shell, shell + 1), repeat=3):
        vectors = spectators + np.asarray(shift) @ cell
        mask = np.linalg.norm(vectors, axis=1) <= bound
        vectors = vectors[mask]
        if not len(vectors):
            continue
        radii = np.clip(-vectors @ direction, radius_min, radius_max)
        minimum = min(minimum, float(np.linalg.norm(vectors + radii[:, None] * direction, axis=1).min()))
    require(math.isfinite(minimum), 'No spectator image in clearance search')
    return minimum


def _geometry(root, run, inventory, archive_root):
    g, spec = run['geometry'], run['specification']
    require(g['passed'] is True and 0 <= number(g['raw_fmax_eV_A'], 'Raw fmax') <= .05,
            'Geometry not actually force-qualified')
    symbols = g['symbols']
    n = 2 if spec['size'] == 'base' else 3
    expected = Counter(Zr=n**3, H=2*n**3)
    if spec['variant'] == 'oxygen':
        expected['O'] = 1
    require(Counter(symbols) == expected, 'Wrong ZrD2/O inventory')
    positions = np.asarray(g['positions_A'], dtype=float)
    cell = np.asarray(g['cell_A'], dtype=float)
    require(positions.shape == (len(symbols), 3) and np.isfinite(positions).all(), 'Invalid actual positions')
    PeriodicNuclearGreen(cell)
    masses = np.asarray(g['masses_u'], dtype=float)
    require(masses.shape == (len(symbols),) and np.isfinite(masses).all() and bool(np.all(masses > 0)),
            'Invalid actual isotope masses')
    close(masses[np.asarray(symbols) == 'H'], [2.01410177812] * symbols.count('H'), 'Actual D masses')
    target = g['target_index']
    require(type(target) is int and 0 <= target < len(symbols) and symbols[target] == 'H', 'Wrong D target')
    direction = np.asarray(g['direction'], dtype=float)
    require(direction.shape == (3,) and np.isfinite(direction).all(), 'Invalid ray')
    close(float(np.linalg.norm(direction)), 1., 'Unit ray')
    design = g['initial_design']
    require(design['pair'] == 'zr-oxygen' and design['variant'] == spec['variant'] and
            design['size'] == spec['size'] and design['target_index'] == target, 'Wrong initial design')
    close(direction, design['direction'], 'Unchanged initial incoming ray')
    close(cell, design['cell_control']['cell_A'], 'Source-set fixed cell')
    minimum = _segment_clearance(positions, cell, target, direction)
    close(g['ray_clearance']['min_clearance_A'], minimum, 'Post-relax all-image segment clearance')
    require(minimum >= 1., 'Actual ray intersects a spectator')
    for key in ('displacements', 'coordination', 'source_geometry'):
        require(key in g, 'Missing actual geometry diagnostics ' + key)
    selected = run['provenance']['selected_geometry']
    selected_path = bound_file(root, selected, 'Selected original geometry JSON')
    original = read_json(selected_path)
    for key in ('symbols', 'positions_A', 'cell_A', 'masses_u', 'target_index', 'direction',
                'raw_fmax_eV_A', 'passed', 'initial_design', 'ray_clearance',
                'displacements', 'coordination'):
        compare_json(g[key], original[key], 'Selected geometry ' + key)
    initial_source = original['source_geometry']
    if initial_source is not None:
        initial_source = dict(initial_source)
        for key in ('local_json', 'local_gpw'):
            item = initial_source[key]
            initial_source[key] = dict(
                item, path=(selected_path.parent.relative_to(root) / item['path']).as_posix())
            bound_file(root, initial_source[key], 'Selected initial geometry ' + key)
    compare_json(g['source_geometry'], initial_source, 'Selected initial geometry provenance')
    source = spec['selected_geometry']
    require({k: selected[k] for k in ('bytes', 'sha256')} ==
            {k: source['binding'][k] for k in ('bytes', 'sha256')}, 'Selected input JSON identity changed')
    gpw = bound_file(root, g['gpw'], 'Final geometry PAW GPW')
    require({k: g['gpw'][k] for k in ('bytes', 'sha256')} ==
            {k: source['gpw_binding'][k] for k in ('bytes', 'sha256')}, 'Selected PAW GPW changed')
    raw = read_json(bound_file(root, g['state_file'], 'Final geometry PAW state'))
    _state(selected_path.parent, raw, 'geometry-PAW', raw['specification'], g, True, inventory, archive_root)
    require(raw['charge'] == 0 and raw['symbols'] == symbols, 'Wrong geometry PAW charge/inventory')
    require(original['schema'] == 'micro-site-geometry-v1' and
            original['specification']['phase'] == 'geometry' and
            original['specification']['variant'] == spec['variant'] and
            original['specification']['size'] == spec['size'], 'Wrong selected geometry phase/variant/size')
    _coordination(g, positions, cell, target)
    require(raw['gpw'] == original['gpw'] and raw['log'] == original['log'],
            'Selected PAW state checkpoint/log identity differs')
    close(raw['positions_A'], positions, 'Final PAW geometry coordinates')
    close(raw['masses_u'], masses, 'Final PAW geometry masses')
    forces = np.asarray(raw['forces_eV_A'], dtype=float)
    require(forces.shape == positions.shape and np.isfinite(forces).all(), 'Missing actual PAW forces')
    close(float(np.linalg.norm(forces, axis=1).max()), g['raw_fmax_eV_A'], 'Actual unconstrained PAW fmax')
    close(raw['max_force_eV_A'], g['raw_fmax_eV_A'], 'Final raw PAW force')
    require(gpw.is_file(), 'Missing final PAW GPW')
    return positions, cell, masses, direction


def _state(root, state, role, spec, geometry, paw, inventory, archive_root, *, valence=VALENCE):
    require(state['schema'] == 'micro-native-state-v1' and state['backend'] == 'gpaw-native' and
            state['scf_converged'] is True, role + ': actual converged Native required')
    if 'state_file' in state:
        standalone = read_json(bound_file(root, state['state_file'], role + ' standalone state'))
        require(standalone == {k: v for k, v in state.items() if k != 'state_file'}, role + ': raw state differs')
    else:
        require(paw, role + ': missing standalone state binding')
    compare_json(state['specification'], spec, role + ': consumed specification')
    _no_extra_electrons(state)
    symbols = state['symbols']
    require(isinstance(symbols, list) and set(symbols) <= set(valence), role + ': species')
    charge = state['charge']
    require(type(charge) is int and charge in (0, 1), role + ': charge')
    ne = sum(valence[s] for s in symbols) - charge
    close(number(state['electrons'], role + ' Ne'), ne, role + ': actual valence/charge Ne')
    require(state['geometry_PAW'] is paw and state['xc'] == 'PBE' and state['spinpol'] is False,
            role + ': wrong actual physics family')
    close(state['occupation_width_eV'], .1, role + ': actual smearing')
    shape = state['realspace_shape']
    require(isinstance(shape, list) and len(shape) == 3 and all(type(x) is int and x > 0 and x % 4 == 0 for x in shape),
            role + ': invalid actual grid')
    cell = np.asarray(geometry['cell_A'], dtype=float)
    lengths = 1 / np.linalg.norm(np.linalg.inv(cell), axis=0)
    expected_grid = np.maximum(4, (lengths / spec['h_A'] / 4 + .5).astype(int) * 4)
    close(shape, expected_grid, role + ': actual grid versus requested h')
    require(type(state['actual_kpoints']) is int and state['actual_kpoints'] == math.prod(spec['kmesh']),
            role + ': actual unreduced k grid')
    require(type(state['nbands']) is int and state['nbands'] == math.ceil(ne / 2) + 8,
            role + ': actual band count')
    require(type(state['scf_iterations']) is int and 0 < state['scf_iterations'] <= 150,
            role + ': actual SCF iterations')
    require(spec['eigensolver'] == {'name': 'davidson', 'niter': 2} and state['eigensolver'] == 'Davidson',
            role + ': actual explicit Davidson required')
    setups = state['actual_setups']
    require(set(setups) == set(symbols), role + ': missing actual setup metadata')
    for symbol, setup in setups.items():
        close(setup['Nv'], valence[symbol], role + ': actual setup valence')
        require(type(setup['ni']) is int and setup['ni'] >= 0 and setup['Nc'] >= 0, role + ': invalid core/projectors')
        if symbol in ('H', 'He'):
            require((setup['ni'] > 0) is (paw and symbol == 'H') and setup['Nc'] == 0,
                    role + ': wrong PAW/all-electron setup')
            if not paw:
                require(setup['type'] == 'all-electron', role + ': wrong actual Hae/He type')
    energy = state['energy']
    require(set(energy) == {'internal_energy_eV', 'free_energy_eV', 'minus_TS_eV'}, role + ': energy fields')
    for key, value in energy.items():
        number(value, role + ':' + key)
    close(energy['internal_energy_eV'], energy['free_energy_eV'] - energy['minus_TS_eV'], role + ': internal energy identity')
    require(energy['minus_TS_eV'] <= 1e-9, role + ': negativeTS sign')
    _log(root, state['log'], inventory, archive_root)
    if 'sanitized_log' in state:
        sanitized = bound_file(root, state['sanitized_log'], role + ': approved sanitized companion')
        if inventory is None:
            expected, changes = sanitize_log((root / state['log']['path']).read_bytes())
            require(sanitized.read_bytes() == expected and state['log_privacy_transforms'] == changes,
                    role + ': sanitized companion is not accepted privacy transform')
    if state.get('gpw') is not None:
        bound_file(root, state['gpw'], role + ': standard GPW')
    else:
        require('He' in symbols and not paw, role + ': missing standard checkpoint')
    close(state['cell_A'], geometry['cell_A'], role + ': frozen cell')
    require(np.asarray(state['positions_A']).shape == (len(symbols), 3), role + ': coordinate shape')
    require(np.asarray(state['masses_u']).shape == (len(symbols),), role + ': mass shape')
    return energy['internal_energy_eV']


def _replay(directory, *, log_inventory=None, archive_root=None, valence=VALENCE,
            required_sources=EPOCH_SOURCES, validate_spec=None, geometry_validator=_geometry):
    """Fail closed on raw/source corruption or changed physical input after rehashing."""
    root = Path(directory).resolve()
    path = root / 'case.json'
    original = binding(path)
    run = read_json(path)
    require(run['schema'] == 'micro-site-case-v1' and isinstance(run['case'], str), 'Wrong case schema')
    spec = run['specification']
    _no_extra_electrons(run)
    if validate_spec is None:
        require(spec['phase'] == 'bo' and spec['variant'] in ('plain', 'oxygen') and
                spec['size'] in ('base', 'large') and spec['case'] == run['case'], 'Wrong commissioned Zr case')
        require(spec['h_A'] in (.2, .16) and spec['kmesh'] in ([1, 1, 1], [2, 2, 2]) and
                spec['radii_A'] in (RADII, [.1]), 'Wrong commissioned numerical control')
    else:
        validate_spec(spec, run)
    epoch, sources = _epoch(root, run, required_sources=required_sources)
    positions, cell, masses, direction = geometry_validator(root, run, log_inventory, archive_root)
    g = run['geometry']
    target = g['target_index']
    symbols = g['symbols']
    vacant_symbols = symbols[:target] + symbols[target+1:]
    vacant_positions = np.delete(positions, target, axis=0)
    vacant_masses = np.delete(masses, target)
    states = run['states']
    require(set(states) == {'occupied-host', 'vacant-host', 'occupied-probes', 'vacant-probes',
                            'occupied-probe-He', 'vacant-probe-H-R0'}, 'Incomplete exact state roles')
    energies = {}
    grid = None

    def admit(state, role, wanted_symbols, wanted_positions, wanted_masses, charge):
        nonlocal grid
        energy = _state(root, state, role, spec, g, False, log_inventory, archive_root, valence=valence)
        require(state['symbols'] == wanted_symbols and state['charge'] == charge, role + ': wrong nuclear/charge role')
        close(state['positions_A'], wanted_positions, role + ': frozen geometry/probe role')
        close(state['masses_u'], wanted_masses, role + ': isotope role')
        if grid is None:
            grid = state['realspace_shape']
        require(state['realspace_shape'] == grid, role + ': unequal quartet grids')
        return energy

    for role, ss, pp, mm in [('occupied-host', symbols, positions, masses),
                             ('vacant-host', vacant_symbols, vacant_positions, vacant_masses)]:
        energies[role] = admit(states[role], role, ss, pp, mm, 0)
    green = PeriodicNuclearGreen(cell)
    radii = spec['radii_A']
    require(len(states['occupied-probes']) == len(states['vacant-probes']) == len(radii), 'Missing actual radii')
    rows = []
    for index, radius in enumerate(radii):
        projectile = positions[target] - radius * direction
        pair = {}
        for label, ss, pp, mm in [('occupied', symbols, positions, masses),
                                 ('vacant', vacant_symbols, vacant_positions, vacant_masses)]:
            state = states[label + '-probes'][index]
            close(state['R_A'], radius, 'Actual descending probe radius')
            close(state['projectile_A'], projectile, 'Actual projectile ray')
            pair[label] = admit(state, label + '-probe-' + str(radius), ss + ['H'],
                                np.vstack((pp, projectile)), np.append(mm, masses[target]), 1)
        q = math.fsum([pair['occupied'], -pair['vacant'], -energies['occupied-host'], energies['vacant-host']])
        correction = green.point_minus_gaussian_eV(radius * direction, math.sqrt(5.) / BOHR_A)
        total = q + correction
        rows.append(dict(R_A=radius, BO_quartet_soft_eV=q, point_minus_gaussian_eV=correction,
                         V_pair_total_eV=total, lowering_eV=COULOMB_EV_A/radius-total))
    require(len(run['rows']) == len(rows) and all(set(row) == ROW_KEYS for row in run['rows']), 'Wrong saved row schema')
    compare_json(run['rows'], rows, 'Independent quartet/Green/lowering replay')
    he_symbols = list(symbols)
    he_symbols[target] = 'He'
    he_masses = masses.copy()
    he_masses[target] *= 2
    energies['occupied-probe-He'] = admit(states['occupied-probe-He'], 'occupied-probe-He',
                                         he_symbols, positions, he_masses, 1)
    energies['vacant-probe-H-R0'] = admit(states['vacant-probe-H-R0'], 'vacant-probe-H-R0',
                                         symbols, positions, masses, 1)
    contact = run['contact']
    convention = contact['absolute_convention']
    for key, expected in dict(alpha1=10., alpha2=300., H_E_reference_Hartree=0., He_E_reference_Hartree=0.,
                              representation='same-alpha all-electron Z2', arbitrary_atomic_offset_eV=0.,
                              energy='internal=free-minus-negativeTS').items():
        require(convention[key] == expected, 'Wrong absolute convention ' + key)
    operators = contact['operator_sources']
    require(set(operators) == {'ae', 'setup', 'basis_data', 'hamiltonian', 'density'}, 'Missing actual operator sources')
    opaths = {key: bound_file(root, value, 'Actual operator ' + key) for key, value in operators.items()}
    mapping_path = next(path for name, path in sources.items() if Path(name).name == 'collision_bare_reference.py')
    proof = verify_mapping(cell, convention['H_M_Hartree'], convention['He_M_Hartree'],
                           convention['ASE_Hartree_eV'], opaths['ae'], mapping_path)
    units_path = bound_file(root, convention['ASE_units_source'], 'Actual ASE unit source')
    units = verify_ase_units(units_path, convention['ASE_Hartree_eV'], convention['ASE_CODATA_version'])
    for actual in (operators['ae'], operators['setup'], convention['ASE_units_source']):
        require(any(all(actual[k] == item[k] for k in ('bytes', 'sha256')) for item in spec['library_bindings']),
                'Consumed operator/unit source not bound by original input')
    close(convention['B_native_eV'], proof['B_native_eV'], 'Actual/source Bnative')
    identity_record = contact['operator_identity']
    require(identity_record['schema'] == 'hydride-contact-operator-identity-v1' and
            identity_record['source'] == operators, 'Wrong constructed operator source')
    coincident = identity_record['coincident_native']
    require(coincident['performed'] is False and isinstance(coincident['reason'], str) and coincident['reason'],
            'Unqualified coincident Native claim')
    guards = {key: bound_file(root, value, 'Unmodified public guard') for key, value in coincident['guard_sources'].items()}
    radial = read_json(bound_file(root, identity_record['radial_file'], 'Actual constructed radial operator'))
    for role, records in states.items():
        for state in records if isinstance(records, list) else [records]:
            for symbol in ('H', 'He'):
                if symbol in state['actual_setups']:
                    for key in ('Z', 'Nv', 'Nc', 'ni', 'E', 'alpha1', 'alpha2', 'M', 'Delta0'):
                        close(state['actual_setups'][symbol][key], radial[symbol][key],
                              role + ': actual loaded operator/radial ' + symbol + '/' + key)
    for role in ('occupied-probe-He', 'vacant-probe-H-R0'):
        close(states[role]['R_A'], 0., role + ': exact contact radius')
        close(states[role]['projectile_A'], positions[target], role + ': exact contact location')
    identity = verify_operator_identity(radial, opaths, guards, convention, 1, states['occupied-probe-He']['electrons'])
    close(identity['constructed_operator_difference_eV'], proof['B_native_eV'], 'Common-density Bnative identity')
    raw = math.fsum([energies['occupied-probe-He'], -energies['vacant-probe-H-R0'],
                     -energies['occupied-host'], energies['vacant-host']])
    signed = math.fsum([raw, proof['unit_constant_residual_eV'], proof['mapping']['regular_images_eV']])
    close(contact['Q_He_eV'], raw, 'Independent contact quartet')
    close(contact['signed_W0_eV'], signed, 'Independent signed W0 after outer hash repair')
    for item in run['logs']:
        _log(root, item, log_inventory, archive_root)
    require(number(run['process_cpu_s'], 'Actual CPU') > 0 and isinstance(run['scope'], str) and run['scope'],
            'Missing actual cost/scope')
    require(binding(path) == original, 'Case changed during replay')
    return dict(schema='micro-site-replay-v1', case=run['case'], specification=spec,
                geometry=g, rows=rows, signedW0_eV=signed, Q_He_eV=raw,
                raw_quartet_energies_eV=energies, states=states, source_epoch=run['source_epoch'],
                source_manifest=epoch, case_binding=original, absolute_convention=convention,
                operator_sources=operators, mapping_proof=proof, unit_source_proof=units,
                functional_identity_replay=identity, scope=SCOPE)


def replay(directory, *, log_inventory=None, archive_root=None):
    """Replay the commissioned Zr family with its strict original defaults."""
    return _replay(directory, log_inventory=log_inventory, archive_root=archive_root)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    print(json.dumps(replay(args.directory), indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
