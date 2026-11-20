"""NumPy-only integrity and physical replay of saved LiD/gas Native quartets.

No electronic-structure package is imported and no missing calculation is inferred.
GPW descriptors are retained as provenance; portable replay needs standalone JSON,
not the private wavefunctions. All derived comparisons use atol=1e-9, rtol=0.
"""
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

_SUPPORT = Path(__file__).resolve().parent
if not (_SUPPORT / 'green.py').is_file():
    _SUPPORT = _SUPPORT.parent / 'delivery'
sys.path.insert(0, str(_SUPPORT))
from green import BOHR_A, COULOMB_EV_A, PeriodicNuclearGreen

RADII = [.1, .2, .3, .5, 1., 1.6]
VOLUME_SHA = '01e4b1b4022b9d540a581f5e270262c1177a6b43f946f97a783cf71444830d4e'
ROW_KEYS = {'R_A', 'BO_quartet_soft_eV', 'point_minus_gaussian_eV',
            'V_pair_total_eV', 'lowering_eV'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'Duplicate JSON key: ' + key)
        result[key] = value
    return result


def read_json(path):
    data = json.loads(Path(path).read_bytes(), object_pairs_hook=_pairs)
    json.dumps(data, allow_nan=False)
    return data


def binding(path, name=None):
    path = Path(path)
    data = path.read_bytes()
    return dict(path=path.name if name is None else name, bytes=len(data),
                sha256=hashlib.sha256(data).hexdigest())


def check_binding(value, label, basename=False):
    require(isinstance(value, dict) and set(value) == {'path', 'bytes', 'sha256'},
            label + ': require exact path/bytes/sha256 binding')
    name = value['path']
    require(isinstance(name, str) and name and '\\' not in name and
            not Path(name).is_absolute() and all(p not in ('', '.', '..') for p in name.split('/')),
            label + ': unsafe relative source path')
    require(not basename or Path(name).name == name, label + ': require basename')
    require(type(value['bytes']) is int and value['bytes'] >= 0, label + ': invalid bytes')
    sha = value['sha256']
    require(isinstance(sha, str) and len(sha) == 64 and
            all(c in '0123456789abcdef' for c in sha), label + ': invalid SHA256')


def bound_file(root, value, label, basename=False):
    check_binding(value, label, basename)
    root = Path(root).resolve()
    path = root / value['path']
    require(path.resolve().is_relative_to(root) and path.is_file() and
            not any(p.is_symlink() for p in (path, *path.parents) if p != root.parent),
            label + ': missing/unsafe source')
    require(binding(path, value['path']) == value, label + ': source byte/hash mismatch')
    return path


def number(value, label):
    require(type(value) in (int, float) and math.isfinite(value), label + ': finite number required')
    return value


def close(actual, expected, label):
    a, b = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(),
            label + ': nonfinite/shape mismatch')
    require(bool(np.all(np.abs(a - b) <= 1e-9)), label + ': changed values (atol=1e-9, rtol=0)')


def compare_json(actual, expected, label='JSON replay'):
    """Compare exact structure/discrete metadata, tolerating only float roundoff."""
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and set(actual) == set(expected), label + ': changed keys')
        for key in expected:
            compare_json(actual[key], expected[key], label + '.' + key)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), label + ': changed list')
        for index, (a, b) in enumerate(zip(actual, expected)):
            compare_json(a, b, label + '[' + str(index) + ']')
    elif type(expected) is float:
        require(type(actual) is float, label + ': changed float type')
        close(number(actual, label), expected, label)
    else:
        require(type(actual) is type(expected) and actual == expected, label + ': changed discrete metadata')


def validate_epoch(directory, run):
    epoch = run['source_epoch']
    require(set(epoch) == {'input', 'source_manifest'}, 'Wrong source epoch fields')
    unit = directory.parent
    ipath = bound_file(unit, epoch['input'], 'Native input', basename=True)
    mpath = bound_file(unit, epoch['source_manifest'], 'Source manifest', basename=True)
    require(ipath.name == 'input.json' and mpath.name == 'source-manifest.json', 'Wrong source epoch names')
    require(read_json(ipath) == run['specification'], 'Input specification changed')
    spec = run['specification']
    original_volume = None
    if spec['case'] == 'F':
        original = read_json(bound_file(unit, spec['volume_source'], 'Original PW relaxed-volume primary', basename=True))
        close(spec['primitive_cell_A'], original['cell_A'], 'Original PW primitive cell')
        close(spec['primitive_positions_A'], original['positions_A'], 'Original PW primitive positions')
        require(spec['primitive_symbols'] == original['initial_geometry']['symbols'] == ['Li', 'H'],
                'Original PW primitive symbols')
        require(original['relaxation']['relaxed'] is True and
                original['settings']['ecut_eV'] == 500 and original['settings']['kmesh'] == [4, 4, 4],
                'Original volume must be the relaxed PW500/k4 primary')
        original_volume = dict(source=spec['volume_source'], volume_A3=original['volume_A3'],
                               settings=original['settings'], relaxation=original['relaxation'],
                               stress_eV_A3=original['stress_eV_A3'],
                               scope='Original PW500/k4 hydrostatic volume/stress; not FD volume reoptimization')
    manifest = read_json(mpath)
    require(set(manifest) == {'files'}, 'Wrong frozen source manifest schema')
    files = manifest['files']
    require(isinstance(files, list) and files, 'Empty frozen executable source epoch')
    names = []
    for item in files:
        bound_file(unit / 'source', item, 'Frozen executable source')
        names.append(item['path'])
    require(len(names) == len(set(names)), 'Duplicate frozen source member')
    return manifest, original_volume


def validate_geometry(geometry, spec):
    case = spec['case']
    require(geometry['case'] == case and geometry['passed'] is True, 'Unqualified/wrong geometry')
    require(0 <= number(geometry['raw_fmax_eV_A'], 'Geometry fmax') <= .1,
            'Geometry must have actual forces <=0.1 eV/A')
    require(type(geometry['steps']) is int and 0 <= geometry['steps'] <= spec['ionic_steps'],
            'Invalid ionic relaxation step count')
    symbols = geometry['symbols']
    require(isinstance(symbols, list) and all(s in ('Li', 'H') for s in symbols), 'Wrong nuclear species')
    positions = np.asarray(geometry['positions_A'], dtype=float)
    require(positions.shape == (len(symbols), 3) and np.isfinite(positions).all(), 'Invalid geometry positions')
    cell = np.asarray(geometry['cell_A'], dtype=float)
    PeriodicNuclearGreen(cell)
    target = geometry['target_index']
    require(type(target) is int and 0 <= target < len(symbols) and symbols[target] == 'H', 'Wrong target nucleus')
    path = geometry['path']
    close(path['target_A'], positions[target], 'Actual relaxed target')
    direction = np.asarray(path['direction'], dtype=float)
    require(direction.shape == (3,), 'Invalid approach direction')
    close(np.linalg.norm(direction), 1., 'Unit approach direction')
    clearance = path['minimum_background_clearance_A']
    if case == 'G':
        require(clearance is None, 'One-atom gas has no background clearance')
    else:
        require(number(clearance, 'Background clearance') > 0,
                'Projectile ray intersects background nucleus')
    if case == 'F':
        repeat = spec['repeat_primitive']
        require(type(repeat) is int and repeat in (2, 3), 'F requires primitive repeat2 or repeat3')
        require(symbols.count('Li') == symbols.count('H') == repeat ** 3, 'Wrong LiD inventory')
        primitive = np.asarray(spec['primitive_cell_A'], dtype=float)
        require(primitive.shape == (3, 3), 'Missing original primitive cell')
        close(cell, repeat * primitive, 'Fixed inherited PW primitive volume')
        original_positions = np.asarray(spec['primitive_positions_A'], dtype=float)
        require(original_positions.shape == (2, 3) and np.isfinite(original_positions).all(),
                'Missing original primitive positions')
        require(geometry['volume_source'] == spec['volume_source'], 'Changed original PW volume provenance')
        check_binding(spec['volume_source'], 'Original PW volume', basename=True)
        require(spec['volume_source']['sha256'] == VOLUME_SHA and
                spec['volume_source']['bytes'] == 3075 and
                spec['primitive_symbols'] == ['Li', 'H'], 'Wrong original PW volume source')
    else:
        require(case in ('G', 'G2') and spec['repeat_primitive'] == 1, 'Wrong gas case')
        require(symbols == ['H'] * (1 if case == 'G' else 2) and target == 0, 'Wrong gas target inventory')
        close(spec['box_A'], 10., 'Gas box specification')
        close(cell, 10. * np.eye(3), 'Gas 10 A periodic box')
        close(direction, np.ones(3) / math.sqrt(3.), 'Gas [111] approach')
        require(geometry['volume_source'] is None and spec['volume_source'] is None,
                'Gas has no inherited solid volume')
    return positions, cell


def validate_spec(spec, case):
    require(spec['schema'] == 'hydride-specification-v1' and spec['case'] == case,
            'Mismatched case specification/schema')
    require(type(spec['repeat_primitive']) is int, 'Noninteger primitive repeat')
    require(spec['kmesh'] == [1, 1, 1] and all(type(n) is int for n in spec['kmesh']),
            'Require discrete Gamma mesh')
    close(spec['force_tolerance_eV_A'], .1, 'Registered force threshold')
    require(type(spec['ionic_steps']) is int and spec['ionic_steps'] == 20, 'Wrong ionic step limit')
    require(number(spec['h_A'], 'Grid spacing') in (.2, .16), 'Unregistered grid spacing')
    require(spec['absolute_stop_utc'] == '2026-11-20T18:00:00Z', 'Wrong absolute stop provenance')
    compare_json(spec['method'], dict(xc='PBE', spinpol=False, occupation_width_eV=.1,
                 geometry_H_setup='standard_PAW', BO_H_setup='ae', Li_setup='default_1valence',
                 Hae_alpha1_Bohr_m2=10, energy_convention='internal=free-minus-negativeTS'),
                 'Registered method')
    expected_padding = 2 * spec['repeat_primitive'] ** 3 if case == 'F' else 4
    require(type(spec['nband_padding']) is int and spec['nband_padding'] == expected_padding,
            'Wrong band-padding specification')


def state_energy(state, directory, spec, role, allow_empty=False):
    require(state['schema'] == 'hydride-native-state-v1', role + ': wrong Native state schema')
    require(state['specification'] == spec, role + ': changed state specification')
    spath = bound_file(directory / 'bo', state['state_file'], role + ' standalone JSON', basename=True)
    require(read_json(spath) == {k: v for k, v in state.items() if k != 'state_file'},
            role + ': standalone state differs from quartet')
    charge = state['charge']
    require(type(charge) is int and charge in (0, 1), role + ': wrong discrete charge')
    symbols = state['symbols']
    require(isinstance(symbols, list) and all(s in ('Li', 'H') for s in symbols), role + ': wrong species')
    close(state['electrons'], len(symbols) - charge, role + ': actual electron inventory (Li/H each one valence electron)')
    e = state['energy']
    require(set(e) == {'internal_energy_eV', 'free_energy_eV', 'minus_TS_eV'}, role + ': wrong energy schema')
    for key, value in e.items():
        number(value, role + ' ' + key)
    close(e['internal_energy_eV'], e['free_energy_eV'] - e['minus_TS_eV'], role + ': internal energy convention')
    require(e['minus_TS_eV'] <= 1e-9, role + ': minus_TS must be negative TS')
    if state['backend'] == 'analytic-empty-vacuum':
        require(allow_empty and role == 'vacant-host' and spec['case'] == 'G' and
                symbols == [] and state['positions_A'] == [] and charge == 0 and
                state['electrons'] == 0 and state['scf_converged'] is False and
                all(value == 0 for value in e.values()) and 'gpw' not in state and
                state['realspace_shape'] is None and state['actual_kpoints'] is None and
                type(state['nbands']) is int and state['nbands'] == 0,
                'Analytic zero allowed ONLY for one-atom G vacant-host, never a failed/missing Native state')
    else:
        require(state['backend'] == 'gpaw-native' and state['scf_converged'] is True and
                state['xc'] == 'PBE' and state['setups'] == {'H': 'ae'}, role + ': require converged PBE/Hae Native')
        require(state['spinpol'] is False, role + ': require nonspin Native')
        close(state['occupation_width_eV'], .1, role + ': actual FD width')
        require(state['setup_valence'] == {symbol: 1 for symbol in set(symbols)} and
                all(type(value) in (int, float) for value in state['setup_valence'].values()),
                role + ': actual one-valence Li/H setups')
        check_binding(state['gpw'], role + ' private GPW', basename=True)
        shape = state['realspace_shape']
        require(isinstance(shape, list) and len(shape) == 3 and
                all(type(n) is int and n > 0 for n in shape), role + ': invalid actual grid')
        require(type(state['actual_kpoints']) is int and state['actual_kpoints'] == 1,
                role + ': actual Gamma k-point count')
        require(type(state['nbands']) is int and state['nbands'] > 0, role + ': invalid bands')
        require(2 * state['nbands'] >= state['electrons'], role + ': insufficient nonspin bands for electron inventory')
        require(type(state['scf_iterations']) is int and state['scf_iterations'] > 0, role + ': missing SCF iterations')
        require('wavefunction_predictor' in state, role + ': missing predictor provenance')
    return e['internal_energy_eV']


def load_run(directory, expected_case=None):
    """Replay a complete saved run (six radii, or R0.1 for either F control)."""
    directory = Path(directory)
    runpath = directory / 'run-result.json'
    run_binding = binding(runpath)
    run = read_json(runpath)
    require(run['schema'] == 'hydride-run-v1', 'Wrong hydride run schema')
    spec = run['specification']
    case = spec['case']
    require(case in ('F', 'G', 'G2') and (expected_case is None or case == expected_case), 'Wrong expected hydride case')
    validate_spec(spec, case)
    manifest, original_volume = validate_epoch(directory, run)
    require(run['geometry_file']['path'] == 'geometry/relaxation.json' and
            run['quartet']['path'] == 'bo/quartet.json', 'Wrong geometry/quartet source names')
    gpath = bound_file(directory, run['geometry_file'], 'Actual geometry')
    qpath = bound_file(directory, run['quartet'], 'Actual quartet')
    geometry, quartet = read_json(gpath), read_json(qpath)
    require(geometry == run['geometry'] == quartet['geometry'], 'Geometry source records disagree')
    require(quartet['schema'] == 'hydride-quartet-v1' and quartet['case'] == case and
            quartet['specification'] == spec, 'Wrong quartet specification/schema')
    positions, cell = validate_geometry(geometry, spec)
    states = quartet['states']
    require(set(states) == {'occupied-host', 'vacant-host', 'q0', 'q1'} and
            set(quartet['rows']) == {'q0', 'q1'}, 'Wrong quartet charge/role branches')
    occupied, vacant = states['occupied-host'], states['vacant-host']
    eo = state_energy(occupied, directory, spec, 'occupied-host')
    ev = state_energy(vacant, directory, spec, 'vacant-host', allow_empty=case == 'G')
    require(occupied['charge'] == vacant['charge'] == 0, 'Both hosts must be neutral')
    require(occupied['symbols'] == geometry['symbols'], 'Occupied host geometry inventory changed')
    close(occupied['positions_A'], positions, 'Frozen occupied-host positions')
    target = geometry['target_index']
    require(vacant['symbols'] == occupied['symbols'][:target] + occupied['symbols'][target + 1:], 'Target-only removal')
    remaining = np.delete(positions, target, axis=0).tolist()
    close(vacant['positions_A'], remaining, 'Frozen vacant-host positions')
    close(occupied['electrons'] - vacant['electrons'], 1., 'Target removal must remove one electron')
    for state in (occupied, vacant):
        close(state['cell_A'], cell, 'Same host cell')
    require(vacant['backend'] == 'analytic-empty-vacuum' or
            vacant['realspace_shape'] == occupied['realspace_shape'], 'Same host realspace grid')
    green = PeriodicNuclearGreen(cell)
    direction = np.asarray(geometry['path']['direction'], dtype=float)
    rows, state_bindings = {}, {'occupied-host': occupied['state_file'], 'vacant-host': vacant['state_file']}
    for charge in (0, 1):
        branch = 'q' + str(charge)
        branch_states = states[branch]
        require(set(branch_states) == {'occupied-probes', 'vacant-probes'}, branch + ': wrong probe roles')
        saved_rows = quartet['rows'][branch]
        require(isinstance(saved_rows, list) and len(saved_rows) in (1, 6), branch + ': incomplete support')
        require(len(branch_states['occupied-probes']) == len(branch_states['vacant-probes']) == len(saved_rows),
                branch + ': unpaired probe support')
        replayed = []
        for o, v, saved in zip(branch_states['occupied-probes'], branch_states['vacant-probes'], saved_rows):
            require(set(saved) == ROW_KEYS, branch + ': wrong row schema')
            radius = number(saved['R_A'], branch + ' radius')
            require(radius > 0 and o['R_A'] == v['R_A'] == radius, branch + ': mismatched radius')
            projectile = np.asarray(geometry['path']['target_A']) - radius * direction
            energies = []
            for role, probe, host in (('occupied', o, occupied), ('vacant', v, vacant)):
                label = branch + '-' + role + '-probe'
                # Check the charge-specific invariant before the generic inventory guard.
                require(probe['charge'] == charge and type(probe['charge']) is int, label + ': wrong probe charge')
                close(probe['electrons'], host['electrons'] + 1 - charge,
                      'Probe electron count ' + label + (' (neutral adds one)' if charge == 0 else ' (bare adds none)'))
                require(probe['symbols'] == host['symbols'] + ['H'], label + ': wrong projectile inventory')
                close(probe['positions_A'][:-1], host['positions_A'], label + ': frozen background')
                close(probe['positions_A'][-1], projectile, label + ': target-relative position')
                close(probe['projectile_A'], projectile, label + ': projectile annotation')
                close(probe['cell_A'], cell, label + ': same quartet cell')
                require(probe['realspace_shape'] == occupied['realspace_shape'], label + ': changed quartet grid')
                energies.append(state_energy(probe, directory, spec, label))
                key = label + '-R' + format(radius, '.3f')
                require(key not in state_bindings, 'Duplicate state radius')
                state_bindings[key] = probe['state_file']
            raw = math.fsum([energies[0], -energies[1], -eo, ev])
            correction = green.point_minus_gaussian_eV(-radius * direction, math.sqrt(5.) / BOHR_A)
            total = raw + correction
            row = dict(R_A=radius, BO_quartet_soft_eV=raw, point_minus_gaussian_eV=correction,
                       V_pair_total_eV=total, lowering_eV=COULOMB_EV_A / radius - total)
            for key, value in row.items():
                close(number(saved[key], branch + ' ' + key), value, branch + ' saved quartet ' + key)
            replayed.append(row)
        replayed.sort(key=lambda row: row['R_A'])
        require([r['R_A'] for r in replayed] in ([.1], RADII), branch + ': wrong sorted registered support')
        rows[branch] = replayed
    require([r['R_A'] for r in rows['q0']] == [r['R_A'] for r in rows['q1']], 'Charge branches have different support')
    require(len(state_bindings) == len({v['path'] for v in state_bindings.values()}), 'Aliased standalone Native states')
    require(binding(runpath) == run_binding, 'Run changed during replay')
    bound_file(directory, run['quartet'], 'Quartet after replay')
    return dict(run=run, quartet=quartet, rows=rows, source_manifest=manifest,
                original_volume_provenance=original_volume,
                run_binding=run_binding, geometry_binding=run['geometry_file'],
                quartet_binding=run['quartet'], state_bindings=state_bindings,
                sources=dict(run=run_binding, geometry=run['geometry_file'], quartet=run['quartet'],
                             source_epoch=run['source_epoch'], states=state_bindings))
