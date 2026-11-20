"""NumPy-only replay of original contact states, never an atomic-offset fit.

A load is an actual integrity/physics check, not a claim made by this source.
The source-bound constructed-operator identity is mandatory before export;
the public coincident-atom SCF guard is never bypassed.
"""
import argparse
import json
import math
import re
from pathlib import Path

import numpy as np

from hydride_replay import (binding, bound_file, check_binding, close, compare_json, load_run,
                            number, read_json, require)
from verify_contact_mapping import verify_ase_units, verify_mapping, verify_operator_identity

ROLES = ('occupied-probe-He', 'vacant-probe-H', 'occupied-host', 'vacant-host')
from log_privacy import LOG_TRANSFORM_COUNTS, sanitize_log


def _log_file(directory, value, inventory, archive_root, label):
    """Raw: verify original bytes. Portable: verify declared ORIGINAL/EXPORTED.

    The export inventory attests the approved User/CWD/library removal and
    dataset/GPW basename transforms; removed original bytes are intentionally
    absent and cannot be reconstructed here.
    Scientific state/source/energy bindings never use this transformation path.
    """
    if inventory is None:
        return bound_file(directory, value, label)
    require(archive_root is not None, 'Portable log inventory requires its archive root')
    check_binding(value, label)
    root = Path(archive_root).resolve()
    name = (Path(directory).resolve().relative_to(root) / value['path']).as_posix()
    require(name in inventory, label + ': missing ORIGINAL/EXPORTED inventory')
    entry = inventory[name]
    require(set(entry) == {'ORIGINAL', 'EXPORTED', 'transformation'}, label + ': wrong inventory schema')
    original = {key: value[key] for key in ('bytes', 'sha256')}
    require(entry['ORIGINAL'] == original, label + ': original scientific log binding changed')
    exported = dict(path=value['path'], **entry['EXPORTED'])
    path = bound_file(directory, exported, label + ' EXPORTED')
    changes = entry['transformation']
    if changes:
        require(path.suffix.lower() in ('.txt', '.log') and type(changes) is list and len(changes) == 1 and
                type(changes[0]) is dict and bool(changes[0]) and set(changes[0]) <= LOG_TRANSFORM_COUNTS and
                all(type(v) is int and v > 0 for v in changes[0].values()) and
                original['bytes'] > entry['EXPORTED']['bytes'] and sanitize_log(path.read_bytes())[1] == [],
                label + ': only declared approved GPAW log identity/path transforms are allowed')
    else:
        require(changes == [] and entry['ORIGINAL'] == entry['EXPORTED'], label + ': undeclared transformation')
    return path


def _method_spec(spec):
    spec = dict(spec)
    spec.pop('absolute_stop_utc', None)
    spec['method'] = dict(spec['method'])
    spec['method'].pop('contact_He_setup', None)
    return spec


def _epoch(directory, epoch):
    require(set(epoch) == {'input', 'source_manifest'}, 'Wrong contact source epoch')
    ipath = bound_file(directory, epoch['input'], 'Contact input')
    mpath = bound_file(directory, epoch['source_manifest'], 'Contact frozen manifest')
    manifest = read_json(mpath)
    require(set(manifest) == {'files'} and isinstance(manifest['files'], list) and manifest['files'],
            'Invalid/empty contact executable epoch')
    names = []
    for item in manifest['files']:
        bound_file(mpath.parent / 'source', item, 'Frozen contact executable')
        names.append(item['path'])
    require(len(names) == len(set(names)), 'Duplicate contact source member')
    for required in ('contact_native.py', 'contact_state.py', 'contact_setup.py'):
        require(sum(Path(name).name == required for name in names) == 1,
                'Missing/ambiguous frozen producer source ' + required)
    return read_json(ipath), manifest


def _state(state, directory, role, spec, operator_sources, convention, log_inventory, archive_root):
    require(state['schema'] == 'hydride-native-state-v1', role + ': wrong raw state schema')
    path = bound_file(directory, state['state_file'], role + ': raw standalone state')
    require(read_json(path) == {k: v for k, v in state.items() if k != 'state_file'},
            role + ': standalone state differs')
    symbols = state['symbols']
    require(isinstance(symbols, list) and set(symbols) <= {'Li', 'H', 'He'}, role + ': wrong species')
    state_spec, actual_spec = _method_spec(state['specification']), _method_spec(spec)
    if state['backend'] == 'analytic-empty-vacuum':
        state_spec['h_A'] = actual_spec['h_A']  # No grid exists for the exact empty Hamiltonian.
    compare_json(state_spec, actual_spec, role + ': physical specification')
    require(type(state['charge']) is int and state['charge'] in (0, 1), role + ': invalid charge')
    expected_ne = sum(2 if s == 'He' else 1 for s in symbols) - state['charge']
    close(number(state['electrons'], role + ' Ne'), expected_ne, role + ': actual Ne / Z / charge')
    energy = state['energy']
    require(set(energy) == {'internal_energy_eV', 'free_energy_eV', 'minus_TS_eV'}, role + ': energy schema')
    for key, value in energy.items():
        number(value, role + ':' + key)
    close(energy['internal_energy_eV'], energy['free_energy_eV'] - energy['minus_TS_eV'],
          role + ': internal=free-minus-negativeTS')
    require(energy['minus_TS_eV'] <= 1e-9, role + ': invalid negative TS')
    if state['backend'] == 'analytic-empty-vacuum':
        require(role == 'vacant-host' and spec['case'] == 'G' and symbols == [] and
                state['positions_A'] == [] and state['charge'] == 0 and state['electrons'] == 0 and
                state['scf_converged'] is False and all(v == 0 for v in energy.values()) and
                state['realspace_shape'] is None and state['actual_kpoints'] is None and state['nbands'] == 0,
                'Analytic empty energy permitted only for exact G empty host')
    else:
        require(state['backend'] == 'gpaw-native' and state['scf_converged'] is True and
                state['xc'] == 'PBE' and state['spinpol'] is False, role + ': actual converged nonspin PBE required')
        close(state['occupation_width_eV'], .1, role + ': actual smearing width')
        shape = state['realspace_shape']
        require(isinstance(shape, list) and len(shape) == 3 and all(type(v) is int and v > 0 for v in shape),
                role + ': invalid actual realspace grid')
        require('gpw' in state or 'log' in state, role + ': no raw Native recovery/log provenance')
        if 'log' in state:
            _log_file(directory, state['log'], log_inventory, archive_root, role + ': original Native log')
        require(type(state['actual_kpoints']) is int and state['actual_kpoints'] == 1,
                role + ': actual Gamma sampling required')
        require(type(state['nbands']) is int and state['nbands'] > 0 and 2 * state['nbands'] >= expected_ne,
                role + ': invalid actual bands')
        require(type(state['scf_iterations']) is int and state['scf_iterations'] > 0,
                role + ': missing actual SCF iterations')
        compare_json(state['setup_valence'], {s: float(2 if s == 'He' else 1) for s in set(symbols)},
                     role + ': actual setup valence')
        setups = state['setups']
        require(setups.get('H') == 'ae' and set(setups) == ({'H', 'He'} if 'He' in symbols else {'H'}),
                role + ': wrong Hae/custom-He setup family')
        if 'He' in symbols:
            he = setups['He']
            for key, expected in dict(name='HaeZ2', type='all-electron', Z=2, Nv=2, Nc=0, ni=0,
                                      E=0., alpha1=10., alpha2=300., occupation=[2.]).items():
                require(he[key] == expected and (type(expected) is not int or type(he[key]) is int),
                        role + ': wrong custom He ' + key)
            close(he['Delta0'], -2 / math.sqrt(4 * math.pi), role + ': He charge coefficient')
            close(he['M'], convention['He_M_Hartree'], role + ': He self reference')
            require(he['source'] == operator_sources, role + ': changed loaded operator source')
    return energy['internal_energy_eV']


def load_contact(directory, *, log_inventory=None, archive_root=None):
    """Validate saved raw/source bytes and return a recomputed signed contact anchor.

    Actual constructed-operator arrays and source-functional equality replace
    the prohibited coincident-two-H SCF, without claiming that SCF occurred.
    """
    directory = Path(directory)
    path = directory / 'contact-native.json'
    original = binding(path)
    run = read_json(path)
    require(run['schema'] == 'hydride-contact-native-v1', 'Wrong contact Native schema')
    case, q, h = run['case'], run['charge'], run['h_A']
    require(case in ('F', 'G') and type(q) is int and q in (0, 1) and type(h) is float and h in (.2, .16),
            'Wrong contact case/charge/registered grid')
    require(run['base_role'] == ('F_h' if case == 'F' and h == .16 else case), 'Wrong finite-case base role')
    inp, manifest = _epoch(directory, run['source_epoch'])
    require(inp == dict(schema='hydride-contact-input-v1', case=case, charge=q, h_A=h,
                        base_role=run['base_role'], absolute_stop_utc='2026-11-20T18:00:00Z'),
            'Frozen contact input differs from actual run')
    refs = run['source_references']
    expected_refs = {'input', 'source_manifest', 'run_result', 'quartet', 'geometry', 'occupied_host', 'vacant_host'}
    if case == 'G' and h == .16 and q == 0:
        expected_refs.add('prior_q1')
    require(set(refs) == expected_refs, 'Wrong original source reference set')
    paths = {key: bound_file(directory, value, 'Original ' + key) for key, value in refs.items()}
    base = load_run(paths['run_result'].parent, case)
    require(paths['quartet'] == paths['run_result'].parent / base['quartet_binding']['path'] and
            paths['geometry'] == paths['run_result'].parent / base['geometry_binding']['path'],
            'Original geometry/quartet binding is not the replayed base')
    require(paths['input'] == paths['run_result'].parent.parent / base['run']['source_epoch']['input']['path'] and
            paths['source_manifest'] == paths['run_result'].parent.parent / base['run']['source_epoch']['source_manifest']['path'],
            'Original executable/input source mismatch')
    for role in ('occupied-host', 'vacant-host'):
        expected = base['quartet']['states'][role]['state_file']
        require(paths[role.replace('-', '_')] == paths['run_result'].parent / 'bo' / expected['path'],
                role + ': original host source mismatch')
    geometry, spec = run['geometry'], run['specification']
    require(geometry == base['run']['geometry'], 'Contact changed finite-R geometry, ray, target or spectators')
    expected_spec = _method_spec(base['run']['specification'])
    expected_spec['h_A'] = h
    compare_json(_method_spec(spec), expected_spec, 'Registered contact physical specification')
    require(spec['method']['contact_He_setup'] == 'HaeZ2', 'Wrong actual contact He operator')
    operators = run['operator_sources']
    require(set(operators) == {'ae', 'setup', 'basis_data', 'hamiltonian', 'density'},
            'Missing loaded operator/functional-assembly source bytes')
    opaths = {key: bound_file(directory, value, 'Loaded operator ' + key) for key, value in operators.items()}
    convention = run['absolute_convention']
    for key, expected in dict(alpha1=10., alpha2=300., H_E_reference_Hartree=0., He_E_reference_Hartree=0.,
                              representation='same-alpha all-electron Z2').items():
        require(convention[key] == expected, 'Wrong absolute operator convention: ' + key)
    mapping_path = paths['source_manifest'].parent / 'source' / 'substrate_framework/collision_bare_reference.py'
    proof = verify_mapping(geometry['cell_A'], convention['H_M_Hartree'], convention['He_M_Hartree'],
                           convention['ASE_Hartree_eV'], opaths['ae'], mapping_path)
    units_path = bound_file(directory, convention['ASE_units_source'], 'Actual ASE units source')
    proof['ASE_units_proof'] = verify_ase_units(units_path, convention['ASE_Hartree_eV'],
                                               convention['ASE_CODATA_version'])
    require(convention['arbitrary_atomic_offset_eV'] == 0. and
            convention['energy'] == 'internal=free-minus-negativeTS', 'Fitted offset/wrong energy convention')
    close(convention['B_native_eV'], proof['B_native_eV'], 'Declared/source B_native')
    states = run['states']
    require(set(states) == set(ROLES), 'Incomplete/unregistered contact quartet; no coincident-H SCF allowed')
    energies = {role: _state(state, directory, role, spec, operators, convention, log_inventory, archive_root)
                for role, state in states.items()}
    positions = np.asarray(geometry['positions_A'], dtype=float)
    target = geometry['target_index']
    symbols = geometry['symbols']
    n = len(symbols)
    desired_symbols = dict(zip(ROLES, [symbols[:target] + ['He'] + symbols[target + 1:], symbols,
                                      symbols, symbols[:target] + symbols[target + 1:]]))
    desired_positions = dict(zip(ROLES, [positions.tolist(), positions.tolist(), positions.tolist(),
                                        np.delete(positions, target, axis=0).tolist()]))
    grid = states['occupied-probe-He']['realspace_shape']
    for role in ROLES:
        state = states[role]
        require(state['symbols'] == desired_symbols[role], role + ': target-only united-atom inventory')
        close(state['positions_A'], desired_positions[role], role + ': frozen target/spectator coordinates')
        close(state['cell_A'], geometry['cell_A'], role + ': frozen cell')
        expected_q = q if role in ROLES[:2] else 0
        require(state['charge'] == expected_q, role + ': wrong branch/host charge')
        close(state['electrons'], (n + 1 - q if role == ROLES[0] else n - q if role == ROLES[1]
                                  else n if role == ROLES[2] else n - 1), role + ': branch Ne')
        require(state['backend'] == 'analytic-empty-vacuum' or state['realspace_shape'] == grid,
                role + ': different quartet grid')
    for role in ('occupied-host', 'vacant-host'):
        if not (case == 'G' and h == .16):
            require(read_json(paths[role.replace('-', '_')]) ==
                    {key: value for key, value in states[role].items() if key != 'state_file'},
                    role + ': original host was not reused byte-exact in content')
            require(refs[role.replace('-', '_')]['sha256'] == states[role]['state_file']['sha256'],
                    role + ': original host bytes changed')
    if 'prior_q1' in paths:
        prior = read_json(paths['prior_q1'])
        require(prior['schema'] == 'hydride-contact-native-v1' and prior['case'] == 'G' and
                prior['charge'] == 1 and prior['h_A'] == .16 and prior['geometry'] == geometry and
                prior['specification'] == spec, 'Wrong source-bound prior q1 fine host')
        require(prior['states']['occupied-host']['state_file']['sha256'] ==
                states['occupied-host']['state_file']['sha256'], 'Prior q1 actual fine host bytes changed')
    if q == 0:
        require(states['vacant-probe-H']['state_file']['sha256'] == states['occupied-host']['state_file']['sha256'],
                'q0 vacant probe must reuse identical full-H neutral host bytes')
    identity_record = run['operator_identity']
    require(identity_record['schema'] == 'hydride-contact-operator-identity-v1' and
            identity_record['source'] == operators, 'Missing/changed actual-operator proof source')
    coincident = identity_record['coincident_native']
    require(coincident['performed'] is False and isinstance(coincident['reason'], str) and coincident['reason'],
            'Prohibited coincident SCF must be honestly recorded as not performed')
    guards = {key: bound_file(directory, value, 'Public guard ' + key)
              for key, value in coincident['guard_sources'].items()}
    radial_path = bound_file(directory, identity_record['radial_file'], 'Actual constructed radial operators')
    identity = verify_operator_identity(read_json(radial_path), opaths, guards, convention, q,
                                        states['occupied-probe-He']['electrons'])
    identity['radial_file'] = identity_record['radial_file']
    identity['coincident_native'] = coincident
    close(identity['constructed_operator_difference_eV'], proof['B_native_eV'],
          'Independent same-density source functional identity')
    raw = math.fsum([energies[ROLES[0]], -energies[ROLES[1]], -energies[ROLES[2]], energies[ROLES[3]]])
    mapping = proof['mapping']
    signed = math.fsum([raw, proof['unit_constant_residual_eV'], mapping['regular_images_eV']])
    require(binding(path) == original, 'Contact Native result changed during replay')
    return dict(schema='hydride-contact-replay-v1', case=case, charge=q, h_A=h, R_A=0.,
                signedW0_eV=signed, Q_He_eV=raw, raw_quartet_energies_eV={r: energies[r] for r in ROLES},
                electrons={r: states[r]['electrons'] for r in ROLES}, geometry=geometry,
                specification=spec, absolute_convention=convention, periodic_mapping=mapping,
                B_native_eV=proof['B_native_eV'], unit_constant_residual_eV=proof['unit_constant_residual_eV'],
                mapping_proof=proof, gas_absolute_identity=identity,
                raw_states=states, native_binding=original, source_epoch=run['source_epoch'],
                source_references=refs, operator_sources=operators, source_manifest=manifest,
                finite_base_run=base['run_binding'], finite_base_sources=base['sources'],
                finite_base_rows=base['rows']['q' + str(q)])


def qualify_pair(coarse, fine):
    require(coarse['case'] == fine['case'] and coarse['charge'] == fine['charge'] and
            coarse['h_A'] == .2 and fine['h_A'] == .16, 'Wrong contact h pair')
    for field in ('case', 'symbols', 'positions_A', 'cell_A', 'masses_u',
                  'target_index', 'path', 'volume_source'):
        compare_json(coarse['geometry'][field], fine['geometry'][field],
                     'Contact h pair physical geometry/ray: ' + field)
    cs, fs = dict(coarse['specification']), dict(fine['specification'])
    cs['h_A'] = fs['h_A']
    compare_json(cs, fs, 'Contact h pair method/charge family')
    require(coarse['source_manifest'] == fine['source_manifest'] and
            coarse['operator_sources'] == fine['operator_sources'], 'Contact h pair changed executable/operator epoch')
    cg = np.asarray(coarse['raw_states']['occupied-probe-He']['realspace_shape'])
    fg = np.asarray(fine['raw_states']['occupied-probe-He']['realspace_shape'])
    require(bool(np.all(fg >= cg) and np.any(fg > cg)), 'Contact h control did not refine the actual grid')
    delta = coarse['signedW0_eV'] - fine['signedW0_eV']
    tolerance = max(5., .1 * abs(fine['signedW0_eV']))
    require(abs(delta) <= tolerance, 'Contact h0.20/0.16 qualification FAILED')
    return dict(coarse_h_A=.2, fine_h_A=.16, coarse_signedW0_eV=coarse['signedW0_eV'],
                fine_signedW0_eV=fine['signedW0_eV'], coarse_minus_fine_eV=delta,
                absolute_delta_eV=abs(delta), tolerance_eV=tolerance, qualified=True,
                scope='Same frozen geometry finite-cell grid sensitivity only; not exact h->0 or isolated-cell convergence')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    args = parser.parse_args()
    print(json.dumps(load_contact(args.directory), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
