"""Replay a hash-bound hydride draft from saved outputs with NumPy and stdlib only."""
import argparse
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(1, str(Path(__file__).resolve().parent.parent))
from delivery_io import (SCHEMA, ROLES, VOLUME, byte_binding, digest, json_bytes,
                         json_value, require, safe_name, same_json, scan_text)
from log_privacy import allowed_log_transform, sanitize_log


def verify_inventory(root):
    root = Path(root)
    require(root.is_dir() and not root.is_symlink(), 'Capsule must be a real directory')
    manifest = json_value(root / 'manifest.json')
    require(manifest['schema'] == SCHEMA and manifest['publishing'] is False and
            type(manifest['native_invocations']) is int and manifest['native_invocations'] == 0,
            'Wrong/private-draft capsule declaration')
    omissions = manifest['omissions']
    require(type(omissions) is list and omissions and
            all(isinstance(item, str) and item.strip() for item in omissions), 'Nonempty omissions required')
    files = manifest['files']
    require(type(files) is dict and files and 'manifest.json' not in files, 'Invalid file inventory')
    require({'readout.py', 'export_pairs.py', 'hydride_replay.py', 'green.py', 'reader.py',
             'delivery_io.py', 'log_privacy.py', 'portability_proof.py', 'capsule.py', 'README.md', 'inputs.json',
             } <= set(files), 'Incomplete portable scientific source or saved readout')
    actual = set()
    for path in root.rglob('*'):
        require(not path.is_symlink() and (path.is_file() or path.is_dir()), 'Link/special file in capsule')
        if path.is_file():
            name = path.relative_to(root).as_posix()
            safe_name(name)
            actual.add(name)
    require(actual == set(files) | {'manifest.json'}, 'Unmanifested or missing capsule files')
    for name, entry in files.items():
        safe_name(name)
        require(set(entry) == {'ORIGINAL', 'EXPORTED', 'transformation'}, 'Wrong inventory entry shape')
        byte_binding(entry['ORIGINAL'])
        byte_binding(entry['EXPORTED'])
        data = (root / name).read_bytes()
        require(digest(data) == entry['EXPORTED'], 'EXPORTED byte integrity failed: ' + name)
        changes = entry['transformation']
        require(type(changes) is list, 'Missing transformation declaration')
        if changes:
            require(Path(name).suffix.lower() in {'.txt', '.log', '.stdout', '.stderr', '.exit'} and
                    allowed_log_transform(changes) and
                    entry['ORIGINAL']['bytes'] > entry['EXPORTED']['bytes'] and
                    sanitize_log(data)[1] == [],
                    'Only approved GPAW log identity/location metadata may be transformed: ' + name)
        else:
            require(entry['ORIGINAL'] == entry['EXPORTED'], 'Undeclared transformation: ' + name)
        scan_text(name, data)
        if name.endswith('.json'):
            json_value(root / name)
    scan_text('manifest.json', (root / 'manifest.json').read_bytes())
    return manifest


def bind(root, manifest, binding, name, expected_path):
    safe_name(name)
    require(type(binding) is dict and set(binding) == {'path', 'bytes', 'sha256'} and
            binding['path'] == expected_path, 'Wrong bound source path: ' + name)
    expected = {key: binding[key] for key in ('bytes', 'sha256')}
    byte_binding(expected)
    require(name in manifest['files'] and expected == manifest['files'][name]['ORIGINAL'] and
            expected == manifest['files'][name]['EXPORTED'] and not manifest['files'][name]['transformation'],
            'ORIGINAL scientific source binding mismatch: ' + name)
    require(digest((root / name).read_bytes()) == expected, 'Scientific source bytes changed: ' + name)


def verify_native(root, manifest, role):
    prefix = 'native/' + role + '/'
    calculation = prefix + 'calculation/'
    directory = root / calculation
    run = json_value(directory / 'run-result.json')
    require(run['schema'] == 'hydride-run-v1', 'Wrong Native run schema')
    bind(root, manifest, run['quartet'], calculation + 'bo/quartet.json', 'bo/quartet.json')
    bind(root, manifest, run['geometry_file'], calculation + 'geometry/relaxation.json', 'geometry/relaxation.json')
    for key, filename in (('input', 'input.json'), ('source_manifest', 'source-manifest.json')):
        bind(root, manifest, run['source_epoch'][key], prefix + filename, filename)
    quartet = json_value(directory / 'bo/quartet.json')
    require(quartet['schema'] == 'hydride-quartet-v1' and
            run['geometry'] == quartet['geometry'] == json_value(directory / 'geometry/relaxation.json') and
            run['specification'] == quartet['specification'] == json_value(root / prefix / 'input.json'),
            'Run/quartet/standalone geometry or specification disagreement: ' + role)
    volume = run['specification']['volume_source']
    require(run['geometry']['volume_source'] == volume, 'Geometry/input volume provenance disagreement')
    if role in {'F', 'F_h', 'F_cell'}:
        name = str(safe_name(volume['path']))
        require('/' not in name and {key: volume[key] for key in ('bytes', 'sha256')} == VOLUME,
                'Wrong original PW-relaxed LiD volume-source binding')
        bind(root, manifest, volume, prefix + name, name)
    else:
        require(volume is None, 'Gas volume-source must be null')
    epoch = json_value(root / prefix / 'source-manifest.json')
    seen = set()
    require(type(epoch['files']) is list and epoch['files'], 'Empty frozen source epoch')
    for entry in epoch['files']:
        name = str(safe_name(entry['path']))
        require(name not in seen, 'Duplicate frozen source path')
        seen.add(name)
        bind(root, manifest, entry, prefix + 'source/' + name, name)
    actual_sources = {name[len(prefix + 'source/'):] for name in manifest['files']
                      if name.startswith(prefix + 'source/')}
    require(actual_sources == seen, 'Incomplete/extra frozen source epoch: ' + role)
    states = quartet['states']
    require(set(states) == {'occupied-host', 'vacant-host', 'q0', 'q1'}, 'Wrong quartet state families')
    expected_states = [states['occupied-host'], states['vacant-host']]
    for charge in ('q0', 'q1'):
        require(set(states[charge]) == {'occupied-probes', 'vacant-probes'}, 'Wrong probe families')
        for family in ('occupied-probes', 'vacant-probes'):
            require(type(states[charge][family]) is list and states[charge][family], 'Missing probe states')
            expected_states.extend(states[charge][family])
    names = set()
    for state in expected_states:
        binding = state['state_file']
        name = str(safe_name(binding['path']))
        require('/' not in name and name.endswith('.json') and name != 'quartet.json' and name not in names,
                'Duplicate/unsafe standalone state name')
        names.add(name)
        bind(root, manifest, binding, calculation + 'bo/' + name, name)
        standalone = json_value(directory / 'bo' / name)
        require(standalone == {key: value for key, value in state.items() if key != 'state_file'},
                'Quartet/standalone state disagreement: ' + role + '/' + name)
    actual_states = {path.name for path in (directory / 'bo').glob('*.json')} - {'quartet.json'}
    require(actual_states == names, 'Missing/unexpected standalone state JSON: ' + role)
    return directory


def replay(root):
    root = Path(root).resolve()
    manifest = verify_inventory(root)
    inputs = json_value(root / 'inputs.json')
    require(set(inputs) == {'runs', 'readout', 'curves', 'extra', 'omissions', 'original_private_configuration'},
            'Wrong portable input shape')
    roles = set(inputs['runs'])
    require(roles in (ROLES, ROLES | {'G2'}), 'Both F controls and same-code gas are required')
    require(inputs['runs'] == {role: 'native/' + role + '/calculation' for role in roles} and
            isinstance(inputs['readout'], str) and Path(inputs['readout']).parent == Path('results') and
            Path(inputs['readout']).suffix == '.json' and inputs['readout'] in manifest['files'] and
            inputs['curves'] == 'curves' and
            inputs['omissions'] == manifest['omissions'], 'Wrong portable layout or omissions')
    require({path.name for path in (root / 'native').iterdir()} == roles, 'Unexpected/missing Native unit')
    byte_binding(inputs['original_private_configuration'])
    require(type(inputs['extra']) is list and all(type(item) is dict and set(item) == {'destination'}
            and item['destination'].startswith('provenance/') for item in inputs['extra']), 'Wrong curated extras')
    extras = {item['destination'] for item in inputs['extra']}
    require(len(extras) == len(inputs['extra']) and
            extras == {name for name in manifest['files'] if name.startswith('provenance/')},
            'Extra provenance inventory disagreement')
    volume_sources = [name for name in extras if manifest['files'][name]['ORIGINAL'] == VOLUME]
    require(volume_sources and all(manifest['files'][name]['EXPORTED'] == VOLUME for name in volume_sources),
            'Original PW-relaxed LiD volume primary missing or transformed')
    runs = {role: verify_native(root, manifest, role) for role in sorted(roles)}
    # Admission must exercise actual all-state scientific guards, not a saved qualified flag.
    from readout import summarize
    from export_pairs import export_curves
    actual = summarize(runs['F'], runs['G'], runs['F_h'], runs['F_cell'], runs.get('G2'))
    require(actual['schema'] == 'hydride-screening-readout-v1' and actual['qualified'] is True,
            'Actual readout is not qualified')
    saved = json_value(root / inputs['readout'])
    require(same_json(actual, saved), 'Replayed readout differs from saved bound readout')
    require(manifest['qualified'] is True, 'Manifest must preserve actual qualification')
    curves = export_curves(root / inputs['readout'], runs['F'], runs['G'],
                           runs['F_h'], runs['F_cell'], runs.get('G2'))
    required_curves = {case + '_' + charge + '_BO_pair.json'
                       for case in ('F', 'G', 'G2') if case in roles for charge in ('q0', 'q1')}
    require(set(curves) == required_curves and
            {path.name for path in (root / 'curves').iterdir()} == required_curves,
            'Missing/unrequested pair curves')
    for name, expected in curves.items():
        require(same_json(json_value(root / 'curves' / name), expected), 'Replayed curve/metadata differs: ' + name)
        require(expected['schema'] == 'coherencev1-target-conditioned-pair-v1' and
                expected['observable'] == 'target_conditioned_BO_total_interaction', 'Wrong pinned pair schema')
    return dict(schema='hydride-screening-capsule-replay-v1', qualified=True, files=len(manifest['files']),
                rows=actual['rows'], F_controls=actual['F_controls'], ratios_R01=actual['ratios_R01'],
                sources=actual['sources'], curves=sorted(curves), omissions=manifest['omissions'], publishing=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('capsule', type=Path)
    args = parser.parse_args()
    try:
        print(json_bytes(replay(args.capsule)).decode(), end='')
    except (ValueError, KeyError, TypeError, OSError, UnicodeError, ImportError) as error:
        print('REFUSED: ' + str(error), file=sys.stderr)
        sys.exit(2)
