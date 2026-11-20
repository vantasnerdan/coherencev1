"""Verify a finite defect capsule and replay saved Native outputs with NumPy only."""
import argparse
import math
import json
from pathlib import Path, PurePosixPath
import re
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from capsule_io import digest, require
from build_capsule import PRIVATE_PATH, GPAW_USER_HEADER

SCHEMA = 'pdd-defect-screening-capsule-v2'
ROLES = ('A', 'C', 'D', 'h016')
EMAIL = re.compile(r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}')
TEXT_SUFFIXES = {'.json', '.py', '.md', '.txt', '.log', '.stdout', '.stderr', '.exit', '.csv', '.tsv'}
B_ORIGINAL_SHA256 = '854c570dc53ac5f0ff792f38393dea18c66dab7aa474d4f03b29d5e4e98a32a5'


def safe_name(name):
    require(isinstance(name, str) and bool(name), 'Empty/nontext capsule path')
    path = PurePosixPath(name)
    require(not path.is_absolute() and str(path) == name and
            all(part not in ('', '.', '..') for part in name.split('/')) and
            not any(character in name for character in ('\\', ':', '\x00')),
            'Unsafe capsule path: ' + repr(name))
    return path


def json_value(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'Duplicate JSON key: ' + key)
            result[key] = value
        return result
    def invalid(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    def floating(value):
        number = float(value)
        require(math.isfinite(number), 'Nonfinite JSON number')
        return number
    return json.loads(Path(path).read_bytes(), object_pairs_hook=pairs, parse_constant=invalid, parse_float=floating)


def scan_text(name, data):
    suffix = PurePosixPath(name).suffix.lower()
    forbidden = {'.gpw', '.pdf', '.npz', '.npy', '.pkl', '.pickle', '.tar', '.gz', '.zip',
                 '.png', '.jpg', '.sqlite', '.db', '.env', '.pem', '.key'}
    require(suffix not in forbidden and (suffix in TEXT_SUFFIXES or name.startswith('provenance/')),
            'Forbidden/nontext capsule payload: ' + name)
    text = data.decode('utf-8')
    require(not (PRIVATE_PATH.search(text) or GPAW_USER_HEADER.search(text) or EMAIL.search(text)),
            'Export privacy scan failed: ' + name)


def verify_inventory(root):
    root = Path(root)
    require(root.is_dir() and not root.is_symlink(), 'Capsule must be a real directory')
    manifest = json_value(root / 'manifest.json')
    require(manifest['schema'] == SCHEMA, 'Wrong defect capsule schema')
    omissions = manifest['omissions']
    require(isinstance(omissions, list) and len(omissions) == 3 and
            all(isinstance(item, str) and item.strip() for item in omissions), 'Exactly three nonempty omissions required')
    files = manifest['files']
    require(isinstance(files, dict) and files and 'manifest.json' not in files, 'Invalid file inventory')
    actual = set()
    for path in root.rglob('*'):
        require(not path.is_symlink(), 'Symlink in capsule')
        require(path.is_dir() or path.is_file(), 'Special file in capsule')
        if path.is_file():
            name = path.relative_to(root).as_posix()
            safe_name(name)
            actual.add(name)
    require(actual == set(files) | {'manifest.json'}, 'Unmanifested or missing capsule files')
    for name, entry in files.items():
        safe_name(name)
        data = (root / name).read_bytes()
        require(digest(data) == entry['EXPORTED'], 'EXPORTED byte integrity failed: ' + name)
        require(isinstance(entry['transformation'], list), 'Missing declared transformation: ' + name)
        for field in ('ORIGINAL', 'EXPORTED'):
            bound = entry[field]
            require(set(bound) == {'bytes', 'sha256'} and type(bound['bytes']) is int and bound['bytes'] >= 0 and
                    isinstance(bound['sha256'], str) and re.fullmatch('[0-9a-f]{64}', bound['sha256']),
                    'Malformed byte binding: ' + name)
        if not entry['transformation']:
            require(entry['ORIGINAL'] == entry['EXPORTED'], 'Undeclared privacy transformation: ' + name)
        if ((name.startswith(('native/', 'sensitivity/', 'results/', 'baseline/', 'curves/')) and name.endswith('.json'))
                or name.endswith('.py') or name == 'README.md'):
            require(entry['ORIGINAL'] == entry['EXPORTED'] and not entry['transformation'],
                    'Scientific/source bytes must be exact: ' + name)
        scan_text(name, data)
    scan_text('manifest.json', (root / 'manifest.json').read_bytes())
    return manifest

def bind(root, manifest, binding, name):
    safe_name(name)
    require(binding['path'] == PurePosixPath(name).name, 'Wrong bound source basename: ' + name)
    require({key: binding[key] for key in ('bytes', 'sha256')} == manifest['files'][name]['ORIGINAL'],
            'ORIGINAL source binding mismatch: ' + name)
    require(digest((root / name).read_bytes()) == manifest['files'][name]['ORIGINAL'],
            'Scientific binding was transformed: ' + name)


def verify_native(root, manifest, role):
    directory = root / 'native' / role
    run = json_value(directory / 'run-result.json')
    quartet = json_value(directory / 'bo/quartet.json')
    bind(root, manifest, run['quartet'], 'native/' + role + '/bo/quartet.json')
    bind(root, manifest, run['input'], 'native/' + role + '/input.json')
    bind(root, manifest, run['frozen_source'], 'native/' + role + '/source-manifest.json')
    require(json_value(directory / 'input.json') == quartet['specification'], 'Changed actual input specification: ' + role)
    epoch = json_value(directory / 'source-manifest.json')
    seen = set()
    for entry in epoch['files']:
        name = str(safe_name(entry['path']))
        require(name not in seen, 'Duplicate frozen source entry')
        seen.add(name)
        # Frozen source bindings carry source-relative paths, not only basenames.
        destination = 'native/' + role + '/source/' + name
        require({key: entry[key] for key in ('bytes', 'sha256')} == manifest['files'][destination]['ORIGINAL'],
                'ORIGINAL frozen executable source mismatch: ' + destination)
    require(seen, 'Empty frozen executable epoch')
    states = quartet['states']
    expected = {}
    for family in ('occupied', 'vacant'):
        expected[family + '-host.json'] = states[family + '-host']
        for probe in states[family + '-probes']:
            stem = family + '-probe-R' + format(probe['R_A'], '.3f').replace('.', 'p')
            require(stem + '.json' not in expected, 'Duplicate Native probe radius')
            expected[stem + '.json'] = {key: value for key, value in probe.items() if key not in ('R_A', 'projectile_A')}
    actual = {path.name for path in (directory / 'bo').glob('*.json')} - {'quartet.json'}
    require(actual == set(expected), 'Missing/unexpected per-state Native JSON: ' + role)
    for name, state in expected.items():
        require(json_value(directory / 'bo' / name) == state, 'Quartet/per-state Native disagreement: ' + role + '/' + name)
    require(run['rows'] == quartet['rows'], 'Native run/quartet table disagreement: ' + role)
    geometry = 'geometry-selection.json' if 'source_frame' in run['relaxation'] else 'relax/relaxation.json'
    require(json_value(directory / geometry) == run['relaxation'], 'Native standalone geometry disagreement: ' + role)
    return directory


def verify_sensitivity(root, manifest):
    prefix = 'sensitivity/'
    required = {prefix + name for name in ('input.json', 'source-manifest.json', 'cell-sensitivity.json')}
    epoch = json_value(root / prefix / 'source-manifest.json')
    seen = set()
    for item in epoch['files']:
        name = str(safe_name(item['path']))
        require(name not in seen, 'Duplicate sensitivity frozen source')
        seen.add(name)
        destination = prefix + 'source/' + name
        required.add(destination)
        require({key: item[key] for key in ('bytes', 'sha256')} == manifest['files'][destination]['ORIGINAL'],
                'ORIGINAL sensitivity executable binding mismatch: ' + destination)
    require(seen, 'Empty sensitivity frozen source epoch')
    optional = set()
    for role in ('occupied-host', 'vacant-host', 'occupied-probe', 'vacant-probe'):
        folder = prefix + 'roles/' + role + '/'
        required.update(folder + name for name in ('state.json', 'cell-role.json'))
        optional.update(folder + name for name in ('state.txt', 'native.log'))
    actual = {name for name in manifest['files'] if name.startswith(prefix)}
    require(required <= actual <= required | optional, 'Missing/unexpected sensitivity payload')
    return root / 'sensitivity'


def replay(root):
    root = Path(root)
    manifest = verify_inventory(root)
    inputs = json_value(root / 'inputs.json')
    expected_roles = set(ROLES)
    saved = json_value(root / 'results/readout.json')
    require(saved['schema'] == 'pd-defect-screening-readout-v2' and saved['point_qualified'] is True,
            'Point-qualified v2 saved readout required')
    require(saved['cell_qualified'] is False, 'Converged cell-size qualification remains omitted')
    require(type(saved['E_required']) is bool, 'Missing point-qualified E decision')
    if saved['E_required']:
        expected_roles.add('E')
    require(set(inputs['runs']) == expected_roles, 'Missing/unrequested actual Native roles')
    require({path.name for path in (root / 'native').iterdir()} == expected_roles,
            'Missing/unrequested Native role directory')
    require(inputs['runs'] == {role: 'native/' + role for role in expected_roles} and
            inputs['readout'] == 'results/readout.json' and inputs['baseline'] == 'baseline/B-reused.json' and
            inputs['curves'] == 'curves', 'Inputs must name the fixed relative capsule layout')
    require(inputs['omissions'] == manifest['omissions'], 'Omission declaration disagreement')
    require('cell_sensitivity' in inputs and
            (inputs['cell_sensitivity'] is None or inputs['cell_sensitivity'] == 'sensitivity'),
            'Sensitivity must be null or the fixed relative sensitivity directory')
    sensitivity = None
    if inputs['cell_sensitivity'] is not None:
        sensitivity = verify_sensitivity(root, manifest)
    else:
        require(not any(name.startswith('sensitivity/') for name in manifest['files']),
                'Undeclared optional sensitivity payload')
    # Replay first: a semantically wrong but fully rebound probe must reach the physical guard.
    from readout_defect import summarize, load_run, RADII
    actual = summarize(root / 'native/A', root / 'native/C', root / 'native/D',
                       root / 'native/h016', root / 'baseline/B-reused.json', sensitivity)
    require(actual == saved, 'Replayed readout differs from saved bound readout')
    require(manifest['point_qualified'] is True and manifest['cell_qualified'] is False and
            manifest['qualification_scope'] == actual['qualification_scope'],
            'Manifest must preserve point-only qualification and omitted converged cell-size check')
    require(actual['point_qualified'] is True and actual['C_controls']['h']['passed'] is True,
            'Unqualified actual C grid control')
    cell = actual['C_controls']['cell']
    require(actual['cell_qualified'] is False and cell['converged'] is False and 'passed' not in cell and
            cell['status'] == ('omitted' if sensitivity is None else 'unrelaxed_sensitivity'),
            'Waived converged cell-size check must remain an omission or unrelaxed sensitivity')
    runs = {role: verify_native(root, manifest, role) for role in expected_roles}
    if saved['E_required']:
        load_run(runs['E'], 'E', RADII)
    baseline = json_value(root / 'baseline/B-reused.json')
    reuse = baseline['g7_reuse']
    require(reuse['original_sha256'] == B_ORIGINAL_SHA256 and reuse['original_bytes'] == 3730 and
            type(reuse['native_runs']) is int and reuse['native_runs'] == 0 and
            reuse['selected_original_indices'] == [0, 1, 2, 4, 5, 6] and reuse['excluded_original_radii_A'] == [.4],
            'B is not the exact registered six-point reuse')
    from export_defect_pair import export_curves
    curves = export_curves(root / 'results/readout.json', root / 'baseline/B-reused.json',
                           runs['A'], runs['C'], runs['D'], runs.get('E'))
    require({path.name for path in (root / 'curves').iterdir()} == set(curves), 'Missing/unrequested pair curves')
    for name, expected in curves.items():
        safe_name(name)
        require(json_value(root / 'curves' / name) == expected, 'Bound numeric curve/metadata mismatch: ' + name)
    require((root / 'curves/B_BO_pair.json').read_bytes() == (root / 'baseline/B-reused.json').read_bytes(),
            'B curve must preserve exact reused bytes')
    return dict(schema='pd-defect-screening-capsule-replay-v2', point_qualified=True, cell_qualified=False,
                files=len(manifest['files']), rows=actual['rows'], contrasts_R01=actual['contrasts_R01'],
                C_controls=actual['C_controls'], qualification_scope=actual['qualification_scope'],
                E_required=actual['E_required'], E_replayed='E' in runs,
                curves=sorted(curves), omissions=manifest['omissions'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('capsule', type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(replay(args.capsule), indent=2, sort_keys=True, allow_nan=False))
    except (ValueError, KeyError, TypeError, OSError) as error:
        print('REFUSED: ' + str(error), file=sys.stderr)
        sys.exit(2)
