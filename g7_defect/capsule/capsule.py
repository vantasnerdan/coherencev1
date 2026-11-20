"""Build a private-path-clean, finite saved-output defect capsule; never launch Native."""
import argparse
import gzip
import io
import json
from pathlib import Path
import re
import sys
import tarfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent
LEGACY = SOURCE.parent / 'delivery'
sys.path.insert(0, str(LEGACY))
from build_capsule import Builder, json_bytes, digest
from capsule_io import require
sys.path.insert(0, str(HERE))
from reader import SCHEMA, ROLES, safe_name, json_value, scan_text

SENSITIVE_NAME = re.compile(r'(?i)(?:^|[-_. /])(?:credentials?|secrets?|tokens?|sessions?|provider|agent-report|full-report)(?:$|[-_. /])')
SENSITIVE_FIELD = re.compile(r'(?i)^(?:password|api_key|access_token|refresh_token|bearer_token|private_key|credential|credentials|provider_session|session_id|complete_agent_report|agent_messages)$')


def extra_source(source, destination):
    safe_name(destination)
    require(destination.startswith('provenance/'), 'Extras belong only under provenance/')
    require(not SENSITIVE_NAME.search(destination) and not SENSITIVE_NAME.search(source.name),
            'Forbidden credential/provider/agent-report extra')
    require(source.is_file() and not source.is_symlink(), 'Extra must be a real text file')
    require(source.suffix in {'.json', '.txt', '.log', '.stdout', '.stderr', '.exit', '.csv', '.tsv', '.md'},
            'Forbidden extra format; no GPW/PDF/provider blobs')
    text = source.read_text(encoding='utf-8')
    if source.suffix == '.json':
        def inspect(value):
            if isinstance(value, dict):
                require(not any(SENSITIVE_FIELD.fullmatch(key) for key in value), 'Private provider/credential/agent content in extra')
                for child in value.values():
                    inspect(child)
            elif isinstance(value, list):
                for child in value:
                    inspect(child)
        inspect(json_value(source))
    require(not re.search(r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----|\bBearer\s+[A-Za-z0-9._-]+', text),
            'Private credential content in extra')


def build(config_path, root, archive, receipt, readme):
    config = json_value(config_path)
    require(set(config) == {'runs', 'readout', 'baseline', 'curves', 'omissions', 'extra', 'cell_sensitivity'},
            'Wrong v2 capsule config keys')
    omissions = config['omissions']
    require(isinstance(omissions, list) and len(omissions) == 3 and
            all(isinstance(item, str) and item.strip() for item in omissions), 'Exactly three Root-written omissions required')
    readout = json_value(Path(config['readout']))
    require(readout['schema'] == 'pd-defect-screening-readout-v2' and readout['point_qualified'] is True and
            readout['cell_qualified'] is False and type(readout['E_required']) is bool,
            'Point-qualified v2 readout with explicitly omitted cell convergence required')
    roles = set(ROLES) | ({'E'} if readout['E_required'] else set())
    require(set(config['runs']) == roles, 'Required completed roles missing, or untriggered E supplied')
    runs = {role: Path(config['runs'][role]) for role in roles}
    baseline = Path(config['baseline'])
    curves = Path(config['curves'])
    require(config['cell_sensitivity'] is None or
            (isinstance(config['cell_sensitivity'], str) and bool(config['cell_sensitivity'].strip())),
            'cell_sensitivity must be null or an actual directory string')
    sensitivity = None if config['cell_sensitivity'] is None else Path(config['cell_sensitivity'])
    if sensitivity is not None:
        require(sensitivity.is_dir() and not sensitivity.is_symlink(), 'Actual sensitivity directory required')
        epoch = json_value(sensitivity / 'source-manifest.json')
        for item in epoch['files']:
            safe_name(item['path'])
    # Admission is actual saved-output replay, not a filename/completed-flag assertion.
    sys.path.insert(0, str(SOURCE))
    from readout_defect import summarize, load_run, RADII
    from export_defect_pair import export_curves
    actual = summarize(runs['A'], runs['C'], runs['D'], runs['h016'], baseline, sensitivity)
    require(actual == readout and actual['point_qualified'] is True and actual['cell_qualified'] is False,
            'Source-bound v2 readout does not replay or misstates qualification')
    if 'E' in runs:
        load_run(runs['E'], 'E', RADII)
    expected_curves = export_curves(Path(config['readout']), baseline, runs['A'], runs['C'], runs['D'], runs.get('E'))
    require(curves.is_dir() and {path.name for path in curves.iterdir()} == set(expected_curves),
            'All and only required pair curves must be supplied')
    for name, expected in expected_curves.items():
        safe_name(name)
        require(json_value(curves / name) == expected, 'Unbound curve or metadata: ' + name)
    require((curves / 'B_BO_pair.json').read_bytes() == baseline.read_bytes(), 'B must be copied byte-for-byte')
    require(isinstance(config['extra'], list), 'Extras must be an explicit source/destination list')
    for item in config['extra']:
        require(set(item) == {'source', 'destination'}, 'Wrong extra metadata keys')
        extra_source(Path(item['source']), item['destination'])
    root, archive, receipt, readme = map(Path, (root, archive, receipt, readme))
    require(not any(path.exists() or path.is_symlink() for path in (root, archive, receipt)), 'Refusing existing delivery output')
    resolved = [path.resolve() for path in (root, archive, receipt)]
    require(len(set(resolved)) == 3 and not any(path.is_relative_to(resolved[0]) for path in resolved[1:]),
            'Archive and receipt must be distinct and outside capsule')
    require(readme.is_file() and not readme.is_symlink(), 'Root-written README required')
    root.mkdir(parents=True, exist_ok=False)
    builder = Builder(root)
    # This wrapper enforces names before the inherited Builder writes anything.
    def add(source, destination, exact=True):
        safe_name(destination)
        require(source.is_file() and not source.is_symlink(), 'Missing/unsafe actual source: ' + destination)
        builder.add(source, destination, exact=exact)
    for role in sorted(roles):
        directory = runs[role]
        prefix = 'native/' + role + '/'
        run = json_value(directory / 'run-result.json')
        add(directory / 'run-result.json', prefix + 'run-result.json')
        for field, filename in (('input', 'input.json'), ('frozen_source', 'source-manifest.json')):
            origin = directory.parent / filename
            binding = run[field]
            require(binding['path'] == filename and digest(origin.read_bytes()) ==
                    {key: binding[key] for key in ('bytes', 'sha256')}, 'Changed original Native ' + filename)
            add(origin, prefix + filename)
        epoch = json_value(directory.parent / 'source-manifest.json')
        for item in epoch['files']:
            name = str(safe_name(item['path']))
            origin = directory.parent / 'source' / name
            require(digest(origin.read_bytes()) == {key: item[key] for key in ('bytes', 'sha256')},
                    'Changed frozen Native executable epoch: ' + name)
            add(origin, prefix + 'source/' + name)
        require((directory / 'bo/quartet.json').is_file(), 'Missing completed Native quartet')
        for origin in sorted((directory / 'bo').iterdir()):
            if origin.suffix in ('.json', '.txt'):
                add(origin, prefix + 'bo/' + origin.name, exact=origin.suffix == '.json')
        # Preserve original optimizer outcomes and selected-frame provenance, without orbitals/trajectories.
        for filename in ('geometry-selection.json', 'relax/relaxation.json'):
            origin = directory / filename
            if origin.exists():
                add(origin, prefix + filename)
        relaxation = directory / 'relax'
        if relaxation.is_dir():
            for origin in sorted(relaxation.iterdir()):
                if origin.suffix in ('.json', '.txt', '.log') and origin.name != 'relaxation.json':
                    add(origin, prefix + 'relax/' + origin.name, exact=origin.suffix == '.json')
        if (directory.parent / 'native.log').is_file():
            add(directory.parent / 'native.log', prefix + 'native.log', exact=False)
    if sensitivity is not None:
        for filename in ('input.json', 'source-manifest.json', 'cell-sensitivity.json'):
            add(sensitivity / filename, 'sensitivity/' + filename)
        epoch = json_value(sensitivity / 'source-manifest.json')
        for item in epoch['files']:
            name = str(safe_name(item['path']))
            origin = sensitivity / 'source' / name
            require(digest(origin.read_bytes()) == {key: item[key] for key in ('bytes', 'sha256')},
                    'Changed original sensitivity executable epoch: ' + name)
            add(origin, 'sensitivity/source/' + name)
        for role in ('occupied-host', 'vacant-host', 'occupied-probe', 'vacant-probe'):
            for filename in ('state.json', 'cell-role.json'):
                add(sensitivity / 'roles' / role / filename, 'sensitivity/roles/' + role + '/' + filename)
            for filename in ('state.txt', 'native.log'):
                origin = sensitivity / 'roles' / role / filename
                if origin.exists():
                    add(origin, 'sensitivity/roles/' + role + '/' + filename, exact=False)
    add(Path(config['readout']), 'results/readout.json')
    add(baseline, 'baseline/B-reused.json')
    for name in sorted(expected_curves):
        add(curves / name, 'curves/' + name)
    for name in ('readout_defect.py', 'export_defect_pair.py', 'reuse_baseline.py', 'cell_sensitivity.py', 'compare_relaxation.py'):
        add(SOURCE / name, name)
    for name in ('green.py', 'capsule_io.py', 'build_capsule.py'):
        add(LEGACY / name, name)
    for name in ('reader.py', 'portability_proof.py', 'capsule.py'):
        add(HERE / name, name)
    add(readme, 'README.md')
    for item in config['extra']:
        add(Path(item['source']), item['destination'], exact=False)
    public_inputs = dict(runs={role: 'native/' + role for role in sorted(roles)},
                         readout='results/readout.json', baseline='baseline/B-reused.json',
                         curves='curves', omissions=omissions,
                         cell_sensitivity=None if sensitivity is None else 'sensitivity',
                         extra=[dict(destination=item['destination']) for item in config['extra']],
                         original_private_configuration=digest(Path(config_path).read_bytes()))
    builder.generated('inputs.json', json_bytes(public_inputs))
    manifest = dict(schema=SCHEMA, files=builder.files, omissions=omissions,
                    point_qualified=actual['point_qualified'], cell_qualified=False,
                    qualification_scope=actual['qualification_scope'],
                    provenance='ORIGINAL binds immutable source bytes; only declared private paths and GPAW User headers may differ in EXPORTED text. Scientific Native JSON and reused B are exact. Manifest is externally bound by archive SHA256, not a signature.',
                    excluded_payloads=['GPW wavefunctions', 'PDFs', 'credentials/provider-session blobs', 'complete agent reports'],
                    reader_command='python reader.py CAPSULE_DIRECTORY',
                    portability_command='python portability_proof.py ARCHIVE OUTPUT_DIRECTORY',
                    native_invocations=0, publishing=False)
    (root / 'manifest.json').write_bytes(json_bytes(manifest))
    for path in sorted(root.rglob('*')):
        if path.is_file():
            scan_text(path.relative_to(root).as_posix(), path.read_bytes())
    # Seal only after the complete staged capsule passes its real reader.
    sys.path.insert(0, str(HERE))
    from reader import replay
    report = replay(root)
    archive.parent.mkdir(parents=True, exist_ok=True)
    with archive.open('xb') as output:
        with gzip.GzipFile(filename='', mode='wb', fileobj=output, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode='w', format=tarfile.PAX_FORMAT) as tar:
                for path in sorted(root.rglob('*')):
                    if path.is_file():
                        data = path.read_bytes()
                        member = tarfile.TarInfo('defect-screening-delivery/' + path.relative_to(root).as_posix())
                        member.size = len(data)
                        member.mtime = 0
                        member.mode = 0o644
                        tar.addfile(member, io.BytesIO(data))
    result = dict(schema='pdd-defect-screening-asset-v2', archive=digest(archive.read_bytes()),
                  manifest=digest((root / 'manifest.json').read_bytes()), files=len(builder.files),
                  point_qualified=report['point_qualified'], cell_qualified=False,
                  qualification_scope=report['qualification_scope'], E_required=report['E_required'],
                  portability_proof='Not run by builder; execute the actual archive portability command.',
                  privacy_scan=dict(passed=True, findings=[], scope='Every exported UTF-8 payload, including sources and any supplied proofs; private paths, GPAW User headers and emails refused.'),
                  publishing=False)
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with receipt.open('xb') as output:
        output.write(json_bytes(result))
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('config', 'root', 'archive', 'receipt', 'readme'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args.config, args.root, args.archive, args.receipt, args.readme)
    except (ValueError, KeyError, TypeError, OSError) as error:
        print('REFUSED: ' + str(error), file=sys.stderr)
        sys.exit(2)
