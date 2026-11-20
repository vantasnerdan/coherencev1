"""Build the private, source-bound hydride draft; never run Native or publish."""
import argparse
import gzip
import io
from pathlib import Path
import sys
import tarfile

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent
SHARED = SOURCE.parent / 'delivery'
sys.path.insert(0, str(SOURCE))
sys.path.insert(0, str(HERE))
from delivery_io import (ARCHIVE_ROOT, SCHEMA, ROLES, VOLUME, Builder, TEXT_SUFFIXES,
                         SENSITIVE_NAME, digest, json_bytes, json_value,
                         require, safe_name, scan_text)
from reader import replay
from log_privacy import sanitize_log


def extra_source(item):
    require(type(item) is dict and set(item) == {'source', 'destination'}, 'Wrong extra metadata keys')
    source, destination = Path(item['source']), item['destination']
    safe_name(destination)
    require(destination.startswith('provenance/') and source.is_file() and not source.is_symlink(),
            'Curated extras must be real text files under provenance/')
    require(source.suffix.lower() in TEXT_SUFFIXES and Path(destination).suffix.lower() == source.suffix.lower(),
            'Extra format must be explicit UTF-8 text, never GPW/PDF or opaque payloads')
    require(not SENSITIVE_NAME.search(source.name), 'Forbidden private extra source name')
    text = source.read_bytes().decode('utf-8')
    if source.suffix.lower() in {'.txt', '.log', '.stdout', '.stderr', '.exit'}:
        text = sanitize_log(text.encode('utf-8'))[0].decode('utf-8')
    scan_text(destination, text.encode('utf-8'))
    if source.suffix.lower() == '.json':
        json_value(source)


def build(config_path, root, archive, receipt, readme):
    config_path = Path(config_path)
    config = json_value(config_path)
    require(set(config) == {'runs', 'readout', 'curves', 'extra', 'omissions'}, 'Wrong hydride capsule config keys')
    require(type(config['runs']) is dict and set(config['runs']) in (ROLES, ROLES | {'G2'}),
            'F/G and both F controls are mandatory; only G2 is optional')
    omissions = config['omissions']
    require(type(omissions) is list and omissions and all(isinstance(item, str) and item.strip() for item in omissions),
            'Nonempty Root-written omissions required')
    require(type(config['extra']) is list, 'Explicit curated extra list required')
    destinations = set()
    for item in config['extra']:
        extra_source(item)
        require(item['destination'] not in destinations, 'Duplicate extra destination')
        destinations.add(item['destination'])
    require(any(digest(Path(item['source']).read_bytes()) == VOLUME for item in config['extra']),
            'Exact original PW-relaxed LiD volume primary required in curated provenance')
    runs = {role: Path(path) for role, path in config['runs'].items()}
    for directory in runs.values():
        require(directory.is_dir() and not directory.is_symlink() and not directory.parent.is_symlink(),
                'Actual regular Native calculation directories required')
    root, archive, receipt, readme = map(Path, (root, archive, receipt, readme))
    outputs = (root, archive, receipt)
    require(not any(path.exists() or path.is_symlink() for path in outputs), 'Refusing existing delivery output')
    resolved = [path.resolve() for path in outputs]
    require(len(set(resolved)) == 3 and not any(path.is_relative_to(resolved[0]) for path in resolved[1:]),
            'Archive and receipt must be distinct and outside capsule')
    sources = [directory.parent.resolve() for directory in runs.values()] + [Path(config['curves']).resolve()]
    require(not any(resolved[0].is_relative_to(source) or source.is_relative_to(resolved[0]) for source in sources),
            'Staging root must not overlap scientific inputs')
    require(readme.is_file() and not readme.is_symlink(), 'Root-written final README required')
    root.mkdir(parents=True, exist_ok=False)
    builder = Builder(root)
    excluded = []
    for role in sorted(runs):
        directory = runs[role]
        prefix = 'native/' + role + '/'
        unit = directory.parent
        for filename in ('input.json', 'source-manifest.json'):
            builder.add(unit / filename, prefix + filename)
        if role in {'F', 'F_h', 'F_cell'}:
            volume = json_value(unit / 'input.json')['volume_source']
            name = str(safe_name(volume['path']))
            require('/' not in name and {key: volume[key] for key in ('bytes', 'sha256')} == VOLUME,
                    'Wrong original LiD volume-source binding')
            builder.add(unit / name, prefix + name)
        epoch = json_value(unit / 'source-manifest.json')
        require(type(epoch['files']) is list and epoch['files'], 'Empty frozen Native source epoch')
        seen = set()
        for binding in epoch['files']:
            name = str(safe_name(binding['path']))
            require(name not in seen, 'Duplicate frozen Native executable source')
            seen.add(name)
            source = unit / 'source' / name
            require(digest(source.read_bytes()) == {key: binding[key] for key in ('bytes', 'sha256')},
                    'Changed frozen Native executable source: ' + name)
            builder.add(source, prefix + 'source/' + name)
        for path in sorted(directory.rglob('*')):
            require(not path.is_symlink() and (path.is_file() or path.is_dir()), 'Unsafe Native saved-output member')
            if path.is_dir():
                continue
            name = path.relative_to(directory).as_posix()
            if path.suffix.lower() in {'.gpw', '.pdf', '.traj'}:
                excluded.append(prefix + 'calculation/' + name)
                continue
            require(path.suffix.lower() in TEXT_SUFFIXES, 'Unknown saved-output format: ' + name)
            builder.add(path, prefix + 'calculation/' + name,
                        header_only=path.suffix.lower() in {'.txt', '.log', '.stdout', '.stderr', '.exit'})
        builder.add(unit / 'native.log', prefix + 'native.log', header_only=True)
    readout_name = 'results/' + Path(config['readout']).name
    builder.add(Path(config['readout']), readout_name)
    curves = Path(config['curves'])
    expected = {case + '_' + charge + '_BO_pair.json' for case in ('F', 'G', 'G2')
                if case in runs for charge in ('q0', 'q1')}
    require(curves.is_dir() and not curves.is_symlink() and {path.name for path in curves.iterdir()} == expected,
            'All and only required pair curves must be supplied')
    for name in sorted(expected):
        builder.add(curves / name, 'curves/' + name)
    for name in ('readout.py', 'export_pairs.py', 'hydride_replay.py', 'log_privacy.py'):
        builder.add(SOURCE / name, name)
    builder.add(SHARED / 'green.py', 'green.py')
    for name in ('capsule.py', 'reader.py', 'portability_proof.py', 'delivery_io.py'):
        builder.add(HERE / name, name)
    builder.add(readme, 'README.md')
    for item in config['extra']:
        builder.add(Path(item['source']), item['destination'],
                    header_only=Path(item['source']).suffix.lower() in {'.txt', '.log', '.stdout', '.stderr', '.exit'})
    portable = dict(runs={role: 'native/' + role + '/calculation' for role in sorted(runs)},
                    readout=readout_name, curves='curves', omissions=omissions,
                    extra=[dict(destination=item['destination']) for item in config['extra']],
                    original_private_configuration=digest(config_path.read_bytes()))
    builder.generated('inputs.json', json_bytes(portable))
    manifest = dict(schema=SCHEMA, files=builder.files, qualified=True, omissions=omissions,
                    excluded_payloads=excluded, native_invocations=0, publishing=False,
                    provenance='ORIGINAL binds immutable inputs. Only declared GPAW User/CWD/library-path metadata and dataset/checkpoint path prefixes may differ in EXPORTED. Scientific JSON and source are exact. Administrative watch output is excluded. Archive SHA256 externally binds this inventory; it is not a signature.',
                    reader_command='python reader.py CAPSULE_DIRECTORY',
                    portability_command='python portability_proof.py ARCHIVE OUTPUT_DIRECTORY')
    data = json_bytes(manifest)
    scan_text('manifest.json', data)
    with (root / 'manifest.json').open('xb') as output:
        output.write(data)
    # These are the very bundled modules later used by the isolated archive proof.
    sys.path.insert(0, str(root))
    report = replay(root)
    archive.parent.mkdir(parents=True, exist_ok=True)
    with archive.open('xb') as output:
        with gzip.GzipFile(filename='', mode='wb', fileobj=output, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode='w', format=tarfile.PAX_FORMAT) as handle:
                for path in sorted(root.rglob('*')):
                    if path.is_file():
                        data = path.read_bytes()
                        member = tarfile.TarInfo(ARCHIVE_ROOT + '/' + path.relative_to(root).as_posix())
                        member.size = len(data)
                        member.mtime = 0
                        member.mode = 0o644
                        handle.addfile(member, io.BytesIO(data))
    transformed = [dict(path=name, ORIGINAL=entry['ORIGINAL'], EXPORTED=entry['EXPORTED'],
                        transformation=entry['transformation']) for name, entry in builder.files.items()
                   if entry['transformation']]
    result = dict(schema='hydride-screening-asset-v1', archive=digest(archive.read_bytes()),
                  manifest=digest((root / 'manifest.json').read_bytes()), files=len(builder.files),
                  archive_members=len(builder.files) + 1, qualified=report['qualified'],
                  replay=report, findings=[dict(check='actual saved-output replay', disposition='accepted'),
                                           dict(check='complete inventory and UTF-8 privacy', disposition='accepted'),
                                           dict(check='actual archive isolated portability proof', disposition='not-executed-by-builder')],
                  privacy_scan=dict(passed=True, findings=[], transformations=transformed),
                  portability_proof='Separate proof must bind this exact archive SHA256.', publishing=False)
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with receipt.open('xb') as output:
        output.write(json_bytes(result))
    print(json_bytes(result).decode(), end='')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('config', 'root', 'archive', 'receipt', 'readme'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    try:
        build(args.config, args.root, args.archive, args.receipt, args.readme)
    except (ValueError, KeyError, TypeError, OSError, UnicodeError, ImportError, tarfile.TarError) as error:
        print('REFUSED: ' + str(error), file=sys.stderr)
        sys.exit(2)
