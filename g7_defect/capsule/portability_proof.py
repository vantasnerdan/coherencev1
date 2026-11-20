"""Extract a real capsule, replay without Native dependencies, and retain two refusals."""
import argparse
import builtins
import hashlib
import importlib.abc
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile

sys.dont_write_bytecode = True
FORBIDDEN = frozenset({'gpaw', '_gpaw', 'ase', 'scipy', 'matplotlib', 'pylab', 'plotly',
                       'seaborn', 'substrate_framework', 'theorem'})
ARCHIVE_ROOT = 'defect-screening-delivery'
MAX_FILES = 20000
MAX_BYTES = 512 * 1024 * 1024


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()


def safe_member(name):
    require(isinstance(name, str) and name and not any(value in name for value in ('\\', ':', '\x00')),
            'Unsafe archive path')
    path = PurePosixPath(name)
    require(not path.is_absolute() and str(path) == name and
            all(part not in ('', '.', '..') for part in name.split('/')) and
            len(path.parts) > 1 and path.parts[0] == ARCHIVE_ROOT, 'Wrong/unsafe archive root')
    return path


def extract(archive, destination):
    destination.mkdir(parents=True, exist_ok=False)
    seen, total = set(), 0
    # No extractall: links, devices, directory aliases and path traversal are never accepted.
    with tarfile.open(archive, mode='r:gz') as handle:
        for member in handle:
            relative = safe_member(member.name)
            require(member.isfile() and not member.issparse(), 'Nonregular/sparse archive member')
            require(member.name not in seen, 'Duplicate archive member')
            seen.add(member.name)
            total += member.size
            require(len(seen) <= MAX_FILES and 0 <= member.size <= MAX_BYTES and total <= MAX_BYTES,
                    'Archive exceeds finite extraction limits')
            target = destination.joinpath(*relative.parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            source = handle.extractfile(member)
            require(source is not None, 'Missing archive payload')
            with source, target.open('xb') as output:
                shutil.copyfileobj(source, output, length=1024 * 1024)
            require(target.stat().st_size == member.size, 'Truncated archive member')
    require(seen, 'Empty capsule archive')
    return destination / ARCHIVE_ROOT


def worker(root):
    class BlockNative(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            if fullname.split('.')[0] in FORBIDDEN:
                raise ImportError('Portability proof blocks ' + fullname)
            return None
    original_import = builtins.__import__
    def guarded(name, *args, **kwargs):
        if name.split('.')[0] in FORBIDDEN:
            raise ImportError('Portability proof blocks ' + name)
        return original_import(name, *args, **kwargs)
    require(not any(name.split('.')[0] in FORBIDDEN for name in sys.modules),
            'Forbidden dependency already loaded')
    sys.meta_path.insert(0, BlockNative())
    builtins.__import__ = guarded
    for name in sorted(FORBIDDEN):
        try:
            __import__(name)
        except ImportError as error:
            require(str(error).startswith('Portability proof blocks '), 'Unexpected import failure')
        else:
            raise ValueError('Import blocker failed: ' + name)
    root = Path(root).resolve()
    sys.path.insert(0, str(root))
    from reader import replay
    result = replay(root)
    require(not any(name.split('.')[0] in FORBIDDEN for name in sys.modules), 'Forbidden dependency loaded')
    result['blocked_imports'] = sorted(FORBIDDEN)
    result['isolated_python'] = sys.flags.isolated == 1
    require(result['isolated_python'], 'Worker requires isolated Python')
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


def semantic_variant(root):
    """Change only an actual A probe's electron count; rebind all earlier integrity gates."""
    manifest_path = root / 'manifest.json'
    manifest = json.loads(manifest_path.read_bytes())
    changes = []
    def replace(name, value):
        path = root / name
        before = digest(path.read_bytes())
        path.write_bytes(json_bytes(value))
        after = digest(path.read_bytes())
        entry = manifest['files'][name]
        entry.update(ORIGINAL=after, EXPORTED=after, transformation=[],
                     provenance='NEGATIVE TEST ONLY: semantic wrong-input mutation of actual archive bytes')
        changes.append(dict(path=name, before=before, after=after))
        return after
    quartet_name = 'native/A/bo/quartet.json'
    quartet = json.loads((root / quartet_name).read_bytes())
    probe = quartet['states']['occupied-probes'][0]
    old_electrons = probe['electrons']
    probe['electrons'] = old_electrons + 1
    stem = 'occupied-probe-R' + format(probe['R_A'], '.3f').replace('.', 'p')
    state_name = 'native/A/bo/' + stem + '.json'
    state = json.loads((root / state_name).read_bytes())
    require(state['electrons'] == old_electrons, 'Actual probe records disagree before negative mutation')
    state['electrons'] = probe['electrons']
    replace(state_name, state)
    new_quartet = replace(quartet_name, quartet)
    run_name = 'native/A/run-result.json'
    run = json.loads((root / run_name).read_bytes())
    run['quartet'].update(new_quartet)
    new_run = replace(run_name, run)
    readout_name = 'results/readout.json'
    readout = json.loads((root / readout_name).read_bytes())
    readout['sources']['A']['run'].update(new_run)
    readout['sources']['A']['quartet'].update(new_quartet)
    replace(readout_name, readout)
    manifest_path.write_bytes(json_bytes(manifest))
    return dict(case='A', probe=stem, electron_change=1, old_electrons=old_electrons,
                wrong_electrons=probe['electrons'], rebinding=changes,
                purpose='Only a negative variant of the real archive; no fabricated positive Native fixture.')


def main(archive, output):
    archive, output = Path(archive), Path(output)
    require(archive.is_file() and not archive.is_symlink(), 'Actual regular capsule archive required')
    require(not output.exists() and not output.is_symlink(), 'Refusing existing proof output')
    output.mkdir(parents=True, exist_ok=False)
    archive_binding = digest(archive.read_bytes())
    actual = extract(archive, output / 'actual')
    # The worker script itself comes from the actual archive, not an external installed reader.
    command = [sys.executable, '-I', str(actual / 'portability_proof.py'), '--worker']
    results = {}
    def run(label, target, expected):
        completed = subprocess.run(command + [str(target)], capture_output=True, check=False)
        (output / (label + '.stdout')).write_bytes(completed.stdout)
        (output / (label + '.stderr')).write_bytes(completed.stderr)
        (output / (label + '.exit')).write_text(str(completed.returncode) + '\n')
        (output / (label + '.argv.json')).write_bytes(json_bytes(
            ['python', '-I', 'ACTUAL_ARCHIVE/portability_proof.py', '--worker', label + '/defect-screening-delivery']))
        results[label] = dict(returncode=completed.returncode, stdout=digest(completed.stdout),
                              stderr=digest(completed.stderr), expected=expected)
        # Save every disposition before evaluating it, including an unexpected refusal.
        (output / 'attempt-results.json').write_bytes(json_bytes(dict(archive=archive_binding, results=results)))
        require((completed.returncode == 0) == (expected == 'accept'), label + ': unexpected CLI disposition')
        return completed
    replay = run('actual-replay', actual, 'accept')
    report = json.loads(replay.stdout)
    byte_root = output / 'byte-tamper' / ARCHIVE_ROOT
    byte_root.parent.mkdir()
    shutil.copytree(actual, byte_root)
    with (byte_root / 'baseline/B-reused.json').open('ab') as handle:
        handle.write(b' ')
    refusal = run('byte-tamper', byte_root, 'refuse')
    expected = 'EXPORTED byte integrity failed: baseline/B-reused.json'
    require(expected in refusal.stderr.decode('utf-8'), 'Byte mutation did not reach byte integrity refusal')
    results['byte-tamper']['refusal_reason'] = expected
    semantic_root = output / 'extra-probe-electron' / ARCHIVE_ROOT
    semantic_root.parent.mkdir()
    shutil.copytree(actual, semantic_root)
    mutation = semantic_variant(semantic_root)
    (output / 'semantic-mutation.json').write_bytes(json_bytes(mutation))
    refusal = run('extra-probe-electron', semantic_root, 'refuse')
    expected = 'Probe must not add an electron: changed values'
    require(expected in refusal.stderr.decode('utf-8'), 'Semantic mutation did not reach the explicit physical electron guard')
    require('EXPORTED byte integrity failed' not in refusal.stderr.decode('utf-8') and
            'Changed Native quartet bytes' not in refusal.stderr.decode('utf-8'),
            'Semantic negative was refused only by an integrity gate')
    results['extra-probe-electron']['refusal_reason'] = expected
    require(digest(archive.read_bytes()) == archive_binding, 'Original archive changed during proof')
    proof = dict(schema='pdd-defect-screening-portability-proof-v2', passed=True,
                 archive=archive_binding, manifest=digest((actual / 'manifest.json').read_bytes()),
                 blocked_imports=sorted(FORBIDDEN), isolated_python=True, actual_positive_fixture='unmodified extracted archive',
                 results=results, semantic_mutation=mutation, point_qualified=report['point_qualified'],
                 cell_qualified=False, qualification_scope=report['qualification_scope'],
                 C_controls=report['C_controls'], E_required=report['E_required'], E_replayed=report['E_replayed'],
                 rows=report['rows'], curves=report['curves'], publishing=False)
    # Raw originals/negative variants/stdout/stderr remain intact. The same inherited scanner
    # checks every exportable proof text; internal extracted variants are checked by the reader.
    sys.path.insert(0, str(actual))
    from reader import scan_text
    for path in sorted(output.iterdir()):
        if path.is_file():
            scan_text(path.name, path.read_bytes())
    data = json_bytes(proof)
    scan_text('portability-proof.json', data)
    (output / 'portability-proof.json').write_bytes(data)
    print(json.dumps(proof, indent=2, sort_keys=True, allow_nan=False))
    return proof


if __name__ == '__main__':
    try:
        if len(sys.argv) == 3 and sys.argv[1] == '--worker':
            worker(sys.argv[2])
        else:
            parser = argparse.ArgumentParser(description=__doc__)
            parser.add_argument('archive', type=Path)
            parser.add_argument('output', type=Path)
            args = parser.parse_args()
            main(args.archive, args.output)
    except (ValueError, KeyError, TypeError, OSError, tarfile.TarError) as error:
        print('REFUSED: ' + str(error), file=sys.stderr)
        sys.exit(2)
