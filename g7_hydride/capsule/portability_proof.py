"""Prove the real archive replays in isolated Python and retain both meaningful refusals."""
import argparse
import builtins
import importlib.abc
import importlib.util
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(1, str(Path(__file__).resolve().parent.parent))
from delivery_io import ARCHIVE_ROOT, digest, json_bytes, json_value, require, scan_text

FORBIDDEN = frozenset({'gpaw', '_gpaw', 'gpaw_native', 'ase', 'scipy', 'matplotlib', 'pylab',
                       'plotly', 'seaborn', 'substrate_framework', 'theorem'})
BUNDLED = frozenset({'reader', 'readout', 'export_pairs', 'hydride_replay', 'green', 'delivery_io', 'portability_proof', 'log_privacy'})
MAX_FILES = 20000
MAX_BYTES = 512 * 1024 * 1024


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
    # No extractall: duplicate paths, links, devices, sparse files and aliases are refused.
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
    require(seen, 'Empty archive')
    return destination / ARCHIVE_ROOT


def worker(root, numpy_path):
    allowed = set(sys.stdlib_module_names) | set(sys.builtin_module_names) | BUNDLED | {'numpy', '__main__'}
    require(Path(numpy_path).is_dir(), 'Installed NumPy search directory required')
    sys.path.insert(0, str(Path(numpy_path).resolve()))
    def check(name):
        head = name.split('.')[0]
        if head in FORBIDDEN or head not in allowed:
            raise ImportError('Portability proof blocks ' + name)
    class NumPyOnly(importlib.abc.MetaPathFinder):
        def find_spec(self, fullname, path=None, target=None):
            check(fullname)
            return None
    original_import = builtins.__import__
    def guarded(name, globals=None, locals=None, fromlist=(), level=0):
        if level == 0 and name:
            check(name)
        return original_import(name, globals, locals, fromlist, level)
    require(sys.flags.isolated == 1 and sys.flags.no_site == 1, 'Worker requires isolated Python without site startup')
    require(not any(name.split('.')[0] in FORBIDDEN for name in sys.modules), 'Forbidden dependency already loaded')
    sys.meta_path.insert(0, NumPyOnly())
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
    # Reject an already-loaded nonstdlib third party as well as attempted imports.
    for name in sys.modules:
        require(name.split('.')[0] in allowed, 'Non-NumPy dependency loaded')
    result.update(blocked_imports=sorted(FORBIDDEN), isolated_python=True, dependency_scope='NumPy and stdlib only')
    print(json_bytes(result).decode(), end='')


def semantic_variant(root):
    """Mutate one actual F q1 probe, then rebind the complete JSON dependency graph."""
    manifest_path = root / 'manifest.json'
    manifest = json_value(manifest_path)
    changes, replacements = [], {}
    def replace(name, value):
        path = root / name
        before = digest(path.read_bytes())
        data = json_bytes(value)
        if path.read_bytes() == data:
            return
        path.write_bytes(data)
        after = digest(data)
        replacements[before['sha256']] = (before, after)
        manifest['files'][name].update(ORIGINAL=after, EXPORTED=after, transformation=[])
        changes.append(dict(path=name, before=before, after=after))
    quartet_name = 'native/F/calculation/bo/quartet.json'
    quartet = json_value(root / quartet_name)
    probe = quartet['states']['q1']['occupied-probes'][0]
    old_electrons = probe['electrons']
    require(type(old_electrons) in (int, float), 'Actual probe has no numeric electron count')
    state_name = 'native/F/calculation/bo/' + probe['state_file']['path']
    state = json_value(root / state_name)
    require(state['electrons'] == old_electrons, 'Actual standalone/quartet probe disagreement')
    probe['electrons'] = old_electrons + 1
    state['electrons'] = probe['electrons']
    replace(state_name, state)
    replace(quartet_name, quartet)
    def rebound(value):
        if isinstance(value, dict):
            result = {key: rebound(child) for key, child in value.items()}
            # Bindings always move their byte count together with their hash.
            original_sha = value.get('sha256')
            if isinstance(original_sha, str) and original_sha in replacements:
                before, after = replacements[original_sha]
                if 'bytes' in value:
                    require(value['bytes'] == before['bytes'], 'Inconsistent binding in semantic negative')
                    result['bytes'] = after['bytes']
                result['sha256'] = after['sha256']
            return result
        if isinstance(value, list):
            return [rebound(child) for child in value]
        if isinstance(value, str) and value in replacements:
            return replacements[value][1]['sha256']
        return value
    json_names = sorted(name for name in manifest['files'] if name.endswith('.json'))
    # The graph contains source -> state -> quartet -> run -> readout -> curve metadata.
    # Iteration also covers copied bindings in future optional G2/curated provenance.
    for unused in range(len(json_names) + 1):
        start = len(changes)
        for name in json_names:
            value = json_value(root / name)
            updated = rebound(value)
            if updated != value:
                replace(name, updated)
        if len(changes) == start:
            break
    else:
        raise ValueError('Binding graph failed to settle')
    manifest_path.write_bytes(json_bytes(manifest))
    # Establish complete byte/source/standalone rebinding independently of the physical guard.
    sys.path.insert(0, str(root))
    from reader import verify_inventory, verify_native
    verified = verify_inventory(root)
    for role in json_value(root / 'inputs.json')['runs']:
        verify_native(root, verified, role)
    return dict(case='F', charge='q1', state=state_name, R_A=probe['R_A'], electron_change=1,
                old_electrons=old_electrons, wrong_electrons=probe['electrons'], rebinding=changes,
                integrity_and_standalone_checks_passed=True,
                purpose='Negative variant of the actual archive only; no fabricated positive Native fixture.')


def main(archive, output):
    archive, output = Path(archive), Path(output)
    require(archive.is_file() and not archive.is_symlink(), 'Actual regular archive required')
    require(not output.exists() and not output.is_symlink(), 'Refusing existing proof output')
    output.mkdir(parents=True, exist_ok=False)
    archive_binding = digest(archive.read_bytes())
    results, findings = {}, []
    proof = dict(schema='hydride-screening-portability-proof-v1', passed=False, archive=archive_binding,
                 results=results, findings=findings, publishing=False)
    def save():
        data = json_bytes(proof)
        scan_text('portability-proof.json', data)
        (output / 'portability-proof.json').write_bytes(data)
    try:
        actual = extract(archive, output / 'actual')
        proof.update(manifest=digest((actual / 'manifest.json').read_bytes()),
                     archive_members=sum(path.is_file() for path in actual.rglob('*')),
                     actual_positive_fixture='unmodified extracted archive')
        numpy_spec = importlib.util.find_spec('numpy')
        require(numpy_spec is not None and numpy_spec.origin is not None, 'Installed NumPy required')
        numpy_path = str(Path(numpy_spec.origin).resolve().parent.parent)
        command = [sys.executable, '-I', '-S', str(actual / 'portability_proof.py'), '--worker']
        def run(label, target, expected):
            try:
                completed = subprocess.run(command + [str(target), numpy_path], capture_output=True, check=False, timeout=180)
            except subprocess.TimeoutExpired as error:
                stdout, stderr = error.stdout or b'', error.stderr or b''
                returncode = None
            else:
                stdout, stderr = completed.stdout, completed.stderr
                returncode = completed.returncode
            for suffix, data in (('.stdout', stdout), ('.stderr', stderr),
                                 ('.exit', (str(returncode) + '\n').encode()),
                                 ('.argv.json', json_bytes(['python', '-I', '-S', 'ACTUAL_ARCHIVE/portability_proof.py',
                                                           '--worker', label + '/' + ARCHIVE_ROOT, 'INSTALLED_NUMPY_DIRECTORY']))):
                (output / (label + suffix)).write_bytes(data)
            accepted = returncode == 0
            matched = returncode is not None and (accepted == (expected == 'accept'))
            results[label] = dict(returncode=returncode, stdout=digest(stdout), stderr=digest(stderr),
                                  expected=expected, disposition='accepted' if accepted else 'refused', matched=matched)
            findings.append(dict(check=label, disposition='passed' if matched else 'failed'))
            save()
            require(matched, label + ': unexpected CLI disposition')
            return stdout, stderr.decode('utf-8')
        stdout, unused = run('actual-replay', actual, 'accept')
        report = json.loads(stdout)
        byte_root = output / 'byte-tamper' / ARCHIVE_ROOT
        byte_root.parent.mkdir()
        shutil.copytree(actual, byte_root)
        target = 'curves/F_q1_BO_pair.json'
        with (byte_root / target).open('ab') as handle:
            handle.write(b' ')
        unused, stderr = run('byte-tamper', byte_root, 'refuse')
        expected = 'EXPORTED byte integrity failed: ' + target
        require(expected in stderr, 'Byte negative did not reach its integrity guard')
        results['byte-tamper']['refusal_reason'] = expected
        semantic_root = output / 'extra-q1-probe-electron' / ARCHIVE_ROOT
        semantic_root.parent.mkdir()
        shutil.copytree(actual, semantic_root)
        mutation = semantic_variant(semantic_root)
        (output / 'semantic-mutation.json').write_bytes(json_bytes(mutation))
        proof['semantic_mutation'] = mutation
        unused, stderr = run('extra-q1-probe-electron', semantic_root, 'refuse')
        expected = 'Probe electron count'
        require(expected in stderr and 'EXPORTED byte integrity failed' not in stderr and
                'source binding mismatch' not in stderr, 'Semantic negative did not reach the physical electron guard')
        results['extra-q1-probe-electron']['refusal_reason'] = expected
        require(digest(archive.read_bytes()) == archive_binding, 'Original archive changed during proof')
        for path in sorted(output.iterdir()):
            if path.is_file():
                scan_text(path.name, path.read_bytes())
        proof.update(passed=True, qualified=report['qualified'], rows=report['rows'],
                     F_controls=report['F_controls'], ratios_R01=report['ratios_R01'], curves=report['curves'],
                     blocked_imports=sorted(FORBIDDEN), isolated_python=True,
                     dependency_scope='NumPy and stdlib only')
        findings.append(dict(check='proof text privacy and unchanged original archive', disposition='passed'))
        save()
    except Exception as error:
        proof['passed'] = False
        findings.append(dict(check='proof completion', disposition='failed', error_type=type(error).__name__))
        save()
        raise
    print(json_bytes(proof).decode(), end='')
    return proof


if __name__ == '__main__':
    try:
        if len(sys.argv) == 4 and sys.argv[1] == '--worker':
            worker(sys.argv[2], sys.argv[3])
        else:
            parser = argparse.ArgumentParser(description=__doc__)
            parser.add_argument('archive', type=Path)
            parser.add_argument('output', type=Path)
            args = parser.parse_args()
            main(args.archive, args.output)
    except (ValueError, KeyError, TypeError, OSError, UnicodeError, ImportError, tarfile.TarError) as error:
        print('REFUSED: ' + str(error), file=sys.stderr)
        sys.exit(2)
