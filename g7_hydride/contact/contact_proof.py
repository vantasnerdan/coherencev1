"""Replay the actual contact archive in NumPy-only isolation and exercise physical refusals."""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

from portability_proof import extract
from delivery_io import digest, json_bytes, json_value, require, scan_text
from hydride_replay import bound_file


def replace_json(path, value):
    path.write_bytes(json_bytes(value))

def verify_bindings(root):
    manifest = json_value(root / 'contact-sidecars.json')
    for section in ('sources', 'curves', 'sidecars'):
        for entry in manifest[section]:
            bound_file(root, entry, 'Repaired physical-variant byte integrity')



def rebind_sources(root, names):
    path = root / 'contact-sidecars.json'
    manifest = json_value(path)
    for entry in manifest['sources']:
        if entry['path'] in names:
            entry.update(digest((root / entry['path']).read_bytes()))
    replace_json(path, manifest)


def extra_electron(root):
    prefix = 'sources/F_q1/coarse/'
    native_path = root / prefix / 'contact-native.json'
    native = json_value(native_path)
    state = native['states']['occupied-probe-He']
    relative = state['state_file']['path']
    state_path = root / prefix / relative
    standalone = json_value(state_path)
    old = standalone['electrons']
    standalone['electrons'] += 1
    state['electrons'] = standalone['electrons']
    replace_json(state_path, standalone)
    state['state_file'].update(digest(state_path.read_bytes()))
    replace_json(native_path, native)
    rebind_sources(root, {prefix + relative, prefix + 'contact-native.json'})
    bound_file(root / prefix, state['state_file'], 'Repaired standalone Native binding')
    return dict(state=prefix + relative, old_electrons=old,
                wrong_electrons=standalone['electrons'], standalone_and_outer_hashes_repaired=True)


def wrong_contact(root):
    path = root / 'F_q1_contact.json'
    value = json_value(path)
    old = value['signedW0_eV']
    value['signedW0_eV'] += 1
    replace_json(path, value)
    manifest_path = root / 'contact-sidecars.json'
    manifest = json_value(manifest_path)
    for entry in manifest['sidecars']:
        if entry['path'] == path.name:
            entry.update(digest(path.read_bytes()))
    replace_json(manifest_path, manifest)
    return dict(old_signedW0_eV=old, wrong_signedW0_eV=value['signedW0_eV'],
                sidecar_hash_repaired=True)


def main(archive, output):
    archive, output = Path(archive).resolve(), Path(output).resolve()
    require(not output.exists(), 'Contact proof output must be new')
    original = digest(archive.read_bytes())
    output.mkdir(parents=True)
    numpy = importlib.util.find_spec('numpy')
    require(numpy is not None and numpy.submodule_search_locations, 'Installed NumPy required')
    numpy_path = str(Path(next(iter(numpy.submodule_search_locations))).parent)
    base = extract(archive, output / 'actual-archive')
    results = {}
    proof = dict(schema='hydride-contact-portability-proof-v1', archive=original,
                 positive_fixture='unmodified actual extracted archive', results=results, passed=False,
                 proof_source=digest(Path(__file__).read_bytes()), publishing=False)

    def save():
        data = json_bytes(proof)
        scan_text('contact-portability-proof.json', data)
        (output / 'contact-portability-proof.json').write_bytes(data)

    def execute(name, root, expected):
        command = [sys.executable, '-I', '-S', str(root / 'portability_proof.py'),
                   '--contact-worker', str(root), numpy_path]
        (output / (name + '.argv.json')).write_bytes(json_bytes(command))
        run = subprocess.run(command, capture_output=True, timeout=90)
        (output / (name + '.stdout')).write_bytes(run.stdout)
        (output / (name + '.stderr')).write_bytes(run.stderr)
        (output / (name + '.exit')).write_text(str(run.returncode) + '\n')
        accepted = run.returncode == 0
        entry = dict(returncode=run.returncode, expected=expected,
                     disposition='accepted' if accepted else 'refused',
                     stdout=digest(run.stdout), stderr=digest(run.stderr))
        results[name] = entry
        save()
        require(accepted == (expected == 'accept'), 'Unexpected actual contact replay disposition: ' + name)
        if accepted:
            value = json.loads(run.stdout)
            require(value['qualified'] and value['isolated_python'] and
                    value['dependency_scope'] == 'NumPy and stdlib only', 'Unqualified isolated replay')
            proof['contacts'] = value['contacts']
            proof['blocked_imports'] = value['blocked_imports']
        else:
            entry['refusal_reason'] = run.stderr.decode('utf-8').strip()
        save()
        return entry

    try:
        execute('actual-replay', base, 'accept')
        bad_byte = output / 'byte-tamper'
        shutil.copytree(base, bad_byte)
        path = bad_byte / 'curves/F_q1_BO_pair.json'
        with path.open('ab') as stream:
            stream.write(b'\n')
        execute('byte-tamper', bad_byte, 'refuse')
        wrong_ne = output / 'extra-q1-electron'
        shutil.copytree(base, wrong_ne)
        proof['electron_mutation'] = extra_electron(wrong_ne)
        verify_bindings(wrong_ne)
        proof['electron_mutation']['all_archive_bindings_checked'] = True
        execute('extra-q1-electron', wrong_ne, 'refuse')
        wrong_w0 = output / 'wrong-contact-W0'
        shutil.copytree(base, wrong_w0)
        proof['contact_mutation'] = wrong_contact(wrong_w0)
        verify_bindings(wrong_w0)
        proof['contact_mutation']['all_archive_bindings_checked'] = True
        execute('wrong-contact-W0', wrong_w0, 'refuse')
        require(digest(archive.read_bytes()) == original, 'Original archive changed during proof')
        proof.update(passed=True, findings=[dict(check=name, disposition='passed') for name in results],
                     dependency_scope='NumPy and stdlib only', isolated_python=True)
        save()
    except Exception as error:
        proof['failure_type'] = type(error).__name__
        save()
        raise
    print(json_bytes(proof).decode(), end='')
    return proof


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    try:
        main(args.archive, args.output)
    except (ValueError, KeyError, TypeError, OSError, ImportError, subprocess.TimeoutExpired) as error:
        print('REFUSED: ' + str(error), file=sys.stderr)
        sys.exit(2)
