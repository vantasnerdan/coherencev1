#!/usr/bin/env python3
"""Actual isolated CLI replay plus byte tamper, wrong-source and wrong-cell refusals."""
import argparse
import builtins
import importlib.abc
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def worker(root):
    original_import = builtins.__import__
    forbidden = {'scipy', 'ase', 'gpaw'}
    def guarded(name, *args, **kwargs):
        if name.split('.')[0] in forbidden:
            raise ImportError('Portability proof blocks ' + name)
        return original_import(name, *args, **kwargs)
    builtins.__import__ = guarded
    for name in sorted(forbidden):
        try:
            __import__(name)
        except ImportError:
            pass
        else:
            raise RuntimeError('Import blocker failed')
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from reader import replay
    print(json.dumps(replay(root), indent=2, sort_keys=True, allow_nan=False))
    require_no_forbidden = not any(name.split('.')[0] in forbidden for name in sys.modules)
    if not require_no_forbidden:
        raise RuntimeError('Forbidden dependency loaded')


def main(root, output):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from capsule_io import digest, require
    output = Path(output); output.mkdir(parents=True, exist_ok=False)
    command = [sys.executable, '-I', str(Path(__file__).resolve()), '--worker']
    results = {}
    def run(label, target, expect_success):
        completed = subprocess.run(command + [str(target)], capture_output=True)
        (output / (label + '.stdout')).write_bytes(completed.stdout)
        (output / (label + '.stderr')).write_bytes(completed.stderr)
        (output / (label + '.exit')).write_text(str(completed.returncode) + '\n')
        # Paths are deliberately logical in exported argv; executed argv is independently
        # recorded by the private immutable delivery report.
        (output / (label + '.argv.json')).write_text(json.dumps(['python', '-I', 'portability_proof.py', '--worker', '<CAPSULE>']) + '\n')
        results[label] = dict(returncode=completed.returncode, stdout=digest(completed.stdout), stderr=digest(completed.stderr), expected='accept' if expect_success else 'refuse')
        require((completed.returncode == 0) == expect_success, label + ': unexpected CLI disposition')
        return completed
    actual = run('actual-replay', root, True)
    report = json.loads(actual.stdout)
    with tempfile.TemporaryDirectory(prefix='collision-proof-') as temporary:
        temporary = Path(temporary)
        for label in ('byte-tamper', 'wrong-source', 'wrong-image-cell'):
            target = temporary / label
            shutil.copytree(root, target)
            manifest_path = target / 'manifest.json'
            manifest = json.loads(manifest_path.read_bytes())
            if label == 'byte-tamper':
                name = 'curves/baseline_BO_pair.json'
                with (target / name).open('ab') as handle:
                    handle.write(b' ')
            elif label == 'wrong-source':
                name = 'runs/base-occ/calculation/bo/excitation.json'
                shutil.copyfile(target / 'runs/gamma-occ/calculation/bo/excitation.json', target / name)
                manifest['files'][name]['EXPORTED'] = digest((target / name).read_bytes())
                manifest_path.write_text(json.dumps(manifest))
            else:
                name = 'controls/image-2x/calculation/occupied-probe/state.json'
                state = json.loads((target / name).read_bytes())
                state['cell_A'][0][0] += .01
                (target / name).write_text(json.dumps(state))
                manifest['files'][name]['EXPORTED'] = digest((target / name).read_bytes())
                manifest_path.write_text(json.dumps(manifest))
            refusal = run(label, target, False)
            text = refusal.stderr.decode()
            expected = {'byte-tamper':'EXPORTED byte integrity failed', 'wrong-source':'ORIGINAL source binding mismatch', 'wrong-image-cell':'Approved doubled cell'}[label]
            require(expected in text, 'Wrong refusal reason: ' + label)
            results[label]['refusal_reason'] = expected
    proof = dict(schema='pdd-actual-numpy-portability-proof-v1', passed=True,
                 blocked_imports=['scipy', 'ase', 'gpaw'], isolated_python=True, results=results,
                 full_trajectory_rows=report['full_trajectory_rows'], replay_differences=report['replay_differences'],
                 controls=report['controls'], encounter_excitation_dt_checked=False, halfdt_drag='pending; may be second asset')
    (output / 'portability-proof.json').write_text(json.dumps(proof, indent=2, sort_keys=True, allow_nan=False) + '\n')
    print(json.dumps(proof, indent=2, sort_keys=True))


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--worker':
        try:
            worker(sys.argv[2])
        except (ValueError, KeyError, OSError) as error:
            print('REFUSED: ' + str(error), file=sys.stderr)
            sys.exit(2)
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('capsule'); parser.add_argument('output')
        args = parser.parse_args()
        main(Path(args.capsule), args.output)
