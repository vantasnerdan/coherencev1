#!/usr/bin/env python3
"""Build finite local capsule, then seal with actual portability proof; no Native/network."""
import argparse
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import tarfile
from capsule_io import digest, require

PRIVATE_PATH = re.compile(r'/(?:home|root|tmp|run/user|mnt|workspace|opt)/[^\s\"\'<>`,;\)\]\}]+')


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()


def sanitize(value, changes, trail=''):
    if isinstance(value, dict):
        return {k: sanitize(v, changes, trail + '/' + k) for k, v in value.items()}
    if isinstance(value, list):
        return [sanitize(v, changes, trail + '/' + str(i)) for i, v in enumerate(value)]
    if isinstance(value, str):
        transformed, count = PRIVATE_PATH.subn('<PRIVATE_PATH>', value)
        if count:
            changes.append(dict(json_pointer=trail, private_path_values_removed=count))
        return transformed
    return value


class Builder:
    def __init__(self, root):
        self.root = root
        self.files = {}

    def add(self, source, destination, exact=False):
        raw = source.read_bytes(); exported = raw; changes = []
        if source.suffix == '.json':
            original = json.loads(raw)
            transformed = sanitize(original, changes)
            if changes:
                exported = json_bytes(transformed)
            require(sanitize(original, []) == json.loads(exported), 'Non-path scientific transformation')
        elif source.suffix != '.npz':
            text = raw.decode('utf-8')
            transformed, count = PRIVATE_PATH.subn('<PRIVATE_PATH>', text)
            if count:
                changes = [dict(private_path_values_removed=count)]
                exported = transformed.encode()
        require(not exact or exported == raw, 'Required exact source contains private paths')
        self.generated(destination, exported, original=digest(raw), changes=changes)

    def generated(self, destination, data, original=None, changes=None):
        require(destination not in self.files, 'Duplicate capsule destination')
        path = self.root / destination
        path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(data)
        self.files[destination] = dict(ORIGINAL=original or digest(data), EXPORTED=digest(data),
            transformation=changes or [], provenance='original-byte inventory; EXPORTED differs only for private path values' if original else 'generated delivery source/evidence')

    def manifest(self, omitted):
        data = dict(schema='pdd-collision-capsule-v1', files=self.files, private_orbitals=omitted,
            provenance='ORIGINAL hashes/sizes bind original scientific epochs, never rewritten to EXPORTED hashes. Original private path values are not distributed. Arrays/numeric values remain unchanged. Manifest integrity is externally bound by tar SHA256; it is not a signature.',
            omissions=['restricted 2keV/head-on frozen snapshot, not a computed locked state', 'interior R0.1–1.6 only; no U0, prepared state or rate', 'encounter excitation NOT dt-checked; finite-path drag with halfdt pending'],
            reader_command='python reader.py .', portability_command='python portability_proof.py CAPSULE_DIRECTORY OUTPUT_DIRECTORY',
            deadline='2026-11-20T05:30:00Z', native_invocations=[], budget_cpu_core_min_cap=20)
        (self.root / 'manifest.json').write_bytes(json_bytes(data))


def build(args):
    root = Path(args.capsule)
    require(not root.exists(), 'Capsule destination already exists')
    root.mkdir(parents=True)
    private = Path(args.private_root); canonical = Path(args.source_root)
    builder = Builder(root); omitted = []
    runs = [(phase + '-' + role, 'p6collision-' + phase + '-' + role + '-001') for phase in ('base', 'gamma') for role in ('occ', 'vac')]
    runs += [('h016', 'p6collision-h016-static-quartet-001'), ('image-2x', 'p6collision-image-2x-quartet-001')]
    for name, run in runs:
        directory = private / 'g6_collision' / run
        prefix = ('runs/' if name not in ('h016', 'image-2x') else 'controls/') + name
        paths = [Path(base) / filename for base, _, names in os.walk(directory, followlinks=True) for filename in names]
        for path in sorted(paths):
            relative = path.relative_to(directory)
            if path.suffix == '.gpw':
                hasher = hashlib.sha256()
                with path.open('rb') as handle:
                    for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b''):
                        hasher.update(chunk)
                omitted.append(dict(relative_original=prefix + '/' + str(relative), bytes=path.stat().st_size,
                                    sha256=hasher.hexdigest(), disposition='raw orbital retained privately; no raw GPW in capsule'))
            elif relative.parts[0] in ('calculation', 'source') or path.name in ('input.json', 'source-bindings.sha256', 'native.log', 'launch.stdout', 'launch.exit', 'launch.argv', 'watch.stdout', 'watch.exit', 'watch.argv', 'observer.stdout', 'observer.exit'):
                builder.add(path, prefix + '/' + str(relative))
    receipt_names = [
        'collision-base-periodic-pair-readout-001', 'collision-gamma-periodic-pair-readout-001',
        'collision-displacement-pair-admission-001', 'collision-gamma-vac-admission-001',
        'collision-base-pair-invariant-admission', 'collision-h016-static-comparison-001',
        'collision-h016-comparison-check-001', 'collision-image-comparison-001',
        'collision-qualified-consumer-proof-001', 'collision-shared-state-h-replay-001',
        'collision-native-nuclear-kernel-001', 'collision-native-nuclear-kernel-002',
        'collision-bare-green-admission', 'collision-gaussian-gauge-probe', 'collision-work-screening-identity',
        'collision-static-field-admission', 'collision-ae-density-admission', 'collision-production-ground-reuse',
        'collision-watch-receipt-retrospective-clarification', 'collision-pilot-cost-disposition',
        'collision-canonical-producer-bindings-20261120', 'collision-unprepared-main-refusal-001',
        'collision-periodic-pair-main-refusal', 'collision-pair-wrong-role-001',
        'collision-baseline-main-consumer-001', 'collision-gamma-main-consumer-001',
        'collision-step-reader-diagnosis', 'collision-readout-repair']
    for stem in receipt_names:
        found = list((private / 'receipts').glob(stem + '.*'))
        require(bool(found), 'Declared receipt unavailable: ' + stem)
        for path in sorted(found):
            builder.add(path, 'receipts/' + path.name)
    for phase in ('baseline', 'gamma'):
        builder.add(private / 'g6_collision/bo-qualified-20261120' / (phase + '_BO_pair.json'), 'curves/' + phase + '_BO_pair.json', exact=True)
    relaxed = private / 'g6_collision/sources/relaxed-pdd-original.json'
    builder.add(relaxed, 'sources/relaxed-pdd-original.json')
    for path in sorted(Path(__file__).parent.glob('*.py')):
        builder.add(path, path.name, exact=True)
    for path in sorted(canonical.glob('*.py')):
        builder.add(path, 'oracle/campaign/' + path.name)
    framework = canonical.parent.parent / 'src/substrate_framework'
    for name in ('collision_bare_reference.py', 'collision_observables.py', 'fixed_projector_contractions.py'):
        builder.add(framework / name, 'oracle/framework/' + name)
    readme = ('The imposed neighbouring-D displacement of +0.0768 A along [100] changes the R=0.1 A screening lowering from 40.5497996 to 39.8626563 eV (−0.6871433 eV), below the 5 eV accuracy target and not a resolved enhancement.\n'
              '- Restricted 2 keV/head-on frozen snapshot; not a computed locked state.\n'
              '- Interior R=0.1–1.6 A only; no U0, prepared state or rate.\n'
              '- Encounter excitation is NOT dt-checked; vacancy drag is finite-path and halfdt is pending, possibly a second asset.\n')
    builder.generated('README.txt', readme.encode())
    builder.manifest(omitted)
    print(json.dumps(dict(capsule=str(root), files=len(builder.files), private_gpw=len(omitted))))


def seal(args):
    root = Path(args.capsule); manifest = json.loads((root / 'manifest.json').read_bytes())
    proof = Path(args.proof)
    for path in sorted(proof.glob('*')):
        if path.is_file():
            name = 'proof/' + path.name; target = root / name; target.parent.mkdir(exist_ok=True)
            raw = path.read_bytes()
            require(not PRIVATE_PATH.search(raw.decode('utf-8')), 'Private path in proof')
            target.write_bytes(raw)
            manifest['files'][name] = dict(ORIGINAL=digest(raw), EXPORTED=digest(raw), transformation=[], provenance='actual isolated portability proof')
    (root / 'manifest.json').write_bytes(json_bytes(manifest))
    findings = []
    for path in sorted(root.rglob('*')):
        if not path.is_file() or path.suffix == '.npz':
            continue
        text = path.read_text()
        if PRIVATE_PATH.search(text) or re.search(r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}', text):
            findings.append(str(path.relative_to(root)))
    require(not findings, 'Privacy scan failed: ' + repr(findings))
    destination = Path(args.tar)
    with destination.open('wb') as output:
        with gzip.GzipFile(filename='', mode='wb', fileobj=output, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode='w') as archive:
                for path in sorted(root.rglob('*')):
                    if path.is_file():
                        raw = path.read_bytes(); info = tarfile.TarInfo('collision-delivery/' + str(path.relative_to(root)))
                        info.size = len(raw); info.mtime = 0; info.mode = 0o644
                        archive.addfile(info, io.BytesIO(raw))
    receipt = dict(schema='pdd-collision-sealed-asset-v1', tar=digest(destination.read_bytes()),
                   manifest=digest((root / 'manifest.json').read_bytes()), privacy_scan=dict(passed=True, findings=[], scope='all exported text including proofs; NPZs are numeric density/operator arrays; no GPW/PDF'),
                   member_count=len(manifest['files']) + 1)
    Path(args.receipt).write_bytes(json_bytes(receipt)); print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--capsule', required=True)
    parser.add_argument('--private-root'); parser.add_argument('--source-root')
    parser.add_argument('--seal', action='store_true'); parser.add_argument('--proof')
    parser.add_argument('--tar'); parser.add_argument('--receipt')
    args = parser.parse_args()
    seal(args) if args.seal else build(args)
