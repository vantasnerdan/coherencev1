"""Verify portable force records against actual compact GPW/geometry banks."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
from force_bank_results import verify_force_record
from surface_io import save


def read_bound(base, reference):
    relative = Path(reference['path'])
    if relative.is_absolute():
        raise ValueError('Portable bank references must be relative')
    raw = (base/relative).read_bytes()
    if len(raw) != reference['bytes'] or sha256(raw).hexdigest() != reference['sha256']:
        raise ValueError('Portable force-record byte identity changed')
    return json.loads(raw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Verification output is immutable')
    path = args.manifest.resolve()
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != args.sha256:
        raise ValueError('Native bank manifest identity changed')
    manifest = json.loads(raw)
    if manifest['schema'] != 'coherence-portable-native-force-banks-v2':
        raise ValueError('Explicit portable force-bank payload schema required')
    count = 0
    for case in manifest['cases'].values():
        for reference in case['records'].values():
            record = read_bound(path.parent, reference)
            if case['payload_key'] is not None:
                record = record[case['payload_key']]
            for field in ('gpw', 'geometry'):
                relative = Path(record[field]['path'])
                if relative.is_absolute():
                    raise ValueError('Portable Native bank path must be relative')
                record[field] = {**record[field], 'path': str(path.parent/relative)}
            verify_force_record(record)
            count += 1
    if count != manifest['record_count']:
        raise ValueError('Native manifest record count differs')
    report = dict(schema='coherence-portable-native-bank-verification-v1',
        Native_launched=False, records_verified=count,
        source_manifest=dict(bytes=len(raw), sha256=args.sha256))
    save(args.output, report)
    print('NATIVE_BANKS_VERIFIED', count, flush=True)


if __name__ == '__main__':
    main()
