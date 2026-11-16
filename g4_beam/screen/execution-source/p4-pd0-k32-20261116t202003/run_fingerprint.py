"""One serial native fingerprint; interrupted work is never consumer-eligible."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time

import gpaw

from screen_ground import ground_state
from screen_io import manifest, require, write_json
from screen_response import response_fingerprint
from screen_saved_ground import admit_saved_ground


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--setup-root', type=Path, required=True)
    parser.add_argument('--phase', required=True)
    parser.add_argument('--reuse-ground', type=Path)
    args = parser.parse_args()
    config = json.loads(args.input.read_bytes())
    require(config['schema'] == 'coherence-screen-fingerprint-input-v1',
            'Wrong fingerprint input schema')
    require(os.sched_getaffinity(0) == {config['cpu']},
            'Native affinity differs from the declared owned single core')
    require(gpaw.mpi.world.size == 1, 'Fingerprint must have exactly one rank')
    gpaw.setup_paths[:] = [str(args.setup_root)]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).resolve().parent
    producer = [manifest(source/name, 'source/'+name) for name in
                ['screen_io.py', 'screen_ground.py', 'screen_response.py',
                 'run_fingerprint.py', 'screen_saved_ground.py']]
    inherited_ground = None
    if args.reuse_ground is not None:
        inherited_ground = admit_saved_ground(config, args.reuse_ground, output)
    require(config.get('reuse_ground') == inherited_ground,
            'Requested reused-ground byte identity differs from actual source')
    initial = dict(phase=args.phase, config=config,
                   started_utc=datetime.now(timezone.utc).isoformat(),
                   producer_source=producer, input=manifest(args.input),
                   source_ground=inherited_ground,
                   setups=[manifest(args.setup_root/name, 'setups/'+name)
                           for name in config['setup_files']],
                   affinity=list(os.sched_getaffinity(0)))
    write_json(output/'producer.json', initial)
    start = time.monotonic()
    print(f'SCREEN_FINGERPRINT_START phase={args.phase}', flush=True)
    if inherited_ground is None:
        ground_state(config, output)
    print('SCREEN_FINGERPRINT_GROUND_COMPLETED', flush=True)
    response_fingerprint(config, output)
    write_json(output/'completed.json', dict(
        phase=args.phase, completed_utc=datetime.now(timezone.utc).isoformat(),
        elapsed_wall_s=time.monotonic()-start,
        ground=manifest(output/'ground.json'),
        response=manifest(output/'response.json'),
        producer=manifest(output/'producer.json'),
        eligible=False,
        reason='A fingerprint does not establish the requested material tolerances'))
    print('SCREEN_FINGERPRINT_NUMERICS_COMPLETED', flush=True)


if __name__ == '__main__':
    main()
