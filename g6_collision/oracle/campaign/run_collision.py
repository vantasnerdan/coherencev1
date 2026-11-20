"""Run one isolated collision stage; the service manager owns runtime and exits."""
import argparse
import json
from pathlib import Path
import time
from screen_io import write_json, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', choices=['ground', 'trajectory', 'bo', 'all'], default='all')
    args = parser.parse_args()
    spec = json.loads(args.input.read_bytes())
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    if args.stage in ('ground', 'all'):
        from collision_ground import prepare
        prepare(spec, output/'ground')
    if args.stage in ('trajectory', 'all'):
        from collision_propagate import propagate
        propagate(spec, output/'ground', output/'trajectory')
    if args.stage in ('bo', 'all'):
        from collision_bo import bo_checkpoints
        bo_checkpoints(spec, output/'ground', output/'trajectory', output/'bo')
    write_json(output/f'completed-{args.stage}.json', dict(stage=args.stage, input=manifest(args.input),
                   elapsed_wall_s=time.monotonic()-start))
    print(f'COLLISION_STAGE_COMPLETE stage={args.stage}', flush=True)


if __name__ == '__main__':
    main()
