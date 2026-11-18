"""Bounded static family driver; all actual DFT/RPA stays in qualified producers."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time

import gpaw
from gpaw.mpi import world
from prepare_static_ground import prepare_saved_ground
from screen_ground import ground_state
from screen_io import manifest, require, write_json
from screen_response import response_fingerprint


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    specification = json.loads(args.input.read_text())
    require(world.size == 1, 'This family is serial')
    for record in specification['bound_sources']:
        actual = manifest(Path(record['path']))
        require(actual['bytes'] == record['bytes'] and actual['sha256'] == record['sha256'],
                'Registered source/setup bytes changed: ' + record['path'])
    gpaw.setup_paths[:] = [specification['setup_directory']]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    config = specification['configuration']
    write_json(output / 'input.json', config)
    start = time.monotonic()
    stages = []
    status = 'failed'
    try:
        begin = time.monotonic()
        print('STATIC_GROUND_START', config['case'], flush=True)
        if specification['ground_mode'] == 'saved':
            prepare_saved_ground(specification, output)
        else:
            require(specification['ground_mode'] == 'cold', 'Unknown ground mode')
            ground_state(config, output)
        elapsed = time.monotonic() - begin
        stages.append(dict(stage='saved_or_cold_ground', elapsed_wall_s=elapsed))
        write_json(output / 'stage-cost-ground.json', stages)
        print('STATIC_GROUND_COMPLETE', config['case'], elapsed, flush=True)
        require(elapsed <= 3 * specification['expected_ground_wall_s'],
                'Ground preparation exceeds3x expected cost; stop before response scaling')
        begin = time.monotonic()
        print('STATIC_RESPONSE_START', config['case'], flush=True)
        response_fingerprint(config, output)
        stages.append(dict(stage='response_family', elapsed_wall_s=time.monotonic() - begin))
        write_json(output / 'stage-cost-final.json', stages)
        status = 'complete'
        print('STATIC_FAMILY_COMPLETE', config['case'], flush=True)
    finally:
        write_json(output / 'terminal.json', dict(status=status, stages=stages,
            elapsed_wall_s=time.monotonic() - start, process_cpu_s=time.process_time(),
            finished_at=datetime.now(timezone.utc).isoformat(), argv=sys.argv,
            recovery='Complete ground/point banks only; no in-flight SCF/response reuse promised'))


if __name__ == '__main__':
    main()
