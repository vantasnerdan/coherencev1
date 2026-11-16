"""Serial, isolated case checkpoints inside one bounded native owning unit."""
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

SOURCE_FILES = ['screen_io.py', 'screen_cases.py', 'screen_dos.py',
                'screen_relax.py', 'screen_ground.py', 'screen_response.py',
                'screen_saved_ground.py', 'run_bulk_batch.py']


def run_case(directory, config, source, setups, reuse):
    require(config['kind'] in ['dos', 'rpa'], 'Unknown commissioned case kind')
    inherited = None
    if reuse is not None:
        inherited = admit_saved_ground(config, reuse, directory)
    require(config.get('reuse_ground') == inherited,
            'Reused ground identity differs from the declared source')
    producer = dict(config=config,
                    started_utc=datetime.now(timezone.utc).isoformat(),
                    input=manifest(directory/'input.json'), source_ground=inherited,
                    producer_source=[manifest(source/name, 'source/'+name)
                                     for name in SOURCE_FILES],
                    setups=[manifest(setups/name, 'setups/'+name)
                            for name in config['setup_files']],
                    affinity=list(os.sched_getaffinity(0)))
    write_json(directory/'producer.json', producer)
    start = time.monotonic()
    if inherited is None:
        ground_state(config, directory)
    print(f'SCREEN_GROUND_COMPLETED case={config["case"]}', flush=True)
    artifacts = dict(ground=manifest(directory/'ground.json'),
                     producer=manifest(directory/'producer.json'))
    if config['kind'] == 'rpa':
        require(len(config['responses']) >= 2, 'A screening-length case needs multiple q')
        response_fingerprint(config, directory)
        artifacts['response'] = manifest(directory/'response.json')
    write_json(directory/'completed.json', dict(
        schema='coherence-screen-bulk-completion-v1', kind=config['kind'],
        case=config['case'], completed_utc=datetime.now(timezone.utc).isoformat(),
        elapsed_wall_s=time.monotonic()-start, **artifacts,
        eligible=False, reason='Requires independent observable-level qualification'))
    print(f'SCREEN_BULK_CASE_COMPLETED case={config["case"]}', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--setup-root', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_bytes())
    require(os.sched_getaffinity(0) == {plan['cpu']}, 'Wrong owned CPU affinity')
    require(gpaw.mpi.world.size == 1, 'Bulk batch requires exactly one rank')
    gpaw.setup_paths[:] = [str(args.setup_root)]
    source = Path(__file__).resolve().parent
    root = args.plan.resolve().parent
    for entry in plan['cases']:
        name = entry['directory']
        require(Path(name).name == name, 'Case path escapes batch')
        directory = root/name
        config = json.loads((directory/'input.json').read_bytes())
        require(config['cpu'] == plan['cpu'], 'Case and batch CPU differ')
        run_case(directory, config, source, args.setup_root, entry.get('reuse'))
    write_json(root/'batch-completed.json', dict(
        phase=plan['phase'], completed_utc=datetime.now(timezone.utc).isoformat(),
        cases=[dict(directory=e['directory'],
                    completed=manifest(root/e['directory']/'completed.json'))
               for e in plan['cases']]))


if __name__ == '__main__':
    main()
