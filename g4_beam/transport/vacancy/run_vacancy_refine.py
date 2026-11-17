"""Execute declared numerical controls on completed source-bound geometries only."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
if os.sched_getaffinity(0) != {7}:
    raise RuntimeError('Vacancy numerical controls require CPU7')
import argparse
import json
from pathlib import Path
import time
from run_vacancy import binding


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    config = json.loads(Path(args.input).read_text())
    output = Path(args.output)
    output.mkdir(exist_ok=False)
    source = Path(config['parent_directory'])
    setup_dir = Path(config['setup_directory']).resolve()
    os.environ['GPAW_SETUP_PATH'] = str(setup_dir)
    sources = [binding(args.input)] + [binding(Path(__file__).with_name(name)) for name in
               ('run_vacancy_refine.py', 'vacancy_refine.py', 'vacancy_geometry.py', 'run_vacancy.py')]
    for name, expected in config['setup_sha256'].items():
        actual = binding(setup_dir / name)
        if actual['sha256'] != expected:
            raise ValueError('Selected PAW byte mismatch: ' + name)
        sources.append(actual)
    for case in config['cases']:
        for name in ('primitive.json', 'final.extxyz', 'final.gpw'):
            sources.append(binding(source / config.get('parent_case', case) / name))
    if 'seed_reference' in config:
        sources.append(binding(config['seed_reference']))
    provenance = dict(schema='coherence-vacancy-execution-v1', sources=sources,
                      affinity=sorted(os.sched_getaffinity(0)),
                      cgroup=Path('/proc/self/cgroup').read_text(),
                      start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                      threads={key: value for key, value in os.environ.items() if key.endswith('NUM_THREADS')})
    (output / 'execution.json').write_text(json.dumps(provenance, indent=2) + '\n')
    from vacancy_refine import refine
    for case in config['cases']:
        print('BEGIN_CASE', case, flush=True)
        start = time.monotonic()
        result = refine(config, case, source / config.get('parent_case', case), output / case)
        print('END_CASE', case, 'seconds', time.monotonic()-start,
              'parent_stationarity_retained', result['converged_relaxation'], flush=True)
    (output / 'completed.json').write_text(json.dumps(dict(cases=config['cases'], completed=True)) + '\n')


if __name__ == '__main__':
    main()
