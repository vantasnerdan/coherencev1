"""Execute the declared cases once, binding raw source/setup bytes before GPAW."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
if os.sched_getaffinity(0) != {7}:
    raise RuntimeError('Vacancy work requires exclusive CPU7 affinity')
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time


def binding(path):
    data = Path(path).read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='inputs.json')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    config = json.loads(Path(args.input).read_text())
    output = Path(args.output)
    output.mkdir(exist_ok=False)
    setup_dir = Path(config['setup_directory']).resolve()
    os.environ['GPAW_SETUP_PATH'] = str(setup_dir)
    sources = [binding(args.input)] + [binding(Path(__file__).with_name(name)) for name in
               ('run_vacancy.py', 'vacancy_geometry.py', 'vacancy_calculate.py')]
    for name, expected in config['setup_sha256'].items():
        actual = binding(setup_dir / name)
        if actual['sha256'] != expected:
            raise ValueError('Selected PAW byte mismatch: ' + name)
        sources.append(actual)
    provenance = dict(schema='coherence-vacancy-execution-v1', sources=sources,
                      affinity=sorted(os.sched_getaffinity(0)),
                      cgroup=Path('/proc/self/cgroup').read_text(),
                      start_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                      threads={key: value for key, value in os.environ.items() if key.endswith('NUM_THREADS')})
    (output / 'execution.json').write_text(json.dumps(provenance, indent=2) + '\n')
    from vacancy_calculate import calculate
    for case in config['cases']:
        print('BEGIN_CASE', case, flush=True)
        start = time.monotonic()
        result = calculate(config, case, output / case)
        print('END_CASE', case, 'seconds', time.monotonic() - start,
              'converged', result['converged_relaxation'], flush=True)
        if not result['converged_relaxation']:
            raise RuntimeError('Relaxation plateau/step cap: no unchanged retry')
    (output / 'completed.json').write_text(json.dumps(dict(cases=config['cases'], completed=True)) + '\n')


if __name__ == '__main__':
    main()
