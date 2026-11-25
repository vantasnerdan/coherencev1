"""Verify the saved TiD2 capsule without the solver or original PAW/PDF bytes."""
import argparse
from math import isclose
from pathlib import Path

from dense_pd import read_json, require
from fluorite_tid2 import (SCHEMA, SOURCE_BYTES, admit_zero, binding, case_specs, check_readout,
                           cuts_readout, load_point, load_reference, pair_readout)


def replay_inputs(root):
    reference = read_json(root/'relax/reference.json')
    model = reference['model']
    sources = model['source_artifacts']
    require(len(sources) == len(SOURCE_BYTES) and
            {r['file']: (r['bytes'], r['sha256']) for r in sources} == SOURCE_BYTES,
            'Reference does not bind the registered inputs')
    require(model['code'] == binding(Path(__file__).with_name('fluorite_tid2.py')),
            'Reference does not bind this producer')
    reference, reference_source = load_reference(root, model)
    trials = [binding(p) for p in sorted((root/'relax').glob('trial-*.json'))]
    require(trials == reference['trial_sources'], 'Relaxation trial bytes differ')
    return reference, reference_source


def same_readout(computed, saved):
    """Allow floating-point evaluation noise, not changed topology or provenance."""
    if type(computed) is not type(saved):
        return False
    if isinstance(computed, float):
        return isclose(computed, saved, rel_tol=0.0, abs_tol=1e-9)
    if isinstance(computed, dict):
        return computed.keys() == saved.keys() and all(same_readout(v, saved[k]) for k, v in computed.items())
    if isinstance(computed, list):
        return len(computed) == len(saved) and all(same_readout(a, b) for a, b in zip(computed, saved))
    return computed == saved


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('native', type=Path)
    args = p.parse_args()
    reference, reference_source = replay_inputs(args.native)
    points = {case: load_point(args.native, reference, reference_source, case) for case in case_specs()}
    zero = points['zero']
    admit_zero(zero)
    for section, result in [('pair', pair_readout(reference, points.__getitem__, zero)),
                            ('setting', check_readout(points.__getitem__, zero)),
                            ('cuts', cuts_readout(reference, points.__getitem__, zero))]:
        rebuilt = dict(result, schema=SCHEMA, model=reference['model'], a_A=reference['a_A'],
                       reference=reference_source)
        require(same_readout(rebuilt, read_json(args.native/f'{section}.json')), f'Replayed {section} differs')
    print(f'Replayed {len(points)} native TiD2 records, three readouts, and all relaxation trials')


if __name__ == '__main__':
    main()
