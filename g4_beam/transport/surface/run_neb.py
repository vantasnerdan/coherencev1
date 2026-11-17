"""Measured-allocation NEB entrypoint; no hidden endpoint replays."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
from surface_consumer import bound_json
from surface_io import manifest, save
from surface_neb import run_neb


def main():
    root = Path(__file__).resolve().parent
    allocation = json.loads((root/'neb-decision.json').read_text())
    if allocation['decision'] != 'continue':
        raise ValueError('Explicit measured whole-cost continuation required')
    for name in ('adsorbed-1ML', 'gas-slab-endpoint'):
        result = bound_json(allocation['endpoints'][name])
        if not result['converged']:
            raise ValueError('Unconverged endpoint cannot initialize MEP')
        final = bound_json(result['final'])
        if final['fmax_eV_A'] > result['fmax_target_eV_A']:
            raise ValueError('Endpoint force criterion inconsistent with completion')
    run_neb(root/'adsorbed-1ML', root/'gas-slab-endpoint', root/'neb',
            interior=allocation['interior_images'], steps=allocation['max_steps'],
            fmax=0.05)
    print('NEB_COMPLETED', flush=True)

if __name__ == '__main__':
    main()
