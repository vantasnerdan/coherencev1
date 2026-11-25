"""Read the exact saved inward selections without rerunning trajectories."""
from pathlib import Path
import json

import numpy as np

import nonlinear as N
from inward import ENERGY, DT, HORIZON, cap_diagnostics, restore
from large_geometry import LargeLattice

HERE = Path(__file__).resolve().parent


def consume():
    rows = []
    for size, m in ((131072, 32), (16384, 16)):
        report = json.loads((HERE/('inward_N'+str(size)+'.json')).read_text())
        assert (report['energy_per_D_eV'], report['dt_fs'], report['horizon_fs']) == (ENERGY, DT, HORIZON)
        lattice = LargeLattice(m)
        with np.load(HERE/('inward_N'+str(size)+'.npz')) as arrays:
            for name in ('original', 'conservative', 'passive_loss', 'loss'):
                final = report[name]['final'] if name in ('conservative', 'loss') else report[name]
                result = restore(arrays, name, final)
                actual = N.summarize(result, 'inward', DT)
                assert actual == final, (size, name)
                launch = arrays[name+'_launch_v']
                check = cap_diagnostics(lattice, launch)
                assert check['valid'], (size, name)
                assert final['selected_event']['signed_inward_A_fs'] > 0
                if name in ('conservative', 'loss'):
                    search = report[name]
                    assert final['score'] >= search['initial']['score']
                    scanned = [t['score'] for t in search['feasible_trials']
                               if t['full_horizon_scanned'] or t['first_original_exit_fs'] is not None]
                    assert not scanned or final['score'] >= max(scanned)*(1-2e-14)
                    assert search['stop_reason'] != 'search_in_progress'
                    assert check == search['launch_constraints'], (size, name, 'saved cap/energy record')
                    if name == 'loss':
                        changed = not np.array_equal(launch, arrays['passive_loss_launch_v'])
                        if size == 16384:
                            assert changed and final['score'] > report['passive_loss']['score']
                        else:
                            assert not changed and final == report['passive_loss']
                rows.append(dict(n_D=size, branch=name,
                    inward_eV=final['score'], event=final['selected_event'],
                    first_exit_fs=final['first_original_exit_fs'],
                    full_horizon_supported=final['full_horizon_supported'], cap=check))
            if size == 131072:
                gradient = arrays['conservative_gradient']
                assert gradient.shape == arrays['original_launch_v'].shape
                assert np.isfinite(gradient).all() and np.linalg.norm(gradient) > 0
                assert report['conservative']['endpoint_support_witness']['passed']
                if len(report['conservative']['accepted_history']) > 1:
                    step = report['conservative']['costs'][0]['step']
                    c = arrays['conservative_trace'][step, 9]
                    assert .25*N.E.MASS*max(0., c)**2 > report['conservative']['initial']['score']
        if size == 131072:
            assert report['conservative']['directional_witness']['passed']
    return dict(schema='pdd-inward-saved-consumer-v1', rows=rows,
        interpretation='finite frozen-model inward selection; no physical quantum probability or material ceiling')


if __name__ == '__main__':
    print(json.dumps(consume(), indent=2, allow_nan=False))
