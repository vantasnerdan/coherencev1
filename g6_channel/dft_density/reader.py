"""Verify and replay this complete density capsule using NumPy only."""
import sys
sys.dont_write_bytecode = True
import contextlib
import io
import json
from pathlib import Path
import tempfile
import numpy as np


def replay(root):
    root = root.resolve()
    sys.path.insert(0, str(root / 'source'))
    from readout import binding, run
    from fit_widths import model_coefficients, residuals
    manifest = json.loads((root / 'manifest.json').read_bytes())
    if manifest['schema'] != 'pdd-channel-density-capsule-v1':
        raise ValueError('Wrong channel-density capsule schema')
    entries = manifest['files']
    actual = {str(path.relative_to(root)) for path in root.rglob('*') if path.is_file()}
    if actual != set(entries) | {'manifest.json'}:
        raise ValueError('Missing or unmanifested capsule file')
    for name, entry in entries.items():
        path = root / name
        if not path.resolve().is_relative_to(root) or path.is_symlink():
            raise ValueError('Unsafe capsule member')
        if binding(path) != entry['EXPORTED']:
            raise ValueError('EXPORTED byte integrity failed: ' + name)
    expected = json.loads((root / 'results/readout.json').read_bytes())
    geometry = json.loads((root / 'results/geometry.json').read_bytes())
    fits = json.loads((root / 'results/width-fit.json').read_bytes())
    errors = []
    with tempfile.TemporaryDirectory(prefix='pdd-density-replay-') as temporary:
        with contextlib.redirect_stdout(io.StringIO()):
            computed = run(root / 'inputs.json', Path(temporary) / 'readout')
        for name in ('h020','h016'):
            errors.append(abs(computed['cases'][name]['integrated_electrons']-expected['cases'][name]['integrated_electrons']))
            errors.append(float(np.max(np.abs(np.asarray(computed['cases'][name]['centres_n_e_A3'])-expected['cases'][name]['centres_n_e_A3']))))
            with np.load(root / 'results' / (name+'-transverse.npz'),allow_pickle=False) as data:
                with np.load(Path(temporary)/'readout'/(name+'-transverse.npz'),allow_pickle=False) as regenerated:
                    errors.append(float(np.max(np.abs(data['n_perp_e_A3']-regenerated['n_perp_e_A3']))))
                coefficients=model_coefficients(fits['cases'][name]['widths_A'],geometry,data['supported_modes'])
                metrics=residuals(coefficients,data['n_perp_e_A3'],data['outside_Pd_augmentation'],geometry)
            for key in ('exterior_RMS_e_A3','exterior_relative_RMS','full_plane_RMS_e_A3','full_plane_relative_RMS'):
                errors.append(abs(metrics[key]-fits['cases'][name]['residuals'][key]))
            errors.append(float(np.max(np.abs(np.asarray(metrics['model_centres_n_e_A3'])-fits['cases'][name]['residuals']['model_centres_n_e_A3']))))
    if max(errors)>1e-10 or computed['centre_h_passed'] != expected['centre_h_passed']:
        raise ValueError('Saved density/projection/fit readout did not replay')
    return dict(schema='pdd-channel-density-replay-v1',passed=True,max_absolute_replay_error=max(errors),
                manifest_files=len(entries),centre_h_passed=computed['centre_h_passed'],
                max_centre_relative_h_change=computed['max_centre_relative_h_change'],
                dependencies='NumPy only; fit optimization is not repeated, saved widths are forward-replayed.')


if __name__=='__main__':
    if len(sys.argv)!=2:
        raise SystemExit('Usage: python reader.py CAPSULE_DIRECTORY')
    print(json.dumps(replay(Path(sys.argv[1])),indent=2),flush=True)
