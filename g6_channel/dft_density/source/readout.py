"""NumPy-only, hash-bound readout of the two admitted Native valence grids."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from geometry import construct, compare_pinned_geometry, outside_mask
from project import project, evaluate


def binding(path):
    data = Path(path).read_bytes()
    return dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def bound_path(record, base):
    path = base / record['path']
    if binding(path) != {key: record[key] for key in ('bytes', 'sha256')}:
        raise ValueError('Changed bound input: ' + Path(record['path']).name)
    return path


def run(config_path, output):
    cpu, wall = time.process_time(), time.perf_counter()
    config = json.loads(config_path.read_bytes())
    base = config_path.parent
    if [case['name'] for case in config['cases']] != ['h020', 'h016']:
        raise ValueError('Require the two admitted h cases in coarse/fine order')
    output.mkdir(parents=True, exist_ok=False)
    results, geometry = {}, None
    for case in config['cases']:
        path = bound_path(case, base)
        with np.load(path, allow_pickle=False) as data:
            density = data['density_valence_e_A3']
            cell, positions, symbols = data['cell_A'], data['positions_A'], data['symbols']
            if not np.array_equal(data['origin_A'], np.zeros(3)):
                raise ValueError('Wrong density grid origin')
            current = construct(cell, positions, symbols, config['pd_augmentation_matching_radius_A'])
            if geometry is None:
                geometry = current
            elif geometry != current:
                raise ValueError('The two pristine cells/strings differ')
            if list(density.shape) != case['grid_shape']:
                raise ValueError('Wrong admitted Native grid shape')
            integral = float(density.sum(dtype=np.float64) * np.linalg.det(cell) / density.size)
            if abs(integral - case['native_integrated_electrons']) > 1e-10:
                raise ValueError('Native integrated count changed')
            surface, coefficients, support, no_nyquist = project(density)
        fractions = np.asarray(geometry['centres_xy_A']) / geometry['lengths_A'][:2]
        centres = evaluate(coefficients, fractions)
        no_nyquist_centres = evaluate(no_nyquist, fractions)
        mask = outside_mask(geometry, surface.shape[0])
        if not mask.any():
            raise ValueError('No transverse samples outside Pd augmentation matching cylinders')
        npz = output / (case['name'] + '-transverse.npz')
        np.savez_compressed(npz, n_perp_e_A3=surface, coefficients_e_A3=coefficients,
                            supported_modes=support, outside_Pd_augmentation=mask,
                            centres_xy_A=geometry['centres_xy_A'], centres_n_e_A3=centres,
                            transverse_lengths_A=geometry['lengths_A'][:2])
        results[case['name']] = dict(native_density_binding={k: case[k] for k in ('bytes', 'sha256')},
            grid_shape=case['grid_shape'], integrated_electrons=integral, expected_electrons=44,
            integral_difference_e=integral - 44, mean_n_e_A3=float(surface.mean()),
            projected_orthocell_integrated_electrons=float(surface.mean() * np.prod(geometry['lengths_A'])),
            centres_n_e_A3=centres.tolist(), centres_n_without_Nyquist_e_A3=no_nyquist_centres.tolist(),
            centre_Nyquist_sensitivity_e_A3=(centres - no_nyquist_centres).tolist(),
            outside_Pd_augmentation_samples=int(mask.sum()), plotting_grid=list(surface.shape),
            transverse_npz=binding(npz))
    pinned = bound_path(config['pinned_channel_source'], base).read_bytes()
    differential = compare_pinned_geometry(geometry, pinned)
    coarse = np.asarray(results['h020']['centres_n_e_A3'])
    fine = np.asarray(results['h016']['centres_n_e_A3'])
    if np.any(fine <= 0):
        raise ValueError('Nonpositive channel-centre density')
    relative = np.abs(coarse - fine) / np.abs(fine)
    result = dict(schema='pdd-channel-density-readout-v1', cases=results,
        geometry_differential=differential, centre_relative_h_change=relative.tolist(),
        max_centre_relative_h_change=float(relative.max()), centre_tolerance=0.05,
        centre_h_passed=bool(relative.max() < 0.05),
        representation='Native skip_core as returned; no correction or renormalization. Canonical real Hermitian Nyquist-split trigonometric interpolant averaged continuously along110.',
        projected_cell='Orthogonal (1-10,001,110), volume twice the original cubic cell; density units remain electrons/A^3, not electrons/A^2.',
        process_cpu_s=time.process_time() - cpu, wall_s=time.perf_counter() - wall)
    (output / 'geometry.json').write_text(json.dumps(geometry, indent=2) + '\n')
    (output / 'readout.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.config, args.output)
