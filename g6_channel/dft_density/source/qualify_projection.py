"""Independent continuous-line quadrature and Nyquist-sensitive projection oracles."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from project import project, evaluate


def direct_line_average(density, fractions):
    n = density.shape[0]
    frequencies = np.fft.fftfreq(n) * n
    coefficients = np.fft.fftn(density) / density.size
    times = np.arange(4 * n) / (4 * n)

    def basis(points):
        values = np.exp(2j * np.pi * np.asarray(points)[:, None] * frequencies)
        values[:, n // 2] = np.cos(np.pi * n * np.asarray(points))
        return values

    averages = []
    for x, y in fractions:
        z_basis = basis([y])[0]
        collapsed = np.einsum('ijk,k->ij', coefficients, z_basis)
        values = np.einsum('ti,ij,tj->t', basis(times + x), collapsed, basis(times - x))
        if np.max(np.abs(values.imag)) > 1e-11:
            raise ValueError('Quadrature oracle produced complex real density')
        averages.append(float(values.real.mean()))
    return np.asarray(averages)


def qualify(config_path, readout_path):
    cpu = time.process_time()
    config = json.loads(config_path.read_bytes())
    geometry = json.loads((readout_path / 'geometry.json').read_bytes())
    n, grid = 8, 128
    x, y, z = np.meshgrid(*[np.arange(n) / n] * 3, indexing='ij')
    X, Y = np.meshgrid(*[np.arange(grid) / grid] * 2, indexing='ij')
    controls = [
        ('zero_longitudinal', 1 + .2 * np.cos(2*np.pi*(2*x-2*y+z)), 1 + .2*np.cos(2*np.pi*(4*X+Y))),
        ('nonzero_longitudinal', 1 + .2 * np.cos(2*np.pi*(x+y+z)), np.ones_like(X)),
        ('Nyquist_product_half_weight', 1 + .3*np.cos(2*np.pi*4*x)*np.cos(2*np.pi*4*y), 1 + .15*np.cos(2*np.pi*8*X)),
        ('transverse_Nyquist', 1 + .2*np.cos(2*np.pi*4*z), 1 + .2*np.cos(2*np.pi*4*Y))]
    errors = {}
    for name, density, expected in controls:
        actual, _, _, _ = project(density, grid)
        errors[name] = float(np.max(np.abs(actual - expected)))
    fractions = np.asarray(geometry['centres_xy_A']) / geometry['lengths_A'][:2]
    for case in config['cases']:
        with np.load(config_path.parent / case['path'], allow_pickle=False) as data:
            density = data['density_valence_e_A3']
        with np.load(readout_path / (case['name'] + '-transverse.npz'), allow_pickle=False) as data:
            projected = evaluate(data['coefficients_e_A3'], fractions)
        independent = direct_line_average(density, fractions)
        errors[case['name'] + '_direct_continuous_line_quadrature'] = float(np.max(np.abs(projected - independent)))
    refused = []
    for name, density, output_grid in [('odd_grid', np.zeros((7,7,7)),128),
                                      ('aliased_output',np.zeros((8,8,8)),16),
                                      ('nonfinite',np.full((8,8,8),np.nan),128)]:
        try:
            project(density, output_grid)
        except ValueError:
            refused.append(name)
        else:
            raise ValueError('Wrong input was not refused: ' + name)
    result = dict(schema='pdd-projection-qualification-v1', errors=errors,
                  tolerance=1e-10, wrong_inputs_refused=refused,
                  process_cpu_s=time.process_time()-cpu,
                  passed=bool(max(errors.values()) < 1e-10),
                  Nyquist_scope='Sampled Nyquist aliases are not distinguishable physical fields; the oracle integrates the explicitly declared canonical real interpolant, not an unseen continuous density.')
    if not result['passed']:
        raise ValueError(json.dumps(result))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--readout',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args = parser.parse_args()
    result = qualify(args.config,args.readout)
    with args.report.open('x') as handle:
        handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
