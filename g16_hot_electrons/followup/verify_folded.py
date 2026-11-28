"""Reconstruct the folded thermal contrasts from immutable native bytes."""
import argparse
from decimal import Decimal as D, localcontext
import hashlib
import json
from pathlib import Path

from verify_table import near


MODES = ('gamma_100', 'x_100', 'x_010')
EXPECTED = {(t, k) for t in (300, 3000, 10000) for k in (12, 16)}


def get_points(root, identities):
    points = {}
    for identity in identities:
        path = Path(identity['path'])
        if path.is_absolute() or '..' in path.parts or path.as_posix() in points:
            raise ValueError('Unsafe or duplicate native member')
        data = (root / path).read_bytes()
        if (len(data), hashlib.sha256(data).hexdigest()) != (identity['bytes'], identity['sha256']):
            raise ValueError(f'Folded native source differs: {path}')
        points[path.as_posix()] = json.loads(data, parse_float=D)
    return points


def stiffness(points, t, k, name, zero):
    folder = name.replace('_', '-')
    m, p = [points[f'k{k}/{t}/{folder}/u{sign}0.0100.json'] for sign in '-+']
    axis = 1 if name == 'x_010' else 0
    parity = -1 if name.startswith('x_') else 1
    force_d0 = m['forces_eV_A'][2][axis] - p['forces_eV_A'][2][axis]
    force_d1 = m['forces_eV_A'][3][axis] - p['forces_eV_A'][3][axis]
    force = (force_d0 + parity * force_d1) * 25
    energy = (m['energy_free_eV'] + p['energy_free_eV'] - 2 * zero['energy_free_eV']) * 5000
    return force, energy


def run(root, response_path, changes_path):
    data = response_path.read_bytes()
    response = json.loads(data, parse_float=D)
    changes = json.loads(changes_path.read_bytes(), parse_float=D)
    if changes['folded_response_sha256'] != hashlib.sha256(data).hexdigest():
        raise ValueError('Thermal contrast does not bind to this response')
    points = get_points(root, response['sources'])
    if len(points) != 42:
        raise ValueError('Incomplete folded native input set')
    results = {}
    with localcontext() as ctx:
        ctx.prec = 34
        for row in response['rows']:
            t, k = row['temperature_K'], row['primitive_mesh_n']
            if (t,k) in results:
                raise ValueError('Repeated temperature and mesh')
            zero = points[f'k{k}/{t}/zero/u+0.0000.json']
            u = (2 * zero['energy_extrapolated_eV'] - zero['energy_free_eV'])/2
            near(row['U_eV_per_D_site'], u, 'per-D-electronic-internal-energy')
            modes = {name: stiffness(points,t,k,name,zero) for name in MODES}
            for name, (f,e) in modes.items():
                near(row[name]['K_force_eV_A2'], f, name + ' force curvature')
                near(row[name]['K_energy_eV_A2'], e, name + ' F curvature', absolute='3e-11')
            results[t,k] = u, modes
        if set(results) != EXPECTED:
            raise ValueError('Incomplete folded response rows')
        checked = set()
        for row in changes['rows']:
            t = row['temperature_K']
            for k in (12, 16):
                if (t,k) in checked:
                    raise ValueError('Repeated folded contrast row')
                checked.add((t,k))
                u, modes = results[t,k]
                cold_u, cold = results[300,k]
                out = row[f'mesh{k}']
                near(out['electronic_excitation_eV_per_D_site'], u-cold_u, 'internal excitation per D')
                d = {name: modes[name][0] - cold[name][0] for name in MODES}
                for name in MODES:
                    near(out['mode_K_change_eV_A2'][name], d[name], name + ' thermal change')
                if 'split_change' in out:
                    g = d['gamma_100']
                    for pol in ('100','010'):
                        x = d[f'x_{pol}']
                        split = out['split_change'][pol]
                        near(split['K_diagonal_folded_eV_A2'], (x+g)/2, 'thermal folded diagonal')
                        near(split['J_cross_folded_eV_A2'], (x-g)/2, 'thermal folded cross J')
                        near(split['J_shell1_plus_shell3_scalar_projection_eV_A2'], (x-g)*D('.75'), 'scalar odd-shell change')
        if checked != EXPECTED:
            raise ValueError('Incomplete folded contrast rows')
    return len(points), len(checked)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('native_root', type=Path)
    p.add_argument('response', type=Path)
    p.add_argument('thermal_changes', type=Path)
    a = p.parse_args()
    points, contrasts = run(a.native_root, a.response, a.thermal_changes)
    print(f'Independently verified {points} folded native points and {contrasts} matched thermal contrasts')


if __name__ == '__main__':
    main()
