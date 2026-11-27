"""Replay key delivered quantities using independent decimal finite differences."""
import argparse
import csv
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

from ase.units import GPa

D = Decimal


def near(found, expected, label, absolute='2e-12'):
    found = D(found)
    expected = D(expected)
    tolerance = max(D(absolute), abs(expected) * D('2e-12'))
    if not (found.is_finite() and expected.is_finite() and abs(found - expected) <= tolerance):
        raise ValueError(f'{label}: found {found}, expected {expected}, tolerance {tolerance}')


def point(root, mesh, t, direction, disp):
    path = root / f'k{mesh}' / str(t) / direction / f'u{disp:+.4f}.json'
    return json.loads(path.read_bytes(), parse_float=D)


def raw(root, mesh, t):
    z = point(root, mesh, t, 'zero', D(0))
    d = D('0.0100')
    lp, rp = (point(root, mesh, t, '100-fd', s * d) for s in (-1, 1))
    f, e, s = [z[k] for k in ('energy_free_eV', 'energy_extrapolated_eV', 'stress_eV_A3')]
    kf = (lp['forces_eV_A'][1][0] - rp['forces_eV_A'][1][0]) / (2*d)
    ke = (lp['energy_free_eV'] + rp['energy_free_eV'] - 2*f) / d**2
    return dict(u=2*e-f, f=f, k=kf, ke=ke, vol=z['volume_A3'],
                force=z['forces_eV_A'][1], p=-sum(s[i][i] for i in range(3))/3/D(str(GPa)))


def check_quartics(root, t, group, csv_row, harmonic):
    for axis, data in group['quartic'].items():
        direction = [D(1), D(0), D(0)] if axis == '100' else [D(1)/D(3).sqrt()] * 3
        forces = []
        energies = []
        steps = (D('.12'), D('.24'))
        for s in steps:
            m, p = (point(root, 12, t, f'{axis}-quartic', j*s) for j in (-1, 1))
            forces.append(sum(a*(b-c) for a,b,c in zip(direction, m['forces_eV_A'][1], p['forces_eV_A'][1]))/2)
            energies.append((m['energy_free_eV'] + p['energy_free_eV'])/2 - harmonic['f'])
        k = harmonic['k']
        k4 = sum((f-k*s)*s**3 for f,s in zip(forces,steps))/sum(s**6 for s in steps)
        near(csv_row[f'K4_{axis}_eV_A4'], k4, 'quartic force coefficient')
        for index, (s, f, e) in enumerate(zip(steps, forces, energies)):
            near(data['fit']['force_residual_eV_A'][index], f-k*s-k4*s**3, 'quartic force residual')
            near(data['fit']['free_energy_holdout_residual_eV'][index], e-k*s**2/2-k4*s**4/4, 'quartic energy holdout')


def run(root, readout_path, csv_path):
    data = readout_path.read_bytes()
    readout = json.loads(data, parse_float=D)
    groups = {(r['temperature_K'], k): r[f'mesh{k}'] for r in readout['rows'] for k in (12, 16)}
    for rel, expected in readout['native_sources'].items():
        b = (root / rel).read_bytes()
        if (hashlib.sha256(b).hexdigest(), len(b)) != (expected['sha256'], expected['bytes']):
            raise ValueError(f'Native source bytes differ: {rel}')
    seen = set()
    reference = {k: raw(root, k, 300) for k in (12,16)}
    digest = hashlib.sha256(data).hexdigest()
    with csv_path.open(newline='') as f, localcontext() as ctx:
        ctx.prec = 34
        for row in csv.DictReader(f):
            t, k = int(row['temperature_K']), int(row['kmesh_axis_count'])
            if row['readout_sha256'] != digest or (t,k) in seen or (t,k) not in groups:
                raise ValueError('Mismatched or repeated response-table identity')
            seen.add((t,k))
            v = raw(root,k,t)
            b = reference[k]
            du = v['u'] - b['u']
            dk = v['k'] - b['k']
            near(row['volume_A3_per_D_site'], v['vol'], 'cell volume')
            near(row['electronic_excitation_eV_per_D_site'], du, 'electronic ΔU')
            near(row['electronic_free_energy_change_eV_per_D_site'], v['f']-b['f'], 'electronic ΔF')
            near(row['K_Gamma_100_eV_A2'], v['k'], 'Γ curvature')
            near(row['K_change_eV_A2'], dk, 'Γ contrast')
            near(row['fractional_K_change'], dk/b['k'], 'fractional Γ contrast')
            near(row['pressure_GPa'], v['p'], 'electronic pressure')
            near(row['pressure_change_GPa'], v['p']-b['p'], 'pressure contrast')
            near(row['D_centered_force_norm_eV_A'], sum(f*f for f in v['force']).sqrt(), 'centered force')
            near(groups[t,k]['harmonic']['K_free_energy_eV_A2'], v['ke'], 'free-energy stencil', absolute='3e-11')
            check_quartics(root, t, groups[t,k], row, v)
            if du > 0:
                near(row['fractional_K_secant_per_eV'], dk/b['k']/du, 'per-eV Γ secant', absolute='1e-6')
            elif row['fractional_K_secant_per_eV']:
                raise ValueError('A non-positive excitation is not a reported secant')
    if seen != set(groups):
        raise ValueError('The exported response table omitted native readout rows')
    return len(seen)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('native_root', type=Path)
    p.add_argument('readout', type=Path)
    p.add_argument('table', type=Path)
    a = p.parse_args()
    print(f'Independently checked {run(a.native_root, a.readout, a.table)} native-derived rows')


if __name__ == '__main__':
    main()
