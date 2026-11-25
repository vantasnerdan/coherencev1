"""Verify saved bulk records, then export measured and explicitly derived signed rows."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

LAYOUT = {'p100': (48, 80), 'p110': (12, 24, 48, 80), 'p111': (48, 80), 'n111': (24, 48, 80)}
DIRECT_MINUS = {('p110', 12), ('p110', 24), ('p110', 48)}


def load(path):
    return json.loads(path.read_bytes())


def sha(path):
    with path.open('rb') as fd:
        return hashlib.file_digest(fd, 'sha256').hexdigest()


def inversion(g, center):
    coordinates = np.asarray(g['reference_positions_A'])
    cell, numbers = np.asarray(g['cell_A']), np.asarray(g['atomic_numbers'])
    frac = (coordinates - center) @ np.linalg.inv(cell)
    image = frac[:, None, :] + frac[None, :, :]
    good = (np.linalg.norm((image - np.rint(image)) @ cell, axis=2) < 1e-9) & (numbers[:, None] == numbers)
    if not np.all(good.sum(axis=1) == 1):
        raise ValueError('No unique species-preserving inversion')
    perm = good.argmax(axis=1)
    if not np.array_equal(perm[perm], np.arange(len(perm))):
        raise ValueError('Inversion is not an involution')
    return perm


def check_point(path, case, u, axis, profile, zero, bindings):
    row = load(path)
    g = profile['geometry']
    i = g['moving_index']
    r0 = np.asarray(g['reference_positions_A'])
    expected = r0.copy()
    expected[i] += u * axis
    force = np.asarray(row['forces_eV_A'])
    if (row['case'] != case or row['material'] != profile['material'] or row['reference_profile_id'] != profile['profile_id']
            or row['zero_provenance'] != bindings or row['moving_index'] != i or row['displacement_A'] != u
            or not np.allclose(row['ray_unit'], axis, atol=1e-13, rtol=0)
            or not np.allclose(row['positions_A'], expected, atol=1e-12, rtol=0) or force.shape != r0.shape
            or not np.isfinite(force).all() or not np.isfinite(row['energy_free_eV'])):
        raise ValueError(f'Native geometry/reference/force mismatch: {case}')
    if (abs(force[i] @ axis - row['native_force_parallel_eV_A']) > 1e-10
            or abs(row['energy_free_eV'] - zero['energy_free_eV'] - row['delta_energy_free_eV']) > 1e-10):
        raise ValueError(f'Native energy or force projection mismatch: {case}')
    if np.linalg.norm(force.sum(axis=0)) > 0.005:
        raise ValueError(f'Raw net force diagnostic exceeded target: {case}')
    row['sample_kind'] = 'native_scf'
    row['native_row_sha256'] = sha(path)
    return row


def derive(row, profile, perm):
    g = profile['geometry']
    r = np.asarray(g['reference_positions_A']).copy()
    i, u = g['moving_index'], -row['displacement_A']
    ray = np.asarray(row['ray_unit'])
    r[i] += u * ray
    force = -np.asarray(row['forces_eV_A'])[perm]
    derived = dict(row, case=row['case'][:4]+'m'+row['case'][5:], displacement_A=u,
                   positions_A=r.tolist(), forces_eV_A=force.tolist(), net_force_eV_A=force.sum(axis=0).tolist(),
                   native_force_parallel_eV_A=float(force[i] @ ray), sample_kind='symmetry_derived',
                   parent_native_case=row['case'])
    return derived


def build(native, refs, partial=False):
    sources = {}
    for letter, material in [('p', 'Pd'), ('n', 'Ni')]:
        pp, zp = refs / f'{material}_profile.json', refs / f'{material}_zero.json'
        p, z = load(pp), load(zp)
        if p != z['profile']:
            raise ValueError(f'Reference differs for {material}')
        binding = dict(profile_sha256=sha(pp), zero_sha256=sha(zp), checkpoint_sha256=z['checkpoint']['sha256'])
        g = p['geometry']
        perm = inversion(g, np.asarray(g['reference_positions_A'][g['moving_index']]))
        sources[letter] = (p, z, binding, perm)
    out, controls, missing = [], [], []
    paired = []
    for label, amps in LAYOUT.items():
        p, z, binding, perm = sources[label[0]]
        axis = np.array(tuple(map(float, label[1:])))
        axis /= np.linalg.norm(axis)
        for a in amps:
            key = (label, a)
            rows = {}
            for sign, scale in [('p', 1), ('m', -1)]:
                case = f'{label}{sign}{a:02d}'
                path = native / f'{case}.json'
                if path.is_file():
                    rows[sign] = check_point(path, case, scale*a/100, axis, p, z, binding)
                    out.append(rows[sign])
                elif sign == 'p' or key in DIRECT_MINUS:
                    missing.append(case)
            if 'p' not in rows:
                continue
            if 'm' in rows:
                plus, minus = rows['p'], rows['m']
                expected = derive(plus, p, perm)
                err_f = float(np.max(np.abs(np.asarray(expected['forces_eV_A']) - minus['forces_eV_A'])))
                err_e = abs(expected['energy_free_eV'] - minus['energy_free_eV'])
                controls.append(dict(case=plus['case'], difference_eV=err_e, force_component_max_difference_eV_A=err_f,
                                     passed=err_e <= .005 and err_f <= .02))
            elif key not in DIRECT_MINUS:
                paired.append(derive(rows['p'], p, perm))
    wide = next((r for r in controls if r['case']=='p110p48'), None)
    if not partial and (missing or not wide or not all(r['passed'] for r in controls)):
        raise ValueError(f'Not a verified complete wide-well set: missing={missing}, controls={controls}')
    if wide and wide['passed'] and all(r['passed'] for r in controls):
        out.extend(paired)
    return dict(schema='g10-bulk-wells-export-v1', rows=out, inversion_controls=controls,
                missing_native_cases=missing, derived_admitted=len(out)-sum(r['sample_kind']=='native_scf' for r in out),
                geometry={p['material']:p['geometry'] for p, _, _, _ in sources.values()})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('native', type=Path)
    p.add_argument('references', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--partial', action='store_true')
    a = p.parse_args()
    result = build(a.native, a.references, a.partial)
    with a.output.open('x') as fd:
        json.dump(result, fd, allow_nan=False, sort_keys=True, indent=2)
        fd.write('\n')
    print(f'Export: {len(result["rows"])} rows; {len(result["missing_native_cases"])} missing native; checks {result["inversion_controls"]}')


if __name__ == '__main__':
    main()
