"""Export source-bound local site cuts and matched two-point FD sensitivities."""
import argparse
import json
from pathlib import Path

import numpy as np

from readout_bulk import load, sha


FINAL_AMPLITUDES = (12, 36)
AXES = ('normal', 'in1', 'in2')


def load_record(path):
    row = load(path)
    if row['schema'] != 'g10-site-point-v1' or not row['scf_converged']:
        raise ValueError(f'Not a converged native site point: {path}')
    row['native_point_sha256'] = sha(path)
    return row


def geometry_keys(row):
    return {k: row[k] for k in ('cell_A','pbc','symbols','masses_u','indices')}


def point_difference(base, row, case, moving):
    if (base['reference']['sha256'] != row['reference']['sha256'] or geometry_keys(base) != geometry_keys(row)
            or base['calculator'] != row['calculator'] or row['point']['case'] != case or row['moving_index'] != moving):
        raise ValueError(f'Site point/zero mismatch: {case}')
    xyz = np.asarray(base['positions_A']).copy()
    move = row['point']
    axis, signed = case.split('-')
    expected_step = int(signed[1:]) / 100 * (1 if signed[0] == 'p' else -1)
    if (move['axis'] != axis or move['displacement_A'] != expected_step
            or row['ray_units'] != base['ray_units'] or move['ray_unit'] != base['ray_units'][axis]
            or move['reference_D_position_A'] != base['positions_A'][moving]):
        raise ValueError(f'Site signed displacement, direction or origin differs: {case}')
    ray = np.asarray(move['ray_unit'])
    xyz[moving] += move['displacement_A']*ray
    forces = np.asarray(row['forces_eV_A'])
    if not np.allclose(xyz, row['positions_A'], atol=1e-12, rtol=0) or forces.shape != xyz.shape or not np.isfinite(forces).all():
        raise ValueError(f'Native position or force vectors differ: {case}')
    return dict(row, relative_free_energy_eV=row['energy_free_eV']-base['energy_free_eV'],
                force_response_eV_A=(forces-base['forces_eV_A']).tolist(),
                projected_force_eV_A=float(forces[moving]@ray))


def check_grid(base_zero, base_point, tight_zero, tight_point):
    for base, tight in [(base_zero,tight_zero),(base_point,tight_point)]:
        if (base['resolution']!='base' or tight['resolution']!='tight' or base['positions_A']!=tight['positions_A']
                or base['point']!=tight['point'] or base['reference']['sha256']!=tight['reference']['sha256']
                or tight['benchmark']['base_point']['sha256']!=base['native_point_sha256']
                or not all(base['calculator'][k]==tight['calculator'][k] for k in base['calculator'] if k not in ('gpts','h_A'))
                or tight['calculator']['h_A']>=base['calculator']['h_A']):
            raise ValueError('Grid sensitivity did not hold geometry and all other settings fixed')
    case = base_point['point']['case']
    baseline = point_difference(base_zero,base_point,case,base_point['moving_index'])
    comparison = point_difference(tight_zero,tight_point,case,tight_point['moving_index'])
    df = np.asarray(comparison['force_response_eV_A']) - baseline['force_response_eV_A']
    de = comparison['relative_free_energy_eV'] - baseline['relative_free_energy_eV']
    return dict(delta_relative_free_energy_eV=de, max_force_response_difference_eV_A=float(np.linalg.norm(df,axis=1).max()),
                force_response_difference_eV_A=df.tolist(), base_grid=base_zero['calculator']['gpts'],
                refined_grid=tight_zero['calculator']['gpts'], tested_displacement=case)


def build(folder):
    result = {}
    for site, movers, bench in [('surface', {'adsorbed':28,'subsurface':26},'normal-p12'),
                                ('vacancy', {'vacancy':7},'in1-p12')]:
        root = folder/site
        base_zero = load_record(root/'zero.json')
        groups = {}
        for label,index in movers.items():
            groups[label] = []
            for axis in AXES:
                for amp in FINAL_AMPLITUDES:
                    for sign in 'mp':
                        case = f'{axis}-{sign}{amp:02d}'
                        groups[label].append(point_difference(base_zero,load_record(root/label/f'{case}.json'),case,index))
        base_point = next(row for row in next(iter(groups.values())) if row['point']['case']==bench)
        comparison = check_grid(base_zero,base_point,load_record(root/'tight-zero.json'),load_record(root/f'tight-{bench}.json'))
        result[site] = dict(zero=base_zero,directions=groups,grid_comparison=comparison)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('saved_root',type=Path)
    p.add_argument('output',type=Path)
    a=p.parse_args()
    result = build(a.saved_root)
    with a.output.open('x') as f:
        json.dump(result,f,allow_nan=False,indent=2,sort_keys=True)
        f.write('\n')


if __name__=='__main__':
    main()
