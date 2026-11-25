"""Check 0/A/B/AB inputs and separate force response from the path gradient."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from pair_geometry import AMPS, displace, paw_overlap, paw_support
from readout_bulk import check_point, load, sha


def native(path, role, u, profile, bind):
    row = load(path)
    pos, i, j, image, axis, distance, _ = displace(profile, u, role)
    force = np.asarray(row['forces_eV_A'])
    if (row['role'] != role or row['displacement_per_site_A'] != u or row['reference_profile_id'] != profile['profile_id']
            or row['zero_provenance'] != bind or row['target_indices'] != [i, j] or row['image_of_second_lattice'] != image.tolist()
            or not np.allclose(row['closing_axis'],axis,atol=1e-12,rtol=0) or abs(row['distance_A']-distance)>1e-12
            or not np.allclose(row['positions_A'],pos,atol=1e-12,rtol=0) or force.shape != pos.shape
            or not np.isfinite(force).all() or not np.isfinite(row['energy_free_eV'])):
        raise ValueError(f'Mismatched native pair row {path}')
    row['native_row_sha256'] = sha(path)
    return row


def combine(zero, a, b, ab, i, j, n):
    f0, fa, fb, fab = (np.asarray(r['forces_eV_A']) for r in (zero, a, b, ab))
    # These two residual forces differentiate the *two-coordinate well* on the
    # diagonal; E_A has no derivative wrt moving B, and vice versa.
    ra, rb = fab[i]-fa[i], fab[j]-fb[j]
    return dict(pair_nonadditive_energy_eV=math.fsum([ab['energy_free_eV'], -a['energy_free_eV'], -b['energy_free_eV'], zero['energy_free_eV']]),
                pair_force_on_A_eV_A=ra.tolist(), pair_force_on_B_eV_A=rb.tolist(), d_pair_energy_du_eV_A=-float((ra-rb) @ n),
                all_atom_force_mixed_difference_eV_A=(fab-fa-fb+f0).tolist())


def join(pair, bulk, refs):
    pp, zp = refs/'Pd_profile.json', refs/'Pd_zero.json'
    p, z = load(pp), load(zp)
    if p != z['profile']:
        raise ValueError('Pair zero/profile mismatch')
    bind = dict(profile_sha256=sha(pp), zero_sha256=sha(zp), checkpoint_sha256=z['checkpoint']['sha256'])
    result, b_check = [], None
    for amp in AMPS:
        u = amp/100
        bpos, i, j, _, ray, _, perm = displace(p, u, 'B')
        old = bulk/f'p110p{amp:02d}.json'
        a = check_point(old, old.stem, u, ray, p, z, bind) if old.is_file() else native(pair/f'A-{amp}.json', 'A', u, p, bind)
        pred_b = -np.asarray(a['forces_eV_A'])[perm]
        b_path = pair/f'B-{amp}.json'
        if b_path.is_file():
            b = native(b_path, 'B', u, p, bind)
            if amp == 24:
                de = abs(a['energy_free_eV'] - b['energy_free_eV'])
                df = float(np.max(np.abs(pred_b - np.asarray(b['forces_eV_A']))))
                b_check = dict(delta_energy_eV=de, max_force_component_difference_eV_A=df, passed=de <= .005 and df <= .02)
        else:
            b = dict(role='B', sample_kind='midpoint_inversion_of_A', energy_free_eV=a['energy_free_eV'],
                     positions_A=bpos.tolist(), forces_eV_A=pred_b.tolist(), parent_native_row_sha256=a['native_row_sha256'])
        ab = native(pair/f'AB-{amp}.json', 'AB', u, p, bind)
        for row in (a, b, ab):
            row.update(paw_overlap(p, row['positions_A']))
        result.append(dict(u_A=u, separation_A=ab['distance_A'],
                           widest_tabulated_Pd_H_partial_wave_overlap=ab['widest_tabulated_Pd_H_partial_wave_overlap'],
                           one_A=a, one_B=b, pair_AB=ab, **combine(z,a,b,ab,i,j,ray)))
    if not b_check or not b_check['passed']:
        raise ValueError(f'Unqualified midpoint inversion control: {b_check}')
    return dict(schema='g10-pair-export-v1', target_indices=[i,j], closing_axis=ray.tolist(), zero_geometry=p['geometry'],
                zero_provenance=bind, zero_free_energy_eV=z['energy_free_eV'], midpoint_control=b_check, rows=result,
                paw_support=paw_support(p), zero_paw_geometry=paw_overlap(p, z['positions_A']),
                force_policy='pair_force_on_A/B is the correction to each moving-site one-body force; all_atom_force_mixed_difference is the two-displacement response at fixed original index, not the derivative on the pair path')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for k in ('pair_native','bulk_native','references','output'):
        p.add_argument(k,type=Path)
    a = p.parse_args()
    result = join(a.pair_native,a.bulk_native,a.references)
    with a.output.open('x') as f:
        json.dump(result,f,allow_nan=False,indent=2,sort_keys=True)
        f.write('\n')


if __name__ == '__main__':
    main()
