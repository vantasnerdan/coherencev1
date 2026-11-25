"""Admit saved pair-setting and tetrahedral-ray points against their exact geometries."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from pair_geometry import H_AE_SHA, H_BASIS_SHA, displace, paw_overlap, paw_support
from readout_bulk import check_point, load, sha
from readout_pair import native


def sensitivity_point(path, treatment, case, p, bind):
    row = load(path)
    u, role = (.8 if case != 'zero' else 0.), ('A' if case == 'one' else 'AB')
    pos, i, j, shift, axis, length, _ = displace(p, u, role)
    forces = np.asarray(row['forces_eV_A'])
    if treatment == 'pair_mesh':
        mesh, hydrogen = 4, dict(hydrogen='standard_PAW')
    else:
        mesh, hydrogen = 3, dict(hydrogen='projector_free_H', h_implementation_sha256=H_AE_SHA,
                                 h_basis_sha256=H_BASIS_SHA)
    if (row['schema'] != 'g10-pair-sensitivity-v1' or row['case'] != case or row['role'] != role
            or row['treatment'] != treatment or row['cutoff_eV'] != p['cutoff_eV']
            or row['kmesh'] != [mesh]*3 or row['bz_kpoints'] != mesh**3
            or any(row.get(key) != value for key, value in hydrogen.items())
            or row['coarse_zero_provenance'] != bind or row['reference_profile_id'] != p['profile_id']
            or row['target_indices'] != [i, j] or row['image_of_second_lattice'] != shift.tolist()
            or row['displacement_per_site_A'] != u or abs(row['distance_A']-length) > 1e-12
            or not np.allclose(row['closing_axis'], axis, rtol=0, atol=1e-12)
            or not np.allclose(row['positions_A'], pos, rtol=0, atol=1e-12)
            or forces.shape != pos.shape or not np.isfinite(forces).all() or not np.isfinite(row['energy_free_eV'])):
        raise ValueError(f'Pair sensitivity geometry or numerical profile mismatch: {path}')
    row['native_row_sha256'] = sha(path)
    return row


def mesh_contrast(z3, p3, z4, p4, indices, axis):
    e3 = p3['energy_free_eV']-z3['energy_free_eV']
    e4 = p4['energy_free_eV']-z4['energy_free_eV']
    i, j = indices
    fz3, f3, fz4, f4 = (np.asarray(r['forces_eV_A']) for r in (z3, p3, z4, p4))
    err = float(np.max(np.abs(f4-f3)))
    return dict(deltaE_coarse_eV=e3, deltaE_fine_eV=e4, fine_minus_coarse_eV=e4-e3,
                zero_max_force_component_difference_eV_A=float(np.max(np.abs(fz4-fz3))),
                closest_max_force_component_difference_eV_A=err,
                closest_path_derivative_coarse_eV_A=-float((f3[i]-f3[j]) @ axis),
                closest_path_derivative_fine_eV_A=-float((f4[i]-f4[j]) @ axis),
                within_targets=abs(e4-e3)<=.005 and err<=.02,
                scope='Total pair approach E(AB)-E(0), not connected E(AB)-E(A)-E(B)+E(0).')


def hard_contrast(standard, projector_free):
    def energies(group):
        zero, one, closest = (group[key]['energy_free_eV'] for key in ('zero', 'one', 'closest'))
        return closest-zero, math.fsum((closest, -one, -one, zero))
    total_base, connected_base = energies(standard)
    total_h, connected_h = energies(projector_free)
    return dict(total_standard_eV=total_base, total_projector_free_eV=total_h,
                total_projector_free_minus_standard_eV=total_h-total_base,
                connected_standard_eV=connected_base, connected_projector_free_eV=connected_h,
                connected_projector_free_minus_standard_eV=connected_h-connected_base,
                scope='The matched B energy equals A by periodic midpoint inversion. '
                      'Same finite plane-wave cutoff; not an isolated overlap-error bound.')


def build(bulk, pair, fine, hard, refs, partial):
    pp, zp = refs/'Pd_profile.json', refs/'Pd_zero.json'
    p, z = load(pp), load(zp)
    if p != z['profile']:
        raise ValueError('Extension zero/profile mismatch')
    bind = dict(profile_sha256=sha(pp), zero_sha256=sha(zp), checkpoint_sha256=z['checkpoint']['sha256'])
    ray = np.ones(3) / np.sqrt(3.)
    tetra, mesh, alternative, missing = [], {}, {}, []
    for case, u in [('p111p80', .8), ('p111p110', 1.1), ('p111p140', 1.4), ('p111tetra', np.sqrt(3.)*p['a_A']/4)]:
        path = bulk/f'{case}.json'
        if path.is_file():
            point = check_point(path, case, u, ray, p, z, bind)
            point.update(paw_overlap(p, point['positions_A']))
            tetra.append(point)
        else:
            missing.append(case)
    for treatment, folder, cases, rows in (
            ('pair_mesh', fine, ('zero', 'closest'), mesh),
            ('pair_hard', hard, ('zero', 'one', 'closest'), alternative)):
        for case in cases:
            path = folder/f'{case}.json'
            if path.is_file():
                point = sensitivity_point(path, treatment, case, p, bind)
                diagnostic = paw_overlap(p, point['positions_A'])
                if treatment == 'pair_hard':
                    # Nominal geometry compared to the *standard* PAW files only.
                    point['standard_PAW_geometry'] = diagnostic
                else:
                    point.update(diagnostic)
                rows[case] = point
            else:
                missing.append(f'{treatment}-{case}')
    contrast, hard_check = None, None
    if len(mesh) == 2 or len(alternative) == 3:
        ab = native(pair/'AB-80.json', 'AB', .8, p, bind)
        if len(mesh) == 2:
            contrast = mesh_contrast(z, ab, mesh['zero'], mesh['closest'], ab['target_indices'], np.asarray(ab['closing_axis']))
        if len(alternative) == 3:
            one = check_point(bulk/'p110p80.json', 'p110p80', .8, np.asarray(ab['closing_axis']), p, z, bind)
            hard_check = hard_contrast(dict(zero=z, one=one, closest=ab), alternative)
    if missing and not partial:
        raise ValueError(f'Missing requested additions: {missing}')
    return dict(schema='g10-pair-additions-v1', tetrahedral_ray=tetra, fine_native=mesh, projector_free_native=alternative,
                pair_mesh_contrast=contrast, pair_projector_free_contrast=hard_check, missing_native_cases=missing,
                paw_support=paw_support(p))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('bulk_native', 'pair_native', 'mesh_native', 'projector_free_native', 'references', 'output'):
        p.add_argument(key, type=Path)
    p.add_argument('--partial', action='store_true')
    a = p.parse_args()
    result = build(a.bulk_native, a.pair_native, a.mesh_native, a.projector_free_native, a.references, a.partial)
    with a.output.open('x') as f:
        json.dump(result, f, allow_nan=False, indent=2, sort_keys=True)
        f.write('\n')


if __name__ == '__main__':
    main()
