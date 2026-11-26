"""Join exact portable D26/D28 slab native records; print the response as JSON."""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path

import numpy as np

REFERENCE = '7aadd280b7a58106a79702832d27c33639d04a2a5e14899065abe6b879a9c6a9'
I, J = 26, 28


def require(ok, message):
    if not ok:
        raise ValueError(message)


def point(root, mesh, role, amp):
    path = root / f'm{mesh}-{role.lower()}-{amp}.json'
    raw = path.read_bytes()
    row = json.loads(raw)
    require(row['schema'] == 'g10-covered-slab-pair-v1' and row['scf_converged'] is True
            and row['reference']['sha256'] == REFERENCE and row['native_grid'] == [24, 24, 136]
            and row['calculator']['kmesh'] == [mesh, mesh, 1]
            and row['electron_count'] == 173, f'Wrong native provenance: {path.name}')
    g = row['point']
    expected = (role, amp/100, [I, J])
    require((g['role'], g['displacement_per_D_A'], g['original_D_indices']) == expected,
            f'Wrong native pair location: {path.name}')
    f = np.asarray(row['forces_eV_A'])
    xyz = np.asarray(row['positions_A'])
    e = float(row['energy_free_eV'])
    require(xyz.shape == f.shape == (29, 3) and np.isfinite(f).all() and np.isfinite(xyz).all()
            and np.isfinite(e), f'Nonfinite or incomplete native array: {path.name}')
    return row, xyz, e, f, dict(file=path.name, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def validate(point, zero, role, amp, n, translation):
    row, xyz = point[:2]
    ref, original = zero[:2]
    require(all(row[k] == ref[k] for k in ('indices', 'symbols', 'masses_u', 'cell_A', 'pbc',
                                           'reference', 'setup_sha256')),
            f'Atom order or setup differs: {point[-1]["file"]}')
    native_settings = row['calculator'].copy()
    native_settings['kmesh'] = ref['calculator']['kmesh']
    require(native_settings == ref['calculator'], 'The mesh comparison also changed other settings')
    moved = original.copy()
    if role in ('A', 'AB'):
        moved[I] += amp/100*n
    if role in ('B', 'AB'):
        moved[J] -= amp/100*n
    stated = row['point']
    actual_r = np.linalg.norm(xyz[J] + translation @ np.asarray(ref['cell_A']) - xyz[I])
    require(np.allclose(moved, xyz, atol=2e-12, rtol=0) and
            np.allclose(stated['axis_A_to_B'], n, atol=2e-12, rtol=0) and
            np.array_equal(stated['B_image_cell_translation'], translation) and
            abs(stated['separation_A']-actual_r) < 2e-12,
            f'Geometry does not match the registered frozen displacement: {point[-1]["file"]}')


def axis(zero):
    row, xyz = zero[:2]
    require(row['symbols'][I] == row['symbols'][J] == 'H' and row['pbc'] == [True, True, False],
            'Unexpected reference atoms or periodicity')
    images = np.asarray(list(product((-1, 0, 1), (-1, 0, 1), (0,))))
    vecs = xyz[J] - xyz[I] + images @ np.asarray(row['cell_A'])
    norms = np.linalg.norm(vecs, axis=1)
    winner = int(norms.argmin())
    r0 = norms[winner]
    require(abs(r0-2.437039831543227) < 2e-12, 'Wrong adopted slab reference')
    n, translation = vecs[winner]/r0, images[winner]
    return n, translation, r0


def curve(root, mesh, amps, reference):
    zero = point(root, mesh, 'zero', 0)
    n, translation, _ = axis(reference)
    validate(zero, reference, 'zero', 0, n, translation)
    e0, f0 = zero[2:4]
    curves, missing = [], []
    for amp in amps:
        paths = [root / f'm{mesh}-{role.lower()}-{amp}.json' for role in ('A', 'B', 'AB')]
        if not all(p.is_file() for p in paths):
            missing.extend(p.name for p in paths if not p.is_file())
            continue
        a, b, ab = [point(root, mesh, role, amp) for role in ('A', 'B', 'AB')]
        for p, role in zip((a, b, ab), ('A', 'B', 'AB'), strict=True):
            validate(p, zero, role, amp, n, translation)
        ea, eb, eab = a[2]-e0, b[2]-e0, ab[2]-e0
        fa, fb, fab = a[3], b[3], ab[3]
        tagged = ab[0]['point']
        curves.append(dict(displacement_per_D_A=amp/100, separation_A=tagged['separation_A'],
                           A_eV=ea, B_eV=eb, AB_eV=eab, connected_eV=eab-ea-eb,
                           AB_dE_du_eV_A=float((fab[J]-fab[I]) @ n),
                           connected_dE_du_eV_A=float(((fab[J]-fb[J])-(fab[I]-fa[I])) @ n),
                           mixed_forces_eV_A=(fab-fa-fb+f0).tolist(),
                           any_tagged_Pd_H_overlap=any(p[0]['point']['nominal_Pd_H_overlap'] for p in (a,b,ab)),
                           AB_min_tagged_Pd_D_A=tagged['min_tagged_Pd_D_A'], sources=[p[-1] for p in (a,b,ab)]))
    return zero, curves, missing


def readout(root):
    zero = point(root, 2, 'zero', 0)
    n, translation, r0 = axis(zero)
    _, curves, missing = curve(root, 2, (6, 12, 24, 36, 48), zero)
    row, _, e0, f0, _ = zero
    mesh = None
    mesh_files = ('m2-ab-36.json', 'm3-zero-0.json', 'm3-ab-36.json')
    if all((root/name).is_file() for name in mesh_files):
        ab = point(root, 2, 'AB', 36)
        fine0 = point(root, 3, 'zero', 0)
        fine = point(root, 3, 'AB', 36)
        for p, role, amp in ((ab,'AB',36), (fine0,'zero',0), (fine,'AB',36)):
            validate(p, zero, role, amp, n, translation)
        de = (fine[2]-fine0[2])-(ab[2]-e0)
        response = (fine[3]-fine0[3])-(ab[3]-f0)
        path_delta = ((fine[3][J]-fine[3][I])-(ab[3][J]-ab[3][I])) @ n
        maxf = np.abs(response).max()
        mesh = dict(displacement_per_D_A=0.36, relative_energy_shift_eV=de,
                    approach_derivative_shift_eV_A=float(path_delta),
                    max_force_response_component_shift_eV_A=float(maxf),
                    energy_5_meV_pass=bool(abs(de)<=0.005), force_0p02_eV_A_pass=bool(maxf<=0.02),
                    sources=[p[-1] for p in (ab,fine0,fine)])
    else:
        missing.extend(name for name in mesh_files if not (root/name).is_file() and name not in missing)
    return dict(schema='g10-covered-slab-pair-readout-v1', zero=zero[-1], axis_A_to_B=n.tolist(),
                sources_reference=row['reference'], rest_separation_A=float(r0),
                triplets=curves, registered_missing=missing, two_point_mesh=mesh)

def refined_readout(root):
    coarse = readout(root)
    reference = point(root, 2, 'zero', 0)
    _, _, r0 = axis(reference)
    zero, rows, missing = curve(root, 3, (6, 12, 24, 36, 48), reference)
    m4_point_names = [f'm4-{role}-36.json' for role in ('a', 'b', 'ab')] + ['m4-zero-0.json']
    m4_missing = [name for name in m4_point_names if not (root/name).is_file()]
    switch = next((row for row in rows if row['displacement_per_D_A'] == 0.36), None)
    control = None
    if not m4_missing:
        fine0, (finer,), _ = curve(root, 4, (36,), reference)
        if switch is not None:
            signed = finer['connected_eV']-switch['connected_eV']
            basis = abs(switch['connected_eV'])
            control = dict(separation_A=finer['separation_A'], mesh_4_zero=fine0[-1], mesh_4_triplet=finer,
                           mesh_3_connected_eV=switch['connected_eV'], connected_shift_eV=signed,
                           abs_percent_shift=100*abs(signed)/basis if basis else None,
                           ten_percent_target_pass=abs(signed)<=0.10*basis,
                           AB_energy_rise_shift_eV=finer['AB_eV']-switch['AB_eV'],
                           connected_derivative_shift_eV_A=(finer['connected_dE_du_eV_A']
                                                           -switch['connected_dE_du_eV_A']))
    return dict(schema='g10-covered-slab-pair-mesh-readout-v2', rest_separation_A=float(r0),
                consumer_mesh=[3, 3, 1], mesh_3_zero=zero[-1], triplets=rows,
                registered_missing=missing+m4_missing, connected_mesh_4_control=control, coarse=coarse)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('native', type=Path)
    parser.add_argument('--refined', action='store_true')
    args = parser.parse_args()
    read = refined_readout if args.refined else readout
    print(json.dumps(read(args.native), indent=2, allow_nan=False))
