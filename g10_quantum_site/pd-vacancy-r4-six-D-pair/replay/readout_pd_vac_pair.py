"""Read the matched D0–D3 vacancy cut without ASE or GPAW."""
import argparse
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path

import numpy as np

REF_SHA = '0bda3d044fdb291542b0283e299bf1d0e93c3bb25642f0def6f11baf8a9ac888'
OUTPUT_SCHEMA = 'g10-pd-vacancy-r4-d6-pair-readout-v1'
I, J = 63, 66
AMPS = (6, 12, 24)
ROLES = ('A', 'B', 'AB')
RADIUS_SUM = 3.4*0.529177210903


def require(value, message):
    if not value:
        raise ValueError(message)


def load(path):
    raw = path.read_bytes()
    row = json.loads(raw)
    x = np.asarray(row['positions_A'])
    f = np.asarray(row['forces_eV_A'])
    e = float(row['energy_free_eV'])
    require(row['scf_converged'] is True and x.shape == f.shape == (69,3)
            and np.isfinite(x).all() and np.isfinite(f).all() and math.isfinite(e),
            f'Bad SCF or array: {path.name}')
    return row, x, f, e, dict(file=path.name, bytes=len(raw), sha256=sha256(raw).hexdigest())


def native_case(path, reference, mesh, role, amp, axis, images):
    p = load(path)
    row, x = p[:2]
    native, old = reference[:2]
    require(row['schema'] == 'g10-pd-vacancy-r4-d6-pair-v1' and row['reference'] == reference[-1],
            f'Wrong reference or schema: {path.name}')
    require(all(row[k] == native[k] for k in ('indices','symbols','masses_u','cell_A','pbc')),
            f'Wrong atom order, cell or isotope: {path.name}')
    expected = native['specification']['calculator'].copy()
    expected['kmesh'] = [mesh]*3
    expected['symmetry'] = {'point_group': False, 'time_reversal': mesh == 3}
    require(row['calculator'] == expected and len(row['ibz_kpoints']) == (8 if mesh == 2 else 14),
            f'Wrong mesh or numerical settings: {path.name}')
    v = row['point']
    pred = old.copy()
    if 'A' in role:
        pred[I] += (amp/100)*axis
    if 'B' in role:
        pred[J] -= (amp/100)*axis
    require((v['role'], v['displacement_per_D_A'], v['tags']) == (role, amp/100, [I,J])
            and np.allclose(x,pred,atol=1e-12,rtol=0) and np.allclose(v['axis_A_to_B'],axis,atol=1e-12,rtol=0),
            f'Wrong frozen displacement: {path.name}')
    r = np.linalg.norm(x[J]-x[I]+images,axis=1).min()
    pd = np.linalg.norm(x[[I,J],None,None,:]-x[None,:I,None,:]+images,axis=-1).min()
    require(abs(v['separation_A']-r)<1e-12 and abs(v['minimum_tagged_Pd_D_A']-pd)<1e-12
            and v['nominal_Pd_H_overlap'] == (pd<RADIUS_SUM), f'Wrong distances or PAW flag: {path.name}')
    return p


def triplet(zero, points, axis):
    a, b, ab = points
    e0 = zero[3]
    ea, eb, et = (p[3]-e0 for p in points)
    f0, fa, fb, fab = (p[2] for p in (zero,*points))
    path_force = (fab[J]-fab[I]) @ axis
    connected_force = ((fab[J]-fb[J])-(fab[I]-fa[I])) @ axis
    g = ab[0]['point']
    return dict(displacement_per_D_A=g['displacement_per_D_A'], separation_A=g['separation_A'],
                A_eV=ea, B_eV=eb, AB_eV=et, connected_eV=math.fsum((ab[3],e0,-a[3],-b[3])),
                AB_dE_du_eV_A=float(path_force), connected_dE_du_eV_A=float(connected_force),
                mixed_forces_eV_A=(fab-fa-fb+f0).tolist(), minimum_AB_Pd_D_A=g['minimum_tagged_Pd_D_A'],
                any_tagged_Pd_H_overlap=any(p[0]['point']['nominal_Pd_H_overlap'] for p in points),
                sources=[p[-1] for p in (zero,*points)])


def readout(root, reference_path):
    z = load(reference_path)
    ref, x, f, *_ = z
    require(z[-1]['sha256']==REF_SHA and ref['ionic_converged'] and ref['distances']['all_trapped_near_centre']
            and np.linalg.norm(f,axis=1).max()<=0.05,'Wrong exact resting reference')
    images = np.array(list(product(range(-2,3),repeat=3))) @ np.asarray(ref['cell_A'])
    pairs = x[J]-x[I]+images
    norms = np.linalg.norm(pairs,axis=1)
    k = norms.argmin()
    r0 = float(norms[k]); axis = pairs[k]/r0
    zero_pd = float(np.linalg.norm(x[[I,J],None,None,:]-x[None,:I,None,:]+images,axis=-1).min())
    missing = []
    def get(mesh, role, amp):
        name=f'm{mesh}-{role.lower()}-{amp}.json'
        p=root/name
        if not p.is_file():
            missing.append(name)
            return None
        return native_case(p,z,mesh,role,amp,axis,images)
    rows = []
    for u in AMPS:
        points=[get(2,r,u) for r in ROLES]
        if all(p is not None for p in points):
            rows.append(triplet(z,points,axis))
    fine0=get(3,'zero',0)
    fines=[get(3,r,max(AMPS)) for r in ROLES]
    coarse_last=next((p for p in rows if p['displacement_per_D_A']==max(AMPS)/100),None)
    check=None
    if fine0 is not None and coarse_last is not None and all(p is not None for p in fines):
        fine=triplet(fine0,fines,axis)
        dw=fine['connected_eV']-coarse_last['connected_eV']
        den=abs(coarse_last['connected_eV'])
        check=dict(fine_triplet=fine, connected_shift_eV=dw, percent_absolute_shift=100*abs(dw)/den if den else None,
                   ten_percent_target_pass=abs(dw)<=.10*den,
                   total_AB_shift_eV=fine['AB_eV']-coarse_last['AB_eV'],
                   connected_derivative_shift_eV_A=fine['connected_dE_du_eV_A']-coarse_last['connected_dE_du_eV_A'])
    return dict(schema=OUTPUT_SCHEMA, zero=z[-1], original_atom_indices=[I,J], rest_separation_A=r0,
                zero_minimum_tagged_Pd_D_A=zero_pd, zero_has_nominal_PAW_overlap=zero_pd<RADIUS_SUM,
                axis_A_to_B=axis.tolist(), base_kmesh=[2,2,2], rows=rows, shortest_point_mesh_check=check,
                registered_missing=missing)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('native',type=Path)
    p.add_argument('reference',type=Path)
    a=p.parse_args()
    print(json.dumps(readout(a.native,a.reference),indent=2,allow_nan=False))
