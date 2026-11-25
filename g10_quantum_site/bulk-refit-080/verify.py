"""One narrow integrated readback; no spectral or sampling replay."""
import csv
import hashlib
import json
import platform
import numpy as np
import scipy
import well

HERE,OUT=well.HERE,well.HERE/'results'

def main():
    source_sha=hashlib.sha256((well.INPUT/'bulk-wells.json').read_bytes()).hexdigest()
    assert source_sha=='18fa6232525243391b84f2fb495cfe9947f7d41eef270931fb3dc8b99857e28d'
    payload=json.loads((OUT/'observables.json').read_text());rows=payload['rows']
    assert {r['material'] for r in rows}=={'PdD','NiH'}
    fields=('material','isotope','temperature_K','variant','msd_A2','kinetic_eV','probability_r_gt_048','probability_r_gt_080')
    old_geometry=well.prior.OLD['profile']['geometry']
    geometry=well.WIDE['geometry']['PdD']
    assert geometry['moving_index']==old_geometry['moving_index']==25
    for key in ('H_indices','atomic_numbers','cell_A','pbc','reference_positions_A'):
        assert np.array_equal(geometry[key],old_geometry[key])
    with (HERE/'tables/observables.csv').open() as f:csvrows=list(csv.DictReader(f))
    assert len(rows)==len(csvrows)
    for row,tab in zip(rows,csvrows):
        assert all(k in row for k in fields)
        for key in fields:
            if row[key] is None:assert tab[key]==''
            elif key in ('material','isotope','variant'):assert row[key]==tab[key]
            else:assert float(tab[key])==row[key]
        assert row['msd_A2'] is None or row['msd_A2']>0
        assert row['kinetic_eV'] is None or row['kinetic_eV']>0
        for key in ('probability_r_gt_048','probability_r_gt_080'):
            assert row[key] is None or 0<=row[key]<=1
    norms=[];marginal=[]
    for T in (300,77):
        for kind in ('position','momentum'):
            data=np.genfromtxt(HERE/f'tables/D_{T}_{kind}.csv',delimiter=',',names=True)
            for v in ('old','quadratic','quartic',*well.VARIANTS):
                a=data['coordinate'];d=data[f'{v}_density'];m=data[f'{v}_measure']
                norm=float(np.trapezoid(m,a));assert abs(norm-1)<2e-4
                assert np.max(abs(m-4*np.pi*a*a*d))<1e-11
                norms.append(dict(temperature_K=T,variant=v,kind=kind,normalization=norm))
        data=np.genfromtxt(HERE/f'tables/D_{T}_signed111.csv',delimiter=',',names=True)
        for v in ('old','quadratic','quartic',*well.VARIANTS):
            norm=float(np.trapezoid(data[v],data['s_A']));sym=float(np.max(abs(data[v]-data[v][::-1])))
            assert abs(norm-1)<2e-4 and sym<1e-10
            marginal.append(dict(temperature_K=T,variant=v,normalization=norm,max_signed_reflection_defect=sym))
    p=json.loads((OUT/'potential.json').read_text())
    assert all(v['min_slope_eV_A']>0 and v['min_energy_eV']>0 for v in p['radial_positivity'].values())
    assert max(p['max_abs_force_plus_gradient_eV_A'].values())<2e-6
    max_anchor_error=0.
    for a in p['anchors']:
        for v in well.VARIANTS:
            pred=a['predictions'][v]
            max_anchor_error=max(max_anchor_error,abs(pred['energy_eV']-a['energy_eV']),abs(pred['outward_slope_eV_A']-a['slope_eV_A']))
    assert max_anchor_error<1e-11
    u=np.array([.71,.43,.56]);u/=np.linalg.norm(u)
    x=np.array([.017,.18,.37,.67,.80])[:,None]*u
    shared_error=max(float(np.max(abs(a-b))) for a,b in zip(well.energy_force(x,'wide_quadratic'),well.energy_force(x,'wide_quartic')))
    assert shared_error==0.
    coupled=json.loads((OUT/'coupled.json').read_text())
    for case in coupled['cases']:
        assert case['original_atom_index']==25 and case['frames']==801
        assert case['old_callback_control_max_abs_eV']<1e-12
        for v in case['variants'].values():
            assert v['overlap']['frame_ess']>=50 and v['overlap']['block_ess']>=4
            assert v['shifted_target_control_abs_error']<1e-12
    for r in coupled['observable_rows']:
        assert r['probability_r_gt_080'] is None and r['kinetic_eV'] is None
    nih_rows=json.loads((HERE/'nih/observable_table.json').read_text())
    joined_nih=[r for r in rows if r['material']=='NiH']
    assert [{k:r[k] for k in fields} for r in joined_nih]==nih_rows
    nih_support=json.loads((HERE/'nih/support.json').read_text())
    assert nih_support['model']['signed_export_sha256']==source_sha
    assert nih_support['model']['original_moving_index']==25
    assert [a['sample_kind'] for a in nih_support['model']['signed_anchors']].count('native_scf')==3
    for evidence in nih_support['isolated'].values():
        for key in ('position_radial_norm','momentum_radial_norm','signed_111_marginal_norm'):
            assert abs(evidence[key]-1)<3e-4
    for case in nih_support['finite_P_tagged'].values():
        assert case['frames_with_sampled_tail_080']==0
        for evidence in case['variants'].values():
            assert evidence['frame_weight_ess']>=50 and evidence['block_weight_ess']>=4
            assert evidence['tagged_virial_fd_max_abs_error_eV']<1e-8
    result=dict(scope='settled joined table/readback, changed PdD distribution and exact-input consumer only; not independent physical validation',
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,source_sha256=source_sha,
        joined_rows=len(rows),materials=sorted({r['material'] for r in rows}),
        csv_readback='exact common scalar/null agreement',normalizations=norms,signed_plane_marginals=marginal,
        exact_anchor_max_abs_error=max_anchor_error,shared_region_variant_max_abs_error=shared_error,
        coupled='consumed recorded graph/callback/direct-shift controls; no sampling or characteristic rerun',
        NiH=dict(consumer='source-linked shared fields read back exactly; recorded owner normalization/virial/overlap controls consumed without a numerical rerun',
                 owner_table_sha256=hashlib.sha256((HERE/'nih/observable_table.json').read_bytes()).hexdigest(),
                 owner_support_sha256=hashlib.sha256((HERE/'nih/support.json').read_bytes()).hexdigest(),
                 owner_distributions_sha256=hashlib.sha256((HERE/'nih/distributions.npz').read_bytes()).hexdigest()))
    (OUT/'verification.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
