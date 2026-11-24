"""One settled changed-consumer check; consumes, never reruns, trajectory evidence."""
from pathlib import Path
import csv
import hashlib
import json
import platform
import numpy as np
import scipy
from well import HERE


def main():
    out=HERE/'results';report=json.loads((out/'report_table.json').read_text())
    numerical=report['numerical_comparison'];coupled=report['coupled']
    assert len(report['isolated'])==18
    assert coupled['well_py_sha256']==hashlib.sha256((HERE/'well.py').read_bytes()).hexdigest()
    checks=[]
    for T in (300,77):
        for var in ('old','quadratic','quartic'):
            z=np.load(out/f'D_{T}_{var}_main_distributions.npz')
            for kind,gridname,dname in [('position','x_positive_A','rho_positive_octant_A_minus3'),('momentum','k_positive_A_inverse','n_positive_octant_A3')]:
                grid=z[gridname];values=z[dname]
                norm=float(8*np.trapezoid(np.trapezoid(np.trapezoid(values,grid,axis=2),grid,axis=1),grid))
                assert abs(norm-1)<2e-6,(T,var,kind,norm)
                assert np.isfinite(values).all() and values.min()>=0
                records=np.genfromtxt(HERE/f'tables/D_{T}_{kind}.csv',delimiter=',',names=True)
                assert np.array_equal(records[f'{var}_density'],z[f'{kind}_density'])
                assert np.array_equal(records[f'{var}_measure'],z[f'{kind}_measure'])
                assert abs(float(np.trapezoid(z[f'{kind}_measure'],z[f'{kind}_radius']))-1)<2e-4
                checks.append(dict(temperature_K=T,variant=var,observable=kind,full_3D_normalization=norm,table_matches_actual_distribution=True))
            signed=z['s111_A'];p=z['marginal111_A_inverse'];norm=float(np.trapezoid(p,signed))
            assert abs(norm-1)<.002
            assert np.max(abs(p-p[::-1]))<1e-10
            checks.append(dict(temperature_K=T,variant=var,observable='signed111 plane marginal',normalization=norm,mean_A=float(np.trapezoid(p*signed,signed))))
    for row in numerical:
        if row['variant']=='old':
            assert abs(row['recreated_reporting_minus_archive_r2_A2'])<1e-9
            assert abs(row['recreated_reporting_minus_archive_kinetic_eV'])<1e-9
        assert abs(row['delta_r2_A2'])<2e-4
        assert abs(row['delta_kinetic_eV'])<1e-4
    potential=json.loads((out/'potential_evidence.json').read_text())
    assert all(v['min_slope_eV_A']>0 and v['min_energy_eV']>0 for v in potential['radial_segment_positivity'].values())
    assert all(v['max_abs_force_plus_energy_gradient_eV_A']<2e-5 for v in potential['gradient'].values())
    for c in coupled['cases']:
        assert c['frames']==801 and c['contiguous_blocks']==9
        assert all(v['overlap']['passes_participation_screen'] for v in c['variants'].values())
        assert c['old_characteristic_field_callback_vs_polynomial_max_abs_block_error']<1e-12
        for v in c['variants'].values():assert v['momentum_characteristic_111']['direct_witness_abs_error']<1e-12
    for filename in ('README.md','replay.md','wells.png','distributions.png','signed111.png','coupled_position.csv','coupled_characteristic_111.csv'):
        assert (HERE/filename).stat().st_size>0
    result=dict(status='passed',scope='settled D old/new quantum distributions, exported tables, signed marginal, archived reporting positive control and actual coupled handoff; no producer or coupled rerun',environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),checks=checks,limitations='empirical basis/quadrature/grid comparison, not rigorous numerical enclosure; no electronic error bound; three rays only; coupled tagged-site finite-P and existing sampled overlap only')
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
