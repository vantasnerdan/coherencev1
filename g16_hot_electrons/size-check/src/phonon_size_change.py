"""Compare the thermal PdD response across supercells at fixed electronic sampling."""
import argparse
from pathlib import Path

import numpy as np
from phonopy import load as load_phonopy

from native import require, write_once
from phonon_geometry import CELLS, TEMPERATURES
from phonon_thermal_changes import input_file, load, stats


def response(root, cell, n, t):
    folder = root / cell / f'k{n}' / str(t)
    record = load(folder,cell,n,t)
    yml = folder/'phonopy.yaml'
    ph = load_phonopy(yml,force_constants_filename=folder/'force_constants.hdf5',produce_fc=False,is_nac=False)
    primitive = ph.primitive
    result = ph.run_qpoints(record['q'],with_dynamical_matrices=True)
    k = result.dynamical_matrices
    require(np.isfinite(k).all() and k.shape==(len(record['q']),6,6)
            and np.max(np.abs(k-k.conj().swapaxes(1,2)))<1e-10
            and np.max(np.abs(result.frequencies-record['w']))<1e-6,'Invalid force-constant basis or spectrum')
    m = np.sqrt(np.repeat(primitive.masses,3))
    k *= m[:,None]*m[None,:]
    basis = dict(symbols=list(primitive.symbols),masses_u=primitive.masses.tolist(),
                 cell_A=primitive.cell.tolist(),positions_fractional=primitive.scaled_positions.tolist())
    return dict(q=record['q'],frequencies=record['w'],k=k,basis=basis,identity=[*record['sources'],input_file(yml)])


def norm(k):
    return np.max(np.abs(np.linalg.eigvalsh(k)),axis=1)


def compare(sc, sh, lc, lh):
    require(all(s['basis']==sc['basis'] and np.array_equal(s['q'],sc['q']) for s in (sh,lc,lh)),
            'Use identical primitive basis, Cartesian axes, isotope masses and q grid')
    a = sh['k']-sc['k']
    b = lh['k']-lc['k']
    err = norm(b-a)
    ref = norm(b)
    tol = np.maximum(.002,.10*ref)
    cold = np.sort(lc['frequencies'],axis=1)-np.sort(sc['frequencies'],axis=1)
    hot = np.sort(lh['frequencies'],axis=1)-np.sort(sh['frequencies'],axis=1)
    maximum = int(err.argmax())
    ratio = int((err/tol).argmax())
    summary = dict(delta_K_operator_norm_error_eV_A2=stats(err),larger_delta_K_operator_norm_eV_A2=stats(ref),
                   cold_absolute_sorted_shift_THz=stats(cold),hot_absolute_sorted_shift_THz=stats(hot),
                   thermal_sorted_shift_difference_THz=stats(hot-cold),q_count=len(err),
                   worst_absolute_error_q=sc['q'][maximum].tolist(),worst_relative_error_q=sc['q'][ratio].tolist(),
                   worst_relative_error_to_tolerance=float(err[ratio]/tol[ratio]),
                   outside_tolerance_q_count=int(np.count_nonzero(err>tol)),
                   size_screen_passed=bool(np.all(err<=tol)))
    return summary, dict(qpoints=sc['q'],difference_operator_norm_eV_A2=err,
                         larger_delta_operator_norm_eV_A2=ref,tolerance_eV_A2=tol)


def run(a):
    require(a.smaller!=a.larger and a.temperature_K>300,'Choose distinct sizes and an elevated temperature')
    combos = [(s,t) for s in (a.smaller,a.larger) for t in (300,a.temperature_K)]
    responses = [response(a.root,s,a.mesh,t) for s,t in combos]
    summary, arrays = compare(*responses)
    a.output.mkdir(parents=True,exist_ok=False)
    detail = a.output/'q_sensitivity.npz'
    np.savez_compressed(detail,**arrays)
    write_once(a.output/'size_sensitivity.json',dict(schema='g16-pdd-size-sensitivity-v1',
               primitive_mesh_n=a.mesh,smaller=a.smaller,larger=a.larger,temperature_K=a.temperature_K,
               sources=[dict(cell=s,temperature_K=t,files=r['identity']) for (s,t),r in zip(combos,responses,strict=True)],
               q_sensitivity_file=input_file(detail),**summary))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--smaller',choices=CELLS,default='primitive2')
    p.add_argument('--larger',choices=CELLS,default='conventional2')
    p.add_argument('--mesh',type=int,choices=(12,16),default=16)
    p.add_argument('--temperature-K',type=int,choices=TEMPERATURES,required=True)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args())


if __name__=='__main__':
    main()
