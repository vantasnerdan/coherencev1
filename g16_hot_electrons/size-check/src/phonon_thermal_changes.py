"""Deliver matrix differences relative to the matched cold PdD crystal, without a shell fit."""
import argparse
from hashlib import file_digest
import json
from pathlib import Path

import numpy as np
from phonopy.file_IO import read_force_constants_hdf5, write_force_constants_to_hdf5

from native import require, write_once
from phonon_geometry import CELLS


def input_file(path):
    with path.open('rb') as stream:
        return dict(file=path.name,bytes=path.stat().st_size,sha256=file_digest(stream,'sha256').hexdigest())


def load(root, cell, mesh, temperature):
    summary=root/'dispersion.json'
    cfg=json.loads(summary.read_bytes())
    require(cfg['schema']=='g16-pdd-full-zone-dispersion-v1'
            and (cfg['cell'],cfg['primitive_mesh_n'],cfg['temperature_K'])==(cell,mesh,temperature),
            'Mismatched supercell, electronic mesh or temperature')
    path=root/'force_constants.hdf5'
    fc=read_force_constants_hdf5(path)
    require(fc.ndim==4 and fc.shape[0]==fc.shape[1] and fc.shape[2:]==(3,3) and np.isfinite(fc).all(),
            'Squared interaction tensor must contain both ionic species and three directions')
    samples=root/'mesh.npz'
    with np.load(samples,allow_pickle=False) as data:
        q, w=data['qpoints'], data['signed_frequencies_THz']
    require(q.shape==(24**3,3) and w.shape==(24**3,6) and np.isfinite(w).all(),'Incomplete phonon mesh')
    return dict(q=q,w=w,fc=fc,sources=[input_file(p) for p in (summary,path,samples)])


def stats(value):
    return dict(max_abs=float(np.abs(value).max()),rms=float(np.sqrt(np.mean(np.square(value)))))


def compare(cold,hot):
    require(np.array_equal(cold['q'],hot['q']) and cold['fc'].shape==hot['fc'].shape,'Incompatible lattice or phonon mesh')
    return hot['fc']-cold['fc'], np.sort(hot['w'],axis=1)-np.sort(cold['w'],axis=1)


def run(a):
    t=a.temperature_K
    require(t>300,'The hot temperature must exceed the common 300 K baseline')
    cold=[load(a.root/a.cell/f'k{n}'/'300',a.cell,n,300) for n in (12,16)]
    hot=[load(a.root/a.cell/f'k{n}'/str(t),a.cell,n,t) for n in (12,16)]
    require(np.array_equal(cold[0]['q'],cold[1]['q']),'The meshes use different phonon q grids')
    ch=[compare(c,h) for c,h in zip(cold,hot,strict=True)]
    a.output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for n, (k,w), c, h in zip((12,16),ch,cold,hot,strict=True):
        fc_path=a.output/f'delta_force_constants_k{n}.hdf5'
        write_force_constants_to_hdf5(k,filename=fc_path)
        rows.append(dict(primitive_mesh_n=n,sorted_frequency_change_THz=stats(w),ifc_change_eV_A2=stats(k),
                         cold_sources=c['sources'],hot_sources=h['sources'],delta_tensor=input_file(fc_path)))
    write_once(a.output/'thermal_change.json',dict(schema='g16-pdd-thermal-ifc-delta-v1',cell=a.cell,
               cold_temperature_K=300,hot_temperature_K=t,measure='holonomic force-constant difference',rows=rows,
               kmesh_difference_of_IFC_change_eV_A2=stats(ch[1][0]-ch[0][0]),
               kmesh_difference_of_sorted_frequency_change_THz=stats(ch[1][1]-ch[0][1]),
               note='Sorted frequency change compares unordered spectra; it does not track branches through crossings.'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--cell',choices=CELLS,required=True)
    p.add_argument('--temperature-K',type=int,required=True)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args())


if __name__=='__main__':
    main()
