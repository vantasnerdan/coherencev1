"""Write full PdD IFCs and Phonopy spectra using saved GPAW displacements."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from phonopy.file_IO import write_force_constants_to_hdf5
import yaml

from native import require, write_once
from phonon_forces import assemble
from phonon_geometry import CELLS, TEMPERATURES

FOLDED_SHA = '34ff5148c979266230f95db3b00b9eb56a552dc0de9e91fe79e0913f58d8750a'
Q_SPECIAL = np.array([[0.,0.,0.],[0.,.25,.25],[0.,.5,.5],[.5,.5,.5]])
LABELS = ('gamma','x_half','x','l')


def matrix_residuals(fc):
    return dict(row_sum_max_eV_A2=float(np.abs(fc.sum(axis=1)).max()),
                permutation_max_eV_A2=float(np.abs(fc-fc.transpose(1,0,3,2)).max()))


def special_points(ph, original=None):
    data=ph.run_qpoints(Q_SPECIAL,with_dynamical_matrices=True)
    mass=ph.primitive.masses
    d=list(ph.primitive.symbols).index('H')
    part=data.dynamical_matrices[:,3*d:3*d+3,3*d:3*d+3]*mass[d]
    require(np.abs(part.imag).max()<1e-8,'A same-sublattice reciprocal diagonal is unexpectedly complex')
    blocks={label:part[i].real.tolist() for i,label in enumerate(LABELS)}
    table=dict(qpoints=Q_SPECIAL.tolist(), signed_frequencies_THz=dict(zip(LABELS,data.frequencies.tolist(),strict=True)),
               clamped_D_matrix_eV_A2=blocks)
    if original is not None:
        observed=np.array([part[0,0,0].real,part[2,0,0].real,part[2,1,1].real])
        expected=np.array([original[k]['K_force_eV_A2'] for k in ('gamma_100','x_100','x_010')])
        err=np.abs(observed-expected)
        lim=np.maximum(.002,.01*np.abs(expected))
        table['old_folded_clamped_D_crosscheck']=dict(order=['gamma_100','x_100','x_010'],
            expected_eV_A2=expected.tolist(), observed_eV_A2=observed.tolist(), abs_error_eV_A2=err.tolist(),
            tolerance_eV_A2=lim.tolist(), matches=(err<=lim).tolist())
    return table


def save_phonopy(ph, path):
    ph.save(path,settings={'force_constants':False})
    # Phonopy 4.11's YAML writer rounds isotope masses to six decimals.
    doc=yaml.safe_load(path.read_text())
    doc['phonopy']['frequency_unit_conversion_factor']=float(ph.unit_conversion_factor)
    for name, structure in (('unit_cell',ph.unitcell),('primitive_cell',ph.primitive),('supercell',ph.supercell)):
        for atom, mass in zip(doc[name]['points'],structure.masses,strict=True):
            atom['mass']=float(mass)
    path.write_text(yaml.safe_dump(doc,sort_keys=False))



def export(ph, out, original=None):
    raw=ph.force_constants.copy()
    before=matrix_residuals(raw)
    raw_points=special_points(ph,original)
    ph.symmetrize_force_constants(level=3,show_drift=False)
    fc=ph.force_constants
    after=matrix_residuals(fc)
    max_projection=float(np.abs(fc-raw).max())
    corrected=special_points(ph,original)
    write_force_constants_to_hdf5(raw,filename=out/'force_constants_raw.hdf5')
    write_force_constants_to_hdf5(fc,filename=out/'force_constants.hdf5')
    save_phonopy(ph,out/'phonopy.yaml')
    ph.auto_band_structure(npoints=81,with_eigenvectors=False,write_yaml=True,filename=out/'band.yaml')
    mesh=ph.run_mesh([24]*3,is_gamma_center=True,is_mesh_symmetry=False,is_time_reversal=False,with_eigenvectors=True)
    require(mesh.qpoints.shape==(24**3,3) and mesh.frequencies.shape==(24**3,6)
            and mesh.eigenvectors.shape==(24**3,6,6),'Incomplete full-zone mesh')
    np.savez_compressed(out/'mesh.npz',qpoints=mesh.qpoints,signed_frequencies_THz=mesh.frequencies,
                        eigenvectors=mesh.eigenvectors)
    return dict(raw_matrix=before, projected_matrix=after, max_projection_eV_A2=max_projection,
                raw_special_points=raw_points, projected_special_points=corrected,
                signed_mode_minimum_THz=float(mesh.frequencies.min()))


def old_comparator(path,n,t):
    content=path.read_bytes()
    require(sha256(content).hexdigest()==FOLDED_SHA,'Wrong independent folded comparison')
    return next((r for r in json.loads(content)['rows'] if r['primitive_mesh_n']==n and r['temperature_K']==t),None)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--cell',choices=CELLS,required=True)
    p.add_argument('--mesh',type=int,choices=(12,16),required=True)
    p.add_argument('--temperature-K',type=int,choices=TEMPERATURES,required=True)
    p.add_argument('--folded',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--compact-bands',action='store_true',help='Replay the seven-band anchor control')
    p.add_argument('--tight-forces',action='store_true',help='Replay the tighter force stopping control')
    p.add_argument('--producer-sources',type=Path,default=Path(__file__).resolve().parent,
                   help='Directory containing the exact source files used for these native points')
    a=p.parse_args()
    old=old_comparator(a.folded,a.mesh,a.temperature_K)
    ph,native=assemble(a.root,a.cell,a.mesh,a.temperature_K,a.producer_sources,
                       compact_bands=a.compact_bands,tight_forces=a.tight_forces)
    a.output.mkdir(parents=True,exist_ok=False)
    summary=export(ph,a.output,old)
    write_once(a.output/'dispersion.json',dict(schema='g16-pdd-full-zone-dispersion-v1',cell=a.cell,primitive_mesh_n=a.mesh,
               temperature_K=a.temperature_K, reference_folded_sha256=FOLDED_SHA,**native,**summary))


if __name__=='__main__':
    main()
