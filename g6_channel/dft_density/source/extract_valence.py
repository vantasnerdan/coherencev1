"""Read a bound pristine GPW; reconstruct valence only without running an SCF."""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import time
import numpy as np
from gpaw import GPAW
from gpaw.mpi import serial_comm
from ase.units import Bohr


def binding(path):
    hasher=hashlib.sha256();size=0
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(4*1024*1024),b''):
            size+=len(chunk);hasher.update(chunk)
    return dict(bytes=size,sha256=hasher.hexdigest())


def extract(record,output):
    start=time.perf_counter();cpu=time.process_time()
    path=Path(record['path']);actual=binding(path)
    if actual!={key:record[key] for key in ('bytes','sha256')}:
        raise ValueError('Changed original pristine GPW')
    if record['charge']!=0 or record['electrons']!=44 or record['symbols'].count('Pd')!=4 or record['symbols'].count('H')!=4:
        raise ValueError('Source is not pristine neutral PdD')
    output.mkdir(parents=True,exist_ok=False)
    calc=GPAW(str(path),legacy_gpaw=True,communicator=serial_comm,txt=str(output/'reader.log'))
    atoms=calc.atoms;cell=atoms.cell.array
    if atoms.get_chemical_symbols()!=record['symbols'] or not np.allclose(atoms.positions,record['positions_A'],rtol=0,atol=1e-12):
        raise ValueError('Stored host geometry changed')
    if not np.array_equal(cell,record['cell_A']) or list(calc.wfs.gd.N_c)!=record['realspace_shape']:
        raise ValueError('Stored host cell/grid changed')
    if calc.get_number_of_electrons()!=44 or calc.parameters.charge!=0 or calc.wfs.kd.nibzkpts!=8:
        raise ValueError('Wrong pristine electron/charge/k support')
    if calc.parameters.xc!='PBE' or calc.parameters.setups!={'Pd':'10','H':'ae'}:
        raise ValueError('Wrong admitted functional/setup')
    density=calc.get_all_electron_density(gridrefinement=2,skip_core=True)
    volume=float(np.linalg.det(cell));electrons=float(np.sum(density,dtype=np.float64)*volume/density.size)
    error=abs(electrons-44.);finite=bool(np.isfinite(density).all())
    radii=[max(setup.rcut_j)*Bohr if setup.ni else 0. for setup in calc.wfs.setups]
    spline_cutoffs=[max(wave.get_cutoff() for wave in setup.get_partial_waves()[0])*Bohr if setup.ni else 0. for setup in calc.wfs.setups]
    npz=output/'valence-density.npz'
    np.savez_compressed(npz,density_valence_e_A3=density,cell_A=cell,positions_A=atoms.positions,
                        symbols=np.asarray(atoms.get_chemical_symbols()),origin_A=np.zeros(3))
    sources={Path(file).name:binding(file) for file in (inspect.getfile(type(calc)),inspect.getfile(type(calc.density)))}
    result=dict(schema='pdd-pristine-valence-density-v1',case=record['name'],original_gpw=actual,
        cell_A=cell.tolist(),positions_A=atoms.positions.tolist(),symbols=atoms.get_chemical_symbols(),
        electrons_expected=44,electrons_integrated=electrons,integration_error_e=error,count_exact=bool(error<=1e-6),
        integral_policy='Native skip_core count is reported, not corrected/renormalized; Axiom accepts its finite-grid PAW count discrepancy for the open-channel use.',
        augmentation_matching_radii_A=radii,partial_wave_spline_cutoff_A=spline_cutoffs,
        grid_shape=list(density.shape),coarse_shape=calc.wfs.gd.N_c.tolist(),actual_kpoints=int(calc.wfs.kd.nibzkpts),
        xc=calc.parameters.xc,setups=calc.parameters.setups,charge=0,skip_core=True,
        minimum_density_e_A3=float(np.min(density)),density_npz=binding(npz),native_reader_sources=sources,
        process_cpu_s=time.process_time()-cpu,wall_s=time.perf_counter()-start,
        passed=finite,scope='Native all-electron VALENCE as returned, frozen Pd core excluded; no charge correction, rescale or new SCF. Centre h-convergence is the deciding accuracy check.')
    (output/'density.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    if not result['passed']:
        raise ValueError('Nonfinite Native valence density')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bindings',type=Path,required=True);parser.add_argument('--case',choices=['h020','h016'],required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    records=json.loads(args.bindings.read_bytes())
    extract(next(record for record in records if record['name']==args.case),args.output)
