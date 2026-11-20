"""Use the accepted P255 constructor; only primitive-site band padding is explicit."""
from pathlib import Path
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from collision_ground import make_calculator, energy_record
from screen_io import manifest, write_json, require


def calculator(atoms,spec,log,charge,*,relaxation=False):
    family=dict(h_A=spec['h_A'],repeat=[spec['repeat_primitive']]*3,kmesh=spec['kmesh'])
    calc=make_calculator(atoms,family,log,charge=charge)
    nvalence=sum(10 if atom.symbol=='Pd' else 1 for atom in atoms)-charge
    # Two empty bands/site = eight/conventional fcc cell, not eight/primitive cell.
    experimental=dict(calc.parameters.experimental)
    experimental['reuse_wfs_method']='paw' if relaxation else 'keep'
    return calc.new(nbands=int(np.ceil(nvalence/2))+2*spec['original_Pd_sites'],
                    txt=str(log),parallel={'domain':1,'band':1,'kpt':1},experimental=experimental,
                    setups={'Pd':'10'} if relaxation else {'Pd':'10','H':'ae'})


def describe(atoms,spec,charge):
    calc=atoms.calc
    electrons=float(calc.get_number_of_electrons())
    expected=sum(10 if atom.symbol=='Pd' else 1 for atom in atoms)-charge
    require(abs(electrons-expected)<1e-7,'Wrong neutral-host/probe electron convention')
    require(bool(calc.scf.converged),'Native SCF did not converge')
    iterations=calc.get_number_of_iterations()
    return dict(schema='pd-defect-native-state-v1',energy=energy_record(calc),
        charge=charge,electrons=electrons,symbols=atoms.get_chemical_symbols(),
        positions_A=atoms.positions.tolist(),cell_A=atoms.cell.array.tolist(),
        realspace_shape=calc.wfs.gd.N_c.tolist(),actual_kpoints=int(calc.wfs.kd.nibzkpts),
        nbands=int(calc.get_number_of_bands()),specification=spec,scf_converged=True,
        xc=calc.parameters.xc,setups=dict(calc.parameters.setups),
        wavefunction_predictor=type(calc.wfs.wfs_mover).__name__,
        scf_iterations=None if iterations is None else int(iterations))


def save_state(atoms,spec,charge,output,stem,extra=None):
    checkpoint=output/(stem+'.gpw')
    temporary=output/(stem+'.writing.gpw')
    atoms.calc.write(str(temporary),mode='all')
    temporary.replace(checkpoint)
    record=describe(atoms,spec,charge)
    record['gpw']=manifest(checkpoint)
    if extra:record.update(extra)
    write_json(output/(stem+'.json'),record)
    return record


def scf(atoms,spec,charge,output,stem):
    start=time.process_time()
    output.mkdir(parents=True,exist_ok=True)
    if atoms.calc is None:
        atoms.calc=calculator(atoms,spec,output/(stem+'.txt'),charge)
    print('DEFECT_SCF_START case='+spec['case']+' state='+stem,flush=True)
    atoms.get_potential_energy()
    record=save_state(atoms,spec,charge,output,stem,dict(process_cpu_s=time.process_time()-start))
    print('DEFECT_SCF_COMPLETE case='+spec['case']+' state='+stem+' cpu_s='+str(record['process_cpu_s']),flush=True)
    return record
