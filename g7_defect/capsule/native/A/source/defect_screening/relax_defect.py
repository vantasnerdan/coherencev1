"""Relax all ions, retaining each completed Native force/GPW boundary."""
import time
import numpy as np
from ase.optimize import BFGS
from defect_geometry import build, vacancy_record, choose_path
from defect_calculator import calculator, save_state, write_json, require


def relax(spec,output):
    output.mkdir(parents=True,exist_ok=False)
    atoms,initial=build(spec)
    write_json(output/'initial.json',initial)
    atoms.calc=calculator(atoms,spec,output/'relax-scf.txt',charge=0,relaxation=True)
    optimizer=BFGS(atoms,restart=str(output/'optimizer.json'),logfile=str(output/'optimizer.log'),trajectory=str(output/'relax.traj'))
    started=time.process_time();frames=[]
    def boundary():
        forces=atoms.get_forces(apply_constraint=False)
        require(np.isfinite(forces).all(),'Nonfinite unconstrained Native force')
        fmax=float(np.linalg.norm(forces,axis=1).max())
        record=save_state(atoms,spec,0,output,'step-'+str(optimizer.nsteps).zfill(3),
            dict(ionic_step=optimizer.nsteps,forces_eV_A=forces.tolist(),raw_fmax_eV_A=fmax,process_cpu_s=time.process_time()-started))
        frames.append(dict(step=optimizer.nsteps,raw_fmax_eV_A=fmax,energy=record['energy'],gpw=record['gpw']))
        print('DEFECT_RELAX_STEP case='+spec['case']+' step='+str(optimizer.nsteps)+' raw_fmax='+str(fmax),flush=True)
    optimizer.attach(boundary,interval=1)
    completed=optimizer.run(fmax=.05,steps=spec.get('max_ionic_steps',40))
    forces=atoms.get_forces(apply_constraint=False)
    fmax=float(np.linalg.norm(forces,axis=1).max())
    result=dict(initial=initial,frames=frames,raw_fmax_eV_A=fmax,passed=bool(completed and fmax<=.05),
                process_cpu_s=time.process_time()-started,positions_A=atoms.positions.tolist(),
                symbols=atoms.get_chemical_symbols(),cell_A=atoms.cell.array.tolist(),target_index=initial['target_index'])
    result['vacancy']=vacancy_record(atoms,initial)
    result['path']=choose_path(atoms,initial,result['vacancy'])
    write_json(output/'relaxation.json',result)
    require(result['passed'],'Ionic relaxation did not reach raw all-ion fmax0.05eV/A')
    return atoms,result


