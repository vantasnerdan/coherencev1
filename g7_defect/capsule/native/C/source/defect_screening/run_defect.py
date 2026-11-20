"""Finite preregistered Native cases; B is reused and never launched here."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from gpaw import GPAW
from gpaw.mpi import serial_comm
from relax_defect import relax
from defect_bo import quartet
from defect_calculator import manifest,write_json,require
from defect_geometry import vacancy_record,choose_path


def load_geometry_frame(frame,spec):
    source=json.loads(frame.read_bytes())
    require(source['scf_converged'] and source['setups']=={'Pd':'10'},'Geometry must be an actual converged standard-H PAW state')
    require(source['raw_fmax_eV_A']<=.1,'Input geometry exceeds Axiom2518 force threshold0.1eV/A')
    physical_source={key:value for key,value in source['specification'].items() if key!='h_A'}
    requested={key:value for key,value in spec.items() if key!='h_A'}
    require(physical_source==requested,'Changed physical settings on a reused geometry')
    path=frame.parent/source['gpw']['path']
    require(manifest(path)==source['gpw'],'Changed original geometry GPW')
    calc=GPAW(str(path),legacy_gpaw=True,communicator=serial_comm,txt=None)
    atoms=calc.get_atoms()
    require(atoms.get_chemical_symbols()==source['symbols'] and
            np.array_equal(atoms.cell.array,source['cell_A']) and
            np.allclose(atoms.positions,source['positions_A'],rtol=0,atol=1e-12),'Wrong saved geometry')
    forces=atoms.get_forces(apply_constraint=False)
    require(np.isfinite(forces).all() and np.allclose(forces,source['forces_eV_A'],rtol=0,atol=1e-5),'Saved force boundary did not reload')
    fmax=float(np.linalg.norm(forces,axis=1).max())
    require(fmax<=.1,'Reloaded geometry exceeds force threshold0.1eV/A')
    initial=json.loads((frame.parent/'initial.json').read_bytes())
    last_path=max(frame.parent.glob('step-*.json'),key=lambda item:int(item.stem.split('-')[-1]))
    last=json.loads(last_path.read_bytes())
    result=dict(initial=initial,frames=[dict(step=source['ionic_step'],raw_fmax_eV_A=fmax,
                energy=source['energy'],gpw=source['gpw'])],passed=True,raw_fmax_eV_A=fmax,
                force_threshold_eV_A=.1,authority='Axiom2518; original0.05 optimizer outcome is unchanged',
                source_frame=manifest(frame),source_last_completed_frame=manifest(last_path),
                source_last_completed_fmax_eV_A=last['raw_fmax_eV_A'],
                free_energy_drift_after_selected_eV=last['energy']['free_energy_eV']-source['energy']['free_energy_eV'],
                positions_A=atoms.positions.tolist(),symbols=atoms.get_chemical_symbols(),
                cell_A=atoms.cell.array.tolist(),target_index=initial['target_index'],
                scope='Selected actual PAW force boundary admitted at0.1eV/A, not a claim the original0.05 optimization passed.')
    result['vacancy']=vacancy_record(atoms,initial)
    result['path']=choose_path(atoms,initial,result['vacancy'])
    return atoms,result


def validate_epoch(input_path):
    source_manifest=input_path.parent/'source-manifest.json'
    require(source_manifest.is_file(),'Missing frozen executable source epoch')
    source=json.loads(source_manifest.read_bytes())
    for item in source['files']:
        path=input_path.parent/'source'/item['path']
        require(manifest(path,name=item['path'])==item,'Frozen executable source changed: '+item['path'])
    return manifest(input_path),manifest(source_manifest)


def run(input_path,output,radii,geometry_frame):
    started,wall=time.process_time(),time.perf_counter()
    input_binding,source_binding=validate_epoch(input_path)
    spec=json.loads(input_path.read_bytes())
    require(spec['case'] in ('A','C','D'),'B must use the qualified existing baseline')
    require(spec['geometry_H_setup']=='standard_PAW','Only the approved PAW-geometry/AE-BO route is current')
    require(spec['original_Pd_sites']==spec['repeat_primitive']**3,'Wrong primitive-site band support')
    allowed=(.1,.2,.3,.5,1.,1.6)
    require(radii and len(set(radii))==len(radii) and all(radius in allowed for radius in radii),'Wrong requested BO support')
    output.mkdir(parents=True,exist_ok=False)
    print('DEFECT_RUN_START case='+spec['case']+' stage=production',flush=True)
    if geometry_frame:
        atoms,relaxation=load_geometry_frame(geometry_frame,spec)
        write_json(output/'geometry-selection.json',relaxation)
    else:
        atoms,relaxation=relax(spec,output/'relax')
    pairs=quartet(atoms,relaxation,spec,radii,output/'bo')
    require(all(state['setups'].get('H')=='ae' for state in
            [pairs['states']['occupied-host'],pairs['states']['vacant-host'],
             *pairs['states']['occupied-probes'],*pairs['states']['vacant-probes']]),
            'Mixed PAW/AE quartet is forbidden')
    result=dict(schema='pd-defect-native-run-v1',stage='production',case=spec['case'],
        input=input_binding,frozen_source=source_binding,relaxation=relaxation,
        quartet=manifest(output/'bo/quartet.json'),
        rows=pairs['rows'],process_cpu_s=time.process_time()-started,wall_s=time.perf_counter()-wall,
        completed=True,scope='Standard-H PAW relaxed geometry, all four BO energies recomputed with Hae; C numerical admission remains separate.')
    write_json(output/'run-result.json',result)
    print('DEFECT_RUN_COMPLETE case='+spec['case']+' stage=production cpu_s='+str(result['process_cpu_s']),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--radii',type=float,nargs='+',required=True)
    parser.add_argument('--geometry-frame',type=Path)
    args=parser.parse_args()
    run(args.input,args.output,args.radii,args.geometry_frame)
