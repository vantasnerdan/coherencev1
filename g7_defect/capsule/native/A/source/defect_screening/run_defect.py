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


def load_relaxed(directory,spec):
    result=json.loads((directory/'relaxation.json').read_bytes())
    require(result['passed'] and result['raw_fmax_eV_A']<=.05,'Input geometry is not force-qualified')
    last=result['frames'][-1]
    path=directory/last['gpw']['path']
    require(manifest(path)==last['gpw'],'Changed original relaxed GPW')
    calc=GPAW(str(path),legacy_gpaw=True,communicator=serial_comm,txt=None)
    atoms=calc.get_atoms()
    require(np.array_equal(atoms.cell.array,result['cell_A']) and np.allclose(atoms.positions,result['positions_A'],rtol=0,atol=1e-12),'Wrong saved relaxed geometry')
    source=json.loads((directory/('step-'+str(last['step']).zfill(3)+'.json')).read_bytes())
    require(source['specification']==spec,'Changed physical settings on a reused relaxed state')
    return atoms,result


def validate_epoch(input_path):
    source_manifest=input_path.parent/'source-manifest.json'
    require(source_manifest.is_file(),'Missing frozen executable source epoch')
    source=json.loads(source_manifest.read_bytes())
    for item in source['files']:
        path=input_path.parent/'source'/item['path']
        require(manifest(path,name=item['path'])==item,'Frozen executable source changed: '+item['path'])
    return manifest(input_path),manifest(source_manifest)


def run(input_path,output,radii,relaxed):
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
    if relaxed:
        atoms,relaxation=load_relaxed(relaxed,spec)
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
    parser.add_argument('--relaxed',type=Path)
    args=parser.parse_args()
    run(args.input,args.output,args.radii,args.relaxed)
