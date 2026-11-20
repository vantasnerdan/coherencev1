"""No-SCF Native initialization, real primitive geometry and wrong-case refusal."""
import argparse
import gc
import json
from pathlib import Path
import numpy as np
from defect_geometry import build,vacancy_record
from defect_calculator import calculator,require,write_json,manifest
from run_defect import validate_epoch


def check(input_path,report):
    input_binding,source_binding=validate_epoch(input_path)
    spec=json.loads(input_path.read_bytes());cases={}
    for name in (spec['case'],):
        current=spec|{'case':name}
        atoms,initial=build(current)
        sites=spec['original_Pd_sites']
        expected=(sites,1) if name=='A' else ((sites-1,1) if name=='C' else (sites-1,sites))
        require((sum(atoms.numbers==46),sum(atoms.numbers==1))==expected,'Wrong Pd/D inventory')
        distances=atoms.get_all_distances(mic=True);np.fill_diagonal(distances,np.inf)
        require(float(distances.min())>1e-6,'Coalesced initial nuclei')
        require(abs(atoms.get_volume()-sites*spec['lattice_A']**3/4)<1e-9,'Wrong primitive fcc volume')
        calc=calculator(atoms,current,report.parent/('preflight-'+name+'.txt'),charge=0,relaxation=True)
        calc.initialize(atoms)
        require(all(setup.ni>0 for atom,setup in zip(atoms,calc.wfs.setups) if atom.symbol=='H'),'Geometry is not standard-H PAW')
        record=dict(Pd=expected[0],D=expected[1],volume_A3=atoms.get_volume(),
            grid_shape=calc.wfs.gd.N_c.tolist(),nbands=int(calc.get_number_of_bands()),
            electrons=float(calc.get_number_of_electrons()),kpoints=int(calc.wfs.kd.nibzkpts),
            minimum_initial_nuclear_distance_A=float(distances.min()),vacancy=vacancy_record(atoms,initial))
        record['geometry_setups']=dict(calc.parameters.setups)
        record['setup_files']=[]
        for setup in calc.wfs.setups:
            filename=getattr(setup,'filename',None)
            if filename and Path(filename).is_file():
                bound=manifest(Path(filename))
                if bound not in record['setup_files']:record['setup_files'].append(bound)
        ae=calculator(atoms,current,report.parent/('preflight-'+name+'-AE.txt'),charge=0)
        ae.initialize(atoms)
        require(all(setup.ni==0 for atom,setup in zip(atoms,ae.wfs.setups) if atom.symbol=='H'),'BO is not Hae')
        record['BO_setups']=dict(ae.parameters.setups)
        del ae
        cases[name]=record
        del calc,atoms;gc.collect()
    try:build(spec|{'case':'B'})
    except ValueError:refused=True
    else:raise ValueError('Forbidden fresh B geometry was accepted')
    result=dict(schema='pd-defect-native-preflight-v1',input=input_binding,frozen_source=source_binding,
                cases=cases,wrong_fresh_B_refused=refused,passed=True,scope='Real Native setup/grid/geometry initialized; no SCF, force or physical defect trend claimed.')
    write_json(report,result)
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True);parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args();check(args.input,args.report)
