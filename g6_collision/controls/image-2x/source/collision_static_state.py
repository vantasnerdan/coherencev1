"""One real fixed-geometry state for the approved P255 representation quartet."""
import argparse
import json
from pathlib import Path
from collections import deque
from time import monotonic
from screen_io import manifest,write_json,require
from collision_geometry import collision_atoms
from collision_ground import make_calculator,energy_record


def attach_iteration_limit(calc,spec):
    limit=spec.get('maximum_five_scf_iterations_wall_s')
    if limit is None:
        return
    ticks=deque(maxlen=6)
    def monitor():
        ticks.append(monotonic())
        if len(ticks)<6:
            return
        elapsed=ticks[-1]-ticks[0]
        print(f'COLLISION_STATIC_COST five_iterations_wall_s={elapsed:.3f} limit_s={limit:.3f}',flush=True)
        if elapsed>limit:
            raise RuntimeError('COLLISION_STATIC_3X_PAUSE: stop investment; prior completed states survive, current SCF resume is unproven; no automatic retry')
    calc.attach(monitor,n=1)


def calculate(spec,role,state,output):
    require(manifest(spec['relaxed_ground_record'])['sha256']==spec['relaxed_ground_sha256'],
            'Changed relaxed host source')
    spec=dict(spec,occupied_target=(role=='occupied'))
    atoms,path=collision_atoms(spec)
    charge=1
    if state=='host':
        del atoms[path['projectile_index']]
        charge=0
    output.mkdir(parents=True,exist_ok=False)
    atoms.calc=make_calculator(atoms,spec,output/'scf.txt',charge=charge)
    attach_iteration_limit(atoms.calc,spec)
    print(f'COLLISION_STATIC_STATE_START role={role} state={state}',flush=True)
    atoms.get_potential_energy()
    gpw=output/'state.gpw';atoms.calc.write(str(gpw),mode='all')
    result=dict(schema='pdd-static-quartet-state-v1',role=role,state=state,charge=charge,
        energy=energy_record(atoms.calc),electrons=float(atoms.calc.get_number_of_electrons()),
        symbols=atoms.get_chemical_symbols(),positions_A=atoms.positions.tolist(),cell_A=atoms.cell.array.tolist(),
        realspace_shape=atoms.calc.wfs.gd.N_c.tolist(),actual_kpoints=int(atoms.calc.wfs.kd.nibzkpts),
        specification=spec,path=path,gpw=manifest(gpw),
        scope='One fixed-geometry SCF state, not a trajectory or an isolated-material qualification')
    write_json(output/'state.json',result)
    print(f'COLLISION_STATIC_STATE_COMPLETE role={role} state={state}',flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--role',choices=['occupied','vacant'],required=True)
    parser.add_argument('--state',choices=['probe','host'],required=True)
    args=parser.parse_args()
    calculate(json.loads(args.input.read_bytes()),args.role,args.state,args.output)
