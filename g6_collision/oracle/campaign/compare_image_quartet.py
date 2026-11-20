"""Consume the Axiom2291 fixed-phase R0.1A neutral-image quartet."""
import argparse
import json
import math
import re
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.units import Bohr
from collision_pair_readout import load_case
from static_state_data import load_state
from screen_io import manifest,require,write_json
from substrate_framework.collision_bare_reference import PeriodicNuclearGreen
from substrate_framework.collision_observables import COULOMB_EV_A


def load_quartet(control,baseline):
    states={};bindings={};base=baseline['ground'];spec=baseline['trajectory']['specification']
    for role in ('occupied','vacant'):
        for state in ('probe','host'):
            key=role+'-'+state
            record,bindings[key]=load_state(control,role,state,h_A=.2,repeat=[2,2,2])
            require(np.array_equal(record['cell_A'],2*np.asarray(base['cell_A'])),'Image cell is not the matched doubled cell')
            require(record['realspace_shape']==[40,40,40] and record['actual_kpoints']==1,'Wrong image grid/k support')
            for name in ('D_displacement_A','direction','impact_A','relaxed_ground_sha256'):
                require(record['specification'][name]==spec[name],'Image also changed '+name)
            text=(control/'calculation'/key/'scf.txt').read_text()
            bands=re.search(r'^Number of bands in calculation:\s+(\d+)',text,re.M)
            require(bands is not None and int(bands[1])==240,'Actual SCF lacks full folded band support')
            states[key]=record
        probe=states[role+'-probe'];host=states[role+'-host'];index=probe['path']['projectile_index']
        require(probe['symbols'][:index]+probe['symbols'][index+1:]==host['symbols'] and
                np.allclose(np.delete(probe['positions_A'],index,axis=0),host['positions_A'],rtol=0,atol=1e-11),
                'Probe/neutral-host geometries do not differ by the projectile alone')
        require(probe['electrons']==host['electrons']==(352 if role=='occupied' else 351),'Wrong neutral electron convention')
    host=states['occupied-host'];vac=states['vacant-host'];target=np.asarray(states['occupied-probe']['path']['target_A'])
    original=baseline['host']
    repeated=Atoms(original['symbols'],positions=original['positions_A'],cell=original['cell_A'],pbc=True).repeat((2,2,2))
    require(host['symbols']==repeated.get_chemical_symbols() and np.allclose(host['positions_A'],repeated.positions,rtol=0,atol=1e-11),
            'Image occupied host is not the original repeated frozen material')
    removed=[i for i,(symbol,pos) in enumerate(zip(host['symbols'],host['positions_A'])) if symbol=='H' and np.linalg.norm(np.asarray(pos)-target)<1e-10]
    require(len(removed)==1,'Image host has no unique original target')
    index=removed[0]
    require(host['symbols'][:index]+host['symbols'][index+1:]==vac['symbols'] and
            np.allclose(np.delete(host['positions_A'],index,axis=0),vac['positions_A'],rtol=0,atol=1e-11),
            'Image vacancy removes more than the one neutral target')
    return states,bindings


def compare(baseline_run,control,reference,bare_path):
    base=load_case(baseline_run)
    require(base['trajectory']['specification']['D_displacement_A']==[0,0,0],'Wrong baseline phase')
    states,bindings=load_quartet(control,base)
    old=json.loads(reference.read_bytes())
    require(old['source_bindings']['occupied']==base['bindings'],'Image baseline source epoch differs')
    oldrow=[r for r in old['rows'] if r['R_A']==.1]
    require(len(oldrow)==1,'No unique originalR0.1 pair point');oldrow=oldrow[0]
    probe=states['occupied-probe'];index=probe['path']['projectile_index']
    point=np.asarray(probe['positions_A'][index]);target=np.asarray(probe['path']['target_A']);vector=point-target
    require(np.allclose(point,base['trajectory']['rows'][-1]['projectile_A'],rtol=0,atol=1e-11) and
            np.allclose(point,states['vacant-probe']['positions_A'][-1],rtol=0,atol=1e-11),'Target/grid phase or projectile point changed')
    require(abs(np.linalg.norm(vector)-.1)<1e-10,'Wrong actual image radius')
    bare=json.loads(bare_path.read_bytes());require(bare['passed'],'Unqualified bare regularizer')
    alpha_a,alpha_b=bare['native_alpha1_Bohr_inverse2'];width=np.sqrt(alpha_a*alpha_b/(alpha_a+alpha_b))/Bohr
    small=PeriodicNuclearGreen(base['ground']['cell_A']);large=PeriodicNuclearGreen(probe['cell_A'])
    signs={'occupied-probe':1,'vacant-probe':-1,'occupied-host':-1,'vacant-host':1}
    raw=math.fsum(signs[key]*record['energy']['internal_energy_eV'] for key,record in states.items())
    correction=large.point_minus_gaussian_eV(vector,width);total=raw+correction
    lowering=COULOMB_EV_A/.1-total;delta=lowering-oldrow['candidate_lowering_vs_isolated_C_over_R_eV']
    tolerance=max(5.,.1*abs(lowering))
    return dict(schema='pdd-neutral-image-quartet-comparison-v1',authority='Axiom2291',R_A=.1,
        raw_neutral_quartet_2x_eV=raw,point_minus_soft_2x_eV=correction,V_pair_periodic_point_2x_eV=total,
        candidate_lowering_base_eV=oldrow['candidate_lowering_vs_isolated_C_over_R_eV'],candidate_lowering_2x_eV=lowering,
        delta_2x_minus_base_lowering_eV=delta,tolerance_eV=tolerance,passed=bool(abs(delta)<=tolerance),
        raw_quartet_delta_eV=raw-oldrow['BO_quartet_soft_eV'],point_minus_soft_delta_eV=correction-oldrow['bare_periodic_point_minus_soft_eV'],
        periodic_bare_point_kernel_delta_eV=large.point_eV(vector)-small.point_eV(vector),
        state_bindings=bindings,baseline=manifest(reference),bare_qualification=manifest(bare_path),
        scope='One baseline-phaseR0.1 representative cell check. A pass admits original-cell BO curves under2291, not a prepared state or dt/drag convergence. The reported bare point-kernel delta is diagnostic and is not subtracted from the total neutral-material comparison.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('baseline_occupied','control','reference','bare','output'):
        parser.add_argument('--'+name.replace('_','-'),type=Path,required=True)
    args=parser.parse_args()
    write_json(args.output,compare(args.baseline_occupied,args.control,args.reference,args.bare))
