"""NumPy-only vacancy time-step control; production BO references are reused, not refitted."""
import argparse
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
import numpy as np
from capsule_io import Capsule,require
from replay import close,density_potential


BASE='runs/base-vac/calculation/'
REFINED='runs/dt1-vac/calculation/'


def read_control(root):
    capsule=Capsule(root)
    ground=capsule.json(BASE+'ground/ground.json')
    production=capsule.json(BASE+'trajectory/trajectory.json')
    refined=capsule.json(REFINED+'trajectory/trajectory.json')
    bo=capsule.json(BASE+'bo/excitation.json')
    capsule.bind(production['ground'],BASE+'ground/ground.json')
    capsule.bind(refined['ground'],BASE+'ground/ground.json')
    capsule.bind(bo['trajectory'],BASE+'trajectory/trajectory.json')
    require(bo['specification']==production['specification'],'Reused BO input epoch differs')
    a=production['specification'];b=refined['specification']
    require(a['dt_as']==2 and b['dt_as']==1,'Require actual half-dt trajectory')
    require({k:v for k,v in a.items() if k!='dt_as'}=={k:v for k,v in b.items() if k!='dt_as'},'Only time step may change')
    require(not a['occupied_target'] and a['D_displacement_A']==[0,0,0] and a['energy_keV']==2,'Require baseline vacancy only')
    require(a['h_A']==.2 and a['kmesh']==[2,2,2] and a['solver_tolerance']==1e-9,'Wrong admitted spatial/solver settings')
    require(ground['electrons']==bo['neutral_host_reference']['electrons']==43 and bo['neutral_host_reference']['charge']==0,'Wrong neutral BO reference')
    close(production['rows'][0]['td_internal_energy_eV'],refined['rows'][0]['td_internal_energy_eV'],'Identical initial GPW energy',1e-7)
    rows=refined['rows'];require(len(rows)>len(production['rows']),'Half-dt trajectory has no refined steps')
    require(all(y['time_as']>x['time_as'] and y['R_A']<x['R_A'] for x,y in zip(rows,rows[1:])),'Nonmonotone refined trajectory')
    positions=np.asarray([r['projectile_A'] for r in rows]);forces=np.asarray([r['projectile_force_eV_A'] for r in rows])
    work=np.r_[0.,np.cumsum(-np.einsum('ij,ij->i',.5*(forces[1:]+forces[:-1]),np.diff(positions,axis=0)))]
    work_difference=close(work,[r['external_work_eV'] for r in rows],'Actual refined force work')
    energy=np.asarray([r['td_internal_energy_eV'] for r in rows]);closure=work-(energy-energy[0])
    require(np.all(abs(closure)<=np.maximum(.5,.01*abs(work))),'Refined independent force-work closure fails')
    particles=np.asarray([r['integrated_all_electrons'] for r in rows]);drift=float(max(abs(particles-particles[0])))
    require(abs(particles[0]-187)<1e-5 and drift<1e-5,'Refined particle count/drift fails')
    operators=refined['fixed_projector_contractions']
    op=max(max(c['addition_max_abs'],c['integration_max_abs']) for c in operators['actual_native_checks'])
    require(operators['moved_projector_refused'] and len(operators['actual_native_checks'])==8 and op<1e-11,'Refined Native projector check/refusal fails')
    require(len(production['checkpoints'])==len(refined['checkpoints'])==len(bo['checkpoints'])==7,'Wrong checkpoint support')
    point_rows=[];density_error=0.;start_R=float(a['start_R_A'])
    for old,new,static in zip(production['checkpoints'],refined['checkpoints'],bo['checkpoints']):
        require(old['R_A']==new['R_A']==static['R_A'],'Unmatched reused BO radius')
        close(old['positions_A'],new['positions_A'],'Reused fixed-geometry BO positions')
        require(old['row_index']==static['row_index'],'Original BO row binding differs')
        oldrow=production['rows'][old['row_index']];newrow=rows[new['row_index']]
        close(oldrow['td_internal_energy_eV'],static['td_internal_energy_eV'],'Original BO/TD epoch')
        close(new['positions_A'][ground['path']['projectile_index']],newrow['projectile_A'],'Refined checkpoint projectile')
        close(newrow['R_A'],new['R_A'],'Refined row radius')
        Ebo=static['bo']['internal_energy_eV']
        Xold=oldrow['td_internal_energy_eV']-Ebo;Xnew=newrow['td_internal_energy_eV']-Ebo
        close(Xold,static['excitation_internal_energy_eV'],'Original vacancy excitation')
        density_file=capsule.bound_file(new['density'],REFINED+'trajectory')
        with np.load(capsule.root/density_file,allow_pickle=False) as data:
            density=data['density_e_A3'];cell=np.asarray(ground['cell_A'])
            require(density.shape==tuple(2*np.asarray(ground['realspace_shape'])),'Wrong refined density grid')
            close(data['cell_A'],cell,'Refined density cell')
            density_error=max(density_error,close(density_potential(density,cell,newrow['projectile_A']),newrow['electron_potential_eV'],'Actual refined density potential'))
            close(float(np.sum(density)*np.linalg.det(cell)/density.size),newrow['integrated_all_electrons'],'Refined density electron count')
        point_rows.append(dict(R_A=new['R_A'],path_length_A=start_R-new['R_A'],production_excitation_eV=Xold,refined_excitation_eV=Xnew))
    first=point_rows[0]
    for point in point_rows:
        length=point['path_length_A']
        point['production_vacancy_mean_loss_eV_A']=(point['production_excitation_eV']-first['production_excitation_eV'])/length if length else None
        point['refined_vacancy_mean_loss_eV_A']=(point['refined_excitation_eV']-first['refined_excitation_eV'])/length if length else None
    final=point_rows[-1];old=final['production_vacancy_mean_loss_eV_A'];new=final['refined_vacancy_mean_loss_eV_A']
    require(old>0 and final['R_A']==.1,'Final positive production drag reference missing')
    relative=abs(new-old)/old
    return dict(schema='pdd-vacancy-halfdt-control-v1',authority='Axiom2315/2345',rows=point_rows,
        full_refined_rows=len(rows),production_mean_loss_eV_A=old,refined_mean_loss_eV_A=new,relative_change=relative,
        acceptance_relative=.1,passed=bool(relative<=.1),max_particle_drift_electrons=drift,max_projector_difference=op,
        max_force_work_energy_residual_eV=float(max(abs(closure))),recomputed_force_work_difference_eV=work_difference,
        max_density_potential_replay_difference_eV=density_error,initial_GPW=ground['ground_gpw'],
        encounter_excitation_dt_checked=False,scope='Only undisplaced2keV head-on vacancy finite-path mean drag is dt-checked. Original BO references reused at verified identical positions; no ground/BO SCF, occupied rerun or mixed-dt encounter inference.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('capsule',type=Path);args=parser.parse_args()
    try:result=read_control(args.capsule)
    except (ValueError,KeyError,TypeError) as err:parser.exit(2,str(err)+'\n')
    print(json.dumps(result,indent=2,allow_nan=False))
