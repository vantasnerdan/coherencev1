"""Join actual P255 quartet/work outputs without claiming isolated material admission."""
import argparse
import json
from pathlib import Path
import numpy as np
from ase.units import Bohr
from screen_io import manifest,require,write_json
from substrate_framework.collision_bare_reference import PeriodicNuclearGreen
from substrate_framework.collision_observables import COULOMB_EV_A


def load_case(run):
    calculation=run/'calculation'
    ground=json.loads((calculation/'ground/ground.json').read_bytes())
    trajectory_path=calculation/'trajectory/trajectory.json'
    trajectory=json.loads(trajectory_path.read_bytes())
    bo_path=calculation/'bo/excitation.json';bo=json.loads(bo_path.read_bytes())
    require(manifest(trajectory_path)==bo['trajectory'],'BO is bound to another trajectory')
    require(trajectory['specification']==bo['specification'],'BO/trajectory input epochs differ')
    host=bo['neutral_host_reference']
    require(host['charge']==0,'Quartet zero requires the actual neutral host')
    return dict(ground=ground,trajectory=trajectory,bo=bo,host=host,
                bindings=dict(ground=manifest(calculation/'ground/ground.json'),
                              trajectory=manifest(trajectory_path),bo=manifest(bo_path)))


def match_cases(occupied,vacant):
    a=occupied['trajectory']['specification'];b=vacant['trajectory']['specification']
    require(a['occupied_target'] and not b['occupied_target'],'Require occupied then vacant roles')
    require({k:v for k,v in a.items() if k!='occupied_target'}==
            {k:v for k,v in b.items() if k!='occupied_target'},'Physical pair settings differ')
    require(np.array_equal(occupied['host']['cell_A'],vacant['host']['cell_A']),'Pair cells differ')
    target=np.asarray(occupied['ground']['path']['target_A'])
    host=occupied['host'];keep=[];removed=0
    for i,(symbol,position) in enumerate(zip(host['symbols'],host['positions_A'])):
        if symbol=='H' and np.linalg.norm(np.asarray(position)-target)<1e-10:
            removed+=1
        else:keep.append(i)
    require(removed==1,'Occupied host lacks its unique target')
    require([host['symbols'][i] for i in keep]==vacant['host']['symbols'] and
            np.allclose(np.asarray(host['positions_A'])[keep],vacant['host']['positions_A'],rtol=0,atol=1e-10),
            'Host zero does not remove only the neutral target')
    require(abs(host['electrons']-vacant['host']['electrons']-1)<1e-7,'Target electron count is not neutral')


def pair_rows(occupied,vacant,bare):
    cell=np.asarray(occupied['host']['cell_A']);green=PeriodicNuclearGreen(cell)
    alpha_a,alpha_b=bare['native_alpha1_Bohr_inverse2']
    width=np.sqrt(alpha_a*alpha_b/(alpha_a+alpha_b))/Bohr
    target=np.asarray(occupied['ground']['path']['target_A'])
    rows=[]
    require(len(occupied['bo']['checkpoints'])==len(vacant['bo']['checkpoints']),'Different BO support')
    for o,v in zip(occupied['bo']['checkpoints'],vacant['bo']['checkpoints']):
        require(o['R_A']==v['R_A'],'Unmatched checkpoint radius')
        to=occupied['trajectory']['rows'][o['row_index']]
        tv=vacant['trajectory']['rows'][v['row_index']]
        point=np.asarray(to['projectile_A'])
        require(np.allclose(point,tv['projectile_A'],rtol=0,atol=1e-11),'Unmatched projectile positions')
        R=float(o['R_A']);require(R>0 and abs(np.linalg.norm(point-target)-R)<1e-9,'Wrong head-on radius')
        raw=(o['bo']['internal_energy_eV']-v['bo']['internal_energy_eV']
             -occupied['host']['energy']['internal_energy_eV']+vacant['host']['energy']['internal_energy_eV'])
        correction=green.point_minus_gaussian_eV(point-target,width)
        total=raw+correction
        rows.append(dict(R_A=R,BO_quartet_soft_eV=raw,bare_periodic_point_minus_soft_eV=correction,
            V_pair_periodic_point_eV=total,candidate_lowering_vs_isolated_C_over_R_eV=COULOMB_EV_A/R-total,
            signed_electronic_candidate_eV=total-COULOMB_EV_A/R,
            U_phi_diagnostic_eV=-(to['electron_potential_eV']-tv['electron_potential_eV']),
            occupied_excitation_eV=o['excitation_internal_energy_eV'],vacant_excitation_eV=v['excitation_internal_energy_eV'],
            encounter_external_work_raw_eV=o['external_work_eV']-v['external_work_eV']))
    first=rows[0]
    initial_encounter=first['occupied_excitation_eV']-first['vacant_excitation_eV']
    for row in rows:
        row['path_length_A']=first['R_A']-row['R_A']
        row['encounter_excitation_change_eV']=row['occupied_excitation_eV']-row['vacant_excitation_eV']-initial_encounter
        row['vacancy_loss_change_eV']=row['vacant_excitation_eV']-first['vacant_excitation_eV']
        row['vacancy_mean_loss_eV_A']=row['vacancy_loss_change_eV']/row['path_length_A'] if row['path_length_A']>0 else None
        row['encounter_external_work_mapped_increment_eV']=(row['encounter_external_work_raw_eV']
            +row['bare_periodic_point_minus_soft_eV']-first['bare_periodic_point_minus_soft_eV'])
        row['anchored_encounter_work_profile_eV']=first['V_pair_periodic_point_eV']+row['encounter_external_work_mapped_increment_eV']
        row['independent_force_work_identity_residual_eV']=(row['anchored_encounter_work_profile_eV']
            -row['V_pair_periodic_point_eV']-row['encounter_excitation_change_eV'])
    return rows


def analyze(occupied_run,vacant_run,bare_path):
    occupied=load_case(occupied_run);vacant=load_case(vacant_run);match_cases(occupied,vacant)
    bare=json.loads(bare_path.read_bytes());require(bare['passed'],'Unqualified bare operator')
    with np.load(bare_path.with_suffix('.npz'),allow_pickle=False) as data:
        require(np.allclose(data['cell_A'],occupied['host']['cell_A'],rtol=0,atol=1e-12),'Bare qualification belongs to another cell')
        require(np.array_equal(data['rho_target_Bohr_inv3'].shape,2*np.asarray(occupied['ground']['realspace_shape'])),
                'Bare qualification belongs to another grid')
    rows=pair_rows(occupied,vacant,bare)
    residual=max(abs(r['independent_force_work_identity_residual_eV']) for r in rows)
    require(residual<.5,'Actual matched force-work identity fails')
    last=rows[-1];lowering=last['candidate_lowering_vs_isolated_C_over_R_eV']
    return dict(schema='pdd-periodic-point-pair-diagnostic-v1',observable='periodic_target_conditioned_BO_interaction',
        rows=rows,source_bindings=dict(occupied=occupied['bindings'],vacant=vacant['bindings'],bare=manifest(bare_path)),
        host_energy_zero='E_probe,o-E_probe,v-E_host,o+E_host,v; actual neutral same-settings host states',
        bare_mapping='Zero-G periodic point minus qualified Gaussian compensation kernel; not an isolated material-image correction',
        material_pair_export_eligible=False,remaining='Neutral-image and h/dt/stopping qualifications; no prepared state or rate is supplied',
        minimum_radius_nominal_encounter_ratio=abs(last['encounter_excitation_change_eV']/lowering) if lowering>0 else None,
        extension_trigger_status='Provisional current-cell/h diagnostic only; no extra energy/direction is authorized',
        max_independent_force_work_identity_residual_eV=residual,
        electronic_state_convention='Internal-energy excess to independently equilibrated BO at Fermi width0.1eV; signed, no clipping')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('occupied','vacant','bare','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    write_json(args.output,analyze(args.occupied,args.vacant,args.bare))
