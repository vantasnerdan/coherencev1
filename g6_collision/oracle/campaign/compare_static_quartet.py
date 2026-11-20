"""Compare only the approved R0.1A,h0.20->0.16 neutral-host quartet."""
import argparse
import json
from pathlib import Path
import numpy as np
from screen_io import require,manifest,write_json
from collision_pair_readout import load_case,match_cases
from static_state_data import load_state




def match_state(record,case,state):
    require(record['actual_kpoints']==case['ground']['actual_kpoints'],'Control k support changed')
    require(abs(record['electrons']-case['ground']['electrons'])<1e-7,'Control electron count changed')
    require(np.array_equal(record['cell_A'],case['ground']['cell_A']),'Control cell changed')
    for name in ('repeat','kmesh','D_displacement_A','direction','impact_A','relaxed_ground_sha256'):
        require(record['specification'][name]==case['trajectory']['specification'][name],
                'Control changed '+name+' as well as h')
    if state=='host':
        symbols=case['host']['symbols'];positions=case['host']['positions_A']
    else:
        symbols=case['ground']['symbols']
        points=[p for p in case['trajectory']['checkpoints'] if abs(p['R_A']-.1)<1e-10]
        require(len(points)==1,'Baseline has no unique actualR0.1 checkpoint')
        positions=points[0]['positions_A']
    require(record['symbols']==symbols and np.allclose(record['positions_A'],positions,rtol=0,atol=1e-10),
            'Control geometry is not the actual matched frozen state')
    require(np.prod(record['realspace_shape'])>np.prod(case['ground']['realspace_shape']),
            'The control grid is not finer')


def compare(occupied_run,vacant_run,control,reference):
    cases=dict(occupied=load_case(occupied_run),vacant=load_case(vacant_run))
    require(all(c['trajectory']['specification']['h_A']==.2 for c in cases.values()),'Wrong coarse h reference')
    match_cases(cases['occupied'],cases['vacant'])
    states={};bindings={}
    for role,case in cases.items():
        for state in ('probe','host'):
            key=role+'-'+state
            states[key],bindings[key]=load_state(control,role,state,h_A=.16,repeat=[1,1,1])
            match_state(states[key],case,state)
        require(abs(states[role+'-probe']['electrons']-states[role+'-host']['electrons'])<1e-7,
                'Neutral host/probe electron-number convention changed')
    for name in ('collision_ground.py','collision_geometry.py'):
        require(manifest(control/'source'/name)==manifest(occupied_run/'source'/name),
                'Control changed the actual factory/geometry source')
    coarse=[]
    for role in ('occupied','vacant'):
        rows=[r for r in cases[role]['bo']['checkpoints'] if r['R_A']==.1]
        require(len(rows)==1,'Wrong baseline BO support')
        coarse.append(rows[0]['bo']['internal_energy_eV']-cases[role]['host']['energy']['internal_energy_eV'])
    raw_coarse=coarse[0]-coarse[1]
    energy=lambda role,state:states[role+'-'+state]['energy']['internal_energy_eV']
    raw_fine=energy('occupied','probe')-energy('vacant','probe')-energy('occupied','host')+energy('vacant','host')
    original=json.loads(reference.read_bytes())
    require(original['source_bindings']['occupied']==cases['occupied']['bindings'] and
            original['source_bindings']['vacant']==cases['vacant']['bindings'],'Pair reference source epoch differs')
    ref=[r for r in original['rows'] if r['R_A']==.1]
    require(len(ref)==1 and abs(ref[0]['BO_quartet_soft_eV']-raw_coarse)<1e-10,
            'Pair reference is not the actual baseline quartet')
    delta=raw_fine-raw_coarse
    lowering_coarse=ref[0]['candidate_lowering_vs_isolated_C_over_R_eV']
    lowering_fine=lowering_coarse-delta
    tolerance=max(5.,.1*abs(lowering_fine))
    return dict(schema='pdd-R0p1-static-h-comparison-v1',authority='Axiom2250',
        raw_quartet_h020_eV=raw_coarse,raw_quartet_h016_eV=raw_fine,
        delta_V_pair_h016_minus_h020_eV=delta,
        candidate_positive_lowering_h020_eV=lowering_coarse,candidate_positive_lowering_h016_eV=lowering_fine,
        delta_positive_lowering_eV=-delta,delta_signed_electronic_correction_eV=delta,
        tolerance_eV=tolerance,passed=bool(abs(delta)<=tolerance),
        fine_grid=states['occupied-probe']['realspace_shape'],coarse_grid=cases['occupied']['ground']['realspace_shape'],
        state_bindings=bindings,baseline_pair_reference=manifest(reference),
        baseline_bindings={role:case['bindings'] for role,case in cases.items()},
        scope='Fixed cell/k/phase,R0.1A neutral-host quartet; common bare correction cancels in Delta_h. Image/dt/drag qualifications remain distinct. No h-TD family or isolated pair export is authorized by this file.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('occupied','vacant','control','reference','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    write_json(args.output,compare(args.occupied,args.vacant,args.control,args.reference))
