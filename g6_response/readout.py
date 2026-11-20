"""Consume complete actual-frequency q families and expose numerical controls."""
import argparse
import json
from pathlib import Path

import numpy as np
from _upstream.dielectric_observables import metallic_screening_length,powder_stopping_q
from screen_io import manifest,require,write_json


def evaluate(directory):
    directory=Path(directory)
    completed=json.loads((directory/'completed.json').read_text())
    cfg=json.loads((directory/'input.json').read_text())
    anchor=next(iter(completed['sources'].values()))['ground']
    expected={'Pd':['Pd'],'PdD':['Pd','H'],'TiD2':['Ti','H','H']}[cfg['host']]
    for record in completed['sources'].values():
        ground=record['ground']
        require(ground['geometry']['symbols']==expected,'Wrong physical composition')
        require(ground['neutral_electrons']==cfg['electrons'],'Wrong PAW electron count')
        for name in ['cell_A','positions_A']:
            require(np.allclose(ground['geometry'][name],anchor['geometry'][name],rtol=0,atol=1e-8),
                    'Numerical controls mix physical geometries')
    groups={}
    heads={}
    curves=[]
    energy=np.array([.5,1,2,3,5,8,13])
    velocity=np.sqrt(2*energy*1000/completed['hartree_eV']/3670.48296788)
    C=completed['coulomb_eV_A']
    for p in completed['points']:
        settings=p['settings']
        path=directory/p['arrays']['path']
        require(manifest(path)==p['arrays'],'Saved matrix bytes differ')
        ground=completed['sources'][settings['bank']]['ground']
        volume=ground['cell_volume_A3']
        group=settings['group']
        with np.load(path,allow_pickle=False) as data:
            Q=data['Q_Ainv'];omega=data['energy_eV']
            inverse=data['inverse_epsilon_coulomb_symmetrized']
            loss=-np.diagonal(inverse,axis1=1,axis2=2).imag
            require(float(loss.min())>=-1e-8*max(1.,float(abs(loss).max())),
                    'Saved loss violates passive equilibrium response; no clipping is permitted')
            eps=1/inverse[0,0,0]
            curves.append(dict(group=group,q_fractional=settings['q_fractional'],
                Q_head_Ainv=float(np.linalg.norm(Q[0])),epsilon_M_real=float(eps.real),
                epsilon_M_imag=float(eps.imag),response_cutoff_eV=settings['ecut_eV'],
                nbands=settings['nbands'],eta_eV=settings['eta_eV'],
                kmesh=ground['ground_state_settings']['kmesh']))
            if group.startswith('lambda'):
                heads.setdefault(group,[]).append((float(np.linalg.norm(Q[0])),1/inverse[0,0,0]))
            else:
                value=powder_stopping_q(Q,omega,loss,velocity*C,volume=volume,
                                       weight=settings['weight'],coulomb=C)
                entry=groups.setdefault(group,dict(S=np.zeros(len(energy)),weight=0))
                entry['S']+=value
                entry['weight']+=settings['weight']
    for group,row in groups.items():
        require(abs(row['weight']-1)<1e-12,'Incomplete normalized BZ family '+group)
    primary=groups['primary-N3']['S']
    n2=groups['primary-N2']['S']
    errors={}
    differences=[]
    for group,row in groups.items():
        if group=='primary-N3':continue
        reference=primary if group=='primary-N2' else n2
        errors[group]=float(np.max(np.abs(row['S']-reference)/np.maximum(np.abs(reference),1e-12)))
        differences.append(abs(row['S']-reference))
    lambdas={}
    for name,rows in heads.items():
        rows.sort()
        q,eps=zip(*rows)
        lambdas[name]=metallic_screening_length(q,eps,coulomb=C)
    lam=lambdas['lambda']['lambda_A']
    screening_errors={name:abs(row['lambda_A']/lam-1) for name,row in lambdas.items() if name!='lambda'}
    coupled_stopping=np.sum(differences,axis=0)/np.maximum(abs(primary),1e-12)
    U=C/lam
    window=max(abs(C/value-U) for value in lambdas['lambda']['q_window_lambda_span_A'])
    screening_absolute=window+sum(abs(C/row['lambda_A']-U) for name,row in lambdas.items() if name!='lambda')
    screening_ok=(all(row['q_window_pass_10percent'] for row in lambdas.values())
                  and all(v<=.1 for v in screening_errors.values())
                  and screening_absolute<=max(5.,.1*U))
    stopping_ok=bool(np.max(coupled_stopping)<=.15)
    first=next(iter(completed['sources'].values()))['ground']
    volume=first['cell_volume_A3']
    symbols=first['geometry']['symbols']
    metals=sum(s!='H' for s in symbols)
    H=symbols.count('H')/metals
    number=metals/volume
    # g4_beam.forward.model Se_atom units:1e15 atoms/cm^2, physical density converts by10*n_A^-3.
    coefficients={'Pd':5.238,'PdD':5.238+H*1.254,'TiD2':4.858+H*1.254}
    benchmark=10*number*coefficients[cfg['host']]*np.sqrt(energy*1.007276/2.01410)
    Aeff=primary/(10*number*np.sqrt(energy*1.007276/2.01410))
    historical={'Pd':.39747864274672146,'PdD':.541,'TiD2':.435}[cfg['host']]
    return dict(schema='g6-lowq-controlled-observables-v1',host=cfg['host'],model=cfg['model'],
        energy_lab_keV=energy.tolist(),velocity_au=velocity.tolist(),stopping_eV_A=primary.tolist(),
        friction_eV_A_per_vau=(primary/velocity).tolist(),
        lowq_lambda_A=lam,lowq_Yukawa_proxy_eV=U,
        lowq_proxy_over_stopping_A=(U/primary).tolist(),
        lowq_proxy_over_friction_A_vau=(U/(primary/velocity)).tolist(),
        stopping_controls_relative=errors,screening_controls_relative=screening_errors,
        screening_details=lambdas,
        coupled_stopping_relative_estimate=coupled_stopping.tolist(),
        coupled_screening_absolute_estimate_eV=float(screening_absolute),
        error_estimate_scope='Conservative sum of observed two-point absolute changes; empirical numerical estimate, not rigorous enclosure.',
        protocol_diagnostics_unmeasured_for_host=[name for name in ['eta','frequency'] if name not in groups],
        numerical_eligibility_scope='Adopted low-Q protocol: k/q/band/G each host, eta/frequency diagnosed on Pd only. Unmeasured host-specific terms are not asserted zero.',
        stopping_groups={name:row['S'].tolist() for name,row in groups.items()},
        epsilon_macro_static_curve=curves,
        friction_definition='Effective S(v)/v at each stated speed; no fitted coefficient or strict v->0 extrapolation.',
        ICRU=dict(expected_formula_unit_A1=coefficients[cfg['host']],effective_A1=Aeff.tolist(),
                  host_A1=coefficients[cfg['host']]-H*1.254,hydrogen_A1=1.254,H_per_metal=H,
                  metal_number_density_Ainv3=number,
                  scope='Untuned Bragg comparison at actual source stoichiometry/density, not transported foil loading',
                  stopping_eV_A=benchmark.tolist(),relative_difference=(primary/benchmark-1).tolist(),
                  within30percent=bool(np.max(abs(primary/benchmark-1))<=.3),tuned=False),
        historical_lowq_screening=dict(
            reference_lambda_A=historical,relative_difference=lam/historical-1,
            source='Pinned Pd-length control; PdD/TiD2 rounded screening README values; comparison is not fitted'),
        eligible=bool(screening_ok and stopping_ok),screening_qualified=bool(screening_ok),stopping_qualified=bool(stopping_ok),
        scope='Authorized low-Q finite-G linear RPA only; C/lambda is an explicitly low-Q Yukawa merit numerator, not short-range screening.',
        omission='FrozenPd4s/4p response is omitted, not assumed zero. No highQ25/Angstrom screening result.',
        source=manifest(directory/'completed.json'),readout_source=manifest(__file__))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--directory',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args()
    record=evaluate(args.directory)
    write_json(args.output,record)
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    main()
