"""Combine observed P257 settings changes and bound the half-stopping window.

This is an L1 envelope of executed controls, not a probabilistic uncertainty
or an infinite-basis proof. No interpolation, fitted bridge or Native run.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from screen_io import manifest,require,write_json
from stopping_threshold import threshold_record,_half_window
from substrate_framework.lithium_fuel_return import lit_icru_electronic


def load_stopping(path,points,cutoff):
    data=json.loads(path.read_bytes())
    quadrature=Path(data['group'])/'quadrature.json'
    require(manifest(quadrature)==data['quadrature'],'Changed actual quadrature bytes')
    q=json.loads(quadrature.read_bytes())
    require(len(q['points'])==points and len(data['source_points'])==points,
            'Wrong complete q-star grid for the named control')
    require(all(p['ecut_eV']==cutoff and p['bands']==96 and p['eta_eV']==.025
                for p in q['points']),'Wrong G/band/eta settings for the named control')
    require(all(np.array_equal(p['energy_eV'],q['points'][0]['energy_eV']) for p in q['points']),
            'Inconsistent per-point frequency grids')
    frequency=np.asarray(q['points'][0]['energy_eV'])
    require(frequency[0]==0 and frequency[-1]==80 and np.all(np.diff(frequency)==.03125),
            'Wrong matched fine frequency grid')
    require(abs(sum(p['q_weight'] for p in q['points'])-1)<1e-12,'Incomplete q weight')
    threshold_record(data)  # Actual saved-output units/domain refusal path.
    return data,q,dict(stopping=manifest(path),quadrature=manifest(quadrature))


def analyze(prior_path,q2_path,g100_path,q4_path):
    prior=json.loads(prior_path.read_bytes())
    base,bq,bb=load_stopping(q2_path,2,200)
    coarse,cq,cb=load_stopping(g100_path,2,100)
    final,fq,fb=load_stopping(q4_path,10,200)
    energy=np.asarray(final['energy_lab_keV'])
    require(all(np.array_equal(energy,x['energy_lab_keV']) for x in (prior,base,coarse)),
            'Controls must use the same actual energy grid')
    require(bq['ground']==cq['ground']==fq['ground'],'Control ground-state bytes differ')
    require(base['volume_A3']==coarse['volume_A3']==final['volume_A3'],'Control formula volumes differ')
    old=np.asarray(prior['partial_L1_control_envelope_eV_A'])
    require(old.shape==energy.shape and np.isfinite(old).all() and np.all(old>=0),
            'Invalid prior observed-control envelope')
    stopping=np.asarray(final['stopping_eV_A'])
    G_delta=np.asarray(base['stopping_eV_A'])-np.asarray(coarse['stopping_eV_A'])
    q_delta=stopping-np.asarray(base['stopping_eV_A'])
    error=old+abs(G_delta)+abs(q_delta)
    tolerance=np.maximum(.15*abs(stopping),.01)
    domain=(energy>=.1)&(energy<=13)
    ratios=error/tolerance
    selected=np.flatnonzero(domain)
    worst=int(selected[np.argmax(ratios[domain])])
    reference=lit_icru_electronic(energy[domain])*10/final['volume_A3']
    thresholds=dict(central=threshold_record(final)['references']['source_ICRU49_velocity_proportional'],
        guaranteed_half_suppression_within_observed_envelope=_half_window(energy[domain],(stopping[domain]+error[domain])/reference),
        possible_half_suppression_within_observed_envelope=_half_window(energy[domain],(stopping[domain]-error[domain])/reference))
    return dict(schema='P257-refined-control-analysis-v1',
        source=manifest(__file__),inputs=dict(prior=manifest(prior_path),q2_G200=bb,q2_G100=cb,q4_G200=fb),
        energy_lab_keV=energy.tolist(),q2_G100_to_G200_delta_eV_A=G_delta.tolist(),
        q2_to_q4_delta_eV_A=q_delta.tolist(),combined_observed_L1_envelope_eV_A=error.tolist(),
        tolerance_eV_A=tolerance.tolist(),max_required_domain_tolerance_ratio=float(ratios[worst]),
        worst_energy_lab_keV=float(energy[worst]),
        observed_15percent_absolute_floor_controls_pass=bool(ratios[worst]<=1),
        threshold=thresholds,
        scope='Observed executed-settings L1 envelope only; not a rigorous infinite-q/G bound, strict zero threshold, model-error allowance or physical energy-return verdict.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('prior','q2','g100','q4','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    write_json(args.output,analyze(args.prior,args.q2,args.g100,args.q4))


if __name__=='__main__':
    main()
