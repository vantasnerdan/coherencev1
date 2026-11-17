"""Keep original estimates; change only the Pd material inputs actually supplied."""
import math


def validate_scenarios(scenarios):
    if not scenarios or len({s['id'] for s in scenarios}) != len(scenarios):
        raise ValueError('Distinct explicit parameter scenarios are required')
    for scenario in scenarios:
        kind=scenario['kind'];overrides=scenario['material_overrides']
        if kind not in ['baseline_estimates','conditional_0K_electronic']:
            raise ValueError('Unknown scientific parameter role')
        if set(overrides)-{'Eb','Eads'}:
            raise ValueError('DFT barrier/prefactor or capacity cannot silently replace sticking/mtrap')
        if not all(math.isfinite(float(v)) for v in overrides.values()):
            raise ValueError('Nonfinite material parameter')
        if kind=='baseline_estimates' and overrides:
            raise ValueError('Baseline estimate scenario cannot hide new parameters')
        if kind!='baseline_estimates':
            if not scenario.get('missing_finite_T_corrections'):
                raise ValueError('Electronic-only rows must label missing finite-T corrections')
            references=scenario['references']
            if 'Eb' in overrides and references.get('Eb')!='matched_bulk_interstitial_D':
                raise ValueError('Vacancy binding must use the matched bulk-interstitial reference')
            if 'Eads' in overrides and references.get('Eads')!='differential_adsorption_half_H2_0K':
                raise ValueError('Adsorption must use the differential half-H2 electronic reference')


def build_jobs(runner,scenarios,pdrh):
    validate_scenarios(scenarios)
    original=runner.knob_jobs();jobs=[]
    for scenario in scenarios:
        fresh=[] if scenario['kind']=='baseline_estimates' else [j for j in original if j[1] in ['Pd','PdD']]
        fresh += [(f'knob_{condition}_{energy}_Pd80Rh20 conditional alloy',condition,
                   dict(E_cm=energy,alloy=pdrh['material_overrides']))
                  for condition in ['Pd','PdD'] for energy in [.25,2.5]]
        for tag,condition,case in fresh:
            case=dict(case);alloy=case.get('alloy')
            material=dict(runner.T.ALLOYS[alloy]) if isinstance(alloy,str) else dict(alloy or {})
            material.update(scenario['material_overrides'])
            if material:case['alloy']=material
            jobs.append(dict(id=scenario['id']+'__'+tag,tag=tag,condition=condition,
                knob=tag.split('_',3)[3],case=case,scenario=scenario,
                transferred_to_alloy=bool(alloy and scenario['material_overrides'])))
    return jobs
