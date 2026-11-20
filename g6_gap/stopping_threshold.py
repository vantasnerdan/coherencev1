"""Read operational half-ICRU stopping windows from saved P257 curves.

No Native computation or fitted normalization. The source ICRU49/Geant4
law is velocity-proportional throughout the requested 0.1--13keV interval;
the inherited forward sqrt law is a distinct comparator. Finite-eta/G
crossings are not exact zero-stopping thresholds or a convergence proof.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from screen_io import manifest,require,write_json
from substrate_framework.lithium_fuel_return import (
    lit_icru_electronic,lit_inherited_sqrt_electronic,stopping_per_formula)


def _half_window(energy,ratio):
    below=ratio<=.5
    crossings=[]
    for i in np.flatnonzero(below[:-1]!=below[1:]):
        left,right=float(energy[i]),float(energy[i+1])
        if ratio[i]==.5:
            right=left
        elif ratio[i+1]==.5:
            left=right
        crossings.append(dict(energy_lab_bracket_keV=[left,right],
                              direction='leaving_suppression' if below[i] else 'entering_suppression'))
    if not below[0]:
        status='not_half_suppressed_at_resolved_lower_edge'
        boundary=None
    elif not crossings:
        status='half_suppression_persists_through_resolved_upper_edge'
        boundary=None
    else:
        status='resolved_low_energy_window_upper_boundary'
        boundary=crossings[0]['energy_lab_bracket_keV']
    return dict(status=status,low_energy_window_upper_boundary_keV=boundary,
                resolved_domain_lab_keV=[float(energy[0]),float(energy[-1])],
                ratio_at_domain_edges=[float(ratio[0]),float(ratio[-1])],
                all_crossings=crossings,
                unique_crossing=bool(len(crossings)==1))


def threshold_record(record):
    energy=np.asarray(record['energy_lab_keV'],dtype=float)
    stopping=np.asarray(record['electronic_cross_section_per_formula'],dtype=float)
    direct=np.asarray(record['stopping_eV_A'],dtype=float)
    require(energy.ndim==1 and energy.size>=2 and np.isfinite(energy).all()
            and np.all(np.diff(energy)>0),'Increasing finite energy grid required')
    require(stopping.shape==energy.shape and direct.shape==energy.shape,
            'Matched stopping and energy arrays required')
    require(np.isfinite(stopping).all() and np.all(stopping>=0),'Passive finite stopping required')
    converted=stopping_per_formula(direct,record['volume_A3'])
    require(np.allclose(converted,stopping,rtol=1e-12,atol=1e-14),
            'Saved cross section disagrees with eV/A and formula volume')
    require(np.any(energy==.1) and np.any(energy==13.),'Explicit0.1/13keV domain edges required')
    keep=(energy>=.1)&(energy<=13)
    e,s=energy[keep],stopping[keep]
    references={}
    for name,function in [('source_ICRU49_velocity_proportional',lit_icru_electronic),
                          ('inherited_forward_sqrt',lit_inherited_sqrt_electronic)]:
        baseline=function(e)
        ratio=s/baseline
        references[name]=dict(**_half_window(e,ratio),
            ratio_to_reference=ratio.tolist(),reference_cross_section_per_formula=baseline.tolist())
    return dict(energy_lab_keV=e.tolist(),references=references,
                factor=.5,units='eV/(1e15 formula/cm2)',
                inference='Sampled brackets, not interpolated exact roots; all reentrant crossings retained. A settings envelope and model omissions must be assessed separately.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stopping',type=Path,nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    results=[]
    for path in args.stopping:
        record=json.loads(path.read_bytes())
        result=threshold_record(record)
        results.append(dict(dataset=path.parent.name+'/'+path.name,
                            configuration=Path(record['group']).name,
                            source=manifest(path),**result))
    write_json(args.output,dict(schema='P257-half-stopping-threshold-v1',
        source=manifest(__file__),records=results,
        eligibility='Operational saved-curve readout only. No strict zero-threshold or full settings-convergence certificate.'))


if __name__=='__main__':
    main()
