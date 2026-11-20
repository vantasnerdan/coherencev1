"""Export BO-only material curves after the explicit P255 h/cell gates."""
import argparse
import json
from pathlib import Path
from screen_io import require,manifest,write_json


def export(baseline,gamma,h_check,image_check,displacement,output):
    base=json.loads(baseline.read_bytes());moved=json.loads(gamma.read_bytes())
    h=json.loads(h_check.read_bytes());image=json.loads(image_check.read_bytes())
    contrast=json.loads(displacement.read_bytes())
    require(h['passed'] and image['passed'] and contrast['passed'],'Unqualified BO h/cell/geometry inputs')
    require(h['baseline_pair_reference']==image['baseline']==contrast['baseline_pair']==manifest(baseline),
            'BO qualification belongs to another baseline')
    require(h['baseline_bindings']['occupied']==base['source_bindings']['occupied'] and
            h['baseline_bindings']['vacant']==base['source_bindings']['vacant'],'BO h source epoch differs')
    require(contrast['gamma_pair']==manifest(gamma),'Displacement admission belongs to another pair')
    require(image['authority']=='Axiom2291' and image['R_A']==.1,'Wrong approved representative image check')
    require(base['schema']==moved['schema']=='pdd-periodic-point-pair-diagnostic-v1','Wrong BO source observable')
    output.mkdir(parents=True,exist_ok=False)
    results={}
    for name,source,path,shift in (('baseline',base,baseline,0.),('gamma',moved,gamma,.0768)):
        rows=sorted(source['rows'],key=lambda row:row['R_A'])
        reference=dict(
            source='P255/issue244 actual neutral-host BO quartet; original paired output SHA256 '+manifest(path)['sha256'],
            configuration=f'Frozen ordered PdD conventional cell; PBE,Pd10,all-electron regularizedH,h0.20A actual20^3,fullGamma-centred k2,Fermi width0.1eV. AllD[100] shift={shift}A,Pd fixed.2keV head-on[110] checkpoint geometry; not a prepared locked state.',
            energy_zero='E_probe,occupied-E_probe,vacant-E_host,occupied+E_host,vacant using actual neutral host-only states; no arbitrary far-point offset.',
            periodic_to_isolated='Native alpha1=10Bohr^-2 nuclear compensation replaced by same-zero-G periodic point-minus-Gaussian Green block. Neutral material image control:2x2x2,matchedGamma k1/240bands,R0.1,lowering40.54979961597->40.79863409606eV. h0.20->0.16 changes pair by0.45104165082eV. Both pass max(5eV,10%); Axiom2291/2314 admits original cell/h for both complete curves. Representative-point qualification, not an exact infinite-cell limit; no bare-only image subtraction or fitted shift.',
            bare='isolated_unit_charge_coulomb')
        data=dict(schema='coherencev1-target-conditioned-pair-v1',observable='target_conditioned_BO_total_interaction',
            R_A=[row['R_A'] for row in rows],V_pair_total_eV=[row['V_pair_periodic_point_eV'] for row in rows],reference=reference,
            qualification=dict(h=manifest(h_check),cell=manifest(image_check),actual_displacement=manifest(displacement),authority=['Axiom2291','Axiom2314']),
            source_bindings=source['source_bindings'],support='Only supplied0.1–1.6A samples; no U(0),prepared state,rate,core continuation or far-tail completion.',
            excluded='Encounter work/excitation and vacancy drag are separate trajectory outputs, not additions to BO. Encounter excitation is not dt-checked (Axiom2315).')
        dest=output/(name+'_BO_pair.json');write_json(dest,data);results[name]=manifest(dest)
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('baseline','gamma','h_check','image_check','displacement','output'):
        parser.add_argument('--'+name.replace('_','-'),type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(export(args.baseline,args.gamma,args.h_check,args.image_check,args.displacement,args.output),indent=2))
