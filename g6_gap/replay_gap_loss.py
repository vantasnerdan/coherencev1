"""Replay the complete P257 delivery with NumPy; no GPAW/ASE/SciPy required."""
import argparse
import copy
import json
import os
from pathlib import Path
import numpy as np
from screen_io import manifest,require,write_json
from substrate_framework.dielectric_observables import powder_stopping_q
from substrate_framework.lithium_fuel_return import stopping_per_formula
from gap_control_analysis import analyze
from gap_delivery_tables import build_tables


def checked_json(entry):
    path=Path(entry['path'])
    require(manifest(path)==entry['binding'],'Changed portable metadata: '+entry['path'])
    return json.loads(path.read_bytes())


def replay_dataset(entry,constants,output):
    expected=checked_json(entry['stopping'])
    energy=np.asarray(expected['energy_lab_keV'])
    hartree=constants['Hartree_eV'];bohr=constants['Bohr_A']
    velocity=np.sqrt(2*energy*1000/hartree/constants['velocity_mass_D_over_e'])
    curves={stride:np.zeros_like(energy) for stride in (1,2,4)}
    total_weight=0.
    for point in expected['source_points']:
        path=Path(point['compact_path'])
        require(manifest(path)==point['compact'],'Changed saved loss bytes: '+str(path))
        with np.load(path,allow_pickle=False) as data:
            w=data['energy_eV'];qg=data['Q_Ainv'];loss=data['loss_wG'];weight=float(data['q_weight'])
        require(loss.min()>=-1e-9,'Nonpassive saved loss')
        total_weight+=weight
        for stride,stopping in curves.items():
            take=np.unique(np.r_[np.arange(0,len(w),stride),len(w)-1])
            stopping+=powder_stopping_q(qg,w[take],loss[take],velocity*hartree*bohr,
                volume=expected['volume_A3'],weight=weight,coulomb=hartree*bohr)
    require(abs(total_weight-1)<1e-12,'Incomplete saved-loss q weights')
    error=max(float(np.max(abs(curves[s]-np.asarray(expected['frequency_decimation_stopping_eV_A'][str(s)])))) for s in curves)
    require(error<1e-11,'Saved loss does not reproduce the original stopping')
    record=copy.deepcopy(expected)
    record['stopping_eV_A']=curves[1].tolist()
    record['frequency_decimation_stopping_eV_A']={str(k):v.tolist() for k,v in curves.items()}
    record['electronic_cross_section_per_formula']=stopping_per_formula(curves[1],record['volume_A3']).tolist()
    destination=output/entry['name'];destination.mkdir()
    write_json(destination/'stopping.json',record)
    return record,error


def prior_controls(records):
    base=records['q2-fine'];fine=np.asarray(base['stopping_eV_A'])
    coarse=np.asarray(base['frequency_decimation_stopping_eV_A']['2'])
    error=abs(coarse-fine)+abs(np.asarray(records['eta-half']['stopping_eV_A'])-fine)
    # These three Native controls share the .0625eV frequency grid, not
    # the .03125eV grid. Match that premise before adding their changes.
    for name in ('bands72','PW700','k16'):
        error+=abs(np.asarray(records[name]['stopping_eV_A'])-coarse)
    return dict(schema='P257-replayed-five-control-envelope-v1',energy_lab_keV=base['energy_lab_keV'],
                partial_L1_control_envelope_eV_A=error.tolist())


def numeric_difference(expected,actual):
    if isinstance(expected,dict):
        require(expected.keys()==actual.keys(),'Changed result keys')
        return max((numeric_difference(v,actual[k]) for k,v in expected.items()),default=0.)
    if isinstance(expected,list):
        require(len(expected)==len(actual),'Changed result support')
        return max((numeric_difference(a,b) for a,b in zip(expected,actual)),default=0.)
    if isinstance(expected,(int,float)):
        return abs(expected-actual)
    require(expected==actual,'Changed result semantics')
    return 0.


def replay(index,output):
    output.mkdir(parents=True,exist_ok=False)
    records={};differences={}
    for entry in index['datasets']:
        records[entry['name']],differences[entry['name']]=replay_dataset(entry,index['constants'],output)
    prior=prior_controls(records)
    expected_prior=checked_json(index['prior_controls'])
    prior_error=numeric_difference(expected_prior['partial_L1_control_envelope_eV_A'],prior['partial_L1_control_envelope_eV_A'])
    require(prior_error<1e-12,'Matched-control reconstruction changed')
    prior_path=output/'prior-controls.json';write_json(prior_path,prior)
    analysis=analyze(prior_path,output/'q2-fine/stopping.json',output/'G100/stopping.json',output/'q4/stopping.json')
    expected_analysis=checked_json(index['control_analysis'])
    envelope_error=numeric_difference(expected_analysis['combined_observed_L1_envelope_eV_A'],analysis['combined_observed_L1_envelope_eV_A'])
    require(envelope_error<1e-12,'Refined settings envelope changed')
    write_json(output/'control-analysis.json',analysis)
    checked_json(index['bridge'])
    decision=build_tables(records['q4'],analysis,Path(index['bridge']['path']),output/'tables')
    reference=checked_json(index['q4_returns'])
    actual=json.loads((output/'tables/q4/return-curves.json').read_bytes())
    G_error=numeric_difference(reference['curves'],actual['curves'])
    require(G_error<1e-10,'Original q4 return consumer changed')
    receipt=dict(schema='P257-portable-replay-verification-v1',passed=True,
        stopping_max_abs_difference_eV_A=differences,prior_envelope_max_abs_difference_eV_A=prior_error,
        refined_envelope_max_abs_difference_eV_A=envelope_error,original_q4_all_curve_max_numeric_difference=G_error,
        numpy_version=np.__version__,scope='Actual saved-loss contraction, matched controls, crossing and three G bands; not a Native RPA rerun or a physical gain verdict')
    write_json(output/'verification.json',receipt)
    return receipt,decision


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index',type=Path,default=Path(__file__).with_name('index.json'))
    parser.add_argument('--output',type=Path,default=Path('replay-results'))
    args=parser.parse_args();source=args.index.resolve();output=args.output.resolve()
    index=json.loads(source.read_bytes());os.chdir(source.parent)
    receipt,_=replay(index,output)
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    main()
