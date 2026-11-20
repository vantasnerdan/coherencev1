"""Exact six-point subset of the already admitted P255 B curve; zero Native calls."""
import argparse
import hashlib
import json
from pathlib import Path

SOURCE_SHA='854c570dc53ac5f0ff792f38393dea18c66dab7aa474d4f03b29d5e4e98a32a5'
RADII=[.1,.2,.3,.5,1.,1.6]


def reuse(source,output):
    raw=source.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:
        raise ValueError('Wrong original qualified B source epoch')
    data=json.loads(raw)
    if data['schema']!='coherencev1-target-conditioned-pair-v1' or data['observable']!='target_conditioned_BO_total_interaction':
        raise ValueError('Wrong B observable')
    selected=[data['R_A'].index(radius) for radius in RADII]
    result=dict(data)
    result['R_A']=RADII
    result['V_pair_total_eV']=[data['V_pair_total_eV'][index] for index in selected]
    result['g7_reuse']=dict(original_bytes=len(raw),original_sha256=SOURCE_SHA,selected_original_indices=selected,
        excluded_original_radii_A=[radius for radius in data['R_A'] if radius not in RADII],
        native_runs=0,scope='Requested exact subset, no interpolation or requalification of the unchanged admitted B curve.')
    with output.open('x') as handle:
        handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output_bytes=output.stat().st_size,sha256=hashlib.sha256(output.read_bytes()).hexdigest(),R_A=RADII,native_runs=0),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();reuse(args.source,args.output)
