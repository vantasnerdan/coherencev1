"""Collect completed Native batches; preserve raw source bindings and model scope."""
import argparse
import json
from pathlib import Path

from screen_io import manifest, require, write_json
from screen_readout import read_case
from screen_report import report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--batch',type=Path,action='append',default=[])
    parser.add_argument('--case',type=Path,action='append',default=[])
    parser.add_argument('--allow-incomplete',action='store_true')
    parser.add_argument('--case-specs',type=Path,required=True)
    parser.add_argument('--saved-dos-supplement',type=Path,required=True)
    parser.add_argument('--pd-dos-control',type=Path,required=True)
    parser.add_argument('--pd-length-control',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--neutral-cycle',type=Path)
    args=parser.parse_args()
    supplement=json.loads(args.saved_dos_supplement.read_bytes())
    rows=[]; ledger=[];folders=list(args.case)
    for batch in args.batch:
        require((batch/'batch-completed.json').is_file(), 'Native batch is not complete')
        folders += [d for d in batch.iterdir() if d.is_dir() and (d/'completed.json').is_file()]
    require(folders,'No completed scientific checkpoints selected')
    for folder in sorted(folders):
        require((folder/'completed.json').is_file(),'Selected scientific case is not complete')
        config=json.loads((folder/'input.json').read_bytes())
        material=config['material_id']
        rows.append(read_case(folder,material,supplement if material=='Pd_hydro_+0pct' else None))
        ledger.append(dict(material_id=material, config=config,
            convergence=dict(DOS='Requested Pd-only k control and declared transferred5% target; occupied projections independently integrated',
                             RPA=(rows[-1]['screening']['eligible_for_Tier1'] if 'screening' in rows[-1] else None),
                             physical_scope='Numerical protocol eligibility, not all-electron physical confidence'),
            consumers=dict(DOS_and_slope=['g4_beam'],lambda_and_ratios=['g1_rank','g3_lattice']),
            source_case=folder.name, source_phase=folder.parent.name,
            checkpoint_from_incomplete_batch=not (folder.parent/'batch-completed.json').is_file(),
            input=manifest(folder/'input.json'),completed=manifest(folder/'completed.json')))
    expected={v['case'] for v in json.loads(args.case_specs.read_bytes())}
    actual={r['material_id'] for r in rows}
    require(actual<=expected and (args.allow_incomplete or actual==expected),
            'Completed materials differ from commissioned Tier1 cases')
    scalars=report(rows,require_complete=not args.allow_incomplete)
    scalars['not_computed_materials']=sorted(expected-actual)
    if args.neutral_cycle is not None:
        neutral=json.loads(args.neutral_cycle.read_bytes())
        require(neutral['schema']=='coherence-neutral-adiabatic-v1','Wrong neutral-cycle observable')
        scalars['neutral_adiabatic']=neutral
    scalars['numerical_controls']=dict(
        Pd_DOS=json.loads(args.pd_dos_control.read_bytes()),
        Pd_lambda=json.loads(args.pd_length_control.read_bytes()),
        control_scope='Requested Pd-only k control; per-case q-window; Pd band/cutoff/eta controls transferred as an explicit model allowance, not independent all-material convergence.')
    args.output.mkdir()
    write_json(args.output/'qualified-cases.json',rows)
    write_json(args.output/'scalars.json',scalars)
    write_json(args.output/'input-ledger.json',dict(cases=ledger,
        requested_cases=manifest(args.case_specs),saved_DOS_supplement=manifest(args.saved_dos_supplement),
        Pd_DOS_control=manifest(args.pd_dos_control),Pd_length_control=manifest(args.pd_length_control),
        neutral_cycle=(manifest(args.neutral_cycle) if args.neutral_cycle is not None else None),
        collector=manifest(Path(__file__))))
    print(json.dumps(dict(cases=len(rows),response_cases=sum('screening' in r for r in rows),
                         output=manifest(args.output/'scalars.json'),
                         headline_slope=scalars['headline_static_Yukawa_proxy_slope']),indent=2))


if __name__=='__main__':
    main()
