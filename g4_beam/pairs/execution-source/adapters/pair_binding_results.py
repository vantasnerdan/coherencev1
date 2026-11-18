"""Consume completed actual cells with matched chemistry and all-mobile forces."""
import argparse
import json
from pathlib import Path
import numpy as np
from substrate_framework.cluster_binding import occupancy_bindings_eV
from force_bank_results import verify_force_record
from surface_io import manifest, save


def load_completed(root):
    root = Path(root)
    execution = json.loads((root/'execution.json').read_text())
    for bound in execution['sources']:
        if manifest(bound['path']) != bound:
            raise ValueError('Original execution source/input bytes changed')
    config = json.loads(Path(execution['sources'][0]['path']).read_text())
    rows = []
    for path in sorted(root.glob('*/primitive.json')):
        record = json.loads(path.read_text())
        case = record['case']
        declared = next(item for item in config['cases'] if item['name'] == case['name'])
        if case != declared or record['settings'] != config['settings']:
            raise ValueError('Completed case differs from its bound execution declaration')
        atoms = verify_force_record(record['bank'])
        counts = {symbol: list(atoms.symbols).count(symbol) for symbol in ('Pd','H')}
        if (len(atoms) != counts['Pd']+counts['H'] or counts['Pd'] != case['n_Pd'] or
                counts['H'] != case['n_D'] or case['host_Pd']-case['n_missing'] != counts['Pd'] or
                atoms.constraints or not all(atoms.pbc)):
            raise ValueError('Actual bank has unmatched chemistry, cell boundary or mobile subspace')
        if record['neutral_electrons'] != 10*counts['Pd']+counts['H']:
            raise ValueError('Wrong neutral Pd10/H1 electron accounting')
        if not np.allclose(atoms.get_masses()[np.asarray(atoms.numbers)==1],2.01410177812,rtol=0,atol=1e-12):
            raise ValueError('Wrong deuterium isotope')
        fmax = float(np.linalg.norm(np.asarray(record['bank']['forces']),axis=1).max())
        if not record['converged_relaxation'] or fmax > config['settings']['fmax_eV_A']:
            raise ValueError('Unrelaxed all-mobile state cannot enter binding ledger')
        rows.append(dict(case=case,energy_eV=record['bank']['energy_eV'],
                         fmax_eV_A=fmax,cell_A=atoms.cell.tolist(),source_primitive=manifest(path),
                         bank=record['bank']))
    return config, execution, rows


def consume(roots):
    groups = [load_completed(root) for root in roots]
    common = groups[0][0]
    producer = [item['sha256'] for item in groups[0][1]['sources'][1:]]
    rows = []
    for config, execution, group in groups:
        if (config['settings'] != common['settings'] or config['setup_sha256'] != common['setup_sha256'] or
                [item['sha256'] for item in execution['sources'][1:]] != producer):
            raise ValueError('Unmatched Hamiltonian, force protocol or actual producer')
        rows.extend(group)
    if len({row['case']['name'] for row in rows}) != len(rows):
        raise ValueError('Duplicate completed case')
    if not rows or any(not np.array_equal(row['cell_A'],rows[0]['cell_A']) for row in rows):
        raise ValueError('Binding subtraction requires the same cell')
    bulk = next(row for row in rows if row['case']['n_missing']==0 and row['case']['n_D']==0)
    bulk_D = next(row for row in rows if row['case']['n_missing']==0 and row['case']['n_D']==1)
    clusters = {}
    for row in rows:
        if row['case']['n_missing'] != 1:
            continue
        n = row['case']['n_D']
        previous = clusters.get(n)
        if previous is None or (row['energy_eV'],row['case']['name']) < (previous['energy_eV'],previous['case']['name']):
            clusters[n] = row
    values = {n:row['energy_eV'] for n,row in clusters.items()}
    result = occupancy_bindings_eV(values,bulk['energy_eV'],bulk_D['energy_eV'])
    declared = [case for config,_,_ in groups for case in config['cases']]
    requested = {case['n_D'] for case in declared if case['n_missing']==1}
    completed_names = {row['case']['name'] for row in rows}
    result.update(schema='coherence-pairs-matched-binding-v1',settings=common['settings'],
                  source_cells=rows,computed_occupancies=sorted(values),
                  selected_cluster_cases={n:row['case']['name'] for n,row in clusters.items()},
                  selection_policy='Lowest extrapolated energy among enumerated local basins at each occupancy; not a global minimum',
                  requested_occupancies=sorted(requested),
                  requested_occupancies_still_missing=sorted(requested-set(values)),
                  unsettled_declared_cases=sorted(case['name'] for case in declared if case['name'] not in completed_names),
                  ionic_free_energies_computed=False,capacity_proven=False,
                  current_transport_Eb_eV=0.23,replace_current_transport=False,
                  replacement_reason='D1/D2 size check and explicit reservoir mapping still required')
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True,action='append')
    parser.add_argument('--reference',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    result=consume([*args.root,args.reference])
    save(args.output,result)
    print(json.dumps({key:result[key] for key in ('blocks','computed_occupancies',
          'requested_occupancies_still_missing','replace_current_transport')}))
