"""Authenticate actual saved primitives and their chemistry before thermodynamic subtraction."""
import hashlib
import json
import math
from pathlib import Path
from vacancy_geometry_audit import audit


def file_binding(path):
    data = Path(path).read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_saved(roots):
    records, geometry_audits, bindings = {}, {}, []
    for root in map(Path, roots):
        execution = json.loads((root / 'execution.json').read_text())
        require(execution['affinity'] == [7], 'Incorrect vacancy CPU provenance')
        for source in execution['sources']:
            actual = file_binding(source['path'])
            require(actual['sha256'] == source['sha256'], 'Execution input/source bytes changed')
        bindings.append(file_binding(root / 'execution.json'))
        for path in sorted(root.glob('*/primitive.json')):
            record = json.loads(path.read_text())
            case = record['case']
            require(case not in records, 'Duplicate primitive case')
            require(record['schema'] == 'coherence-vacancy-primitive-v1', 'Unsupported primitive schema')
            require(record['converged_relaxation'], 'Unrelaxed primitive cannot enter binding ledger')
            require(math.isfinite(record['energy_eV']), 'Nonfinite primitive energy')
            cfg = record['config']
            n_pd = cfg['repeat'] ** 3
            vacancy = case.startswith('vacancy')
            n_d = int(case.split('_D')[1]) if vacancy and '_D' in case else int(case == 'bulk_D')
            require(record['symbols'].count('Pd') == n_pd - int(vacancy),
                    'Metal vacancy requires exactly one missing Pd atom')
            require(record['symbols'].count('H') == n_d, 'Incorrect deuterium stoichiometry')
            require(record['geometry']['missing_Pd'] == vacancy, 'Vacancy definition mismatch')
            require(record['neutral_electrons'] == (n_pd-int(vacancy))*10+n_d,
                    'Neutral Pd10/H1 electron count mismatch')
            require(record['max_constrained_force_eV_A'] <= cfg['fmax_eV_A'],
                    'Force threshold not satisfied')
            require(cfg['setups'] == {'Pd': '10', 'H': 'paw'}, 'Wrong PAW selection')
            require(cfg['setup_sha256'] == {
                'Pd.10.PBE.gz': '3c6cfd1de7e6d6815d58c9521b2bc7cec27ede8e83554ed1f0151f979761bbb4',
                'H.PBE.gz': '3136ad73c02d471f99dd5f82e451b2061952970556a260d38b53b5020dc94ae7'},
                'Wrong PAW bytes')
            records[case] = record
            geometry_audits[case] = audit(path.parent / 'final.extxyz', record)
            bindings.extend(file_binding(path.parent / name) for name in
                            ('primitive.json', 'initial.extxyz', 'final.extxyz', 'final.gpw', 'relax.traj'))
    require(all(case in records for case in ('bulk', 'bulk_D', 'vacancy', 'vacancy_D1')),
            'Incomplete primitive ledger: bulk, bulk_D, vacancy, vacancy_D1 required')
    common = records['bulk']['config']
    matched_keys = ('lattice_A', 'repeat', 'ecut_eV', 'kmesh', 'occupation_width_eV',
                    'spinpol', 'fmax_eV_A', 'mobile_radius_A', 'setups', 'setup_sha256')
    require(all(all(record['config'][key] == common[key] for key in matched_keys)
                for record in records.values()), 'Unmatched primitive Hamiltonian/cell/relaxation protocol')
    return records, geometry_audits, bindings
