"""Four actual AE energies; unrelaxed embedded-cell sensitivity, never a cell-convergence gate."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

SUPPORT = Path(__file__).resolve().parent
if not (SUPPORT/'green.py').is_file():
    SUPPORT = SUPPORT.parent/'delivery'
sys.path.insert(0, str(SUPPORT))
from green import BOHR_A, COULOMB_EV_A, PeriodicNuclearGreen
from capsule_io import require

ROLES = ('occupied-host', 'vacant-host', 'occupied-probe', 'vacant-probe')


def binding(path):
    raw = path.read_bytes()
    return dict(path=path.name, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def calculate(directory):
    from readout_defect import energy, close
    directory = Path(directory)
    spec = json.loads((directory/'input.json').read_bytes())
    require(spec['sensitivity_scope'] == 'unrelaxed_embedded_local_geometry', 'Wrong cell sensitivity scope')
    target_index = spec['target_index']; symbols = spec['symbols']
    positions = np.asarray(spec['initial_positions_A']); cell = np.asarray(spec['cell_A'])
    direction = np.asarray(spec['probe_direction']); target = positions[target_index]
    require(symbols.count('Pd') == 26 and symbols.count('H') == 1 and symbols[target_index] == 'H', 'Wrong embedded inventory')
    close(np.linalg.norm(direction), 1., 'Frozen unit ray')
    states = {}; wrappers = {}; source = json.loads((directory/'source-manifest.json').read_bytes())
    for item in source['files']:
        require(binding(directory/'source'/item['path']) == dict(item, path=Path(item['path']).name), 'Changed original role executable')
    for role in ROLES:
        folder = directory/'roles'/role
        wrapper = json.loads((folder/'cell-role.json').read_bytes())
        require(wrapper['completed'] is True and wrapper['role'] == role and wrapper['force_converged'] is False,
                'Incomplete or mislabelled sensitivity role')
        require(wrapper['input'] == binding(directory/'input.json') and
                wrapper['frozen_source'] == binding(directory/'source-manifest.json') and
                wrapper['state'] == binding(folder/'state.json'), 'Changed sensitivity role epoch')
        state = json.loads((folder/'state.json').read_bytes()); states[role] = state
        require(state['specification'] == spec, 'Different sensitivity physical input')
        expected_symbols = list(symbols); expected_positions = positions.copy()
        if role.startswith('vacant'):
            del expected_symbols[target_index]; expected_positions = np.delete(expected_positions, target_index, axis=0)
        charge = int(role.endswith('probe'))
        if charge:
            expected_symbols.append('H'); expected_positions = np.vstack((expected_positions, target-.1*direction))
        require(state['symbols'] == expected_symbols and state['charge'] == charge, 'Wrong sensitivity role inventory/charge')
        close(state['positions_A'], expected_positions, 'Frozen embedded geometry')
        close(state['cell_A'], cell, 'Same sensitivity cell')
        energy(state)
        wrappers[role] = dict(state=binding(folder/'state.json'), role=binding(folder/'cell-role.json'))
    for role in ('occupied', 'vacant'):
        close(states[role+'-probe']['electrons'], states[role+'-host']['electrons'], 'Probe must not add an electron')
    close(states['occupied-host']['electrons']-states['vacant-host']['electrons'], 1., 'Neutral target removal')
    raw = math.fsum([energy(states['occupied-probe']), -energy(states['vacant-probe']),
                     -energy(states['occupied-host']), energy(states['vacant-host'])])
    correction = PeriodicNuclearGreen(cell).point_minus_gaussian_eV(-.1*direction, math.sqrt(5.)/BOHR_A)
    total = raw+correction
    return dict(schema='pd-defect-cell-sensitivity-v1', R_A=.1,
                row=dict(R_A=.1, BO_quartet_soft_eV=raw, point_minus_gaussian_eV=correction,
                         V_pair_total_eV=total, lowering_eV=COULOMB_EV_A/.1-total),
                input=binding(directory/'input.json'), frozen_source=binding(directory/'source-manifest.json'),
                states=wrappers, base_force_frame=spec['seed_source']['force_frame'],
                cell_A=cell.tolist(), direction=direction.tolist(), source_geometry=spec['seed_source'],
                electronic_states_converged=True, forces_computed=False, force_converged=False, cell_converged=False,
                scope='8-site relaxed local Pd/D patch embedded in a 27-site cell with the far field unrelaxed. One R0.1 electronic sensitivity only; NOT a relaxed or converged cell-size check.')


def replay(directory):
    from readout_defect import close
    actual = calculate(directory)
    saved = json.loads((Path(directory)/'cell-sensitivity.json').read_bytes())
    for key, value in actual['row'].items():
        close(saved['row'][key], value, 'Saved unrelaxed sensitivity '+key)
    a, s = dict(actual), dict(saved); del a['row']; del s['row']
    require(a == s, 'Changed sensitivity provenance/status')
    return actual


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    args = parser.parse_args()
    value = calculate(args.directory)
    payload = json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n'
    with (args.directory/'cell-sensitivity.json').open('x') as stream:
        stream.write(payload)
    print(payload, end='')
