"""Independent saved-state replay of representative h and neutral image gates."""
import itertools
import math
import re
import numpy as np
from capsule_io import require
from green import COULOMB_EV_A, PeriodicNuclearGreen
from replay import close


def controls(capsule, cases, baseline, width):
    result = {}
    base = cases['base-occ']
    for control in ('h016', 'image-2x'):
        states = {}
        receipt_name = 'receipts/collision-' + ('h016-static' if control == 'h016' else 'image') + '-comparison-001.json'
        receipt = capsule.json(receipt_name)
        for role, short in (('occupied', 'occ'), ('vacant', 'vac')):
            case = cases['base-' + short]
            for state in ('probe', 'host'):
                key = role + '-' + state
                name = 'controls/' + control + '/calculation/' + key + '/state.json'
                record = capsule.json(name)
                capsule.bind(receipt['state_bindings'][key], name)
                require(record['role'] == role and record['state'] == state, 'Wrong static control role')
                require(record['charge'] == (1 if state == 'probe' else 0), 'Wrong control charge')
                spec = record['specification']; original = case['td']['specification']
                for field in ('D_displacement_A', 'direction', 'impact_A', 'relaxed_ground_sha256', 'solver_tolerance'):
                    require(spec[field] == original[field], 'Control changed ' + field)
                require(spec['h_A'] == (.16 if control == 'h016' else .2), 'Wrong control h')
                require(spec['repeat'] == ([1, 1, 1] if control == 'h016' else [2, 2, 2]), 'Wrong image repetition')
                require(spec['kmesh'] == ([2, 2, 2] if control == 'h016' else [1, 1, 1]), 'Wrong control kmesh')
                if control == 'h016':
                    close(record['cell_A'], case['ground']['cell_A'], 'h cell')
                    require(record['actual_kpoints'] == 8 and record['realspace_shape'] == [24, 24, 24], 'Wrong fine-grid support')
                    symbols = case['host']['symbols'] if state == 'host' else case['ground']['symbols']
                    positions = case['host']['positions_A'] if state == 'host' else case['td']['checkpoints'][-1]['positions_A']
                    require(record['symbols'] == symbols and record['electrons'] == case['ground']['electrons'], 'h species/electrons differ')
                    close(record['positions_A'], positions, 'h matched geometry')
                else:
                    close(record['cell_A'], 2 * np.array(case['ground']['cell_A']), 'Approved doubled cell')
                    require(record['actual_kpoints'] == 1 and record['realspace_shape'] == [40, 40, 40], 'Wrong image grid/k support')
                    require(record['electrons'] == (352 if role == 'occupied' else 351), 'Wrong neutral image electrons')
                    text = (capsule.root / ('controls/' + control + '/calculation/' + key + '/scf.txt')).read_text()
                    bands = re.search(r'^Number of bands in calculation:\s+(\d+)', text, re.M)
                    require(bands is not None and int(bands[1]) == 240, 'Missing folded image band support')
                states[key] = record
            probe = states[role + '-probe']; host = states[role + '-host']; index = probe['path']['projectile_index']
            require(probe['symbols'][:index] + probe['symbols'][index + 1:] == host['symbols'], 'Control host removes wrong atom')
            close(np.delete(probe['positions_A'], index, axis=0), host['positions_A'], 'Control host geometry')
        energy = lambda role, state: states[role + '-' + state]['energy']['internal_energy_eV']
        if control == 'h016':
            raw = energy('occupied', 'probe') - energy('vacant', 'probe') - energy('occupied', 'host') + energy('vacant', 'host')
            delta = raw - baseline[-1]['BO_quartet_soft_eV']
            lowering = baseline[-1]['candidate_lowering_vs_isolated_C_over_R_eV'] - delta
            computed = dict(raw_quartet_h020_eV=baseline[-1]['BO_quartet_soft_eV'], raw_quartet_h016_eV=raw,
                            delta_V_pair_h016_minus_h020_eV=delta, candidate_positive_lowering_h020_eV=baseline[-1]['candidate_lowering_vs_isolated_C_over_R_eV'],
                            candidate_positive_lowering_h016_eV=lowering, delta_positive_lowering_eV=-delta,
                            delta_signed_electronic_correction_eV=delta, tolerance_eV=max(5., .1 * abs(lowering)))
        else:
            host = states['occupied-host']; vacant = states['vacant-host']; probe = states['occupied-probe']
            original = base['host']; cell = np.array(original['cell_A'])
            repeated = np.concatenate([np.array(original['positions_A']) + np.array(shift) @ cell for shift in itertools.product(range(2), repeat=3)])
            require(host['symbols'] == original['symbols'] * 8, 'Image species are not repeated original material')
            close(host['positions_A'], repeated, 'Repeated image host')
            target = np.array(probe['path']['target_A'])
            indices = [i for i, (s, p) in enumerate(zip(host['symbols'], host['positions_A'])) if s == 'H' and np.linalg.norm(np.array(p) - target) < 1e-10]
            require(len(indices) == 1, 'Image target not unique')
            i = indices[0]
            require(host['symbols'][:i] + host['symbols'][i+1:] == vacant['symbols'], 'Wrong image vacancy species')
            close(np.delete(host['positions_A'], i, axis=0), vacant['positions_A'], 'Image target-only removal')
            point = np.array(probe['positions_A'][probe['path']['projectile_index']])
            close(point, base['td']['rows'][-1]['projectile_A'], 'Image projectile phase')
            close(point, states['vacant-probe']['positions_A'][-1], 'Image paired projectile')
            vector = point - target
            close(np.linalg.norm(vector), .1, 'Image radius')
            raw = math.fsum(sign * states[key]['energy']['internal_energy_eV'] for key, sign in (('occupied-probe', 1), ('vacant-probe', -1), ('occupied-host', -1), ('vacant-host', 1)))
            small = PeriodicNuclearGreen(base['ground']['cell_A']); large = PeriodicNuclearGreen(probe['cell_A'])
            correction = large.point_minus_gaussian_eV(vector, width)
            total = raw + correction; lowering = COULOMB_EV_A / .1 - total
            computed = dict(raw_neutral_quartet_2x_eV=raw, point_minus_soft_2x_eV=correction,
                            V_pair_periodic_point_2x_eV=total, candidate_lowering_base_eV=baseline[-1]['candidate_lowering_vs_isolated_C_over_R_eV'],
                            candidate_lowering_2x_eV=lowering, delta_2x_minus_base_lowering_eV=lowering-baseline[-1]['candidate_lowering_vs_isolated_C_over_R_eV'],
                            tolerance_eV=max(5., .1 * abs(lowering)), raw_quartet_delta_eV=raw-baseline[-1]['BO_quartet_soft_eV'],
                            point_minus_soft_delta_eV=correction-baseline[-1]['bare_periodic_point_minus_soft_eV'],
                            periodic_bare_point_kernel_delta_eV=large.point_eV(vector)-small.point_eV(vector))
        differences = {key: close(value, receipt[key], control + '/' + key) for key, value in computed.items()}
        delta_key = 'delta_V_pair_h016_minus_h020_eV' if control == 'h016' else 'delta_2x_minus_base_lowering_eV'
        require(abs(computed[delta_key]) <= computed['tolerance_eV'] and receipt['passed'], 'Static gate failed')
        result[control] = dict(values=computed, differences=differences, passed=True)
    return result
