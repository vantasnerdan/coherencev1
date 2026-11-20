"""Recompute physical readouts from full trajectories, forces, densities and BO energies."""
import numpy as np
from capsule_io import require
from green import COULOMB_EV_A, PeriodicNuclearGreen


def close(a, b, label, tolerance=1e-9):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(), label + ': nonfinite/shape')
    difference = float(np.max(np.abs(a - b))) if a.size else 0.
    require(difference <= tolerance, label + ': replay difference ' + str(difference))
    return difference


def density_potential(density, cell, point):
    reciprocal = 2 * np.pi * np.linalg.inv(cell).T
    axes = [np.fft.fftfreq(n) * n for n in density.shape]
    vectors = (np.stack(np.meshgrid(*axes, indexing='ij'), axis=-1) @ reciprocal).reshape(-1, 3)
    q2 = np.einsum('ij,ij->i', vectors, vectors)
    kernel = np.zeros(q2.shape)
    kernel[q2 > 0] = -4 * np.pi * COULOMB_EV_A / q2[q2 > 0]
    coefficients = np.fft.fftn(density).ravel() / density.size * kernel
    return float((coefficients * np.exp(1j * (vectors @ point))).sum().real)


def load_case(capsule, phase, role):
    directory = 'runs/' + phase + '-' + role + '/calculation/'
    ground_name = directory + 'ground/ground.json'
    trajectory_name = directory + 'trajectory/trajectory.json'
    bo_name = directory + 'bo/excitation.json'
    ground = capsule.json(ground_name)
    td = capsule.json(trajectory_name)
    bo = capsule.json(bo_name)
    capsule.bind(td['ground'], ground_name)
    capsule.bind(bo['trajectory'], trajectory_name)
    require(td['specification'] == bo['specification'], 'BO/trajectory input epochs differ')
    spec = td['specification']
    for key in ('h_A', 'kmesh', 'repeat', 'D_displacement_A', 'occupied_target', 'direction', 'impact_A', 'start_R_A', 'relaxed_ground_sha256'):
        require(ground['specification'][key] == spec[key], 'Reused ground changed ' + key)
    if phase == 'base':
        reuse = capsule.json('receipts/collision-production-ground-reuse.json')
        admitted = [state for state in reuse['states'] if state['role'] == role]
        require(len(admitted) == 1 and admitted[0]['gpw_sha256'] == ground['ground_gpw']['sha256']
                and admitted[0]['electrons'] == ground['electrons'], 'Ground reuse source differs')
    require(spec['occupied_target'] == (role == 'occ'), 'Wrong occupied/vacant source role')
    require(spec['D_displacement_A'] == ([0, 0, 0] if phase == 'base' else [.0768, 0, 0]), 'Wrong displacement source')
    require(spec['energy_keV'] == 2 and spec['h_A'] == .2 and spec['dt_as'] == 2 and spec['impact_A'] == 0
            and spec['kmesh'] == [2, 2, 2] and spec['solver_tolerance'] == 1e-9, 'Wrong production settings')
    rows = td['rows']
    require(len(rows) == 176 and len(td['checkpoints']) == len(bo['checkpoints']) == 7, 'Incomplete trajectory/support')
    require(all(b['time_as'] > a['time_as'] and b['R_A'] < a['R_A'] for a, b in zip(rows, rows[1:])), 'Nonmonotone trajectory')
    positions = np.array([r['projectile_A'] for r in rows])
    forces = np.array([r['projectile_force_eV_A'] for r in rows])
    work = np.concatenate(([0.], np.cumsum(-np.einsum('ij,ij->i', .5 * (forces[1:] + forces[:-1]), np.diff(positions, axis=0)))))
    diffs = {'force_work_eV': close(work, [r['external_work_eV'] for r in rows], 'Full trajectory force work')}
    particle = np.array([r['integrated_all_electrons'] for r in rows])
    require(np.max(np.abs(particle - particle[0])) < 1e-5, 'Particle conservation failed')
    checks = td['fixed_projector_contractions']
    require(checks['moved_projector_refused'] and max(max(r['addition_max_abs'], r['integration_max_abs']) for r in checks['actual_native_checks']) < 1e-10,
            'Unqualified fixed projector operator')
    host = bo['neutral_host_reference']
    require(host['charge'] == 0 and host['electrons'] == ground['electrons'], 'Wrong neutral host')
    close(host['cell_A'], ground['cell_A'], 'Host cell')
    index = ground['path']['projectile_index']
    require(ground['symbols'][:index] + ground['symbols'][index + 1:] == host['symbols'], 'Host species mismatch')
    close(np.delete(ground['positions_A'], index, axis=0), host['positions_A'], 'Host geometry')
    excitation = []
    checkpoint_work = []
    for point, checkpoint in zip(bo['checkpoints'], td['checkpoints']):
        require(point['row_index'] == checkpoint['row_index'] and point['R_A'] == checkpoint['R_A'], 'BO/checkpoint epoch differs')
        capsule.bound_file(point['actual_checkpoint'], directory + 'trajectory')
        actual = capsule.json(directory + 'trajectory/' + point['actual_checkpoint']['path'])
        require(actual == checkpoint, 'Checkpoint metadata mismatch')
        row = rows[point['row_index']]
        close(actual['positions_A'][index], row['projectile_A'], 'Checkpoint projectile')
        close(row['R_A'], np.linalg.norm(np.array(row['projectile_A']) - ground['path']['target_A']), 'Actual radius')
        close(point['td_internal_energy_eV'], row['td_internal_energy_eV'], 'TD energy epoch')
        value = row['td_internal_energy_eV'] - point['bo']['internal_energy_eV']
        close(value, point['excitation_internal_energy_eV'], 'Excitation from TD minus BO')
        excitation.append(value)
        checkpoint_work.append(float(work[point['row_index']]))
        close(work[point['row_index']], point['external_work_eV'], 'BO work binding')
        density_name = capsule.bound_file(checkpoint['density'], directory + 'trajectory')
        with np.load(capsule.root / density_name, allow_pickle=False) as density:
            data = density['density_e_A3']
            require(data.shape == tuple(2 * np.array(ground['realspace_shape'])), 'Wrong density grid')
            close(density['cell_A'], ground['cell_A'], 'Density cell')
            phi = density_potential(data, np.asarray(ground['cell_A']), row['projectile_A'])
            diffs['density_potential_eV'] = max(diffs.get('density_potential_eV', 0), close(phi, row['electron_potential_eV'], 'Density potential'))
            close(np.sum(data) * np.linalg.det(ground['cell_A']) / data.size, row['integrated_all_electrons'], 'Density electron count')
    closure = np.asarray(checkpoint_work) - (np.array([rows[p['row_index']]['td_internal_energy_eV'] for p in bo['checkpoints']]) - rows[0]['td_internal_energy_eV'])
    require(np.max(np.abs(closure)) < .5, 'Force-work energy closure failed')
    return dict(ground=ground, td=td, bo=bo, host=host, excitation=excitation, work=checkpoint_work,
                differences=diffs, max_force_work_energy_residual_eV=float(np.max(np.abs(closure))),
                bindings=dict(ground=ground_name, trajectory=trajectory_name, bo=bo_name))


def pair(capsule, occupied, vacant, width):
    a = occupied['td']['specification']; b = vacant['td']['specification']
    require({k:v for k,v in a.items() if k != 'occupied_target'} == {k:v for k,v in b.items() if k != 'occupied_target'}, 'Pair settings differ')
    target = np.array(occupied['ground']['path']['target_A'])
    host = occupied['host']; other = vacant['host']
    removed = [i for i, (s, p) in enumerate(zip(host['symbols'], host['positions_A'])) if s == 'H' and np.linalg.norm(np.array(p) - target) < 1e-10]
    require(len(removed) == 1, 'No unique neutral target')
    i = removed[0]
    require(host['symbols'][:i] + host['symbols'][i + 1:] == other['symbols'] and host['electrons'] - other['electrons'] == 1, 'Vacancy species/electron mismatch')
    close(np.delete(host['positions_A'], i, axis=0), other['positions_A'], 'Neutral vacancy geometry')
    close(host['cell_A'], other['cell_A'], 'Pair cells')
    green = PeriodicNuclearGreen(host['cell_A'])
    rows = []
    for j, (o, v) in enumerate(zip(occupied['bo']['checkpoints'], vacant['bo']['checkpoints'])):
        require(o['R_A'] == v['R_A'], 'Pair radii differ')
        to = occupied['td']['rows'][o['row_index']]; tv = vacant['td']['rows'][v['row_index']]
        close(to['projectile_A'], tv['projectile_A'], 'Paired projectile')
        raw = o['bo']['internal_energy_eV'] - v['bo']['internal_energy_eV'] - host['energy']['internal_energy_eV'] + other['energy']['internal_energy_eV']
        correction = green.point_minus_gaussian_eV(np.array(to['projectile_A']) - target, width)
        total = raw + correction; R = o['R_A']
        rows.append(dict(R_A=R, BO_quartet_soft_eV=raw, bare_periodic_point_minus_soft_eV=correction,
                         V_pair_periodic_point_eV=total, candidate_lowering_vs_isolated_C_over_R_eV=COULOMB_EV_A / R - total,
                         signed_electronic_candidate_eV=total - COULOMB_EV_A / R,
                         U_phi_diagnostic_eV=-(to['electron_potential_eV'] - tv['electron_potential_eV']),
                         occupied_excitation_eV=occupied['excitation'][j], vacant_excitation_eV=vacant['excitation'][j],
                         encounter_external_work_raw_eV=occupied['work'][j] - vacant['work'][j]))
    first = rows[0]
    initial = first['occupied_excitation_eV'] - first['vacant_excitation_eV']
    for row in rows:
        length = first['R_A'] - row['R_A']
        row['path_length_A'] = length
        row['encounter_excitation_change_eV'] = row['occupied_excitation_eV'] - row['vacant_excitation_eV'] - initial
        row['vacancy_loss_change_eV'] = row['vacant_excitation_eV'] - first['vacant_excitation_eV']
        row['vacancy_mean_loss_eV_A'] = row['vacancy_loss_change_eV'] / length if length else None
        row['encounter_external_work_mapped_increment_eV'] = row['encounter_external_work_raw_eV'] + row['bare_periodic_point_minus_soft_eV'] - first['bare_periodic_point_minus_soft_eV']
        row['anchored_encounter_work_profile_eV'] = first['V_pair_periodic_point_eV'] + row['encounter_external_work_mapped_increment_eV']
        row['independent_force_work_identity_residual_eV'] = row['anchored_encounter_work_profile_eV'] - row['V_pair_periodic_point_eV'] - row['encounter_excitation_change_eV']
    require(max(abs(r['independent_force_work_identity_residual_eV']) for r in rows) < .5, 'Independent paired work closure failed')
    return rows


def compare_rows(actual, expected, label):
    require(len(actual) == len(expected), label + ': support differs')
    differences = {}
    for row, old in zip(actual, expected):
        require(set(row) == set(old), label + ': readout columns differ')
        for key, value in row.items():
            if value is None:
                require(old[key] is None, label + ': invalid zero-length drag')
            else:
                differences[key] = max(differences.get(key, 0), close(value, old[key], label + '/' + key))
    return differences
