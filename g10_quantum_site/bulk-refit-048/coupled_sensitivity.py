"""One changed PdD site in the saved joint old primitive path ensembles.

Run from any directory with candidate variant names from well.py as arguments.
Only this site's on-site term changes; its co-beaded neighbor field is retained.
"""

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
G10 = HERE.parent
sys.path.insert(0, str(G10))
from coupled_momentum import HBAR2_OVER_2AMU_EV_A2, KB_EV_PER_K, displaced_path
from model import MASS, full, potential

BLOCK_SIZE = 89  # 9 equal contiguous time blocks using all 801 frames
EDGES = np.arange(41) / 50  # 0.02 Å shells; the overflow is reported separately
EXTERIORS = (0.24, 0.48, 0.60)
DIRECTION = np.ones(3) / np.sqrt(3)
DISPLACEMENTS = np.arange(13) / 20
CUTS = 4
WITNESS_A = 0.25


def ratio_stats(values, weights):
    """Columns are per-frame observables; return ratio and block influence rows."""
    den = weights.reshape(-1, BLOCK_SIZE).sum(axis=1)
    num = (weights[:, None] * values).reshape(len(den), BLOCK_SIZE, -1).sum(axis=1)
    value = num.sum(axis=0) / den.sum()
    return value, value + (num - den[:, None] * value) / den.mean()


def se(blocks):
    return blocks.std(axis=0, ddof=1) / np.sqrt(len(blocks))


def ess(weights):
    return float(weights.sum() ** 2 / (weights @ weights))


def path_observables(x):
    radii = np.linalg.norm(x, axis=-1)
    num_bins = len(EDGES)  # len(EDGES)-1 finite shells, then overflow
    bins = np.searchsorted(EDGES[1:], radii, side='right')
    by_frame = bins + num_bins * np.arange(len(x))[:, None]
    histogram = np.bincount(by_frame.ravel(), minlength=num_bins * len(x))
    histogram = histogram.reshape(len(x), num_bins) / x.shape[1]
    columns = {'MSD_A2': (radii**2).mean(axis=1)}
    for cut in EXTERIORS:
        outside = radii > cut
        columns[f'bead_fraction_gt_{cut:.2f}_A'] = outside.mean(axis=1)
        columns[f'frame_any_gt_{cut:.2f}_A'] = outside.any(axis=1)
    columns['histogram_overflow_at_or_above_0.80_A'] = histogram[:, -1]
    return columns, histogram[:, :-1]


def energy_callback(local):
    # The second slot carries the field in eV/Å, not a dynamical particle or
    # an extra spring. The borrowed operator shifts only slot 0. All other
    # terms in the full bead PES are constant in this shift and cancel.
    def energy(packed):
        position, field = packed[..., 0, :], packed[..., 1, :]
        return local(position) - np.einsum('fbi,fbi->fb', position, field)
    return energy


def characteristic(packed, temperature, local, logw, coefficients=None):
    energy = ({'energy_eV': energy_callback(local)} if coefficients is None
              else {'site_coefficients': coefficients})
    return displaced_path(packed, MASS['D'], temperature, 0, DISPLACEMENTS,
                          direction=DIRECTION, block_size=BLOCK_SIZE, cuts=CUTS,
                          log_frame_weights=logw, **energy)


def direct_witness(x, field, old_bead, local, logw, beta):
    """Independent old-base to shifted-target numerator at a nonzero separation."""
    p = x.shape[1]
    h = (np.arange(p) + 0.5) / p - 0.5
    free_action = WITNESS_A**2 * MASS['D'] / (4 * beta * HBAR2_OVER_2AMU_EV_A2)
    values, exterior = [], []
    for cut in range(CUTS):
        d = np.roll(h, cut * p // CUTS)[:, None] * (WITNESS_A * DIRECTION)
        for shift in (-d, d):
            shifted = x + shift
            delta = local(shifted) - old_bead - np.einsum('fbi,bi->fb', field, shift)
            wl = -beta * delta.mean(axis=1) - free_action
            values.append(wl)
            exterior.append((np.linalg.norm(shifted, axis=-1) > 0.48).any(axis=1))
    vals = np.asarray(values)
    total = logsumexp(vals)
    outside = np.asarray(exterior)
    return (float(np.exp(total - np.log(2 * CUTS) - logsumexp(logw))),
            float(np.exp(logsumexp(vals[outside]) - total)))


def checked_energy(energy, x):
    e = np.asarray(energy(x), dtype=float)
    if e.shape != x.shape[:-1] or not np.isfinite(e).all():
        raise ValueError('Expected one finite local energy per position from well.energy_force')
    return e


def analyze(temperature, local_models, old, well_old, saved_row):
    case = f'PdD_{temperature}K'
    path = G10 / 'output' / 'coupled' / f'{case}_paths.npz'
    with np.load(path) as data:
        x = data['beads_A']
        field = data['neighbor_field_eV_A']
    expected = (801, {300: 32, 77: 128}[temperature], 3)
    if x.shape != expected or field.shape != expected or not np.isfinite(x).all() or not np.isfinite(field).all():
        raise ValueError(f'{case}: unexpected retained path or co-beaded field')
    packed = np.stack((x, field), axis=2)
    site = old['site']
    beta = 1 / (KB_EV_PER_K * temperature)
    old_local = dict(old, K=old['K'][site:site + 1], f=old['f'][site:site + 1], J=np.zeros((1, 1)))

    def old_energy(q):
        return potential(q[..., None, :], old_local)[0]

    old_bead = old_energy(x)
    old_model_error = float(np.max(np.abs(checked_energy(well_old, x) - old_bead)))
    if old_model_error > 1e-12:
        raise ValueError(f'{case}: well.old does not match the archived source PES')
    columns, histogram = path_observables(x)
    names = tuple(columns)
    scalars = np.column_stack(tuple(columns.values()))
    old_n = (len(x) // 12) * 12  # archived position block_stats retained the last 792
    archive_msd_error = abs(float(scalars[-old_n:, 0].mean()) - saved_row['mean_square_radius_A2'])
    archive_exterior_error = abs(float(columns['bead_fraction_gt_0.24_A'].mean()) - saved_row['position_outside_support_fraction'])
    if max(archive_msd_error, archive_exterior_error) > 1e-12:
        raise ValueError(f'{case}: archived old raw bead controls disagree')

    observables = np.column_stack((scalars, histogram))
    old_v, old_blocks = ratio_stats(observables, np.ones(len(x)))
    old_poly = characteristic(packed, temperature, old_energy, None,
                              (old['K'][site], old['f'][site], old['a'], old['b'], field))
    case_info = {
        'temperature_K': temperature, 'beads': x.shape[1], 'frames': len(x),
        'case': case,
        'contiguous_blocks': len(x) // BLOCK_SIZE, 'frames_per_block': BLOCK_SIZE,
        'path_file': str(path.relative_to(G10)), 'path_sha256': sha256(path.read_bytes()).hexdigest(),
        'max_saved_bead_radius_A': float(np.linalg.norm(x, axis=-1).max()),
        'neighbor_field_rms_eV_A': float(np.sqrt(np.mean(np.sum(field**2, axis=-1)))),
        'archived_old_MSD_control_abs_error_A2': archive_msd_error,
        'archived_old_exterior_0.24_A_control_abs_error': archive_exterior_error,
        'well_old_vs_archive_max_abs_bead_energy_error_eV': old_model_error,
        'variants': {},
    }
    position_rows, characteristic_rows = [], []
    estimates = {}
    for name, local in [('old', old_energy), *local_models.items()]:
        delta_bead = checked_energy(local, x) - old_bead
        logw = -beta * delta_bead.mean(axis=1)
        weights = np.exp(logw - logw.max())
        w_block = weights.reshape(-1, BLOCK_SIZE).sum(axis=1)
        frame_ess, block_ess = ess(weights), ess(w_block)
        # This is the old G10 participation rule, not proof against unseen modes.
        sampled_overlap = frame_ess >= 50 and block_ess >= 4
        val, blocks = ratio_stats(observables, weights)
        estimates[name] = val, blocks
        errors, paired = se(blocks), se(blocks - old_blocks)
        metrics = {}
        for j, key in enumerate(names):
            positive_blocks = int(np.count_nonzero((weights * scalars[:, j]).reshape(-1, BLOCK_SIZE).sum(axis=1)))
            metrics[key] = {'value': val[j], 'block_se': errors[j], 'change_vs_old': val[j] - old_v[j],
                            'paired_change_block_se': paired[j], 'nonzero_numerator_blocks': positive_blocks}
        info = {
            'overlap': {'frame_weight_ess': frame_ess, 'block_weight_ess': block_ess,
                        'largest_frame_weight_fraction': float(weights.max() / weights.sum()),
                        'largest_block_weight_fraction': float(w_block.max() / w_block.sum()),
                        'log_weight_span': float(np.ptp(logw)), 'passes_participation_screen': sampled_overlap},
            'position': metrics,
            'position_normalization': float(val[len(names):].sum() + metrics[names[-1]]['value']),
        }
        for k, (lo, hi) in enumerate(zip(EDGES[:-1], EDGES[1:])):
            j = len(names) + k
            width = hi - lo
            volume = (4 * np.pi / 3) * (hi**3 - lo**3)
            nblocks = np.count_nonzero((weights * histogram[:, k]).reshape(-1, BLOCK_SIZE).sum(axis=1))
            position_rows.append([temperature, name, lo, hi, val[j], errors[j], val[j] / width,
                                  errors[j] / width, val[j] / volume, errors[j] / volume,
                                  (val[j] - old_v[j]) / width, paired[j] / width, nblocks, lo >= 0.48])
        if not sampled_overlap:
            info['momentum_characteristic_111'] = {'status': 'not evaluated: closed-path weight participation too low'}
            case_info['variants'][name] = info
            continue
        char = characteristic(packed, temperature, local, logw)
        if name == 'old':
            poly_error = float(np.max(np.abs(char.blocks - old_poly.blocks)))
            case_info['old_characteristic_field_callback_vs_polynomial_max_abs_block_error'] = poly_error
            if poly_error > 1e-12:
                raise ValueError(f'{case}: old field-aware callback and archived polynomial estimator disagree')
            old_char = char
        anchor = int(np.flatnonzero(DISPLACEMENTS == WITNESS_A)[0])
        direct, outside_share = direct_witness(x, field, old_bead, local, logw, beta)
        witness_error = float(abs(direct - char.value[anchor]))
        if char.frames_used != len(x) or abs(char.value[0] - 1) > 1e-14 or witness_error > 1e-12:
            raise ValueError(f'{case}/{name}: displaced path normalization or shifted-target witness disagrees')
        char_paired = se(char.blocks - old_char.blocks)
        for idx, radius in enumerate(DISPLACEMENTS):
            numerator_ess = char.numerator_block_ess[idx]
            characteristic_rows.append([temperature, name, radius, char.value[idx], char.se[idx],
                                        char.value[idx] - old_char.value[idx], char_paired[idx],
                                        numerator_ess, numerator_ess >= 4])
        info['momentum_characteristic_111'] = {
            'status': 'estimated at finite P; no Fourier inversion or far-tail claim',
            'min_numerator_block_ess': float(char.numerator_block_ess.min()),
            'all_displacements_pass_numerator_block_ess_4': bool(np.all(char.numerator_block_ess >= 4)),
            'beta_over_beads_eV_inverse': char.beta_over_beads_eV_inv,
            'kinetic_factor_eV_A2': char.kinetic_factor_eV_A2,
            'direct_witness_displacement_A': WITNESS_A, 'direct_witness_abs_error': witness_error,
            'witness_numerator_share_from_any_shifted_bead_gt_0.48_A': outside_share,
        }
        case_info['variants'][name] = info
    if 'quadratic' in estimates and 'quartic' in estimates:
        primary, primary_blocks = estimates['quadratic']
        other, other_blocks = estimates['quartic']
        paired = se(other_blocks - primary_blocks)
        case_info['quartic_minus_quadratic_position'] = {
            key: {'change': float(other[j] - primary[j]), 'paired_change_block_se': paired[j]}
            for j, key in enumerate(names)
        }
    return case_info, position_rows, characteristic_rows


def table(name, columns, rows):
    with (HERE / name).open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('variants', nargs='+', choices=('quadratic', 'quartic'),
                        help='candidate well.energy_force variants; old is always included')
    args = parser.parse_args()
    if len(set(args.variants)) != len(args.variants):
        parser.error('give distinct candidate variant names')
    import well

    def local_for(variant):
        return lambda x: well.energy_force(x, variant)[0]

    local_models = {name: local_for(name) for name in args.variants}
    old = full('Pd', 'joint')
    cases = {r['case']: r for r in json.loads((G10 / 'output' / 'coupled' / 'trajectory_summary.json').read_text())}
    metadata = {
        'replay_from_repo_root': (
            'OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 '
            'python -B '
            f'g10_quantum_site/bulk-refit-048/coupled_sensitivity.py {" ".join(args.variants)}'
        ),
        'well_py_sha256': sha256(Path(well.__file__).read_bytes()).hexdigest(),
        'scope': 'tagged-site-only sensitivity in the full PdD cell; NOT an all-cell refit',
        'sampling': 'source joint old homogeneous PdD primitive paths; mass, graph, box and host Hamiltonian held fixed; host distribution may respond through joint reweighting',
        'weight': 'one frame weight exp[-beta mean_beads(V_new(x)-V_old(x))] on all retained beads; neighbor and spring terms cancel; archived all-cell alternate weights are not used',
        'position': 'direct joint bead marginal, f(r)=dP/dr; empty bins mean unobserved, not physical zeros; shell-averaged spatial density is shell probability divided by shell volume',
        'support': '0.24 Å is the old nominal directional extent; 0.48 Å is the new maximum signed ray radius, NOT a validated ball; even some interior off-axis points and every exterior prediction rely on continuation. 0.60/0.80 Å are diagnostics, not data limits.',
        'uncertainty': 'paired ratio influence SE across nine equal contiguous time blocks; not a formal confidence interval or a correction for remaining correlations. Frame ESS and block ESS measure sampled weight participation only. No P extrapolation, PES uncertainty or unseen-overlap guarantee.',
        'momentum': 'existing exact finite-P midpoint displaced-path characteristic in signed [111] component, with the actual co-beaded graph field, deuteron mass, four cuts and antithetic signs. Separation is conjugate to k=p/hbar; no coordinate-to-k cutoff, Gaussian momentum or dynamics.',
        'site_original_zero_based_atom_index': int(old['original_indices'][old['site']]),
        'site_slot_in_original_occupied_graph': old['site'],
        'occupied_sites': len(old['K']), 'mass_amu': MASS['D'], 'old_source': old['source'],
        'old_stiffness_eV_A2': old['K'][old['site']].tolist(),
        'old_F0_eV_A': old['f'][old['site']].tolist(),
        'old_a_b_eV_A4': [old['a'], old['b']], 'cases': [],
    }
    position, momentum = [], []
    for temperature in (300, 77):
        result, p, m = analyze(temperature, local_models, old, local_for('old'), cases[f'PdD_{temperature}K'])
        metadata['cases'].append(result)
        position.extend(p)
        momentum.extend(m)
    warm, cold = (case['variants'] for case in metadata['cases'])
    metadata['warm_minus_cold_position'] = {}
    for name in warm:
        by_metric = {}
        for metric in ('MSD_A2', 'bead_fraction_gt_0.24_A', 'bead_fraction_gt_0.48_A'):
            w, c = (case[name]['position'][metric] for case in (warm, cold))
            change = float(w['value'] - c['value'])
            error_scale = float(w['block_se'] + c['block_se'])
            by_metric[metric] = {'change': change, 'sum_of_block_se': error_scale,
                                 'warm_gt_cold_by_3_times_sum_of_se': change > 3 * error_scale}
        metadata['warm_minus_cold_position'][name] = by_metric
    table('coupled_position.csv', ['temperature_K', 'variant', 'r_lo_A', 'r_hi_A', 'shell_probability',
          'shell_probability_block_se', 'f_per_A', 'f_block_se_per_A', 'shell_mean_density_per_A3',
          'shell_mean_density_block_se_per_A3', 'f_change_vs_old_per_A', 'paired_f_change_block_se_per_A',
          'nonzero_numerator_blocks', 'wholly_beyond_max_measured_ray'], position)
    table('coupled_characteristic_111.csv', ['temperature_K', 'variant', 'separation_A', 'characteristic',
          'block_se', 'change_vs_old', 'paired_change_block_se', 'numerator_block_ess',
          'passes_numerator_block_ess_4'], momentum)
    (HERE / 'coupled_sensitivity.json').write_text(json.dumps(metadata, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
