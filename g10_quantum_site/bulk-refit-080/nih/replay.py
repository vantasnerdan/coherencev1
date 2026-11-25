"""Replay selected NiH isolated/finite-P observables without rerunning trajectories."""
import argparse
from hashlib import sha256
import importlib
import json
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from nih_well import DIRECTIONS, NiHWell, SIGNED_SHA256, VARIANTS, load_anchors

HERE = Path(__file__).resolve().parent
TEMPERATURES = (77, 300)
MODES = ('msd_A2', 'kinetic_eV', 'probability_r_gt_048', 'probability_r_gt_080')
CUTOFF_EV = 1.05  # same finite-basis energy ceiling as the archived Galerkin consumer


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def json_write(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def common(temperature, variant, values):
    return dict(material='NiH', isotope='H', temperature_K=temperature, variant=variant,
                **dict(zip(MODES, values)))


def load_reference(root):
    g10 = root / 'g10_quantum_site'
    sys.path[:0] = [str(g10 / 'bulk-refit-048'), str(g10)]
    return tuple(importlib.import_module(name) for name in ('model', 'exact', 'solve'))


def old_spectrum(old_exact, well, mass, npar):
    raw = old_exact.solve(well.k, well.a, well.b, mass, npar)
    return dict(raw, isotope='H', variant='old', order=None, legacy_spectrum=raw,
                blocks=[dict(par=par, mult=mult, energies=e[e <= CUTOFF_EV], vectors=v[:, e <= CUTOFF_EV])
                        for par, mult, e, v in raw['blocks']])


def matched_spectrum(solver, well, mass, length, variant, npar, order):
    # The read-only .48 solver has a fixed well and oscillator scale at module level.
    # Substitute both providers only during this call; nothing in NiH uses the PdD scale.
    def energy(points, name):
        return well.energy_force(points, name, gradient=False)

    with patch.multiple(solver, energy_force=energy, length_scale=lambda m: length, MASS={'H': mass}):
        return solver.solve('H', variant, npar, order)


def signed_plane(grid, rho):
    n = DIRECTIONS[2]
    a = np.array([1., -1., 0.]) / np.sqrt(2)
    b = np.cross(n, a)
    s = np.linspace(-1.0, 1.0, 201)
    t = np.linspace(-1.2, 1.2, 121)
    aa, bb = np.meshgrid(t, t, indexing='ij')
    plane = aa[..., None] * a + bb[..., None] * b
    interp = RegularGridInterpolator((grid,) * 3, rho, bounds_error=False, fill_value=0.)
    density = [np.trapezoid(np.trapezoid(interp(abs(z * n + plane)), t, axis=1), t) for z in s]
    return s, np.asarray(density)


def distribution_data(solver, spectrum, temperature, fine, full=True):
    ngrid, angles = (101, 24) if fine else (71, 16)
    r = np.linspace(0., 1.2, 301)
    xgrid = np.linspace(0., r[-1], ngrid)
    rho, omitted = solver.density(spectrum, temperature, xgrid)
    pos = solver.radial_distribution(spectrum, temperature, r, xgrid, rho, angular=angles)
    msd, kinetic = solver.moments(spectrum, temperature)
    values = [msd, kinetic, *np.interp([.48, .80], r, pos['tail']).tolist()]
    e0 = float(spectrum['blocks'][0]['energies'][0])
    support = dict(npar=spectrum['npar'], quadrature_order=spectrum['order'], length_A=spectrum['length'],
                   ground_energy_eV=e0, position_grid_pitch_A=float(xgrid[1]), angular_GL_order=angles,
                   position_window_A=float(r[-1]), position_radial_norm=float(pos['tail'][0]),
                   omitted_active_thermal_weight=omitted,
                   finite_basis_energy_ceiling_thermal_weight_bound=float(
                       8 * spectrum['npar'] ** 3 * np.exp(-(CUTOFF_EV - e0) / (solver.KB * temperature))))
    if abs(support['position_radial_norm'] - 1) > 3e-4 or np.min(pos['density']) < -1e-16:
        raise ValueError(f'Isolated position integration failed: {spectrum["variant"]}/{temperature}')
    if not full:
        return values, support, {}
    k = np.linspace(0., 48., 241)
    kgrid = np.linspace(0., k[-1], ngrid)
    nk, _ = solver.density(spectrum, temperature, kgrid, momentum=True)
    mom = solver.radial_distribution(spectrum, temperature, k, kgrid, nk, momentum=True, angular=angles)
    s, plane = signed_plane(xgrid, rho)
    support.update(momentum_grid_pitch_A_inverse=float(kgrid[1]), momentum_window_A_inverse=float(k[-1]),
                   momentum_radial_norm=float(mom['tail'][0]), momentum_min_radial_density_A3=float(mom['density'].min()),
                   signed_111_marginal_norm=float(np.trapezoid(plane, s)),
                   position_radial_msd_A2=float(np.trapezoid(r * r * pos['measure'], r)),
                   radial_momentum_kinetic_eV=float(solver.kinetic(spectrum['mass']) * np.trapezoid(k * k * mom['measure'], k)))
    if (abs(support['momentum_radial_norm'] - 1) > 3e-4 or np.min(mom['density']) < -1e-16
            or abs(support['signed_111_marginal_norm'] - 1) > 3e-4):
        raise ValueError(f'Isolated momentum or plane integration failed: {spectrum["variant"]}/{temperature}')
    arrays = dict(position_r_A=r, position_spherical_density_A_inverse3=pos['density'],
                  position_radial_measure_A_inverse=pos['measure'], position_tail_to_window=pos['tail'],
                  position_line_density_A_inverse3=pos['lines'], momentum_k_A_inverse=k,
                  momentum_spherical_density_A3=mom['density'], momentum_radial_measure_A=mom['measure'],
                  momentum_tail_to_window=mom['tail'], momentum_line_density_A3=mom['lines'],
                  position_signed_111_s_A=s, position_signed_111_plane_marginal_A_inverse=plane)
    return values, support, arrays


def paired_tagged(root, well, model):
    # The neighbour field is from the same bead/frame. The unchanged springs,
    # other-site on-sites and all pair terms cancel in the closed-path ratio.
    g10 = root / 'g10_quantum_site'
    pathdir = g10 / 'output' / 'coupled'
    archive = {r['case']: r for r in json.loads((pathdir / 'trajectory_summary.json').read_text())}
    rows, support = [], {}
    for t in TEMPERATURES:
        case = f'NiH_{t}K'
        path = pathdir / f'{case}_paths.npz'
        with np.load(path) as data:
            x, field = data['beads_A'], data['neighbor_field_eV_A']
        beads = {77: 128, 300: 32}[t]
        if x.shape != (801, beads, 3) or field.shape != x.shape or not (np.isfinite(x).all() and np.isfinite(field).all()):
            raise ValueError(f'{case}: invalid saved primitive paths or field')
        r2 = np.einsum('fbi,fbi->fb', x, x)
        msd = r2.mean(axis=1)
        archive_msd_error = abs(msd[-792:].mean() - archive[case]['mean_square_radius_A2'])
        archive_p024_error = abs(np.mean(r2 > .24**2) - archive[case]['position_outside_support_fraction'])
        if max(archive_msd_error, archive_p024_error) > 1e-12:
            raise ValueError(f'{case}: wrong path population or archived comparison')
        tail48, tail80 = [(r2 > r**2).mean(axis=1) for r in (.48, .80)]
        thermal_eV = model.KB * t
        xc = x.mean(axis=1, keepdims=True)
        fluct = x - xc
        old_e = well.energy_force(x, 'old', gradient=False)[0]
        angles = {}
        block_count, block_size = 9, 89
        info = dict(path_relative_to_reference=str(path.relative_to(root)), path_sha256=digest(path),
                    beads=beads, frames=len(x), contiguous_blocks=block_count, frames_per_block=block_size,
                    max_bead_radius_A=float(np.sqrt(r2.max())),
                    archive_old_msd_abs_error_A2=float(archive_msd_error), archive_old_p024_abs_error=float(archive_p024_error),
                    frames_with_sampled_tail_048=int(np.count_nonzero(tail48)),
                    frames_with_sampled_tail_080=int(np.count_nonzero(tail80)),
                    kinetic_estimator='tagged primitive centroid virial with the full saved co-beaded neighbor field',
                    momentum_distribution=None, variants={})
        for variant in VARIANTS[:-1]:
            new_e, new_f = well.energy_force(x, variant)
            energy_diff = new_e - old_e
            logw = -energy_diff.mean(axis=1) / thermal_eV
            w = np.exp(logw - logw.max())
            den = w.reshape(block_count, block_size).sum(axis=1)
            kinetic = 1.5 * thermal_eV - .5 * np.einsum('fbi,fbi->f', fluct, new_f + field) / beads
            observables = np.column_stack([msd, kinetic, tail48])
            num = (w[:, None] * observables).reshape(block_count, block_size, -1).sum(axis=1)
            estimates = num.sum(axis=0) / den.sum()
            influence = estimates + (num - den[:, None] * estimates) / den.mean()
            angles[variant] = estimates, influence
            se = influence.std(axis=0, ddof=1) / np.sqrt(block_count)
            paired = (influence - angles['old'][1]).std(axis=0, ddof=1) / np.sqrt(block_count)
            frame_ess, block_ess = w.sum()**2 / (w @ w), den.sum()**2 / (den @ den)
            tail_blocks = int(np.count_nonzero(num[:, 2]))
            ok = frame_ess >= 50 and block_ess >= 4
            # Verify the force-based virial against the actual differentiated site action,
            # scaling only tagged bead fluctuations and retaining the full pair weight once.
            x0, d0, f0 = x[:3], fluct[:3], field[:3]
            h = 1e-5
            ep = well.energy_force(x0 + h * d0, variant, gradient=False)[0]
            em = well.energy_force(x0 - h * d0, variant, gradient=False)[0]
            fd = (ep - em) / (2 * h) - np.einsum('fbi,fbi->fb', f0, d0)
            direct = -np.einsum('fbi,fbi->fb', (new_f + field)[:3], d0)
            virial_error = float(np.max(abs(fd - direct)))
            if virial_error > 1e-8:
                raise ValueError(f'{case}/{variant}: tagged virial differentiation disagrees')
            supported = [ok, ok, ok and tail_blocks >= 4]
            published = [float(v) if use else None for v, use in zip(estimates, supported)]
            rows.append(common(t, 'finite_P_tagged_' + variant, published + [None]))
            info['variants'][variant] = dict(frame_weight_ess=float(frame_ess), block_weight_ess=float(block_ess),
                                              largest_frame_weight_fraction=float(w.max() / w.sum()),
                                              log_weight_span=float(np.ptp(logw)), tail048_contributing_blocks=tail_blocks,
                                              tagged_virial_fd_max_abs_error_eV=virial_error,
                                              paired=dict(zip(MODES, [dict(block_se=float(err),
                                                  change_vs_old=float(v - base), paired_change_block_se=float(p))
                                                  if use else None for err, v, base, p, use in
                                                  zip(se, estimates, angles['old'][0], paired, supported)] + [None])))
        if all(row['probability_r_gt_048'] is not None for row in rows[-3:]):
            diff = angles['w'][0] - angles['w2'][0]
            err = (angles['w'][1] - angles['w2'][1]).std(axis=0, ddof=1) / np.sqrt(block_count)
            info['w_minus_w2'] = dict(zip(MODES, [dict(change=float(d), paired_block_se=float(e))
                                                   for d, e in zip(diff, err)] + [None]))
        support[str(t)] = info
    return rows, support


def archived_grid_control(root, old_exact, reference, model):
    # Compare like with like: the archive uses a coarser Cartesian grid and
    # linear interpolation. It is not the exact operator-moment oracle.
    support = {}
    for t in TEMPERATURES:
        path = root / 'g10_quantum_site' / 'output' / 'exact' / f'NiH_{t}K_joint.npz'
        with np.load(path) as saved:
            x, r = saved['x_A'], saved['radial_A']
            k, kr = saved['k_A_inverse'], saved['momentum_k_A_inverse']
            rho, _ = old_exact.density_grid(reference['legacy_spectrum'], t, x)
            nk, _ = old_exact.density_grid(reference['legacy_spectrum'], t, k, momentum=True)
            position, _ = old_exact.distributions(x, rho, .24, r)
            momentum, _ = old_exact.distributions(k, nk, None, kr)
            dx = float(np.max(abs(position['density'] - saved['radial_density_A_minus3'])))
            dk = float(np.max(abs(momentum['density'] - saved['momentum_density_A3'])))
            old_saved = [float(np.trapezoid(r**2 * saved['radial_measure_A_minus1'], r)),
                         float(model.kinetic_coefficient(reference['mass'])
                               * np.trapezoid(kr**2 * saved['momentum_measure_A'], kr)),
                         *np.interp([.48, .80], r, saved['radial_tail']).tolist()]
        if max(dx, dk) > 3e-11:
            raise ValueError(f'{path.name}: prior NiH Hamiltonian does not reproduce its archived grid')
        support[t] = dict(sha256=digest(path), same_grid_position_density_max_abs_error_A_inverse3=dx,
                          same_grid_momentum_density_max_abs_error_A3=dk, archived=dict(zip(MODES, old_saved)))
    return support


def run(root, signed):
    model, exact, solver = load_reference(root)
    old = model.full('Ni', 'joint')
    site = old['site']
    moving_index = int(old['original_indices'][site])
    anchors = load_anchors(signed, moving_index)
    well = NiHWell(float(old['K'][site, 0, 0]), old['a'], old['b'], anchors)
    probes = well.probe()
    coupled, path_support = paired_tagged(root, well, model)
    print('NiH signed-well and tagged-path checks complete', flush=True)
    archived_reference = old_spectrum(exact, well, model.MASS['H'], 10)
    archive_support = archived_grid_control(root, exact, archived_reference, model)
    reference = old_spectrum(exact, well, model.MASS['H'], 12)
    rows, support, arrays, results = [], {}, {}, {}
    for variant in VARIANTS:
        if variant == 'old':
            spectrum = reference
        else:
            spectrum = matched_spectrum(solver, well, reference['mass'], reference['length'], variant, 12, 52)
        for t in TEMPERATURES:
            values, meta, table = distribution_data(solver, spectrum, t, fine=True)
            results[variant, t] = values
            key = f'{variant}_{t}K'
            support[key] = meta
            arrays.update({key + '_' + k: array for k, array in table.items()})
            print('NiH isolated', key, values, flush=True)
    # One bounded combined basis/quadrature/readout change per angular closure.
    # No new exterior solve here: the transport compares both exteriors on the same fine rule.
    for variant in ('w', 'w2'):
        coarse_spectrum = matched_spectrum(solver, well, reference['mass'], reference['length'], variant, 10, 42)
        for t in TEMPERATURES:
            coarse, coarse_meta, _ = distribution_data(solver, coarse_spectrum, t, fine=False, full=False)
            fine = results[variant, t]
            support[f'{variant}_{t}K']['numeric_change'] = dict(coarse_settings=coarse_meta,
                                                                fine_minus_coarse=dict(zip(MODES, np.subtract(fine, coarse).tolist())))
    for t in TEMPERATURES:
        current = results['old', t]
        archived = archive_support[t]
        archived['operator_and_fine_readout_minus_archived'] = dict(zip(MODES, np.subtract(
            current, list(archived['archived'].values())).tolist()))
        old_msd, old_ke = solver.moments(archived_reference, t)
        archived['old_npar_12_minus_10_operator_moments'] = dict(zip(MODES[:2], np.subtract(
            current[:2], [old_msd, old_ke]).tolist()))
        support[f'old_{t}K']['archived_grid_control'] = archived
        for variant in VARIANTS:
            values = results[variant, t].copy()
            if variant != 'old':
                compare = 'w' if variant == 'w_stiff_exterior' else variant
                det = abs(support[f'{compare}_{t}K']['numeric_change']['fine_minus_coarse'][MODES[-1]])
                support[f'{variant}_{t}K']['p080_numeric_resolution_probe'] = dict(
                    raw_finite_window_value=values[-1], reference_comparison_variant=compare,
                    abs_fine_coarse_difference=det, retained=values[-1] > 3 * det)
                if values[-1] <= 3 * det:
                    values[-1] = None
            rows.append(common(t, 'isolated_' + variant, values))
    sources = [model.NI, model.NIQ, 'g2_refit/computed_inputs.json',
               'g10_quantum_site/model.py', 'g10_quantum_site/exact.py', 'g10_quantum_site/bulk-refit-048/solve.py']
    e_old, f_old = well.energy_force(np.array([r['signed_displacement_A'] for r in anchors])[:, None] * DIRECTIONS[2], 'old')
    for anchor, e, f in zip(anchors, e_old, f_old @ DIRECTIONS[2]):
        anchor.update(old_energy_eV=float(e), old_force_parallel_eV_A=float(f))
    model_info = dict(on_site_k_eV_A2=well.k, old_q100_eV_A4=well.a, old_q110_eV_A4=well.q110,
                      old_q111_eV_A4=well.q111, old_cubic_cross_coefficient_eV_A4=well.b,
                      min_spline_dE_du_eV_A2=well.min_dE_du, q110_measured=False,
                      extrapolation_curvature_eV_A2=dict(w=well.k, w2=well.k, w_stiff_exterior=8 * well.k),
                      original_moving_index=moving_index, folded_site_index=site, graph_sites=len(old['J']),
                      graph_signed_neighbor_sum_eV_A2=float(old['J'][site].sum()),
                      graph_J_binary64_le_sha256=sha256(old['J'].astype('<f8').tobytes()).hexdigest(),
                      graph_max_eigenvalue_eV_A2=float(np.linalg.eigvalsh(old['J'])[-1]),
                      physical_H_mass_amu=reference['mass'], kinetic_coefficient_eV_A2=model.kinetic_coefficient(reference['mass']),
                      matched_NiH_length_A=reference['length'],
                      source_sha256={name: digest(root / name) for name in sources},
                      signed_export_sha256=SIGNED_SHA256, signed_anchors=anchors,
                      signed_force_and_splice_probe=probes)
    contrasts = {}
    for t in TEMPERATURES:
        contrasts[str(t)] = {name: dict(zip(MODES, np.subtract(results[left, t], results[right, t]).tolist()))
                             for name, left, right in [('w_minus_old', 'w', 'old'), ('w2_minus_old', 'w2', 'old'),
                                                       ('w_minus_w2', 'w', 'w2'),
                                                       ('stiff_exterior_minus_w', 'w_stiff_exterior', 'w')]}
    json_write('support.json', dict(model=model_info, isolated=support, isolated_contrasts=contrasts, finite_P_tagged=path_support))
    json_write('observable_table.json', rows + coupled)
    np.savez_compressed(HERE / 'distributions.npz', **arrays)


if __name__ == '__main__':
    argp = argparse.ArgumentParser(description=__doc__)
    argp.add_argument('--reference-root', type=Path, default=HERE.parents[2],
                      help='checkout or root of the pinned prior G10 reference assets')
    argp.add_argument('--signed-export', type=Path,
                      default=HERE.parents[1] / '.input-080/g10_quantum_site/bulk-wide-wells/bulk-wells.json',
                      help='verified provider bulk-wells.json')
    args = argp.parse_args()
    run(args.reference_root.resolve(), args.signed_export.resolve())
