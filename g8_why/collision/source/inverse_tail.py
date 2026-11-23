"""Finite-shell positive inverse for the fixed-profile free-target diagnostic.

Run with one CPU/thread. Reported minima are only over the explicitly stored
finite shell set, NOT universal bounds for continuous momentum distributions.
"""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
import uuid
import warnings

import numpy as np
import scipy
from scipy.optimize import linprog, OptimizeWarning

from yield_kernel import DEFAULT_BASE, FixedProfileShellKernel


def cpu_seconds():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_utime+usage.ru_stime


def reporting_coefficients(kernel, material, energies, shells, current_power, p=None, n=None):
    common = kernel.common
    p = kernel.matrix(material, energies, shells, 'p') if p is None else p
    n = kernel.matrix(material, [6.5], shells, 'n')[0] if n is None else n
    join = common.yield_C(material, 4.5, 'n')/common.yield_C(material, 4.5, 'p')*1.0413535612834366
    collection = common.DD.collected_fraction(np.asarray(energies), 1., 0.)
    a = p*join/(collection*common.DD.I_cup(np.asarray(energies))**current_power)[:, None]
    b = n/common.DD.collected_fraction(6.5, 1., 0.)/common.DD.I_cup(6.5)**current_power
    return a, b


def conditional_shell_weights(a, b, observed, shells):
    rows = []
    for j, target in enumerate(shells):
        if target == 0:
            continue
        for i, measured in enumerate(observed):
            denominator = (a[i, j]-a[i, 0])-measured*(b[j]-b[0])
            f = (measured*b[0]-a[i, 0])/denominator if denominator != 0 else None
            rows.append(dict(energy_index=i, shell_target_eV=float(target),
                             conditional_shell_fraction=f,
                             disposition='no_solution' if f is None else
                             'feasible_single_energy_only' if 0 <= f <= 1 else
                             'impossible_fraction_above_one' if f > 1 else 'impossible_negative_fraction'))
    return rows

def pairwise_infeasibility(a, observed, sigma_ln, band_sigma):
    """A normalized positive mixture cannot leave any shell ratio envelope."""
    lower = observed*np.exp(-band_sigma*sigma_ln)
    upper = observed*np.exp(band_sigma*sigma_ln)
    certificates = []
    for i in range(len(observed)):
        for j in range(i+1, len(observed)):
            ratios = a[i]/a[j]
            possible = [float(ratios.min()), float(ratios.max())]
            required = [float(lower[i]/upper[j]), float(upper[i]/lower[j])]
            if required[0] > possible[1] or required[1] < possible[0]:
                certificates.append(dict(numerator_energy_index=i, denominator_energy_index=j,
                                         shell_ratio_range=possible, required_ratio_range=required))
    return certificates



def positive_inverse(a, b, observed, sigma_ln, shells, cap_eV, band_sigma, threshold_eV):
    keep = shells <= cap_eV
    support = shells[keep]
    aa, bb = a[:, keep], b[keep]
    certificates = pairwise_infeasibility(aa, observed, sigma_ln, band_sigma)
    if certificates:
        return dict(cap_target_eV=float(cap_eV), band_sigma=float(band_sigma),
                    strict_tail_threshold_eV=float(threshold_eV), support_count=int(len(support)),
                    status=2, disposition='finite_grid_infeasible',
                    solver_message='Positive-mixture shell-ratio envelope excludes requested ratios; no solver required',
                    pairwise_certificates=certificates,
                    domain='nonnegative normalized mixtures on the exact stored shell grid <= cap; strict T > threshold')
    center = observed[:, None]*bb
    options = dict(primal_feasibility_tolerance=1e-9, dual_feasibility_tolerance=1e-9, threads=1)
    if band_sigma == 0:
        difference = aa-center
        scale = np.max(np.abs(difference), axis=1)
        a_eq = np.vstack([np.ones(len(support)), difference/scale[:, None]])
        b_eq = np.r_[1., np.zeros(len(observed))]
        a_ub = b_ub = None
    else:
        lower = observed*np.exp(-band_sigma*sigma_ln)
        upper = observed*np.exp(band_sigma*sigma_ln)
        difference = np.vstack([aa-upper[:, None]*bb, lower[:, None]*bb-aa])
        scale = np.max(np.abs(difference), axis=1)
        a_ub = difference/scale[:, None]
        b_ub = np.zeros(len(difference))
        a_eq, b_eq = np.ones((1, len(support))), np.ones(1)
    objective = (support > threshold_eV).astype(float)
    # SciPy forwards the native HiGHS threads=1 option; suppress only its
    # expected Python-wrapper warning about that native option.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Unrecognized options detected.*threads', category=OptimizeWarning)
        fit = linprog(objective, A_ub=a_ub, b_ub=b_ub, A_eq=a_eq, b_eq=b_eq,
                      bounds=(0., None), method='highs', options=options)
    result = dict(cap_target_eV=float(cap_eV), band_sigma=float(band_sigma),
                  strict_tail_threshold_eV=float(threshold_eV), support_count=int(len(support)),
                  status=int(fit.status), solver_message=fit.message,
                  domain='nonnegative normalized mixtures on the exact stored shell grid <= cap; strict T > threshold')
    if fit.success:
        weights = np.maximum(fit.x, 0.)
        weights /= weights.sum()
        ratios = aa@weights/(bb@weights)
        residual = np.log(ratios/observed)
        violation = float(np.max(np.maximum(0., np.abs(residual)-band_sigma*sigma_ln)))
        result.update(disposition='finite_grid_minimum' if violation < 1e-5 else 'numerically_unresolved_witness',
                      finite_grid_minimum_tail_weight=float(objective@weights),
                      witness=[dict(target_eV=float(t), probability=float(w)) for t, w in zip(support, weights) if w > 0],
                      witness_reported_R=ratios.tolist(), witness_log_residual=residual.tolist(),
                      max_log_band_violation=violation)
    else:
        result['disposition'] = 'finite_grid_infeasible' if fit.status == 2 else 'solver_unresolved'
    return result


def run(output_dir, base_root=DEFAULT_BASE):
    start_cpu = cpu_seconds()
    start_wall = time.monotonic()
    resource.setrlimit(resource.RLIMIT_CPU, (1500, 1500))
    kernel = FixedProfileShellKernel(base_root, angular_order=64)
    common = kernel.common
    energies = common.ENERGIES
    # The 1eV geometric grid resolves cold-to-fast shells; named exact supports
    # include threshold and cap boundaries so the strict > definition is explicit.
    named = np.array([0., 1., 10., 30., 100., 300., 1000., 3000., 10000., 30000., 100000.])
    shells = np.unique(np.r_[named, np.geomspace(1., 100000., 121)])
    caps = [100., 300., 1000., 3000., 10000., 30000., 100000.]
    thresholds = [100., 300., 1000.]
    report = dict(schema='g8-free-target-shell-inverse-v1',
                  applicability='ON-SHELL FREE-TARGET KINEMATIC DIAGNOSTIC ONLY; neither bound-state prediction nor rigorous bound-state upper bound',
                  energies_nominal_CM_keV=energies.tolist(), stationary_target_incident_D1_LAB_keV=(2*energies).tolist(),
                  target_shell_LAB_kinetic_eV=shells.tolist(), angular_order=kernel.angular_order,
                  normalization=dict(physical='reactions per true charged C, fixed transported nD/path/sheet/stopping; Y/Y40 divides physical same-channel stationary yield',
                                     paper='common.readout: own 6.5keV neutron reference; fixed stationary smooth40 n/p4.5 stitch times 1.0413535612834366; b_hi=1,b_lo=0; k0 primary,k1 current-power alternative',
                                     beam='85% D+ and 15% D2+ electrical current; multiplicity 1,2; 5% FWHM five-Hermite-node spread; unchanged stopping/profile at nominal CM key',
                                     relative_flux='sigma(E_relative_CM,+0.040keV)*sqrt(2*E_relative_CM/LAB_beam_keV)',
                                     target='same normalized isotropic target distribution in every layer, sheet, and nominal energy; no density or S-factor adjustment'),
                  source_bindings=kernel.source_bindings(), materials={}, findings=[], final_checks={})
    for material in common.MATERIALS:
        p = kernel.matrix(material, energies, shells, 'p')
        n = kernel.matrix(material, energies, shells, 'n')
        obs, err = common.observed(material, energies)
        material_result = dict(kernel_p_reactions_per_true_C=p.tolist(), kernel_n_reactions_per_true_C=n.tolist(),
                               physical_p_Y_over_Y40=(p/p[:, :1]).tolist(), physical_n_Y_over_Y40=(n/n[:, :1]).tolist(),
                               observed_R=obs.tolist(), observed_sigma_lnR=err.tolist(), inverse={})
        for power in (0, 1):
            a, b = reporting_coefficients(kernel, material, energies[:4], shells, power, p[:4], n[-1])
            solutions = [positive_inverse(a, b, obs[:4], err[:4], shells, cap, band, threshold)
                         for cap in caps for band in (0., 1., 2.) for threshold in thresholds]
            physical_shape = obs[0]/obs[3]*common.DD.I_cup(.25)**power/common.DD.I_cup(2.5)**power
            material_result['inverse'][f'k{power}'] = dict(
                reported_numerator_coefficients=a.tolist(), reported_neutron_reference_coefficients=b.tolist(),
                conditional_cold_plus_one_shell=conditional_shell_weights(a, b, obs[:4], shells),
                positive_distribution_solutions=solutions,
                physical_p_shape_025_over_25_required=float(physical_shape),
                finite_grid_physical_p_shape_range=[float(np.min(p[0]/p[3])), float(np.max(p[0]/p[3]))])
        moments = []
        for e in energies:
            for t in [0., 100., 300., 1000., 3000., 10000., 30000., 100000.]:
                for angle in (0., 90.):
                    moments.append(dict(nominal_CM_keV=float(e), target_kinetic_eV=t,
                                        **kernel.shell_moments(material, e, t, detector_deg=angle)))
        material_result['reaction_weighted_COM_and_nominal_proton_moments'] = moments
        report['materials'][material] = material_result
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    persisted_path = output_dir/f'inverse-distributions-{uuid.uuid4().hex}.json'
    persisted_bytes = (json.dumps(dict(target_shell_LAB_kinetic_eV=shells.tolist(),
                                      materials=report['materials']), allow_nan=False)+'\n').encode()
    with persisted_path.open('xb') as stream:
        stream.write(persisted_bytes)
    persisted = json.loads(persisted_path.read_bytes())
    report['persisted_inverse'] = dict(path=str(persisted_path), bytes=len(persisted_bytes),
                                       sha256=hashlib.sha256(persisted_bytes).hexdigest())
    # One narrow final check: actual persisted inverse witness goes back through
    # its public physical-yield and unchanged paper readout consumer, plus T=0
    # through the authoritative forward_yields consumer (not a fixture).
    stationary = 0.
    for material in common.MATERIALS:
        for channel in ('p', 'n'):
            actual = common.forward_yields(energies, common.reference(material).profiles,
                Ue_keV=common.U, measurement={'channel': channel}, reaction='eq7',
                stopping_fn=common.transport_stopping, reference_E_cm_keV=None)['physical_yield_per_C']
            own = kernel.matrix(material, energies, [0.], channel)[:, 0]
            stationary = max(stationary, float(np.max(np.abs(own/actual-1))))
    report['final_checks']['stationary_shared_forward_consumer_max_relative_error'] = stationary
    assert stationary < 3e-13
    witness_checks = []
    for material, values in persisted['materials'].items():
        for power in (0, 1):
            candidates = [r for r in values['inverse'][f'k{power}']['positive_distribution_solutions']
                          if r['disposition'] == 'finite_grid_minimum']
            if not candidates:
                witness_checks.append(dict(material=material, current_power=power, disposition='no_feasible_grid_witness'))
                continue
            chosen = candidates[0]
            support = [r['target_eV'] for r in chosen['witness']]
            weights = [r['probability'] for r in chosen['witness']]
            direct = kernel.reported(material, energies, support, weights, power)
            error = float(np.max(np.abs(direct[:4]/chosen['witness_reported_R']-1)))
            assert error < 3e-12
            witness_checks.append(dict(material=material, current_power=power,
                max_relative_error=error, reported_R_consumer=direct.tolist(),
                cap_target_eV=chosen['cap_target_eV'], band_sigma=chosen['band_sigma'], witness=chosen['witness']))
    report['final_checks']['positive_inverse_public_consumer'] = witness_checks
    # A single selected angular quadrature sensitivity, not refinement of any
    # well/harmonic model. Include near cancellation and high-energy shell.
    fine = FixedProfileShellKernel(base_root, angular_order=128)
    quadrature = []
    for e in energies:
        for t in [100., 300., 1000., 3000., 100000.]:
            low = kernel.shell('PdD', e, t)
            high = fine.shell('PdD', e, t)
            quadrature.append(dict(material='PdD', nominal_CM_keV=float(e), target_eV=t,
                                   relative_change_128_over_64=float(high/low-1)))
    report['final_checks']['angular_64_vs_128_selected_PdD'] = quadrature
    report['final_checks']['angular_max_absolute_relative_change_checked'] = max(abs(r['relative_change_128_over_64']) for r in quadrature)
    report['findings'] = [
        dict(finding='Free-target mapping and relative flux are explicit; stationary limit reproduces actual persisted transported consumer.', disposition='implemented_and_checked'),
        dict(finding='The full physical p/n shell matrix and Y/Y40 at all five nominal CM settings are retained independently of the paper normalization.', disposition='delivered'),
        dict(finding='Each conditional cold-plus-shell fraction solves ONE reported ratio including the moving neutron denominator; different energy requirements need not share a fraction.', disposition='diagnostic_not_unique_distribution'),
        dict(finding='Positive-distribution inverses use only the four low-energy proton ratios and normalized 6.5 neutron anchor, independently for each material and k convention. Log-error bands are pointwise boxes, not confidence regions.', disposition='finite_grid_identifiability_only'),
        dict(finding='Minima above100/300/1000eV apply ONLY to nonnegative normalized weights on the stored finite grid under each stated cap. They are not continuum-distribution lower bounds; feasibility can improve with additional support.', disposition='no_unqualified_lower_bound'),
        dict(finding='A normalized measured shape does not establish an absolute physical excess or its density; no detector efficiency, target density, site fraction or S(E) was fitted.', disposition='normalizations_separated'),
        dict(finding='A bound-state n(p) cannot automatically be interpreted as on-shell kinetic energy; binding, final excitations, recoil and energy conservation are absent from this diagnostic.', disposition='root_bound_collision_operator_required'),
        dict(finding='Doppler moments assume isotropic target azimuth, CM proton direction and a fixed nominal3.020MeV proton; they omit aperture, angular nuclear distribution, detector response and energy-dependent product kinematics.', disposition='conditional_line_signature_only'),
        dict(finding='No Native, SCF, TDDFT, harmonic calculation, profile optimization or nuclear fitting was performed.', disposition='scope_preserved')]
    for material, values in report['materials'].items():
        for power in (0, 1):
            inverse = values['inverse'][f'k{power}']
            cap_solutions = [r for r in inverse['positive_distribution_solutions']
                             if r['cap_target_eV'] == 100000. and r['strict_tail_threshold_eV'] == 100.]
            report['findings'].append(dict(
                material=material, current_power=power,
                finding='Finite-grid positive-distribution feasibility through100keV target kinetic energy',
                central=next(r['disposition'] for r in cap_solutions if r['band_sigma'] == 0.),
                one_sigma_box=next(r['disposition'] for r in cap_solutions if r['band_sigma'] == 1.),
                two_sigma_box=next(r['disposition'] for r in cap_solutions if r['band_sigma'] == 2.),
                required_physical_p_shape_025_over_25=inverse['physical_p_shape_025_over_25_required'],
                shell_physical_p_shape_range=inverse['finite_grid_physical_p_shape_range'],
                disposition='finite_support_only_not_bound_state'))
    report['findings'].append(dict(
        finding='Experimental values and log-error estimates are the pinned errbars.csv digitization. The original experimental paper was not independently re-digitized or consumed in this inverse slice.',
        disposition='explicit_primary_experimental_source_gap'))
    report['Doppler_formula'] = dict(
        reaction_average='<F>_R = sum_j path_weight_j integral dOmega_t/(4pi) n(T) sigma(Erel) vrel/vbeam F / Y',
        COM='V=(vb+vt)/2; <Vx>=<Vy>=0; <Vx^2>=<Vy^2>=<vt^2(1-mu^2)/8>_R',
        exact_nonrel_symbolic='Ep_lab=Ep_star(Erel)+mp*V^2/2+sqrt(2*mp*Ep_star(Erel)) n_star dot V; average this expression and its square with the SAME reaction weights',
        fixed_Ep_linear='mean shift = mp<V^2>_R/2 + p_star cos(theta)<Vz>_R; leading Doppler variance = p_star^2[cos(theta)^2 Var_R(Vz)+sin(theta)^2<Vx^2>_R]',
        full_fixed_Ep_variance='Var_R(mp*V^2/2+p_star*cos(theta)*Vz) + p_star^2*sin(theta)^2 <Vperp^2>_R/2; this full expression is used numerically')
    report['software'] = dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                              constants='scipy.constants atomic_mass,electron_volt,proton_mass; deuteron mass is unchanged DD.M_D atomic-mass units')
    for name in ['yield_kernel.py', 'inverse_tail.py']:
        path = Path(__file__).with_name(name)
        data = path.read_bytes()
        report['source_bindings']['files'].append(dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
    archive = Path(base_root).parent/'sources/coherencev1-2318593.tar.gz'
    data = archive.read_bytes()
    report['source_bindings']['archive'] = dict(path=str(archive), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    profile = json.loads((Path(base_root)/'g4_beam/merit/profiles.json').read_bytes())
    report['source_bindings']['profile_record_count'] = len(profile['records'])
    report['source_bindings']['profile_top_level_metadata'] = {key: value for key, value in profile.items() if key != 'records'}
    report['resource'] = dict(measured_process_CPU_seconds=cpu_seconds()-start_cpu,
                              measured_process_elapsed_seconds=time.monotonic()-start_wall,
                              process_CPU_seconds_at_report=cpu_seconds(),
                              assigned_CPU_cap_seconds=1500, affinity='taskset -c7', numerical_threads=1)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir/f'inverse-report-{time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())}-{uuid.uuid4().hex}.json'
    data = (json.dumps(report, indent=2, allow_nan=False)+'\n').encode()
    with path.open('xb') as stream:
        stream.write(data)
    digest = hashlib.sha256(data).hexdigest()
    manifest = dict(report_path=str(path), bytes=len(data), sha256=digest,
                    measured_process_CPU_seconds=cpu_seconds()-start_cpu)
    with path.with_suffix('.manifest.json').open('x') as stream:
        json.dump(manifest, stream, indent=2)
        stream.write('\n')
    print(json.dumps(manifest))
    return path


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('<PRIVATE_PATH>'))
    parser.add_argument('--base-root', type=Path, default=DEFAULT_BASE)
    args = parser.parse_args()
    run(args.output_dir, args.base_root)
