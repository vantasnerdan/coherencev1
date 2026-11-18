"""Compare declared static profiles without conflating fits or PAW activation."""
import numpy as np

from g4_beam.screen.dynamic.largeq_bank import require


def _grid(case, n):
    matches = [q for q in case['quadratures'] if q['n'] == n]
    require(len(matches) == 1, 'Expected exactly one complete declared q quadrature')
    return matches[0]


def _profile(cases, key, n):
    return np.asarray(_grid(cases[key], n)['primary_U_eV'])


def _difference(cases, lhs, rhs, n):
    return _profile(cases, lhs, n) - _profile(cases, rhs, n)


def _reference_difference(cases, lhs, rhs, n):
    left = _grid(cases[lhs], n)['reference_profiles']
    right = _grid(cases[rhs], n)['reference_profiles']
    require([(p['kappa_Ainv'], p['gaussian_width_Ainv']) for p in left] ==
            [(p['kappa_Ainv'], p['gaussian_width_Ainv']) for p in right],
            'Paired reference parameters differ')
    curves = np.asarray([np.asarray(a['U_eV']) - b['U_eV'] for a, b in zip(left, right)])
    return np.max(np.abs(curves - curves[0]), axis=0)


def _assessment(reference, changes, absolute_tolerance):
    tolerance = np.maximum(absolute_tolerance, .1 * np.abs(reference))
    envelope = np.max(np.abs(np.asarray(list(changes.values()))), axis=0)
    return dict(tolerance_eV=tolerance.tolist(), sampled_change_envelope_eV=envelope.tolist(),
                all_sampled_changes_within_tolerance=bool(np.all(envelope <= tolerance)),
                changes_eV={key: value.tolist() for key, value in changes.items()},
                interpretation='Empirical sampled sensitivities, not a rigorous bound or physical model error.')


def _require_matched_geometry(left, right):
    require(left['material'] == right['material'] and
            left['geometry']['symbols'] == right['geometry']['symbols'] and
            np.allclose(left['geometry']['cell_A'], right['geometry']['cell_A'], rtol=0, atol=1e-9) and
            np.allclose(left['geometry']['positions_A'], right['geometry']['positions_A'], rtol=0, atol=1e-8) and
            np.allclose(left['source_A'], right['source_A'], rtol=0, atol=1e-8),
            'Activation or numerical control uses an unmatched physical geometry/site')


def compare_static(result):
    cases = {case['id']: case for case in result['cases']}
    required = {material + suffix for material in ('Pd', 'PdD') for suffix in (
        '18_1200', '10_500', '18_500', '18_k8', '18_1200_nb100',
        '10_500_nb96', '18_1200_G400', '10_500_G400')}
    require(set(cases) == required, 'Final comparison requires the complete registered two-host/control set')
    absolute, activation = {}, {}
    for material in ('Pd', 'PdD'):
        primary = material + '18_1200'
        for suffix in ('18_500', '18_k8', '18_1200_nb100', '18_1200_G400',
                       '10_500', '10_500_nb96', '10_500_G400'):
            _require_matched_geometry(cases[primary], cases[material + suffix])
        eighteen, ten = cases[material + '18_500'], cases[material + '10_500']
        require(eighteen['ground_state']['ecut_eV'] == ten['ground_state']['ecut_eV'] == 500 and
                eighteen['ground_state']['kmesh'] == ten['ground_state']['kmesh'] == [12, 12, 12] and
                eighteen['response']['ecut_eV'] == ten['response']['ecut_eV'] == 200 and
                eighteen['response']['eta_eV'] == ten['response']['eta_eV'] and
                eighteen['response']['nbands'] - ten['response']['nbands'] == 4,
                'Semicore activation requires matched PW/q/G/eta and nominal empty-state range')
        main3 = _profile(cases, primary, 3)
        changes = dict(
            q3_minus_q2=main3 - _profile(cases, primary, 2),
            reference_parameters=np.asarray(_grid(cases[primary], 3)['reference_parameter_max_difference_eV']),
            bands160_minus100=_difference(cases, primary, primary + '_nb100', 2),
            G200_minus400=_difference(cases, primary, primary + '_G400', 2),
            k12_minus8=_difference(cases, primary, material + '18_k8', 2),
            PW1200_minus500=_difference(cases, primary, material + '18_500', 3))
        absolute[material] = dict(U18_N3_eV=main3.tolist(),
            U18_N2_eV=_profile(cases, primary, 2).tolist(),
            empirical_numerics=_assessment(main3, changes, 5.))
        shell3 = _difference(cases, material + '18_500', material + '10_500', 3)
        shell2 = _difference(cases, material + '18_500', material + '10_500', 2)
        activation[material] = dict(
            matched_PW_eV=500, response_G_eV=200, response_bands_18_10=[160, 156],
            kmesh=[12, 12, 12], matched_U18_minus_U10_N3_eV=shell3.tolist(),
            matched_U18_minus_U10_N2_eV=shell2.tolist(),
            q3_minus_q2_eV=(shell3 - shell2).tolist(),
            paired_reference_sensitivity_eV=_reference_difference(
                cases, material + '18_500', material + '10_500', 3).tolist(),
            interpretation='Response change on explicit4s/4p activation with matched geometry, '
                'PW500,q/G and nominal empty-band count. Includes PAW-representation and '
                'ground-density changes; not an isolated independent-shell susceptibility. '
                'Band/G/k/basis sensitivities are reported separately; no ten-valence k control is invented.')
    loading3 = _difference(cases, 'PdD18_1200', 'Pd18_1200', 3)
    loading2 = _difference(cases, 'PdD18_1200', 'Pd18_1200', 2)
    loading_changes = dict(
        q3_minus_q2=loading3 - loading2,
        paired_reference_parameters=_reference_difference(cases, 'PdD18_1200', 'Pd18_1200', 3),
        paired_band_sensitivity=loading2 - _difference(cases, 'PdD18_1200_nb100', 'Pd18_1200_nb100', 2),
        paired_G_sensitivity=loading2 - _difference(cases, 'PdD18_1200_G400', 'Pd18_1200_G400', 2),
        paired_k_sensitivity=loading2 - _difference(cases, 'PdD18_k8', 'Pd18_k8', 2),
        paired_basis_sensitivity=loading3 - _difference(cases, 'PdD18_500', 'Pd18_500', 3))
    loading_assessment = _assessment(loading3, loading_changes, 1.)
    contact_margin = max(1., loading_assessment['sampled_change_envelope_eV'][0])
    loading = dict(direction='PdD minus Pd', N3_eV=loading3.tolist(), N2_eV=loading2.tolist(),
        empirical_numerics=loading_assessment,
        contact_sign_outside_sampled_envelope=bool(abs(loading3[0]) > contact_margin),
        contact_sign=int(np.sign(loading3[0])), fitted_loading_step_eV=[550., 790.],
        comparison_scope='Computed fixed-lattice static contrast beside a phenomenological fitted loading step; '
                         'the fit is not a microscopic observable or calibration to add to U.')
    ten_controls = {}
    for material in ('Pd', 'PdD'):
        key = material + '10_500'
        main3 = _profile(cases, key, 3)
        ten_controls[material] = dict(U10_largeq_N3_eV=main3.tolist(),
            q3_minus_q2_eV=(main3 - _profile(cases, key, 2)).tolist(),
            bands156_minus96_eV=_difference(cases, key, key + '_nb96', 2).tolist(),
            G200_minus400_eV=_difference(cases, key, key + '_G400', 2).tolist(),
            qualification='Retained PW500 ground potential expanded to156 states without newSCF; '
                          'no new or assumed ten-valence k-grid convergence.')
    maximum = max(max(values['U18_N3_eV']) for values in absolute.values())
    return dict(schema='two-host-static-comparison-v1', radii_A=result['radii_A'],
        U18=absolute, semicore_activation=activation, loading18=loading,
        matched_largeq_U10=ten_controls,
        step3_trigger=dict(metric='Maximum sampled primary N3 U18 over both hosts/radii',
                          maximum_eV=maximum, above100eV=bool(maximum > 100),
                          real_time_started=False,
                          authority='Any above100eV result still needs explicit setup/cost approval before real-time.'),
        physical_scope='Neutral nonmagnetic PBE/RPA, fixed Pd/PdD lattice geometries, admitted PAW setups, '
                       'FD.1eV,eta.005eV and finite virtual/q/G bases. Not a global physical screening exclusion.')
