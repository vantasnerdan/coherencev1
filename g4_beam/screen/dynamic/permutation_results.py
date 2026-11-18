"""Admit the actual generic-q permutation pair before reducing a family."""
import argparse
import json
from pathlib import Path

import numpy as np
from g4_beam.screen.dynamic._upstream.microscopic_screening import q_residual_eV
from g4_beam.screen.dynamic._upstream.screening_quadrature import validate_site_group
from g4_beam.screen.dynamic.largeq_bank import binding, load_ground, require, response_points
from g4_beam.screen.dynamic.largeq_results import COULOMB_EV_A, RADII_A


def compare_permutation(directory):
    bank = load_ground(directory, 'PdD', 18)
    points = list(response_points(directory, bank))
    require(len(points) == 2, 'The production gate requires two published complete points')
    q0, q1 = [np.asarray(p['settings']['q_fractional']) for p in points]
    require(not np.allclose(q0, q1, rtol=0, atol=1e-12),
            'Repeated identical q is not a permutation discriminator')
    require(np.allclose(np.sort(q0), np.sort(q1), rtol=0, atol=1e-12),
            'Actual pair is not related by the proposed coordinate permutation')
    operations = validate_site_group(
        bank['cell'], bank['ground']['geometry']['positions_A'],
        bank['ground']['geometry']['symbols'], bank['source_A'],
        bank['config']['ground_state']['kmesh'])
    witnesses = [operation for operation in operations
                 if np.allclose(q0 @ operation, q1, rtol=0, atol=1e-12)]
    require(witnesses, 'The actual ionic/source site group does not relate the q pair')
    settings0, settings1 = [dict(point['settings']) for point in points]
    settings0.pop('q_fractional')
    settings1.pop('q_fractional')
    require(settings0 == settings1, 'The symmetry pair has unmatched response settings')
    profiles = [q_residual_eV(
        point['Q_Ainv'], point['inverse'], bank['source_A'], RADII_A,
        bank['volume_A3'], 1., 2., 4., COULOMB_EV_A) for point in points]
    error = float(np.max(np.abs(profiles[0] - profiles[1])))
    scale = max(1., *(float(np.max(np.abs(profile))) for profile in profiles))
    require(error < 1e-8 * scale,
            'Independent generic-q site responses disagree; do not reduce the family')
    return dict(schema='actual-generic-permutation-admission-v1',
        admitted=True, q_fractional=[q0.tolist(), q1.tolist()],
        operations=witnesses, radii_A=RADII_A.tolist(),
        unit_weight_subtracted_integrands_eV=[profile.tolist() for profile in profiles],
        max_abs_difference_eV=error, scale_eV=scale,
        actual_point_seconds=[point['metadata']['elapsed_wall_s'] for point in points],
        checks=[point['checks'] for point in points],
        interpretation='Independent RPA covariance check, not a total BZ screening potential.',
        source_bindings=bank['bindings'] + [record for point in points for record in point['bindings']],
        runtime_bindings=[binding(Path(__file__)),
            binding(Path(__file__).parent / '_upstream/microscopic_screening.py'),
            binding(Path(__file__).parent / '_upstream/screening_quadrature.py')])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bank', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = compare_permutation(args.bank)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(output=str(args.output), admitted=result['admitted'],
                          max_abs_difference_eV=result['max_abs_difference_eV'],
                          actual_point_seconds=result['actual_point_seconds'],
                          binding=binding(args.output))))


if __name__ == '__main__':
    main()
