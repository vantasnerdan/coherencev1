"""Read the two-q cost discriminator; this is not a normalized BZ integral."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.constants import elementary_charge, epsilon_0

from g4_beam.screen.dynamic._upstream.microscopic_screening import q_residual_eV
from g4_beam.screen.dynamic.largeq_bank import binding, load_ground, require, response_points


RADII_A = np.asarray([0, .001, .002, .005, .01, .02, .05, .1, .2, .5, 1.])
COULOMB_EV_A = elementary_charge / (4 * np.pi * epsilon_0) * 1e10


def pilot_result(directory, material, pd_valence):
    bank = load_ground(directory, material, pd_valence)
    rows = []
    bindings = bank['bindings'][:]
    for point in response_points(directory, bank):
        profile = q_residual_eV(
            point['Q_Ainv'], point['inverse'], bank['source_A'], RADII_A,
            bank['volume_A3'], 1., 2., 4., COULOMB_EV_A)
        rows.append(dict(index=point['index'], settings=point['settings'],
                         elapsed_wall_s=point['metadata']['elapsed_wall_s'],
                         checks=point['checks'],
                         epsilon_macro_local_fields=point['metadata']['epsilon_macro_local_fields'],
                         epsilon_head_no_local_fields=point['metadata']['epsilon_head_no_local_fields'],
                         unit_weight_subtracted_q_integrand_eV=profile.tolist()))
        bindings.extend(point['bindings'])
    require(rows, 'No producer-published complete response point exists')
    inversion = None
    if len(rows) == 2:
        q0, q1 = [np.asarray(row['settings']['q_fractional']) for row in rows]
        require(np.allclose(q0 + q1, 0, rtol=0, atol=1e-12),
                'The discriminator does not contain actual opposite q points')
        difference = np.asarray(rows[0]['unit_weight_subtracted_q_integrand_eV']) - np.asarray(
            rows[1]['unit_weight_subtracted_q_integrand_eV'])
        error = float(np.max(np.abs(difference)))
        scale = max(1., *(abs(value) for row in rows
                          for value in row['unit_weight_subtracted_q_integrand_eV']))
        require(error < 1e-8 * scale,
                'Opposite-q site kernels disagree; do not reduce the production grid')
        inversion = dict(max_abs_difference_eV=error, relative_scale=scale,
                         admitted_for_this_test_pair=True)
    return dict(
        schema='coherence-static-largeq-discriminator-readout-v1',
        material=material, pd_valence=pd_valence,
        neutral_electrons=bank['ground']['neutral_electrons'],
        ground_state_settings=bank['config']['ground_state'],
        cell_volume_A3=bank['volume_A3'], source_A=bank['source_A'].tolist(),
        radii_A=RADII_A.tolist(), coulomb_eV_A=COULOMB_EV_A,
        intermediate_reference=dict(kappa_Ainv=2., gaussian_width_Ainv=4.),
        interpretation='Unit-weight Gaussian-reference-subtracted q integrands only. '
                       'Two opposite q points are not a normalized BZ quadrature. '
                       'No continuum addback or total U18 is claimed.',
        completed_points=rows, inversion_check=inversion,
        source_bindings=bindings,
        canonical_kernel=binding(Path(__file__).parent / '_upstream/microscopic_screening.py'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bank', required=True, type=Path)
    parser.add_argument('--material', choices=('Pd', 'PdD'), required=True)
    parser.add_argument('--pd-valence', choices=(10, 18), type=int, required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = pilot_result(args.bank, args.material, args.pd_valence)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(output=str(args.output), sha256=binding(args.output)['sha256'],
                          completed_points=len(result['completed_points']),
                          inversion_check=result['inversion_check'])))


if __name__ == '__main__':
    main()
