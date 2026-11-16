"""Read actual bound native diagnostics; refuse wrong cases or physical eligibility."""
import argparse
import json
from pathlib import Path

import numpy as np

from screen_io import manifest, require, write_json


def bound_file(root, reference):
    relative = Path(reference['path'])
    require(not relative.is_absolute() and '..' not in relative.parts,
            'Artifact path escapes the saved result')
    path = root/relative
    require(manifest(path, reference['path']) == reference,
            f'Artifact byte binding failed: {reference["path"]}')
    return path


def read_fingerprint(root, expected_case, require_eligible=False):
    root = Path(root)
    completed = json.loads((root/'completed.json').read_bytes())
    require(completed['eligible'] is False,
            'A fingerprint must not claim physical eligibility')
    producer = json.loads(bound_file(root, completed['producer']).read_bytes())
    config = producer['config']
    require(config['case'] == expected_case, 'Saved case differs from requested case')
    bound_file(root, producer['input'])
    require(json.loads((root/'input.json').read_bytes()) == config,
            'Producer configuration differs from saved input')
    for reference in producer['producer_source']:
        bound_file(root, reference)
    ground = json.loads(bound_file(root, completed['ground']).read_bytes())
    response = json.loads(bound_file(root, completed['response']).read_bytes())
    require(ground['case'] == response['case'] == expected_case,
            'Mixed ground/response case identities')
    require(ground['gpw'] == response['gpw'], 'Response used a different saved state')
    require(response['eligible'] is False and ground['dos']['eligible'] is False,
            'Unconverged source incorrectly claims eligibility')
    require(not require_eligible,
            'Fingerprint has no k/eta/angular/tail convergence: consumer use refused')
    reciprocal = 2*np.pi*np.linalg.inv(np.asarray(ground['geometry']['cell_A'])).T
    require(len(response['points']) == len(config['responses']),
            'Incomplete response settings family')
    points = []
    for row, settings in zip(response['points'], config['responses']):
        require(row['settings'] == settings, 'Response settings mismatch')
        with np.load(bound_file(root, row['arrays']), allow_pickle=False) as data:
            q = data['Q_Ainv']
            epsilon = data['epsilon_coulomb_symmetrized']
            inverse = data['inverse_epsilon_coulomb_symmetrized']
            require(epsilon.shape == inverse.shape == (len(q), len(q)),
                    'Saved dielectric matrix shape mismatch')
            require(np.isfinite(epsilon).all() and np.isfinite(inverse).all(),
                    'Saved dielectric matrix is nonfinite')
            require(np.allclose(q[0], np.asarray(settings['q_fractional'])@reciprocal,
                                rtol=1e-10, atol=1e-10),
                    'Saved response head is not the requested physical q')
            residual = float(np.linalg.norm(epsilon@inverse-np.eye(len(q)), ord=np.inf))
            require(residual < 1e-8, 'Saved dielectric-inverse identity failed')
            require(np.allclose(1/np.diag(epsilon),
                                data['inverse_diagonal_no_local_fields'],
                                rtol=1e-10, atol=1e-12), 'No-local-field convention mismatch')
            require(np.allclose(np.diag(inverse),
                                data['inverse_diagonal_local_fields'],
                                rtol=1e-10, atol=1e-12), 'Local-field convention mismatch')
            head_nlf = epsilon[0, 0]
            head_lf = 1/inverse[0, 0]
            require(np.allclose([head_nlf.real, head_nlf.imag],
                                row['epsilon_head_no_local_fields'],
                                rtol=1e-10, atol=1e-12), 'No-local-field head metadata mismatch')
            require(np.allclose([head_lf.real, head_lf.imag],
                                row['epsilon_macro_local_fields'],
                                rtol=1e-10, atol=1e-12), 'Local-field head metadata mismatch')
            points.append(dict(settings=settings, Q_head_Ainv=float(np.linalg.norm(q[0])),
                               matrix_dimension=len(q), inverse_residual=residual,
                               epsilon_head_NLF=[float(head_nlf.real), float(head_nlf.imag)],
                               epsilon_macro_LF=[float(head_lf.real), float(head_lf.imag)]))
    return dict(schema='coherence-screen-fingerprint-readout-v1', case=expected_case,
                completed=manifest(root/'completed.json'), dos=ground['dos'],
                points=points, eligible=False,
                reason='Diagnostics only; this is not a material screening result')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-case', required=True)
    parser.add_argument('--require-eligible', action='store_true')
    args = parser.parse_args()
    result = read_fingerprint(args.input, args.expected_case, args.require_eligible)
    write_json(args.output, result)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
