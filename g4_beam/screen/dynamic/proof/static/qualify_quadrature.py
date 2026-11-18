"""Qualify midpoint coverage/site covariance and actual q-specific k costs."""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path

import gpaw
from gpaw import GPAW
from gpaw.mpi import serial_comm
from gpaw.response.df import DielectricFunction
from gpaw.response.kpoints import KPointDomainGenerator, ResponseKPointGrid
from gpaw.response.symmetry import QSymmetryAnalyzer
import numpy as np
from scipy.constants import elementary_charge, epsilon_0
from scipy.spatial import cKDTree

from g4_beam.screen.dynamic._upstream import microscopic_screening, screening_quadrature
from g4_beam.screen.dynamic._upstream.microscopic_screening import q_residual_eV
from g4_beam.screen.dynamic._upstream.screening_quadrature import midpoint_classes, validate_site_group

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bank', required=True, type=Path)
parser.add_argument('--setups', required=True)
parser.add_argument('--output', required=True, type=Path)
args = parser.parse_args()
config = json.loads((args.bank / 'input.json').read_text())
ground = json.loads((args.bank / 'ground.json').read_text())
cell = np.asarray(ground['geometry']['cell_A'])
positions = np.asarray(ground['geometry']['positions_A'])
source = positions[1]
mesh = config['ground_state']['kmesh']
operations = validate_site_group(cell, positions, ground['geometry']['symbols'], source, mesh)
try:
    validate_site_group(cell, positions, ground['geometry']['symbols'], source + [.031, .047, .059], mesh)
except ValueError as error:
    wrong_site_refusal = str(error)
else:
    raise AssertionError('An inequivalent source falsely admitted the site group')

coverage = []
for n in (2, 3, 6):
    classes = midpoint_classes(n)
    coordinates = (np.arange(n) + .5) / n
    coordinates = np.where(coordinates >= .5, coordinates - 1, coordinates)
    full = np.asarray(list(product(coordinates, repeat=3)))
    integer_members = np.asarray([member for row in classes for member in row['members_numerators']])
    full_numerators = np.rint(full * (2 * n)).astype(int)
    assert np.max(np.abs(full * (2 * n) - full_numerators)) < 1e-12
    assert {tuple(q) for q in full_numerators} == {tuple(q) for q in integer_members}
    assert len(integer_members) == len({tuple(q) for q in integer_members}) == n**3
    flattened = integer_members / (2 * n)
    assert sum(row['multiplicity'] for row in classes) == n**3
    assert np.min(np.linalg.norm(flattened, axis=1)) > 0
    representatives = np.asarray([row['q_fractional'] for row in classes])
    assert np.max(np.abs(representatives * mesh - np.rint(representatives * mesh))) < 1e-12
    weights = np.asarray([row['weight'] for row in classes])
    assert abs(weights.sum() - 1) < 1e-14
    def functions(q):
        c = np.cos(2 * np.pi * q)
        s = np.sin(2 * np.pi * q)
        cross = s[:, 0] * s[:, 1] + s[:, 0] * s[:, 2] + s[:, 1] * s[:, 2]
        return np.column_stack((np.exp(c.sum(axis=1) + .1 * cross),
                                (1 + c.sum(axis=1))**2 + .3 * cross**2,
                                c.prod(axis=1)**2, cross))
    direct = functions(full).mean(axis=0)
    reduced = weights @ functions(representatives)
    error = float(np.max(np.abs(direct - reduced)))
    assert error < 2e-13
    wrong_weights = weights.copy()
    wrong_weights[0] += .001
    wrong_weights[1] -= .001
    mutation_error = float(np.max(np.abs(wrong_weights @ functions(representatives) - direct)))
    assert mutation_error > 1e-6
    coverage.append(dict(n=n, classes=classes, full_points=n**3,
                         independent_sum_max_abs_error=error, wrong_weight_error=mutation_error))

gpaw.setup_paths[:] = [args.setups]
calc = GPAW(str(args.bank / 'ground.gpw'), legacy_gpaw=False, communicator=serial_comm, txt=None)
assert calc.get_number_of_electrons() == ground['neutral_electrons'] == 19
# Constructor/adapters only: no SCF, chi0 or inverse-dielectric calculation.
df = DielectricFunction(calc, frequencies=[0.0], ecut=200, eta=.005,
                        nbands=160, hilbert=False, intraband=True,
                        nblocks=1, world=serial_comm, txt=None)
kd = df.gs.kd
for operation in operations:
    assert any(np.array_equal(operation, stored) for stored in kd.symmetry.op_scc)
analyzer = QSymmetryAnalyzer()
kpoints = ResponseKPointGrid(kd)
def groups(q):
    symmetries = analyzer.analyze_symmetries(np.asarray(q), kd)
    generator = KPointDomainGenerator(symmetries, kpoints)
    return len(generator.group_kpoints()), len(symmetries)
assert groups([.25, .25, .25]) == (364, 6)
for plan in coverage:
    for row in plan['classes']:
        row['native_k_groups'], row['native_q_symmetries'] = groups(row['q_fractional'])
    plan['cost_equivalent_q111_points'] = sum(row['native_k_groups'] / 364 for row in plan['classes'])

with np.load(args.bank / 'response-00.npz', allow_pickle=False) as data:
    Q = data['Q_Ainv']
    inverse = data['inverse_epsilon_coulomb_symmetrized']
tree = cKDTree(Q)
radii = np.asarray([0, .001, .01, .1, .5, 1.])
coulomb = elementary_charge / (4 * np.pi * epsilon_0) * 1e10
profile = q_residual_eV(Q, inverse, source, radii, ground['cell_volume_A3'], 1., 2., 4., coulomb)
covariance = []
for operation in operations:
    operation = np.asarray(operation)
    rotated = Q @ operation
    if not np.allclose(rotated[0], Q[0], atol=1e-10):
        continue
    distance, permutation = tree.query(rotated)
    assert distance.max() < 1e-9 and len(set(permutation)) == len(Q)
    transformed = inverse[np.ix_(permutation, permutation)]
    matrix_error = float(np.max(np.abs(transformed - inverse)))
    assert matrix_error < 1e-8
    rotated_profile = q_residual_eV(rotated, transformed, source, radii,
                                    ground['cell_volume_A3'], 1., 2., 4., coulomb)
    profile_error = float(np.max(np.abs(rotated_profile - profile)))
    assert profile_error < 1e-8 * max(1., np.max(np.abs(profile)))
    covariance.append(dict(operation=operation.tolist(), matrix_max_abs_error=matrix_error,
                           profile_max_abs_error_eV=profile_error))
assert len(covariance) == 6

def bind(path):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest)
result = dict(schema='qualified-midpoint-site-cost-v1', coverage=coverage,
              actual_site_operations=operations, wrong_site_refusal=wrong_site_refusal,
              actual_matrix_covariance=covariance,
              sources=[bind(path) for path in [__file__, microscopic_screening.__file__,
                       screening_quadrature.__file__, args.bank / 'input.json',
                       args.bank / 'ground.json', args.bank / 'ground.gpw',
                       args.bank / 'response-00.npz']], no_scf_or_response_calculated=True)
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(dict(output=str(args.output), grids=[dict(n=p['n'], classes=len(p['classes']),
                           equivalent_q111=p['cost_equivalent_q111_points'],
                           independent_error=p['independent_sum_max_abs_error']) for p in coverage],
                      covariance_max=max(row['matrix_max_abs_error'] for row in covariance),
                      source_phase_max_eV=max(row['profile_max_abs_error_eV'] for row in covariance),
                      wrong_site_refusal=wrong_site_refusal)))
