"""Read actual pass-4-format saved response banks without running SCF/RPA."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def binding(path):
    path = Path(path)
    digest = hashlib.sha256()
    size = 0
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
            size += len(chunk)
    return dict(path=path.name, bytes=size, sha256=digest.hexdigest())


def bound_file(directory, record):
    # The saved path is provenance, not authority to open another directory.
    path = directory / Path(record['path']).name
    actual = binding(path)
    require(actual['bytes'] == record['bytes'] and
            actual['sha256'] == record['sha256'],
            'Saved bank bytes differ: ' + path.name)
    return path, actual


def load_ground(directory, material, pd_valence):
    directory = Path(directory)
    config = json.loads((directory / 'input.json').read_text())
    ground = json.loads((directory / 'ground.json').read_text())
    expected_symbols = {'Pd': ['Pd'], 'PdD': ['Pd', 'H']}[material]
    require(pd_valence in (10, 18), 'Only registered Pd10/Pd18 hosts are supported')
    for geometry in (config['geometry'], ground['geometry']):
        require(geometry['symbols'] == expected_symbols,
                'Physical target composition differs from ' + material)
    expected_ne = pd_valence + (material == 'PdD')
    require(ground['neutral_electrons'] == config['valence_electrons'] == expected_ne,
            'Actual neutral-electron count differs from the selected PAW host')
    require(ground['case'] == config['case'] and ground['model'] == config['model'],
            'Ground and registered case/model differ')
    require(ground['ground_state_settings'] == config['ground_state'],
            'Ground settings differ from the registered input')
    require(not config['ground_state']['spinpol'] and
            abs(ground['magnetic_moment_muB']) < 1e-10,
            'The admitted static reduction requires the nonmagnetic host')
    seed, actual = config['geometry'], ground['geometry']
    cell = np.asarray(actual['cell_A'], float)
    require(np.allclose(cell, seed['cell_A'], rtol=0, atol=1e-9),
            'Fixed registered cell changed')
    displacement = (np.asarray(actual['positions_A']) -
                    np.asarray(seed['positions_A'])) @ np.linalg.inv(cell)
    displacement -= np.rint(displacement)
    require(np.max(np.abs(displacement @ cell)) < 1e-8,
            'Fixed registered atomic geometry changed')
    volume = float(abs(np.linalg.det(cell)))
    require(np.isclose(volume, ground['cell_volume_A3'], rtol=1e-10, atol=1e-10),
            'Saved physical volume differs from the actual cell')
    _, gpw_binding = bound_file(directory, ground['gpw'])
    records = [binding(directory / name) for name in ('input.json', 'ground.json')]
    records.append(gpw_binding)
    source = (np.asarray(actual['positions_A'][1], float) if material == 'PdD'
              else np.asarray([-.5, .5, .5]) @ cell)
    return dict(config=config, ground=ground, cell=cell, volume_A3=volume,
                source_A=source, bindings=records, material=material,
                pd_valence=pd_valence)


def response_points(directory, bank):
    directory = Path(directory)
    reciprocal = 2 * np.pi * np.linalg.inv(bank['cell']).T
    for index, settings in enumerate(bank['config']['responses']):
        record_path = directory / f'response-{index:02d}.json'
        if not record_path.exists():
            continue  # Only producer-published complete points are recoverable.
        record = json.loads(record_path.read_text())
        require(record['settings'] == settings, 'Response settings differ from input')
        path, arrays_binding = bound_file(directory, record['arrays'])
        with np.load(path, allow_pickle=False) as arrays:
            Q = arrays['Q_Ainv']
            epsilon = arrays['epsilon_coulomb_symmetrized']
            inverse = arrays['inverse_epsilon_coulomb_symmetrized']
        q = np.asarray(settings['q_fractional'], float)
        n = len(Q)
        require(Q.shape == (n, 3) and epsilon.shape == inverse.shape == (n, n),
                'Saved physical vector/matrix dimensions differ')
        require(n == record['matrix_dimension'] and np.isfinite(Q).all() and
                np.isfinite(epsilon).all() and np.isfinite(inverse).all(),
                'Nonfinite or incorrectly recorded response bank')
        require(np.allclose(Q[0], q @ reciprocal, rtol=1e-10, atol=1e-10),
                'Saved head has incorrect physical q or length units')
        fractional_G = Q @ np.linalg.inv(reciprocal) - q
        require(np.max(np.abs(fractional_G - np.rint(fractional_G))) < 1e-8,
                'Saved vectors are not q plus this cell reciprocal lattice')
        require(np.min(np.linalg.norm(Q, axis=1)) > 0,
                'This finite-q bank contains a singular zero vector')
        identity_error = float(np.linalg.norm(epsilon @ inverse - np.eye(n), ord=np.inf))
        hermitian_error = float(np.linalg.norm(inverse - inverse.conj().T, ord=np.inf))
        require(identity_error < 1e-8, 'Saved dielectric inverse identity failed')
        require(hermitian_error < 1e-8, 'Static inverse is not Hermitian at the admitted tolerance')
        eigenvalues = np.linalg.eigvalsh(inverse)
        require(eigenvalues[0] > -1e-8 and eigenvalues[-1] < 1 + 1e-8,
                'Static RPA inverse violates the admitted passive range')
        checks = dict(inverse_identity_infinity_residual=identity_error,
                      inverse_hermitian_infinity_residual=hermitian_error,
                      inverse_eigenvalue_range=[float(eigenvalues[0]), float(eigenvalues[-1])],
                      matrix_dimension=n, max_Q_Ainv=float(np.linalg.norm(Q, axis=1).max()))
        yield dict(index=index, settings=settings, metadata=record, Q_Ainv=Q,
                   inverse=inverse, checks=checks,
                   bindings=[binding(record_path), arrays_binding])
