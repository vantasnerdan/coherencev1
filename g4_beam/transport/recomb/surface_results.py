"""Constrained endpoint convergence from actual saved GPW energy/force banks."""
import json
from pathlib import Path
import numpy as np
from ase.io.trajectory import read_atoms
from ase.io.ulm import Reader
from surface_io import manifest


def verify_endpoint(point, checkpoint, fmax_target):
    if manifest(Path(checkpoint['path'])) != checkpoint:
        raise ValueError('Endpoint GPW byte binding changed')
    with Reader(checkpoint['path']) as reader:
        atoms = read_atoms(reader.atoms)
        actual = reader.results.asdict()
    if (list(atoms.symbols) != point['geometry']['symbols']
            or not np.array_equal(atoms.positions, point['geometry']['positions_A'])
            or not np.array_equal(atoms.cell, point['geometry']['cell_A'])):
        raise ValueError('Endpoint GPW coordinate binding changed')
    if not np.isclose(actual['energy'], point['energy_extrapolated_eV'],
                      rtol=0, atol=1e-12):
        raise ValueError('Endpoint energy differs from actual GPW')
    raw = np.asarray(point['raw_forces_eV_A'])
    if (not np.all(np.isfinite(raw))
            or not np.allclose(raw, actual['forces'], rtol=0, atol=1e-12)):
        raise ValueError('Endpoint raw forces differ from actual GPW')
    constrained = raw.copy()
    for constraint in atoms.constraints:
        constraint.adjust_forces(atoms, constrained)
    if not np.allclose(constrained, point['forces_eV_A'], rtol=0, atol=1e-12):
        raise ValueError('Endpoint declared constrained forces do not match GPW')
    fmax = float(np.linalg.norm(constrained, axis=1).max())
    if not np.isclose(fmax, point['fmax_eV_A'], rtol=0, atol=1e-12) or fmax > fmax_target:
        raise ValueError('Endpoint has not reached its actual mobile-force target')
    return atoms


def read_converged(directory):
    directory = Path(directory)
    result = json.loads((directory / 'result.json').read_bytes())
    if not result['converged']:
        raise ValueError('No converged endpoint result')
    saved = result['final']
    if manifest(Path(saved['path'])) != saved:
        raise ValueError('Endpoint force-point bytes changed')
    point = json.loads(Path(saved['path']).read_bytes())
    atoms = verify_endpoint(point, result['checkpoint'], result['fmax_target_eV_A'])
    return result, point, atoms
