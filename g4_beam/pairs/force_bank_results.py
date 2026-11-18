"""Authenticate recorded raw forces against actual compact GPW result banks."""
from pathlib import Path
import numpy as np
from ase.io import read
from ase.io.trajectory import read_atoms
from ase.io.ulm import Reader
from g4_beam.pairs.surface_io import manifest


def verify_force_record(record):
    for name in ('gpw', 'geometry'):
        bound = record[name]
        if manifest(Path(bound['path'])) != bound:
            raise ValueError(f'Force-record {name} byte binding changed')
    geometry = read(record['geometry']['path'])
    with Reader(record['gpw']['path']) as reader:
        atoms = read_atoms(reader.atoms)
        actual = reader.results.asdict()
    if (list(atoms.symbols) != list(geometry.symbols)
            or not np.array_equal(atoms.positions, geometry.positions)
            or not np.array_equal(atoms.positions, record['positions_A'])
            or not np.array_equal(atoms.cell, geometry.cell)
            or not np.array_equal(atoms.pbc, geometry.pbc)
            or not np.array_equal(atoms.get_masses(), geometry.get_masses())):
        raise ValueError('Recorded coordinate/cell/PBC/isotope differs from actual GPW')
    forces = np.asarray(record['forces'], dtype=float)
    if (forces.shape != (len(atoms), 3) or not np.all(np.isfinite(forces))
            or not np.allclose(forces, actual['forces'], rtol=0, atol=1e-12)):
        raise ValueError('Recorded raw forces differ from actual compact GPW')
    energy = float(record['energy_eV'])
    if (not np.isfinite(energy)
            or not np.isclose(energy, actual['energy'], rtol=0, atol=1e-12)):
        raise ValueError('Recorded energy differs from actual compact GPW')
    return atoms
