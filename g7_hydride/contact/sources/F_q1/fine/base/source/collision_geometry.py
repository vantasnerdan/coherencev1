"""Matched frozen PdD/target-vacancy geometry from the actual relaxed bulk record."""
from pathlib import Path
import json
import numpy as np
from ase import Atom
from ase.build import bulk


def collision_atoms(spec):
    original = json.loads(Path(spec['relaxed_ground_record']).read_bytes())
    cell = np.asarray(original['geometry']['cell_A'])
    a = float(np.linalg.norm(cell[0]) * np.sqrt(2.0))
    atoms = bulk('PdH', 'rocksalt', a=a, cubic=True).repeat(tuple(spec['repeat']))
    target = a * np.array([0.5, 0.5, 0.5])
    # ASE's conventional rocksalt cell may choose an equivalent origin.
    distances = np.linalg.norm(atoms.positions - target, axis=1)
    target_index = int(np.argmin(distances))
    if atoms[target_index].symbol != 'H' or distances[target_index] > 1e-8:
        raise ValueError('The relaxed conventional body octahedron was not reconstructed')
    displacement = np.asarray(spec.get('D_displacement_A', [0.0, 0.0, 0.0]))
    for atom in atoms:
        if atom.symbol == 'H':
            atom.position += displacement
            atom.mass = 2.01410177812
    target = target + displacement
    if not spec['occupied_target']:
        del atoms[target_index]
    direction = np.asarray(spec['direction'], dtype=float)
    direction /= np.linalg.norm(direction)
    transverse = np.cross(direction, [0.0, 0.0, 1.0])
    if np.linalg.norm(transverse) < 1e-8:
        transverse = np.cross(direction, [0.0, 1.0, 0.0])
    transverse /= np.linalg.norm(transverse)
    point = target - spec['start_R_A'] * direction + spec['impact_A'] * transverse
    atoms += Atom('H', position=point, mass=2.01410177812)
    atoms.set_pbc(True)
    return atoms, dict(target_A=target.tolist(), direction=direction.tolist(),
                       transverse=transverse.tolist(), projectile_index=len(atoms)-1,
                       lattice_A=a, target_present=spec['occupied_target'])


def projectile_point(path, longitudinal_R_A, impact_A):
    return (np.asarray(path['target_A']) - longitudinal_R_A * np.asarray(path['direction'])
            + impact_A * np.asarray(path['transverse']))
