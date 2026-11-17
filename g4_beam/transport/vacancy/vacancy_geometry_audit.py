"""Interpret saved relaxed atoms rather than equating energy arithmetic with trapping."""
import itertools
import numpy as np
from ase.io import read
from ase.geometry import find_mic


def audit(path, record):
    atoms = read(path)
    if atoms.get_chemical_symbols() != record['symbols']:
        raise ValueError('Saved relaxed structure species differ from primitive ledger')
    # ASE extxyz stores positions to eight decimal places (rounding <=5e-9 A).
    if not np.allclose(atoms.positions, record['final_positions_A'], atol=5.1e-9, rtol=0):
        raise ValueError('Saved relaxed positions differ from primitive ledger')
    if not np.allclose(atoms.cell, record['cell_A'], atol=1e-10, rtol=0):
        raise ValueError('Saved lattice differs from primitive ledger')
    vectors = np.array([shift for shift in itertools.product(range(-2, 3), repeat=3)
                        if shift != (0, 0, 0)]) @ atoms.cell.array
    pd = [i for i, atom in enumerate(atoms) if atom.symbol == 'Pd']
    d = [i for i, atom in enumerate(atoms) if atom.symbol == 'H']
    sites = []
    for i in d:
        _, distances = find_mic(atoms.positions[pd] - atoms.positions[i], atoms.cell, pbc=True)
        _, vacancy_distance = find_mic(atoms.positions[i] - np.array(record['geometry']['missing_position_A']),
                                       atoms.cell, pbc=True)
        sites.append(dict(index=i, nearest_Pd_A=float(distances.min()),
                          distance_from_removed_Pd_A=float(vacancy_distance),
                          final_position_A=atoms.positions[i].tolist()))
    pair_distances = [float(atoms.get_distance(i, j, mic=True)) for i in d for j in d if i < j]
    return dict(shortest_periodic_translation_A=float(np.linalg.norm(vectors, axis=1).min()),
                D_sites=sites, D_D_distances_A=pair_distances,
                qualification='Geometry distances are observations; local-minimum and image adequacy need separate controls')
