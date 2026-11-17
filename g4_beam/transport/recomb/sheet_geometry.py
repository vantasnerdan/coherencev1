"""Minimal ordered adsorption phases for sheet free energies, not reaction pairs.

0.5ML uses one D/1x2;0.25ML one D/2x2. All top Pd and D relax; bottom2
layers remain fixed. Subsurface D, chamber species and larger reconstructions
are omitted. Clean cells use the qualified uniform clean primitive geometry.
"""
import numpy as np
from surface_geometry import slab


def sheet_cells(coverage, clean_primitive, filled_primitive):
    if coverage not in (.5, .25):
        raise ValueError('Declared half or quarter monolayer required')
    if (list(clean_primitive.symbols) != ['Pd'] * 3
            or list(filled_primitive.symbols) != ['Pd'] * 3 + ['H']):
        raise ValueError('Qualified clean/filled primitive ordering required')
    size = (1, 2) if coverage == .5 else (2, 2)
    nsites = size[0] * size[1]
    ads = slab(hydrogens=1, layers=3, size=size)
    clean = ads[:-1]
    ideal = slab(hydrogens=0, layers=3, size=(1, 1))
    shift = clean_primitive.positions - ideal.positions
    for layer in range(3):
        ids = slice(layer * nsites, (layer + 1) * nsites)
        clean.positions[ids] += shift[layer]
    ads.positions[:-1] = clean.positions
    ads.positions[2 * nsites:3 * nsites, 2] += coverage * (
        filled_primitive.positions[2, 2] - clean_primitive.positions[2, 2])
    top_z = float(ads.positions[2 * nsites:3 * nsites, 2].mean())
    ads.positions[-1, 2] = top_z + (
        filled_primitive.positions[3, 2] - filled_primitive.positions[2, 2])
    if (len(ads.constraints[0].get_indices()) != 2 * nsites
            or len(ads) - 3 * nsites != 1 or 1 / nsites != coverage):
        raise ValueError('Sheet stoichiometry or fixed-layer subspace changed')
    return clean, ads
