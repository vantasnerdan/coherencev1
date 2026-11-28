"""Replay the retained G16 size and band checks from the diagnostic capsule."""
import csv
import json
from pathlib import Path
import sys

import h5py
import numpy as np

from native import require
from phonon_forces import assemble
from phonon_geometry import phonon
from phonon_size_change import compare, norm, response


CASES = (
    ('dispersion', 'primitive2', 300, '010', False, False),
    ('dispersion', 'primitive2', 3000, '010', False, False),
    ('dispersion', 'conventional2', 300, '014', False, False),
    ('dispersion', 'conventional2', 3000, '014', False, False),
    ('dispersion-band-control', 'primitive2', 3000, '016', True, False),
    ('dispersion-tight-reference', 'primitive2', 3000, '017', False, True),
    ('dispersion-tight-band-control', 'primitive2', 3000, '017', True, True),
)


def saved_matrix(path):
    with h5py.File(path) as f:
        return f['force_constants'][...]


def replay_native(root, folder, cell, t, source, compact, tight):
    lane = f'primitive2{folder.removeprefix("dispersion")}' if folder != 'dispersion' else cell
    path = root / 'native-phonons' / lane / 'k16' / str(t)
    ph, points = assemble(path, cell, 16, t, root / 'sources' / f'full-zone-producer-{source}',
                          compact_bands=compact, tight_forces=tight)
    products = root / folder / cell / 'k16' / str(t)
    record = json.loads((products / 'dispersion.json').read_bytes())
    require(points['sources'] == record['sources'], f'Mismatched native output in {products}')
    np.testing.assert_allclose(ph.force_constants, saved_matrix(products / 'force_constants_raw.hdf5'),
                               rtol=0, atol=1e-10)
    ph.symmetrize_force_constants(level=3, show_drift=False)
    np.testing.assert_allclose(ph.force_constants, saved_matrix(products / 'force_constants.hdf5'), rtol=0, atol=1e-10)
    return response(root / folder, cell, 16, t)


def gamma(root, record):
    path = root / record['original_sources'][0]['file']
    with path.open() as f:
        rows = {int(r['temperature_K']): r for r in csv.DictReader(f) if r['kmesh_axis_count'] == '16'}
    anchor_eV = float(rows[3000]['electronic_excitation_eV_per_D_site'])
    require(anchor_eV == record['original_primitive_excitation_rows']['3000']['electronic_excitation_eV_per_D_site'],
            'Different original electronic excitation anchor')
    anchors = {'original_primitive': [float(rows[t]['K_Gamma_100_eV_A2']) for t in (300, 3000)]}
    for cell in ('primitive2', 'conventional2'):
        ph = phonon(cell)
        d = list(ph.primitive.symbols).index('H')
        mover = ph.primitive.p2s_map[d]
        mask = np.asarray(ph.supercell.symbols) == 'H'
        anchors[cell] = []
        for t in (300, 3000):
            row = record['supercell_rows'][cell][str(t)]
            matrix = saved_matrix(root / row['inputs'][0]['file'])[mover]
            g = matrix[mask].sum(axis=0)
            onsite = matrix[mover].trace() / 3
            np.testing.assert_allclose([g[0, 0], onsite, g.trace() / 3 - onsite],
                                       [row[k] for k in ('K_Gamma_eV_A2', 'K_onsite_eV_A2', 'K_sum_other_D_eV_A2')],
                                       rtol=0, atol=1e-10)
            anchors[cell].append(g[0, 0])
    for label, (cold, hot) in anchors.items():
        row = record['contrast'][label]
        np.testing.assert_allclose([hot - cold, 100 * (hot - cold) / cold],
                                   [row['delta_K_Gamma_eV_A2'], row['percent_change']], rtol=0, atol=1e-9)


def main(root):
    sc, sh, lc, lh, old_compact, full, compact = (replay_native(root, *row) for row in CASES)
    summary, arrays = compare(sc, sh, lc, lh)
    with (root / 'records' / 'size_sensitivity.json').open() as f:
        size = json.load(f)
    require(not summary['size_screen_passed'] and summary['outside_tolerance_q_count'] == size['outside_tolerance_q_count'],
            'Full-zone size conclusion differs')
    with np.load(root / 'records' / 'q_sensitivity.npz', allow_pickle=False) as expected:
        for key, array in arrays.items():
            np.testing.assert_allclose(array, expected[key], rtol=0, atol=1e-10)
    band = json.loads((root / 'records' / 'matched-tight-band-control.json').read_bytes())
    comparisons = ((sh, old_compact), (sh, full), (old_compact, compact), (full, compact))
    for key, (a, b) in zip(band['summary'], comparisons, strict=True):
        require(a['basis'] == b['basis'] and np.array_equal(a['q'], b['q']), 'Band-control basis differs')
        delta = norm(b['k'] - a['k'])
        np.testing.assert_allclose([delta.max(), np.sqrt(np.mean(delta**2))],
                                   [band['summary'][key][k] for k in ('max_eV_A2', 'rms_eV_A2')], rtol=0, atol=1e-10)
    gamma(root, json.loads((root / 'records' / 'collective-gamma-size-anchor.json').read_bytes()))
    print('Seven native IFC replays; collective Γ, failed size screen and both band/force controls verified.')


if __name__ == '__main__':
    main(Path(sys.argv[1]))
