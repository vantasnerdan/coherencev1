"""Material-order matched 0/one-A/one-B/two-ion geometries for bulk PdD."""
import numpy as np

from readout_bulk import inversion

AMPS = (12, 24, 36, 48, 80)
PAW_PARTIAL_WAVE = {
    'Pd': {'max_rc_bohr': 2.5, 'setup': 'setups/Pd.10.PBE.gz'},
    'H': {'max_rc_bohr': 0.9, 'setup': 'setups/H.PBE.gz'},
}
PD_H_RC_SUM_A = sum(v['max_rc_bohr'] for v in PAW_PARTIAL_WAVE.values()) * 0.529177210903
H_AE_SHA = 'd6daba8aa25bbeeb716c7c9a3a3c7ec67fd1b56673550b523f308400ad8064cc'
H_BASIS_SHA = '4d5e69dc06a8cd80a8875c874b44dd79e96d8e8beb6b8501e8f922c0c4184db3'
SENSITIVITY_GEOMETRY = {'zero': (0., 'AB'), 'one': (.8, 'A'),
                        'closest': (.8, 'AB'), 'r1928': (.48, 'AB')}


def paw_support(profile):
    shas = {r['path']: r['sha256'] for r in profile['sources']}
    return dict(setups={element: dict(spec, sha256=shas[spec['setup']]) for element, spec in PAW_PARTIAL_WAVE.items()},
                tabulated_rc_sum_A=PD_H_RC_SUM_A, bound='Nominal valence partial-wave radii; not an error bound.')


def paw_overlap(profile, positions):
    """Minimum Pd–D distance over all pairs in the adopted orthogonal cell."""
    g = profile['geometry']
    numbers, cell, r = np.asarray(g['atomic_numbers']), np.asarray(g['cell_A']), np.asarray(positions)
    gram = cell @ cell.T
    if not np.allclose(gram, np.diag(np.diag(gram)), atol=1e-12, rtol=0):
        raise ValueError('Pd–D overlap diagnostic requires the adopted orthogonal cell')
    delta = r[numbers == 46, None, :] - r[numbers == 1][None, :, :]
    frac = delta @ np.linalg.inv(cell)
    closest = float(np.min(np.linalg.norm((frac-np.rint(frac)) @ cell, axis=-1)))
    return dict(minimum_periodic_Pd_D_A=closest, widest_tabulated_Pd_H_partial_wave_overlap=closest < PD_H_RC_SUM_A)


def geometry(profile):
    g = profile['geometry']
    r0, a = np.asarray(g['reference_positions_A']), np.asarray(g['cell_A'])
    d = np.asarray(g['H_indices'])
    i = g['moving_index']
    dx = r0[d] - r0[i]
    frac = dx @ np.linalg.inv(a)
    shifts = -np.rint(frac).astype(int)
    mic = dx + shifts @ a
    expect = np.array([1., 1., 0.]) * profile['a_A'] / 2
    matches = np.flatnonzero(np.linalg.norm(mic - expect, axis=1) < 1e-9)
    if len(matches) != 1 or i != 25 or int(d[matches[0]]) != 63:
        raise ValueError('Nearest-neighbour tagged pair differs from registered original indices')
    k = matches[0]
    j, shift = int(d[k]), shifts[k]
    distance = np.linalg.norm(mic[k])
    nonzero = np.linalg.norm(mic, axis=1)
    if not np.isclose(distance, min(nonzero[nonzero > 1e-9]), atol=1e-9, rtol=0):
        raise ValueError('Selected image is not a D nearest neighbour')
    perm = inversion(g, r0[i]+mic[k]/2)
    if perm[i] != j or perm[j] != i:
        raise ValueError('Midpoint inversion does not exchange this pair')
    return i, j, shift, mic[k]/distance, distance, perm


def displace(profile, u, role):
    i, j, shift, n, rzero, perm = geometry(profile)
    g = profile['geometry']
    pos = np.asarray(g['reference_positions_A']).copy()
    if role in ('A', 'AB'):
        pos[i] += u*n
    if role in ('B', 'AB'):
        pos[j] -= u*n
    cell = np.asarray(g['cell_A'])
    distance = float(np.linalg.norm(pos[j] + shift @ cell - pos[i]))
    if role == 'AB' and not np.isclose(distance, rzero-2*u, atol=1e-12, rtol=0):
        raise ValueError('Pair displacement is not equal/opposite along the original image')
    return pos, i, j, shift, n, distance, perm
