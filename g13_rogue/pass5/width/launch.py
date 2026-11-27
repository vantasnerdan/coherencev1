"""Launch-only heat-width and electron handover sensitivity of accepted Pass 5."""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import source as s


t = s.t
# First row is the accepted baseline, also used in the handover comparison.
CASES = [(w, 0.) for w in (2., .5, 1., 5., 10.)] + [(2., tau) for tau in (.1, 1., 10.)]
WIDTH, TAU = np.array(CASES).T * np.array([[1e-9], [1e-12]])


def thermal_variants(total, electronic, k, omega):
    """Keep recoil heat prompt; electron response is causal with mean delay TAU."""
    transfer_minus_one = (1j * TAU * omega[..., None]) / (1 - 1j * TAU * omega[..., None])
    form = np.exp(-.5 * (k[..., None] * WIDTH)**2)
    return (total[..., None] + electronic[..., None] * transfer_minus_one) * form


def terms(mechanical, variance, thermal):
    mechanical = mechanical[..., None]
    power = np.broadcast_to(abs(mechanical)**2 + variance[..., None], thermal.shape)
    return np.stack((power, abs(thermal)**2, 2 * np.real(mechanical * thermal.conj())), axis=-1)


def guided(data, values):
    """Pass 5's free-face thickness projection, with electron and other heat separated."""
    ks = np.geomspace(.1, 1.2e5, 120)
    modes = np.arange(64) * np.pi / t.H
    mode_weights = np.r_[1., np.full(63, 2.)] / t.H
    basis = [np.array([t.legval(2 * tr['z'] / t.H - 1, np.eye(15)[j]) *
                       np.sqrt((2*j+1)/(t.RHO*t.H)) for j in range(15)]) for tr in data]
    mode_at_row = [np.cos(np.outer(modes, tr['z'])) for tr in data]
    for ia, grain in enumerate((9, 12, 18, 23, 40, 45, 51, 54)):
        theta = 2 * np.pi * (ia + .5) / 8
        for k, dk in zip(ks, np.gradient(ks)):
            omega, vec, _, _, z, ww, strain = t.plate_matrices(k, theta, grain)
            freq = omega / (2 * np.pi)
            band = np.searchsorted(t.BANDS, freq, side='right') - 1
            use = np.flatnonzero((band >= 0) & (band < 2) & (freq < 1e7))
            if not len(use):
                continue
            om, vec = omega[use], vec[:, use]
            div = (strain[:, :3].sum(axis=1) @ vec).conj()
            proj = ((div.T * ww) @ np.cos(np.outer((z + 1)*t.H/2, modes))) * mode_weights
            proj *= (t.BETA/(t.RHO*t.CP)) / (t.DTH*(k*k + modes*modes)[None, :] - 1j*om[:, None])
            for species, tr in enumerate(data):
                phase = np.exp(1j * tr['t'][:, None] * om)
                up = basis[species].T @ vec[30:]
                ap = t.PARENT_P[species] * np.sum(tr['momentum'][:, None] * up.conj() * phase, axis=0)
                ap *= np.exp(-.5 * (k * s.BATH)**2)
                total = mode_at_row[species] @ (tr['deposit'].sum(axis=1)[:, None] * phase)
                total += (np.cos(modes*tr['terminal_z'])[:, None] * tr['energy'][3] *
                          np.exp(1j * om * tr['terminal_t']))
                electron = mode_at_row[species] @ (tr['deposit'][:, 0, None] * phase)
                at, ae = (np.sum(proj * h.T, axis=1) for h in (total, electron))
                thermal = thermal_variants(at, ae, np.full_like(om, k), om)
                energy = terms(ap, np.zeros_like(om), thermal) * (.5*k*dk/(2*np.pi*8))
                for pos, j in enumerate(use):
                    values[species, band[j]] += energy[pos]


def bulk(data, values):
    n, dw = s.angular_grid()
    intervals = [(max(t.BANDS[b], 1e7), t.BANDS[b+1], b) for b in range(1, 7)]
    intervals += [(2.5e12, 3.399e12, 7), (3.399e12, 4.699e12, 8), (4.699e12, 6.799e12, 9)]
    _, polarization = t.christoffel(n)
    longitudinal = np.einsum('ni,nib->nb', n, polarization)
    nodes, weights = leggauss(8)
    bath, caps = s.BATH, t.FCAP
    for species, tr in enumerate(data):
        mark = s.marks(tr)
        for lo, hi, band in intervals:
            span = np.log(hi/lo)
            freqs = np.exp(np.log(lo) + (nodes+1)*span/2)
            for f, df in zip(freqs, weights*span/2*freqs):
                # Get the common unshaped phasors only once. Reapply the fixed
                # mechanical form separately; it is not part of this heat scan.
                s.BATH, t.FCAP = 0., np.maximum(caps, f*1.0001)
                try:
                    ap, vv, at, k, vg = s.amplitudes(f, n, tr, mark)
                finally:
                    s.BATH, t.FCAP = bath, caps
                omega = 2*np.pi*f
                ep = np.exp(1j*omega*tr['t'] - 1j*k[:, :, None]*n[:, 2, None, None]*tr['z'])
                electron_heat = ep @ tr['deposit'][:, 0]
                ae = (-1j*t.BETA/(t.RHO*t.CP) * k * longitudinal * electron_heat /
                      (t.DTH*k*k - 1j*omega))
                form = np.exp(-.5 * (k*bath)**2)
                thermal = thermal_variants(at, ae, k, np.asarray(omega))
                dos = k*k/(2*t.RHO*(2*np.pi)**3*vg)*dw[:, None]*(2*np.pi*df)
                dos[:, f >= caps] = 0.
                values[species, band] += np.einsum('nb,nbic->ic', dos, terms(ap*form, vv*form**2, thermal))


def launch(reference):
    with Path(reference).open('rb') as fh:
        if t.hashlib.file_digest(fh, 'sha256').hexdigest() != t.EXPECTED['signed_r2.npz']:
            raise ValueError('Expected the same signed Rung 2 input as Pass 5')
    with np.load(reference, allow_pickle=False) as archive:
        data = s.tracks({'heatenergy': archive['heatenergy']})
    values = np.zeros((2, 10, len(CASES), 3))
    guided(data, values)
    bulk(data, values)
    for species, tr in enumerate(data):
        q = float(tr['energy'].sum())
        for case in range(len(CASES)):
            mech, therm, cross = values[species, :, case].sum(axis=0)
            h = brentq(lambda h: h + (mech+h*h*therm+h*cross)/q - 1, 0., 1., xtol=1e-15)
            values[species, :, case, 1] *= h*h
            values[species, :, case, 2] *= h
    with np.load(HERE.parent/'source.npz', allow_pickle=False) as accepted:
        # Compare physical per-species and per-band terms, including source
        # above 2.5 THz that contributes to the heat budget but not the table.
        ref = accepted['terms'][0].sum(axis=2)
        np.testing.assert_allclose(values[:, :, 0], ref, rtol=2e-10, atol=1e-42)
    power = t.RATE * np.einsum('s,sbic->ib', t.FRACTION, values)
    return power


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--signed-reference', type=Path,
                        default=s.ROOT/'g13_rogue/pass3/rung2/signed_r2.npz')
    power = launch(parser.parse_args().signed_reference)
    fields = ['width_nm', 'electron_handover_ps', 'seven_band_total_W', 'above_ladder_W', 'all_modeled_total_W']
    fields += [f'{lo:g}-{hi:g}_Hz_W' for lo, hi in zip(t.BANDS[:-1], t.BANDS[1:])]
    with (HERE/'launch.csv').open('w', newline='') as fh:
        csv_writer = csv.writer(fh)
        csv_writer.writerow(fields)
        for (width, tau), row in zip(CASES, power):
            csv_writer.writerow((width, tau, row[:7].sum(), row[7:].sum(), row.sum(), *row[:7]))
    print('Accepted baseline matches; eight launch-only cases written to', HERE/'launch.csv')


if __name__ == '__main__':
    main()
