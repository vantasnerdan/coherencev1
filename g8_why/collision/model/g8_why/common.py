"""Fixed-profile g8 response and explicit paper readout (NumPy/SciPy only)."""
from functools import lru_cache
from pathlib import Path
import csv
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from g4_beam.forward import model as DD
from g4_beam.forward.consumer import QE, csda_nodes, forward_yields
from g4_beam.forward.profile_thermal import transport_stopping
from g4_beam.merit.api import reference, paper_rows

OUT = Path(__file__).resolve().parent
ENERGIES = np.array([.25, .5, 1., 2.5, 6.5])
MATERIALS = ('Pd', 'PdD', 'Ti', 'TiD')
GX, GW = np.polynomial.hermite_e.hermegauss(5)
GW /= GW.sum()
U = .040


@lru_cache(None)
def response(material, nominal, fractions=(.85, .15, 0.), incident_lab=None):
    """Return local lab energies and nD*path*mb weights per true charged C.

    Profiles, sheet and stopping are immutable. Only incident transport changes
    for an explicitly different incident species or neutral energy. The profile
    key always remains the actual nominal run energy, never a molecular energy.
    """
    layers, sheet, _ = reference(material).profiles(float(nominal))
    energy, weight = [], []
    species = ((1, 1.),) if incident_lab is not None else tuple(enumerate(fractions, 1))
    for multiplicity, fraction in species:
        if not fraction:
            continue
        for g, w in zip(GX, GW):
            lab = (incident_lab if incident_lab is not None else 2*nominal/multiplicity)*(1+.05/2.3548*g)
            en, _, ww, _ = csda_nodes(float(lab), layers, transport_stopping)
            share = fraction*multiplicity*w/QE
            energy.extend(en)
            weight.extend(ww*share)
            energy.append(lab)
            weight.append(sheet*1e-31*share)
    return np.asarray(energy), np.asarray(weight)


def yield_C(material, nominal, channel='p', multiplier=None, fractions=(.85, .15, 0.), incident_lab=None):
    lab, weights = response(material, float(nominal), fractions, incident_lab)
    ec = np.maximum(lab/2, .05)
    sigma = DD.sigma(ec, U, channel)
    if multiplier is not None:
        sigma = sigma*multiplier(ec, channel)
    return float(weights@sigma)


def readout(material, energies, yields, *, current_power=0, b_hi=1., b_lo=0., stitch=1.0413535612834366, rejoin=False):
    """Own-6.5 neutron reference, fixed smooth40 p/n calibration at 4.5.

    yields(e,ch) is physical reactions/true C. Cup focusing c=Icup/Itrue;
    signal=Y/c/Icup**current_power. The explicit common stitch is fixed, not
    refitted by foil or hypothesis. rejoin=True instead recalibrates the p/n
    join for each reaction law and is only a named reporting sensitivity.
    """
    def signal(e, ch):
        return yields(float(e), ch)/DD.collected_fraction(e, b_hi, b_lo)/DD.I_cup(e)**current_power
    ref = signal(6.5, 'n')
    if rejoin:
        join = signal(4.5, 'n')/signal(4.5, 'p')*stitch
    else:
        join = yield_C(material, 4.5, 'n')/yield_C(material, 4.5, 'p')*stitch
    return np.array([signal(e, 'p')*join/ref if e <= 4. else signal(e, 'n')/ref for e in energies])

def observed(material, energies):
    table = {float(r['E_cm']): (float(r['R']), float(r['sigma_lnR']))
             for r in paper_rows() if r['material'] == material}
    return np.array([table[float(e)][0] for e in energies]), np.array([table[float(e)][1] for e in energies])


def baseline_rows():
    rows = []
    for material in MATERIALS:
        raw = lambda e, ch: yield_C(material, e, ch)
        obs, err = observed(material, ENERGIES)
        r0 = readout(material, ENERGIES, raw)
        r1 = readout(material, ENERGIES, raw, current_power=1)
        for i, e in enumerate(ENERGIES):
            rows.append(dict(material=material, E_cm_keV=float(e),
                smooth40_p_reactions_per_true_C=raw(e, 'p'),
                smooth40_n_reactions_per_true_C=raw(e, 'n'),
                observed_own65_R=obs[i], sigma_lnR=err[i],
                baseline_reported_R_k0_focus=r0[i], observed_normalized_excess_k0=obs[i]/r0[i],
                baseline_reported_R_k1_focus=r1[i], observed_normalized_excess_k1=obs[i]/r1[i],
                cup_current_uA=float(DD.I_cup(e)),
                physical_factor_identifiable=False))
    return rows


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def source_hashes():
    paths = ('g4_beam/merit/profiles.json', 'g4_beam/merit/api.py',
             'g4_beam/forward/consumer.py', 'g4_beam/forward/model.py',
             'g4_beam/forward/profile_thermal.py', 'g4_beam/forward/errbars.csv',
             'g8_why/common.py')
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}


def baseline_proof():
    """One real unchanged-consumer positive check, including the native sheet."""
    worst = 0.
    for material in MATERIALS:
        for ch in ('p', 'n'):
            direct = forward_yields(ENERGIES, reference(material).profiles,
                Ue_keV=U, measurement={'channel': ch}, reaction='eq7',
                stopping_fn=transport_stopping, reference_E_cm_keV=None)['physical_yield_per_C']
            own = np.array([yield_C(material, e, ch) for e in ENERGIES])
            worst = max(worst, float(np.max(np.abs(own/direct-1))))
    assert worst < 2e-13
    rows = baseline_rows()
    assert all(r['observed_normalized_excess_k0'] == 1. for r in rows if r['E_cm_keV'] == 6.5)
    return {'shared_smooth40_consumer_max_relative': worst,
            'reference_observation_is_one': True, 'profile_source_hashes': source_hashes()}
