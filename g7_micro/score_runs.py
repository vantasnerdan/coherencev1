"""Add a screening submodel to the existing calibrated per-run shape scorer.

No data ingest, masking, raw-count simulation, or fifth-cause registration.
CSV/NPZ are the existing g6_plateau contract plus explicit beam/profile fields.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from g6_plateau.score_runs import score
from g7_micro.model import raw_yield, raw_depth_yield


def predict_run(row, profile, law):
    """Return raw reactions/true C with the actual supplied frozen run profile.

    profile: {nominal_DD_cm_key_keV,layers,nads_m2,state,source_id}.
    law: {family,a,b,sheet}; a=f and b=U_high or exponential scale, in keV.
    No profile interpolation, nominal-foil replacement, implicit species or spread.
    """
    if row['yield_convention'] != 'reactions_per_true_C':
        raise ValueError('Join must use calibrated reactions per true charged C, not detector counts')
    lab = float(row['lab_keV'])
    if lab <= 0 or not np.isfinite(lab) or float(profile['nominal_DD_cm_key_keV']) != lab/2:
        raise ValueError('Exact source key must be nominal lab/2, never DT CM')
    if not profile['source_id'] or not profile['state']:
        raise ValueError('Actual profile source and state/history are required')
    isotope = row['isotope']
    if isotope not in ('DD', 'DT'):
        raise ValueError('Explicit DD or DT isotope pair required')
    channel = row['channel']
    if channel not in (('p', 'n', 'total') if isotope == 'DD' else ('fusion',)):
        raise ValueError('Explicit calibrated DD channel or DT fusion convention required')
    measurement = {key: float(row[key]) for key in ('D1_fraction', 'D2_fraction', 'spread_FWHM')}
    frozen = (profile['layers'], profile['nads_m2'], profile['state'])
    if law['family'] == 'depth':
        return raw_depth_yield(lab/2, frozen, f_s=float(law['f_s']), f_b=float(law['f_b']),
                               d_nm=float(law['d_nm']), U_high_keV=float(law['U_high_keV']),
                               isotope=isotope, channel=None if channel in ('total', 'fusion') else channel,
                               sheet=law['sheet'], measurement=measurement)
    return raw_yield(lab/2, {lab/2: frozen}, law['family'], float(law['a']), float(law['b']),
                     isotope=isotope, channel=None if channel in ('total', 'fusion') else channel,
                     sheet=law['sheet'], measurement=measurement)


def compare_runs(rows, covariance, profiles, laws):
    """Identical rows/covariance/nuisances for old hypotheses and each screening law."""
    result = score(rows, covariance)
    for name, material_laws in laws.items():
        replaced = []
        for row in rows:
            law = material_laws[row['material']]
            predicted = predict_run(row, profiles[row['profile_id']], law)
            replaced.append(dict(row, screeningrise=predicted))
        # Reuse the real scorer without changing its frozen four-cause registry.
        item = score(replaced, covariance)[0]
        result.append(dict(item, cause=name, interpretation='screening_submodel_not_new_cause'))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', type=Path)
    parser.add_argument('log_covariance', type=Path)
    parser.add_argument('profiles', type=Path)
    parser.add_argument('laws', type=Path)
    args = parser.parse_args()
    with args.runs.open() as stream:
        rows = list(csv.DictReader(stream))
    with np.load(args.log_covariance, allow_pickle=False) as data:
        ids = [str(value) for value in data['run_id'].tolist()]
        if len(ids) != len(rows) or len(set(ids)) != len(ids) or set(ids) != {r['run_id'] for r in rows}:
            raise ValueError('Same complete, consistently masked run IDs required')
        index = {run_id: i for i, run_id in enumerate(ids)}
        permutation = [index[row['run_id']] for row in rows]
        cov = data['covariance']
        if cov.shape != (len(ids), len(ids)):
            raise ValueError('Covariance dimensions must match its run IDs')
        covariance = cov[np.ix_(permutation, permutation)]
    profiles = json.loads(args.profiles.read_text())
    laws = json.loads(args.laws.read_text())
    print(json.dumps(compare_runs(rows, covariance, profiles, laws), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
