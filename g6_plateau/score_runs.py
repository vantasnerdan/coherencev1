"""Compare the four frozen-U hypotheses on a normalized g4_zenodo run join.

Input CSV and keyed log covariance are defined in README.md. This is the
calibrated positive-count shape comparison, not a raw-count/Poisson decoder.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky, solve_triangular

CAUSES = ('screeningrise', 'surface', 'neutrals', 'lattice')


def score(runs, covariance):
    ids = [row['run_id'] for row in runs]
    if not ids or any(not run_id for run_id in ids) or len(set(ids)) != len(ids):
        raise ValueError('run_id must identify every scored observation exactly once')
    groups = [row['normalization_group'] for row in runs]
    free = sorted(set(groups)-{'absolute'})
    if any(not group for group in groups):
        raise ValueError('Use absolute or an explicit shared normalization group')
    df = len(runs)-len(free)
    if df <= 0:
        raise ValueError('Group normalizations absorb all per-run shape information')
    values = np.array([[float(row[k]) for k in ('observed_yield_per_true_C', *CAUSES)] for row in runs])
    if not np.all(np.isfinite(values) & (values > 0)):
        raise ValueError('Require positive finite calibrated yields; use a raw-count likelihood otherwise')
    if covariance.shape != (len(ids), len(ids)) or not np.isfinite(covariance).all() or not np.allclose(covariance, covariance.T, rtol=1e-12, atol=0.):
        raise ValueError('Supply one finite, symmetric log-yield covariance on the scored rows')
    ll = cholesky(covariance, lower=True)
    log_ratio = np.log(values[:, [0]]) - np.log(values[:, 1:])
    residuals = solve_triangular(ll, log_ratio, lower=True)
    design = np.array([[group == name for name in free] for group in groups], dtype=float)
    weighted = solve_triangular(ll, design, lower=True)
    if free:
        normalization, _, rank, _ = np.linalg.lstsq(weighted, residuals, rcond=None)
        if rank != len(free):
            raise ValueError('Rank-deficient normalization design')
        residuals -= weighted@normalization
    else:
        normalization = np.empty((0, len(CAUSES)))
    return [dict(cause=cause, chi2=float(residuals[:, j]@residuals[:, j]), dof=df,
                 scales={name: float(np.exp(normalization[k, j])) for k, name in enumerate(free)})
            for j, cause in enumerate(CAUSES)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', type=Path)
    parser.add_argument('log_covariance', type=Path)
    args = parser.parse_args()
    with args.runs.open() as f:
        rows = list(csv.DictReader(f))
    with np.load(args.log_covariance, allow_pickle=False) as data:
        ids = [str(value) for value in data['run_id'].tolist()]
        if len(ids) != len(rows) or len(set(ids)) != len(ids) or set(ids) != {r['run_id'] for r in rows}:
            raise ValueError('Covariance and consistently masked run join must cover the same run IDs')
        index = {run_id: i for i, run_id in enumerate(ids)}
        permutation = [index[row['run_id']] for row in rows]
        cov = data['covariance']
        if cov.shape != (len(ids), len(ids)):
            raise ValueError('Covariance dimensions must match its run IDs')
        scores = score(rows, cov[np.ix_(permutation, permutation)])
    print(json.dumps(scores, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
