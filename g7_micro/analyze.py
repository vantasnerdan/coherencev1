"""Bounded comparative fits to digitized normalized Figure 2, never raw counts.

Run from checkout root: python -m g7_micro.analyze
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import subprocess
import sys
import tempfile

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
from scipy.sparse import lil_matrix

from g4_beam.forward import model as DD
from g4_beam.forward.consumer import forward_yields
from g4_beam.forward.profile_thermal import transport_stopping
from g4_beam.merit.api import ENERGIES, bestfits, reference
from g6_plateau.forecast import projected_yield
from g7_micro.model import BULK, FrozenResponse, U_GRID, LOG_GRID, population, raw_depth_yield

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FAMILIES = ('sigmoid', 'anchored_rise', 'two_site', 'exponential')


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class Benchmark:
    def __init__(self):
        self.data = DD.load_data(str(ROOT/'g4_beam/forward/errbars.csv'))
        self.baseline = bestfits()['screeningrise']
        self.E0 = self.baseline['params']['E0']
        self.responses = {m: [FrozenResponse(e, reference(m).profiles(float(e))) for e in ENERGIES]
                          for m in DD.MATS}
        self.indices = {m: np.array([np.flatnonzero(ENERGIES == e)[0] for e in self.data[m][0]])
                        for m in DD.MATS}
        self.sparsity = lil_matrix((sum(len(self.data[m][0]) for m in DD.MATS), 9), dtype=int)
        row = 0
        for i, m in enumerate(DD.MATS):
            n = len(self.data[m][0])
            self.sparsity[row:row+n, 0] = 1
            self.sparsity[row:row+n, 1+2*i:3+2*i] = 1
            row += n
        self.population_tables = {}

    def population_values(self, material, family, a, b, sheet, order):
        if material not in self.population_tables:
            self.population_tables[material] = [
                PchipInterpolator(LOG_GRID, np.log(np.array([
                    response.components(U_GRID, ch) for response in self.responses[material]])), axis=-1)
                for ch in ('p', 'n')]
        u, probability = population(family, a, b, order)
        values = np.array([np.exp(fn(np.log(u+.005)))@probability
                           for fn in self.population_tables[material]])
        if sheet == 'bulk40':
            values[:, :, 1] = np.array([np.exp(fn(np.log(BULK+.005)))[:, 1]
                                       for fn in self.population_tables[material]])
        return np.sum(values, axis=2)

    def curve(self, material, family, a, b, power, stitch, sheet, exact=False, order=96):
        if family in ('two_site', 'exponential') and not exact:
            values = self.population_values(material, family, a, b, sheet, order)
        else:
            values = np.array([[response.yield_per_C(family, a, b, channel=ch, sheet=sheet,
                                                    exact=exact, E0=self.E0, order=order)
                                for response in self.responses[material]] for ch in ('p', 'n')])
        correction = DD.collected_fraction(ENERGIES, 1, 0)*(DD.I_cup(ENERGIES)*1e-6)**power
        values /= correction
        ip = int(np.flatnonzero(ENERGIES == 4.5)[0])
        low = ENERGIES <= 4.
        prediction = values[1]/values[1, -1]
        prediction[low] = values[0, low]*values[1, ip]/values[0, ip]*stitch/values[1, -1]
        return prediction

    def decode(self, family, x):
        if family in ('two_site', 'exponential'):
            return 10.**x[0], 10.**x[1]
        if family == 'sigmoid':
            return x[0], x[1]
        return x[0], x[1]

    def residual(self, x, family, power, sheet, *, high_only=False):
        result = []
        for i, m in enumerate(DD.MATS):
            a, b = self.decode(family, x[1+2*i:3+2*i])
            pred = self.curve(m, family, a, b, power, np.exp(x[0]), sheet)
            e, observed, error = self.data[m]
            residual = (np.log(pred[self.indices[m]])-np.log(observed))/error
            result.extend(residual[e >= 2.5] if high_only else residual)
        return np.array(result)

    def fit(self, family, power, sheet, *, high_only=False):
        params = self.baseline['params']
        if family == 'sigmoid':
            pairs = [[params['Ue_'+m], params['Ues_'+m]] for m in DD.MATS]
            bounds = ([0., 0.], [15., 30.])
        elif family == 'anchored_rise':
            pairs = [[params['Ues_'+m]-BULK, self.E0] for m in DD.MATS]
            bounds = ([0., .05], [30., 3.])
        else:
            pairs = [[-2., np.log10(4.)] for _ in DD.MATS]
            bounds = ([-12., np.log10(BULK) if family == 'two_site' else -4.], [0., 6.])
        lower = np.r_[np.log(.3), np.tile(bounds[0], 4)]
        upper = np.r_[np.log(3.), np.tile(bounds[1], 4)]
        x0 = np.r_[np.log(params['s']), np.ravel(pairs)]
        candidates = []
        starts = [x0]
        if family in ('two_site', 'exponential'):
            high = x0.copy(); high[1::2] = -4.; high[2::2] = 2.
            starts.append(high)
        for start in starts:
            fit = least_squares(self.residual, np.clip(start, lower+1e-10, upper-1e-10),
                                args=(family, power, sheet), kwargs={'high_only': high_only},
                                bounds=(lower, upper), jac_sparsity=None if high_only else self.sparsity,
                                x_scale='jac', max_nfev=180, ftol=2e-7, xtol=2e-7, gtol=2e-6)
            candidates.append(fit)
        fit = min(candidates, key=lambda fit: fit.fun@fit.fun)
        material = {}
        predictions = []
        for i, m in enumerate(DD.MATS):
            a, b = self.decode(family, fit.x[1+2*i:3+2*i])
            curve = self.curve(m, family, a, b, power, np.exp(fit.x[0]), sheet, exact=True)
            e, obs, err = self.data[m]
            residual = (np.log(curve[self.indices[m]])-np.log(obs))/err
            material[m] = dict(a=float(a), b=float(b), chi2=float(residual@residual),
                               low_chi2=float(residual[e < 2.5]@residual[e < 2.5]),
                               high_chi2=float(residual[e >= 2.5]@residual[e >= 2.5]),
                               n=len(e), low_n=int(sum(e < 2.5)))
            for j, energy in enumerate(e):
                predictions.append(dict(material=m, E_dd_cm_keV=float(energy),
                                        observed_R=float(obs[j]), sigma_lnR=float(err[j]),
                                        predicted_R=float(curve[self.indices[m][j]])))
        return dict(family=family, power=power, sheet=sheet, high_only=high_only,
                    E0_external_keV=self.E0 if family == 'sigmoid' else None,
                    stitch=float(np.exp(fit.x[0])), k=9, n=sum(len(v[0]) for v in self.data.values()),
                    chi2=sum(v['chi2'] for v in material.values()),
                    surrogate_chi2=float(fit.fun@fit.fun), success=bool(fit.success),
                    nfev=int(fit.nfev), message=fit.message, materials=material,
                    predictions=predictions, encoded_x=fit.x.tolist())

    def profile(self, fitted):
        """2D (f,U) region, profiling common stitch and every other material.

        Conditional-on-comparand fixed E0/power/sheet; descriptive delta only.
        At each grid cell other materials retain their already-optimized shapes
        while common log stitch is analytically reprofiled. Additionally, the
        one-dimensional U profile reoptimizes every nuisance shape and stitch.
        """
        power, sheet = fitted['power'], fitted['sheet']
        base = np.array(fitted['encoded_x'])
        ugrid = np.unique(np.r_[np.geomspace(.06, 30., 46), np.geomspace(40., 1e6, 25),
                               [v['b'] for v in fitted['materials'].values()], DD.BG**2])
        fgrid = np.unique(np.r_[0., np.geomspace(1e-10, 1., 81),
                               [v['a'] for v in fitted['materials'].values()]])
        rows, regions = [], []
        for mat in ('PdD', 'TiD'):
            mi = DD.MATS.index(mat)
            start = base.copy()
            one_profile = []
            for u in ugrid:
                fixed = 2+2*mi
                free = np.array([j for j in range(9) if j != fixed])
                full = start.copy(); full[fixed] = np.log10(u)
                lo = np.r_[np.log(.3), np.tile([-12., np.log10(BULK)], 4)][free]
                hi = np.r_[np.log(3.), np.tile([0., 6.], 4)][free]
                def res(x):
                    point = full.copy(); point[free] = x
                    return self.residual(point, 'two_site', power, sheet)
                fit = least_squares(res, np.clip(full[free], lo+1e-10, hi-1e-10), bounds=(lo, hi),
                                    x_scale='jac', max_nfev=70, ftol=2e-6, xtol=2e-6, gtol=1e-5)
                point = full.copy(); point[free] = fit.x
                start = point
                one_profile.append((u, float(fit.fun@fit.fun), 10**point[1+2*mi], point.copy()))
            minimum = min(fitted['chi2'], min(row[1] for row in one_profile))
            for u, chi2, f, point in one_profile:
                rows.append(dict(material=mat, kind='profile_U_all_nuisance', U_high_keV=float(u),
                                 f=float(f), chi2=chi2, delta_chi2=chi2-minimum,
                                 stitch=float(np.exp(point[0]))))
            # Conditional grid is useful for occupancy translation, not a mislabeled full profile.
            other_res = self.residual(base, 'two_site', power, sheet)
            offsets = np.cumsum([0]+[len(self.data[m][0]) for m in DD.MATS])
            keep = np.ones(len(other_res), bool); keep[offsets[mi]:offsets[mi+1]] = False
            e, observed, error = self.data[mat]
            lows = np.concatenate([self.data[m][0] <= 4. for m in DD.MATS])
            errors = np.concatenate([self.data[m][2] for m in DD.MATS])
            design = lows/errors
            fixed_other = other_res.copy()-design*base[0]
            for u in ugrid:
                for f in fgrid:
                    prediction = self.curve(mat, 'two_site', float(f), float(u), power, 1., sheet)
                    res = fixed_other.copy()
                    res[~keep] = (np.log(prediction[self.indices[mat]])-np.log(observed))/error
                    logs = np.clip(-float(design@res)/float(design@design), np.log(.3), np.log(3.))
                    chi2 = float((res+logs*design)@(res+logs*design))
                    rows.append(dict(material=mat, kind='grid_other_shapes_fixed_stitch_profiled',
                                     U_high_keV=float(u), f=float(f), chi2=chi2,
                                     delta_chi2=chi2-minimum, stitch=float(np.exp(logs))))
            for threshold in (2.30, 5.99, 11.83):
                selected = [r for r in rows if r['material'] == mat and
                            r['kind'] == 'grid_other_shapes_fixed_stitch_profiled' and
                            r['delta_chi2'] <= threshold]
                prof = [r for r in rows if r['material'] == mat and r['kind'] == 'profile_U_all_nuisance'
                        and r['delta_chi2'] <= threshold]
                regions.append(dict(material=mat, descriptive_delta_chi2=threshold,
                                    grid_points=len(selected),
                                    conditional_f_min=min((r['f'] for r in selected), default=None),
                                    conditional_f_max=max((r['f'] for r in selected), default=None),
                                    conditional_U_min=min((r['U_high_keV'] for r in selected), default=None),
                                    conditional_U_max=max((r['U_high_keV'] for r in selected), default=None),
                                    full_profile_U_min=min((r['U_high_keV'] for r in prof), default=None),
                                    full_profile_U_max=max((r['U_high_keV'] for r in prof), default=None),
                                    upper_grid_boundary_allowed=any(r['U_high_keV'] == 1e6 for r in prof)))
        return rows, regions


def verify_adapter():
    """Deterministic arithmetic fixture, not observations or simulated raw data."""
    from g6_plateau.score_runs import score
    from g7_micro.score_runs import predict_run
    rows, profiles = [], {}
    for i, e in enumerate((.5, 1.)):
        layers, nads, state = reference('PdD').profiles(e)
        profile_id = 'arithmetic_profile_'+str(i)
        profiles[profile_id] = dict(nominal_DD_cm_key_keV=e, layers=layers, nads_m2=nads,
                                    state=state, source_id='g4_beam/merit/profiles.json')
        y = float(forward_yields([e], reference('PdD').profiles, Ue_keV=BULK,
                                reference_E_cm_keV=None, stopping_fn=transport_stopping)
                  ['physical_yield_per_C'][0])
        rows.append(dict(run_id='arithmetic_'+str(i), normalization_group='shared_fixture',
                         observed_yield_per_true_C=y, screeningrise=y, surface=y, neutrals=y, lattice=y,
                         material='PdD', lab_keV=2*e, isotope='DD', channel='p', D1_fraction=.85,
                         D2_fraction=.15, spread_FWHM=.05, profile_id=profile_id,
                         yield_convention='reactions_per_true_C'))
    laws = {'fixture_bulk': {'PdD': dict(family='two_site', a=0., b=3., sheet='shared')},
            'fixture_point3': {'PdD': dict(family='two_site', a=1., b=3., sheet='shared')},
            'fixture_depth_bulk': {'PdD': dict(family='depth', f_s=0., f_b=0., d_nm=3.,
                                              U_high_keV=3., sheet='shared')}}
    covariance = np.array([[.04, .01], [.01, .09]])
    with tempfile.TemporaryDirectory(prefix='g7_micro_arithmetic_') as directory:
        path = Path(directory)
        write_csv(path/'runs.csv', rows)
        (path/'profiles.json').write_text(json.dumps(profiles))
        (path/'laws.json').write_text(json.dumps(laws))
        # Deliberately reversed covariance IDs exercises the actual CLI join.
        np.savez(path/'covariance.npz', run_id=np.array([r['run_id'] for r in rows][::-1]),
                 covariance=covariance[::-1, ::-1])
        process = subprocess.run([sys.executable, '-m', 'g7_micro.score_runs',
                                  str(path/'runs.csv'), str(path/'covariance.npz'),
                                  str(path/'profiles.json'), str(path/'laws.json')],
                                 cwd=ROOT, capture_output=True, text=True, check=True)
        result = {r['cause']: r for r in json.loads(process.stdout)}
    paired = [dict(row, screeningrise=predict_run(row, profiles[row['profile_id']],
                                                 laws['fixture_point3']['PdD'])) for row in rows]
    expected = score(paired, covariance)[0]['chi2']
    observed = result['fixture_point3']['chi2']
    if abs(expected-observed) > 1e-10 or result['fixture_bulk']['chi2'] > 1e-20 or expected <= 0:
        raise AssertionError(result)
    return dict(real_cli='python -m g7_micro.score_runs', fixture_only=True,
                covariance_ID_permutation=True, original_score_chi2=score(rows, covariance)[0]['chi2'],
                identity_distribution_chi2=result['fixture_bulk']['chi2'],
                depth_identity_chi2=result['fixture_depth_bulk']['chi2'],
                paired_point3_chi2=observed, paired_shared_scorer_abs_error=abs(expected-observed))


def verify(benchmark, fits):
    """One settled consumer check: scalar identities, paired ranking and quadrature."""
    checks = {}
    target = reference('PdD')
    e = .5
    response = benchmark.responses['PdD'][int(np.flatnonzero(ENERGIES == e)[0])]
    for u in (BULK, 3.):
        expected = forward_yields([e], target.profiles, Ue_keV=u, reference_E_cm_keV=None,
                                  stopping_fn=transport_stopping)['physical_yield_per_C'][0]
        actual = response.yield_per_C('two_site', 1., u, exact=True)
        checks['DD_point_'+str(u)+'_relative'] = abs(actual/expected-1)
    zero = response.yield_per_C('two_site', 0., 8., exact=True)
    bulk = response.yield_per_C('two_site', 1., BULK, exact=True)
    checks['f_zero_identity_relative'] = abs(zero/bulk-1)
    dt = FrozenResponse(e, target.profiles(e), isotope='DT')
    expected = projected_yield(e, target.profiles, 3., 'none', {}, {}, isotope='DT')
    actual = dt.yield_per_C('two_site', 1., 3., exact=True)
    checks['DT_point_relative'] = abs(actual/expected-1)
    interpolation, tail_error = 0., 0.
    for fit in fits:
        if fit['high_only']:
            continue
        interpolation = max(interpolation, abs(fit['chi2']-fit['surrogate_chi2']))
        if fit['family'] == 'exponential':
            for mat in ('PdD', 'TiD'):
                a, b = (fit['materials'][mat][k] for k in ('a', 'b'))
                for response in benchmark.responses[mat]:
                    low = response.yield_per_C('exponential', a, b, sheet=fit['sheet'], exact=True)
                    high = response.yield_per_C('exponential', a, b, sheet=fit['sheet'], exact=True, order=192)
                    tail_error = max(tail_error, abs(low/high-1))
    checks.update(max_abs_chi2_interpolation_error=interpolation,
                  max_fitted_tail_relative_96_vs_192=tail_error,
                  raw_observations_used=False, synthetic_counts_used=False)
    checks['run_adapter'] = verify_adapter()
    if max(checks[k] for k in ('DD_point_0.04_relative', 'DD_point_3.0_relative',
                               'f_zero_identity_relative', 'DT_point_relative')) > 2e-12:
        raise AssertionError(checks)
    if interpolation > .05 or tail_error > .002:
        raise AssertionError(checks)
    return checks


def dt_table(fits, depth=None):
    rows = []
    for fit in fits:
        if fit['high_only'] or fit['power'] != 1 or fit['family'] not in ('two_site', 'exponential'):
            continue
        for mat in ('PdD', 'TiD'):
            target = reference(mat)
            a, b = (fit['materials'][mat][k] for k in ('a', 'b'))
            for lab in (1., 2., 5., 13.):
                key = lab/2
                response = FrozenResponse(key, target.profiles(key), isotope='DT')
                mixed = response.yield_per_C(fit['family'], a, b, sheet=fit['sheet'], exact=True)
                moment = response.yield_per_C(fit['family'], a, b, sheet=fit['sheet'], exact=True, moment=True)
                bare = sum(response.components([0.], exact=True)[:, 0])/1.602176634e-19
                rows.append(dict(material=mat, family=fit['family'], sheet=fit['sheet'], lab_keV=lab,
                                 profile_nominal_DD_cm_key_keV=key, mixed_reactions_per_true_C=mixed,
                                 bare_reactions_per_true_C=bare, screened_over_bare=mixed/bare,
                                 yield_weighted_U_keV=moment/mixed, f=a, tail_scale_or_U_high_keV=b,
                                 f_s=None, f_b=None, d_nm=None))
    if depth is not None:
        for selection in depth['best']:
            if selection['power'] != 1 or selection['sheet'] != 'shared':
                continue
            fit = selection['best']
            mat = selection['material']
            target = reference(mat)
            for lab in (1., 2., 5., 13.):
                key = lab/2
                response = FrozenResponse(key, target.profiles(key), isotope='DT')
                parameters = dict(f_s=fit['f_s'], f_b=fit['f_b'], d_nm=fit['d_nm'],
                                  U_high_keV=fit['U_high_keV'], isotope='DT', sheet='shared')
                mixed = raw_depth_yield(key, response.profile, **parameters)
                moment = raw_depth_yield(key, response.profile, **parameters, moment=True)
                bare = float(np.sum(response.components([0.], exact=True))/1.602176634e-19)
                # A uniform depth law must match the stationary yield moment.
                uniform = raw_depth_yield(key, response.profile, **dict(parameters, f_s=.03, f_b=.03),
                                          moment=True)
                expected = response.yield_per_C('two_site', .03, fit['U_high_keV'], exact=True, moment=True)
                if abs(uniform/expected-1) > 2e-6:
                    raise AssertionError('DT depth moment fails uniform-population identity')
                rows.append(dict(material=mat, family='depth', sheet='shared', lab_keV=lab,
                                 profile_nominal_DD_cm_key_keV=key, mixed_reactions_per_true_C=mixed,
                                 bare_reactions_per_true_C=bare, screened_over_bare=mixed/bare,
                                 yield_weighted_U_keV=moment/mixed, f=None,
                                 tail_scale_or_U_high_keV=fit['U_high_keV'],
                                 f_s=fit['f_s'], f_b=fit['f_b'], d_nm=fit['d_nm']))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=HERE)
    parser.add_argument('--reuse-fits', action='store_true',
                        help='Reuse this output directory results.json after a focused correction')
    parser.add_argument('--dt-only', action='store_true',
                        help='Rebuild only conditional DT table from saved stationary and depth fits')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    if args.dt_only:
        fits = json.loads((args.out/'results.json').read_text())['fits']
        depth = json.loads((args.out/'depth_results.json').read_text())
        rows = dt_table(fits, depth)
        write_csv(args.out/'dt.csv', rows)
        print(json.dumps({'DT_rows': len(rows), 'depth_rows': sum(r['family'] == 'depth' for r in rows),
                          'DT_uniform_depth_moment_identity': 'passed', 'elapsed_seconds': time.monotonic()-start}))
        return
    benchmark = Benchmark()
    if args.reuse_fits:
        fits = json.loads((args.out/'results.json').read_text())['fits']
    else:
        fits = []
        for sheet in ('shared', 'bulk40'):
            for power in (0, 1, 2):
                for family in FAMILIES:
                    fit = benchmark.fit(family, power, sheet)
                    fits.append(fit)
                    print(f"{sheet} power={power} {family} chi2={fit['chi2']:.6f} success={fit['success']}", flush=True)
        # Diagnostic holdout only: external E0 used these curves previously.
        for family in FAMILIES:
            fits.append(benchmark.fit(family, 1, 'shared', high_only=True))
        (args.out/'results.json').write_text(json.dumps({'fits': fits}, indent=2, allow_nan=False)+'\n')
    checks = verify(benchmark, fits)
    two = next(f for f in fits if f['family'] == 'two_site' and f['power'] == 1 and
               f['sheet'] == 'shared' and not f['high_only'])
    profiles, regions = benchmark.profile(two)
    from g7_micro.depth import evaluate as evaluate_depth
    depth = evaluate_depth(benchmark, args.out)
    source_paths = ['g4_beam/merit/api.py', 'g4_beam/merit/refit.py', 'g4_beam/merit/profiles.json',
                    'g4_beam/merit/fits_transported.json', 'g4_beam/forward/errbars.csv',
                    'g4_beam/forward/model.py', 'g4_beam/forward/consumer.py',
                    'g4_beam/forward/profile_thermal.py', 'g6_plateau/forecast.py', 'g6_plateau/score_runs.py']
    summary = dict(source_main='2ca008af1ad630056e65d2e7c19452c668bd3c87',
                   source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths},
                   baseline=benchmark.baseline, external_E0_keV=benchmark.E0, fits=fits, regions=regions,
                   verification=checks, depth=depth, elapsed_seconds=time.monotonic()-start)
    (args.out/'results.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    table = []
    for fit in fits:
        for mat, value in fit['materials'].items():
            table.append(dict(family=fit['family'], power=fit['power'], sheet=fit['sheet'],
                              high_only=fit['high_only'], material=mat, a=value['a'], b=value['b'],
                              stitch=fit['stitch'], total_chi2=fit['chi2'], material_chi2=value['chi2'],
                              low_chi2=value['low_chi2'], high_chi2=value['high_chi2'],
                              n=value['n'], low_n=value['low_n'], k=fit['k'], success=fit['success']))
    write_csv(args.out/'fits.csv', table)
    write_csv(args.out/'profile.csv', profiles)
    write_csv(args.out/'regions.csv', regions)
    write_csv(args.out/'dt.csv', dt_table(fits, depth))
    print(json.dumps(checks, indent=2), flush=True)
    print('wall_seconds', time.monotonic()-start, flush=True)


if __name__ == '__main__':
    main()
