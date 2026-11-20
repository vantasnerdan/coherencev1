"""Exact two-zone populations: ratio-feasible region and full-curve cost.

Depth is a spatial population assumption, not an invented phase assignment.
The two low-energy central values define one linear constraint, not a likelihood.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar, brentq

from g4_beam.forward import model as DD
from g4_beam.merit.api import ENERGIES
from g4_beam.forward.consumer import csda_nodes
from g4_beam.forward.profile_thermal import transport_stopping
from g7_micro.analyze import Benchmark, write_csv
from g7_micro.model import BULK, depth_components, raw_depth_yield

HERE = Path(__file__).resolve().parent
DEPTH_NM = np.array([1., 2., 3., 5., 10., 20., 30., 50.])
SHIFTS = np.unique(np.r_[np.geomspace(.08, 30., 54), np.geomspace(40., 1e6, 25), DD.BG**2])


def segment(c, a, b):
    """Intersect c+a*f_s+b*f_b=0 with triangle 0<=f_b<=f_s<=1."""
    points = []
    scale = max(abs(c), abs(a), abs(b))
    if scale == 0:
        return [(0., 0.), (1., 1.)]
    for origin, direction in ((np.array([0., 0.]), np.array([1., 0.])),
                              (np.array([1., 0.]), np.array([0., 1.])),
                              (np.array([0., 0.]), np.array([1., 1.]))):
        denominator = a*direction[0]+b*direction[1]
        if abs(denominator) < 1e-14*scale:
            continue
        t = -(c+a*origin[0]+b*origin[1])/denominator
        if -1e-12 <= t <= 1+1e-12:
            point = origin+np.clip(t, 0., 1.)*direction
            if not any(np.linalg.norm(point-old) < 1e-10 for old in points):
                points.append(point)
    return points


def inventory(response, depth_nm):
    """Encountered number inventories per area, with no nuclear/yield weighting."""
    fuel = top = metal = sheet = 0.
    for lab, share, _, _ in response.tracks:
        _, _, _, end = csda_nodes(lab, response.layers, transport_stopping)
        for zt, zb, _, nd, _, _, nm in response.layers:
            length = max(0., min(zb, end)-zt)
            local = max(0., min(zb, end, depth_nm*1e-9)-zt)
            fuel += share*length*nd
            top += share*local*nd
            metal += share*length*nm
        sheet += share*response.nads
    return dict(fuel_per_metal=fuel/metal, top_fuel_per_metal=top/metal,
                sheet_fuel_per_metal=sheet/metal, sheet_fraction=sheet/(fuel+sheet))


def allowed_segment(cost, optimum_t, threshold):
    grid = np.unique(np.r_[np.linspace(0., 1., 81), optimum_t])
    objective = np.array([cost(t)-threshold for t in grid])
    allowed = grid[objective <= 0].tolist()
    for i in range(len(grid)-1):
        if objective[i]*objective[i+1] < 0:
            allowed.append(brentq(lambda t: cost(t)-threshold, grid[i], grid[i+1], xtol=1e-12))
    return allowed


def curve_band(cells, material, out_rows):
    """Profile depth, deeper population, sheet convention and one-foil stitch.

    The retained metal-normalized quantity is occupied special D, a necessary
    lower bound on available single-occupancy sites, not an occupancy inversion.
    Delta chi-square is descriptive under the digitized approximate errors.
    """
    minimum = min(cell['row']['full_chi2_one_foil'] for cell in cells)
    threshold = minimum+4.
    for shift in SHIFTS:
        selected = [cell for cell in cells if cell['row']['U_high_keV'] == shift]
        if not selected:
            continue
        ratio_low, ratio_high, band, volume_band, depths, sheets = [], [], [], [], [], []
        for cell in selected:
            start, end, inv, row, cost = (cell[k] for k in ('start', 'end', 'inventory', 'row', 'cost'))
            def occupied(t, include_sheet=True):
                fs, fb = start+t*(end-start)
                return float(fs*(inv['top_fuel_per_metal']+
                                 (inv['sheet_fuel_per_metal'] if include_sheet and row['sheet'] == 'shared' else 0.))+
                             fb*(inv['fuel_per_metal']-inv['top_fuel_per_metal']))
            limits = [occupied(0.), occupied(1.)]
            ratio_low.append(min(limits)); ratio_high.append(max(limits))
            if row['full_chi2_one_foil'] > threshold:
                continue
            allowed = allowed_segment(cost, cell['optimum_t'], threshold)
            band.extend(occupied(t) for t in allowed)
            volume_band.extend(occupied(t, False) for t in allowed)
            depths.append(row['d_nm']); sheets.append(row['sheet'])
        out_rows.append(dict(material=material, U_high_keV=float(shift), power=1,
                             profile_nominal_DD_cm_key_keV=1.,
                             ratio_only_occupied_min_per_metal=min(ratio_low),
                             ratio_only_occupied_max_per_metal=max(ratio_high),
                             full_min_chi2_one_foil=minimum, descriptive_delta_chi2=4.,
                             full_band_occupied_min_per_metal=min(band) if band else None,
                             full_band_occupied_max_per_metal=max(band) if band else None,
                             full_band_volume_occupied_min_per_metal=min(volume_band) if volume_band else None,
                             full_band_volume_occupied_max_per_metal=max(volume_band) if volume_band else None,
                             full_band_d_min_nm=min(depths) if depths else None,
                             full_band_d_max_nm=max(depths) if depths else None,
                             full_band_sheet_choices=';'.join(sorted(set(sheets)))))


def evaluate(benchmark, out):
    rows, best, verification, curves, inventories = [], [], {}, [], []
    for material in ('PdD', 'TiD'):
        energies, observed, errors = benchmark.data[material]
        ratio = float(observed[energies == .25][0]/observed[energies == 1.][0])
        i0, i1 = [int(np.flatnonzero(ENERGIES == e)[0]) for e in (.25, 1.)]
        ip = int(np.flatnonzero(ENERGIES == 4.5)[0])
        low = ENERGIES <= 4.
        idx = benchmark.indices[material]
        design = (energies <= 4.)/errors
        cells = []
        for depth in DEPTH_NM:
            inv = inventory(benchmark.responses[material][i1], depth)
            inventories.append(dict(material=material, profile_nominal_DD_cm_key_keV=1.,
                                    d_nm=float(depth), **inv))
            components = np.array([[depth_components(response, depth, SHIFTS, channel)
                                    for response in benchmark.responses[material]]
                                   for channel in ('p', 'n')])
            for sheet in ('shared', 'bulk40'):
                for power in (0, 1, 2):
                    correction = DD.collected_fraction(ENERGIES, 1, 0)*(DD.I_cup(ENERGIES)*1e-6)**power
                    values = components/correction[None, :, None, None]
                    for ui, shift in enumerate(SHIFTS):
                        base = values[:, :, ui, 0]
                        top = values[:, :, ui, 1]+(values[:, :, ui, 3] if sheet == 'shared' else 0.)
                        deep = values[:, :, ui, 2]
                        feasible = segment(base[0, i0]-ratio*base[0, i1],
                                           top[0, i0]-ratio*top[0, i1],
                                           deep[0, i0]-ratio*deep[0, i1])
                        if not feasible:
                            continue
                        start = np.array(feasible[0]); end = np.array(feasible[-1])
                        def cost(t, details=False, base=base, top=top, deep=deep, start=start, end=end):
                            fs, fb = start+t*(end-start)
                            yields = base+fs*top+fb*deep
                            predicted = yields[1]/yields[1, -1]
                            predicted[low] = yields[0, low]*yields[1, ip]/yields[0, ip]/yields[1, -1]
                            residual = (np.log(predicted[idx])-np.log(observed))/errors
                            logs = np.clip(-float(design@residual)/float(design@design), np.log(.3), np.log(3.))
                            residual += logs*design
                            return (float(residual@residual), float(fs), float(fb), float(np.exp(logs)),
                                    float(yields[0, i0]/yields[0, i1])) if details else float(residual@residual)
                        fit = minimize_scalar(cost, bounds=(0., 1.), method='bounded', options={'xatol': 1e-7})
                        t = min((0., 1., float(fit.x)), key=cost)
                        chi2, fs, fb, stitch, achieved = cost(t, True)
                        rows.append(dict(material=material, power=power, sheet=sheet, d_nm=float(depth),
                                         U_high_keV=float(shift), f_s=fs, f_b=fb,
                                         feasible_f_s_min=float(min(p[0] for p in feasible)),
                                         feasible_f_s_max=float(max(p[0] for p in feasible)),
                                         feasible_f_b_min=float(min(p[1] for p in feasible)),
                                         feasible_f_b_max=float(max(p[1] for p in feasible)),
                                         R025_over_R1=achieved, observed_ratio=ratio,
                                         full_chi2_one_foil=chi2, stitch_one_foil=stitch,
                                         n=len(energies)))
                        if power == 1:
                            cells.append(dict(start=start, end=end, inventory=inv, row=rows[-1],
                                              cost=cost, optimum_t=t))
        curve_band(cells, material, curves)
        for sheet in ('shared', 'bulk40'):
            for power in (0, 1, 2):
                candidates = [r for r in rows if r['material'] == material and r['sheet'] == sheet and
                              r['power'] == power]
                best.append(dict(material=material, sheet=sheet, power=power, feasible_cells=len(candidates),
                                 best=min(candidates, key=lambda r: r['full_chi2_one_foil']) if candidates else None))
        # A uniform two-zone law must reduce exactly to the stationary scalar law.
        response = benchmark.responses[material][i1]
        reference_value = response.yield_per_C('two_site', .03, 3., exact=True)
        depth_value = raw_depth_yield(response.e, response.profile, f_s=.03, f_b=.03,
                                      d_nm=3., U_high_keV=3.)
        verification[material+'_uniform_depth_relative'] = abs(depth_value/reference_value-1)
    if max(verification.values()) > 2e-6:
        raise AssertionError(verification)
    if rows:
        write_csv(out/'depth_region.csv', rows)
    write_csv(out/'depth_inventory_curves.csv', curves)
    write_csv(out/'inventories.csv', inventories)
    summary = dict(best=best, verification=verification, depths_nm=DEPTH_NM.tolist(),
                   U_high_scan_keV=SHIFTS.tolist(),
                   interpretation='Central-ratio feasible segments; not confidence or posterior regions. '
                                  'Full chi2 profiles one independent stitch per foil, not the original shared stitch. '
                                  'No phase fractions were inferred from unmeasured phase metadata.')
    (out/'depth_results.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    return summary


def add_volume_bands(benchmark, out):
    """Postprocess only saved full-band cells; never repeat the depth fit scan."""
    with (out/'depth_region.csv').open() as stream:
        region = list(csv.DictReader(stream))
    with (out/'depth_inventory_curves.csv').open() as stream:
        curves = list(csv.DictReader(stream))
    volumes = {}
    for material in ('PdD', 'TiD'):
        e, observed, errors = benchmark.data[material]
        idx = benchmark.indices[material]
        i0, i1, ip = [int(np.flatnonzero(ENERGIES == v)[0]) for v in (.25, 1., 4.5)]
        low = ENERGIES <= 4.
        design = (e <= 4.)/errors
        correction = DD.collected_fraction(ENERGIES, 1, 0)*(DD.I_cup(ENERGIES)*1e-6)
        threshold = next(float(r['full_min_chi2_one_foil'])+4. for r in curves if r['material'] == material)
        ratio = float(observed[e == .25][0]/observed[e == 1.][0])
        accepted = [r for r in region if r['material'] == material and int(r['power']) == 1 and
                    float(r['full_chi2_one_foil']) <= threshold]
        for row in accepted:
            shift, depth = float(row['U_high_keV']), float(row['d_nm'])
            values = np.array([[depth_components(response, depth, [shift], channel)[0]
                                for response in benchmark.responses[material]] for channel in ('p', 'n')])
            values /= correction[None, :, None]
            base, deep = values[:, :, 0], values[:, :, 2]
            top = values[:, :, 1]+(values[:, :, 3] if row['sheet'] == 'shared' else 0.)
            points = segment(base[0, i0]-ratio*base[0, i1], top[0, i0]-ratio*top[0, i1],
                             deep[0, i0]-ratio*deep[0, i1])
            start, end = np.array(points[0]), np.array(points[-1])
            direction = end-start
            best_point = np.array([float(row['f_s']), float(row['f_b'])])
            tbest = float((best_point-start)@direction/(direction@direction)) if direction@direction else 0.
            def cost(t):
                fs, fb = start+t*direction
                yields = base+fs*top+fb*deep
                predicted = yields[1]/yields[1, -1]
                predicted[low] = yields[0, low]*yields[1, ip]/yields[0, ip]/yields[1, -1]
                residual = (np.log(predicted[idx])-np.log(observed))/errors
                logs = np.clip(-float(design@residual)/float(design@design), np.log(.3), np.log(3.))
                residual += logs*design
                return float(residual@residual)
            inv = inventory(benchmark.responses[material][i1], depth)
            for t in allowed_segment(cost, tbest, threshold):
                fs, fb = start+t*direction
                occupied = fs*inv['top_fuel_per_metal']+fb*(inv['fuel_per_metal']-inv['top_fuel_per_metal'])
                volumes.setdefault((material, shift), []).append(float(occupied))
    for row in curves:
        values = volumes.get((row['material'], float(row['U_high_keV'])), [])
        row['full_band_volume_occupied_min_per_metal'] = min(values) if values else None
        row['full_band_volume_occupied_max_per_metal'] = max(values) if values else None
    write_csv(out/'depth_inventory_curves.csv', curves)
    return dict(postprocessed_accepted_cells=True, volume_band_U_cells=len(volumes))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=HERE)
    parser.add_argument('--volume-bands-only', action='store_true',
                        help='Add volume-only extrema using saved accepted depth cells')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    benchmark = Benchmark()
    result = add_volume_bands(benchmark, args.out) if args.volume_bands_only else evaluate(benchmark, args.out)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
