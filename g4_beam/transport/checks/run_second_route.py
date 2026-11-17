"""Source-bound checks of the actual baseline quantities; no baseline replay."""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'screen'))
from screen_io import manifest, write_json
from lss_range import projected_range
from profile_checks import check_profiles
from thermo_checks import dilute_lag, dilute_sieverts, film_lifetimes


def ranges(cache, old_rust):
    accepted = {}
    for line in old_rust.read_text().splitlines()[1:]:
        fields = line.split()
        if len(fields) > 4 and float(fields[2]) == 1.53:
            accepted[(fields[0], int(fields[1]))] = float(
                next(v for v in fields if v.startswith('Rp=')).split('=')[1])
    rows = []
    for target, Z, mass, density, scale in [('Pd', 46, 106.42, 68., 1.54),
                                           ('Ti', 22, 47.87, 56.7, 1.52)]:
        component = [(Z, mass, density, scale)]
        coarse = projected_range(component, energy_nodes=256, impact_order=64, radial_order=64)
        fine = projected_range(component, energy_nodes=512, impact_order=96, radial_order=96)
        for E in [500, 13000]:
            projected = float(np.interp(E, fine['energy_eV'], fine['projected_nm']))
            before = float(np.interp(E, coarse['energy_eV'], coarse['projected_nm']))
            row = dict(target=target, E_lab_eV=E,
                       analytic_infinite_medium_projected_mean_nm=projected,
                       analytic_mean_path_nm=float(np.interp(E, fine['energy_eV'], fine['path_nm'])),
                       quadrature_relative_change=abs(projected-before)/projected,
                       quadrature_pass_5percent=bool(abs(projected-before)/projected < .05),
                       cached_half_space_implanted_mean_nm=float(cache[f'{target}_{E}_meta'][2]),
                       earlier_RustBCA_KrC_scaledLS_implanted_mean_nm=accepted[(target, E)],
                       LS_scale=scale, density_nm3=density)
            if target == 'Ti':
                row['cached_TiD2_implanted_mean_nm'] = float(cache[f'TiD2_{E}_meta'][2])
            rows.append(row)
    return dict(rows=rows, method='LSS1963 Eq4.1 with ZBL scattering and continuous LS electronic stopping; independently discretised Volterra equation',
                boundary_warning='LSS p24 explicitly defines an infinite homogeneous-medium projected moment. Cached/RustBCA means are conditional on implantation in a half-space. Different observables: do not overwrite the cache or claim a fresh4% agreement.',
                earlier_comparison='Earlier accepted RustBCA comparison is retained as source data, not rerun or newly certified')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--transport-module', type=Path, required=True)
    parser.add_argument('--profiles', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--earlier-rust-ranges', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.monotonic(); cpu = time.process_time_ns()
    sys.path.insert(0, str(args.transport_module.resolve().parent))
    spec = importlib.util.spec_from_file_location('transport_reference', args.transport_module)
    model = importlib.util.module_from_spec(spec); spec.loader.exec_module(model)
    profiles = json.loads(args.profiles.read_bytes())
    if len(profiles['records']) != 16:
        raise ValueError('Expected the sixteen original full profiles')
    with np.load(args.cache, allow_pickle=False) as cache:
        range_result = ranges(cache, args.earlier_rust_ranges)
    sieverts = dilute_sieverts(model.PD, 300., model.P_GAS)
    direct = model.Isotherm(model.PD, 300.).x_of_mu(.5*model.KB*300*np.log(model.P_GAS))
    sieverts['baseline_numerical_inverse_D_per_Pd'] = direct
    sieverts['relative_difference'] = abs(direct/sieverts['D_per_Pd']-1)
    result = dict(schema='coherence-transport-second-route-v1',
                  sources=dict(transport=manifest(args.transport_module, 'g4_beam/transport/transport.py'),
                               profiles=manifest(args.profiles, 'baseline/profiles_full.json'),
                               cache=manifest(args.cache), earlier_rust=manifest(args.earlier_rust_ranges),
                               code=[manifest(Path(__file__).parent/name) for name in
                                     ['lss_range.py', 'profile_checks.py', 'thermo_checks.py', 'run_second_route.py']]),
                  ranges=range_result, dilute_Sieverts=sieverts,
                  Pd_dilute_backside_lag=dilute_lag(model.PD, model.L_FOIL, 300.),
                  overlayer=film_lifetimes(model),
                  profiles=check_profiles(profiles['records'], model.stopping, model.sigma),
                  process_cpu_ns=str(time.process_time_ns()-cpu),
                  elapsed_wall_s=time.monotonic()-start,
                  scope='Independent numerical routes under the explicitly supplied baseline physics; not experimental validation of all estimated parameters')
    write_json(args.output, result)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
