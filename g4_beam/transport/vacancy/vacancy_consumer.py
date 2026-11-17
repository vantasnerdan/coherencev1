"""Consume authentic completed primitive JSON; never promote an untested trap capacity."""
import argparse
import json
import math
from pathlib import Path
from vacancy_saved import load_saved, require


def load_group(root, references=()):
    records, geometry_audits, bindings = load_saved([root, *references])
    common = records['bulk']['config']
    e = {case: record['energy_eV'] for case, record in records.items()}
    free_e = {case: record['electronic_free_energy_eV'] for case, record in records.items()}
    insertion = e['bulk_D'] - e['bulk']
    n_values = sorted(int(case.split('_D')[1]) for case in records if case.startswith('vacancy_D'))
    require(n_values == list(range(1, max(n_values) + 1)), 'Missing incremental occupancy')
    increments = []
    free_increments = []
    relative = {0: 0.0}
    for n in n_values:
        previous = e['vacancy'] if n == 1 else e[f'vacancy_D{n-1}']
        increments.append(dict(n=n, binding_eV_per_D=previous + insertion - e[f'vacancy_D{n}']))
        free_previous = free_e['vacancy'] if n == 1 else free_e[f'vacancy_D{n-1}']
        free_increments.append(dict(n=n, binding_eV_per_D=free_previous
                                    + free_e['bulk_D'] - free_e['bulk'] - free_e[f'vacancy_D{n}']))
        relative[n] = e[f'vacancy_D{n}'] - e['vacancy'] - n * insertion
    # Exact pairwise crossings define the lower envelope among actually tested n.
    crossings = sorted({(relative[n] - relative[m]) / (n-m)
                        for n in relative for m in relative if n > m})
    boundaries = [-math.inf] + crossings + [math.inf]
    intervals = []
    for low, high in zip(boundaries[:-1], boundaries[1:]):
        sample = (low + high) / 2 if math.isfinite(low) and math.isfinite(high) else (
            high - 1 if not math.isfinite(low) else low + 1)
        preferred = min(relative, key=lambda n: relative[n] - n * sample)
        intervals.append(dict(delta_mu_lower_eV=None if not math.isfinite(low) else low,
                              delta_mu_upper_eV=None if not math.isfinite(high) else high,
                              preferred_n_among_tested=preferred))
    return dict(primitive_energies_eV=e, bulk_D_insertion_eV=insertion,
                vacancy_formation_eV=e['vacancy']-(common['repeat']**3-1)*e['bulk']/common['repeat']**3,
                incremental_trapping=increments, relative_vacancy_energies_eV=relative,
                finite_electronic_smearing_binding=free_increments,
                electronic_free_energies_eV=free_e,
                zero_temperature_reservoir_intervals=intervals,
                maximum_n_tested=max(n_values), capacity_proven=False,
                relaxed_geometry_audit=geometry_audits,
                config=common, source_bindings=bindings,
                reference='delta_mu=mu_D-[E(bulk+D)-E(bulk)]; positive B favors transfer from bulk interstitial to vacancy',
                units=dict(binding='eV/D', occupancy='D per missing Pd atom'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('root')
    parser.add_argument('--output', required=True)
    parser.add_argument('--reference', action='append', default=[])
    args = parser.parse_args()
    result = load_group(args.root, args.reference)
    result.update(eligible=False,
                  uncertainty_eV=None,
                  original_Eb_replacement_eligible=False, original_mtrap_replacement_eligible=False,
                  reasons=['No dominant image/k/cutoff sensitivity supplied',
                           'No ionic ZPE/vibrational free energies; H/D share electronic surface only',
                           'No finite-temperature mapping of bulk electronic insertion energy to original isotherm mu0',
                           'Maximum occupancy tested is not proven capacity; no unique reservoir-independent n'])
    Path(args.output).write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({key: result[key] for key in ('incremental_trapping', 'eligible', 'maximum_n_tested')}))


if __name__ == '__main__':
    main()
