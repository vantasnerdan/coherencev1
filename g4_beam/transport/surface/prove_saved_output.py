"""One narrow actual saved-output intake and scientific wrong-input refusal."""
import argparse
import copy
import json
from pathlib import Path
from surface_consumer import primitive, bound_json, rate_range, admit_original_transport, adsorption
from surface_io import manifest, save


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replacement', required=True)
    parser.add_argument('--primitive', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    replacement = bound_json(manifest(args.replacement))
    binding = manifest(args.primitive)
    value = primitive(binding)
    outcomes = {'primitive_binding': binding,
                'consumed_energy_extrapolated_eV': value['energy_extrapolated_eV'],
                'consumed_free_energy_eV': value['free_energy_eV'],
                'consumed_fmax_eV_A': value['fmax_eV_A'],
                'consumed_atom_count': len(value['geometry']['symbols']),
                'replacement_binding': manifest(args.replacement), 'refusals': []}
    if replacement['eligibility']['adsorption']:
        refs = replacement['primitive_bindings']
        actual = {name: primitive(refs[name]) for name in ('clean', 'adsorbed', 'H_atom', 'H2_gas')}
        outcomes['consumed_adsorption'] = adsorption(
            actual['clean'], actual['adsorbed'], actual['H_atom'], actual['H2_gas'],
            replacement['adsorbate_count'])
        wrong_atom = copy.deepcopy(actual['H_atom'])
        wrong_atom['spin_polarized'] = False
        try:
            adsorption(actual['clean'], actual['adsorbed'], wrong_atom, actual['H2_gas'],
                       replacement['adsorbate_count'])
        except ValueError as exc:
            outcomes['refusals'].append({'case': 'unpolarized_atomic_H_reservoir',
                                        'refused': True, 'reason': str(exc)})
        else:
            raise AssertionError('Wrong atomic reservoir accepted')
    if replacement['eligibility']['barrier']:
        outcomes['conditional_rate_range'] = rate_range(replacement, 300., 0.5)
    missing_barrier = copy.deepcopy(replacement)
    missing_barrier['eligibility']['barrier'] = False
    for name, operation in [
        ('incomplete_barrier_is_not_rate', lambda: rate_range(missing_barrier, 300., 0.5)),
        ('ineligible_original_transport_replacement', lambda: admit_original_transport(replacement, None))]:
        try:
            result = operation()
        except ValueError as exc:
            outcomes['refusals'].append({'case': name, 'refused': True, 'reason': str(exc)})
        else:
            outcomes['refusals'].append({'case': name, 'refused': False, 'actual_output': result})
    # Changing only the saved byte identity detects a wrong physical primitive.
    wrong = copy.deepcopy(binding)
    wrong['sha256'] = ('0' if wrong['sha256'][0] != '0' else '1')+wrong['sha256'][1:]
    try:
        primitive(wrong)
    except ValueError as exc:
        outcomes['refusals'].append({'case': 'wrong_saved_primitive_identity', 'refused': True, 'reason': str(exc)})
    else:
        raise AssertionError('Wrong input was accepted')
    if not all(x['refused'] for x in outcomes['refusals']):
        raise AssertionError('A scientific eligibility gate was bypassed')
    outcomes['claim'] = 'Actual saved primitive numerical intake and refusal; NOT completed surface-rate calibration'
    save(Path(args.output), outcomes)
    print(json.dumps(outcomes, indent=2))

if __name__ == '__main__':
    main()
