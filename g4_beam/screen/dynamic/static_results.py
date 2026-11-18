"""Evaluate complete registered microscopic static quadratures; no Native work."""
import argparse
import csv
import json
from pathlib import Path

from g4_beam.screen.dynamic.largeq_bank import binding
from g4_beam.screen.dynamic.static_quadrature import evaluate_registry


def write_csv(path, result):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            'model', 'material', 'pd_valence', 'ground_PW_eV', 'kmesh',
            'response_G_eV', 'response_bands', 'q_midpoint_n', 'r_A',
            'U_eV', 'reference_parameter_sensitivity_eV'])
        writer.writeheader()
        for case in result['cases']:
            for quadrature in case['quadratures']:
                for radius, value, sensitivity in zip(
                        result['radii_A'], quadrature['primary_U_eV'],
                        quadrature['reference_parameter_max_difference_eV']):
                    writer.writerow(dict(model=case['id'], material=case['material'],
                        pd_valence=case['pd_valence'],
                        ground_PW_eV=case['ground_state']['ecut_eV'],
                        kmesh='x'.join(map(str, case['ground_state']['kmesh'])),
                        response_G_eV=case['response']['ecut_eV'],
                        response_bands=case['response']['nbands'],
                        q_midpoint_n=quadrature['n'], r_A=radius, U_eV=value,
                        reference_parameter_sensitivity_eV=sensitivity))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registry', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--csv', required=True, type=Path)
    args = parser.parse_args()
    registry = json.loads(args.registry.read_text())
    result = evaluate_registry(registry, args.registry.parent)
    result['registry_binding'] = binding(args.registry)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    write_csv(args.csv, result)
    print(json.dumps(dict(output=binding(args.output), csv=binding(args.csv),
        cases=[dict(id=case['id'], contact_eV=[dict(n=q['n'], U=q['primary_U_eV'][0],
                     reference_sensitivity=q['reference_parameter_max_difference_eV'][0])
                     for q in case['quadratures']]) for case in result['cases']])))


if __name__ == '__main__':
    main()
