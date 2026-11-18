"""Authenticate every admitted final saved bank without starting a calculator."""
import argparse
import json
from pathlib import Path

from g4_beam.pairs.force_bank_results import verify_force_record
from g4_beam.pairs.pair_binding_results import load_completed
from g4_beam.pairs.surface_io import manifest, save


def verify():
    report = json.loads(Path('scientific-result.json').read_text())
    rows = []
    for family in ('reference', 'axis', 'tetra'):
        _, _, completed = load_completed(Path('data') / family)
        rows.extend(completed)
    expected = report['material_cells']
    if len(rows) != 12 or len(expected) != 12:
        raise ValueError('The complete twelve-cell final bank set is required')
    checks = []
    for cell in expected:
        primitive = cell['source_primitive']
        if manifest(primitive['path']) != primitive:
            raise ValueError('Portable primitive bytes changed')
        point = json.loads(Path(primitive['path']).read_text())
        atoms = verify_force_record(point['bank'])
        if (point['bank']['gpw'] != cell['gpw'] or
                point['bank']['geometry'] != cell['geometry'] or
                point['bank']['energy_eV'] != cell['energy_eV'] or
                point['step'] != cell['step'] or
                point['case']['name'] != cell['case'] or
                point['neutral_electrons'] != cell['neutral_electrons'] or
                point['electronic_free_energy_eV'] != cell['free_energy_eV']):
            raise ValueError('Final bank differs from the frozen report')
        checks.append(dict(case=cell['case'], family=cell['family'],
                           atoms=len(atoms), gpw=cell['gpw'], geometry=cell['geometry'],
                           original_report_energy_eV=cell['energy_eV']))
    for model in report['canonical_model_sources']:
        if manifest(model['path']) != model:
            raise ValueError('Exact-byte canonical model snapshot changed')
    return dict(status='passed', final_banks=12, canonical_models=4,
                raw_force_energy_geometry_isotope_authenticated=True,
                no_calculator_or_scf_started=True, checks=checks)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = verify()
    save(args.output, result)
    print(json.dumps(dict(status=result['status'], final_banks=result['final_banks'],
                          canonical_models=result['canonical_models'],
                          no_calculator_or_scf_started=True)))
