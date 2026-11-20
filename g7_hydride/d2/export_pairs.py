"""Export hash-bound total BO pair samples for the pinned g6_pairlaw consumer.

Only NumPy and the standard library are needed. A saved readout is admitted only
when a fresh physical replay reproduces its complete structure and provenance.
"""
import argparse
import json
import math
from pathlib import Path

from hydride_replay import (BOHR_A, RADII, binding, compare_json,
                            number, read_json, require)
from readout import summarize

SCHEMA = 'coherencev1-target-conditioned-pair-v1'
OBSERVABLE = 'target_conditioned_BO_total_interaction'
CONSUMER = 'g6_pairlaw/script/material_pair.py@2645b19a'


def check_pair(curve):
    require(curve['schema'] == SCHEMA and curve['observable'] == OBSERVABLE,
            'Require pinned target-conditioned total BO pair schema')
    require('prepared_state' not in curve, 'BO export must not invent a prepared state')
    require(curve['R_A'] == RADII and len(curve['V_pair_total_eV']) == len(RADII),
            'Require exact sorted six-R support')
    for radius, value in zip(curve['R_A'], curve['V_pair_total_eV']):
        require(number(radius, 'Radius') > 0, 'Positive radius required')
        number(value, 'Total BO interaction')
    reference = curve['reference']
    for key in ('source', 'configuration', 'energy_zero', 'periodic_to_isolated'):
        require(isinstance(reference[key], str) and reference[key].strip(), 'Missing reference ' + key)
    require(reference['bare'] == 'isolated_unit_charge_coulomb', 'Wrong bare Coulomb reference')
    json.dumps(curve, allow_nan=False)


def material_pair(label, readout, readout_binding):
    case, branch = label.split('_')
    charge = int(branch[1])
    rows = readout['rows'][label]
    source = readout['sources'][case]
    provenance = readout['provenance'][case]
    charge_data = provenance['charges'][branch]
    curve = dict(
        schema=SCHEMA, observable=OBSERVABLE,
        R_A=[row['R_A'] for row in rows],
        V_pair_total_eV=[row['V_pair_total_eV'] for row in rows],
        reference=dict(
            bare='isolated_unit_charge_coulomb',
            source='P255/issue248 same-code actual Native ' + label + '; quartet SHA256 ' + source['quartet']['sha256'],
            configuration=(
                ('Force-relaxed rocksalt LiD finite FD cell at inherited original PW volume. ' if case == 'F' else
                 'D atom in a periodic 10 A cubic box. ' if case == 'G' else
                 'D2 molecular target in a periodic 10 A cubic box; approach outside target bond. ') +
                'Neutral occupied and target-removed hosts; one target nucleus and one electron removed. ' +
                ('Neutral probe adds one electron.' if charge == 0 else 'Bare unit-charge probe adds no electron.')),
            energy_zero=(
                'E_probe,occupied - E_probe,vacant - E_host,occupied + E_host,vacant, each actual '
                'PBE/Hae nonspin Native internal energy = free energy - negative TS, FD width 0.1 eV. '
                'Only one-atom G vacant-host is the analytic zero-nucleus/zero-electron empty Hamiltonian; '
                'all one-H zero-electron q1 probe states are actual converged Native. '
                'No fitted offset or far-point subtraction; no geometry PAW energy enters the quartet.'),
            periodic_to_isolated=(
                'Replace the Hae alpha1=10 Bohr^-2 Gaussian nuclear block with the independently '
                'replayed periodic same-zero-G point-minus-Gaussian block, width sqrt(5)/Bohr. '
                'Both F charges pass actual h0.20/0.16 and relaxed repeat2/repeat3 comparisons at R0.1 '
                'within max(5 eV,10% of control lowering). No bare-only image subtraction and no exact '
                'isolated electronic limit is claimed; gas remains a finite-box method reference.')),
        source_bindings=dict(readout=readout_binding, native=source, all_controls=readout['sources']),
        configuration=dict(case=case, isotope='D', probe_charge=charge,
                           probe_added_electrons=1 - charge,
                           occupied_host_electrons=provenance['occupied_host_electrons'],
                           vacant_host_electrons=provenance['vacant_host_electrons'],
                           occupied_probe_electrons=charge_data['occupied_probe_electrons'],
                           vacant_probe_electrons=charge_data['vacant_probe_electrons'],
                           specification=provenance['specification'],
                           volume_A3=provenance['volume_A3'], volume_scope=provenance['volume_scope'],
                           path=provenance['geometry']['path'],
                           realspace_shape=provenance['realspace_shape'],
                           actual_kpoints=provenance['actual_kpoints'], nbands=provenance['nbands'],
                           BO_setups=provenance['BO_setups']),
        geometry_provenance=provenance['geometry'],
        original_volume_provenance=provenance['original_volume_provenance'],
        charge_gauge_scope=readout['charge_gauge_scope'],
        mapping_provenance=dict(rows=rows, Hae_alpha1_Bohr_minus2=10,
                                gaussian_width_inverse_A=math.sqrt(5.) / BOHR_A,
                                zero_G='same zero-G periodic point and Gaussian convention',
                                lowering_definition='14.399645468667815/R_A - V_pair_total_eV'),
        qualification=dict(qualified=True, both_F_charge_h_and_cell_controls=readout['F_controls'],
                           ratios_R01=readout['ratios_R01'], PdD_reference=readout['PdD_reference']),
        consumer_contract=CONSUMER,
        support=readout['scope'])
    check_pair(curve)
    return curve


def export_curves(readout_path, f, g, fine, large, g2=None):
    """Return filename->JSON only after complete independent Native-source replay."""
    readout_path = Path(readout_path)
    original = binding(readout_path)
    saved = read_json(readout_path)
    replayed = summarize(f, g, fine, large, g2)
    compare_json(saved, replayed, 'Complete hash-bound hydride readout')
    curves = {label + '_BO_pair.json': material_pair(label, replayed, original)
              for label in replayed['rows']}
    expected = {case + '_' + charge + '_BO_pair.json'
                for case in (('F', 'G') if g2 is None else ('F', 'G', 'G2')) for charge in ('q0', 'q1')}
    require(set(curves) == expected, 'Wrong exported case/charge set')
    require(binding(readout_path) == original, 'Readout changed during export')
    return curves


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('readout', 'f', 'g', 'fine', 'large', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--g2', type=Path)
    args = parser.parse_args()
    require(not args.output.exists(), 'Output directory must be new (no stale/overwritten curves)')
    curves = export_curves(args.readout, args.f, args.g, args.fine, args.large, args.g2)
    encoded = {name: json.dumps(curve, indent=2, sort_keys=True, allow_nan=False) + '\n'
               for name, curve in curves.items()}
    args.output.mkdir()
    for name, content in sorted(encoded.items()):
        with (args.output / name).open('x', encoding='utf-8') as stream:
            stream.write(content)
    print(json.dumps(dict(output=args.output.name, files=sorted(curves), qualified=True), allow_nan=False))


if __name__ == '__main__':
    main()
