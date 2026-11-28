"""Report matched-temperature changes in identifiable PdD mode combinations."""
import argparse
import hashlib
import json
from pathlib import Path

from native import require, write_once

MODES = ('gamma_100', 'x_100', 'x_010')
SPLIT_NAMES = ('K_diagonal_folded_eV_A2', 'J_cross_folded_eV_A2', 'J_shell1_plus_shell3_scalar_projection_eV_A2')


def diff(row, base):
    require(row['primitive_mesh_n'] == base['primitive_mesh_n'] and base['temperature_K'] == 300,
            'Nonmatching temperature reference')
    modes = {name: row[name]['K_force_eV_A2'] - base[name]['K_force_eV_A2'] for name in MODES}
    alternate = {name: row[name]['K_energy_eV_A2'] - base[name]['K_energy_eV_A2'] for name in MODES}
    out = dict(electronic_excitation_eV_per_D_site=row['electronic_excitation_eV_per_D_site'],
               mode_K_change_eV_A2=modes, energy_stencil_mode_K_change_eV_A2=alternate,
               all_Gamma_agree_with_primitive=all(q['gamma_matches_primitive'] for q in (row,base)))
    if out['all_Gamma_agree_with_primitive']:
        out['split_change'] = {pol: {key: row[f'folded_split_{pol}'][key] - base[f'folded_split_{pol}'][key]
                                    for key in SPLIT_NAMES} for pol in ('100', '010')}
    return out


def build(response):
    require(response['schema'] == 'g16-pdd-folded-mode-v1', 'Not the G16 folded response')
    raw = {(x['temperature_K'],x['primitive_mesh_n']):x for x in response['rows']}
    require(len(response['rows']) == len(raw) == 6 and set(raw) == {(t,k) for t in (300,3000,10000) for k in (12,16)},
            'Need all two-mesh folded-mode rows')
    result = []
    for t in (300,3000,10000):
        a,b = (diff(raw[t,k],raw[300,k]) for k in (12,16))
        compare = {name: abs(a['mode_K_change_eV_A2'][name] - b['mode_K_change_eV_A2'][name]) for name in MODES}
        row = dict(temperature_K=t, mesh12=a, mesh16=b, delta_K_mesh_difference_eV_A2=compare)
        if 'split_change' in a and 'split_change' in b:
            row['split_mesh_difference'] = {pol: {key: abs(a['split_change'][pol][key] - b['split_change'][pol][key])
                                                   for key in SPLIT_NAMES} for pol in ('100','010')}
        result.append(row)
    return dict(schema='g16-pdd-folded-thermal-contrast-v1', electronic_temperature_reference_K=300,
                sign_and_units='F=-Kii*u+sum(J*u_neighbour), projected per D site', rows=result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('folded_response', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    data = a.folded_response.read_bytes()
    result = build(json.loads(data))
    result['folded_response_sha256'] = hashlib.sha256(data).hexdigest()
    write_once(a.output, result)


if __name__ == '__main__':
    main()
