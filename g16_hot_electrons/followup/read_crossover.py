"""Read the primitive crossover against the integrated item-one endpoints."""
import argparse
import hashlib
import json
from pathlib import Path

from native import write_once
from readout import NativeTable
from thermal_response import mesh_check, thermal_delta

T = (1000, 1500, 2000, 3000)
ITEM1_SHA = 'f19aaffe9b34b331b2031a8b62713a14e290e2a155d5fc9a04bf567c431e6383'


def status(mesh12, mesh16, diagnostics):
    a = mesh12['contrast']['K_change_eV_A2']
    b = mesh16['contrast']['K_change_eV_A2']
    floor = max(0.0002, diagnostics['delta_K_mesh_sensitivity_eV_A2'],
                diagnostics['delta_K_force_energy_difference_eV_A2'])
    if a > floor and b > floor:
        return 'stiffer_than_300K'
    if a < -floor and b < -floor:
        return 'softer_than_300K'
    return 'not_resolved'


def sign_brackets(rows, label):
    result = []
    prior = None
    for row in rows:
        current = row[label]
        if current == 'not_resolved':
            continue
        if prior is not None and current != prior[label]:
            result.append([prior['temperature_K'], row['temperature_K']])
        prior = row
    return result


def raw_sign(a, b):
    values = [v['contrast']['K_change_eV_A2'] for v in (a,b)]
    if min(values) > 0:
        return 'stiffer_than_300K'
    if max(values) < 0:
        return 'softer_than_300K'
    return 'not_resolved'


def build(root, producer, original):
    raw = original.read_bytes()
    if hashlib.sha256(raw).hexdigest() != ITEM1_SHA:
        raise ValueError('Integrated item-one response does not match')
    accepted = {r['temperature_K']: r for r in json.loads(raw)['rows']}
    reader = NativeTable(root, producer)
    result = []
    for t in T:
        row = accepted.get(t)
        if row is None:
            values = {}
            for k in (12, 16):
                item = reader.one_temperature(k, t, True)
                item['contrast'] = thermal_delta(item, accepted[300][f'mesh{k}'])
                values[f'mesh{k}'] = item
            row = dict(temperature_K=t, mesh_check=mesh_check(*values.values()), **values)
        result.append(dict(temperature_K=t, mesh12=row['mesh12'], mesh16=row['mesh16'], mesh_check=row['mesh_check'],
                           sign_status=status(row['mesh12'], row['mesh16'], row['mesh_check']),
                           raw_both_mesh_sign=raw_sign(row['mesh12'], row['mesh16']),
                           origin='accepted_item_1' if t in accepted else 'new_native'))
    brackets = sign_brackets(result, 'sign_status')
    return dict(schema='g16-pdd-crossover-v1', comparison_temperature_K=300, item1_response_sha256=ITEM1_SHA,
                native_sources=reader.sources, rows=result, sampled_sign_change_brackets_K=brackets,
                raw_same_sign_brackets_K=sign_brackets(result, 'raw_both_mesh_sign'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native-root', type=Path, required=True)
    p.add_argument('--item1', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    write_once(a.output, build(a.native_root, Path(__file__).with_name('native.py'), a.item1))


if __name__ == '__main__':
    main()
