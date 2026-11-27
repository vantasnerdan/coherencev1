"""Compare PdD response at the same track excitation on both electronic meshes."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from native import write_once

TRACK_SHA256 = '73fa1794047daf82b3dc25247e280b281dba8c6a6421c74ce87694926df18bed'
TRACK_BYTES = 299176
LOW_CURVATURE_RESOLUTION_EV_A2 = 2e-4


def track_levels(path, volume_A3):
    raw = path.read_bytes()
    if len(raw) != TRACK_BYTES or hashlib.sha256(raw).hexdigest() != TRACK_SHA256:
        raise ValueError('Track density does not match the original G15 pass-1 source')
    samples = {}
    for row in csv.DictReader(raw.decode().splitlines()):
        channel = row['channel']
        if channel not in ('thermal_electrons', 'nonthermal_electrons'):
            continue
        key = tuple(row[name] for name in ('scenario', 'species', 'time_s', 'radius_m'))
        density = float(row['excess_eV_m3'])
        group = samples.setdefault(key, {})
        if channel in group or not math.isfinite(density) or density < 0:
            raise ValueError('Conflicting track density samples')
        group[channel] = density
    result = {}
    for label, combined in [('thermal_only', False), ('thermal_plus_nonthermal', True)]:
        energy = (lambda g: sum(g.values())) if combined else (lambda g: g['thermal_electrons'])
        key, values = max(samples.items(), key=lambda pair: energy(pair[1]))
        density = energy(values)
        result[label] = dict(excitation_eV_per_D_site=density * volume_A3 * 1e-30, excess_eV_m3=density,
                             peak_sample=dict(zip(('scenario', 'species', 'time_s', 'radius_m'), key)))
    return result


def bracket(rows, mesh, target):
    values = [r[f'mesh{mesh}'] for r in rows]
    exc = [v['contrast']['electronic_excitation_eV_per_D_site'] for v in values]
    if not all(math.isfinite(u) for u in exc) or exc[0] != 0:
        raise ValueError('Invalid excitation origin or value')
    for i, (lo, hi) in enumerate(zip(exc, exc[1:])):
        if not lo < hi:
            return dict(status='not_monotonic_up_to_excitation')
        if target > hi:
            continue
        parts = values[i:i+2]
        weight = (target - lo) / (hi - lo)
        interp = lambda key: (1-weight) * parts[0]['contrast'][key] + weight * parts[1]['contrast'][key]
        k = interp('K_change_eV_A2')
        fd = abs(k - interp('K_free_energy_change_eV_A2'))
        return dict(status='bracketed', temperature_K=[r['temperature_K'] for r in rows[i:i+2]],
                    electronic_excitation_bracket_eV_per_D_site=[lo, hi], interpolated_delta_K_eV_A2=k,
                    interpolated_fractional_K_change=interp('fractional_K_change'),
                    interpolated_force_energy_difference_eV_A2=fd)
    return dict(status='outside_sampled_excitation_range')


def join(readout, path):
    if readout['schema'] != 'g16-pdd-thermal-material-v1' or readout['reference_temperature_K'] != 300:
        raise ValueError('Expected a 300 K-referenced PdD response')
    rows = readout['rows']
    volume = rows[0]['mesh12']['zero']['volume_A3']
    result = track_levels(path, volume)
    for level in result.values():
        x = level['excitation_eV_per_D_site']
        a, b = (bracket(rows, k, x) for k in (12, 16))
        level.update(mesh12=a, mesh16=b)
        if a['status'] == b['status'] == 'bracketed':
            ka, kb = (v['interpolated_delta_K_eV_A2'] for v in (a,b))
            spread = abs(ka-kb)
            uncertainty = max(spread, LOW_CURVATURE_RESOLUTION_EV_A2,
                              a['interpolated_force_energy_difference_eV_A2'], b['interpolated_force_energy_difference_eV_A2'])
            level.update(mesh_spread_eV_A2=spread, numerical_comparison_floor_eV_A2=uncertainty,
                         sign_above_floor_on_both_meshes=min(abs(ka),abs(kb))>uncertainty and ka*kb>0)
    return dict(schema='g16-pdd-track-join-v1', track_source=dict(bytes=TRACK_BYTES, sha256=TRACK_SHA256),
                electron_energy_reference_K=300, volume_A3_per_D_site=volume, modelled_peaks=result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('readout', type=Path)
    p.add_argument('track_density', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    data = a.readout.read_bytes()
    result = join(json.loads(data), a.track_density)
    result['readout_sha256'] = hashlib.sha256(data).hexdigest()
    write_once(a.output, result)


if __name__ == '__main__':
    main()
