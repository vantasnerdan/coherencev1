"""Actual relaxed density/geometry -> two alternative conditional yield terms."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from g4_beam.merit.api import Target, compare, reference
from substrate_framework.cluster_screening import local_tf_shift_eV, spectator_shift_eV


def binding(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    return dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def completed_density(path):
    data = json.loads(Path(path).read_text())
    source = data['source_point']
    if binding(source['path']) != source:
        raise ValueError('The actual source primitive changed')
    primitive = json.loads(Path(source['path']).read_text())
    if (not primitive.get('converged_relaxation') or
            primitive['bank'] != data['source_bank'] or
            primitive['case'] != data['case'] or
            data['raw_fmax_eV_A'] > data['stationarity_target_eV_A'] or
            not data['midpoint_grid_target_met']):
        raise ValueError('A completed, stationary, grid-qualified actual pair is required')
    if data['closest_pair'] is None:
        raise ValueError('This saved target has no distinct D-D pair')
    fine = next(row for row in data['density_rows'] if row['grid_refinement'] == 2)
    return data, fine['PAW_valence']['midpoint_e_per_A3']


def consume(cluster_path, reference_path, selection_label):
    cluster, density = completed_density(cluster_path)
    ref, ref_density = completed_density(reference_path)
    distance = cluster['closest_pair']['periodic_distance_A']
    ref_distance = ref['closest_pair']['periodic_distance_A']
    shifts = {
        'frozen_positive_spectator': spectator_shift_eV(distance, ref_distance),
        'local_valence_TF_proxy': local_tf_shift_eV(density, ref_density),
    }
    base = reference('PdD')
    rows = []
    for closure, shift in shifts.items():
        target = Target(name=f"{cluster['case']['name']}:{closure}",
                        profiles=base.profiles, material=base.material,
                        stopping_fn=base.stopping_fn, columns=base.columns,
                        measurement=dict(base.measurement), Ue_shift_keV=shift/1000)
        for row in compare(target):
            row.update(closure=closure, additional_shift_eV=shift)
            rows.append(row)
    return dict(schema='coherence-pairs-conditional-yield-v1',
                selection_label=selection_label, case=cluster['case']['name'],
                inputs=[binding(cluster_path), binding(reference_path)],
                closest_D_D_A=distance, reference_closest_D_D_A=ref_distance,
                valence_midpoint_e_A3=density, reference_valence_midpoint_e_A3=ref_density,
                shifts_eV=shifts, rows=rows,
                baseline_commit='7ec4b53c1a9608d7d5f9ae1fd266a89625d93f8c',
                assumptions=[
                    'Every transported target D is assigned to one pair:100percent-pairing sensitivity, not measured pair fraction',
                    'Frozen positive spectator and local TF are alternatives, never summed',
                    'Changes are relative to isovolume defect-free rocksalt Pd8D8, not absolute screening added twice',
                    'Positive-spectator formula is the frozen leading coalescence limit, not a solved three-body trajectory',
                    'PAW-valence midpoint density is only a uniform-electron TF proxy, not a measured conduction susceptibility',
                    'Transported run-end profiles, target D count, stopping, adsorbate and fitted cause parameters are unchanged',
                    'All three charge readings and paper-weighted Eq7 fits are the existing compare(target) outputs',
                ], physical_pair_fraction_computed=False, nonlocal_dielectric_computed=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cluster-density', required=True)
    parser.add_argument('--reference-density', required=True)
    parser.add_argument('--selection-label', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    started = time.process_time()
    result = consume(args.cluster_density, args.reference_density, args.selection_label)
    result['process_cpu_s'] = time.process_time()-started
    Path(args.output).write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: result[key] for key in
                     ('case', 'selection_label', 'shifts_eV', 'process_cpu_s')}))
