"""Consume actual GPW banks and expose missed peaks before a path claim."""
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from ase.io import read
from scipy.interpolate import CubicHermiteSpline
from force_bank_results import verify_force_record
from surface_io import manifest


def bound_json(reference):
    raw = Path(reference['path']).read_bytes()
    if sha256(raw).hexdigest() != reference['sha256']:
        raise ValueError('Native reference byte identity changed')
    return json.loads(raw)


def reference_energy(clean, gas):
    clean_record, gas_record = bound_json(clean), bound_json(gas)
    multiplier = clean.get('multiplier', 1)
    if multiplier not in (1, 2):
        raise ValueError('Explicit one-cell or qualified two-cell folding required')
    value = multiplier*float(clean_record[clean['energy_field']])+float(gas_record['energy_extrapolated_eV'])
    if not np.isfinite(value):
        raise ValueError('Finite authenticated clean plus gas energy required')
    return value


def interval_peak(left, right):
    """Force-consistent Hermite diagnostic, not a new electronic calculation.

    Raw forces differentiate free energy. E-F is interpolated linearly, not
    silently treated as force-consistent extrapolated energy.
    """
    displacement = np.asarray(right['bank']['positions_A'])-np.asarray(left['bank']['positions_A'])
    slopes = [-float(np.sum(np.asarray(point['bank']['forces'])*displacement))
              for point in (left, right)]
    energies = [float(point['bank']['energy_eV']) for point in (left, right)]
    free = [float(point['free_energy_eV']) for point in (left, right)]
    entropy_slope = (energies[1]-free[1])-(energies[0]-free[0])
    curve = CubicHermiteSpline([0., 1.], energies, np.asarray(slopes)+entropy_slope)
    interior = [float(value) for value in curve.derivative().roots(extrapolate=False)
                if np.isfinite(value) and 0 < value < 1]
    fractions = [0., 1.]+interior
    values = [float(curve(value)) for value in fractions]
    peak = int(np.argmax(values))
    return dict(fraction=fractions[peak], predicted_energy_eV=values[peak],
                sampled_endpoint_max_eV=max(energies),
                predicted_unsampled_rise_eV=max(values)-max(energies),
                interpretation='Interpolation diagnostic only; no interpolated energy is admitted as DFT')


def _bound_reference_path(reference, base):
    path = base/Path(reference['path'])
    actual = manifest(path)
    if any(actual[key] != reference[key] for key in ('bytes', 'sha256')):
        raise ValueError('Saved Native reference byte identity changed')
    return path


def read_completed(directory, route_reference, reference_base):
    directory = Path(directory)
    base = Path(reference_base).resolve()
    route_path = _bound_reference_path(route_reference, base)
    frames = read(route_path, index=':')
    inputs = json.loads((directory/'input.json').read_bytes())
    if inputs['route'] != route_reference:
        raise ValueError('Native run binds a different route')
    terminal = json.loads((directory/'terminal.json').read_bytes())
    points = []
    for reference in terminal['completed_points']:
        point_path = _bound_reference_path(reference, base)
        point = json.loads(point_path.read_bytes())
        if point['stage'] != 'completed':
            raise ValueError('Unfinished Native frame is not an energy')
        bank = dict(point['bank'])
        for field in ('gpw', 'geometry'):
            bank[field] = {**bank[field], 'path': str(base/Path(bank[field]['path']))}
        atoms = verify_force_record(bank)
        point = {**point, 'bank': bank}
        index = point['route_index']
        if index not in inputs['selected_indices'] or not 0 <= index < len(frames):
            raise ValueError('Native point outside declared route selection')
        expected = frames[index]
        if (list(atoms.symbols) != list(expected.symbols)
                or not np.array_equal(atoms.positions, expected.positions)
                or not np.array_equal(atoms.cell, expected.cell)
                or not np.array_equal(atoms.pbc, expected.pbc)
                or not np.array_equal(atoms.get_masses(), expected.get_masses())):
            raise ValueError('Actual GPW geometry differs from selected route frame')
        points.append(point)
    indices = [point['route_index'] for point in points]
    if len(indices) != len(set(indices)):
        raise ValueError('Duplicate Native point cannot complete a route')
    return inputs, terminal, sorted(points, key=lambda point: point['route_index']), frames


def _clamp(frame):
    fixed = set()
    for constraint in frame.constraints:
        if constraint.todict()['name'] != 'FixAtoms':
            raise ValueError('Only the declared fixed-atom clamp is supported')
        fixed.update(int(index) for index in constraint.get_indices())
    return fixed


def compose_completed(segments, route_reference):
    """Select genuine records by exact target pose, preserving source runs."""
    if not segments or isinstance(segments, (str, Path)):
        raise ValueError('Explicit Native source-route selections required')
    if manifest(Path(route_reference['path'])) != route_reference:
        raise ValueError('Target route byte identity changed')
    frames = read(route_reference['path'], index=':')
    inputs, terminals, points, selections = [], [], [], []
    for segment in segments:
        source, terminal, available, source_frames = read_completed(
            segment['directory'], segment['route'], segment['reference_base'])
        requested = segment['indices']
        if not requested or len(requested) != len(set(requested)):
            raise ValueError('A nonempty disjoint source selection is required')
        by_index = {point['route_index']: point for point in available}
        if any(index not in by_index or index >= len(frames) for index in requested):
            raise ValueError('Selected actual source point is unavailable')
        for index in requested:
            left, right = source_frames[index], frames[index]
            arrays = ((left.numbers, right.numbers), (left.positions, right.positions),
                      (left.cell, right.cell), (left.pbc, right.pbc),
                      (left.get_masses(), right.get_masses()))
            if any(not np.array_equal(a, b) for a, b in arrays) or _clamp(left) != _clamp(right):
                raise ValueError('Actual source pose differs from the target route')
            points.append(by_index[index])
        inputs.append(source)
        terminals.append(terminal)
        selections.append(dict(directory=str(segment['directory']),
                               route=segment['route'], indices=list(requested),
                               reference_base=str(segment['reference_base'])))
    settings = ('kmesh', 'cutoff_eV', 'xc', 'Pd_valence', 'H_valence', 'convergence', 'symmetry')
    if any(any(item[key] != inputs[0][key] for key in settings) for item in inputs[1:]):
        raise ValueError('Native path segments use different electronic settings')
    points.sort(key=lambda point: point['route_index'])
    indices = [point['route_index'] for point in points]
    if len(indices) != len(set(indices)):
        raise ValueError('Overlapping Native selections cannot complete a route')
    return inputs, terminals, points, frames, selections


def analyze_path(segments, route_reference, clean_reference, gas_reference,
                 *, peak_target_eV=.03, refinement_target_eV=.01):
    inputs, terminals, points, frames, selections = compose_completed(segments, route_reference)
    baseline = reference_energy(clean_reference, gas_reference)
    indices = [point['route_index'] for point in points]
    complete = indices == list(range(len(frames))) and len(frames) >= 3
    result = dict(schema='coherence-connected-elbow-readout-v1',
        Native_inputs=inputs, Native_terminals=terminals, Native_selections=selections,
        references=dict(route=route_reference, clean=clean_reference, gas=gas_reference),
        clean_plus_gas_energy_eV=baseline, complete_route=complete,
        peak_target_eV=peak_target_eV, refinement_target_eV=refinement_target_eV,
        curve=[dict(route_index=point['route_index'], energy_eV=point['bank']['energy_eV'],
            E_minus_gas_eV=point['bank']['energy_eV']-baseline,
            separation_A=point['separation_A'], height_A=point['height_A'],
            actual_bank=point['bank']['gpw']) for point in points],
        capture_probability_computed=False)
    if not complete:
        result['status'] = 'incomplete_route_no_classification'
        return result
    if (len(frames) % 2 != 1
            or any(frame.info.get('nested_coarse_point') != (index % 2 == 0)
                   for index, frame in enumerate(frames))):
        raise ValueError('Declared nested midpoint route required for refinement classification')
    peaks = [interval_peak(left, right) for left, right in zip(points[:-1], points[1:])]
    sampled = max(point['bank']['energy_eV'] for point in points)
    coarse = max(point['bank']['energy_eV'] for point in points[::2])
    interpolation = max(point['predicted_energy_eV'] for point in peaks)
    change = sampled-coarse
    missed = interpolation-sampled
    asymptote = abs(points[0]['bank']['energy_eV']-points[2]['bank']['energy_eV'])
    result.update(sampled_peak_excess_eV=sampled-baseline,
        coarse_to_refined_peak_change_eV=change,
        force_interpolated_peak_excess_eV=interpolation-baseline,
        maximum_predicted_unsampled_rise_eV=missed,
        far_height_change_eV=asymptote, interval_diagnostics=peaks,
        sampled_peak_index=indices[int(np.argmax([p['bank']['energy_eV'] for p in points]))])
    passed = (change <= refinement_target_eV and missed <= refinement_target_eV
              and max(sampled, interpolation)-baseline <= peak_target_eV and asymptote <= .01)
    result['status'] = ('nonactivated_route_candidate_pending_k_slab_controls' if passed
                        else 'route_does_not_yet_decide_activation')
    return result
