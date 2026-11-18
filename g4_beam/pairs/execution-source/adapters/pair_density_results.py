"""Read actual compact GPW densities; never infer density from pair distance."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from gpaw import GPAW
from gpaw.mpi import serial_comm
from substrate_framework.cluster_observables import sample_periodic_density
from force_bank_results import verify_force_record
from surface_io import manifest, save


def pair_geometry(atoms):
    indices = [i for i, symbol in enumerate(atoms.symbols) if symbol == 'H']
    pairs = []
    for offset, i in enumerate(indices):
        for j in indices[offset + 1:]:
            vector = atoms.get_distance(i, j, mic=True, vector=True)
            direct = atoms.positions[j] - atoms.positions[i]
            midpoint = atoms.positions[i] + 0.5 * vector
            pairs.append(dict(atom_indices=[i, j], periodic_distance_A=float(np.linalg.norm(vector)),
                              same_cell_distance_A=float(np.linalg.norm(direct)),
                              periodic_displacement_A=vector.tolist(),
                              uses_periodic_image=bool(np.linalg.norm(vector-direct) > 1e-8),
                              midpoint_A=midpoint.tolist()))
    return sorted(pairs, key=lambda pair: pair['periodic_distance_A'])


def read_point(path):
    path = Path(path).resolve()
    point = json.loads(path.read_text())
    execution_path = path.parent.parent / 'execution.json'
    execution = json.loads(execution_path.read_text())
    for bound in execution['sources']:
        if manifest(bound['path']) != bound:
            raise ValueError('Original execution source/input bytes changed')
    config = json.loads(Path(execution['sources'][0]['path']).read_text())
    case = next(case for case in config['cases'] if case['name'] == path.parent.name)
    atoms = verify_force_record(point['bank'])
    if (not all(atoms.pbc) or atoms.constraints or
            list(atoms.symbols).count('Pd') != case['n_Pd'] or
            list(atoms.symbols).count('H') != case['n_D']):
        raise ValueError('Actual saved bank does not match the requested cluster target')
    pairs = pair_geometry(atoms)
    target = pairs[0]['midpoint_A'] if pairs else None
    calc = GPAW(point['bank']['gpw']['path'], legacy_gpaw=False,
                communicator=serial_comm, txt=None)
    expected_valence = case['n_Pd'] * 10 + case['n_D']
    if abs(calc.get_number_of_electrons() - expected_valence) > 1e-8:
        raise ValueError('Saved GPW has the wrong neutral Pd10/H1 electron count')
    expected_total = int(atoms.numbers.sum())
    rows = []
    volume = atoms.get_volume()
    for refinement in (1, 2):
        row = dict(grid_refinement=refinement)
        for skip_core, label, expected in ((False, 'all_electron', expected_total),
                                           (True, 'PAW_valence', expected_valence)):
            density = calc.get_all_electron_density(gridrefinement=refinement, skip_core=skip_core)
            integral = float(density.mean() * volume)
            if not np.isfinite(integral):
                raise ValueError('Nonfinite actual density')
            if not skip_core and abs(integral - expected) > 1e-6 * expected:
                raise ValueError('Actual all-electron grid integral/units do not match the saved target')
            row[label] = dict(shape=list(density.shape), integral_electrons=integral,
                              expected_electrons=expected,
                              relative_integral_residual=abs(integral-expected)/expected,
                              midpoint_e_per_A3=None if target is None else
                              sample_periodic_density(density, atoms.cell, target))
        rows.append(row)
    change = {}
    if target is not None:
        for label in ('all_electron', 'PAW_valence'):
            coarse, fine = [row[label]['midpoint_e_per_A3'] for row in rows]
            change[label] = abs(fine-coarse)/max(abs(fine), 1e-30)
    return dict(schema='coherence-pairs-actual-density-v1', case=case,
                source_point=manifest(path), source_execution=manifest(execution_path),
                source_bank=point['bank'], periodic_pairs=pairs,
                closest_pair=None if not pairs else pairs[0], density_rows=rows,
                midpoint_relative_grid_change=change,
                midpoint_grid_target_met=bool(change) and max(change.values()) <= 0.05,
                raw_fmax_eV_A=point['max_unconstrained_force_eV_A'],
                stationarity_target_eV_A=config['settings']['fmax_eV_A'],
                saved_point_only=True, density_units='electrons/angstrom^3',
                density_convention='GPAW PAW reconstructed all-electron and frozen-core-excluded valence; spin sum',
                density_is_not_screening_energy=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--point', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    start = time.process_time()
    result = read_point(args.point)
    result['readout_process_cpu_s'] = time.process_time()-start
    save(args.output, result)
    print(json.dumps({key: result[key] for key in ('case', 'closest_pair',
          'midpoint_relative_grid_change', 'readout_process_cpu_s')}))


if __name__ == '__main__':
    main()
