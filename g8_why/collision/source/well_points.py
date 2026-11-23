"""Cold BO points and the same-electron-count field for one frozen PdD well."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from well_native import atoms_from, bound, point
from collision_ground import make_calculator
from screen_io import write_json
from well_density import publish_field


def run(spec, output):
    output.mkdir(parents=True, exist_ok=False)
    baseline = json.loads(bound(spec['baseline_record']).read_bytes())
    target = int(baseline['target_index'])
    started = time.process_time()
    rows = []
    vacancy_field = bool(spec.get('vacancy_field', False))
    for i, offset in enumerate(spec['offsets_A']):
        atoms = atoms_from(baseline)
        positions = atoms.positions.copy()
        if vacancy_field:
            if np.any(offset):
                raise ValueError('Vacancy field has no target nucleus to displace')
            del atoms[target]
            positions = atoms.positions.copy()
        else:
            positions[target] += np.asarray(offset)
        positions += np.asarray(spec.get('rigid_translation_A', [0., 0., 0.]))
        atoms.set_positions(positions)
        calc = make_calculator(atoms, spec['scf_specification'], output/f'point-{i:02d}.txt',
                               charge=-1 if vacancy_field else 0)
        if spec['hydrogen_setup'] == 'default_PAW':
            calc = calc.new(setups={'Pd': '10'}, txt=str(output/f'point-{i:02d}.txt'))
        elif spec['hydrogen_setup'] != 'ae':
            raise ValueError('Unspecified physical H representation')
        atoms.calc = calc
        row = point(calc, atoms, f'point-{i:02d}', output)
        row['offset_A'] = offset
        row['target_force_eV_A'] = None if vacancy_field else row['forces_eV_A'][target]
        if vacancy_field:
            if abs(row['electrons_valence']-44) > 1e-8:
                raise ValueError('The auxiliary vacancy must retain44 valence electrons')
            field = publish_field(calc, atoms, output/f'field-{i:02d}.npz')
            row['valence_field'] = field['binding']
            row['density_valence_electrons'] = field['electron_counts']['valence']
            row['density_decomposition'] = field
        rows.append(row)
        del calc, atoms
    write_json(output/'points.json', dict(schema='g8-well-points-v1', specification=spec,
        target_index=target, points=rows, process_cpu_s=time.process_time()-started,
        scope='Same44-electron frozen PdD; target removed at charge-1 only for the auxiliary field'
              if vacancy_field else 'Frozen neutral PdD well; no density/profile fit or pair campaign'))
    print('G8_WELL_'+'POINTS_COMPLETE', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(json.loads(args.input.read_bytes()), args.output)
