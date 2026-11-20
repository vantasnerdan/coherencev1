"""One actual AE energy of an explicitly unrelaxed embedded-cell sensitivity."""
import argparse
import json
from pathlib import Path
import numpy as np
from ase import Atom, Atoms
from defect_calculator import scf, manifest, write_json, require
from defect_geometry import MASS_D
from run_defect import validate_epoch

ROLES = ('occupied-host', 'vacant-host', 'occupied-probe', 'vacant-probe')


def run(input_path, output, role):
    input_binding, source_binding = validate_epoch(input_path)
    spec = json.loads(input_path.read_bytes())
    require(spec['case'] == 'C' and spec['repeat_primitive'] == 3 and spec['original_Pd_sites'] == 27,
            'Wrong embedded-cell sensitivity')
    require(spec['sensitivity_scope'] == 'unrelaxed_embedded_local_geometry', 'No converged-cell claim is allowed')
    atoms = Atoms(spec['symbols'], positions=spec['initial_positions_A'], cell=spec['cell_A'], pbc=True)
    atoms.set_masses([MASS_D if atom.symbol == 'H' else atom.mass for atom in atoms])
    target_index = spec['target_index']
    require(atoms[target_index].symbol == 'H' and atoms.get_chemical_symbols().count('Pd') == 26 and
            atoms.get_chemical_symbols().count('H') == 1, 'Wrong target/inventory')
    target = atoms.positions[target_index].copy()
    direction = np.asarray(spec['probe_direction'])
    require(abs(np.linalg.norm(direction)-1) < 1e-12, 'Nonunit frozen ray')
    if role.startswith('vacant'):
        del atoms[target_index]
    charge = int(role.endswith('probe'))
    if charge:
        atoms += Atom('H', position=target-.1*direction, mass=MASS_D)
    output.mkdir(parents=True, exist_ok=False)
    print('DEFECT_RUN_START case=C stage=production role='+role, flush=True)
    state = scf(atoms, spec, charge, output, 'state')
    result = dict(schema='pd-defect-cell-role-v1', role=role, completed=True,
                  input=input_binding, frozen_source=source_binding, state=manifest(output/'state.json'),
                  forces_computed=False, force_converged=False,
                  scope='One converged electronic energy at fixed embedded geometry; no ionic relaxation or cell convergence claim.')
    write_json(output/'cell-role.json', result)
    print('CELL_ROLE_COMPLETE role='+role+' cpu_s='+str(state['process_cpu_s']), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--role', choices=ROLES, required=True)
    args = parser.parse_args()
    run(args.input, args.output, args.role)
