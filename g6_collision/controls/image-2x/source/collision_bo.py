"""Independent frozen-geometry BO energies for actual saved TD checkpoints."""
import json
import gc
from ase import Atoms
from screen_io import manifest, write_json
from substrate_framework.collision_observables import excitation_and_closure
from collision_ground import make_calculator, energy_record


def host_reference(spec, ground, output):
    atoms = Atoms(ground['symbols'], positions=ground['positions_A'],
                  cell=ground['cell_A'], pbc=True)
    del atoms[ground['path']['projectile_index']]
    atoms.calc = make_calculator(atoms, spec, output/'host-reference.txt', charge=0)
    print('COLLISION_NEUTRAL_HOST_START', flush=True)
    atoms.get_potential_energy()
    saved = output/'host-reference.gpw'
    atoms.calc.write(str(saved), mode='all')
    result = dict(charge=0, energy=energy_record(atoms.calc), gpw=manifest(saved),
                  electrons=float(atoms.calc.get_number_of_electrons()),
                  symbols=atoms.get_chemical_symbols(), positions_A=atoms.positions.tolist(),
                  cell_A=atoms.cell.array.tolist(),
                  scope='Neutral host-only state; same frozen host/settings, projectile removed')
    write_json(output/'host-reference.json', result)
    print('COLLISION_NEUTRAL_HOST_COMPLETE', flush=True)
    del atoms
    gc.collect()
    return result


def bo_checkpoints(spec, ground_dir, trajectory_dir, output):
    output.mkdir(parents=True, exist_ok=True)
    ground = json.loads((ground_dir/'ground.json').read_bytes())
    trajectory = json.loads((trajectory_dir/'trajectory.json').read_bytes())
    selected = spec['bo_checkpoints_A']
    host = host_reference(spec, ground, output)
    rows = []
    for checkpoint in trajectory['checkpoints']:
        if not any(abs(checkpoint['R_A']-r) < 1e-9 for r in selected):
            continue
        index = checkpoint['index']
        td = trajectory['rows'][checkpoint['row_index']]
        destination = output/f'bo-{index:02d}.json'
        if index == 0:
            energy = ground['energy']
            gpw = ground['ground_gpw']
        else:
            atoms = Atoms(ground['symbols'], positions=checkpoint['positions_A'],
                          cell=ground['cell_A'], pbc=True)
            atoms.calc = make_calculator(atoms, spec, output/f'bo-{index:02d}.txt', charge=1)
            print(f'COLLISION_BO_START index={index} R_A={checkpoint["R_A"]}', flush=True)
            atoms.get_potential_energy()
            energy = energy_record(atoms.calc)
            saved = output/f'bo-{index:02d}.gpw'
            atoms.calc.write(str(saved), mode='all')
            gpw = manifest(saved)
            del atoms
            gc.collect()
        row = dict(R_A=checkpoint['R_A'], row_index=checkpoint['row_index'],
                   bo=energy, td_internal_energy_eV=td['td_internal_energy_eV'],
                   external_work_eV=td['external_work_eV'], gpw=gpw,
                   actual_checkpoint=manifest(trajectory_dir/f'checkpoint-{index:02d}.json'))
        write_json(destination, row)
        rows.append(row)
        print(f'COLLISION_BO_COMPLETE index={index}', flush=True)
    excitations, residuals = excitation_and_closure(
        [r['td_internal_energy_eV'] for r in rows],
        [r['bo']['internal_energy_eV'] for r in rows],
        [r['external_work_eV'] for r in rows],
        initial_excitation_eV=rows[0]['td_internal_energy_eV']-rows[0]['bo']['internal_energy_eV'])
    for row, excitation, residual in zip(rows, excitations, residuals):
        row['excitation_internal_energy_eV'] = float(excitation)
        row['closure_residual_eV'] = float(residual)
    result = dict(schema='pdd-collision-bo-excitation-v1', specification=spec, checkpoints=rows,
                  interpretation='Internal-energy excess relative to independently equilibrated BO at Fermi width0.1eV; signed, no clipping',
                  neutral_host_reference=host,
                  trajectory=manifest(trajectory_dir/'trajectory.json'))
    write_json(output/'excitation.json', result)
    return result
