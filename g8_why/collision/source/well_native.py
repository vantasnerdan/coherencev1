"""Collect actual local-D BO forces and same-electron-count vacancy density.

The pilot contains two new SCFs only: D displaced 0.1 A, then vacancy charge -1.
All Pd and other D coordinates remain fixed. Cached neutral host is not rerun.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from ase import Atoms
from ase.units import Hartree
from screen_io import manifest, write_json, write_npz
from collision_ground import make_calculator


def bound(item):
    path = Path(item['path'])
    data = path.read_bytes()
    if len(data) != item['bytes'] or hashlib.sha256(data).hexdigest() != item['sha256']:
        raise ValueError(f'Original input changed: {path.name}')
    return path


def energy(calc):
    h = calc.hamiltonian
    return dict(free_energy_eV=float(h.e_total_free * Hartree),
                internal_energy_eV=float((h.e_total_free-h.e_entropy) * Hartree),
                minus_TS_eV=float(h.e_entropy * Hartree))


def point(calc, atoms, label, output):
    started = time.process_time()
    print(f'G8_WELL_SCF_START label={label}', flush=True)
    atoms.get_potential_energy()
    forces = atoms.get_forces()
    gpw = output / (label + '.gpw')
    calc.write(str(gpw), mode='all')
    result = dict(label=label, energy=energy(calc), forces_eV_A=forces.tolist(),
                  positions_A=atoms.positions.tolist(), cell_A=atoms.cell.array.tolist(),
                  symbols=atoms.get_chemical_symbols(), electrons_valence=float(calc.get_number_of_electrons()),
                  gpw=manifest(gpw), process_cpu_s=time.process_time()-started)
    write_json(output / (label + '.json'), result)
    print(f'G8_WELL_POINT_COMPLETE label={label} cpu_s={result["process_cpu_s"]:.6g}', flush=True)
    return result


def atoms_from(record):
    atoms = Atoms(record['symbols'], positions=record['positions_A'],
                  cell=record['cell_A'], pbc=True)
    atoms.set_masses([2.014101778 if a.symbol == 'H' else a.mass for a in atoms])
    return atoms


def pilot(spec, output):
    started = time.process_time()
    output.mkdir(parents=True, exist_ok=False)
    baseline = json.loads(bound(spec['baseline_record']).read_bytes())
    atoms = atoms_from(baseline)
    target = int(baseline['target_index'])
    if abs(baseline['electrons_valence']-44) > 1e-8:
        raise ValueError('Changed neutral host electron count')
    write_json(output/'baseline.json', baseline)
    print('G8_WELL_CACHED_HOST_READ', flush=True)
    positions = atoms.positions.copy()
    positions[target] += [0.1, 0., 0.]
    atoms.set_positions(positions)
    # Hae has no movable PAW partial waves: cold SCF preserves the Hamiltonian.
    occ = make_calculator(atoms, spec['scf_specification'], output/'occupied.txt', charge=0)
    atoms.calc = occ
    displaced = point(occ, atoms, 'axis-010', output)
    f0 = np.asarray(baseline['forces_eV_A'])[target]
    f1 = np.asarray(displaced['forces_eV_A'])[target]
    stiffness = float(-(f1[0]-f0[0])/.1)
    write_json(output/'local-curvature.json', dict(target_index=target,
        force_secant_eV_A2=stiffness, displacement_A=.1,
        delta_free_energy_eV=displaced['energy']['free_energy_eV']-baseline['energy']['free_energy_eV'],
        delta_internal_energy_eV=displaced['energy']['internal_energy_eV']-baseline['energy']['internal_energy_eV'],
        convention='Forces differentiate electronic free energy at saved smearing; internal also retained',
        scope='One finite-displacement local force secant, not a converged phonon Hessian'))
    del occ, atoms
    host = json.loads(bound(spec['vacant_host_record']).read_bytes())
    vatoms = atoms_from(host)
    vac = make_calculator(vatoms, spec['scf_specification'], output/'vacancy-44.txt', charge=-1)
    vatoms.calc = vac
    vacancy = point(vac, vatoms, 'vacancy-44', output)
    if abs(vacancy['electrons_valence']-44) > 1e-8:
        raise ValueError('Vacancy must retain all44 host valence electrons')
    density = vac.get_all_electron_density(gridrefinement=2)
    integrated = float(density.sum()*abs(np.linalg.det(vatoms.cell.array))/density.size)
    receipt = write_npz(output/'vacancy-density.npz', density_e_A3=density,
        cell_A=vatoms.cell.array, positions_A=vatoms.positions,
        nuclear_Z=np.asarray(vatoms.get_atomic_numbers(), dtype=int))
    radial = {}
    setup = vac.wfs.setups[next(i for i,a in enumerate(vatoms) if a.symbol == 'Pd')]
    data = getattr(setup, 'data', None)
    if data is not None and hasattr(data, 'nc_g') and hasattr(data, 'rgd'):
        radial['raw_nc_g'] = np.asarray(data.nc_g)
        radial['raw_r_g_bohr'] = np.asarray(data.rgd.r_g)
    radial_receipt = write_npz(output/'pd-core-raw.npz', **radial) if radial else None
    result = dict(schema='g8-local-well-pilot-v1', authority=spec['authority'],
        specification=spec, baseline='baseline.json', displaced='axis-010.json',
        curvature='local-curvature.json', vacancy='vacancy-44.json', density=receipt,
        integrated_all_electrons=integrated, expected_all_electrons=188,
        vacant_nuclear_charge=int(vatoms.numbers.sum()), net_cell_charge=-1,
        raw_core_radial=radial_receipt,
        core_radial_scope='Raw GPAW arrays; their radial normalization is not inferred here',
        setup_attributes=sorted(vars(setup)),
        data_attributes=sorted(vars(data)) if data is not None else [],
        process_cpu_s=time.process_time()-started,
        scope='Local neutral-PdD BO secant and frozen vacancy-electron field, not a completed nuclear well or tail')
    write_json(output/'pilot.json', result)
    print('G8_WELL_'+'PILOT_COMPLETE', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    pilot(json.loads(args.input.read_bytes()), args.output)
