"""Saved-geometry numerical controls and explicitly changed, source-bound seeds."""
import json
import math
from pathlib import Path
import numpy as np
from ase import Atom
from ase.io import read, write
from ase.optimize import BFGS
from vacancy_geometry import make_geometry
from gpaw import GPAW, PW, FermiDirac
from gpaw.mpi import serial_comm


def refine(config, case, source, output):
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    parent = json.loads((source / 'primitive.json').read_text())
    if not parent['converged_relaxation']:
        raise ValueError('Numerical control requires actually relaxed parent primitive')
    atoms = read(source / 'final.extxyz')
    # Extxyz attaches parent-sized SinglePointCalculator arrays.  The new
    # calculation must not carry those results across an atom-count change.
    atoms.calc = None
    moments = list(parent['final_local_moments_muB'])
    geometry = dict(parent['geometry'])
    if 'seed_D_positions_A' in config:
        if parent['symbols'].count('H') or case != 'vacancy_D1':
            raise ValueError('Changed seed requires an actual empty metal vacancy')
        for position in config['seed_D_positions_A']:
            atoms.append(Atom('H', position, mass=2.01410177812))
            moments.append(0.1)
        geometry = make_geometry(config, case)[1]
    atoms.set_initial_magnetic_moments(moments)
    geometry['initial_moments_muB'] = moments
    write(output / 'initial.extxyz', atoms)
    electrons = sum(10 if atom.symbol == 'Pd' else 1 for atom in atoms)
    nbands = math.ceil(electrons / 2) + max(8, math.ceil(electrons * 0.08))
    calc = GPAW(legacy_gpaw=False, communicator=serial_comm,
                mode=PW(config['ecut_eV']), xc='PBE', setups=config['setups'],
                kpts={'size': config['kmesh'], 'gamma': True},
                spinpol=config['spinpol'], charge=0, nbands=nbands,
                occupations=FermiDirac(config['occupation_width_eV']),
                convergence={'energy': 1e-5, 'density': 1e-4, 'eigenstates': 1e-7},
                maxiter=100, parallel={'domain': 1, 'band': 1, 'kpt': 1},
                txt=str(output / 'gpaw.txt'))
    atoms.calc = calc
    optimizer = None
    if config['relax_steps']:
        optimizer = BFGS(atoms, trajectory=str(output / 'relax.traj'),
                         logfile=str(output / 'relax.log'), maxstep=0.15)
        optimizer.run(fmax=config['fmax_eV_A'], steps=config['relax_steps'])
    energy = float(atoms.get_potential_energy())
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces()
    max_force = float(np.linalg.norm(forces, axis=1).max())
    write(output / 'final.extxyz', atoms)
    if optimizer is None:
        write(output / 'relax.traj', atoms)
    calc.write(str(output / 'final.gpw'))
    result = dict(schema='coherence-vacancy-primitive-v1', case=case,
                  config=config, geometry=geometry, nbands=nbands,
                  neutral_electrons=float(calc.get_number_of_electrons()),
                  energy_eV=energy, electronic_free_energy_eV=free_energy,
                  energy_convention='GPAW zero-width extrapolated electronic energy; no ionic ZPE/entropy',
                  converged_relaxation=bool(max_force <= config['fmax_eV_A']),
                  relaxation_steps=0 if optimizer is None else optimizer.nsteps,
                  numerical_control_at_saved_parent_coordinates=optimizer is None,
                  max_constrained_force_eV_A=max_force,
                  unconstrained_forces_eV_A=atoms.get_forces(apply_constraint=False).tolist(),
                  final_total_moment_muB=float(atoms.get_magnetic_moment()),
                  final_local_moments_muB=atoms.get_magnetic_moments().tolist(),
                  final_positions_A=atoms.positions.tolist(), cell_A=atoms.cell.array.tolist(),
                  symbols=atoms.get_chemical_symbols(), eligible=False)
    (output / 'primitive.json').write_text(json.dumps(result, indent=2) + '\n')
    return result
