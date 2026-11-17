"""One bounded force relaxation; retain actual completed primitives and trajectories."""
import json
import math
from pathlib import Path
import numpy as np
from ase.io import write
from ase.optimize import BFGS
from gpaw import GPAW, PW, FermiDirac
from gpaw.mpi import serial_comm
from vacancy_geometry import make_geometry


def calculate(config, case, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    atoms, geometry = make_geometry(config, case)
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
    optimizer = BFGS(atoms, trajectory=str(output / 'relax.traj'),
                     logfile=str(output / 'relax.log'), maxstep=0.15)
    converged = bool(optimizer.run(fmax=config['fmax_eV_A'], steps=config['relax_steps']))
    energy = float(atoms.get_potential_energy())
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces()
    write(output / 'final.extxyz', atoms)
    calc.write(str(output / 'final.gpw'))
    result = dict(schema='coherence-vacancy-primitive-v1', case=case,
                  config=config, geometry=geometry, nbands=nbands,
                  neutral_electrons=float(calc.get_number_of_electrons()),
                  energy_eV=energy, electronic_free_energy_eV=free_energy,
                  energy_convention='GPAW zero-width extrapolated electronic energy; no ionic ZPE/entropy',
                  converged_relaxation=converged, relaxation_steps=optimizer.nsteps,
                  max_constrained_force_eV_A=float(np.linalg.norm(forces, axis=1).max()),
                  unconstrained_forces_eV_A=atoms.get_forces(apply_constraint=False).tolist(),
                  final_total_moment_muB=float(atoms.get_magnetic_moment()),
                  final_local_moments_muB=atoms.get_magnetic_moments().tolist(),
                  final_positions_A=atoms.positions.tolist(), cell_A=atoms.cell.array.tolist(),
                  symbols=atoms.get_chemical_symbols(), eligible=False)
    (output / 'primitive.json').write_text(json.dumps(result, indent=2) + '\n')
    return result
