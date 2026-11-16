"""Neutral periodic PW production, with explicitly constrained bulk relaxation."""

import gpaw
from ase import Atoms
from gpaw import FermiDirac, GPAW, PW
from gpaw.mpi import serial_comm

from screen_cases import geometry as geometry_record
from screen_dos import analyse_dos
from screen_io import manifest, require, write_json
from screen_relax import relax_bulk


def ground_state(config, output):
    seed = config['geometry']
    atoms = Atoms(symbols=seed['symbols'], positions=seed['positions_A'],
                  cell=seed['cell_A'], pbc=True)
    atoms.set_initial_magnetic_moments(seed['initial_magnetic_moments_muB'])
    settings = config['ground_state']
    calc = GPAW(legacy_gpaw=False, communicator=serial_comm,
                mode=PW(settings['ecut_eV']), xc='PBE',
                kpts={'size': settings['kmesh'], 'gamma': True},
                occupations=FermiDirac(settings['occupation_width_eV']),
                setups=settings['setups'], nbands=settings['scf_nbands'],
                spinpol=settings['spinpol'], symmetry={'symmorphic': True},
                convergence=settings['convergence'], maxiter=150,
                parallel={'domain': 1, 'band': 1, 'kpt': 1},
                txt=str(output/'ground.txt'))
    atoms.calc = calc
    relaxation = relax_bulk(atoms, config['relaxation'], output)
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    electrons = float(calc.get_number_of_electrons())
    require(electrons == config['valence_electrons'],
            'Actual PAW electrons differ from the registered per-cell count')
    forces = atoms.get_forces()
    metadata = dict(schema='coherence-screen-ground-v1',
                    case=config['case'], model=config['model'],
                    ground_state_settings=settings,
                    initial_geometry=seed, geometry=geometry_record(atoms),
                    relaxation=relaxation,
                    cell_volume_A3=float(atoms.get_volume()),
                    neutral_electrons=electrons, free_energy_eV=free_energy,
                    forces_eV_A=forces.tolist(),
                    stress_eV_A3=atoms.get_stress().tolist(),
                    magnetic_moment_muB=float(calc.get_magnetic_moment()),
                    gpaw_version=gpaw.__version__)
    calc.diagonalize_full_hamiltonian(nbands=settings['saved_nbands'])
    require(calc.get_number_of_bands() == settings['saved_nbands'],
            'Requested empty states were not saved')
    gpw = output/'ground.gpw'
    calc.write(str(gpw), mode='all')
    metadata['gpw'] = manifest(gpw)
    atoms.calc = None
    del calc
    metadata['dos'] = analyse_dos(gpw, atoms, config, output, electrons)
    write_json(output/'ground.json', metadata)
    return metadata
