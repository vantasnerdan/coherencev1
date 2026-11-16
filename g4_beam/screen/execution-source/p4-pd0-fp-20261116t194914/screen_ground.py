"""Saved neutral PW ground state, exact empty states, and projector-defined DOS."""

import gpaw
import numpy as np
from ase import Atoms
from gpaw import FermiDirac, GPAW, PW
from gpaw.dos import DOSCalculator
from gpaw.mpi import serial_comm
from scipy.special import expit

from screen_io import manifest, require, write_json, write_npz


def ground_state(config, output):
    geometry = config['geometry']
    atoms = Atoms(symbols=geometry['symbols'], positions=geometry['positions_A'],
                  cell=geometry['cell_A'], pbc=True)
    settings = config['ground_state']
    calc = GPAW(legacy_gpaw=False, communicator=serial_comm,
                mode=PW(settings['ecut_eV']), xc='PBE',
                kpts={'size': settings['kmesh'], 'gamma': True},
                occupations=FermiDirac(settings['occupation_width_eV']),
                setups=settings['setups'], nbands=settings['scf_nbands'],
                convergence=settings['convergence'], maxiter=150,
                parallel={'domain': 1, 'band': 1, 'kpt': 1},
                txt=str(output/'ground.txt'))
    atoms.calc = calc
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces()
    calc.diagonalize_full_hamiltonian(nbands=settings['saved_nbands'])
    require(calc.get_number_of_bands() == settings['saved_nbands'],
            'Requested empty states were not saved')
    gpw = output/'ground.gpw'
    calc.write(str(gpw), mode='all')
    metadata = dict(schema='coherence-screen-ground-v1',
                    case=config['case'], model=config['model'],
                    ground_state_settings=settings, geometry=geometry,
                    cell_volume_A3=float(atoms.get_volume()),
                    neutral_electrons=float(calc.get_number_of_electrons()),
                    free_energy_eV=free_energy,
                    forces_eV_A=forces.tolist(),
                    gpaw_version=gpaw.__version__, gpw=manifest(gpw))
    atoms.calc = None
    del calc
    # Reload the saved producer. DOS construction may shift its eigenvalue array;
    # it must never mutate the calculator subsequently used for response.
    dos_calc = GPAW(str(gpw), legacy_gpaw=False,
                    communicator=serial_comm, txt=None)
    dos = DOSCalculator.from_calculator(dos_calc)
    require(dos.collinear and dos.nspins == 1,
            'This fingerprint requires a nonmagnetic collinear state')
    eigen = dos.eig_skn
    weights = dos.weight_k
    count = float(dos.degeneracy*np.einsum(
        'k,skn->', weights,
        expit(-eigen/settings['occupation_width_eV'])))
    require(abs(count-metadata['neutral_electrons']) < 1e-5,
            'Saved-band occupations do not reproduce the electron count')
    n_ef = {}
    for width in config['dos']['widths_eV']:
        value = float(dos.raw_dos([0.0], width=width)[0])
        independent = float(dos.degeneracy*np.einsum(
            'k,skn->', weights, np.exp(-(eigen/width)**2))
            / (np.sqrt(np.pi)*width))
        require(np.isclose(value, independent, rtol=1e-10, atol=1e-12),
                'DOS spin, energy, or k-weight normalization mismatch')
        n_ef[str(width)] = value
    energies = np.arange(*config['dos']['occupied_window_eV'])
    pdos = sum(dos.raw_pdos(energies, a=a, l=2,
                           width=config['dos']['moment_width_eV'])
               for a in config['dos']['metal_indices'])
    occupied_d = float(np.trapezoid(pdos, energies))
    require(occupied_d > 0 and np.min(pdos) >= -1e-12,
            'Invalid occupied metal d projection')
    centre = float(np.trapezoid(energies*pdos, energies)/occupied_d)
    width = float(np.sqrt(np.trapezoid((energies-centre)**2*pdos,
                                     energies)/occupied_d))
    metadata['dos'] = dict(
        settings=config['dos'], fermi_energy_eV=float(dos.fermi_level),
        n_EF_states_per_eV_per_cell=n_ef,
        occupied_d_states_per_cell=occupied_d,
        occupied_d_centre_relative_EF_eV=centre,
        occupied_d_rms_width_eV=width,
        D_projection_present=any(s == 'H' for s in atoms.symbols),
        numerical_electron_count=count,
        eligible=False,
        reason='Fingerprint only; k-mesh and observable convergence not established',
        spectrum=write_npz(output/'dos.npz', energies_eV=energies,
                           occupied_metal_d_states_per_eV=pdos))
    write_json(output/'ground.json', metadata)
    return metadata
