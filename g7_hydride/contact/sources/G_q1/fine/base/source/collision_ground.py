"""Prepare the charged moving-proton electronic initial state; keep host ions frozen."""
import numpy as np
from ase.units import Hartree
from gpaw import GPAW, FermiDirac, Mixer
from gpaw.mpi import serial_comm
from gpaw.utilities import h2gpts
from screen_io import write_json, manifest
from collision_geometry import collision_atoms


def make_calculator(atoms, spec, txt, *, charge):
    nvalence = sum(10 if a.symbol == 'Pd' else 1 for a in atoms) - charge
    return GPAW(legacy_gpaw=True, mode='fd',
                gpts=h2gpts(spec['h_A'], atoms.cell, idiv=4),
                setups={'Pd': '10', 'H': 'ae'}, xc='PBE', charge=charge,
                nbands=int(np.ceil(nvalence / 2.0)) + 8 * int(np.prod(spec['repeat'])),
                kpts={'size': spec['kmesh'], 'gamma': True},
                symmetry={'point_group': False, 'time_reversal': False},
                occupations=FermiDirac(0.1), mixer=Mixer(0.05, 5, 100),
                convergence={'energy': 1e-6, 'density': 1e-5, 'eigenstates': 1e-8},
                maxiter=150, spinpol=False,
                parallel={'domain': 1, 'band': 1, 'kpt': 1},
                communicator=serial_comm, txt=str(txt))


def energy_record(calc):
    hamiltonian = calc.hamiltonian
    # GPAW's internal energy omits the negative T*S term in free energy.
    internal = (hamiltonian.e_total_free - hamiltonian.e_entropy) * Hartree
    return dict(internal_energy_eV=float(internal),
                free_energy_eV=float(hamiltonian.e_total_free * Hartree),
                minus_TS_eV=float(hamiltonian.e_entropy * Hartree))


def prepare(spec, output):
    output.mkdir(parents=True, exist_ok=True)
    atoms, path = collision_atoms(spec)
    atoms.calc = make_calculator(atoms, spec, output / 'ground.txt', charge=1)
    print('COLLISION_GROUND_START', flush=True)
    atoms.get_potential_energy()
    calc = atoms.calc
    calc.write(str(output / 'ground.gpw'), mode='all')
    actual = dict(specification=spec, path=path, energy=energy_record(calc),
                  positions_A=atoms.positions.tolist(), cell_A=atoms.cell.array.tolist(),
                  symbols=atoms.get_chemical_symbols(), electrons=float(calc.get_number_of_electrons()),
                  realspace_shape=calc.wfs.gd.N_c.tolist(),
                  actual_kpoints=int(calc.wfs.kd.nibzkpts),
                  ground_gpw=manifest(output / 'ground.gpw'))
    write_json(output / 'ground.json', actual)
    print('COLLISION_GROUND_COMPLETE', flush=True)
    return actual
