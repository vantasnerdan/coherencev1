"""Bounded Native states, raw energies and byte-exact original host reuse."""
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import time

import numpy as np
from ase.units import Hartree
from gpaw import GPAW, FermiDirac, Mixer
from gpaw.mpi import serial_comm
from gpaw.utilities import h2gpts

from contact_setup import ContactHeliumSetup


def require(condition, message):
    if not condition:
        raise ValueError(message)


def binding(path, root):
    path, root = Path(path), Path(root)
    data = path.read_bytes()
    return dict(path=path.relative_to(root).as_posix(), bytes=len(data),
                sha256=hashlib.sha256(data).hexdigest())


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')


def check_stop(stop):
    require(datetime.now(timezone.utc) < datetime.fromisoformat(stop.replace('Z', '+00:00')),
            'Commissioned absolute stop reached; no further Native work admitted')


def validate_state(record, atoms, charge, h_A):
    require(record['schema'] == 'hydride-native-state-v1', 'Unknown state schema')
    require(record['charge'] == charge and record['specification']['h_A'] == h_A,
            'State charge/grid does not match requested branch')
    require(record['symbols'] == atoms.get_chemical_symbols(), 'State nuclear inventory changed')
    require(np.allclose(np.asarray(record['positions_A']).reshape(-1, 3), atoms.positions,
                        rtol=0, atol=1e-12) and
            np.allclose(record['cell_A'], atoms.cell.array, rtol=0, atol=1e-12),
            'Frozen geometry changed')
    expected = sum(2 if symbol == 'He' else 1 for symbol in record['symbols']) - charge
    require(abs(record['electrons'] - expected) < 1e-7, 'Wrong electron inventory')
    if len(atoms):
        require(record['backend'] == 'gpaw-native' and record['scf_converged'],
                'Host/state is not an actually converged Native state')
    else:
        require(record['backend'] == 'analytic-empty-vacuum' and charge == 0 and
                all(value == 0 for value in record['energy'].values()),
                'Only the exact empty neutral reference is analytic')
    energy = record['energy']
    require(all(math.isfinite(float(x)) for x in energy.values()), 'Nonfinite raw energy')
    require(abs(energy['internal_energy_eV'] - energy['free_energy_eV'] +
                energy['minus_TS_eV']) < 1e-9, 'Wrong internal/free/entropy convention')


def reuse_state(source, atoms, charge, h_A, output, stem):
    source, output = Path(source), Path(output)
    record = json.loads(source.read_bytes())
    validate_state(record, atoms, charge, h_A)
    require(record['setups'] == {'H': 'ae'}, 'Original host has a different H operator')
    destination = output/'states'/(stem + '.json')
    destination.parent.mkdir(parents=True, exist_ok=True)
    require(not destination.exists(), 'Refusing to replace an existing state')
    shutil.copyfile(source, destination)
    # Original record bytes, including original GPW/source descriptors, survive.
    # GPW files are intentionally not copied or represented as local recovery.
    return record | {'state_file': binding(destination, output)}


def calculator(atoms, spec, charge, log, helium):
    symbols = atoms.get_chemical_symbols()
    require(set(symbols) <= {'H', 'He', 'Li'}, 'Unexpected contact species')
    electrons = sum(2 if s == 'He' else 1 for s in symbols) - charge
    require(electrons >= 0 and charge in (0, 1), 'Unsupported contact electron branch')
    setups = {'H': 'ae'}
    if 'He' in symbols:
        require(symbols.count('He') == 1, 'Contact must have exactly one Z=2 nucleus')
        setups['He'] = helium
    return GPAW(legacy_gpaw=True, mode='fd',
                gpts=h2gpts(spec['h_A'], atoms.cell, idiv=4), setups=setups,
                xc='PBE', charge=charge,
                nbands=max(4, int(np.ceil(electrons / 2)) + spec['nband_padding']),
                kpts={'size': spec['kmesh'], 'gamma': True},
                symmetry={'point_group': False, 'time_reversal': False},
                occupations=FermiDirac(0.1), mixer=Mixer(0.05, 5, 100),
                convergence={'energy': 1e-6, 'density': 1e-5, 'eigenstates': 1e-8},
                maxiter=150, spinpol=False,
                experimental={'reuse_wfs_method': 'keep'},
                parallel={'domain': 1, 'band': 1, 'kpt': 1},
                communicator=serial_comm, txt=str(log))


def scf(atoms, spec, charge, output, stem, sources, stop):
    output = Path(output)
    require(len(atoms) > 0, 'Empty vacuum must reuse the original exact zero record')
    check_stop(stop)
    state_dir = output/'states'
    state_dir.mkdir(parents=True, exist_ok=True)
    log = state_dir/(stem + '.txt')
    require(not log.exists(), 'Refusing to overwrite a Native log')
    helium = ContactHeliumSetup() if 'He' in atoms.get_chemical_symbols() else None
    atoms.calc = calculator(atoms, spec, charge, log, helium)
    atoms.calc.attach(lambda: check_stop(stop), 1)
    started = time.process_time()
    print('CONTACT_SCF_START state=' + stem + ' q=' + str(charge), flush=True)
    atoms.get_potential_energy()
    calc = atoms.calc
    require(bool(calc.scf.converged), 'Native SCF did not converge')
    valence = {}
    for atom, setup in zip(atoms, calc.wfs.setups):
        expected = 2 if atom.symbol == 'He' else 1
        require(setup.Nv == expected, 'Wrong actual setup valence')
        if atom.symbol in ('H', 'He'):
            require(setup.Z == expected and setup.Nc == 0 and setup.ni == 0 and
                    setup.type == 'all-electron' and setup.E == 0 and
                    setup.alpha1 == 10 and setup.alpha2 == 300,
                    'Wrong actual all-electron contact operator')
        if atom.symbol == 'He':
            require(setup is helium and not hasattr(setup, 'build') and
                    setup.f_j == [2.0] and
                    abs(setup.Delta0 + 2 / math.sqrt(4 * math.pi)) < 1e-15 and
                    abs(setup.M + 4 * math.sqrt(10 / (2 * math.pi))) < 1e-15,
                    'Z=2 operator was wrapped, replaced or rescaled')
        valence[atom.symbol] = float(setup.Nv)
    hamiltonian = calc.hamiltonian
    energy = dict(internal_energy_eV=float((hamiltonian.e_total_free -
                                           hamiltonian.e_entropy) * Hartree),
                  free_energy_eV=float(hamiltonian.e_total_free * Hartree),
                  minus_TS_eV=float(hamiltonian.e_entropy * Hartree))
    setups = {'H': 'ae'}
    if helium is not None:
        setups['He'] = helium.descriptor(sources)
    count = calc.get_number_of_iterations()
    record = dict(schema='hydride-native-state-v1', backend='gpaw-native',
                  energy=energy, charge=charge, electrons=float(calc.get_number_of_electrons()),
                  symbols=atoms.get_chemical_symbols(), positions_A=atoms.positions.tolist(),
                  cell_A=atoms.cell.array.tolist(), masses_u=atoms.get_masses().tolist(),
                  realspace_shape=calc.wfs.gd.N_c.tolist(), actual_kpoints=int(calc.wfs.kd.nibzkpts),
                  nbands=int(calc.get_number_of_bands()), specification=spec,
                  scf_converged=True, scf_iterations=None if count is None else int(count),
                  xc=calc.parameters.xc, setups=setups, setup_valence=valence,
                  spinpol=bool(calc.parameters.spinpol), occupation_width_eV=.1,
                  wavefunction_predictor=type(calc.wfs.wfs_mover).__name__,
                  process_cpu_s=time.process_time() - started,
                  recovery='raw converged energies/geometry/log only; no GPW serialization or wavefunction restart claimed')
    validate_state(record, atoms, charge, spec['h_A'])
    calc.close()
    calc.log.close()
    record['log'] = binding(log, output)
    path = state_dir/(stem + '.json')
    write_json(path, record)
    atoms.calc = None
    print('CONTACT_SCF_COMPLETE state=' + stem, flush=True)
    return record | {'state_file': binding(path, output)}
