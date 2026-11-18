"""Retained-wavefunction moment producer, never SCF; source identities are in ledger.json."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.special import expit
from ase.units import Bohr, Hartree
import gpaw
from gpaw import GPAW
from gpaw.mpi import serial_comm
from gpaw.response.df import DielectricFunction
from g4_beam.screen.dynamic._upstream.electronic_speed import fermi_rms_velocity_au


def binding(path):
    data = Path(path).read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def require_binding(record):
    actual = binding(record['path'])
    if actual['bytes'] != record['bytes'] or actual['sha256'] != record['sha256']:
        raise ValueError('Actual retained source byte binding mismatch: '+record['path'])


def one_case(case, output):
    for record in case['source_bindings']:
        require_binding(record)
    calc = GPAW(case['ground_gpw'], legacy_gpaw=False, communicator=serial_comm, txt=None)
    atoms = calc.get_atoms()
    if atoms.get_chemical_symbols() != case['geometry']['symbols'] or not np.allclose(atoms.cell.array, case['geometry']['cell_A'], rtol=0, atol=1e-10) or not np.allclose(atoms.positions, case['geometry']['positions_A'], rtol=0, atol=1e-10):
        raise ValueError('Retained ground is not the bound material geometry')
    nspins = calc.get_number_of_spins()
    if nspins != 1 or abs(calc.get_number_of_electrons()-case['expected_valence_electrons']) > 1e-8:
        raise ValueError('Retained ground is not the neutral nonmagnetic ten-valence target')
    if len(calc.get_bz_k_points()) != int(np.prod(case['kmesh'])):
        raise ValueError('Retained ground k grid differs from its declared source')
    bands = calc.get_number_of_bands()
    df = DielectricFunction(calc, frequencies=[0.], ecut=50, hilbert=False,
                            eta=.005, nbands=bands, intraband=True, rate='eta',
                            nblocks=1, world=serial_comm, integrationmode='point integration',
                            txt=str(output/'drude.txt'))
    optical = df.chi0calc.chi0_opt_ext_calc
    if optical.drude_calc is None:
        raise ValueError('No metallic intraband moment in the actual saved ground')
    drude = optical.drude_calc.calculate(optical.wd, optical.rate)
    tensor = drude.plasmafreq_vv
    width_Ha = df.gs.get_occupations_width()
    if abs(width_Ha*Hartree-case['FD_width_eV']) > 1e-12:
        raise ValueError('Saved finite-width occupations differ from the bound source')
    weights = np.asarray(calc.get_k_point_weights())
    if abs(weights.sum()-1) > 1e-12:
        raise ValueError('Expected normalized irreducible k weights')
    eigenvalues_eV = np.array([calc.get_eigenvalues(kpt=k, spin=0) for k in range(len(weights))])
    fermi_eV = calc.get_fermi_level()
    occupation = expit((fermi_eV-eigenvalues_eV)/(width_Ha*Hartree))
    electrons = float(2*np.sum(weights[:, None]*occupation))
    if abs(electrons-case['expected_valence_electrons']) > 1e-5:
        raise ValueError('Spin-counted finite-width eigenvalue occupations do not reproduce neutral valence')
    dos_per_Ha = float(2*np.sum(weights[:, None]*occupation*(1-occupation))/width_Ha)
    volume_Bohr3 = atoms.get_volume()/Bohr**3
    velocity = fermi_rms_velocity_au(tensor, dos_per_Ha, volume_Bohr3)
    spectrum = output/'occupation-spectrum.npz'
    np.savez_compressed(spectrum, k_weights=weights, eigenvalues_eV=eigenvalues_eV,
                        occupations_per_spin=occupation, plasma_squared_au=tensor)
    return dict(schema='actual-pass4-tenvalence-intraband-v1', case=case['case'],
                material=case['material'], geometry=case['geometry'], source_bindings=case['source_bindings'],
                no_new_SCF=True, no_new_solid=True, gpaw_version=gpaw.__version__,
                kmesh=case['kmesh'], saved_bands=bands, expected_electrons=case['expected_valence_electrons'],
                reconstructed_electrons=electrons, FD_width_eV=width_Ha*Hartree,
                D_T_states_per_Hartree_cell=dos_per_Ha,
                D_T_states_per_eV_cell=dos_per_Ha/Hartree, volume_Bohr3=volume_Bohr3,
                plasma_squared_au_real=tensor.real.tolist(), plasma_squared_au_imag=tensor.imag.tolist(),
                band_RMS_velocity_au=velocity, spectrum=binding(spectrum),
                U10_static_Yukawa_proxy_eV=case['U10_static_Yukawa_proxy_eV'],
                screening_scope='Accepted pass4 low-q ten-valence C/lambda Yukawa proxy, NOT a full large-q/contact dielectric inversion',
                velocity_scope='Actual intraband band-RMS finite-width moment; one isotropic Huke factor remains a model closure, not full dynamical response')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    specification = json.loads(args.input.read_text())
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    for record in specification['setup_bindings']:
        require_binding(record)
    gpaw.setup_paths[:] = [specification['setup_directory']]
    started = time.monotonic()
    results = []
    status = 'failed'
    try:
        for case in specification['cases']:
            print('TENVALENCE_DRUDE_START', case['case'], flush=True)
            directory = output/case['case']
            directory.mkdir()
            result = one_case(case, directory)
            (directory/'result.json').write_text(json.dumps(result, indent=2)+'\n')
            results.append(result)
            print('TENVALENCE_DRUDE_COMPLETE', case['case'], result['band_RMS_velocity_au'], flush=True)
        (output/'result.json').write_text(json.dumps(dict(cases=results), indent=2)+'\n')
        status = 'complete'
    finally:
        terminal = dict(status=status, elapsed_wall_s=time.monotonic()-started,
                        process_cpu_s=time.process_time(), finished_at=datetime.now(timezone.utc).isoformat(),
                        argv=sys.argv, no_new_SCF=True)
        (output/'terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')


if __name__ == '__main__':
    main()
