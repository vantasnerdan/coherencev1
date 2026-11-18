"""Reuse an authenticated saved ground state; optional empty-state expansion."""
from copy import deepcopy
import json
from pathlib import Path
import shutil

import numpy as np

from gpaw import GPAW
from gpaw.mpi import serial_comm
from screen_cases import geometry as geometry_record
from screen_io import manifest, require, write_json


def _same_geometry(seed, actual):
    if seed['symbols'] != actual['symbols']:
        return False
    cell = np.asarray(actual['cell_A'])
    if not np.allclose(seed['cell_A'], cell, rtol=0, atol=1e-9):
        return False
    delta = (np.asarray(seed['positions_A']) - actual['positions_A']) @ np.linalg.inv(cell)
    return np.max(np.abs(delta - np.rint(delta))) < 1e-9


def _saved_settings(original_config, calc):
    settings = deepcopy(original_config['ground_state'])
    require(calc.get_number_of_spins() == 1,
            'This registered family requires the actual spin-paired saved state')
    actual_spinpol = False
    if 'spinpol' in settings:
        require(settings['spinpol'] == actual_spinpol,
                'Original spin declaration differs from the actual saved state')
    settings['spinpol'] = actual_spinpol
    return settings


def prepare_saved_ground(specification, output):
    config = specification['configuration']
    source = Path(specification['saved_ground_directory'])
    original_config = json.loads((source / 'input.json').read_text())
    original = json.loads((source / 'ground.json').read_text())
    require(_same_geometry(config['geometry'], original['geometry']) and
            config['case'] == original_config['case'] and
            config['model'] == original_config['model'],
            'Saved actual endpoint differs from the registered fixed response host')
    calc = GPAW(str(source / 'ground.gpw'), legacy_gpaw=False,
                communicator=serial_comm, txt=str(output / 'saved-state-read.txt'))
    require(calc.atoms.pbc.all() and
            _same_geometry(config['geometry'], geometry_record(calc.atoms)),
            'Actual saved GPW does not contain the registered periodic endpoint')
    require(calc.get_number_of_electrons() == original['neutral_electrons'],
            'Actual saved-state electron count differs')
    source_settings = _saved_settings(original_config, calc)
    expected = deepcopy(config['ground_state'])
    expected['saved_nbands'] = source_settings['saved_nbands']
    require(expected == source_settings,
            'Saved ground state cannot supply changed SCF settings')
    actual = manifest(source / 'ground.gpw')
    require(actual['bytes'] == original['gpw']['bytes'] and
            actual['sha256'] == original['gpw']['sha256'],
            'Original saved ground bytes changed')
    require(original['neutral_electrons'] == config['valence_electrons'],
            'Saved PAW valence differs from registered host')
    target = output / 'ground.gpw'
    old_bands = original_config['ground_state']['saved_nbands']
    new_bands = config['ground_state']['saved_nbands']
    require(new_bands >= old_bands, 'Saved expansion cannot discard registered states')
    require(calc.get_number_of_bands() == old_bands,
            'Actual saved-state band count differs from original declaration')
    if new_bands != old_bands:
        calc.diagonalize_full_hamiltonian(nbands=new_bands)
        require(calc.get_number_of_bands() == new_bands,
                'Requested empty-state expansion did not occur')
        calc.write(str(target), mode='all')
        method = 'Full saved-potential Hamiltonian diagonalization only; no new SCF'
    else:
        shutil.copyfile(source / 'ground.gpw', target)
        method = 'Byte-identical saved-state copy; no new SCF or diagonalization'
    del calc
    ground = deepcopy(original)
    ground['ground_state_settings'] = config['ground_state']
    if 'magnetic_moment_muB' in original:
        require(abs(original['magnetic_moment_muB']) < 1e-10,
                'Original magnetic record contradicts the observed spin-paired state')
    ground['magnetic_moment_muB'] = 0.0
    ground['gpw'] = manifest(target)
    if new_bands == old_bands:
        require(ground['gpw']['bytes'] == actual['bytes'] and
                ground['gpw']['sha256'] == actual['sha256'],
                'Saved-state copy changed the authenticated ground bytes')
    ground['saved_state_reuse'] = dict(method=method, observed_spin_channels=1,
        original_input=manifest(source / 'input.json'),
        original_ground_record=manifest(source / 'ground.json'), original_gpw=actual)
    write_json(output / 'ground.json', ground)
    return ground
