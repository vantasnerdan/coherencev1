"""Consume complete Native phonons only after authenticating all actual GPWs."""
import json
from pathlib import Path
import numpy as np
from ase.io import read
from ase.utils.filecache import MultiFileJSONCache
from force_bank_results import verify_force_record
from phonon_core import RecordedPhonons, SUPERCELL, consume, layout
from surface_io import manifest


def _bound_json(bound):
    path = Path(bound['path'])
    if manifest(path) != bound:
        raise ValueError('Phonon saved-output byte binding changed')
    return json.loads(path.read_bytes())


def read_native_spectrum(directory, *, primitive_energy_eV,
                         mobile_fmax_target=.01, folding_energy_tolerance=.001):
    directory = Path(directory)
    saved = json.loads((directory / 'hessian.json').read_bytes())
    inputs = _bound_json(saved['inputs'])
    if not saved['accepted']:
        raise ValueError('Native withholding did not admit a spectrum')
    source = inputs['primitive']
    if manifest(Path(source['path'])) != source:
        raise ValueError('Settled primitive binding changed')
    primitive = read(source['path'])
    active, full = layout(primitive)
    for bound in saved['force_cache']:
        _bound_json(bound)
    caches = {}
    for name in ('cache', 'withheld-cache'):
        if not (directory / name).is_dir():
            raise ValueError('A complete force cache is missing')
        cache = MultiFileJSONCache(directory / name)
        for record in cache.values():
            verify_force_record(record)
        caches[name] = cache
    expected = {'eq'} | {f'{i}{a}{s}' for i in active
                         for a in 'xyz' for s in ('-', '+')}
    withheld_expected = {f'{i}x{s}' for i in full[len(active):]
                          for s in ('-', '+')}
    if (set(caches['cache']) != expected
            or set(caches['withheld-cache']) != withheld_expected):
        raise ValueError('Native displacement coverage differs from declared full-cell reduction')
    ph = RecordedPhonons(primitive, None, supercell=SUPERCELL,
                         delta=inputs['delta_A'], name=str(directory / 'cache'))
    ph.set_atoms(active)
    actual = consume(ph, caches['withheld-cache'])
    if not actual['accepted']:
        raise ValueError('Actual translated-cell forces refuse this reconstruction')
    for key in ('raw_C_N_eV_A2',):
        if not np.allclose(actual[key], saved[key], rtol=0, atol=1e-12):
            raise ValueError('Saved raw force constants differ from actual GPW forces')
    if not np.allclose(actual['spectrum']['signed_quanta_eV'],
                       saved['spectrum']['signed_quanta_eV'], rtol=0, atol=1e-12):
        raise ValueError('Saved signed modes differ from actual mass-weighted force constants')
    center = caches['cache']['eq']
    fmax = float(np.linalg.norm(np.asarray(center['forces'])[full], axis=1).max())
    fold = float(center['energy_eV'] - 2 * primitive_energy_eV)
    if fmax > mobile_fmax_target or abs(fold) > folding_energy_tolerance:
        raise ValueError('Repeated-cell center fails stationarity or primitive energy folding')
    actual.update(center_energy_eV=center['energy_eV'], center_mobile_fmax_eV_A=fmax,
                  two_cell_energy_folding_difference_eV=fold,
                  saved_result=manifest(directory / 'hessian.json'))
    return actual
