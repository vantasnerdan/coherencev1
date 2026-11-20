"""Matched H/Li Native family; exact empty vacuum is explicit, never a fallback."""
import json
from pathlib import Path
import shutil
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from collision_ground import make_calculator, energy_record
from screen_io import manifest, require, write_json


def calculator(atoms, spec, log, charge=0, *, relaxation=False):
    require(charge in (0, 1), 'Only commissioned probe charges zero and one')
    require(set(atoms.get_chemical_symbols()) <= {'H', 'Li'}, 'Unexpected hydride species')
    family = dict(h_A=spec['h_A'], repeat=[spec['repeat_primitive']]*3, kmesh=spec['kmesh'])
    base = make_calculator(atoms, family, log, charge=charge)
    electrons = len(atoms)-charge  # Actual default Li.PBE and H/Hae both have Nv=1.
    require(electrons >= 0, 'Negative electron inventory')
    experimental = dict(base.parameters.experimental)
    experimental['reuse_wfs_method'] = 'paw' if relaxation else 'keep'
    return base.new(nbands=max(4, int(np.ceil(electrons/2))+spec['nband_padding']),
                    setups={} if relaxation else {'H': 'ae'},
                    experimental=experimental, txt=str(log))


def describe(atoms, spec, charge):
    calc = atoms.calc
    electrons = float(calc.get_number_of_electrons())
    require(abs(electrons-(len(atoms)-charge)) < 1e-7, 'Wrong H/Li electron inventory')
    require(bool(calc.scf.converged), 'Native SCF did not converge')
    valence = {}
    for atom, setup in zip(atoms, calc.wfs.setups):
        require(abs(float(setup.Nv)-1) < 1e-12, 'Unexpected PAW valence')
        valence[atom.symbol] = float(setup.Nv)
        if atom.symbol == 'H':
            expected_ae = calc.parameters.setups.get('H') == 'ae'
            require((setup.ni == 0) == expected_ae, 'Wrong actual H setup family')
    count = calc.get_number_of_iterations()
    return dict(schema='hydride-native-state-v1', backend='gpaw-native',
                energy=energy_record(calc), charge=charge, electrons=electrons,
                symbols=atoms.get_chemical_symbols(), positions_A=atoms.positions.tolist(),
                cell_A=atoms.cell.array.tolist(), masses_u=atoms.get_masses().tolist(),
                realspace_shape=calc.wfs.gd.N_c.tolist(), actual_kpoints=int(calc.wfs.kd.nibzkpts),
                nbands=int(calc.get_number_of_bands()), specification=spec, scf_converged=True,
                xc=calc.parameters.xc, setups=dict(calc.parameters.setups),
                spinpol=bool(calc.parameters.spinpol), occupation_width_eV=.1,
                setup_valence=valence, wavefunction_predictor=type(calc.wfs.wfs_mover).__name__,
                scf_iterations=None if count is None else int(count))


def save_state(atoms, spec, charge, output, stem, extra=None):
    output = Path(output)
    temporary = output/(stem+'.writing.gpw')
    atoms.calc.write(str(temporary), mode='all')
    checkpoint = output/(stem+'.gpw')
    temporary.replace(checkpoint)
    record = describe(atoms, spec, charge)
    record['gpw'] = manifest(checkpoint)
    if extra:
        record.update(extra)
    path = output/(stem+'.json')
    write_json(path, record)
    return record | {'state_file': manifest(path)}


def reuse_state(atoms, spec, charge, output, stem, source, extra):
    source = Path(source)
    record = json.loads(source.read_bytes())
    require(record['backend'] == 'gpaw-native' and record['scf_converged'], 'Cannot reuse non-Native/failed state')
    require(record['specification'] == spec and record['charge'] == charge, 'Changed reused physical input')
    require(record['setups'] == {'H': 'ae'}, 'Cannot reuse PAW geometry energy in quartet')
    require(record['symbols'] == atoms.get_chemical_symbols(), 'Changed reused state inventory')
    require(np.allclose(record['positions_A'], atoms.positions, rtol=0, atol=1e-12) and
            np.allclose(record['cell_A'], atoms.cell.array, rtol=0, atol=1e-12), 'Changed reused geometry')
    for key, value in (extra or {}).items():
        require(record[key] == value, 'Changed reused state annotation')
    checkpoint = source.parent/record['gpw']['path']
    require(manifest(checkpoint) == record['gpw'], 'Changed reused GPW')
    require(record['gpw']['path'] == stem+'.gpw', 'Reused state stem mismatch')
    shutil.copyfile(checkpoint, output/(stem+'.gpw'))
    shutil.copyfile(source, output/(stem+'.json'))
    return record | {'state_file': manifest(output/(stem+'.json'))}


def scf(atoms, spec, charge, output, stem, extra=None, reuse=None):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if reuse is not None:
        return reuse_state(atoms, spec, charge, output, stem, reuse, extra)
    if len(atoms) == 0:
        require(spec['case'] == 'G' and charge == 0 and stem == 'vacant-host',
                'Empty reference is only the exact atomic-G vacant host')
        record = dict(schema='hydride-native-state-v1', backend='analytic-empty-vacuum',
                      energy=dict(internal_energy_eV=0., free_energy_eV=0., minus_TS_eV=0.),
                      charge=0, electrons=0., symbols=[], positions_A=[], cell_A=atoms.cell.array.tolist(),
                      realspace_shape=None, actual_kpoints=None, nbands=0, specification=spec,
                      scf_converged=False, xc='PBE', setups={'H': 'ae'}, spinpol=False,
                      occupation_width_eV=.1, setup_valence={},
                      reference='Zero nuclei/electrons: kinetic, Hartree, XC, ionic and atomic-reference sums are empty.')
        write_json(output/(stem+'.json'), record)
        return record | {'state_file': manifest(output/(stem+'.json'))}
    started = time.process_time()
    if atoms.calc is None:
        atoms.calc = calculator(atoms, spec, output/(stem+'.txt'), charge)
    print('HYDRIDE_SCF_START case='+spec['case']+' state='+stem, flush=True)
    atoms.get_potential_energy()
    record = save_state(atoms, spec, charge, output, stem,
                        (extra or {}) | {'process_cpu_s': time.process_time()-started})
    print('HYDRIDE_SCF_COMPLETE case='+spec['case']+' state='+stem+
          ' cpu_s='+str(record['process_cpu_s']), flush=True)
    return record
