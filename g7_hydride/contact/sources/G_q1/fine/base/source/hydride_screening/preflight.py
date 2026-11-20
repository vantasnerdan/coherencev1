"""Actual no-SCF geometry/setup/electron initialization before admission."""
import argparse
import gc
from pathlib import Path

from ase import Atom
import numpy as np

from geometry import build, DIRECTION, MASS_D
from native_state import calculator, require, write_json
from run_hydride import validate_epoch


def check(input_path, report):
    spec, epoch = validate_epoch(input_path)
    atoms = build(spec)
    cases = []
    for relaxation, charge in [(True, 0), (False, 0), (False, 1)]:
        current = atoms.copy()
        if not relaxation:
            target = current.positions[1 if spec['case'] == 'F' else 0]
            current += Atom('H', position=target-.1*DIRECTION, mass=MASS_D)
        calc = calculator(current, spec, Path(report).parent/('preflight-'+str(relaxation)+'-'+str(charge)+'.txt'),
                          charge, relaxation=relaxation)
        calc.initialize(current)
        expected = len(current)-charge
        require(abs(float(calc.get_number_of_electrons())-expected) < 1e-7, 'Wrong initialized electron count')
        for atom, setup in zip(current, calc.wfs.setups):
            require(abs(float(setup.Nv)-1) < 1e-12, 'Wrong actual Li/H valence')
            if atom.symbol == 'H':
                require((setup.ni > 0) == relaxation, 'Wrong actual geometry/BO H setup')
        cases.append(dict(relaxation=relaxation, charge=charge, electrons=expected,
                          setups=dict(calc.parameters.setups), grid_shape=calc.wfs.gd.N_c.tolist(),
                          kpoints=int(calc.wfs.kd.nibzkpts), nbands=int(calc.get_number_of_bands())))
        del calc, current
        gc.collect()
    require(np.isfinite(atoms.positions).all() and atoms.get_volume() > 0, 'Invalid physical geometry')
    result = dict(schema='hydride-native-preflight-v1', passed=True, specification=spec,
                  source_epoch=epoch, cases=cases, scope='Actual initialization only; no SCF/force/BO claim.')
    write_json(report, result)
    print('HYDRIDE_PREFLIGHT_COMPLETE case='+spec['case'], flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    check(args.input, args.report)
