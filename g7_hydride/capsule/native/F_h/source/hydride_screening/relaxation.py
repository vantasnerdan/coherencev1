"""Bounded standard-H relaxation and byte-bound completed-geometry reuse."""
import json
from pathlib import Path
import shutil

import numpy as np
from ase import Atoms
from ase.optimize import BFGS
from gpaw import GPAW

from geometry import build, geometry_record
from native_state import calculator, manifest, require, save_state, write_json


def relax(spec, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    atoms = build(spec)
    atoms.calc = calculator(atoms, spec, output/'standard-H.txt', relaxation=True)
    optimizer = BFGS(atoms, maxstep=.1, logfile=str(output/'optimizer.log'),
                     trajectory=str(output/'optimizer.traj'))
    latest = {}

    def record_boundary():
        forces = atoms.get_forces()
        stem = 'step-'+format(optimizer.nsteps, '03d')
        state = save_state(atoms, spec, 0, output, stem,
                           dict(forces_eV_A=forces.tolist(), ionic_step=int(optimizer.nsteps)))
        frame = geometry_record(atoms, spec, forces, optimizer.nsteps)
        frame.update(state=state['state_file'], gpw=state['gpw'])
        write_json(output/(stem+'-geometry.json'), frame)
        latest.clear()
        latest.update(frame)
        print('HYDRIDE_FORCE case='+spec['case']+' step='+str(optimizer.nsteps)+
              ' fmax='+str(frame['raw_fmax_eV_A']), flush=True)

    optimizer.attach(record_boundary, interval=1)
    converged = optimizer.run(fmax=spec['force_tolerance_eV_A'], steps=spec['ionic_steps'])
    write_json(output/'relaxation.json', dict(latest))
    require(bool(converged) and latest.get('passed') is True,
            'Standard-H geometry reached the finite force/step cap')
    atoms.calc = None
    return atoms, dict(latest)


def load_geometry(spec, source, output):
    source, output = Path(source), Path(output)
    frame = json.loads(source.read_bytes())
    require(frame['case'] == spec['case'] and frame['passed'] is True, 'Wrong or failed selected geometry')
    require(frame['volume_source'] == spec['volume_source'], 'Changed selected volume source')
    require(frame['raw_fmax_eV_A'] <= spec['force_tolerance_eV_A'], 'Selected force exceeds gate')
    checkpoint = source.parent/frame['gpw']['path']
    require(manifest(checkpoint) == frame['gpw'], 'Changed selected GPW bytes')
    state_path = source.parent/frame['state']['path']
    require(manifest(state_path) == frame['state'], 'Changed selected force state')
    original_state = json.loads(state_path.read_bytes())
    require(original_state['setups'] == {} and original_state['xc'] == 'PBE', 'Selected geometry is not standard-H PBE')
    calc = GPAW(str(checkpoint), txt=None)
    atoms = calc.get_atoms()
    require(atoms.get_chemical_symbols() == frame['symbols'], 'Selected GPW inventory differs')
    require(np.allclose(atoms.positions, frame['positions_A'], atol=1e-12, rtol=0) and
            np.allclose(atoms.cell.array, frame['cell_A'], atol=1e-12, rtol=0), 'Selected GPW geometry differs')
    forces = atoms.get_forces()
    actual_fmax = float(np.linalg.norm(forces, axis=1).max())
    require(abs(actual_fmax-frame['raw_fmax_eV_A']) <= 1e-5 and
            actual_fmax <= spec['force_tolerance_eV_A'], 'Selected GPW force is not admitted')
    atoms.calc = None
    shutil.copytree(source.parent, output)
    return atoms, frame
