"""Endpoint stages, runnable only after measured pilot allocation decision."""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import json
from pathlib import Path
import numpy as np
from ase.constraints import FixSubsetCom
from ase.io import read
from gpaw import GPAW
from gpaw.mpi import serial_comm
from surface_geometry import slab, gas
from surface_io import save, manifest
from surface_optimize import optimize


def main():
    root = Path(__file__).resolve().parent
    decision = json.loads((root/'production-decision.json').read_text())
    if decision['decision'] != 'continue' or decision['measured_pilot_binding'] is None:
        raise ValueError('Measured pilot continuation decision required')
    pilot = json.loads((root/'pilot/pilot.json').read_text())
    last = json.loads(Path(pilot['points'][-1]['path']).read_text())
    calc = GPAW(str(root/'pilot/point-002.gpw'), legacy_gpaw=False,
                communicator=serial_comm, txt=str(root/'restart-proof-scf.txt'))
    atoms = read(root/'pilot/relax.traj', index=-1)
    atoms.calc = calc
    # Actual saved native state is read, then its concrete primitive is checked.
    e = float(atoms.get_potential_energy())
    f = atoms.get_forces()
    if abs(e-last['energy_extrapolated_eV']) > 1e-8 or not np.allclose(f, last['forces_eV_A'], atol=1e-8):
        raise ValueError('Saved native state does not reproduce completed pilot point')
    save(root/'restart-proof.json', {'reused_completed_force_point': manifest(root/'pilot/point-002.json'),
         'gpw': manifest(root/'pilot/point-002.gpw'), 'energy_eV': e,
         'max_force_difference_eV_A': float(np.max(np.abs(f-np.array(last['forces_eV_A'])))),
         'claim': 'Completed calculator state reuse only; no FIRE/SCF in-flight-history claim'})
    adsorbed = optimize(atoms, root/'adsorbed-1ML', restart_calc=calc)
    clean = optimize(slab(hydrogens=0), root/'clean')
    optimize(gas('H'), root/'H-atom', steps=0)
    optimize(gas('H2'), root/'H2-gas', steps=20)
    final = adsorbed.copy()
    # Clean endpoint retains matched cell and atom ordering. Molecular COM is an
    # explicit gas-reservoir coordinate, not a claimed unconstrained translation minimum.
    final.positions[:-2] = clean.positions
    vector = adsorbed.positions[-1]-adsorbed.positions[-2]
    vector[2] = 0
    vector /= np.linalg.norm(vector)
    centre = adsorbed.positions[-2:].mean(axis=0)
    centre[2] = clean.positions[:, 2].max()+5.0
    final.positions[-2:] = [centre-0.375*vector, centre+0.375*vector]
    final.set_constraint([*clean.constraints, FixSubsetCom(indices=[len(final)-2, len(final)-1])])
    optimize(final, root/'gas-slab-endpoint', steps=25)
    save(root/'endpoints-completed.json', {'stages': ['adsorbed-1ML', 'clean', 'H-atom', 'H2-gas', 'gas-slab-endpoint'],
         'stage_results': [manifest(root/stage/'result.json') for stage in
                           ['adsorbed-1ML', 'clean', 'H-atom', 'H2-gas', 'gas-slab-endpoint']]})
    print('ENDPOINTS_COMPLETED', flush=True)

if __name__ == '__main__':
    main()
