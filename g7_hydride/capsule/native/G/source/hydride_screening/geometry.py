"""Source-bound LiD repeats and isolated atomic/molecular D geometry."""
import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from native_state import require

MASS_D = 2.01410177812
MASS_LI6 = 6.0151228874
RADII = [.1, .2, .3, .5, 1., 1.6]
DIRECTION = np.ones(3)/np.sqrt(3.)


def build(spec):
    case = spec['case']
    if case == 'F':
        require(spec['primitive_symbols'] == ['Li', 'H'], 'Wrong LiD primitive inventory')
        n = spec['repeat_primitive']
        require(n in (2, 3), 'Only commissioned finite LiD cells')
        atoms = Atoms(spec['primitive_symbols'], positions=spec['primitive_positions_A'],
                      cell=spec['primitive_cell_A'], pbc=True, masses=[MASS_LI6, MASS_D]).repeat((n, n, n))
        require(atoms.get_chemical_symbols().count('Li') == n**3 and
                atoms.get_chemical_symbols().count('H') == n**3, 'Wrong repeated inventory')
    elif case in ('G', 'G2'):
        length = float(spec['box_A'])
        require(length >= 10, 'Gas box smaller than commissioned large-box reference')
        center = np.full(3, length/2)
        positions = [center] if case == 'G' else [center-.37*DIRECTION, center+.37*DIRECTION]
        atoms = Atoms('H' if case == 'G' else 'H2', positions=positions,
                      cell=np.eye(3)*length, pbc=True, masses=[MASS_D]*len(positions))
    else:
        raise ValueError('Uncommissioned hydride/gas case')
    require(not atoms.constraints, 'Host ions may not be constrained')
    return atoms


def geometry_record(atoms, spec, forces, steps):
    target_index = 1 if spec['case'] == 'F' else 0
    target = atoms.positions[target_index]
    background = np.delete(atoms.positions, target_index, axis=0)
    clearance = None
    if len(background):
        probes = target-np.linspace(.1, 1.6, 301)[:, None]*DIRECTION
        vectors = (probes[:, None, :]-background[None, :, :]).reshape(-1, 3)
        _, distances = find_mic(vectors, atoms.cell, pbc=True)
        clearance = float(np.min(distances))
    fmax = float(np.linalg.norm(forces, axis=1).max())
    return dict(case=spec['case'], passed=fmax <= spec['force_tolerance_eV_A'],
                symbols=atoms.get_chemical_symbols(), positions_A=atoms.positions.tolist(),
                masses_u=atoms.get_masses().tolist(), cell_A=atoms.cell.array.tolist(),
                target_index=target_index, forces_eV_A=np.asarray(forces).tolist(),
                raw_fmax_eV_A=fmax, steps=int(steps), relaxation_h_A=spec['h_A'],
                volume_source=spec['volume_source'],
                path=dict(target_A=target.tolist(), direction=DIRECTION.tolist(),
                          minimum_background_clearance_A=clearance,
                          description='Head-on [111]; open tetrahedral side for F, outside the bond for G2.'))
