"""Fixed spectator/endpoint ions for an explicitly signed neutral electronic cycle."""
import itertools
import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell
from ase.geometry import find_mic
from screen_io import require


def geometry(atoms):
    return dict(symbols=atoms.get_chemical_symbols(), positions_A=atoms.positions.tolist(),
                cell_A=atoms.cell.array.tolist(), pbc=atoms.pbc.tolist(),
                initial_moments_muB=atoms.get_initial_magnetic_moments().tolist())


def host_cases(row):
    material = row['material_id']
    seed = row['geometry']
    if material == 'Pd_hydro_+0pct':
        a = np.linalg.norm(seed['cell_A'][0])*np.sqrt(2)
        host = bulk('Pd', 'fcc', a=a, cubic=True)
    else:
        require(material == 'Pd3Rh_L12', 'Wrong neutral-cycle host')
        host = Atoms(seed['symbols'], positions=seed['positions_A'],
                     cell=seed['cell_A'], pbc=True)
        a = host.cell[0, 0]
    require(len(host) == 4 and np.allclose(host.cell.array, np.eye(3)*a),
            'Expected the actual cubic four-metal host')
    configurations = []
    for size, transform in [(4, np.eye(3, dtype=int)),
                             (8, np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))]:
        spectators = make_supercell(host, transform)
        first = np.array([.5, .5, .5])*a
        second = np.array([.5, 0., 0.])*a
        distance = float(find_mic(first-second, spectators.cell, pbc=True)[1])
        shortest = min(np.linalg.norm(np.array(v)@spectators.cell.array)
                       for v in itertools.product(range(-2, 3), repeat=3) if any(v))
        for endpoint in ['DDfar', 'He']:
            atoms = spectators.copy()
            if endpoint == 'DDfar':
                atoms += Atoms('HH', positions=[first, second])
            else:
                atoms += Atoms('He', positions=[first])
            electrons = sum({'Pd':10, 'Rh':9, 'H':1, 'He':2}[s] for s in atoms.symbols)
            configurations.append(dict(case=f'{material}-{size}-{endpoint}',
                host=material, endpoint=endpoint, metal_atoms=size, geometry=geometry(atoms),
                valence_electrons=electrons, kmesh=[4]*3, spinpol=False,
                minimum_host_image_A=float(shortest), DD_nearest_image_A=distance,
                site_definition='A body octahedron (six Pd); B edge octahedron (six Pd in Pd, four Pd/two Rh in L12). Same A/B sites in both cells; He at A.',
                expected_wall_s=120 if size==4 else 300))
    return configurations


def vacuum_cases():
    cases = []
    for length in [8., 10.]:
        for symbol, spin, electrons in [('H', 1., 1), ('He', 0., 2)]:
            atoms = Atoms(symbol, positions=[[length/2]*3], cell=[length]*3, pbc=True)
            atoms.set_initial_magnetic_moments([spin])
            cases.append(dict(case=f'vacuum-{symbol}-{int(length)}A', host='vacuum',
                endpoint=symbol, metal_atoms=0, geometry=geometry(atoms),
                valence_electrons=electrons, kmesh=[1]*3, spinpol=bool(spin),
                expected_wall_s=180, vacuum_box_A=length))
    return cases
