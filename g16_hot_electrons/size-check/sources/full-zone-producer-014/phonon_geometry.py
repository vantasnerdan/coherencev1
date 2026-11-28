"""One real-space Pd/D displacement basis and exact primitive-mesh folds."""
import numpy as np
from ase import Atoms
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

from native import reference, require

TEMPERATURES = (300, 3000, 10000, 30000, 1000, 1500, 2000, 310, 350)
CELLS = {'primitive2': 2*np.eye(3, dtype=int), 'conventional2': 2*(np.ones((3,3), dtype=int)-2*np.eye(3,dtype=int))}
STEP_A = .01


def phonon(cell):
    atoms = reference()
    unit = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.cell.array,
                        positions=atoms.positions, masses=atoms.get_masses())
    ph = Phonopy(unit, supercell_matrix=CELLS[cell].T, primitive_matrix=np.eye(3))
    ph.generate_displacements(distance=STEP_A, is_diagonal=False, is_plusminus=True)
    if cell == 'conventional2':
        # Phonopy chooses [100] in the cubic box; preserve the original
        # Cartesian [011] finite step so size is the only changed variable.
        direction = np.asarray(phonon('primitive2').displacements[0][1:])
        data = ph.dataset
        for row in data['first_atoms']:
            sign = np.sign(row['displacement'][0])
            require(sign != 0, 'Unregistered cubic displacement direction')
            row['displacement'] = sign * direction
        ph.dataset = data
    require(len(ph.displacements)==4 and len(ph.primitive)==2, 'Unexpected Pd/D site symmetry or displacement count')
    return ph


def ase_atoms(cell):
    return Atoms(cell.symbols, positions=cell.positions, cell=cell.cell, masses=cell.masses, pbc=True)


def inputs(ph):
    ideal = ase_atoms(ph.supercell)
    cases = {'zero': (ideal, None, None)}
    for (index, *u), cell in zip(ph.displacements, ph.supercells_with_displacements, strict=True):
        step = np.asarray(u)
        sign = 'm' if step[np.flatnonzero(np.abs(step)>1e-12)[0]] < 0 else 'p'
        element = {'Pd':'pd','H':'d'}[ideal.symbols[index]]
        key = f'{element}-{sign}'
        require(key not in cases, 'Phonopy generated more than one direction at a cubic site')
        cases[key] = (ase_atoms(cell), int(index), step)
    require(set(cases)=={'zero','pd-p','pd-m','d-p','d-m'}, 'Both signed species are required')
    return cases


def electronic_mesh(n, cell):
    require(n in (12,16), 'Use a registered primitive electronic mesh')
    s = CELLS[cell]
    numerator = np.indices((n,)*3).reshape(3,-1).T
    fold = (numerator@s.T + n//2)%n - n//2
    points, count = np.unique(fold,axis=0,return_counts=True)
    require(len(points)*abs(round(np.linalg.det(s)))==n**3 and np.all(count==count[0]), 'Unequal electronic folding preimages')
    return points / n


def inversion(ideal, mover):
    """Native force transport from the positive structure to its negative partner."""
    x = ideal.get_scaled_positions(wrap=False)
    delta = (2*x[mover]-x)[:,None,:] - x[None,:,:]
    delta -= delta.round()
    distances = np.linalg.norm(delta @ ideal.cell.array, axis=2)
    ids = distances.argmin(axis=1)
    require(np.array_equal(ids[ids],np.arange(len(x))) and distances[np.arange(len(x)),ids].max()<1e-9
            and np.array_equal(ideal.numbers[ids],ideal.numbers), 'Inversion is not a periodic same-species permutation')
    return ids
