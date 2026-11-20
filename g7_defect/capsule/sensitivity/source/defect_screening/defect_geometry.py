"""Fixed-volume A/C/D, real ionic positions, explicit vacancy occupancy and ray."""
import numpy as np
from ase import Atom, Atoms
from ase.build import bulk
from ase.geometry import find_mic

MASS_D = 2.01410177812


def build(spec):
    case = spec['case']
    if case not in ('A', 'C', 'D'):
        raise ValueError('Only new A/C/D geometries; B must be reused')
    n, a = int(spec['repeat_primitive']), float(spec['lattice_A'])
    if n < 2 or a <= 0:
        raise ValueError('Invalid admitted fcc repeat/lattice')
    ordered = bulk('PdH', 'rocksalt', a=a).repeat((n, n, n))
    pd = np.flatnonzero(ordered.numbers == 46)
    hydrogen = np.flatnonzero(ordered.numbers == 1)
    missing = ordered.positions[pd[0]].copy()
    seeded = []
    if case == 'A':
        atoms = Atoms('Pd' * len(pd), positions=ordered.positions[pd], cell=ordered.cell, pbc=True)
        atoms += Atom('H', position=ordered.positions[hydrogen[0]], tag=1, mass=MASS_D)
        centre = ordered.positions[hydrogen[0]].copy()
    elif case == 'C':
        atoms = Atoms('Pd' * (len(pd)-1), positions=ordered.positions[pd[1:]], cell=ordered.cell, pbc=True)
        atoms += Atom('H', position=missing, tag=1, mass=MASS_D)
        centre = missing.copy()
    else:
        atoms = ordered.copy()
        atoms.set_tags(np.zeros(len(atoms), int))
        vectors, distances = find_mic(atoms.positions[hydrogen]-missing, atoms.cell, True)
        target = int(hydrogen[np.argmin(distances)])
        atoms.positions[target] = missing
        tags = atoms.get_tags(); tags[target] = 1; atoms.set_tags(tags)
        for local in np.argsort(distances):
            index = int(hydrogen[local])
            if index == target or distances[local] > a*.51:
                continue
            atoms.positions[index] = missing + vectors[local] * (1.2 / distances[local])
            seeded.append(dict(index_before_Pd_removal=index, original_distance_A=float(distances[local]), seed_radius_A=1.2))
        del atoms[int(pd[0])]
        centre = missing.copy()
    shift = .5*np.sum(atoms.cell.array, axis=0)-centre
    atoms.translate(shift)
    vacancy = missing+shift if case != 'A' else None
    # Break artificial Pd symmetry, never pin Pd or force a target to stay trapped.
    indices = np.flatnonzero(atoms.numbers == 46)
    displacement = .02*np.sin((np.arange(len(indices))[:,None]+1)*(np.arange(3)[None,:]+1)*np.sqrt(2.))
    displacement -= displacement.mean(axis=0)
    atoms.positions[indices] += displacement
    atoms.wrap()
    explicit = spec.get('initial_positions_A')
    if explicit is not None:
        positions = np.asarray(explicit, dtype=float)
        if case not in ('C', 'D') or positions.shape != atoms.positions.shape or not np.isfinite(positions).all():
            raise ValueError('C/D starting positions must preserve the finite atomic inventory')
        atoms.positions[:] = positions
        vacancy = np.asarray(spec['vacancy_initial_A'], dtype=float)
        if vacancy.shape != (3,) or not np.isfinite(vacancy).all():
            raise ValueError('Explicit geometry requires its finite registered vacancy centre')
    atoms.set_masses([MASS_D if atom.symbol == 'H' else atom.mass for atom in atoms])
    target_index = int(np.flatnonzero(atoms.get_tags() == 1)[0])
    return atoms, dict(case=case, target_index=target_index, original_Pd_sites=n**3,
        probe_direction=spec.get('probe_direction'),
        seed_source=spec.get('seed_source'),
        vacancy_initial_A=None if vacancy is None else vacancy.tolist(), seeded_neighbour_D=seeded,
        initial_positions_A=atoms.positions.tolist(), symbols=atoms.get_chemical_symbols(), cell_A=atoms.cell.array.tolist(),
        constraint='No fixed ions; raw all-ion forces. Fixed lattice only.',
        seed=('Explicit source-bound starting positions; all ions remain free.' if explicit is not None else
              'Deterministic zero-mean Pd displacement0.02A; D inventory fixed, D target initially at vacancy for C/D.'))


def vacancy_record(atoms, initial):
    if initial['vacancy_initial_A'] is None:
        return None
    pd = np.asarray(initial['symbols']) == 'Pd'
    vectors, _ = find_mic(atoms.positions[pd]-np.asarray(initial['initial_positions_A'])[pd], atoms.cell, True)
    drift = vectors.mean(axis=0)
    centre = np.asarray(initial['vacancy_initial_A'])+drift
    h = np.flatnonzero(atoms.numbers == 1)
    _, distances = find_mic(atoms.positions[h]-centre, atoms.cell, True)
    target_vector, target_distance = find_mic(atoms.positions[initial['target_index']]-centre, atoms.cell, True)
    return dict(registered_centre_A=centre.tolist(), mean_Pd_translation_A=drift.tolist(),
        hydrogen_indices=h.tolist(), hydrogen_distances_A=distances.tolist(),
        occupancy_radius_A=1.5, trapped_n=int(np.sum(distances < 1.5)),
        occupancy_sensitivity={str(radius): int(np.sum(distances < radius)) for radius in (1.3,1.5,1.7)},
        target_offset_A=target_vector.tolist(), target_distance_A=float(target_distance),
        scope='Geometric occupancy of this relaxed fixed-inventory local basin, not a saturation capacity.')


def choose_path(atoms, initial, vacancy):
    target = atoms.positions[initial['target_index']]
    if initial['case'] == 'A':
        return dict(target_A=target.tolist(), direction=(np.array([1.,1.,0.])/np.sqrt(2)).tolist(),
                    scope='Same Cartesian[110] approach as the reused B baseline')
    crystal = np.array(np.meshgrid(*[np.arange(-1,2)]*3)).reshape(3,-1).T
    crystal = crystal[np.any(crystal != 0, axis=1)]
    crystal = crystal/np.linalg.norm(crystal,axis=1)[:,None]
    z = 1-2*(np.arange(128)+.5)/128
    angle = np.arange(128)*np.pi*(3-np.sqrt(5))
    sphere = np.column_stack((np.sqrt(1-z*z)*np.cos(angle),np.sqrt(1-z*z)*np.sin(angle),z))
    candidates = np.concatenate((crystal,sphere))
    fixed = initial.get('probe_direction')
    if fixed is not None:
        direction = np.asarray(fixed, dtype=float)
        if direction.shape != (3,) or not np.isfinite(direction).all() or abs(np.linalg.norm(direction)-1) > 1e-12:
            raise ValueError('Numerical control requires the frozen unit probe direction')
        candidates = direction[None,:]
    offset = np.asarray(vacancy['target_offset_A'])
    if fixed is None and np.linalg.norm(offset) > 1e-3:
        candidates = candidates[(candidates@offset) > 0]
    background = np.delete(atoms.positions,initial['target_index'],axis=0)
    radii = np.linspace(.1,1.6,61)
    best = None
    for direction in candidates:
        points = target-radii[:,None]*direction
        _, distance = find_mic((points[:,None,:]-background).reshape(-1,3),atoms.cell,True)
        clearance = float(distance.min())
        if best is None or clearance > best[0]:
            best = (clearance,direction)
    return dict(target_A=target.tolist(),direction=best[1].tolist(),
                minimum_background_clearance_A=best[0],candidate_directions=len(candidates),
                sampled_radii_A=radii.tolist(),vacancy=vacancy,
                matched_control_direction=fixed is not None,
                vacancy_side_projection_A=float(best[1]@offset),
                scope=('Frozen base-C Cartesian ray for the cell comparison; no independently optimized path.' if fixed is not None else
                       'Best sampled head-on ray by background-nucleus clearance; vacancy-side hemisphere when target is displaced, all directions at the vacancy centre. Not a global angular optimum.'))
