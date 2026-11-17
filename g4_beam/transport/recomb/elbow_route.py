"""One continuous, explicit periodic D2-to-fcc route: an upper-bound path.

A low connected route can establish non-activation within the stated tolerance;
a positive sampled maximum cannot establish the global activation barrier.
"""
import numpy as np
from ase.geometry import find_mic
from elbow_geometry import elbow_pair


def uniform_two_cell(primitive):
    symbols = list(primitive.symbols)
    npd = symbols.count('Pd')
    if npd not in (3, 4) or symbols not in (['Pd']*npd, ['Pd']*npd+['H']):
        raise ValueError('Ordered clean or filled Pd3/Pd4 primitive required')
    fixed = set()
    for constraint in primitive.constraints:
        if constraint.todict()['name'] != 'FixAtoms':
            raise ValueError('Only the declared bottom-two-layer clamp supported')
        fixed.update(int(i) for i in constraint.get_indices())
    if fixed != {0, 1} or tuple(primitive.pbc) != (True, True, False):
        raise ValueError('Actual bottom-two clamp and xy periodicity required')
    repeated = primitive.repeat((1, 2, 1))
    count = len(primitive)
    order = [index for site in range(count) for index in (site, site+count)]
    repeated = repeated[order]
    repeated.calc = None
    if 'H' in symbols and not np.allclose(repeated.get_masses()[-2:], 2.014101778, rtol=0, atol=1e-12):
        raise ValueError('Explicit D isotope masses required')
    return repeated


def connected_knots(clean_primitive, filled_primitive, gas_bond_A):
    clean = uniform_two_cell(clean_primitive)
    filled = uniform_two_cell(filled_primitive)
    if len(clean) != 6 or list(filled.symbols) != ['Pd']*6+['H']*2:
        raise ValueError('This selected route is the three-layer1x2 model')
    if (not np.allclose(clean.cell, filled.cell, rtol=0, atol=1e-12)
            or not np.allclose(clean.positions[:4], filled.positions[:4], rtol=0, atol=1e-12)):
        raise ValueError('Source endpoints disagree on cell or clamped positions')
    if not .5 < gas_bond_A < 1.:
        raise ValueError('Actual molecular bond in angstrom required')
    charts = [(gas_bond_A, 6.), (gas_bond_A, 4.), (.85, 1.8),
              (1.2, 1.5), (1.3423676, 1.4525441), (1.8, 1.3),
              (float(np.linalg.norm(clean.cell[1]/2)), 1.1)]
    knots = [elbow_pair(clean, separation, height, 'top')
             for separation, height in charts]
    # Choose equivalent periodic images of the qualified filled endpoint nearest
    # the bridge chart. The interpolated final leg is then continuous in R^(3N).
    endpoint = filled.copy()
    displacement, _ = find_mic(filled.positions-knots[-1].positions,
                              filled.cell, pbc=filled.pbc)
    endpoint.positions[:] = knots[-1].positions+displacement
    middle = knots[-1].copy()
    middle.positions[:] = .5*(middle.positions+endpoint.positions)
    knots.extend((middle, endpoint))
    for knot in knots:
        knot.calc = None
        if not np.allclose(knot.positions[:4], clean.positions[:4], rtol=0, atol=1e-12):
            raise ValueError('Route moved a clamped source coordinate')
    return knots


def subdivide(knots):
    """Nested midpoint mesh of the SAME piecewise-linear Cartesian route."""
    points = []
    for left, right in zip(knots[:-1], knots[1:]):
        if (list(left.symbols) != list(right.symbols)
                or not np.array_equal(left.cell, right.cell)):
            raise ValueError('A continuous route needs one cell and atom ordering')
        points.append(left.copy())
        middle = left.copy()
        middle.positions[:] = .5*(left.positions+right.positions)
        middle.calc = None
        points.append(middle)
    points.append(knots[-1].copy())
    return points


def refine_selected_knots(frames, relaxed_ridge_positions_A):
    """The final six-frame correction; retain every unchanged actual pose."""
    if len(frames) != 17:
        raise ValueError('The admitted nine-knot/17-frame route is required')
    ridge = np.asarray(relaxed_ridge_positions_A, dtype=float)
    if (ridge.shape != frames[8].positions.shape or not np.all(np.isfinite(ridge))
            or not np.allclose(ridge[:4], frames[8].positions[:4], rtol=0, atol=1e-12)):
        raise ValueError('Saved relaxed ridge must preserve the source clamp')
    knots = [frame.copy() for frame in frames[::2]]
    bond = np.linalg.norm(frames[0].positions[-1]-frames[0].positions[-2])
    pair = knots[2].positions[-2:]
    center = pair.mean(axis=0)
    direction = (pair[1]-pair[0])/np.linalg.norm(pair[1]-pair[0])
    knots[2].positions[-2:] = center+np.array([-.5, .5])[:, None]*bond*direction
    knots[4].positions[4:] = ridge[4:]
    updated = subdivide(knots)
    result = [frame.copy() for frame in frames]
    changed = (3, 4, 5, 7, 8, 9)
    for index in changed:
        result[index].positions[:] = updated[index].positions
    for index, frame in enumerate(result):
        frame.calc = None
        frame.info['nested_coarse_point'] = index % 2 == 0
        if not np.array_equal(frame.positions[:4], frames[index].positions[:4]):
            raise ValueError('Final refinement moved a clamped atom')
    if any(not np.array_equal(result[i].positions, frames[i].positions)
           for i in range(17) if i not in changed):
        raise ValueError('An unselected source pose changed')
    return result, changed


def replace_clean_host(point, clean_host):
    """Map a selected chart point to a matched clean3L/4L top-layer frame."""
    npd = len(point)-2
    if (list(point.symbols) != ['Pd']*npd+['H']*2
            or len(clean_host) not in (6, 8) or set(clean_host.symbols) != {'Pd'}
            or not np.allclose(point.cell[:2], clean_host.cell[:2], rtol=0, atol=1e-12)
            or tuple(clean_host.pbc) != tuple(point.pbc)):
        raise ValueError('Matched clean host and two-D chart required')
    pair = point[-2:].copy()
    pair.positions[:] += clean_host.positions[-2]-point.positions[npd-2]
    result = clean_host.copy()
    result.calc = None
    result += pair
    return result
