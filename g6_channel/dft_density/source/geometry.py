"""Actual cubic PdD rows in the pinned channel consumer's (1-10,001,110) frame."""
import ast
import numpy as np

VECTORS = np.array([[1., -1., 0.], [0., 0., 1.], [1., 1., 0.]])


def periodic_distances(points, rows, lengths):
    delta = np.asarray(points)[:, None, :] - np.asarray(rows)[None, :, :]
    delta -= np.rint(delta / lengths) * lengths
    return np.linalg.norm(delta, axis=2).min(axis=1)


def construct(cell, positions, symbols, pd_radius):
    a = float(cell[0, 0])
    if a <= 0 or not np.allclose(cell, np.eye(3) * a, rtol=0, atol=1e-12):
        raise ValueError('Channel projection requires the admitted cubic cell')
    if list(symbols).count('Pd') != 4 or list(symbols).count('H') != 4:
        raise ValueError('Channel geometry is not pristine PdD')
    norms = np.linalg.norm(VECTORS, axis=1)
    lengths = a * norms
    unit = VECTORS / norms[:, None]
    translations = np.array(np.meshgrid(*[np.arange(-3, 4)] * 3)).reshape(3, -1).T
    rows = {}
    for symbol, label in [('Pd', 'Pd'), ('H', 'D')]:
        basis = positions[np.asarray(symbols) == symbol]
        points = (translations[:, None, :] * a + basis).reshape(-1, 3)
        coords = points @ unit.T
        selected = coords[np.all((coords >= -1e-9) & (coords < lengths - 1e-9), axis=1)]
        xy, counts = np.unique(np.round(selected[:, :2], 10), axis=0, return_counts=True)
        rows[label] = dict(xy_A=xy.tolist(), line_number_A_inv=(counts / lengths[2]).tolist())
        if not np.isclose(np.sum(counts) / np.prod(lengths), 4 / a**3, rtol=0, atol=1e-12):
            raise ValueError('String multiplicities do not preserve actual atomic density')
    all_rows = np.concatenate([np.asarray(row['xy_A']) for row in rows.values()])
    x, y = np.unique(all_rows[:, 0]), np.unique(all_rows[:, 1])
    if len(all_rows) != len(x) * len(y):
        raise ValueError('Expected rectangular pristine PdD string lattice')
    mx = (x + np.r_[x[1:], x[0] + lengths[0]]) / 2
    my = (y + np.r_[y[1:], y[0] + lengths[1]]) / 2
    centres = np.array(np.meshgrid(mx, my, indexing='ij')).reshape(2, -1).T
    pd_distance = periodic_distances(centres, rows['Pd']['xy_A'], lengths[:2])
    return dict(a_A=a, lengths_A=lengths.tolist(), unit_basis_rows=unit.tolist(), rows=rows,
                centres_xy_A=centres.tolist(), centres_nearest_Pd_A=pd_distance.tolist(),
                pd_augmentation_matching_radius_A=float(pd_radius),
                centres_outside_Pd_augmentation=(pd_distance > pd_radius).tolist(),
                centre_definition='All periodic rectangular-hole midpoints of the complete Pd+D string lattice')


def compare_pinned_geometry(geometry, source):
    """Execute only the pinned geometry and its three required constants, not the app."""
    tree = ast.parse(source)
    nodes = []
    for node in tree.body:
        if isinstance(node, ast.Assign):
            names = {child.id for target in node.targets for child in ast.walk(target) if isinstance(child, ast.Name)}
            if names & {'FCC', 'OCT', 'AXES'}:
                nodes.append(node)
        elif isinstance(node, ast.FunctionDef) and node.name == 'geometry':
            nodes.append(node)
    namespace = {'np': np}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), '<bound-channel-geometry>', 'exec'), namespace)
    lengths, rows = namespace['geometry']({'a': geometry['a_A'] / 10, 'x': 1}, '110')
    errors = [float(np.max(np.abs(lengths * 10 - geometry['lengths_A'])))]
    for label, (xy, line) in zip(('Pd', 'D'), rows):
        ours = geometry['rows'][label]
        if xy.shape != np.asarray(ours['xy_A']).shape:
            raise ValueError('Pinned consumer string multiplicity differs')
        errors.extend([float(np.max(np.abs(xy * 10 - ours['xy_A']))),
                       float(np.max(np.abs(line / 10 - ours['line_number_A_inv'])))])
    if max(errors) > 1.1e-9:
        raise ValueError('Actual string geometry differs from pinned consumer convention')
    return dict(max_absolute_geometry_difference=max(errors), tolerance=1.1e-9)


def outside_mask(geometry, grid):
    lengths = np.asarray(geometry['lengths_A'][:2])
    x, y = np.meshgrid(np.arange(grid) / grid * lengths[0], np.arange(grid) / grid * lengths[1], indexing='ij')
    points = np.column_stack((x.ravel(), y.ravel()))
    distance = periodic_distances(points, geometry['rows']['Pd']['xy_A'], lengths)
    return (distance > geometry['pd_augmentation_matching_radius_A']).reshape(grid, grid)
