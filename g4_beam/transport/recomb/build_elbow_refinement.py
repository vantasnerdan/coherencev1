"""Build the final six-frame correction from authenticated saved geometries."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import numpy as np
from ase.io import read, write
from elbow_route import refine_selected_knots
from surface_io import manifest


def bound_bytes(path, expected_sha256):
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Source byte identity changed')
    return raw


def build(args):
    bound_bytes(args.source_route, args.source_route_sha256)
    ridge = json.loads(bound_bytes(args.ridge_point, args.ridge_point_sha256))
    frames = read(args.source_route, index=':')
    geometry = ridge['geometry']
    if (list(frames[8].symbols) != geometry['symbols']
            or not np.array_equal(frames[8].cell, geometry['cell_A'])
            or not np.array_equal(frames[8].pbc, geometry['pbc'])):
        raise ValueError('Relaxed ridge belongs to a different atomic model')
    updated, changed = refine_selected_knots(frames, geometry['positions_A'])
    args.output.mkdir(parents=True, exist_ok=False)
    route = args.output/'route.traj'
    write(route, updated)
    saved = read(route, index=':')
    for index, (before, after) in enumerate(zip(updated, saved)):
        if (not np.array_equal(before.positions, after.positions)
                or not np.array_equal(before.get_masses(), after.get_masses())
                or not np.array_equal(before.cell, after.cell)
                or not np.array_equal(before.pbc, after.pbc)):
            raise ValueError('Saved refinement geometry changed')
        if index % 2 and not np.allclose(after.positions,
                .5*(saved[index-1].positions+saved[index+1].positions), rtol=0, atol=1e-14):
            raise ValueError('Refined route is not the same nested midpoint curve')
    result = dict(schema='coherence-final-elbow-refinement-v1',
        source_route=manifest(args.source_route), source_ridge_point=manifest(args.ridge_point),
        route=manifest(route), changed_indices=list(changed),
        unchanged_indices=[i for i in range(len(frames)) if i not in changed],
        inherited_energy_for_changed_frames=False, saved_roundtrip_exact=True,
        further_path_refinement_authorized=False,
        code=[manifest(Path(__file__)), manifest(Path(__file__).with_name('elbow_route.py'))])
    (args.output/'manifest.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-route', type=Path, required=True)
    parser.add_argument('--source-route-sha256', required=True)
    parser.add_argument('--ridge-point', type=Path, required=True)
    parser.add_argument('--ridge-point-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    build(parser.parse_args())


if __name__ == '__main__':
    main()
