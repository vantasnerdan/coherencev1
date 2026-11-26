"""Matched frozen D26/D28 points on the admitted covered PdD(111) slab."""
import argparse
import json
from pathlib import Path
import runpy

import numpy as np
from ase.geometry import find_mic

from surface_vacancy import (DEFAULT_ROOT, MANIFEST_PATH, calculator, file_identity,
                             load_surface_checkpoint, native, profile, publish, require,
                             root_action, slab_geometry, source_files)

REF_SHA = '7aadd280b7a58106a79702832d27c33639d04a2a5e14899065abe6b879a9c6a9'
ZERO_SHA = '526ec8d8d70a61ea4cb2d20f907ed757bb121050c598539a3ec3b0ad365fc33b'
I, J = 26, 28
AMPS = (6, 12, 24, 36, 48)
MESH_CHECK = 36
SCHEMA = 'g10-covered-slab-pair-v1'


def geometry(atoms, role, amp):
    require((role == 'zero' and amp == 0) or (role in ('A', 'B', 'AB') and amp in AMPS), 'Unregistered slab-pair case')
    r = atoms.positions
    raw = r[J] - r[I]
    v, distance = find_mic(raw, atoms.cell, atoms.pbc)
    require(atoms.numbers[[I, J]].tolist() == [1, 1] and abs(distance-2.437039831543227) < 1e-12,
            'Wrong original covered-slab pair')
    shift = np.rint(np.linalg.solve(atoms.cell.array.T, v-raw)).astype(int)
    axis = v/distance
    u = amp/100
    moved = atoms.copy()
    if role in ('A', 'AB'):
        moved.positions[I] += u*axis
    if role in ('B', 'AB'):
        moved.positions[J] -= u*axis
    pair_distance = float(np.linalg.norm(moved.positions[J] + shift @ atoms.cell.array - moved.positions[I]))
    movers = np.array([I, J]); metal = np.flatnonzero(atoms.numbers == 46)
    clearance = moved.get_distances(I, metal, mic=True).tolist(), moved.get_distances(J, metal, mic=True).tolist()
    minimum = min(min(clearance[0]), min(clearance[1]))
    return moved, dict(role=role, displacement_per_D_A=u, original_D_indices=movers.tolist(),
                       axis_A_to_B=axis.tolist(), B_image_cell_translation=shift.tolist(), separation_A=pair_distance,
                       min_tagged_Pd_D_A=minimum, nominal_Pd_H_overlap=bool(minimum < 3.4*0.529177210903))


def run(args, site, atoms, writer, binding):
    import gpaw
    from gpaw.mpi import world
    from loading_pd import bind_ranks

    bind_ranks(world, args.ranks)
    gpaw.setup_paths[:] = [str(Path(site.sources['H_PAW']['path']).parent)]
    moved, path = geometry(atoms, args.role, args.amp)
    cfg = profile(site, 'base')
    if args.mesh == 3:
        require((args.role, args.amp) in (('zero', 0), ('AB', MESH_CHECK)), 'Unregistered fine-mesh location')
        cfg['kmesh'] = [3, 3, 1]
    out = args.output
    log = out.with_suffix('.scf')

    def prepare():
        out.parent.mkdir(parents=True, exist_ok=True)
        require(not out.exists() and not log.exists(), 'Existing attempt must be retained separately')

    root_action(world, prepare)
    moved.calc = calculator(cfg, log, world)
    if world.rank == 0:
        print(f'G10_SLAB_PAIR_START {args.mesh} {args.role}-{args.amp}', flush=True)
    result = native(moved, cfg)
    if args.mesh == 2 and args.role == 'zero':
        old = DEFAULT_ROOT/'coherence-g10-wide/sites/points/surface/zero.json'
        require(file_identity(old)['sha256'] == ZERO_SHA, 'Original coarse zero differs')
        anchor = json.loads(old.read_bytes())
        require(abs(anchor['energy_free_eV']-result['energy_free_eV']) <= 1e-4
                and np.max(np.abs(np.asarray(anchor['forces_eV_A'])-result['forces_eV_A'])) <= 5e-4,
                'Independent zero does not reproduce the admitted coarse zero')
    if world.rank == 0:
        setup_hashes = {k: site.sources[k]['sha256'] for k in ('Pd_PAW', 'H_PAW')}
        result.update(schema=SCHEMA, point=path, reference=binding, setup_sha256=setup_hashes, scf_log=log.name)
    publish(world, writer, out, result)
    if world.rank == 0:
        print(f'G10_SLAB_PAIR_DONE {args.mesh} {args.role}-{args.amp}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    sub = p.add_subparsers(dest='action', required=True)
    sub.add_parser('plan')
    calc = sub.add_parser('point')
    calc.add_argument('role', choices=('zero', 'A', 'B', 'AB'))
    calc.add_argument('amp', type=int)
    calc.add_argument('--mesh', type=int, choices=(2, 3), default=2)
    calc.add_argument('--ranks', type=int, required=True)
    calc.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    manifest, sources = source_files(DEFAULT_ROOT, DEFAULT_ROOT/MANIFEST_PATH, 'surface')
    site = slab_geometry(manifest, sources)
    atoms, ref = load_surface_checkpoint(site, args.reference, profile(site, 'base'), final=True)
    require(ref['sha256'] == REF_SHA, 'Slab pair requires the exact adopted reference bytes')
    binding = dict(file=args.reference.name, bytes=ref['bytes'], sha256=ref['sha256'])
    if args.action == 'plan':
        print(json.dumps([geometry(atoms, 'AB', amp)[1] for amp in AMPS], indent=2))
        return
    writer = runpy.run_path(sources['atomic_writer']['path'])['write_json']
    run(args, site, atoms, writer, binding)


if __name__ == '__main__':
    main()
