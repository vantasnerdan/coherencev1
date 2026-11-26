"""Neutral Pd vacancy with two or six D, at matched primitive-cell sizes."""
import argparse
from importlib.metadata import version
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np
from ase.geometry import find_mic

from dense_pd import commit_json
from surface_vacancy import file_identity, geometry, isotope_atoms, publish, read_record, require, root_action

A = 3.940
PRIMITIVE = A/2*(np.ones((3, 3)) - np.eye(3))
FMAX = 0.05
SCHEMA = 'g10-pd-mono-vacancy-multi-D-v1'
PAW = {
    'Pd.10.PBE.gz': (106793, '3c6cfd1de7e6d6815d58c9521b2bc7cec27ede8e83554ed1f0151f979761bbb4'),
    'H.PBE.gz': (12501, '3136ad73c02d471f99dd5f82e451b2061952970556a260d38b53b5020dc94ae7'),
}
LITERATURE = {
    'vek2009.pdf': (193128, '7eb62053b6dcfda508a85eef370b675a6161d82b2cf00de75cf955755339ea8d'),
    'pd2020.html': (191269, 'e1374ddd0837bb3312fdbae1483bb076b75d0ef9b133d4f9b3b82da63a9ef7d8'),
    'pd2020-fig2.jpg': (65353, 'f154adb524dc4f3fdefff9d66c098b40f6359459244b5216617256af33e3159b'),
}


def verify(paw, articles):
    sources = {}
    for dirpath, bindings in ((paw, PAW), (articles, LITERATURE)):
        for name, expected in bindings.items():
            found = file_identity(dirpath/name)
            require((found['bytes'], found['sha256']) == expected, f'Actual source bytes do not match {name}')
            sources[name] = {'file': name, 'sha256': found['sha256'], 'bytes': found['bytes']}
    require(version('gpaw') == '26.7.0' and version('ase') == '3.29.0', 'Unexpected software version')
    return sources


def seed(repeat, occupancy):
    cell = repeat*PRIMITIVE
    centre = 0.5*cell.sum(axis=0)
    grid = range(-(repeat//2), repeat-repeat//2)
    keys = [p for p in product(grid, repeat=3) if p != (0, 0, 0)]
    corners = [(-1, -1, -1), (-1, 1, -1)]
    if occupancy == 6:
        corners = [p for p in product((-1, 1), repeat=3) if p not in ((-1, 1, -1), (1, 1, 1))]
    pd_xyz = np.asarray(keys) @ PRIMITIVE
    d_xyz = 1.40/np.sqrt(3)*np.asarray(corners)
    atoms = isotope_atoms(['Pd']*len(keys) + ['H']*occupancy, np.vstack((pd_xyz, d_xyz))+centre, cell, True)
    require(len(atoms) == repeat**3 - 1 + occupancy, 'Wrong vacancy site count')
    return atoms, centre, corners


def distances(atoms, centre, occupancy):
    indices = tuple(range(len(atoms)-occupancy, len(atoms)))
    xyz, radii = find_mic(atoms.positions[list(indices)] - centre, atoms.cell, atoms.pbc)
    i, j = np.asarray(list(combinations(range(occupancy), 2))).T
    _, mic = find_mic(xyz[j]-xyz[i], atoms.cell, atoms.pbc)
    direct = np.linalg.norm(xyz[j]-xyz[i], axis=1)
    images = np.array([p for p in product(range(-2, 3), repeat=3) if p != (0,0,0)]) @ atoms.cell.array
    cross = np.linalg.norm((xyz[:,None] - xyz[None,:])[:,:,None,:] + images, axis=-1).min()
    return dict(D_radii_A=np.asarray(radii).tolist(),
                D_D=[dict(i=indices[k], j=indices[l], D_label_i=int(k), D_label_j=int(l),
                          separation_A=float(s), centre_image_pair=bool(abs(s-d)<1e-10))
                     for k,l,s,d in zip(i,j,mic,direct,strict=True)],
                nearest_other_vacancy_image_D_D_A=float(cross),
                all_trapped_near_centre=bool(np.max(radii) < A/2))


def profile(atoms):
    ne = int(sum(10 if n == 46 else 1 for n in atoms.numbers))
    return dict(expected_electrons=ne, cutoff_eV=500, kmesh=[2,2,2], gamma=True, xc='PBE',
                convergence=dict(energy=1e-6, density=1e-5, forces=1e-4, eigenstates=4e-8, bands='occupied'),
                setups={'Pd': '10', 'H': 'paw'}, charge=0, spinpol=False, eigensolver='davidson',
                symmetry={'point_group': False, 'time_reversal': False},
                smear_eV=0.10, maxiter=175, nbands=ne//2+24)


def case_input(args):
    sources = verify(args.paw, args.articles)
    atoms, centre, tags = seed(args.repeat, args.occupancy)
    metrics = distances(atoms, centre, args.occupancy)
    require(all(p['centre_image_pair'] for p in metrics['D_D']) and
            metrics['all_trapped_near_centre'] and metrics['nearest_other_vacancy_image_D_D_A'] > 5.0,
            'Invalid registration or image isolation at the seed')
    spec = dict(schema=SCHEMA, repeat=args.repeat, occupancy=args.occupancy, seed_corner_signs=tags,
                vacancy_centre_A=centre.tolist(), calculator=profile(atoms), source_artifacts=sources)
    return atoms, centre, spec, metrics


def saved_seed(atoms, spec, resume):
    record, identity = read_record(resume)
    g = geometry(atoms)
    xyz = np.asarray(record['positions_A'])
    require(record['specification'] == json.loads(json.dumps(spec)) and record['scf_converged'] is True
            and all(record[k] == g[k] for k in ('cell_A','pbc','symbols','masses_u','indices'))
            and xyz.shape == atoms.positions.shape and np.isfinite(xyz).all(), 'Wrong resume geometry or setting')
    atoms.positions[:] = xyz
    return dict(file=resume.name, sha256=identity['sha256'], bytes=identity['bytes'])


def relax(args, atoms, centre, spec):
    import gpaw
    from gpaw import FermiDirac, Mixer, PW
    from gpaw.mpi import world
    from gpaw.new.ase_interface import GPAW
    from ase.optimize import LBFGS
    from loading_pd import bind_ranks

    bind_ranks(world, args.ranks)
    require(8 % world.size == 0 and args.steps > 0, 'Use k-only rank count dividing eight and bounded ionic steps')
    gpaw.setup_paths[:] = [str(args.paw.resolve())]
    cfg = spec['calculator']
    out = args.output
    resume = saved_seed(atoms, spec, args.resume_from) if args.resume_from else None
    root_action(world, lambda: out.mkdir(parents=True, exist_ok=False))
    keys = ('xc','convergence','setups','charge','spinpol','eigensolver','symmetry','maxiter','nbands')
    atoms.calc = calc = GPAW(**{k:cfg[k] for k in keys}, mode=PW(cfg['cutoff_eV']),
                            kpts={'size':cfg['kmesh'],'gamma':cfg['gamma']}, occupations=FermiDirac(cfg['smear_eV']),
                            mixer=Mixer(0.04, 7, 80), parallel={'kpt':world.size,'domain':1,'band':1},
                            communicator=world, txt=str(out/'relax.scf'))
    progress_force, progress_step = float('inf'), 0
    outside_steps = np.zeros(args.occupancy, dtype=int)
    with LBFGS(atoms, logfile=out/'relax.log', maxstep=0.12, comm=world) as opt:
        for converged in opt.irun(fmax=FMAX, steps=args.steps):
            force = atoms.get_forces(apply_constraint=False)
            free = float(atoms.get_potential_energy(force_consistent=True))
            max_force = float(np.linalg.norm(force, axis=1).max())
            stats = calc.dft.ibzwfs
            iterations = int(calc.get_number_of_iterations())
            ne = float(calc.get_number_of_electrons())
            kp = calc.get_ibz_k_points()
            require(np.isfinite(free) and np.isfinite(force).all() and np.array_equal(calc.atoms.positions, atoms.positions)
                    and 0 < iterations <= cfg['maxiter'] and ne == cfg['expected_electrons']
                    and len(kp) == 8 and calc.get_number_of_bands() == cfg['nbands']
                    and stats.kpt_comm.size == world.size and stats.domain_comm.size == stats.band_comm.size == 1
                    and calc.dft.scf_loop.update_density_and_potential, 'Native result or topology failed qualification')
            complete = bool(converged and max_force <= FMAX)
            row = dict(**geometry(atoms), specification=spec, ionic_step=opt.nsteps, energy_free_eV=free,
                       forces_eV_A=force.tolist(), raw_fmax_eV_A=max_force, scf_iterations=iterations,
                       scf_converged=True, ionic_converged=complete, electron_count=ne, ibz_kpoints=kp.tolist(),
                       mpi_ranks=world.size, resume=resume, distances=distances(atoms, centre, args.occupancy))
            radii = np.asarray(row['distances']['D_radii_A'])
            # Do not stop on a single displacement near the neighbouring octahedral shell.
            outside_steps = np.where(radii > A/2 + 0.25, outside_steps + 1, 0)
            if max_force <= 0.90 * progress_force:
                progress_force, progress_step = max_force, opt.nsteps
            publish(world, commit_json, out/f'step-{opt.nsteps:04}.json', row)
            if world.rank == 0:
                print(f'G10_PD_VACANCY r{args.repeat} d{args.occupancy} step={opt.nsteps} fmax={max_force:.8f}', flush=True)
            if np.any(outside_steps >= 3) or (complete and not row['distances']['all_trapped_near_centre']):
                if world.rank == 0:
                    print('G10_PD_VACANCY_STOP cage-retention', flush=True)
                return 4
            if complete:
                publish(world, commit_json, out/'reference.json', row)
                return 0
            if opt.nsteps - progress_step >= 12:
                if world.rank == 0:
                    print('G10_PD_VACANCY_STOP force-plateau', flush=True)
                return 3
    return 2  # Completed raw steps are retained, but not a resting geometry.


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--paw', type=Path, required=True)
    p.add_argument('--articles', type=Path, required=True)
    p.add_argument('--repeat', type=int, choices=(3,4), required=True)
    p.add_argument('--occupancy', type=int, choices=(2,6), required=True)
    sub = p.add_subparsers(dest='action', required=True)
    sub.add_parser('plan')
    run = sub.add_parser('relax')
    run.add_argument('--ranks', type=int, required=True)
    run.add_argument('--steps', type=int, required=True)
    run.add_argument('--resume-from', type=Path)
    run.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    atoms, centre, spec, metrics = case_input(args)
    if args.action == 'plan':
        print(json.dumps(dict(**geometry(atoms), specification=spec, distances=metrics), indent=2))
        return 0
    return relax(args, atoms, centre, spec)


if __name__ == '__main__':
    raise SystemExit(main())
