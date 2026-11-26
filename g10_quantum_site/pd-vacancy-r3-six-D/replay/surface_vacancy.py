"""G10 neutral PdD(111) and G7 C vacancy: relax a surface, then sample frozen wells."""
import argparse
from dataclasses import dataclass, replace
import hashlib
from importlib.metadata import version
from itertools import product
import json
from pathlib import Path
import runpy

import numpy as np
from ase import Atom, Atoms
from ase.constraints import FixAtoms
from ase.geometry import find_mic

MANIFEST_SHA = 'ee38e402445cbeb5afceceea2e4345c2980b16ecca78df05b9d965bae0f537f5'
DEFAULT_ROOT = Path('.')
MANIFEST_PATH = Path('coherence-g10-wide/sources/slab-vacancy/source-manifest.json')
DEUTERIUM_U = 2.014101778
SURFACE_FMAX = 0.05
SITE_KEYS = {
    'surface': ('g9_pilot_journal', 'g9_host_log', 'g9_slab_recipe', 'g9_calculator'),
    'vacancy': ('g7_selected', 'g7_standard_record', 'g7_standard_gpw', 'g7_scf_log',
                'g7_calculator', 'g7_calculator_base', 'g7_relaxation', 'g7_charged_probe'),
}
AXES = ('normal', 'in1', 'in2')
POINTS = {'zero': (None, 0.0), **{
    f'{axis}-{sign}{amp:02d}': (axis, (1 if sign == 'p' else -1) * amp / 100)
    for axis, sign, amp in product(AXES, ('p', 'm'), (12, 24, 36))
}}
BENCHMARKS = {'surface': 'normal-p12', 'vacancy': 'in1-p12'}
RELAX_SCHEMA = 'g10-surface-relaxation-v1'
POINT_SCHEMA = 'g10-site-point-v1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def file_identity(path):
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return dict(path=str(path.resolve()), bytes=path.stat().st_size, sha256=digest)


def read_record(path):
    path = path.resolve()
    data = path.read_bytes()
    return json.loads(data), dict(path=str(path), bytes=len(data), sha256=sha_bytes(data))


def source_files(root, manifest_path, kind):
    manifest, identity = read_record(manifest_path)
    require(identity['sha256'] == MANIFEST_SHA, 'Source manifest must be the pinned G10 inventory')
    require(all(version(key) == val for key, val in manifest['software'].items()),
            'Installed GPAW/ASE versions do not match the source manifest')
    sources = {'manifest': identity}
    for key in (*SITE_KEYS[kind], 'atomic_writer', 'Pd_PAW', 'H_PAW'):
        relpath, size, digest = manifest['files'][key]
        bound = file_identity(root / relpath)
        require((bound['bytes'], bound['sha256']) == (size, digest), f'Source bytes differ: {key}')
        sources[key] = bound
    return manifest, sources


def unit(vector):
    return vector / np.linalg.norm(vector)


def frame(normal, tangent):
    n = unit(normal)
    t = unit(tangent - (n @ tangent) * n)
    return dict(zip(AXES, (n, t, np.cross(n, t))))


def isotope_atoms(symbols, positions, cell, pbc):
    atoms = Atoms(symbols, positions=positions, cell=cell, pbc=pbc)
    masses = atoms.get_masses()
    masses[atoms.numbers == 1] = DEUTERIUM_U
    atoms.set_masses(masses)
    return atoms


def geometry(atoms):
    return dict(cell_A=atoms.cell.array.tolist(), pbc=atoms.pbc.tolist(),
                positions_A=atoms.positions.tolist(), symbols=atoms.get_chemical_symbols(),
                masses_u=atoms.get_masses().tolist(), indices=list(range(len(atoms))))


@dataclass
class Site:
    kind: str
    atoms: Atoms
    moving: int
    axes: dict
    fixed: list
    origin: dict
    sources: dict

    def site_record(self):
        return dict(moving_index=self.moving, ray_units={k: v.tolist() for k, v in self.axes.items()},
                    origin=self.origin, sources=self.sources)


def slab_geometry(manifest, sources):
    line_info = manifest['g9_geometry_line']
    with Path(sources['g9_pilot_journal']['path']).open('rb') as stream:
        line = stream.readline()
    require(line_info['line'] == 1 and len(line) == line_info['bytes_including_lf']
            and sha_bytes(line) == line_info['sha256'], 'G9 geometry line differs')
    row = json.loads(line)
    require(row['phase'] == 'geometry' and row['pd_D_symbols'] == ['Pd'] * 16 + ['H'] * 12,
            'Expected G9 16 Pd + 12 D ideal slab')
    slab = isotope_atoms(row['pd_D_symbols'], row['host_positions_A'], row['cell_A'], (True, True, False))
    receiver = row['receiving_indices'][0]
    require(receiver == 26 and slab.numbers[receiver] == 1, 'G9 receiving D order changed')
    z = slab.positions[:16, 2]
    levels, counts = np.unique(z, return_counts=True)
    require(counts.tolist() == [4, 4, 4, 4]
            and np.isclose(levels[-1], row['surface_z_A'], atol=1e-12, rtol=0), 'G9 Pd layers differ')
    layers = [np.flatnonzero(z == level) for level in levels]
    top, below, third = layers[-1], layers[-2], layers[-3]
    xy = slab.positions[receiver].copy()
    xy[2] = 0.0

    def lateral_distances(indices):
        offsets = slab.positions[indices] - xy
        offsets[:, 2] = 0.0
        return find_mic(offsets, slab.cell, slab.pbc)[1]

    # An fcc hollow in ABC stacking is directly above a third-layer Pd, not
    # a second-layer Pd (hcp hollow). Pick the hollow above receiving D 26.
    nns = np.argsort(lateral_distances(top))
    top_r = lateral_distances(top)[nns]
    third_r, below_r = lateral_distances(third), lateral_distances(below)
    d_nn = np.linalg.norm(slab.cell[0]) / 2
    require(np.count_nonzero(third_r < 1e-10) == 1 and min(below_r) > 1.0
            and np.allclose(top_r[:3], d_nn / np.sqrt(3), atol=1e-10, rtol=0)
            and top_r[3] > top_r[2] + 0.5, 'Projected receiving ray is not an fcc hollow')
    ideal = geometry(slab)
    axes = frame(np.cross(slab.cell[0], slab.cell[1]), slab.cell[0])
    require(np.allclose(axes['normal'], (0, 0, 1), atol=1e-12, rtol=0), 'Unexpected slab frame')
    height = 1.10
    ads = xy.copy()
    ads[2] = levels[-1] + height
    moving = len(slab)
    slab.append(Atom('H', ads, mass=DEUTERIUM_U))
    require(np.array_equal(slab.positions[:-1], row['host_positions_A']), 'Ideal host was repositioned')
    registry = dict(name='fcc hollow', adsorbate_index=moving, coverage_ML=1 / len(top),
                    initial_height_A=height, top_Pd_indices=sorted(top[nns[:3]].tolist()),
                    all_top_Pd_indices=top.tolist(), third_layer_Pd_index=int(third[np.argmin(third_r)]),
                    nearest_second_layer_lateral_distance_A=float(min(below_r)),
                    top_lateral_distances_A=top_r.tolist())
    origin = dict(g9_geometry_line=line_info, ideal_host=ideal, ideal_receiving_index=receiver,
                  ideal_receiving_indices=row['receiving_indices'], initial_registry=registry)
    return Site('surface', slab, moving, axes, layers[0].tolist(), origin, sources)


def vacancy_geometry(sources):
    selection = json.loads(Path(sources['g7_selected']['path']).read_bytes())
    record = json.loads(Path(sources['g7_standard_record']['path']).read_bytes())
    for name in ('g7_standard_record', 'g7_standard_gpw'):
        receipt = selection['source_frame'] if name.endswith('record') else selection['frames'][0]['gpw']
        require(all(sources[name][key] == receipt[key] for key in ('sha256', 'bytes')),
                f'G7 selected state does not bind {name}')
    spec = record['specification']
    require(selection['passed'] and spec['case'] == 'C' and spec['geometry_H_setup'] == 'standard_PAW'
            and record['ionic_step'] == 36 and record['scf_converged'] and record['setups'] == {'Pd': '10'}
            and record['charge'] == 0 and record['electrons'] == 71 and record['nbands'] == 52,
            'Expected the neutral G7 C standard-PAW step 36')
    require(all(record[key] == selection[key] for key in ('cell_A', 'positions_A', 'symbols', 'raw_fmax_eV_A'))
            and selection['symbols'] == ['Pd'] * 7 + ['H'] and selection['target_index'] == 7,
            'G7 selected source coordinates, force gate or target order differ')
    atoms = isotope_atoms(record['symbols'], record['positions_A'], record['cell_A'], True)
    moving = selection['target_index']
    registered = selection['vacancy']
    centre = np.array(registered['registered_centre_A'])
    offset, radius = find_mic(atoms.positions[moving] - centre, atoms.cell, atoms.pbc)
    source_forces = np.array(record['forces_eV_A'])
    fmax = np.linalg.norm(source_forces, axis=1).max()
    require(registered['hydrogen_indices'] == [moving]
            and np.isclose(radius, registered['hydrogen_distances_A'][0], atol=1e-12, rtol=0)
            and np.isclose(fmax, record['raw_fmax_eV_A'], atol=1e-12, rtol=0)
            and fmax <= selection['force_threshold_eV_A'], 'G7 trap coordinates or selected force gate differ')
    # There is no intrinsic surface normal. Use the cavity-to-D MIC ray and
    # two orthogonal tangents, in a fixed right-handed Cartesian frame.
    axes = frame(offset, np.array([1.0, 0.0, 0.0]))
    origin = dict(g7_case='C', selected_step=36, bound_index=moving, registered_vacancy_centre_A=centre.tolist(),
                  bound_D_offset_from_centre_A=offset.tolist(), bound_D_radius_A=float(radius),
                  source_energy_free_eV=record['energy']['free_energy_eV'], source_forces_eV_A=source_forces.tolist(),
                  source_raw_fmax_eV_A=float(fmax), source_gate_eV_A=selection['force_threshold_eV_A'])
    return Site('vacancy', atoms, moving, axes, [], origin, sources)


def profile(site, resolution):
    from gpaw.utilities import h2gpts

    coarse = resolution == 'base'
    h = 0.20 if coarse else 0.16
    surface = site.kind == 'surface'
    gpts = h2gpts(h, site.atoms.cell.array, idiv=8 if surface else 4).tolist()
    if coarse:
        require(gpts == ([24, 24, 136] if surface else [24, 24, 24]), 'Source grid mismatch')
    expected_electrons = int(sum(10 if n == 46 else 1 for n in site.atoms.numbers))
    return dict(engine='GPAW 26.7.0 new ASE API', mode='fd', xc='PBE', charge=0, spinpol=False,
                setups={'Pd': '10', 'H': 'paw'}, expected_electrons=expected_electrons,
                h_A=h, gpts=gpts, kmesh=[2, 2, 1] if surface else [1, 1, 1], gamma=not surface,
                nbands=(expected_electrons + 1) // 2 + (8 if surface else 16),
                symmetry={'point_group': False, 'time_reversal': False},
                occupations={'name': 'fermi-dirac', 'width': 0.1},
                poissonsolver={'name': 'fast', **({'metallic_electrodes': 'both'} if surface else {})},
                mixer={'backend': 'pulay', 'method': 'separate', 'beta': 0.10 if surface else 0.05,
                       'nmaxold': 5, 'weight': 50 if surface else 100},
                convergence={'energy': 1e-6, 'density': 1e-5, 'eigenstates': 1e-8},
                maxiter=160 if surface else 150)


def calculator(config, log, world):
    from gpaw.new.ase_interface import GPAW

    nk = int(np.prod(config['kmesh']))
    require(nk % world.size == 0, f'Use a rank count dividing {nk} source k-points')
    keys = ('mode', 'xc', 'charge', 'spinpol', 'setups', 'gpts', 'nbands', 'symmetry',
            'occupations', 'poissonsolver', 'mixer', 'convergence', 'maxiter')
    return GPAW(**{key: config[key] for key in keys}, kpts={'size': config['kmesh'], 'gamma': config['gamma']},
                parallel={'kpt': world.size, 'domain': 1, 'band': 1}, communicator=world, txt=str(log))


def native(atoms, config):
    energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces(apply_constraint=False)
    calc = atoms.calc
    dft = calc.dft
    wave = dft.ibzwfs
    iters = int(calc.get_number_of_iterations())
    grid = [int(n) for n in calc.get_number_of_grid_points()]
    electrons = float(calc.get_number_of_electrons())
    kpts = calc.get_ibz_k_points()
    require(np.isfinite(energy) and forces.shape == (len(atoms), 3) and np.isfinite(forces).all()
            and np.array_equal(calc.atoms.positions, atoms.positions), 'Invalid or wrong-position native quantities')
    require(0 < iters <= config['maxiter'] and dft.scf_loop.update_density_and_potential
            and grid == config['gpts'] and electrons == config['expected_electrons']
            and calc.get_number_of_bands() == config['nbands'] and len(kpts) == int(np.prod(config['kmesh']))
            and wave.kpt_comm.size == calc.comm.size and wave.domain_comm.size == wave.band_comm.size == 1
            and len(wave.ibz.symmetries) == 1, 'Native density, grid, electron, band or k-point topology differs')
    return dict(**geometry(atoms), energy_free_eV=energy, forces_eV_A=forces.tolist(),
                net_force_eV_A=forces.sum(axis=0).tolist(), raw_fmax_eV_A=float(np.linalg.norm(forces, axis=1).max()),
                calculator=config, scf_converged=True, scf_iterations=iters, native_grid=grid,
                electron_count=electrons, ibz_kpoints=kpts.tolist(), ibz_weights=calc.get_k_point_weights().tolist(),
                mpi_ranks=calc.comm.size)


def root_action(world, fn):
    from ase.parallel import broadcast

    result = None
    if world.rank == 0:
        try:
            result = (True, fn())
        except Exception as exc:
            result = (False, f'{type(exc).__name__}: {exc}')
    ok, value = broadcast(result, comm=world)
    if not ok:
        raise RuntimeError(value)
    return value


def publish(world, writer, path, record):
    return root_action(world, lambda: writer(path, record))


def adsorbate_geometry(site, atoms):
    registry = site.origin['initial_registry']
    top_ids = np.array(registry['all_top_Pd_indices'])
    vectors, lengths = find_mic(atoms.positions[site.moving] - atoms.positions[top_ids], atoms.cell, atoms.pbc)
    receiver = site.origin['ideal_receiving_index']
    shift, distance = find_mic(atoms.positions[receiver] - site.atoms.positions[receiver], atoms.cell, atoms.pbc)
    nearest = np.argsort(lengths)[:3]
    return dict(height_above_nearest_three_top_Pd_A=float(vectors[nearest, 2].mean()),
                nearest_three_top_Pd_indices=top_ids[nearest].tolist(),
                nearest_three_top_Pd_distances_A=lengths[nearest].tolist(),
                seed_to_current_D_MIC_A=find_mic(atoms.positions[site.moving] - site.atoms.positions[site.moving],
                                                 atoms.cell, atoms.pbc)[0].tolist(),
                subsurface_index=receiver, subsurface_displacement_from_G9_ideal_A=shift.tolist(),
                subsurface_displacement_norm_A=float(distance))


def load_surface_checkpoint(site, path, config, *, final):
    record, identity = read_record(path)
    g = geometry(site.atoms)
    require(record['schema'] == RELAX_SCHEMA and record['site'] == 'surface'
            and record['moving_index'] == site.moving and record['fixed_indices'] == site.fixed
            and record['calculator'] == config and record['scf_converged'], 'Wrong surface relaxation checkpoint')
    require(all(record['sources'][k]['sha256'] == v['sha256'] for k, v in site.sources.items())
            and all(record[k] == g[k] for k in ('cell_A', 'pbc', 'symbols', 'masses_u', 'indices')),
            'Surface checkpoint source or atom order differs')
    positions = np.asarray(record['positions_A'])
    require(positions.shape == site.atoms.positions.shape and np.isfinite(positions).all()
            and np.array_equal(positions[site.fixed], site.atoms.positions[site.fixed]),
            'Surface checkpoint positions or fixed-Pd coordinates differ')
    atoms = site.atoms.copy()
    atoms.positions[:] = positions
    if final:
        forces = np.asarray(record['forces_eV_A'])
        mobile = np.setdiff1d(g['indices'], site.fixed)
        require(record['ionic_converged'] and forces.shape == positions.shape and np.isfinite(forces).all()
                and np.linalg.norm(forces[mobile], axis=1).max() <= SURFACE_FMAX
                and adsorbate_geometry(site, atoms)['height_above_nearest_three_top_Pd_A'] > 0,
                'Reference is not an admitted above-surface relaxation endpoint')
    return atoms, identity


def relax(site, args, world, writer):
    from ase.optimize import BFGS, LBFGS

    require(args.steps > 0, 'Ionic step limit must be positive')
    cfg = profile(site, 'base')
    atoms, resume = (load_surface_checkpoint(site, args.resume_from, cfg, final=False)
                     if args.resume_from else (site.atoms.copy(), None))
    root_action(world, lambda: args.output.mkdir(parents=True, exist_ok=False))
    atoms.set_constraint(FixAtoms(indices=site.fixed))
    atoms.calc = calculator(cfg, args.output / 'relax.scf', world)
    if world.rank == 0:
        print('G10_SURFACE_RELAX_START', flush=True)
    optimizer = LBFGS if resume else BFGS
    maxstep = 0.05 if resume else 0.10
    with optimizer(atoms, logfile=args.output / 'optimizer.log', maxstep=maxstep, comm=world) as opt:
        for converged in opt.irun(fmax=SURFACE_FMAX, steps=args.steps):
            state = native(atoms, cfg)
            mobile = np.setdiff1d(state['indices'], site.fixed)
            fmax = float(np.linalg.norm(np.asarray(state['forces_eV_A'])[mobile], axis=1).max())
            location = adsorbate_geometry(site, atoms)
            state.update(schema=RELAX_SCHEMA, site=site.kind, **site.site_record(), fixed_indices=site.fixed,
                         attempt_ionic_step=opt.nsteps, ionic_converged=bool(converged),
                         ionic_fmax_eV_A=fmax, ionic_threshold_eV_A=SURFACE_FMAX,
                         ionic_optimizer=optimizer.__name__, ionic_max_step_A=maxstep,
                         adsorption_geometry=location, resume_source=resume)
            receipt = publish(world, writer, args.output / f'step-{opt.nsteps:04d}.json', state)
            if world.rank == 0:
                print(f'G10_SURFACE_RELAX step={opt.nsteps} ionic_fmax_eV_A={fmax:.8f}', flush=True)
            if converged:
                require(location['height_above_nearest_three_top_Pd_A'] > 0,
                        'Tracked D is not above the terminal Pd face; raw endpoint is saved')
                state['selected_step'] = receipt
                publish(world, writer, args.output / 'reference.json', state)
                return
    raise RuntimeError('Surface ionic stop reached without convergence; completed step checkpoints are retained')


def displaced(site, atoms, case):
    axis, step = POINTS[case]
    ray = site.axes[axis] if axis is not None else np.zeros(3)
    centre = atoms.positions[site.moving].copy()
    moved = atoms.copy()
    moved.positions[site.moving] += step * ray
    others = np.arange(len(atoms)) != site.moving
    require(np.array_equal(atoms.positions[others], moved.positions[others])
            and np.allclose(moved.positions[site.moving] - centre, step * ray, atol=1e-12, rtol=0),
            'Frozen displacement changed unexpected degrees of freedom')
    return moved, dict(case=case, axis=axis, ray_unit=ray.tolist(), displacement_A=step,
                       reference_D_position_A=centre.tolist())


def frozen_point(site, args, world, writer, benchmark):
    base_cfg = profile(site, 'base')
    if site.kind == 'surface':
        require(args.reference is not None, 'Surface points require the converged relaxation reference.json')
        reference, identity = load_surface_checkpoint(site, args.reference, base_cfg, final=True)
        if args.target == 'subsurface':
            site = replace(site, moving=site.origin['ideal_receiving_index'])
    else:
        require(args.reference is None and args.target is None,
                'Vacancy uses the G7 C step-36 geometry and its bound D; omit surface selectors')
        reference, identity = site.atoms, site.sources['g7_selected']
    case = args.case or BENCHMARKS[site.kind]
    atoms, point = displaced(site, reference, case)
    cfg = profile(site, 'tight' if benchmark else 'base')
    if benchmark:
        base, base_id = read_record(args.baseline)
        require(base['schema'] == POINT_SCHEMA and base['site'] == site.kind and base['point'] == point
                and base['calculator'] == base_cfg and base['resolution'] == 'base' and base['scf_converged']
                and base['reference']['sha256'] == identity['sha256'] and base['moving_index'] == site.moving
                and all(base[k] == v for k, v in geometry(atoms).items())
                and all(base['sources'][k]['sha256'] == v['sha256'] for k, v in site.sources.items()),
                'Benchmark needs the matched native baseline point on the exact same reference and geometry')
        baseline_forces = np.asarray(base['forces_eV_A'])
        require(np.isfinite(base['energy_free_eV']) and baseline_forces.shape == atoms.positions.shape
                and np.isfinite(baseline_forces).all(), 'Invalid baseline observables')
    path = args.output
    log = path.with_name(path.name + '.scf')

    def prepare():
        path.parent.mkdir(parents=True, exist_ok=True)
        require(not path.exists() and not log.exists(), 'An existing point or SCF attempt cannot be overwritten')

    root_action(world, prepare)
    atoms.calc = calculator(cfg, log, world)
    if world.rank == 0:
        print(f'G10_SITE_START site={site.kind} case={case} grid={cfg["gpts"]}', flush=True)
    row = native(atoms, cfg)
    row.update(schema=POINT_SCHEMA, site=site.kind, **site.site_record(), reference=identity,
               resolution='tight' if benchmark else 'base', point=point,
               frozen_indices=[idx for idx in range(len(atoms)) if idx != site.moving],
               force_parallel_eV_A=(float(np.dot(row['forces_eV_A'][site.moving], point['ray_unit']))
                                    if point['axis'] is not None else None))
    if benchmark:
        row['benchmark'] = dict(base_point=base_id, changed_setting='FD grid spacing only',
                                energy_free_delta_tight_minus_base_eV=row['energy_free_eV'] - base['energy_free_eV'],
                                force_delta_tight_minus_base_eV_A=(np.asarray(row['forces_eV_A']) - baseline_forces).tolist())
    publish(world, writer, path, row)
    if world.rank == 0:
        print(f'G10_SITE_DONE site={site.kind} case={case} energy_free_eV={row["energy_free_eV"]:.10f}', flush=True)


def inspect(site):
    atoms = site.atoms
    rays = np.stack(list(site.axes.values()))
    require(np.allclose(rays @ rays.T, np.eye(3), atol=1e-12, rtol=0) and np.linalg.det(rays) > 0,
            'Local frame is not right-handed orthonormal')
    cases = {}
    for case in POINTS:
        moved, point = displaced(site, atoms, case)
        indices = np.delete(np.arange(len(atoms)), site.moving)
        distances = find_mic(moved.positions[indices] - moved.positions[site.moving], moved.cell, moved.pbc)[1]
        cases[case] = dict(**point, displaced_D_position_A=moved.positions[site.moving].tolist(),
                           nearest_ion_index=int(indices[np.argmin(distances)]), nearest_ion_distance_A=float(min(distances)))
    base, tight = (profile(site, tier) for tier in ('base', 'tight'))
    require({k for k in base if base[k] != tight[k]} == {'h_A', 'gpts'}
            and np.all(np.asarray(tight['gpts']) > base['gpts']), 'Benchmark does not refine only the actual grid')
    return dict(schema='g10-site-geometry-inspection-v1', site=site.kind, **site.site_record(),
                initial=geometry(atoms), fixed_indices_during_relaxation=site.fixed,
                proposed_displacements_from_initial=cases, actual_point_reference=('surface reference.json' if site.fixed
                                                                                   else 'G7 C step 36 as given'),
                base_calculator=base, tight_calculator=tight, benchmark_case=BENCHMARKS[site.kind])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=DEFAULT_ROOT, help='Root of the pinned original inputs')
    p.add_argument('--manifest', type=Path, help='Optional relocated copy of the exact pinned manifest')
    p.add_argument('--ranks', type=int, help='Bind each rank to a core in the admitted manager affinity mask')
    sub = p.add_subparsers(dest='action', required=True)
    info = sub.add_parser('inspect', help='Geometry and source inspection only; no SCF')
    info.add_argument('site', choices=SITE_KEYS)
    info.add_argument('--output', type=Path, help='Omit to write inspection JSON to stdout')
    opt = sub.add_parser('relax', help='Independently relax the covered surface; output is a new directory')
    opt.set_defaults(site='surface')
    opt.add_argument('--output', type=Path, required=True)
    opt.add_argument('--steps', type=int, default=80)
    opt.add_argument('--resume-from', type=Path, help='Restart from a completed geometry using LBFGS with 0.05 Å steps; no SCF or optimizer history is reused')
    for action in ('point', 'benchmark'):
        q = sub.add_parser(action, help='Evaluate one frozen geometry with a new native SCF')
        q.add_argument('site', choices=SITE_KEYS)
        q.add_argument('case', choices=POINTS, nargs=None if action == 'point' else '?')
        q.add_argument('--reference', type=Path, help='Converged surface reference.json; omit for vacancy')
        q.add_argument('--target', choices=('adsorbed', 'subsurface'),
                       help='Select original D 28 (default) or G9 receiving D 26 on the same relaxed surface')
        q.add_argument('--output', type=Path, required=True, help='New result JSON path')
        if action == 'benchmark':
            q.add_argument('--baseline', type=Path, required=True, help='Existing exact same base-grid point JSON')
    args = p.parse_args()
    manifest, sources = source_files(args.root, args.manifest or args.root / MANIFEST_PATH, args.site)
    site = slab_geometry(manifest, sources) if args.site == 'surface' else vacancy_geometry(sources)
    writer = runpy.run_path(sources['atomic_writer']['path'])['write_json']
    import gpaw
    from gpaw.mpi import world
    if args.ranks is not None:
        from loading_pd import bind_ranks
        bind_ranks(world, args.ranks)

    gpaw.setup_paths[:] = [str(Path(sources['H_PAW']['path']).parent)]
    if args.action == 'inspect':
        evidence = inspect(site)
        if args.output:
            publish(world, writer, args.output, evidence)
        elif world.rank == 0:
            print(json.dumps(evidence, indent=2, allow_nan=False))
    elif args.action == 'relax':
        relax(site, args, world, writer)
    else:
        frozen_point(site, args, world, writer, args.action == 'benchmark')


if __name__ == '__main__':
    main()
