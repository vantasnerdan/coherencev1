"""Prepare and evaluate a constrained bulk fluorite TiD2 static fingerprint.

plan does not call GPAW. relax, point and readout are separate, explicit stages;
use a managed launch for relax/point. No launch, retry or rate model lives here.
"""

import argparse
import hashlib
from importlib.metadata import version
from itertools import product
import json
import math
from pathlib import Path

from ase import Atoms
import numpy as np
from scipy.optimize import brentq

from dense_pd import commit_json, read_json, require

SCHEMA = 'g10-fluorite-tid2-v1'
SEED_A = 4.424
SOURCE_BYTES = {
    'er-et-al-0910.1003v1.pdf': (963946, '79339285761f7ab2fd1e3262d068ac684da53f0d2514dedb8e9a52992f4db16c'),
    'Ti.PBE.gz': (101840, '64143c102107cae1ac7766e2a3292fb5ab6d0fa8122b7f8be5b25c29a4a40a2c'),
    'H.PBE.gz': (12501, '3136ad73c02d471f99dd5f82e451b2061952970556a260d38b53b5020dc94ae7'),
}
SOFTWARE = {'gpaw': '26.7.0', 'ase': '3.29.0'}
D_MASS_U = 2.01410177812
SCF = {'energy': 1e-7, 'density': 1e-4, 'forces': 1e-4, 'eigenstates': 4e-8, 'bands': 'occupied'}
RELAX_STRESS = 5e-4  # eV/angstrom^3, each Cartesian stress component at zero external pressure
MAX_FORCE = 5e-3     # eV/angstrom, largest atomic norm at the undisturbed reference
MAX_NET = 5e-3       # eV/angstrom, diagnostic for displaced points
CHECK_E = 5e-3      # eV, nominal case-specific setting diagnostic
CHECK_F = 2e-2      # eV/angstrom, maximum vector norm across sites

MILLANGSTROM = (120, 240, 480, 800)
DIRECTIONS = {'100': (1, 0, 0), '110': (1, 1, 0), '111': (1, 1, 1)}
CELL_PER_A = 0.5 * (np.ones((3, 3)) - np.eye(3))
SITE_A, SITE_B = 8, 19
PAIR_MILLIANGSTROM = (120, 240, 480, 600)
CHECK_CASE = '111m480'


def binding(path):
    with path.open('rb') as stream:
        return {'file': path.name, 'bytes': path.stat().st_size,
                'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}


def source_bindings(source):
    artifacts = [binding(source / name) for name in SOURCE_BYTES]
    for artifact in artifacts:
        expected = SOURCE_BYTES[artifact['file']]
        require((artifact['bytes'], artifact['sha256']) == expected,
                f"Source bytes differ: {artifact['file']}")
    return artifacts


def fluoritic_atoms(a, repeat):
    """Ordered Ti, D+, D− blocks; cell rows in conventional cubic Cartesian axes."""
    require(np.isfinite(a) and a > 0, 'Invalid conventional cubic lattice constant')
    t = np.array(list(product(range(repeat), repeat=3))) @ CELL_PER_A
    sites = np.concatenate((t, t + 0.25, t - 0.25))
    atoms = Atoms(['Ti'] * len(t) + ['H'] * (2 * len(t)), positions=a * sites,
                  cell=(a * repeat) * CELL_PER_A, pbc=True)
    masses = atoms.get_masses()
    masses[len(t):] = D_MASS_U
    atoms.set_masses(masses)
    return atoms


def geometry(atoms):
    return dict(pbc=atoms.pbc.tolist(), cell_A=atoms.cell.array.tolist(),
                ase_symbols=atoms.get_chemical_symbols(), atom_index=list(range(len(atoms))),
                physical_species=['D' if n == 1 else 'Ti' for n in atoms.numbers],
                masses_u=atoms.get_masses().tolist(), positions_A=atoms.positions.tolist())


def pair_geometry(atoms, a):
    ray = np.array(DIRECTIONS['100'], dtype=float)
    direct = atoms.positions[SITE_B] - atoms.positions[SITE_A]
    mic = atoms.get_distance(SITE_A, SITE_B, mic=True, vector=True)
    distances = atoms.get_distances(SITE_A, [i for i in range(8, 24) if i != SITE_A], mic=True)
    require(atoms.numbers[[SITE_A, SITE_B]].tolist() == [1, 1] and
            np.allclose(direct, (a / 2) * ray, atol=1e-12, rtol=0) and
            np.allclose(mic, direct, atol=1e-12, rtol=0) and
            np.isclose(distances.min(), a / 2, atol=1e-12, rtol=0),
            'Pair is not the registered +[100] nearest-neighbour image')
    return dict(A_index=SITE_A, B_index=SITE_B, A_to_B_A=direct.tolist(),
                B_image_cell_translation=[0, 0, 0], d0_A=float(np.linalg.norm(direct)),
                closing_amplitudes_per_site_A=[u/1000 for u in PAIR_MILLIANGSTROM])


def settings(label):
    mesh = {'relax': 8, 'base': 4, 'check': 6}[label]
    return dict(label=label, mode='PW', cutoff_eV=500, xc='PBE',
                setups={'Ti': 'paw', 'H': 'paw'}, kmesh=[mesh] * 3, gamma=True,
                nbands=11 if label == 'relax' else 76,
                charge=0, spinpol=False, smearing_eV=0.1,
                symmetry={'point_group': label == 'relax', 'time_reversal': True},
                eigensolver={'name': 'davidson', 'niter': 2}, mixer=[0.05, 5, 50],
                maxiter=160, convergence=SCF, zero_external_pressure=True)


def profile(source):
    software = {name: version(name) for name in (*SOFTWARE, 'numpy', 'scipy')}
    require(all(software[k] == v for k, v in SOFTWARE.items()), 'Wrong GPAW or ASE version')
    return dict(schema=SCHEMA, material='cubic-fluorite-TiD2', software=software,
                seed_source=dict(a_A=SEED_A, citation='Er et al., arXiv:0910.1003v1, Table 1, p. 7',
                                 role='documented f-TiH2 electronic-structure geometry, not measured TiD'),
                source_artifacts=source_bindings(source), code=binding(Path(__file__)),
                energy_convention='GPAW force-consistent electronic free energy, eV per simulated cell',
                force_convention='native full Cartesian electronic forces, eV/angstrom, zero-based atom_index',
                stress_convention='native ASE Cartesian tensile-positive stress, eV/angstrom^3',
                scope='periodic frozen-ion Born–Oppenheimer fingerprint; no isotope-dependent statics or foil phase assignment')


def case_specs():
    cases = {'zero': dict(setting='base', direction=None, moves_A=[])}
    for label, millia, sign in product(DIRECTIONS, MILLANGSTROM, ('m', 'p')):
        signed = millia / 1000 * (-1 if sign == 'm' else 1)
        cases[f'{label}{sign}{millia:03d}'] = dict(setting='base', direction=label,
                                                   moves_A=[[SITE_A, signed]])
    for millia in PAIR_MILLIANGSTROM:
        u = millia/1000
        cases.setdefault(f'100p{millia:03d}', dict(setting='base', direction='100', moves_A=[[SITE_A,u]]))
        cases[f'pair-b-{millia:03d}'] = dict(setting='base', direction='100', moves_A=[[SITE_B,-u]])
        cases[f'pair-ab-{millia:03d}'] = dict(setting='base', direction='100', moves_A=[[SITE_A,u],[SITE_B,-u]])
    for name in ('zero', CHECK_CASE):
        cases['check-' + name] = dict(cases[name], setting='check')
    return cases


def ray_for(direction):
    ray = np.array(DIRECTIONS[direction], dtype=float)
    return ray / np.linalg.norm(ray)


def point_geometry(a, spec):
    atoms = fluoritic_atoms(a, 2)
    if spec['moves_A']:
        ray = ray_for(spec['direction'])
        for i, shift in spec['moves_A']:
            atoms.positions[i] += shift * ray
    return atoms


def runtime(source):
    import gpaw
    from gpaw.mpi import world

    require(world.size == 1, 'Use one GPAW process per case; run independent cases in separate processes')
    gpaw.setup_paths[:] = [str(source.resolve())]


def calculate(atoms, setting, log, stress=False):
    from gpaw import FermiDirac, GPAW, Mixer, PW

    with log.open('x') as stream:
        atoms.calc = calc = GPAW(
            legacy_gpaw=False, xc=setting['xc'], mode=PW(setting['cutoff_eV']),
            kpts={'size': setting['kmesh'], 'gamma': setting['gamma']}, setups=setting['setups'],
            charge=setting['charge'], spinpol=setting['spinpol'], nbands=setting['nbands'],
            occupations=FermiDirac(setting['smearing_eV']), symmetry=setting['symmetry'],
            mixer=Mixer(*setting['mixer']), eigensolver=setting['eigensolver'],
            convergence=setting['convergence'], maxiter=setting['maxiter'], txt=stream,
        )
        row = native_snapshot(atoms, calc, setting, stress)
        atoms.calc = None
        del calc
    return dict(row, scf_log=binding(log))


def native_snapshot(atoms, calc, setting, stress):
    initial = atoms.positions.copy()
    energy = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces(apply_constraint=False)
    sigma = atoms.get_stress(voigt=False, apply_constraint=False) if stress else None
    wfs = calc.dft.ibzwfs
    electrons = float(calc.get_number_of_electrons())
    require(electrons == 14 * len(atoms) // 3 and calc.dft.scf_loop.update_density_and_potential,
            'Unexpected valence count or fixed-density calculation')
    require(setting['symmetry']['point_group'] or len(wfs.ibz.symmetries) == 1,
            'Point-group reduction was not disabled for the frozen geometry')
    require(np.array_equal(atoms.positions, initial) and forces.shape == initial.shape and
            np.isfinite(energy) and np.isfinite(forces).all() and
            (sigma is None or (sigma.shape == (3, 3) and np.isfinite(sigma).all())),
            'Nonfinite native result or changed physical coordinates')
    eigensolver = type(calc.dft.scf_loop.eigensolver).__name__
    return dict(geometry(atoms), energy_free_eV=energy, forces_eV_A=forces.tolist(),
                stress_eV_A3=None if sigma is None else sigma.tolist(),
                mean_stress_eV_A3=None if sigma is None else float(np.trace(sigma) / 3),
                net_force_eV_A=forces.sum(axis=0).tolist(),
                max_force_eV_A=max_vector_norm(forces), electron_count=electrons,
                actual_nbands=int(calc.get_number_of_bands()), scf_iterations=calc.dft.scf_loop.niter,
                actual_eigensolver=eigensolver, ibz_kpoints=len(calc.get_ibz_k_points()))


def max_vector_norm(vectors):
    return float(np.linalg.norm(vectors, axis=1).max())


def stress_max(row):
    return float(np.abs(row['stress_eV_A3']).max())


def relax(output, model):
    folder = output / 'relax'
    folder.mkdir()
    setting = settings('relax')
    trials = {}

    def stress_at(a):
        if a not in trials:
            trial = folder / f'trial-{len(trials):02d}.json'
            print(f'G10_FLUORITE_RELAX_TRIAL a_A={a:.9f}', flush=True)
            row = calculate(fluoritic_atoms(a, 1), setting, trial.with_suffix('.scf'), stress=True)
            trials[a] = row = dict(row, a_A=a, settings=setting, model=model)
            commit_json(trial, row)
        return trials[a]['mean_stress_eV_A3']

    stress_at(SEED_A)
    lo, hi = SEED_A * 0.96, SEED_A * 1.04
    require(stress_at(lo) < 0 < stress_at(hi), 'No compressive-to-tensile bracket in the registered cubic window')
    root = brentq(stress_at, lo, hi, xtol=5e-5, rtol=1e-8, maxiter=12)
    stress_at(root)
    a = min(trials, key=lambda x: abs(trials[x]['mean_stress_eV_A3']))
    row = trials[a]
    require(stress_max(row) <= RELAX_STRESS and row['max_force_eV_A'] <= MAX_FORCE,
            'Cubic stress/force did not meet the reference tolerances; only raw trials are retained')
    reference = fluoritic_atoms(a, 2)
    commit_json(folder / 'reference.json', dict(
        model=model, a_A=a, primitive_final=row, sampling_reference=geometry(reference),
        pair=pair_geometry(reference, a), bracket_A=[lo, hi], root_xtol_A=5e-5,
        trial_sources=[binding(path) for path in sorted(folder.glob('trial-*.json'))],
        stress_tolerance_eV_A3=RELAX_STRESS, force_tolerance_eV_A=MAX_FORCE,
        relaxation='Only isotropic cubic a varied; Ti 4a and D 8c stay at exact fluorite Wyckoff positions.'))
    print(f'G10_FLUORITE_RELAX_DONE a_A={a:.9f}', flush=True)


def load_reference(output, model):
    path = output / 'relax' / 'reference.json'
    row = read_json(path)
    a = row['a_A']
    primitive = row['primitive_final']
    reference = fluoritic_atoms(a, 2)
    require(row['model'] == primitive['model'] == model and primitive['a_A'] == a and
            primitive['settings'] == settings('relax') and primitive['electron_count'] == 14 and
            stress_max(primitive) <= RELAX_STRESS and max_vector_norm(primitive['forces_eV_A']) <= MAX_FORCE and
            all(primitive[k] == v for k, v in geometry(fluoritic_atoms(a, 1)).items()) and
            row['sampling_reference'] == geometry(reference) and row['pair'] == pair_geometry(reference, a),
            'Saved cubic relaxation does not match this model and geometry')
    return row, binding(path)


def load_point(output, ref, ref_source, name):
    path = output / name / 'point.json'
    row = read_json(path)
    spec = case_specs()[name]
    expect = geometry(point_geometry(ref['a_A'], spec))
    require(row['case'] == name and row['spec'] == spec and row['reference'] == ref_source and
            row['settings'] == settings(spec['setting']) and all(row[k] == v for k, v in expect.items()),
            f'Saved point does not match the registered fixed geometry/settings: {name}')
    forces = np.asarray(row['forces_eV_A'])
    sigma = np.asarray(row['stress_eV_A3'])
    require(row['a_A'] == ref['a_A'] and row['electron_count'] == 112 and
            forces.shape == (24, 3) and np.isfinite(forces).all() and math.isfinite(row['energy_free_eV']) and
            (name not in ('zero', 'check-zero') or (sigma.shape == (3, 3) and np.isfinite(sigma).all())),
            f'Invalid saved native energy/force/stress vectors: {name}')
    return dict(row, source=binding(path))


def admit_zero(row):
    forces = np.asarray(row['forces_eV_A'])
    require(max_vector_norm(forces) <= MAX_FORCE and stress_max(row) <= 2 * RELAX_STRESS and
            np.linalg.norm(forces.sum(axis=0)) <= MAX_NET,
            'Point-group-free 24-atom zero failed the frozen reference force/stress diagnostic')


def run_point(output, ref, ref_source, name):
    spec = case_specs()[name]
    if name != 'zero':
        admit_zero(load_point(output, ref, ref_source, 'zero'))
    folder = output / name
    folder.mkdir()
    setting = settings(spec['setting'])
    print(f'G10_FLUORITE_POINT_START case={name}', flush=True)
    atoms = point_geometry(ref['a_A'], spec)
    row = calculate(atoms, setting, folder / 'point.scf', stress=name in ('zero', 'check-zero'))
    ray = ray_for(spec['direction']) if spec['moves_A'] else None
    projections = [] if ray is None else [[i, float(np.asarray(row['forces_eV_A'][i]) @ ray)]
                                           for i, _ in spec['moves_A']]
    commit_json(folder / 'point.json', dict(
        row, schema=SCHEMA, case=name, a_A=ref['a_A'], spec=spec, reference=ref_source,
        settings=setting, moved_site_force_parallel_eV_A=projections,
        net_force_within_target=bool(np.linalg.norm(row['net_force_eV_A']) <= MAX_NET)))
    if name == 'zero':
        admit_zero(row)
    print(f'G10_FLUORITE_POINT_DONE case={name}', flush=True)


def free_delta(point, zero):
    return math.fsum((point['energy_free_eV'], -zero['energy_free_eV']))


def pair_row(ref, load, zero, millia):
    a_row, b_row, ab_row = (load(name) for name in (f'100p{millia:03d}', f'pair-b-{millia:03d}', f'pair-ab-{millia:03d}'))
    f0, fa, fb, fab = (np.asarray(row['forces_eV_A']) for row in (zero, a_row, b_row, ab_row))
    pair = point_geometry(ref['a_A'], ab_row['spec'])
    separation = pair.get_distance(SITE_A, SITE_B, mic=True)
    require(abs(separation - (ref['pair']['d0_A'] - 2*millia/1000)) < 1e-12,
            'Measured pair geometry does not close along the registered nearest image')
    h = (np.arange(len(pair)) >= 8)
    min_ti_d = min(pair.get_distances(i, np.flatnonzero(~h), mic=True).min() for i in (SITE_A,SITE_B))
    ray = ray_for('100')
    dforces = (fab - fa) - (fb - f0)
    # The single-A term differentiates only on A; single-B only on B.
    dq = float(ray @ (fab[SITE_A] - fa[SITE_A] - fab[SITE_B] + fb[SITE_B]))
    return dict(displacement_per_site_A=millia/1000, separation_A=float(separation),
                min_moving_D_to_Ti_A=float(min_ti_d), atom_index=zero['atom_index'],
                controls={label: row['source'] for label, row in zip(('0', 'A', 'B', 'AB'),
                                                                    (zero, a_row, b_row, ab_row))},
                delta_energy_free_eV=math.fsum((ab_row['energy_free_eV'], -a_row['energy_free_eV'],
                                                -b_row['energy_free_eV'], zero['energy_free_eV'])),
                delta_forces_eV_A=dforces.tolist(),
                force_definition='F_AB(i) - F_A(i) - F_B(i) + F_0(i) at identical physical indices',
                delta_compression_path_force_eV_A=dq,
                path_force_definition='-d(delta energy)/du; [100] dot (F_AB(A) - F_A(A) - F_AB(B) + F_B(B))',
                all_control_net_forces_within_target=all(r['net_force_within_target'] for r in (zero, a_row, b_row, ab_row)))


def pair_readout(ref, load, zero):
    return dict(pair=ref['pair'], rows=[pair_row(ref,load,zero,u) for u in PAIR_MILLIANGSTROM])


def check_readout(load, zero):
    coarse, fine_zero, fine = (load(name) for name in (CHECK_CASE, 'check-zero', 'check-' + CHECK_CASE))
    coarse_delta, fine_delta = free_delta(coarse, zero), free_delta(fine, fine_zero)
    fs0, fs1, ft0, ft1 = (np.asarray(row['forces_eV_A']) for row in (zero, coarse, fine_zero, fine))
    df = (ft1 - ft0) - (fs1 - fs0)
    de = fine_delta - coarse_delta
    return dict(setting_change='Gamma-centered kmesh 4x4x4 -> 6x6x6 only, exactly fixed reference cell and positions',
                point=CHECK_CASE, coarse_delta_eV=coarse_delta, check_delta_eV=fine_delta,
                change_delta_eV=de, max_change_force_response_eV_A=max_vector_norm(df),
                max_change_raw_force_at_point_eV_A=max_vector_norm(ft1 - fs1),
                change_zero_stress_eV_A3=(np.asarray(fine_zero['stress_eV_A3']) - zero['stress_eV_A3']).tolist(),
                energy_target_eV=CHECK_E, force_target_eV_A=CHECK_F,
                tested_pair_within_targets=abs(de) <= CHECK_E and max_vector_norm(df) <= CHECK_F,
                limits='One two-point k-mesh comparison, not a cutoff, cell-size, stress, fit or whole-curve error bound.',
                point_sources=[r['source'] for r in (zero, coarse, fine_zero, fine)])


def cut_symmetries(a, direction):
    # Two-fold site rotations send +ray to -ray for these directions only.
    matrix = np.diag({'100': (-1, 1, -1), '110': (-1, -1, 1)}[direction])
    atoms = fluoritic_atoms(a, 2)
    positions = atoms.positions
    centre = positions[SITE_A]
    rotated = (positions - centre) @ matrix.T + centre
    delta = (rotated[:, None, :] - positions[None, :, :]) @ np.linalg.inv(atoms.cell.array)
    matches = (np.abs(delta - np.rint(delta)) < 1e-10).all(axis=2)
    matches &= atoms.numbers[:, None] == atoms.numbers[None, :]
    require((matches.sum(axis=0) == 1).all() and (matches.sum(axis=1) == 1).all() and
            np.array_equal(matrix @ ray_for(direction), -ray_for(direction)),
            'Registered site symmetry failed its geometry/permutation check')
    return matrix, matches.argmax(axis=1)


def cuts_readout(ref, load, zero):
    output = {}
    for direction in DIRECTIONS:
        ray = ray_for(direction)
        control = None if direction == '111' else cut_symmetries(ref['a_A'], direction)
        points = []
        for amp in MILLANGSTROM:
            minus, plus = (load(f'{direction}{sign}{amp:03d}') for sign in ('m', 'p'))
            signed = [dict(source=row['source'], displacement_A=row['spec']['moves_A'][0][1],
                           delta_energy_free_eV=free_delta(row, zero),
                           target_force_parallel_eV_A=float(np.asarray(row['forces_eV_A'][SITE_A]) @ ray),
                           net_force_within_target=row['net_force_within_target']) for row in (minus, plus)]
            diagnostic = None
            if control is not None:
                rotation, perm = control
                df = np.asarray(minus['forces_eV_A'])[perm] - np.asarray(plus['forces_eV_A']) @ rotation.T
                de = free_delta(plus, minus)
                mismatch = max_vector_norm(df)
                diagnostic = dict(plus_minus_energy_free_eV=de, maximum_rotated_force_mismatch_eV_A=mismatch,
                                  within_targets=abs(de) <= CHECK_E and mismatch <= CHECK_F)
            points.append(dict(signed_points=signed, twofold_symmetry_diagnostic=diagnostic))
        output[direction] = points
    return dict(zero=zero['source'], directions=output,
                signed_111='Independent inequivalent tetrahedral directions: +[111] faces three nearest Ti; -[111] points toward one Ti. No averaging.')


def plan(model):
    return dict(model=model, settings={key: settings(key) for key in ('relax', 'base', 'check')},
                relaxation_seed=geometry(fluoritic_atoms(SEED_A, 1)),
                sampling_seed=geometry(fluoritic_atoms(SEED_A, 2)),
                pair_at_seed=pair_geometry(fluoritic_atoms(SEED_A, 2), SEED_A),
                pair_controls=[{'0':'zero', 'A':f'100p{u:03d}', 'B':f'pair-b-{u:03d}', 'AB':f'pair-ab-{u:03d}'}
                               for u in PAIR_MILLIANGSTROM],
                cases=case_specs(),
                limits='Only volume is optimized within constrained cubic fluorite. Sample 2x2x2 fcc primitive replicas (24 atoms), not the dilute limit. Both signs are native SCFs. Change k mesh at exactly zero and one displaced point only.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True, help='Directory of exact pinned PDF and PAW bytes')
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('plan', help='Print registered plan without starting any GPAW calculation')
    p_relax = sub.add_parser('relax', help='Solve zero isotropic stress on the fluorite manifold')
    p_point = sub.add_parser('point', help='Fresh fixed-geometry SCF for each named case')
    p_readout = sub.add_parser('readout', help='Process saved native point JSON; no solver')
    for p in (p_relax, p_point, p_readout):
        p.add_argument('--output', type=Path, required=True, help='Dedicated calculation directory')
    p_point.add_argument('case', nargs='+', choices=tuple(case_specs()))
    p_readout.add_argument('section', choices=('pair', 'setting', 'cuts'))
    args = parser.parse_args()
    model = profile(args.source)
    if args.command == 'plan':
        print(json.dumps(plan(model), indent=2, allow_nan=False))
        return
    if args.command == 'relax':
        runtime(args.source)
        args.output.mkdir(parents=True, exist_ok=True)
        relax(args.output, model)
        return
    ref, ref_source = load_reference(args.output, model)
    if args.command == 'point':
        runtime(args.source)
        for name in args.case:
            run_point(args.output, ref, ref_source, name)
        return

    def load(name):
        return load_point(args.output, ref, ref_source, name)

    zero = load('zero')
    admit_zero(zero)
    if args.section == 'pair':
        result = pair_readout(ref, load, zero)
    elif args.section == 'setting':
        result = check_readout(load, zero)
    else:
        result = cuts_readout(ref, load, zero)
    commit_json(args.output / f'{args.section}.json', dict(result, schema=SCHEMA, model=model,
                                                           a_A=ref['a_A'], reference=ref_source))


if __name__ == '__main__':
    main()
