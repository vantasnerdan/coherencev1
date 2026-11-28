"""Native fixed-ion thermal electronic PdD points on the registered ladder."""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path

import numpy as np
from ase.build import bulk
from ase.units import GPa, Hartree, kB

A = 4.084
D_MASS_U = 2.01410177812
TEMPERATURES_K = (300, 310, 350, 1000, 3000, 10000, 30000)
AXES = {'zero': None, '100': (1., 0., 0.), '111': (1 / np.sqrt(3),) * 3}
PAW = {
    'Pd.10.PBE.gz': (106793, '3c6cfd1de7e6d6815d58c9521b2bc7cec27ede8e83554ed1f0151f979761bbb4'),
    'H.PBE.gz': (12501, '3136ad73c02d471f99dd5f82e451b2061952970556a260d38b53b5020dc94ae7'),
}
CONVERGENCE = dict(energy=1e-8, density=2e-7, eigenstates=1e-8, forces=2e-6, bands='occupied')
CUTOFF_EV = 500
BANDS = 32


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def write_once(path, result):
    temp = path.with_suffix('.tmp')
    with temp.open('x') as out:
        json.dump(result, out, indent=2, allow_nan=False)
        out.write('\n')
    temp.rename(path)


def reference():
    atoms = bulk('PdH', 'rocksalt', a=A)
    require(atoms.get_chemical_symbols() == ['Pd', 'H'], 'Unexpected atom order')
    atoms.set_masses([atoms.get_masses()[0], D_MASS_U])
    return atoms


def sources(paw):
    pins = {}
    for name, (size, dig) in PAW.items():
        path = paw / name
        with path.open('rb') as stream:
            sha = hashlib.file_digest(stream, 'sha256').hexdigest()
        require((path.stat().st_size, sha) == (size, dig), f'Incorrect PAW bytes for {name}')
        pins[name] = dict(bytes=size, sha256=sha)
    return pins


def native_record(atoms, calc, spec, direction, u):
    free = float(atoms.get_potential_energy(force_consistent=True))
    forces = atoms.get_forces()
    # GPAW's `energy` is its smearing extrapolation, not U at the given Te.
    neg_ts = float(calc.hamiltonian.e_entropy * Hartree)
    extrap = float(atoms.get_potential_energy())
    internal = free - neg_ts
    require(neg_ts <= 0 and abs(extrap - (free - 0.5 * neg_ts)) < 1e-9, 'Finite-temperature energy convention differs')
    w = calc.wfs
    top = max(float(q.f_n[-1] / (2 * w.kd.weight_k[q.k])) for q in w.kpt_u)
    sum_occ = sum(float(q.f_n.sum()) for q in w.kpt_u)
    require(calc.scf.converged and w.world.size == 1 and calc.get_number_of_electrons() == 11
            and abs(sum_occ - 11) < 1e-8 and len(w.kpt_u) == len(w.kd.weight_k) and calc.get_number_of_bands() == BANDS
            and top <= 1e-6 and np.isfinite(forces).all(), f'Native convergence, occupations or geometry failed (edge f={top})')
    stress = atoms.get_stress(voigt=False) if direction == 'zero' else None
    return dict(schema='g16-pdd-thermal-point-v1', specification=spec,
                cell_A=atoms.cell.array.tolist(), symbols=atoms.get_chemical_symbols(), masses_u=atoms.get_masses().tolist(),
                positions_A=atoms.positions.tolist(), volume_A3=float(atoms.get_volume()), direction=direction, displacement_A=u,
                energy_free_eV=free, minus_T_S_eV=neg_ts, energy_internal_eV=internal, energy_extrapolated_eV=extrap,
                forces_eV_A=forces.tolist(), stress_eV_A3=None if stress is None else stress.tolist(),
                pressure_GPa=None if stress is None else float(-np.trace(stress) / (3 * GPa)),
                electron_count=sum_occ, maximum_last_band_occupation=top, n_ibz_kpoints=len(w.kpt_u),
                scf_iterations=int(calc.get_number_of_iterations()))


def run(args):
    import gpaw
    from gpaw import FermiDirac, Mixer, PW
    from gpaw.old.calculator import GPAW
    from gpaw.mpi import world

    require(world.size == 1 and version('gpaw') == '26.7.0' and version('ase') == '3.29.0', 'Wrong GPAW/ASE execution profile')
    paw = sources(args.paw)
    gpaw.setup_paths[:] = [str(args.paw.resolve())]
    atoms = reference()
    xyz = atoms.positions.copy()
    axis = AXES[args.direction]
    if axis is None:
        require(args.displacements == [0.], 'Zero geometry must be isolated from displaced symmetry')
    else:
        require(all(0 < abs(b) <= 0.25 for b in args.displacements) and len(set(args.displacements)) == len(args.displacements),
                'Displaced direction needs nonzero unique shifts')
        xyz[1] += args.displacements[0] * np.asarray(axis)
        atoms.positions[:] = xyz
        xyz = reference().positions
    base = args.output
    base.mkdir(parents=True, exist_ok=False)
    width = float(kB * args.temperature_K)
    spec = dict(temperature_K=args.temperature_K, smearing_eV=width, kmesh=[args.mesh] * 3, gamma=True,
                cutoff_eV=CUTOFF_EV, PAW=paw, xc='PBE', nbands=BANDS, charge=0, spinpol=False,
                symmetry='actual displaced structure', convergence=CONVERGENCE,
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    write_once(base / 'specification.json', spec)
    atoms.calc = GPAW(mode=PW(CUTOFF_EV), setups={'Pd': '10'}, xc='PBE', charge=0, spinpol=False,
                      kpts={'size': [args.mesh] * 3, 'gamma': True}, nbands=BANDS,
                      occupations=FermiDirac(width), symmetry={'point_group': True}, mixer=Mixer(0.05, 5, 50),
                      convergence=CONVERGENCE, maxiter=120, txt=str(base / 'scf.log'))
    for u in args.displacements:
        if axis is not None:
            atoms.positions[1] = xyz[1] + u * np.asarray(axis)
        label = f'{u:+.4f}'
        print(f'G16_NATIVE_START k={args.mesh} K={args.temperature_K} {args.direction} {label}', flush=True)
        result = native_record(atoms, atoms.calc, spec, args.direction, u)
        write_once(base / f'u{label}.json', result)
        print(f'G16_NATIVE_DONE k={args.mesh} K={args.temperature_K} {args.direction} {label}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--paw', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--temperature-K', type=int, choices=TEMPERATURES_K, required=True)
    p.add_argument('--mesh', type=int, choices=(12, 16), required=True)
    p.add_argument('--direction', choices=AXES, required=True)
    p.add_argument('displacements', type=float, nargs='+')
    run(p.parse_args())


if __name__ == '__main__':
    main()
