"""Native signed Γ or X PdD force probes on a matched folded electronic grid."""
import argparse
import hashlib
from importlib.metadata import version
from pathlib import Path

import numpy as np
from ase.units import Hartree, kB

from fold import SUPERCELL, geometry, mesh
from native import BANDS, CONVERGENCE, CUTOFF_EV, require, sources, write_once


def native_record(atoms, calc, spec, u):
    f = float(atoms.get_potential_energy(force_consistent=True))
    force = atoms.get_forces()
    negative_ts = float(calc.hamiltonian.e_entropy * Hartree)
    extrapolated = float(atoms.get_potential_energy())
    u_electronic = f - negative_ts
    require(negative_ts <= 0 and abs(extrapolated - (f - 0.5 * negative_ts)) < 1e-9,
            'Incorrect electronic energy decomposition')
    w = calc.wfs
    tail = max(float(q.f_n[-1] / (2 * w.kd.weight_k[q.k])) for q in w.kpt_u)
    electron_count = sum(float(q.f_n.sum()) for q in w.kpt_u)
    require(calc.scf.converged and w.world.size == 1 and calc.get_number_of_electrons() == 22
            and abs(electron_count - 22) < 1e-8 and len(w.kpt_u) == len(w.kd.weight_k)
            and calc.get_number_of_bands() == 2 * BANDS and tail <= 1e-6
            and np.isfinite(force).all() and np.array_equal(atoms.positions, calc.atoms.positions),
            'Folded electronic state differs or has an occupied band edge')
    return dict(schema='g16-pdd-folded-point-v1', specification=spec, displacement_A=u,
                positions_A=atoms.positions.tolist(), cell_A=atoms.cell.array.tolist(), symbols=atoms.get_chemical_symbols(),
                masses_u=atoms.get_masses().tolist(), volume_A3=float(atoms.get_volume()),
                energy_free_eV=f, energy_internal_eV=u_electronic, minus_T_S_eV=negative_ts,
                energy_extrapolated_eV=extrapolated, forces_eV_A=force.tolist(), scf_iterations=calc.get_number_of_iterations(),
                electron_count=electron_count, maximum_last_band_occupation=tail,
                bz_kpoints_count=w.kd.nbzkpts, ibz_kpoints_count=len(w.kpt_u))


def run(args):
    import gpaw
    from gpaw import FermiDirac, Mixer, PW
    from gpaw.mpi import world
    from gpaw.old.calculator import GPAW

    require(world.size == 1 and version('gpaw') == '26.7.0' and version('ase') == '3.29.0', 'Unexpected native environment')
    gpaw.setup_paths[:] = [str(args.paw.resolve())]
    paw = sources(args.paw)
    atoms, sites, signs, axis = geometry(args.mode, args.polarization)
    grid = mesh(args.mesh)
    ideal = atoms.positions.copy()
    width = float(kB * args.temperature_K)
    spec = dict(temperature_K=args.temperature_K, smearing_eV=width, mode=args.mode, polarization=args.polarization,
                primitive_gamma_mesh=[args.mesh] * 3, folded_electronic_kpoints_fractional=grid.tolist(),
                primitive_to_four_atom_integer_rows=SUPERCELL.tolist(), D_indices=sites.tolist(), D_displacement_signs=signs.tolist(),
                polarization_cartesian=axis.tolist(), cutoff_eV=CUTOFF_EV, nbands=2 * BANDS, xc='PBE',
                PAW=paw, net_charge=0, spinpol=False, convergence=CONVERGENCE,
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    write_once(out / 'specification.json', spec)
    steps = [0.] if args.zero else [-args.step_A, args.step_A]
    if not args.zero:
        atoms.positions[sites] = ideal[sites] + (steps[0] * signs)[:, None] * axis
    atoms.calc = GPAW(mode=PW(CUTOFF_EV), setups={'Pd': '10'}, xc='PBE', charge=0, spinpol=False,
                      kpts=grid, nbands=2 * BANDS, occupations=FermiDirac(width), symmetry={'point_group': True},
                      mixer=Mixer(0.05, 5, 50), convergence=CONVERGENCE, maxiter=120, txt=str(out / 'scf.log'))
    for u in steps:
        atoms.positions[sites] = ideal[sites] + (u * signs)[:, None] * axis
        label = f'{u:+.4f}'
        print(f'G16_FOLD_START k={args.mesh} K={args.temperature_K} {args.mode}/{args.polarization} {label}', flush=True)
        write_once(out / f'u{label}.json', native_record(atoms, atoms.calc, spec, u))
        print(f'G16_FOLD_DONE k={args.mesh} K={args.temperature_K} {args.mode}/{args.polarization} {label}', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--paw', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--temperature-K', type=int, choices=(300, 3000, 10000), required=True)
    p.add_argument('--mesh', type=int, choices=(12, 16), required=True)
    p.add_argument('--mode', choices=('gamma', 'x'), required=True)
    p.add_argument('--polarization', choices=('100', '010'), required=True)
    p.add_argument('--step-A', type=float, choices=(0.01, 0.02), default=0.01)
    p.add_argument('--zero', action='store_true')
    a = p.parse_args()
    require(not a.zero or (a.mode, a.polarization) == ('gamma', '100'), 'Only one zero geometry is used')
    run(a)


if __name__ == '__main__':
    main()
