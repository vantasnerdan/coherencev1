"""Matched D0/D3 displacements from the exact relaxed 64-site Pd vacancy."""
import argparse
import json
from pathlib import Path

from ase.geometry import find_mic
import numpy as np

from pd_multivac import case_input, commit_json, geometry, read_record, require, root_action, publish

SCHEMA = 'g10-pd-vacancy-r4-d6-pair-v1'
REF_SHA = '0bda3d044fdb291542b0283e299bf1d0e93c3bb25642f0def6f11baf8a9ac888'
I, J = 63, 66
AMPS = (6, 12, 24)
ROLES = ('A', 'B', 'AB')


def source(args):
    atoms, _, spec, _ = case_input(args)
    ref, bind = read_record(args.reference)
    identity = dict(file=args.reference.name, bytes=bind['bytes'], sha256=bind['sha256'])
    require(bind['sha256'] == REF_SHA and ref['specification'] == json.loads(json.dumps(spec))
            and ref['ionic_converged'] and ref['distances']['all_trapped_near_centre']
            and ref['raw_fmax_eV_A'] <= 0.05 and ref['scf_converged'], 'Wrong resting vacancy reference')
    atoms.positions[:] = ref['positions_A']
    return atoms, spec['calculator'], identity


def displaced(atoms, role, amp):
    require((role == 'zero' and amp == 0) or (role in ROLES and amp in AMPS), 'Unregistered displacement')
    vec, r0 = find_mic(atoms.positions[J] - atoms.positions[I], atoms.cell, atoms.pbc)
    require(abs(r0 - 1.966074631104379) < 1e-12, 'Wrong registered pair/image')
    n = vec/r0
    moved = atoms.copy()
    if 'A' in role:
        moved.positions[I] += (amp/100)*n
    if 'B' in role:
        moved.positions[J] -= (amp/100)*n
    min_pd = min(min(moved.get_distances(i, range(I), mic=True)) for i in (I, J))
    r = moved.get_distance(I, J, mic=True)
    require(abs(r - (r0 - len(role)*(amp/100))) < 1e-12, 'Pair separation does not close')
    path = dict(role=role, displacement_per_D_A=amp/100, tags=[I, J], axis_A_to_B=n.tolist(),
                separation_A=float(r), minimum_tagged_Pd_D_A=float(min_pd),
                nominal_Pd_H_overlap=bool(min_pd < 3.4*0.529177210903))
    return moved, path


def calculate(args, atoms, cfg, identity):
    import gpaw
    from gpaw import FermiDirac, Mixer, PW
    from gpaw.mpi import world
    from gpaw.new.ase_interface import GPAW
    from loading_pd import bind_ranks

    bind_ranks(world, args.ranks)
    gpaw.setup_paths[:] = [str(args.paw.resolve())]
    moved, path = displaced(atoms, args.role, args.amp)
    nk = 8
    if args.mesh == 3:
        require(args.amp in (0, max(AMPS)), 'Unregistered fine mesh check location')
        cfg['kmesh'] = [3, 3, 3]
        # Time reversal only: 27 full, 14 irreducible k points. No spatial symmetry.
        cfg['symmetry'] = {'point_group': False, 'time_reversal': True}
        nk = 14
    require(nk % world.size == 0, 'Use ranks dividing the actual irreducible mesh')
    out = args.output
    log = out.with_suffix('.scf')
    def prepare():
        out.parent.mkdir(parents=True, exist_ok=True)
        require(not out.exists() and not log.exists(), 'Keep previous attempt at its original path')
    root_action(world, prepare)
    keys = ('xc','convergence','setups','charge','spinpol','eigensolver','symmetry','maxiter','nbands')
    moved.calc = calc = GPAW(**{k:cfg[k] for k in keys}, mode=PW(cfg['cutoff_eV']),
                            kpts={'size':cfg['kmesh'],'gamma':cfg['gamma']}, occupations=FermiDirac(cfg['smear_eV']),
                            mixer=Mixer(0.04, 7, 80), parallel={'kpt':world.size,'domain':1,'band':1},
                            communicator=world, txt=str(log))
    energy = float(moved.get_potential_energy(force_consistent=True))
    force = moved.get_forces(apply_constraint=False)
    stats = calc.dft.ibzwfs
    kpts = calc.get_ibz_k_points()
    require(np.isfinite(energy) and force.shape == (69,3) and np.isfinite(force).all()
            and np.array_equal(calc.atoms.positions, moved.positions) and len(kpts) == nk
            and calc.get_number_of_electrons() == cfg['expected_electrons']
            and calc.get_number_of_bands() == cfg['nbands'] and stats.kpt_comm.size == world.size
            and stats.domain_comm.size == stats.band_comm.size == 1
            and calc.dft.scf_loop.update_density_and_potential, 'Wrong native geometry, SCF or topology')
    row = dict(**geometry(moved), schema=SCHEMA, reference=identity, calculator=cfg, point=path,
               energy_free_eV=energy, forces_eV_A=force.tolist(), scf_converged=True,
               scf_iterations=calc.get_number_of_iterations(), ibz_kpoints=kpts.tolist(), mpi_ranks=world.size)
    publish(world, commit_json, out, row)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--paw', type=Path, required=True)
    p.add_argument('--articles', type=Path, required=True)
    p.set_defaults(repeat=4, occupancy=6)
    sub = p.add_subparsers(dest='action', required=True)
    sub.add_parser('plan')
    run = sub.add_parser('point')
    run.add_argument('role', choices=('zero', *ROLES))
    run.add_argument('amp', type=int)
    run.add_argument('--mesh', type=int, choices=(2,3), required=True)
    run.add_argument('--ranks', type=int, required=True)
    run.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    atoms, cfg, identity = source(args)
    if args.action == 'plan':
        print(json.dumps([displaced(atoms, r, u)[1] for u in AMPS for r in ROLES], indent=2))
    else:
        calculate(args, atoms, cfg, identity)


if __name__ == '__main__':
    main()
