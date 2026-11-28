"""GPAW forces for one Phonopy PdD displacement at a registered electron temperature."""
import argparse
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
import os

import numpy as np
from ase.units import kB
from scipy.special import xlogy

from native import CONVERGENCE, CUTOFF_EV, require, sources, write_once
from phonon_geometry import CELLS, TEMPERATURES, electronic_mesh, inputs, phonon

SCHEMA = 'g16-pdd-phonon-native-v1'
CASES = ('zero','pd-p','d-p','pd-m','d-m')


def profile(natoms, temperature, *, compact_bands=False, tight_forces=False):
    per_formula = 10 if temperature <= 3000 else 16 if temperature == 10000 else 32
    compact_bands = compact_bands or natoms == 128
    tight_forces = tight_forces or natoms == 128
    convergence = CONVERGENCE
    if tight_forces:
        require(temperature in (300, 3000), 'Tight force provision is registered only for the size anchors')
        convergence = {**CONVERGENCE, 'forces': 2e-8, 'density': 2e-8}
    if compact_bands:
        require(temperature in (300, 3000), 'Compact band provision is registered only for the size anchors')
        per_formula = 7
    cfg = dict(cutoff_eV=CUTOFF_EV, xc='PBE', charge=0, spinpol=False,
               nbands=natoms//2 * per_formula, valence_electrons=11*natoms//2,
               maxiter=120, convergence=convergence, symmetry='actual displaced crystal',
               smearing_eV=float(kB*temperature))
    if compact_bands:
        cfg['occupation_tail_max'] = 1e-10
    return cfg


def capture(atoms, calc, spec):
    forces = atoms.get_forces()
    free = float(atoms.get_potential_energy(force_consistent=True))
    extrap = float(atoms.get_potential_energy())
    cfg = spec['calculator']
    w = calc.dft.ibzwfs
    _, occ = w.get_all_eigs_and_occs(broadcast=True)
    weight = np.asarray(calc.get_k_point_weights())
    f = occ[0]
    n_e = float(2*np.dot(f.sum(axis=1),weight))
    tail = float(f[:,-1].max())
    neg_ts = float(2*cfg['smearing_eV']*np.dot((xlogy(f,f)+xlogy(1-f,1-f)).sum(axis=1),weight))
    bz = calc.get_bz_k_points()
    ibz = calc.get_ibz_k_points()
    require(w.nspins==1 and w.spin_degeneracy==2 and np.array_equal(bz,spec['electronic_grid'])
            and np.isfinite(forces).all() and -1e-14<=f.min()<=f.max()<=1
            and 0<=tail<=cfg.get('occupation_tail_max', 1e-6)
            and abs(n_e-cfg['valence_electrons'])<1e-8 and calc.get_number_of_electrons()==cfg['valence_electrons']
            and calc.get_number_of_bands()==cfg['nbands'] and cfg['maxiter']>=calc.get_number_of_iterations()>0
            and abs(extrap-(free-.5*neg_ts))<1e-8 and np.array_equal(calc.atoms.positions,atoms.positions)
            and calc.dft.scf_loop.update_density_and_potential,
            f'Native PdD result unqualified (band tail {tail:.3g}, electron count {n_e:.9f})')
    return dict(schema=SCHEMA, specification=spec, positions_A=atoms.positions.tolist(), cell_A=atoms.cell.array.tolist(),
                symbols=atoms.get_chemical_symbols(), masses_u=atoms.get_masses().tolist(),
                forces_eV_A=forces.tolist(), energy_free_eV=free, energy_internal_eV=free-neg_ts,
                minus_T_S_eV=neg_ts, energy_extrapolated_eV=extrap, occupation_tail=tail,
                electron_count=n_e, scf_iterations=int(calc.get_number_of_iterations()),
                ibz_points=ibz.tolist(), ibz_weights=weight.tolist())


def run(a):
    import gpaw
    from gpaw import FermiDirac, Mixer, PW
    from gpaw.mpi import world
    from gpaw.new.ase_interface import GPAW

    require(version('gpaw')=='26.7.0' and version('ase')=='3.29.0' and version('phonopy')=='4.11.0', 'Wrong scientific runtime')
    cores = tuple(int(s) for s in a.cores.split(','))
    require(len(cores)==len(set(cores))==world.size and set(cores)<=os.sched_getaffinity(0), 'Invalid requested affinity')
    os.sched_setaffinity(0,{cores[world.rank]})
    gpaw.setup_paths[:] = [str(a.paw.resolve())]
    atoms, atom, vector = inputs(phonon(a.cell))[a.case]
    cfg = profile(len(atoms),a.temperature_K,compact_bands=a.compact_bands,tight_forces=a.tight_forces)
    grid = electronic_mesh(a.mesh,a.cell)
    source_dir = Path(__file__).resolve().parent
    spec = dict(case=a.case, temperature_K=a.temperature_K, primitive_mesh_n=a.mesh, cell=a.cell,
                supercell_rows=CELLS[a.cell].tolist(), displaced_atom=atom,
                displacement_A=None if vector is None else vector.tolist(), electronic_grid=grid.tolist(),
                calculator=cfg, PAW=sources(a.paw),
                source_sha256={f:sha256((source_dir/f).read_bytes()).hexdigest()
                               for f in ('native.py','phonon_geometry.py','phonon_native.py')})
    out = a.output
    if world.rank==0:
        out.mkdir(parents=True,exist_ok=False)
        write_once(out/'specification.json',spec)
        print(f'G16_PHONON_START {out.name}',flush=True)
    world.barrier()
    atoms.calc = calc = GPAW(mode=PW(cfg['cutoff_eV']), xc=cfg['xc'], setups={'Pd':'10'}, charge=0, spinpol=False,
                            kpts=grid, nbands=cfg['nbands'], occupations=FermiDirac(cfg['smearing_eV']),
                            symmetry={'point_group':True}, eigensolver='davidson', convergence=cfg['convergence'],
                            parallel={'kpt':world.size,'domain':1,'band':1}, mixer=Mixer(.08,6,50),
                            maxiter=cfg['maxiter'], txt=str(out/'scf.log'), communicator=world)
    row = capture(atoms,calc,spec)
    if world.rank==0:
        write_once(out/'point.json',row)
        print(f'G16_PHONON_DONE {out.name}',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--temperature-K',type=int,choices=TEMPERATURES,required=True)
    p.add_argument('--mesh',type=int,choices=(12,16),required=True)
    p.add_argument('--cell',choices=CELLS,required=True)
    p.add_argument('--case',choices=CASES,required=True)
    p.add_argument('--cores',required=True)
    p.add_argument('--compact-bands',action='store_true',help='Use the seven-band anchor control')
    p.add_argument('--tight-forces',action='store_true',help='Use the tighter force stopping control')
    p.add_argument('--paw',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args())


if __name__=='__main__':
    main()
