"""Conservative PAW valence/core separation; never rescale an electron density.

GPAW old/density.py explicitly warns that skip_core=True need not integrate
to the valence count: it omits the final per-atom monopole corrections. Here
we use the ordinary conserved full AE density, subtract a separately conserved
frozen core using the SAME upstream localized-function correction, and retain
original radial Pd core data for analytic electrostatics. No SCF is repeated.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from ase.units import Bohr
from gpaw import GPAW
from gpaw.lfc import LFC
from screen_io import write_json, write_npz
from well_native import bound


def conserved_valence(calc, atoms, gridrefinement=2):
    full_s, gd = calc.density.get_all_electron_density(
        atoms=atoms, gridrefinement=gridrefinement, skip_core=False)
    if gd.comm.size != 1:
        raise ValueError('The commissioned single-process density publisher is required')
    setups = calc.density.setups
    splines = []
    for setup in setups:
        nc = setup.nct if setup.type == 'all-electron' and setup.ni == 0 else setup.get_partial_waves()[2]
        splines.append([nc])
    positions = atoms.get_scaled_positions() % 1.
    localized = LFC(gd, splines)
    localized.set_positions(positions)
    atom_for_function = np.empty(len(localized.M_W), np.intc)
    start = 0
    for atom in localized.atom_indices:
        size = len(localized.sphere_a[atom].M_w)
        atom_for_function[start:start+size] = atom
        start += size
    core = gd.zeros()
    # This is the core-only instance of upstream density.py:574-590.
    residual_charge = -np.asarray([setup.Nc for setup in setups], float)
    localized.lfc.ae_core_density_correction(1., core, atom_for_function, residual_charge)
    gd.comm.sum(residual_charge)
    grid_positions = np.around(gd.N_c*positions).astype(int) % gd.N_c-gd.beg_c
    for residual, index in zip(residual_charge, grid_positions):
        if (index >= 0).all() and (index < gd.n_c).all():
            core[tuple(index)] -= residual/gd.dv
    full = full_s.sum(axis=0)
    valence = full-core
    counts = dict(full_AE=float(gd.integrate(full)), core=float(gd.integrate(core)),
                  valence=float(gd.integrate(valence)), solver_valence=float(calc.get_number_of_electrons()))
    if abs(counts['valence']-counts['solver_valence']) > 1e-7:
        raise ValueError(f'Conservative decomposition failed:{counts}')
    return valence/Bohr**3, counts


def publish_field(calc, atoms, output, gridrefinement=2):
    density, counts = conserved_valence(calc, atoms, gridrefinement)
    setup = calc.wfs.setups[next(i for i,a in enumerate(atoms) if a.symbol == 'Pd')]
    core = setup.data
    binding = write_npz(output, density_valence_e_A3=density,
        cell_A=atoms.cell.array, positions_A=atoms.positions, nuclear_Z=atoms.numbers,
        pd_nc_g=np.asarray(core.nc_g), pd_r_g_bohr=np.asarray(core.rgd.r_g),
        pd_dr_g_bohr=np.asarray(core.rgd.dr_g))
    return dict(binding=binding, electron_counts=counts, gridrefinement=gridrefinement,
        method='Conserved full AE minus independently conserved frozen cores using upstream PAW monopole correction',
        scope='Conservative grid distribution; pointwise density at underresolved core grid points is not asserted')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.process_time()
    spec = json.loads(args.input.read_bytes())
    args.output.mkdir(parents=True, exist_ok=False)
    calc = GPAW(str(bound(spec['gpw'])), txt=str(args.output/'restart-read.txt'), legacy_gpaw=True)
    atoms = calc.atoms.copy()
    rows = [publish_field(calc, atoms, args.output/f'field-g{ref}.npz', ref)
            for ref in spec['gridrefinements']]
    write_json(args.output/'fields.json', dict(schema='g8-conserved-valence-fields-v1',
        specification=spec, fields=rows, process_cpu_s=time.process_time()-started,
        scope='Cached exact same44-electron state; density reconstruction only, zero new SCFs'))
    print('G8_CONSERVED_'+'FIELDS_COMPLETE', flush=True)
