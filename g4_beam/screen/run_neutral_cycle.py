"""Neutral DD/He fixed-ion endpoint energies; no relaxation or potential surrogate."""
import argparse
import gc
import json
from pathlib import Path
import time

import gpaw
import numpy as np
from ase import Atoms
from gpaw import FermiDirac, GPAW, PW
from gpaw.mpi import serial_comm
from screen_io import manifest, require, write_json


def endpoint(config, output):
    start = time.monotonic()
    seed = config['geometry']
    atoms = Atoms(seed['symbols'], positions=seed['positions_A'], cell=seed['cell_A'], pbc=seed['pbc'])
    atoms.set_initial_magnetic_moments(seed['initial_moments_muB'])
    spin = float(sum(seed['initial_moments_muB']))
    bands = int(np.ceil((config['valence_electrons']+abs(spin))/2))+8
    width = .01 if config['host']=='vacuum' else .1
    calc = GPAW(legacy_gpaw=False, communicator=serial_comm, mode=PW(500), xc='PBE',
                setups={'Pd':'10', 'Rh':'9'},
                kpts={'size':config['kmesh'], 'gamma':True},
                occupations=FermiDirac(width, fixmagmom=config['spinpol']),
                spinpol=config['spinpol'], nbands=bands, symmetry={'point_group':False},
                convergence={'energy':1e-6, 'density':1e-5, 'eigenstates':1e-8},
                maxiter=150, parallel={'domain':1, 'band':1, 'kpt':1}, txt=str(output/'ground.txt'))
    atoms.calc = calc
    energy = float(atoms.get_potential_energy())
    free_energy = float(atoms.get_potential_energy(force_consistent=True))
    actual_electrons = float(calc.get_number_of_electrons())
    require(actual_electrons == config['valence_electrons'], 'Wrong neutral electron count')
    require(np.isfinite(energy) and np.isfinite(free_energy), 'Nonfinite endpoint energy')
    forces = atoms.get_forces().tolist()
    gpw = output/'ground.gpw'; calc.write(str(gpw), mode='all')
    result = dict(case=config['case'], config=config,
                  energy_extrapolated_zero_temperature_eV=energy, free_energy_eV=free_energy,
                  actual_neutral_valence_electrons=actual_electrons,
                  occupation_width_eV=width, PW_ecut_eV=500., scf_nbands=bands,
                  forces_eV_A=forces, elapsed_wall_s=time.monotonic()-start,
                  gpaw_version=gpaw.__version__, gpw=manifest(gpw),
                  atomic_reference_energy_eV=float(calc.get_reference_energy()),
                  ion_scope='All ions fixed; only electrons relax. D represented by H in the Born-Oppenheimer electronic Hamiltonian.')
    write_json(output/'energy.json', result)
    atoms.calc=None; del calc; gc.collect()
    return result


def reference_calibration(atoms, sources):
    require(set(atoms)=={'H','He'}, 'Neutral calibration needs both isolated atoms')
    refs={s:r['atomic_reference_energy_eV'] for s,r in atoms.items()}
    absolute={s:r['energy_extrapolated_zero_temperature_eV']+refs[s] for s,r in atoms.items()}
    return dict(schema='coherence-neutral-PAW-reference-calibration-v1',
                atomic_reference_energy_eV=refs,
                reference_offset_2H_minus_He_eV=2*refs['H']-refs['He'],
                neutral_atom_absolute_energy_eV=absolute,
                neutral_vacuum_adiabatic_difference_eV=2*absolute['H']-absolute['He'],
                sources=sources, gpaw_version=gpaw.__version__,
                convention='Actual GPAW get_reference_energy() added to reported relative-to-atom energies; no experimental constant substituted')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    config=json.loads(args.input.read_bytes())
    require(config['schema']=='coherence-neutral-cycle-plan-v1','Wrong scientific input role')
    sources=[manifest(path, 'source/'+path.name) for path in sorted((args.output/'source').glob('*.py'))]
    write_json(args.output/'producer.json',dict(input=manifest(args.input), producer_source=sources,
              setup_bytes=config['setup_bytes'], observable=config['observable']))
    records=[];atomic={};atomic_sources=[]
    for case in config['cases']:
        folder=args.output/case['case'];folder.mkdir()
        write_json(folder/'input.json',case)
        result=endpoint(case,folder)
        records.append(dict(case=case['case'],energy=manifest(folder/'energy.json',case['case']+'/energy.json')))
        if case['host']=='vacuum' and case['vacuum_box_A']==10.:
            atomic[case['endpoint']]=result
            atomic_sources.append(dict(symbol=case['endpoint'],gpw=result['gpw'],
                                       energy=manifest(folder/'energy.json')))
        print('NEUTRAL_ENDPOINT_COMPLETE '+case['case'],flush=True)
        require(result['elapsed_wall_s'] <= 3*case['expected_wall_s'],
                'Threefold endpoint-cost deviation; stop before further production')
    calibration=None
    if atomic:
        calibration=write_json(args.output/'reference-calibration.json',
                               reference_calibration(atomic,atomic_sources))
    write_json(args.output/'completed.json',dict(schema='coherence-neutral-cycle-completed-v1',
               producer=manifest(args.output/'producer.json'),records=records,
               reference_calibration=calibration))
    print('NEUTRAL_CYCLE_COMPLETE',flush=True)


if __name__=='__main__':
    main()
