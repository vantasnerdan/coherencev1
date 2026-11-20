"""Propagate the electronic state at prescribed midpoint ionic positions.

The moving regularized all-electron H has no PAW projectors. All Pd projector
centres are frozen, so S is constant and the moving-projector term is zero.
"""
import json
import time
import numpy as np
from ase.units import Hartree, Bohr
from gpaw.tddft import TDDFT
from gpaw.tddft.units import attosec_to_autime
from gpaw.old.forces import calculate_forces
from screen_io import manifest, write_json, write_npz
from substrate_framework.collision_observables import (
    PeriodicElectronPotential, deuteron_speed_A_as, external_work_eV)
from collision_geometry import projectile_point
from substrate_framework.fixed_projector_contractions import install_fixed_projectors


def _position(calc, projectile, point):
    calc.atoms.positions[projectile] = point
    calc.set_positions(calc.atoms)


def _sample(calc, potential, projectile, R, wall_start):
    energy = float(calc.get_td_energy() * Hartree)
    forces = calculate_forces(calc.wfs, calc.density, calc.hamiltonian, calc.log)
    density = calc.get_all_electron_density(gridrefinement=2)
    if potential is None:
        potential = PeriodicElectronPotential(density.shape, calc.atoms.cell.array)
    point = calc.atoms.positions[projectile].copy()
    phi = float(potential.samples(density, [point])[0])
    row = dict(R_A=float(R), time_as=float(calc.time / attosec_to_autime),
               projectile_A=point.tolist(), projectile_force_eV_A=(forces[projectile]*Hartree/Bohr).tolist(),
               td_internal_energy_eV=energy, electron_potential_eV=phi,
               integrated_all_electrons=potential.electrons(density),
               elapsed_wall_s=time.monotonic()-wall_start)
    return row, density, potential


def propagate(spec, ground_dir, output):
    output.mkdir(parents=True, exist_ok=True)
    ground = json.loads((ground_dir/'ground.json').read_bytes())
    calc = TDDFT(str(ground_dir/'ground.gpw'), propagator={'name': 'SICN'},
                 solver={'name': 'BiCGStab', 'tolerance': spec['solver_tolerance']},
                 parallel={'domain': 1, 'band': 1, 'kpt': 1},
                 txt=str(output/'td.txt'))
    calc.tddft_init()
    projectile = ground['path']['projectile_index']
    if calc.wfs.setups[projectile].ni != 0:
        raise ValueError('Prescribed mover has PAW projectors; constant-overlap propagation is not applicable')
    contractions = install_fixed_projectors(calc)
    write_json(output/'projector-equivalence.json', contractions)
    speed = deuteron_speed_A_as(spec['energy_keV'])
    R = float(spec['start_R_A'])
    end = float(spec['end_R_A'])
    knots = sorted(set([float(x) for x in spec['checkpoints_A']] + [R, end]), reverse=True)
    knots = [x for x in knots if end <= x <= R]
    rows, checkpoints = [], []
    potential = None
    start = time.monotonic()
    step = 0
    for knot_index, knot in enumerate(knots):
        if knot_index:
            while R > knot + 1e-12:
                dt = min(spec['dt_as'], (R-knot)/speed)
                midpoint = R-0.5*speed*dt
                endpoint = R-speed*dt
                _position(calc, projectile, projectile_point(ground['path'], midpoint, spec['impact_A']))
                calc.propagate(dt, 1)
                R = endpoint
                _position(calc, projectile, projectile_point(ground['path'], R, spec['impact_A']))
                row, density, potential = _sample(calc, potential, projectile, R, start)
                rows.append(row)
                step += 1
                if step % 10 == 0:
                    print(f'COLLISION_TD_POINT step={step} R_A={R:.8g} wall_s={row["elapsed_wall_s"]:.4g}', flush=True)
                if len(rows) > 11:
                    recent = rows[-1]['elapsed_wall_s']-rows[-11]['elapsed_wall_s']
                    expected = spec.get('maximum_ten_step_wall_s')
                    if expected is not None and recent > expected:
                        raise RuntimeError('Propagation cost trigger reached; stop before scaling the trajectory')
        else:
            row, density, potential = _sample(calc, potential, projectile, R, start)
            rows.append(row)
            mismatch = abs(row['td_internal_energy_eV']-ground['energy']['internal_energy_eV'])
            if mismatch > 1e-4:
                raise ValueError(f'Initial BO/TD internal-energy convention mismatch: {mismatch} eV')
        tag = f'checkpoint-{knot_index:02d}'
        density_receipt = write_npz(output/(tag+'-density.npz'), density_e_A3=density,
                                    cell_A=calc.atoms.cell.array, projectile_A=calc.atoms.positions[projectile])
        gpw = output/(tag+'.gpw')
        calc.write(str(gpw), mode='all')
        checkpoint = dict(index=knot_index, row_index=len(rows)-1, R_A=float(knot),
                          gpw=manifest(gpw), density=density_receipt,
                          positions_A=calc.atoms.positions.tolist())
        write_json(output/(tag+'.json'), checkpoint)
        checkpoints.append(checkpoint)
        write_json(output/f'progress-{knot_index:02d}.json', dict(rows=rows, checkpoints=checkpoints))
        print(f'COLLISION_CHECKPOINT index={knot_index} R_A={knot:.8g}', flush=True)
    work = external_work_eV([r['projectile_A'] for r in rows],
                            [r['projectile_force_eV_A'] for r in rows])
    for row, value in zip(rows, work):
        row['external_work_eV'] = float(value)
    result = dict(schema='pdd-collision-trajectory-v1', specification=spec,
                  ground=manifest(ground_dir/'ground.json'), rows=rows, checkpoints=checkpoints,
                  gauge='Electron-only periodic Poisson potential: zero G=0 mode; units eV on positive unit charge',
                  projector_motion='Only projector-free all-electron H moves; all PAW projectors remain fixed',
                  fixed_projector_contractions=contractions,
                  elapsed_wall_s=time.monotonic()-start)
    write_json(output/'trajectory.json', result)
    print('COLLISION_TRAJECTORY_COMPLETE', flush=True)
    return result
