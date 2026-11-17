"""Serial climbing-image NEB, explicit endpoints, complete band-force history."""
import time
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from ase.mep import NEB
from ase.optimize import FIRE
from surface_calculator import calculator
from surface_io import save, geometry, manifest


def run_neb(initial_dir, final_dir, out, interior=5, steps=45, fmax=0.05):
    out.mkdir(parents=True, exist_ok=False)
    initial = read(initial_dir/'final.traj')
    final = read(final_dir/'final.traj')
    # ASE trajectories carry actual converged endpoint SinglePointCalculator data.
    for endpoint in (initial, final):
        energy = float(endpoint.get_potential_energy())
        forces = endpoint.get_forces()
        endpoint.calc = SinglePointCalculator(endpoint, energy=energy, forces=forces)
    images = [initial] + [initial.copy() for _ in range(interior)] + [final]
    neb = NEB(images, climb=True, method='improvedtangent', parallel=False)
    neb.interpolate(method='idpp', apply_constraint=False)
    for i, image in enumerate(images[1:-1], start=1):
        image.calc = calculator(image, out/f'image-{i}-scf.txt')
    started = time.monotonic()
    opt = FIRE(neb, logfile=str(out/'neb.log'), dt=0.08, maxstep=0.12)
    records = []
    previous_wall = [0.0]
    def checkpoint():
        band_forces = neb.get_forces()
        record = {'step': opt.nsteps, 'wall_s': time.monotonic()-started,
                  'band_fmax_eV_A': float(np.linalg.norm(band_forces.reshape(-1, 3), axis=1).max()),
                  'band_forces_eV_A': band_forces.tolist(), 'images': []}
        for image in images:
            record['images'].append({'energy_extrapolated_eV': float(image.get_potential_energy()),
                                     'forces_eV_A': image.get_forces().tolist(),
                                     'geometry': geometry(image)})
        p = out/f'band-{opt.nsteps:03d}.json'
        save(p, record)
        write(out/f'band-{opt.nsteps:03d}.traj', images)
        records.append(manifest(p))
        print('NEB_FORCE_POINT', opt.nsteps, record['band_fmax_eV_A'],
              [x['energy_extrapolated_eV'] for x in record['images']], flush=True)
        work_wall = record['wall_s']-previous_wall[0]
        expected_wall = interior*(156.34 if opt.nsteps == 0 else 64.0)
        previous_wall[0] = record['wall_s']
        if work_wall > 3*expected_wall:
            save(out/'cost-obstruction.json', {
                'actual_band_work_wall_s': work_wall,
                'expected_from_real_pilot_wall_s': expected_wall,
                'completed_geometry_record': records[-1],
                'disposition': 'Stop at3x measured work cost; no larger local cap or automatic retry'})
            raise RuntimeError('NEB exceeded3x pilot-predicted band work cost')
    opt.attach(checkpoint)
    converged = opt.run(fmax=fmax, steps=steps)
    for i, image in enumerate(images[1:-1], start=1):
        image.calc.write(str(out/f'image-{i}.gpw'), mode='all')
    save(out/'result.json', {'converged': bool(converged),
         'method': 'ASE climbing image improvedtangent; IDPP initial path',
         'fmax_target_eV_A': fmax, 'interior_images': interior,
         'wall_s': time.monotonic()-started, 'records': records,
         'final': records[-1], 'coverage_ML': 1.0,
         'reaction': '2D(fcc)/two surface Pd -> D2 gas above clean matched slab'})
    if not converged:
        raise RuntimeError('Band force criterion unmet within original step/cost cap')
