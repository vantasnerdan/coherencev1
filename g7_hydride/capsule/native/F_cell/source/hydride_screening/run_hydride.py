"""Finite pilot or complete two-charge hydride/gas case from a frozen epoch."""
import argparse
import gc
import json
from pathlib import Path
import time

from ase import Atom

from geometry import build, DIRECTION, MASS_D
from native_state import manifest, require, scf, write_json
from quartet import quartet, probe_stem
from relaxation import load_geometry, relax

METHOD = dict(xc='PBE', spinpol=False, occupation_width_eV=.1,
              geometry_H_setup='standard_PAW', BO_H_setup='ae',
              Li_setup='default_1valence', Hae_alpha1_Bohr_m2=10,
              energy_convention='internal=free-minus-negativeTS')


def validate_epoch(input_path):
    input_path = Path(input_path)
    specification = json.loads(input_path.read_bytes())
    require(specification['schema'] == 'hydride-specification-v1' and
            specification['method'] == METHOD, 'Changed commissioned Native family')
    epoch_path = input_path.parent/'source-manifest.json'
    epoch = json.loads(epoch_path.read_bytes())
    for item in epoch['files']:
        require((manifest(input_path.parent/'source'/item['path']) | {'path': item['path']}) == item,
                'Changed frozen executable source')
    if specification['volume_source'] is not None:
        require(manifest(input_path.parent/specification['volume_source']['path']) == specification['volume_source'],
                'Changed original relaxed PW volume source')
    return specification, dict(input=manifest(input_path), source_manifest=manifest(epoch_path))


def pilot(spec, gas_spec, output, epoch):
    atoms, geometry = relax(spec, output/'geometry')
    host = scf(atoms, spec, 0, output/'bo', 'occupied-host')
    atoms.calc = None
    del atoms
    gc.collect()
    gas = build(gas_spec)
    target = gas.positions[0].copy()
    del gas[0]
    gas += Atom('H', position=target-.1*DIRECTION, mass=MASS_D)
    stem = probe_stem('vacant', 1, .1)
    probe = scf(gas, gas_spec, 1, output/'gas-probe', stem,
                extra=dict(R_A=.1, projectile_A=gas.positions[-1].tolist()))
    require(probe['electrons'] == 0 and probe['backend'] == 'gpaw-native',
            'Zero-electron gas probe must actually complete in Native')
    result = dict(schema='hydride-feasibility-pilot-v1', specification=spec,
                  geometry=geometry, F_occupied_host=host, G_zero_electron_probe=probe,
                  source_epoch=epoch, scope='Feasibility/cost only; no complete lowering or control qualification.')
    return result


def run(input_path, output, radii, geometry_source=None, reuse=None, stage='production', gas_input=None):
    spec, epoch = validate_epoch(input_path)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    print('HYDRIDE_RUN_START case='+spec['case']+' stage='+stage, flush=True)
    started = time.process_time()
    if stage == 'pilot':
        require(spec['case'] == 'F' and gas_input is not None, 'Pilot needs F and declared G input')
        gas_spec, gas_epoch = validate_epoch(gas_input)
        require(gas_spec['case'] == 'G', 'Pilot gas input must be atomic G')
        result = pilot(spec, gas_spec, output, epoch)
        result['gas_source_epoch'] = gas_epoch
        result['process_cpu_s'] = time.process_time()-started
        write_json(output/'pilot.json', result)
        print('HYDRIDE_PILOT_COMPLETE cpu_s='+str(result['process_cpu_s']), flush=True)
        return result
    if geometry_source:
        atoms, geometry = load_geometry(spec, geometry_source, output/'geometry')
    else:
        atoms, geometry = relax(spec, output/'geometry')
    result = quartet(atoms, geometry, spec, radii, output/'bo', reuse=reuse)
    record = dict(schema='hydride-run-v1', specification=spec, geometry=geometry,
                  geometry_file=manifest(output/'geometry/relaxation.json') | {'path': 'geometry/relaxation.json'},
                  quartet=manifest(output/'bo/quartet.json') | {'path': 'bo/quartet.json'},
                  source_epoch=epoch, process_cpu_s=time.process_time()-started)
    write_json(output/'run-result.json', record)
    print('HYDRIDE_RUN_COMPLETE case='+spec['case']+' cpu_s='+str(record['process_cpu_s']), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', choices=['pilot', 'production'], default='production')
    parser.add_argument('--gas-input', type=Path)
    parser.add_argument('--geometry-frame', type=Path)
    parser.add_argument('--reuse-state', action='append', default=[])
    parser.add_argument('--radii', type=float, nargs='+', default=[1.6, 1., .5, .3, .2, .1])
    args = parser.parse_args()
    reuse = dict(value.split('=', 1) for value in args.reuse_state)
    run(args.input, args.output, args.radii, args.geometry_frame, reuse, args.stage, args.gas_input)
