"""Source-bound frozen-geometry contact with a functional operator identity.

Ordinary GPAW rejects coincident atoms. The exact co-located two-H operator
is compared algebraically and by sampled built radial functions, never
represented as a completed coincident-atom SCF.
"""
import argparse
from copy import deepcopy
import gc
import inspect
import json
from pathlib import Path
import shutil
import time

import numpy as np
from ase import Atoms
from ase.units import Hartree
import gpaw.ae
import gpaw.basis_data
import gpaw.setup
import gpaw.utilities
import ase.units

from contact_setup import ContactHeliumSetup
from contact_state import binding, check_stop, require, reuse_state, scf, write_json


METHOD = dict(xc='PBE', spinpol=False, occupation_width_eV=.1,
              geometry_H_setup='standard_PAW', BO_H_setup='ae',
              Li_setup='default_1valence', Hae_alpha1_Bohr_m2=10,
              energy_convention='internal=free-minus-negativeTS')


def load(path):
    return json.loads(Path(path).read_bytes())


def verify_binding(path, item):
    require(binding(path, Path(path).parent) | {'path': item['path']} == item,
            'Changed source-bound file: ' + item['path'])


def snapshot_epoch(input_path, destination, output):
    input_path = Path(input_path)
    source_root = input_path.parent
    epoch_path = source_root/'source-manifest.json'
    epoch = load(epoch_path)
    destination.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(input_path, destination/'input.json')
    shutil.copyfile(epoch_path, destination/'source-manifest.json')
    for item in epoch['files']:
        relative = Path(item['path'])
        require(not relative.is_absolute() and '..' not in relative.parts,
                'Unsafe source manifest path')
        source = source_root/'source'/relative
        verify_binding(source, item)
        target = destination/'source'/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    return dict(input=binding(destination/'input.json', output),
                source_manifest=binding(destination/'source-manifest.json', output))


def snapshot_baseline(base_dir, output):
    base_dir = Path(base_dir)
    run = load(base_dir/'run-result.json')
    require(run['schema'] == 'hydride-run-v1', 'Unknown original finite-R run')
    original_input = base_dir.parent/run['source_epoch']['input']['path']
    original_manifest = base_dir.parent/run['source_epoch']['source_manifest']['path']
    verify_binding(original_input, run['source_epoch']['input'])
    verify_binding(original_manifest, run['source_epoch']['source_manifest'])
    require(load(original_input) == run['specification'], 'Original specification changed')
    epoch = snapshot_epoch(original_input, output/'base', output)
    references = dict(epoch)
    for key, relative in [('run_result', 'run-result.json'),
                          ('quartet', run['quartet']['path']),
                          ('geometry', run['geometry_file']['path']),
                          ('occupied_host', 'bo/occupied-host.json'),
                          ('vacant_host', 'bo/vacant-host.json')]:
        source = base_dir/relative
        if key in ('quartet', 'geometry'):
            verify_binding(source, run['quartet' if key == 'quartet' else 'geometry_file'])
        target = output/'base'/'calculation'/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        references[key] = binding(target, output)
    quartet = load(base_dir/run['quartet']['path'])
    require(quartet['geometry'] == run['geometry'] and
            quartet['specification'] == run['specification'], 'Original quartet changed')
    require(load(base_dir/run['geometry_file']['path']) == run['geometry'],
            'Original frozen geometry file changed')
    for role in ('occupied-host', 'vacant-host'):
        raw = load(base_dir/'bo'/(role + '.json'))
        entry = quartet['states'][role]
        require(raw == {k: v for k, v in entry.items() if k != 'state_file'},
                'Original host record differs from accepted quartet')
        verify_binding(base_dir/'bo'/(role + '.json'), entry['state_file'])
    # Preserve every standalone positive-R probe JSON used by independent
    # baseline replay; GPWs are not necessary for its raw-energy provenance.
    for branch in ('q0', 'q1'):
        for role in ('occupied-probes', 'vacant-probes'):
            for state in quartet['states'][branch][role]:
                item = state['state_file']
                relative = Path(item['path'])
                require(relative.name == str(relative), 'Unsafe baseline probe binding')
                source = base_dir/'bo'/relative
                verify_binding(source, item)
                require(load(source) == {k: v for k, v in state.items() if k != 'state_file'},
                        'Original finite probe differs from accepted quartet')
                shutil.copyfile(source, output/'base'/'calculation'/'bo'/relative)
    volume = run['specification'].get('volume_source')
    if volume is not None:
        relative = Path(volume['path'])
        require(relative.name == str(relative), 'Unsafe original volume binding')
        source = base_dir.parent/relative
        verify_binding(source, volume)
        shutil.copyfile(source, output/'base'/relative)
    return run, references


def snapshot_operators(output):
    sources = {}
    from gpaw.old import density, hamiltonian
    modules = dict(ae=gpaw.ae, setup=gpaw.setup, basis_data=gpaw.basis_data,
                   hamiltonian=hamiltonian, density=density)
    for name, module in modules.items():
        target = output/'operator-source'/(name + '.py')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(inspect.getsourcefile(module), target)
        sources[name] = binding(target, output)
    # These exact bytes establish the public coincident-atom obstruction too.
    from gpaw.old import calculator as old_calculator
    for name, module in [('legacy_calculator', old_calculator), ('utilities', gpaw.utilities)]:
        target = output/'operator-source'/(name + '.py')
        shutil.copyfile(inspect.getsourcefile(module), target)
    shutil.copyfile(inspect.getsourcefile(ase.units), output/'operator-source'/'ase_units.py')
    return sources


def operator_identity(output, sources):
    hydrogen = gpaw.ae.HydrogenAllElectronSetup(alpha1=10., alpha2=300.)
    hydrogen.build(None)
    helium = ContactHeliumSetup()
    r = np.linspace(0., 2., 801)
    empty_names = ('Delta_pL', 'Delta_iiL', 'dO_ii', 'dC_ii', 'M_pp',
                   'K_p', 'M_p', 'MB_p', 'X_p')

    def sample(setup):
        ghat = setup.ghat_l[0].map(r)
        return dict(Z=setup.Z, Nv=setup.Nv, Nc=setup.Nc, ni=setup.ni,
                    E=float(setup.E), M=float(setup.M), Delta0=float(setup.Delta0),
                    alpha1=setup.alpha1, alpha2=setup.alpha2,
                    occupation=list(setup.f_j), ghat0_radial=ghat.tolist(),
                    compensation_radial=(setup.Delta0 * ghat).tolist(),
                    vbar_radial=setup.vbar.map(r).tolist(),
                    nct_radial=setup.nct.map(r).tolist(),
                    tauct_radial=setup.tauct.map(r).tolist(),
                    empty_shapes={name: list(getattr(setup, name).shape)
                                  for name in empty_names},
                    Kc=float(setup.Kc), MB=float(setup.MB), ExxC=float(setup.ExxC),
                    Nct=float(setup.Nct), xc_correction=setup.xc_correction,
                    pt_count=len(setup.pt_j), type=setup.type)

    radial = dict(schema='hydride-contact-radial-v1', r_Bohr=r.tolist(),
                  H=sample(hydrogen), He=sample(helium))
    require(np.allclose(radial['He']['compensation_radial'],
                        2 * np.asarray(radial['H']['compensation_radial']),
                        rtol=0, atol=1e-12) and
            np.allclose(radial['He']['vbar_radial'],
                        2 * np.asarray(radial['H']['vbar_radial']),
                        rtol=0, atol=1e-12), 'Built Z2 radial operator is not two Hae operators')
    path = output/'operator-source'/'radial-operators.json'
    write_json(path, radial)
    return dict(schema='hydride-contact-operator-identity-v1',
                radial_file=binding(path, output), source=sources,
                coincident_native=dict(performed=False,
                    reason='public coincident-atom guard; source-functional identity instead',
                    guard_sources={name: binding(output/'operator-source'/(name+'.py'), output)
                                   for name in ('legacy_calculator', 'utilities')}),
                algebra=dict(
                    domain='every common trial density/orbital and commissioned charge/electron branch on the same cell/grid',
                    equal_terms=['kinetic', 'xc', 'Hartree of identical density+compensation',
                                 'linear local vbar', 'zero core/projector/overlap corrections',
                                 'zero atomic Eref'],
                    compensation='2*Delta0_H*ghat_H=Delta0_He*ghat_He',
                    local_potential='2*vbar_H=vbar_He',
                    difference='E_twoH[n]-E_He[n]=2*M_H-M_He=-2*M_H',
                    expected_difference_Hartree=2*hydrogen.M-helium.M,
                    expected_difference_eV=(2*hydrogen.M-helium.M)*Hartree,
                    wrong_M_half_scaling_residual_eV=(helium.M-2*hydrogen.M)*Hartree,
                    wrong_vbar_unscaled_max_abs_Hartree=float(np.max(np.abs(hydrogen.vbar.map(r)))),
                    source_entry_points={
                        'M':'hamiltonian.py Hamiltonian.update_corrections: e_coulomb += setup.M',
                        'vbar':'hamiltonian.py RealSpaceHamiltonian.update_pseudo_potential: LFC sum and finegd.integrate(vbar_g, nt_g)',
                        'Delta0':'density.py CompensationChargeExpansionCoefficients.calculate: Q_L[0] += setup.Delta0',
                        'Eref':'setup.py Setups.__init__: self.Eref += n * setup.E'}))


def empty_fine_reference(original, spec, output):
    record = deepcopy(original)
    require(record['backend'] == 'analytic-empty-vacuum' and record['symbols'] == [] and
            record['electrons'] == 0 and all(x == 0 for x in record['energy'].values()),
            'Original gas vacuum is not the exact empty zero')
    record['specification'] = spec
    record['reference'] = 'Exact empty neutral vacuum at this grid: no Native SCF or GPW claimed.'
    path = output/'states'/'vacant-host.json'
    write_json(path, record)
    return record | {'state_file': binding(path, output)}


def prior_fine_host(prior_path, geometry, atoms, spec, output):
    prior_path = Path(prior_path)
    if prior_path.is_dir():
        prior_path /= 'contact-native.json'
    prior = load(prior_path)
    require(prior['schema'] == 'hydride-contact-native-v1' and prior['case'] == 'G' and
            prior['charge'] == 1 and prior['h_A'] == .16 and
            prior['geometry'] == geometry and prior['specification'] == spec,
            'Prior q1 is not the same Gfine source/geometry/method')
    entry = prior['states']['occupied-host']
    source = prior_path.parent/entry['state_file']['path']
    verify_binding(source, entry['state_file'])
    require(load(source) == {k: v for k, v in entry.items() if k != 'state_file'},
            'Prior q1 host record changed')
    record = reuse_state(source, atoms, 0, .16, output, 'occupied-host')
    if 'log' in record:
        log = prior_path.parent/record['log']['path']
        verify_binding(log, record['log'])
        destination = output/record['log']['path']
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(log, destination)
    target = output/'provenance'/'prior-q1-contact-native.json'
    shutil.copyfile(prior_path, target)
    return record, binding(target, output)


def run(input_path, base_dir, output, prior_q1=None):
    started = time.process_time()
    config = load(input_path)
    require(set(config) == {'schema', 'case', 'charge', 'h_A', 'base_role', 'absolute_stop_utc'} and
            config['schema'] == 'hydride-contact-input-v1' and config['case'] in ('F', 'G') and
            type(config['charge']) is int and config['charge'] in (0, 1) and
            config['h_A'] in (.2, .16), 'Uncommissioned contact input')
    require(config['absolute_stop_utc'] == '2026-11-20T18:00:00Z', 'Changed absolute stop')
    expected_role = 'F_h' if config['case'] == 'F' and config['h_A'] == .16 else config['case']
    require(config['base_role'] == expected_role, 'Wrong original reference role')
    check_stop(config['absolute_stop_utc'])
    print('CONTACT_NATIVE_START case=' + config['case'] +
          ' q=' + str(config['charge']) + ' h_A=' + str(config['h_A']), flush=True)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    epoch = snapshot_epoch(input_path, output/'provenance', output)
    run_record, references = snapshot_baseline(base_dir, output)
    # Bind all three actually imported producer modules to the staged epoch.
    import contact_setup
    import contact_state
    executable_manifest = load(output/'provenance'/'source-manifest.json')
    for module_path in (Path(__file__), Path(contact_setup.__file__), Path(contact_state.__file__)):
        matches = [item for item in executable_manifest['files']
                   if Path(item['path']).name == module_path.name]
        require(len(matches) == 1, 'Producer module absent/ambiguous in source epoch')
        verify_binding(module_path, matches[0])
    operators = snapshot_operators(output)
    proof = operator_identity(output, operators)
    original_spec = run_record['specification']
    require(original_spec['case'] == config['case'] and original_spec['method'] == METHOD,
            'Original method/case not commissioned')
    require(original_spec['h_A'] == (config['h_A'] if config['case'] == 'F' else .2),
            'Wrong original grid reference')
    spec = deepcopy(original_spec)
    spec.pop('absolute_stop_utc')
    spec['h_A'] = config['h_A']
    spec['method']['contact_He_setup'] = 'HaeZ2'
    geometry = run_record['geometry']
    target_index = geometry['target_index']
    atoms = Atoms(geometry['symbols'], positions=geometry['positions_A'],
                  masses=geometry['masses_u'], cell=geometry['cell_A'], pbc=True)
    require(geometry['passed'] and atoms[target_index].symbol == 'H' and
            np.allclose(atoms.positions[target_index], geometry['path']['target_A'],
                        rtol=0, atol=1e-12) and
            np.allclose(geometry['path']['direction'], np.ones(3)/np.sqrt(3),
                        rtol=0, atol=1e-12), 'Wrong original target or [111] ray')
    vacant = atoms.copy()
    del vacant[target_index]
    helium_atoms = atoms.copy()
    helium_atoms[target_index].symbol = 'He'
    helium_atoms[target_index].mass = 2 * atoms[target_index].mass
    states = {}
    fine_g = config['case'] == 'G' and config['h_A'] == .16
    if fine_g:
        if config['charge'] == 0:
            require(prior_q1 is not None, 'Gfine q0 requires actual completed q1 fine host')
            states['occupied-host'], references['prior_q1'] = prior_fine_host(
                prior_q1, geometry, atoms, spec, output)
        else:
            states['occupied-host'] = scf(atoms.copy(), spec, 0, output, 'occupied-host',
                                          operators, config['absolute_stop_utc'])
        states['vacant-host'] = empty_fine_reference(
            load(output/references['vacant_host']['path']), spec, output)
    else:
        for role, frame in [('occupied-host', atoms), ('vacant-host', vacant)]:
            key = role.replace('-', '_')
            states[role] = reuse_state(output/references[key]['path'], frame, 0,
                                       config['h_A'], output, role)
    states['occupied-probe-He'] = scf(helium_atoms, spec, config['charge'], output,
                                      'occupied-probe-He', operators, config['absolute_stop_utc'])
    gc.collect()
    if config['charge'] == 0:
        states['vacant-probe-H'] = reuse_state(
            output/states['occupied-host']['state_file']['path'], atoms, 0,
            config['h_A'], output, 'vacant-probe-H')
    else:
        states['vacant-probe-H'] = scf(atoms.copy(), spec, 1, output, 'vacant-probe-H',
                                       operators, config['absolute_stop_utc'])
    host_ne = states['occupied-host']['electrons']
    require(abs(states['occupied-probe-He']['electrons'] - host_ne - 1 + config['charge']) < 1e-7 and
            abs(states['vacant-probe-H']['electrons'] - host_ne + config['charge']) < 1e-7,
            'Probe legs do not preserve commissioned charge/electron branch')
    setup = ContactHeliumSetup()
    result = dict(schema='hydride-contact-native-v1', case=config['case'], charge=config['charge'],
                  h_A=config['h_A'], base_role=config['base_role'], geometry=geometry,
                  specification=spec, states=states, source_epoch=epoch,
                  source_references=references, operator_sources=operators,
                  absolute_convention=dict(alpha1=10., alpha2=300., H_M_Hartree=setup.M/4,
                                           He_M_Hartree=setup.M, H_E_reference_Hartree=0.,
                                           He_E_reference_Hartree=0., ASE_Hartree_eV=Hartree,
                                           ASE_CODATA_version=ase.units.__codata_version__,
                                           ASE_units_source=binding(output/'operator-source'/'ase_units.py', output),
                                           B_native_eV=-setup.M/2*Hartree,
                                           representation='same-alpha all-electron Z2',
                                           arbitrary_atomic_offset_eV=0.,
                                           energy='internal=free-minus-negativeTS',
                                           formula='W0=Q_He+B_native-B_source+regular_images'),
                  operator_identity=proof,
                  recovery={key: (
                      'raw converged JSON/geometry/energy and bound log only; no GPW or wavefunction restart claimed'
                      if 'log' in value else
                      'original JSON/geometry/energy retained byte-exact, or declared exact empty zero; old GPW descriptors are provenance only, no local GPW/log claimed')
                            for key, value in states.items()},
                  process_cpu_s=time.process_time()-started,
                  scope='Signed contact candidate only; independent operator replay and h/cell qualification required; no derivative or isolated-limit claim')
    write_json(output/'contact-native.json', result)
    print('CONTACT_NATIVE_COMPLETE case=' + config['case'] +
          ' q=' + str(config['charge']) + ' h_A=' + str(config['h_A']), flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--base-dir', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--prior-q1')
    args = parser.parse_args()
    run(args.input, args.base_dir, args.output, args.prior_q1)
