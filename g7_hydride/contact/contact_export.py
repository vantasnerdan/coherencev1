"""Export separate signed R=0 regular-part sidecars without changing positive-R curves.

Export requires an actual complete q1 F/G coarse/fine quartet first. q0 is
optional only as a complete second F/G coarse/fine group. Replay uses NumPy and
stdlib, original saved states and frozen sources; it never reruns Native.
"""
import argparse
import json
from pathlib import Path
import shutil

from contact_replay import _log_file, load_contact, qualify_pair
from log_privacy import sanitize_log
from export_pairs import check_pair
from hydride_replay import binding, bound_file, close, compare_json, read_json, require

SCHEMA = 'hydride-signed-contact-sidecar-v1'
MANIFEST_SCHEMA = 'hydride-contact-sidecars-v1'
SLOPE_REASON = ('Unqualified internal-energy-force derivative: internal=free-minus-negativeTS; '
                'no matched derivative/entropy-force calculation was commissioned or qualified. '
                'No zero-slope or symmetry assumption is made.')
SCOPE = ('Signed regular contact W0=lim_R->0(Vtotal-C/R), not the singular Vtotal(0), '
         'not a positive screening magnitude. Negative W0 denotes lowering. Same frozen '
         'finite periodic cell, target, spectators and [111] ray as the SHA-bound finite-R curve. '
         'The Gaussian nuclear zero-G mapping does not remove periodic electronic or multipole '
         'images. Contact h sensitivity is qualified separately; existing finite-R F cell '
         'qualification is not a contact-cell convergence test. No exact isolated limit is claimed.')


def _sidecar(coarse, fine, curve, curve_binding, coarse_path, fine_path, gas_identity):
    case, q = coarse['case'], coarse['charge']
    check_pair(curve)
    require(coarse['h_A'] == .2 and fine['h_A'] == .16, 'Missing registered contact h pair')
    cfg = curve['configuration']
    require(cfg['case'] == case and type(cfg['probe_charge']) is int and cfg['probe_charge'] == q,
            'Wrong curve host/charge branch')
    require(cfg['probe_added_electrons'] == 1-q, 'Wrong finite-curve electron branch')
    for role, field in (('occupied-host', 'occupied_host_electrons'), ('vacant-host', 'vacant_host_electrons')):
        close(cfg[field], coarse['electrons'][role], 'Curve/contact host Ne')
    for role, field in (('occupied-probe-He', 'occupied_probe_electrons'), ('vacant-probe-H', 'vacant_probe_electrons')):
        close(cfg[field], [coarse['electrons'][role]] * len(curve['R_A']), 'Curve/contact probe Ne')
    require(curve['geometry_provenance'] == coarse['geometry'], 'Curve/contact frozen geometry or ray differs')
    compare_json(curve['source_bindings']['native'], coarse['finite_base_sources'],
                 'Curve/contact original Native source identity')
    close(curve['R_A'], [r['R_A'] for r in coarse['finite_base_rows']], 'Unchanged positive-R support')
    close(curve['V_pair_total_eV'], [r['V_pair_total_eV'] for r in coarse['finite_base_rows']],
          'Unchanged positive-R energies')
    require(curve['qualification']['qualified'] is True, 'Unqualified finite-R curve')
    require(gas_identity is not None, 'Missing same-charge source-bound gas absolute-convention proof')
    require(gas_identity['charge'] == q, 'Wrong absolute-convention proof charge')
    require(gas_identity['sources'] == {k: dict(bytes=v['bytes'], sha256=v['sha256'])
                                      for k, v in coarse['operator_sources'].items()},
            'F/G absolute identity belongs to different loaded operators')
    h = qualify_pair(coarse, fine)
    return dict(schema=SCHEMA, case=case, charge=q, electronic_branch='neutral-added-electron' if q == 0 else 'bare-no-added-electron',
                observable='signed_regular_contact_limit_of_target_conditioned_BO_total',
                R_A=0., signedW0_eV=coarse['signedW0_eV'],
                definition='lim_R->0(V_pair_total_eV-C_eV_A/R_A)',
                Q_He_eV=coarse['Q_He_eV'], raw_quartet_energies_eV=coarse['raw_quartet_energies_eV'],
                electrons=coarse['electrons'], geometry=coarse['geometry'], ray=coarse['geometry']['path'],
                specification=coarse['specification'], absolute_convention=coarse['absolute_convention'],
                periodic_mapping=coarse['periodic_mapping'], B_native_eV=coarse['B_native_eV'],
                B_native_minus_B_source_eV=coarse['unit_constant_residual_eV'],
                mapping_proof=coarse['mapping_proof'], gas_absolute_identity=gas_identity,
                uncertainty=dict(h_sensitivity=h, contact_cell_convergence='not measured',
                                 finite_R_cell_qualification=curve['qualification']['both_F_charge_h_and_cell_controls'],
                                 isolation_claim=False),
                derivative_eV_per_A=None, derivative_status=SLOPE_REASON,
                finite_curve=curve_binding, positive_R_arrays_untouched=True,
                native_sources=dict(coarse=coarse_path, fine=fine_path),
                raw_replay=dict(coarse=coarse, fine=fine), scope=SCOPE)


def _copy_tree(source, destination, archive_root, log_inventory):
    source, destination = Path(source), Path(destination)
    require(source.is_dir() and not source.is_symlink(), 'Missing/unsafe contact source directory')
    destination.mkdir(parents=True)
    for path in sorted(source.rglob('*')):
        require(not path.is_symlink(), 'Symlink not allowed in contact source tree')
        relative = path.relative_to(source)
        target = destination / relative
        if path.is_dir():
            target.mkdir(exist_ok=True)
        else:
            require(path.is_file(), 'Nonregular contact source member')
            target.parent.mkdir(parents=True, exist_ok=True)
            original = binding(path)
            if path.suffix.lower() in ('.txt', '.log'):
                raw = path.read_bytes()
                exported, changes = sanitize_log(raw)
                with target.open('xb') as stream:
                    stream.write(exported)
                after = binding(target)
                name = target.relative_to(archive_root).as_posix()
                log_inventory[name] = dict(ORIGINAL={k: original[k] for k in ('bytes', 'sha256')},
                                           EXPORTED={k: after[k] for k in ('bytes', 'sha256')},
                                           transformation=changes)
                require(target.read_bytes() == exported and binding(path) == original,
                        'Approved log identity/path export integrity failed')
            else:
                shutil.copyfile(path, target)
                require(binding(path) == original and binding(target) == original, 'Source bytes changed during copy')


def _json_new(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def export_sidecars(groups, curves, output):
    """groups maps integer q to {'F': (coarse,fine), 'G': (coarse,fine)}.

    The output is a new self-contained raw/source subtree plus untouched curve
    copies, sidecars and a binding manifest. A partial failed export has no
    complete manifest and cannot pass replay_sidecars.
    """
    require(set(groups) in ({1}, {0, 1}) and all(type(q) is int for q in groups), 'Complete q1 pair must precede optional q0')
    for group in groups.values():
        require(set(group) == {'F', 'G'} and all(len(v) == 2 for v in group.values()), 'Require both hosts and both h grids')
    output, curves = Path(output), Path(curves)
    require(not output.exists(), 'Contact export output must be new')
    require(all(not output.resolve().is_relative_to(Path(p).resolve()) for g in groups.values()
                for pair in g.values() for p in pair), 'Output cannot be inside an input tree')
    # Complete physical/source replay before any output is qualified.
    loaded = {}
    curve_records = {}
    for q in (1, 0):
        if q not in groups:
            continue
        loaded[q] = {}
        for case in ('F', 'G'):
            coarse, fine = (load_contact(p) for p in groups[q][case])
            require(coarse['case'] == fine['case'] == case and coarse['charge'] == fine['charge'] == q,
                    'Contact CLI input host/charge mismatch')
            qualify_pair(coarse, fine)
            loaded[q][case] = (coarse, fine)
            name = case + '_q' + str(q) + '_BO_pair.json'
            path = curves / name
            curve_records[case, q] = (path, read_json(path), binding(path, 'curves/' + name))
        gas_identity = loaded[q]['G'][0]['gas_absolute_identity']
        require(loaded[q]['F'][0]['source_manifest'] == loaded[q]['G'][0]['source_manifest'],
                'F/G contact executable epochs differ')
        if q == 0:
            require(loaded[q]['G'][1]['source_references']['prior_q1']['sha256'] ==
                    loaded[1]['G'][1]['native_binding']['sha256'], 'q0 fine gas did not reuse the archived actual q1 host')
        for case in ('F', 'G'):
            coarse, fine = loaded[q][case]
            _, curve, cb = curve_records[case, q]
            stem = case + '_q' + str(q)
            _sidecar(coarse, fine, curve, cb, 'sources/' + stem + '/coarse',
                     'sources/' + stem + '/fine', gas_identity)
    output.mkdir(parents=True)
    (output / 'curves').mkdir()
    sidecars = []
    log_inventory = {}
    for q in (1, 0):
        if q not in groups:
            continue
        gas_identity = loaded[q]['G'][0]['gas_absolute_identity']
        for case in ('F', 'G'):
            stem = case + '_q' + str(q)
            coarse_path, fine_path = 'sources/' + stem + '/coarse', 'sources/' + stem + '/fine'
            for source, relative in zip(groups[q][case], (coarse_path, fine_path)):
                _copy_tree(source, output / relative, output, log_inventory)
            cpath, curve, cb = curve_records[case, q]
            require(binding(cpath, cb['path']) == cb, 'Finite curve changed during export')
            shutil.copyfile(cpath, output / cb['path'])
            require(binding(output / cb['path'], cb['path']) == cb, 'Copied finite curve bytes differ')
            record = _sidecar(*loaded[q][case], curve, cb, coarse_path, fine_path, gas_identity)
            name = stem + '_contact.json'
            _json_new(output / name, record)
            sidecars.append(binding(output / name))
    sources = [binding(path, path.relative_to(output).as_posix()) for path in sorted((output / 'sources').rglob('*')) if path.is_file()]
    manifest = dict(schema=MANIFEST_SCHEMA, charges=[q for q in (1, 0) if q in groups],
                    sidecars=sidecars, sources=sources, log_inventory=log_inventory,
                    curves=[curve_records[case, q][2] for q in (1, 0) if q in groups for case in ('F', 'G')])
    _json_new(output / 'contact-sidecars.json', manifest)
    return manifest


def replay_sidecars(directory):
    """Recompute all sidecars from archived raw/source bytes; no Native/SciPy/GPAW.

    Every structure and discrete flag is exact; floating values use atol=1e-9,
    rtol=0. Positive-R curve bytes and original source manifests are mandatory.
    """
    directory = Path(directory)
    mpath = directory / 'contact-sidecars.json'
    original = binding(mpath)
    manifest = read_json(mpath)
    require(set(manifest) == {'schema', 'charges', 'sidecars', 'sources', 'curves', 'log_inventory'} and
            manifest['schema'] == MANIFEST_SCHEMA and manifest['charges'] in ([1], [1, 0]) and
            all(type(q) is int for q in manifest['charges']), 'Wrong contact sidecar manifest')
    require(isinstance(manifest['sources'], list) and manifest['sources'], 'Missing raw/source archive')
    for key in ('sources', 'curves', 'sidecars'):
        names = [item['path'] for item in manifest[key]]
        require(len(names) == len(set(names)), 'Duplicate ' + key + ' archive member')
        for item in manifest[key]:
            bound_file(directory, item, 'Archived ' + key)
    actual_sources = {path.relative_to(directory).as_posix() for path in (directory / 'sources').rglob('*') if path.is_file()}
    require(actual_sources == {item['path'] for item in manifest['sources']}, 'Unbound/missing raw archive source')
    inventory = manifest['log_inventory']
    require(type(inventory) is dict, 'Missing ORIGINAL/EXPORTED log accounting')
    expected_logs = {name for name in actual_sources if Path(name).suffix.lower() in ('.txt', '.log')}
    require(set(inventory) == expected_logs, 'Log ORIGINAL/EXPORTED accounting differs from archived logs')
    for name, entry in inventory.items():
        _log_file(directory, dict(path=name, **entry['ORIGINAL']), inventory, directory, 'Portable log')
    expected_sidecars = [case + '_q' + str(q) + '_contact.json' for q in manifest['charges'] for case in ('F', 'G')]
    expected_curves = ['curves/' + case + '_q' + str(q) + '_BO_pair.json' for q in manifest['charges'] for case in ('F', 'G')]
    require([i['path'] for i in manifest['sidecars']] == expected_sidecars and
            [i['path'] for i in manifest['curves']] == expected_curves, 'Wrong complete charge/host archive support')
    results = {}
    q1_fine_g = None
    for q in manifest['charges']:
        loaded = {}
        for case in ('F', 'G'):
            stem = case + '_q' + str(q)
            loaded[case] = tuple(load_contact(directory / 'sources' / stem / level,
                                              log_inventory=inventory, archive_root=directory)
                                 for level in ('coarse', 'fine'))
            require(all(v['case'] == case and v['charge'] == q for v in loaded[case]), 'Archived source branch mismatch')
        gas_identity = loaded['G'][0]['gas_absolute_identity']
        require(loaded['F'][0]['source_manifest'] == loaded['G'][0]['source_manifest'],
                'Archived F/G executable epochs differ')
        if q == 1:
            q1_fine_g = loaded['G'][1]
        else:
            require(loaded['G'][1]['source_references']['prior_q1']['sha256'] ==
                    q1_fine_g['native_binding']['sha256'], 'Archived q0 fine host is not the actual q1 source')
        for case in ('F', 'G'):
            stem = case + '_q' + str(q)
            cb = next(v for v in manifest['curves'] if v['path'] == 'curves/' + stem + '_BO_pair.json')
            curve = read_json(bound_file(directory, cb, 'Unchanged finite curve'))
            expected = _sidecar(*loaded[case], curve, cb, 'sources/' + stem + '/coarse',
                                'sources/' + stem + '/fine', gas_identity)
            saved = read_json(directory / (stem + '_contact.json'))
            compare_json(saved, expected, 'Complete source-bound ' + stem + ' contact sidecar')
            results[stem] = dict(signedW0_eV=expected['signedW0_eV'],
                                 h_sensitivity=expected['uncertainty']['h_sensitivity'])
    require(binding(mpath) == original, 'Contact manifest changed during replay')
    return dict(schema='hydride-contact-sidecar-replay-v1', qualified=True,
                charges=manifest['charges'], contacts=results, manifest=original,
                scope=SCOPE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', type=Path, help='Replay an existing contact archive instead of exporting')
    for q in (1, 0):
        for name in ('f', 'g', 'f-fine', 'g-fine'):
            parser.add_argument('--' + name + '-q' + str(q), type=Path)
    parser.add_argument('--curves', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.replay is not None:
        require(args.curves is None and args.output is None and
                all(getattr(args, name.replace('-', '_') + '_q' + str(q)) is None
                    for q in (1, 0) for name in ('f', 'g', 'f-fine', 'g-fine')),
                '--replay cannot be combined with export inputs')
        result = replay_sidecars(args.replay)
    else:
        require(args.curves is not None and args.output is not None, 'Export requires --curves and --output')
        groups = {}
        for q in (1, 0):
            values = [getattr(args, name + '_q' + str(q)) for name in ('f', 'g', 'f_fine', 'g_fine')]
            if q == 1 or any(v is not None for v in values):
                require(all(v is not None for v in values), 'Each claimed q requires F/G and both fine controls')
                groups[q] = dict(F=(values[0], values[2]), G=(values[1], values[3]))
        result = export_sidecars(groups, args.curves, args.output)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
