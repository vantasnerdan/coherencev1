"""Small, path-redacted Pd delivery after private saved-output replay; no Native."""
import argparse
import json
import math
from pathlib import Path
import re
import shutil

import reader_pd
from reader import (HERE, SUPPORT, BOHR_A, COULOMB_EV_A, PeriodicNuclearGreen,
                    binding, bound_file, close, compare_json, read_json, require)
from verify_contact_mapping import regular_mapping

# Absolute path tokens, including ones embedded in prose, are not delivery data.
_ABSOLUTE = re.compile(r"(?<![\w:/])/(?!/)[^\s\"'<>;,]+")
_WINDOWS = re.compile(r"[A-Za-z]:[\\/][^\s\"'<>;,]+")


def portable(value):
    """Retain provider basenames/hashes, never absolute machine paths."""
    if isinstance(value, dict):
        return {portable(k): portable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [portable(v) for v in value]
    if isinstance(value, str):
        value = _WINDOWS.sub(lambda m: re.split(r'[\\/]', m.group())[-1], value)
        return _ABSOLUTE.sub(lambda m: Path(m.group()).name, value)
    return value


def _new_json(path, value):
    with path.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')



def _common_support(cases):
    radii = {role: sorted(row['R_A'] for row in run['rows']) for role, run in cases.items()}
    common = sorted(set(radii.get('plain', [])) & set(radii.get('oxygen', [])))
    values = {role: {row['R_A']: row for row in run['rows']} for role, run in cases.items()}
    rows = [dict(R_A=radius,
                 plain_V_pair_total_eV=values['plain'][radius]['V_pair_total_eV'],
                 oxygen_V_pair_total_eV=values['oxygen'][radius]['V_pair_total_eV'],
                 oxygen_minus_plain_lowering_eV=values['oxygen'][radius]['lowering_eV'] -
                                                 values['plain'][radius]['lowering_eV'])
            for radius in common]
    return dict(schema='micro-pd-actual-common-support-v1', actual_radii_A=radii,
                missing_radii_A={role: [r for r in reader_pd.RADII if r not in actual]
                                 for role, actual in radii.items()}, common_radii_A=common, rows=rows,
                endpoint_R_A=.1 if .1 in common else None, interpolation=False,
                geometry_error_bound_eV=None, full_scientific_qualification=False)

def _curve(run, comparisons):
    rows = sorted(run['rows'], key=lambda row: row['R_A'])
    return dict(schema='coherencev1-target-conditioned-pair-v1', material='Pd',
                observable='target_conditioned_BO_total_interaction', **reader_pd._qualification(run['geometry']),
                paired_hk_qualified=comparisons['paired_hk_qualified'],
                geometry_force_qualified=comparisons['geometry_force_qualified'],
                full_scientific_qualification=False,
                ray_policy=reader_pd._ray_status(run['geometry']),
                R_A=[r['R_A'] for r in rows], V_pair_total_eV=[r['V_pair_total_eV'] for r in rows],
                completed_radii_A=[r['R_A'] for r in rows],
                missing_radii_A=[r for r in reader_pd.RADII if r not in {row['R_A'] for row in rows}],
                reference=dict(bare='isolated_unit_charge_coulomb',
                    energy_zero='Eop-Evp-Eho+Ehv; internal=free-minus-negativeTS; no fitted offset or far-point subtraction',
                    periodic_to_isolated='Accepted same-zero-G point-minus-Hae-Gaussian nuclear block only; electronic periodic images remain'),
                configuration=dict(specification=run['specification'], geometry=run['geometry'],
                    isotope='D', electronic_symbol='H', probe_charge=1, probe_added_electrons=0,
                    core_scope='Pd10/O6 frozen-core PAW; projector-free Hae/H1; no core-response correction'),
                source_bindings=dict(original_case=run['case_binding'], source_epoch=run['source_epoch']),
                qualification=dict(qualified=comparisons['qualified'], controls=comparisons,
                    paired_hk_qualified=comparisons['paired_hk_qualified'],
                    geometry_force_qualified=comparisons['geometry_force_qualified'],
                    geometry_error_bound_eV=None, full_scientific_qualification=False,
                    full_curve_higher_k_qualified=False, paired_cell_qualified=False, isolation_claim=False),
                scope=reader_pd.SCOPE)


def _common_curve(run, comparisons, radii):
    result = _curve(dict(run, rows=[row for row in run['rows'] if row['R_A'] in radii]), comparisons)
    result.update(curve_support='actual measured bare/coated common support only; no interpolation',
                  actual_completed_radii_A=sorted(row['R_A'] for row in run['rows']))
    return result


def _contact(run, curve_binding, comparisons):
    return dict(schema='hydride-signed-contact-sidecar-v1', material='Pd', case=run['case'],
                **reader_pd._qualification(run['geometry']),
                ray_policy=reader_pd._ray_status(run['geometry']),
                charge=1, electronic_branch='bare-no-added-electron', R_A=0.,
                observable='signed_regular_contact_limit_of_target_conditioned_BO_total',
                signedW0_eV=run['signedW0_eV'], Q_He_eV=run['Q_He_eV'],
                definition='lim_R->0(V_pair_total_eV-C_eV_A/R_A)',
                absolute_convention=run['absolute_convention'], mapping_proof=run['mapping_proof'],
                unit_source_proof=run['unit_source_proof'], functional_identity_replay=run['functional_identity_replay'],
                source_bindings=dict(original_case=run['case_binding'], source_epoch=run['source_epoch']),
                finite_curve=curve_binding, numerical_controls=comparisons,
                paired_hk_qualified=comparisons['paired_hk_qualified'],
                geometry_force_qualified=comparisons['geometry_force_qualified'],
                full_scientific_qualification=False,
                derivative_eV_per_A=None,
                derivative_status='Not qualified: no matched internal-energy/entropy-force derivative calculation; no zero-slope assumption',
                scope='Signed regular W0, not singular Vtotal(0), positive screening magnitude or Ueq. ' + reader_pd.SCOPE)


def readme(comparisons, geometries=None, cases=None):
    geometries = {} if geometries is None else geometries
    cases = {} if cases is None else cases
    support = _common_support(cases)
    coverage = ['| BO role | Actual completed R (A) | Missing commissioned R (A) |',
                '| --- | --- | --- |']
    for role, radii in sorted(support['actual_radii_A'].items()):
        coverage.append(f"| {role} | {radii} | {support['missing_radii_A'][role]} |")
    hk = 'PASS' if comparisons['paired_hk_qualified'] else 'NOT QUALIFIED'
    force = 'PASS' if comparisons['geometry_force_qualified'] else 'NOT FORCE-QUALIFIED'
    missing = ', '.join(comparisons['missing_controls']) or 'none'
    stops = ', '.join(sorted(set(comparisons['recorded_deadlines_utc'].values()))) or 'no BO input'
    statuses = dict(comparisons['geometry_by_role'])
    statuses.update({'geometry-only:' + role: reader_pd._qualification(run['geometry'])
                     for role, run in geometries.items()})
    lines = ['| Geometry | Status | Actual free fmax (eV/A) | Actual raw fmax (eV/A) |',
             '| --- | --- | ---: | ---: |']
    for role in dict.fromkeys(('plain', 'oxygen', *sorted(statuses))):
        q = statuses.get(role)
        lines.append(f'| {role} | missing | not available | not available |' if q is None else
                     f"| {role} | {q['relaxation_status']} | {q['free_fmax_eV_A']:.17g} | {q['raw_fmax_eV_A']:.17g} |")
    partial = any(not q['force_qualified'] for q in statuses.values())
    warning = ('**PARTLY RELAXED GEOMETRY: authority-labelled/conditional calibration only, NOT force-qualified.**\n\n'
               if partial else '')
    rays = ['| BO role | Ray status | Final coated admitted |',
            '| --- | --- | --- |']
    for role, ray in sorted(comparisons['ray_by_role'].items()):
        rays.append(f"| {role} | {ray['status']} | {str(ray['final_coated_admitted']).lower()} |")
    return (f'# Compact Pd/PdO saved-output delivery\n\n{warning}'
            '**GEOMETRY-CONVERGENCE ERROR IS UNBOUNDED: geometry_error_bound_eV=null. '
            'Residual forces and h/k tests do not bound pair-lowering/contrast error.**\n\n'
            f'Paired h/k endpoint criterion: {hk}. Geometry force criterion: {force}. Missing case files: {missing}.\n\n'
            + '\n'.join(lines) + '\n\n'
            + '\n'.join(rays) + '\n\n'
            + '\n'.join(coverage) + '\n\n'
            + f"Actual common positive-R support: {support['common_radii_A']} A. Only these exact measured radii are contrasted; no interpolation or invented far points. R0.1 is selected by value, never by execution ordinal.\n\n"
            'The `*_BO_pair_common.json` files retain only that measured overlap; the ordinary per-variant curves retain every actual completed point. `common-support.json` records the paired contrast and endpoint support. Signed contact remains a separate sidecar, not a fabricated positive-R point.\n\n'
            f'Recorded BO stop(s): {stops}; not an authorization assertion. '
            '`--allow-partial` permits missing case files only, never missing geometry/cap/site evidence or false force labels.\n\n'
            'Bare 15Pd+D versus coated 19Pd+4O+D in the same fixed cell; the added sheet is 4Pd+4O, not oxygen insertion at fixed total Pd. '
            'Only bottom substrate indices 10–14 are fixed. D and every other atom are free. '
            'The actual relaxed target position, migration and raw fixed-atom forces are reported, not an inferred retained-octahedral site.\n\n'
            'Entry_2980/2904 retains strict <=0.30 eV/A partial eligibility. The distinct entry_2994/3001 route permits the specifically authorized ORIGINAL step5 after the actual startup-failed resume: **partly relaxed, fmax 0.42 eV/A, no geometry error bound**. '
            'Its original Native state, specification, accepted step and process clock remain unchanged. The actual failed attempt is separate closed evidence and cost: zero new resumed SCFs/steps, no qualified checkpoint import or geometry advance. A startup failure is not misrepresented as a new cap. '
            'Any second-frame spread is sensitivity only, never a rigorous geometry-error bound.\n\n'
            'Commissioned six base radii: 1.6, 1.2, 0.8, 0.4, 0.2, 0.1 A; only the actual completed sets above are delivered. Near branches execute 0.1, 0.2, 0.4 first; optional far branches start independent cold sequences and reuse original fixed host/contact records, not recomputed hosts. '
            'Parallel container metadata is a checked derived view, not an executed monolithic specification; each Raw state retains its actual branch input, original bytes and execution/sequence indices. '
            'Paired Gamma versus 2x2x1 and h0.20 versus h0.16 checks apply only at R0.1/contact. '
            'Each contrast-change tolerance is max(5 eV, 10% of the absolute base coated-minus-bare contrast). '
            'Endpoint h/k qualification is separate from geometry force qualification; no full scientific qualification is claimed. '
            'Slab-thickness, cell/vacuum and full-radius h/k controls were omitted; no corresponding uncertainty bound or isolated-limit claim.\n\n'
            'Under entry_2971 the immutable precommit selects on force-qualified bare plus the actual accepted coated boundary: vacuum-side normal first if its paired clearance is >=1 A, otherwise the most-open common declared candidate on those CURRENT frames. '
            'The coated boundary is a geometric witness only: not final, not force-qualified, not partly-relaxed BO admission and not BO-eligible. Bare data retain final-coated admission pending. '
            'Actual final coated data require strict force-qualified, explicit entry_2904 partial eligibility, or distinct exact entry_2994/3001 labelled-original admission and continuous clearance >1e-6 A along the EXACT same precommitted direction; their label is **precommitted shared-open**, never maximum-common-final. The 1 A threshold is normal preference only, not a final-admission cutoff. '
            'A blocked final ray is refused without changing the bare ray. No global or final-frame optimum is claimed. Geometry-only records have no selected common ray. '
            'Positive-R data are not spliced to signed regular contact. No Ueq, action/Eq7, or experimental phase claim.\n\n'
            'Private inputs were checked by reader_pd.replay with the accepted Native-state, epoch, contact-operator, Green and energy helpers. '
            'This small package contains no GPWs, raw logs or original private source trees. Provider basenames and byte/SHA bindings are retained. '
            'Run `python3 export_pd.py verify DIRECTORY` for archive hashes, numerical replay from the saved energy ledger, and comparisons. '
            'That small-package check is not a fresh raw-log/checkpoint/operator-source/cap audit; use reader_pd.py on the original private case roots for that.\n')


def export_cases(roots, output, *, allow_partial=False, geometry_roots=None):
    geometry_roots = {} if geometry_roots is None else geometry_roots
    require((roots or geometry_roots) and set(roots) <= set(reader_pd.ROLES) and
            set(geometry_roots) <= {'plain', 'oxygen'}, 'Unknown/empty Pd case or geometry roles')
    roots = {role: Path(path).resolve() for role, path in roots.items()}
    geometry_roots = {role: Path(path).resolve() for role, path in geometry_roots.items()}
    output = Path(output).resolve()
    require(not output.exists() and all(not output.is_relative_to(root)
            for root in (*roots.values(), *geometry_roots.values())),
            'Require a new export outside every private input tree')
    cases = {role: reader_pd.replay(root) for role, root in roots.items()}
    require(all(run['schema'] == 'micro-pd-replay-v1' for run in cases.values()),
            'BO case roles require complete cases; use --geometry for geometry-only inputs')
    geometries = {role: reader_pd.replay_geometry(root) for role, root in geometry_roots.items()}
    for role, run in geometries.items():
        require(run['specification']['variant'] == role, 'Wrong named geometry variant')
        if role in cases:
            require(all(run['geometry_binding'][key] == cases[role]['specification']['selected_geometry']['binding'][key]
                        for key in ('bytes', 'sha256')), 'Separately supplied geometry is not the BO-selected original')
    comparisons = reader_pd.paired_comparisons(cases, allow_partial=allow_partial)
    output.mkdir(parents=True)
    # Saved state records are small JSON metadata, never checkpoint/log bytes.
    evidence = portable(cases)
    _new_json(output / 'replay-evidence.json', evidence)
    _new_json(output / 'geometry-evidence.json', portable(geometries))
    _new_json(output / 'paired-comparisons.json', portable(comparisons))
    _new_json(output / 'common-support.json', portable(_common_support(cases)))
    for variant in ('plain', 'oxygen'):
        if variant not in evidence:
            continue
        run = evidence[variant]
        name = variant + '_q1_BO_pair.json'
        _new_json(output / name, _curve(run, portable(comparisons)))
        _new_json(output / (variant + '_q1_contact.json'),
                  _contact(run, binding(output / name, name), portable(comparisons)))
        common = _common_support(cases)['common_radii_A']
        if common:
            _new_json(output / (variant + '_q1_BO_pair_common.json'),
                      _common_curve(run, portable(comparisons), common))
    with (output / 'README.md').open('x', encoding='utf-8') as stream:
        stream.write(readme(comparisons, geometries, cases))
    (output / 'support').mkdir()
    for name in ('hydride_replay.py', 'verify_contact_mapping.py', 'log_privacy.py'):
        shutil.copyfile(SUPPORT / name, output / 'support' / name)
    green = SUPPORT / 'green.py'
    if not green.is_file():
        green = SUPPORT.parent / 'delivery' / 'green.py'
    shutil.copyfile(green, output / 'support/green.py')
    for name in ('reader.py', 'reader_pd.py', 'reader_pd_chain.py', 'reader_pd_label.py',
                 'reader_pd_parallel.py', 'reader_pd_terminal.py', 'export_pd.py'):
        shutil.copyfile(HERE / name, output / name)
    manifest = dict(schema='micro-pd-small-delivery-v1', publishing=False, native_invocations=0,
                    roles=[r for r in reader_pd.ROLES if r in roots], allow_partial=allow_partial,
                    geometry_roles=sorted(geometries),
                    qualified=comparisons['qualified'],
                    paired_hk_qualified=comparisons['paired_hk_qualified'],
                    geometry_force_qualified=comparisons['geometry_force_qualified'],
                    geometry_error_bound_eV=None, full_scientific_qualification=False,
                    replay_scope='Saved-energy arithmetic and archive integrity; raw sources remain private',
                    files=[binding(p, p.relative_to(output).as_posix()) for p in sorted(output.rglob('*')) if p.is_file()])
    _new_json(output / 'manifest.json', manifest)
    return manifest


def _energy(state):
    energy = state['energy']
    close(energy['internal_energy_eV'], energy['free_energy_eV'] - energy['minus_TS_eV'], 'Small-package internal energy')
    return energy['internal_energy_eV']


def _replay_arithmetic(run):
    """No Native/source admission here: the original private replay is separate."""
    g, states = run['geometry'], run['states']
    green = PeriodicNuclearGreen(g['cell_A'])
    direction = reader_pd.np.asarray(g['direction'], dtype=float)
    rows = []
    occupied, vacant = _energy(states['occupied-host']), _energy(states['vacant-host'])
    radii = run['specification']['radii_A']
    require(len(states['occupied-probes']) == len(states['vacant-probes']) == len(radii), 'Small-package radius count')
    for radius, op, vp in zip(radii, states['occupied-probes'], states['vacant-probes']):
        q = math.fsum([_energy(op), -_energy(vp), -occupied, vacant])
        correction = green.point_minus_gaussian_eV(radius * direction, math.sqrt(5.) / BOHR_A)
        total = q + correction
        rows.append(dict(R_A=radius, BO_quartet_soft_eV=q, point_minus_gaussian_eV=correction,
                         V_pair_total_eV=total, lowering_eV=COULOMB_EV_A/radius-total))
    compare_json(run['rows'], rows, 'Small-package quartet/Green replay')
    q = math.fsum([_energy(states['occupied-probe-He']), -_energy(states['vacant-probe-H-R0']), -occupied, vacant])
    proof = run['mapping_proof']
    compare_json(proof['mapping'], regular_mapping(g['cell_A']), 'Accepted regular contact mapping replay')
    close(proof['unit_constant_residual_eV'], proof['B_native_eV'] - proof['mapping']['B_source_eV'],
          'Source/native unit residual arithmetic')
    signed = math.fsum([q, proof['unit_constant_residual_eV'], proof['mapping']['regular_images_eV']])
    close(run['Q_He_eV'], q, 'Small-package contact quartet')
    close(run['signedW0_eV'], signed, 'Small-package signed contact')


def replay_export(directory):
    root = Path(directory).resolve()
    manifest = read_json(root / 'manifest.json')
    require(manifest['schema'] == 'micro-pd-small-delivery-v1' and manifest['publishing'] is False and
            manifest['native_invocations'] == 0 and type(manifest['allow_partial']) is bool,
            'Wrong Pd small-package schema')
    names = [item['path'] for item in manifest['files']]
    require(len(names) == len(set(names)), 'Duplicate small-package member')
    actual = set()
    for path in root.rglob('*'):
        require(not path.is_symlink() and (path.is_dir() or path.is_file()), 'Unsafe small-package member')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
    require(actual == set(names) | {'manifest.json'}, 'Unbound/missing small-package member')
    for item in manifest['files']:
        bound_file(root, item, 'Small-package member')
    cases = read_json(root / 'replay-evidence.json')
    require(set(cases) == set(manifest['roles']), 'Small-package role inventory')
    compare_json(cases, portable(cases), 'No absolute private machine paths')
    for run in cases.values():
        _replay_arithmetic(run)
        compare_json(run['ray_policy'], reader_pd._ray_status(run['geometry']), 'Small-package typed ray status')
    geometries = read_json(root / 'geometry-evidence.json')
    require(set(geometries) == set(manifest['geometry_roles']) <= {'plain', 'oxygen'}, 'Small-package geometry role inventory')
    compare_json(geometries, portable(geometries), 'No private geometry paths')
    for role, run in geometries.items():
        require(run['schema'] == 'micro-pd-geometry-replay-v1' and run['qualified'] is False and
                run['specification']['variant'] == role and run['BO_available'] is False,
                'Geometry-only result misrepresented as BO')
        qualification = reader_pd._qualification(run['geometry'])
        for key, value in qualification.items():
            compare_json(run[key], value, 'Small-package geometry qualification ' + key)
    comparisons = reader_pd.paired_comparisons(cases, allow_partial=manifest['allow_partial'])
    compare_json(read_json(root / 'paired-comparisons.json'), comparisons, 'Small-package paired comparisons')
    compare_json(read_json(root / 'common-support.json'), _common_support(cases), 'Actual common-support curve')
    require(manifest['qualified'] is comparisons['qualified'], 'Small-package qualification changed')
    for key in ('paired_hk_qualified', 'geometry_force_qualified', 'geometry_error_bound_eV', 'full_scientific_qualification'):
        compare_json(manifest[key], comparisons[key], 'Small-package separate qualification ' + key)
    for variant in ('plain', 'oxygen'):
        if variant in cases:
            name = variant + '_q1_BO_pair.json'
            compare_json(read_json(root / name), _curve(cases[variant], comparisons), 'Small-package curve')
            compare_json(read_json(root / (variant + '_q1_contact.json')),
                         _contact(cases[variant], binding(root / name, name), comparisons), 'Small-package contact')
            common = _common_support(cases)['common_radii_A']
            if common:
                compare_json(read_json(root / (variant + '_q1_BO_pair_common.json')),
                             _common_curve(cases[variant], comparisons, common), 'Measured common-support curve')
    require((root / 'README.md').read_text() == readme(comparisons, geometries, cases),
            'Small-package scope/qualification README')
    return dict(schema='micro-pd-small-replay-v1', replayed=True, qualified=comparisons['qualified'],
                paired_hk_qualified=comparisons['paired_hk_qualified'],
                geometry_force_qualified=comparisons['geometry_force_qualified'], geometry_error_bound_eV=None,
                full_scientific_qualification=False, comparisons=comparisons, raw_source_audit=False, scope=reader_pd.SCOPE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    build = sub.add_parser('create')
    build.add_argument('--case', action='append', default=[], metavar='ROLE=DIRECTORY')
    build.add_argument('--geometry', action='append', default=[], metavar='VARIANT=DIRECTORY')
    build.add_argument('--output', required=True)
    build.add_argument('--allow-partial', action='store_true', help='Allow missing case files only; never waive geometry admission')
    verify = sub.add_parser('verify')
    verify.add_argument('directory')
    args = parser.parse_args()
    if args.command == 'verify':
        result = replay_export(args.directory)
    else:
        roots = {}
        for assignment in args.case:
            role, path = assignment.split('=', 1)
            require(role not in roots, 'Duplicate case role')
            roots[role] = path
        geometries = {}
        for assignment in args.geometry:
            role, path = assignment.split('=', 1)
            require(role not in geometries, 'Duplicate geometry role')
            geometries[role] = path
        result = export_cases(roots, args.output, allow_partial=args.allow_partial, geometry_roots=geometries)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
