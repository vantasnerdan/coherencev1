"""Readout of the base frames from Pd draft 453374618, integrated at a31c95e."""
import argparse
import json
from pathlib import Path
from math import exp, isfinite

import numpy as np

from g7_micro.zr_calibration.action import convert as dd

INPUT_COMMIT = 'a31c95ebb9e722fce5c0937d43a2cafab7172133'
KERNEL_SHA = 'f6b690ca2862426c27b73e4d3b9d3d97f2bc6ece213c76e97978c0a468d71df0'
MANIFEST_SHA = '1512c2f774da850bae09f27e33bc4e8ce122ee60dae01bf5b25503f917941208'


def load(capsule):
    manifest_raw = (capsule / 'manifest.json').read_bytes()
    if dd.digest(manifest_raw) != MANIFEST_SHA or dd.digest(Path(dd.__file__).read_bytes()) != KERNEL_SHA:
        raise ValueError('require the pinned Pd input manifest and unmodified DD action kernel')
    manifest = json.loads(manifest_raw)
    identity = {v['path']: v for v in manifest['files']}

    def bound(rel):
        data = (capsule / rel).read_bytes()
        dd.bind(data, identity[rel])
        return data

    support = json.loads(bound('common-support.json'))
    pair = json.loads(bound('paired-comparisons.json'))
    if (support['interpolation'] or pair['geometry_error_bound_eV'] is not None or
            pair['paired_hk_qualified'] or pair['paired_cell_qualified']):
        raise ValueError('this readout is for the measured common-R, unqualified base snapshot')
    base = pair['pairs']['base']
    sites, metadata = {}, {}
    for role in ('plain', 'oxygen'):
        original_rel = f'{role}_q1_BO_pair.json'
        common_rel = f'{role}_q1_BO_pair_common.json'
        contact_rel = f'{role}_q1_contact.json'
        raw = bound(original_rel)
        c = json.loads(bound(contact_rel))
        # In this export the sidecar binds the original curve, not the common
        # projection. Both measured R sets are complete; verify their identity.
        dd.bind(raw, c['finite_curve'])
        original, v = json.loads(raw), json.loads(bound(common_rel))
        cfg = v['configuration']
        ray = v['ray_policy']['direction_n_hat']
        if (c['finite_curve']['path'] != original_rel or c['R_A'] != 0 or c['charge'] != 1 or
                c['schema'] != 'hydride-signed-contact-sidecar-v1' or
                c['observable'] != 'signed_regular_contact_limit_of_target_conditioned_BO_total' or
                c['derivative_eV_per_A'] is not None or not isfinite(c['signedW0_eV']) or
                c['electronic_branch'] != 'bare-no-added-electron' or ray != c['ray_policy']['direction_n_hat'] or
                cfg['probe_charge'] != 1 or cfg['probe_added_electrons'] != 0 or cfg['isotope'] != 'D' or
                v['R_A'] != original['R_A'] or v['R_A'] != support['common_radii_A'] or
                v['V_pair_total_eV'] != original['V_pair_total_eV'] or
                v['reference'] != original['reference'] or v['source_bindings'] != original['source_bindings'] or
                c['geometry_error_bound_eV'] is not None or v['geometry_error_bound_eV'] is not None):
            raise ValueError('original, measurement projection and signed contact do not agree')
        geometry = cfg['geometry']
        source = v['source_bindings']['original_case']
        # The Pd export gives source and configuration at top level, while the
        # older reader requires text inside reference. Adapt only metadata in
        # memory, from the already bound Pd facts. No arrays or gauge change.
        ref = dict(v['reference'], source=f"{source['path']} SHA-256 {source['sha256']}",
                   configuration=f"{c['case']}: {v['relaxation_status']}; q+1 D probe, frozen saved frame")
        r, total, _ = dd.read_interaction(dict(v, reference=ref))
        w = total-dd.C/r
        if (r[0] != .1 or abs(w[0]+base[f'{role}_R01_lowering_eV']) > 1e-6 or
                c['signedW0_eV'] != base[f'{role}_signedW0_eV'] or
                total.tolist() != [x[f'{role}_V_pair_total_eV'] for x in support['rows']]):
            raise ValueError('producer comparison is not the measured base pair')
        sites[role] = r, w, c, raw, original_rel
        metadata[role] = dict(source=source, finite_curve=identity[original_rel], common_curve=identity[common_rel],
                              contact=identity[contact_rel], W0_eV=c['signedW0_eV'], R_A=r.tolist(), W_eV=w.tolist(),
                              cell_A=geometry['cell_A'], target_A=geometry['positions_A'][geometry['target_index']],
                              geometry=pair['geometry_by_role'][role], ray=ray,
                              precommit_sha256=v['ray_policy']['precommit']['sha256'], raw_reference=v['reference'])
    p, o = metadata.values()
    if any(p[k] != o[k] for k in ('cell_A', 'ray', 'precommit_sha256', 'raw_reference')):
        raise ValueError('common comparison needs the same cell, precommitted ray and gauge')
    return sites, metadata, pair


def evaluate(capsule):
    sites, metadata, pair = load(capsule)
    w_curves = {name: {shape: dd.prepare(dd.potential(r, w, c['signedW0_eV'], shape)) for shape in dd.SHAPES}
                for name, (r, w, c, _, _) in sites.items()}
    energies = sorted({*np.geomspace(.05, 6.5, 121).tolist(), *dd.ANCHORS})
    rows = []
    for e in energies:
        for shape in dd.SHAPES:
            site = {}
            for name, samples in w_curves.items():
                a = dd.gain(samples[shape], e*1000.)
                status = dd.eq7(e*1000., a['G'] if a is not None else None)
                if a is None or status['U_eV'] is None:
                    raise ValueError(f'not supported for paired comparison: {e} {shape} {name}: {status}')
                site[name] = dict(a, **status)
            plain, oxygen = site['plain'], site['oxygen']
            dg = oxygen['G']-plain['G']
            rows.append(dict(E_cm_keV=e, shape=shape, **site, delta_G_O_minus_clean=dg,
                             delta_U_O_minus_clean_eV=oxygen['U_eV']-plain['U_eV'], coated_over_bare=exp(dg)))
    guards = dd.check(sites['plain'][0], sites, rows)
    return dict(schema='pd-DD-base-frame-action-readout-v1', input_commit=INPUT_COMMIT, producer_draft=453374618,
                manifest_sha256=MANIFEST_SHA, kernel_sha256=KERNEL_SHA, isotope='DD', Rmin_A=dd.R0,
                comparison='coated minus bare; WKB actions at the same finite nuclear radius, same Eq.7 S(E)',
                scope='base frames only, no geometry/electronic numerical bound; 10 keV CM extension separate',
                anchors_cm_keV=list(dd.ANCHORS), source=metadata, endpoint_comparison=pair, guards=guards, rows=rows)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('capsule', type=Path)
    p.add_argument('result', type=Path)
    args = p.parse_args()
    args.result.write_text(json.dumps(evaluate(args.capsule), indent=2, allow_nan=False)+'\n')
