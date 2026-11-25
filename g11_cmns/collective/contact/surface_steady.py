"""Covered PdD D28–D26 steady sensitivity. No lattice propagation or fit to yields."""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from contact import BULK_PAIR, KB, SitePair, attach_readout, solve
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
WALLS = ('terminal', 'near', 'stiff')


def descriptor(slab, bulk, wall):
    directions = slab['directions']
    site = {}
    for label in ('adsorbed', 'subsurface'):
        samples = [r for r in directions[label] if r['point']['axis'] == 'normal']
        site[label] = {r['point']['displacement_A']:r for r in samples}
    locations = [np.array(site[name][.12]['point']['reference_D_position_A']) for name in site]
    axis = locations[1]-locations[0]
    r0 = float(np.linalg.norm(axis))
    axis /= r0
    nz = abs(float(axis[2]))
    result = deepcopy(bulk)
    result.update(pair_id=f'covered-PdD-D28-D26-{wall}', R0_A=r0, axis=axis.tolist(),
                  preparation={'thermal_K':[300.]}, sites=[], wells={})
    for label, samples in site.items():
        index = samples[.12]['moving_index']
        f0 = slab['zero']['forces_eV_A'][index][2]
        near = (samples[-.12]['force_parallel_eV_A']-samples[.12]['force_parallel_eV_A'])/(.24)*nz*nz
        result['sites'].append({'id':f'D{index}', 'well':label})
        result['wells'][label] = {'source':'native signed slab normal; projected onto near-normal pair line'}
        for name, sign in (('inward', math.copysign(1., axis[2]) if label == 'adsorbed' else -math.copysign(1., axis[2])),
                           ('outward', -math.copysign(1., axis[2]) if label == 'adsorbed' else math.copysign(1., axis[2]))):
            q = [.12/nz, .36/nz]
            e = [samples[s*sign]['relative_free_energy_eV'] for s in (.12,.36)]
            slopes = [-sign*nz*samples[s*sign]['force_parallel_eV_A'] for s in (.12,.36)]
            slope0 = -sign*nz*f0
            terminal = float(CubicHermiteSpline([0.,*q], [0.,*e], [slope0,*slopes])(q[-1],2))
            # At origin_cut=0 K enters only the unmeasured exterior. All three
            # choices have exactly the same measured signed E and force fit.
            curvatures = {'terminal':max(.25*near, terminal), 'near':near, 'stiff':4*near}
            result['wells'][label][name] = dict(K_eV_A2=curvatures[wall], K4_eV_A4=0.,
                origin_cut_A=0., q_A=q, E_eV=e, dE_dq_eV_A=slopes, dE_dq_at_zero_eV_A=slope0,
                last_sample_curvature_eV_A2=terminal, near_secant_curvature_eV_A2=near)
    return result


def run_variant(pair, ecut):
    rows = []
    for potential in pair.potentials():
        out: dict[str, Any] = solve(potential, ecut=ecut)
        attach_readout(out)
        cold = next(row for row in out['observables'] if row['preparation']=='cold_ground')
        hot = next(row for row in out['observables'] if row['preparation']=='thermal_300K')
        e = np.asarray(out['spectrum']['energy_eV'])
        lc = np.asarray(out['spectrum']['log10_contact_A_minus3'])
        lw = -(e-e[0])/(300*KB)
        subset = e < ecut-1
        truncated_contact = (float(logsumexp((lw+lc*np.log(10))[subset]))
                             -float(logsumexp(lw[subset])))/np.log(10)
        rows.append(dict(join_A=potential.join, include_flagged=potential.include_flagged,
                         join_lattice_eV=float(potential.lattice(potential.join)),
                         ground_energy_eV=float(e[0]), ground=cold, thermal_300K=hot,
                         thermal_cut_check_change_dex=hot['log10_contact_A_minus3']-float(truncated_contact),
                         thermal_contact_peak_energy_eV=hot['contact_peak_energy_eV']))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('site_wells', type=Path, help='Exact delivered site-wells.json')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--ecut', type=float, default=3.)
    args = p.parse_args()
    source = args.site_wells.read_bytes()
    slab = json.loads(source)['surface']
    bulk = json.loads(BULK_PAIR.read_text())
    slab_descriptors = {wall:descriptor(slab,bulk,wall) for wall in WALLS}
    representative = slab_descriptors['near']
    distance = deepcopy(bulk)
    distance.update(R0_A=representative['R0_A'], pair_id='bulk-wells-short-distance-counterfactual',
                    preparation={'thermal_K':[300.]})
    cases = {}
    output = dict(site_wells_sha256=sha256(source).hexdigest(), bulk_pair_sha256=sha256(BULK_PAIR.read_bytes()).hexdigest(),
                  energy_cutoff_eV=args.ecut, cases=cases)
    with TemporaryDirectory() as d:
        path = Path(d)/'pair.json'
        for case, data in {'distance_only':distance, **slab_descriptors}.items():
            path.write_text(json.dumps(data))
            pair = SitePair(path)
            rows = run_variant(pair,args.ecut)
            cases[case] = dict(descriptor=data, variants=rows)
            print(case, [(r['join_A'],r['include_flagged'],r['ground']['log10_contact_A_minus3'],
                          r['thermal_300K']['log10_contact_A_minus3']) for r in rows], flush=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
