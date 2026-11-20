"""Project the delivered PdD 2-keV head-on snapshot without extending its support."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from g4_beam.forward.consumer import forward_yields
from g4_beam.merit.api import reference
from g6_pairlaw.script.material_pair import read_interaction
from g6_pairlaw.script.pairlaw import C

HERE = Path(__file__).resolve().parent
LAB_KEV = np.array([.5, 1., 1.5, 2., 2.5])
AMPLITUDE_A = .0768


def collision_rows(baseline, displaced, baseline_work, displaced_work):
    """Matched signed differences. Positive work raises the effective barrier."""
    r, v, base = read_interaction(baseline)
    rg, vg, gamma = read_interaction(displaced)
    keys = ('bare', 'energy_zero', 'periodic_to_isolated')
    if not np.array_equal(r, rg) or any(base['reference'][k] != gamma['reference'][k] for k in keys):
        raise ValueError('BO exports must share the target-relative grid and reference')
    cases = []
    for work, curve, bo in zip((baseline_work, displaced_work), (baseline, displaced), (v, vg)):
        if work['schema'] != 'pdd-periodic-point-pair-diagnostic-v1':
            raise ValueError('Require the anchored-work readout, not an excitation or drag table')
        if work['source_bindings'] != curve['source_bindings']:
            raise ValueError('Work and mapped BO source bindings describe different snapshots')
        case = sorted(work['rows'], key=lambda row: row['R_A'])
        if not np.array_equal([row['R_A'] for row in case], r):
            raise ValueError('Encounter and BO supports differ')
        if not np.allclose([row['V_pair_periodic_point_eV'] for row in case], bo, rtol=0., atol=1e-8):
            raise ValueError('Work readout and mapped BO energies describe different snapshots')
        cases.append(case)
    return [dict(R_A=float(x), base_BO_lowering_eV=ub, displaced_BO_lowering_eV=ug,
                 delta_BO_lowering_eV=float(v0-v1),
                 base_anchored_encounter_work_eV=b['anchored_encounter_work_profile_eV'],
                 displaced_anchored_encounter_work_eV=g['anchored_encounter_work_profile_eV'],
                 delta_encounter_work_eV=g['anchored_encounter_work_profile_eV']-b['anchored_encounter_work_profile_eV'],
                 delta_encounter_excitation_eV=g['encounter_excitation_change_eV']-b['encounter_excitation_change_eV'],
                 base_vacant_mean_loss_eV_A=b['vacancy_mean_loss_eV_A'],
                 displaced_vacant_mean_loss_eV_A=g['vacancy_mean_loss_eV_A'])
            for x, ub, ug, v0, v1, b, g in zip(r, base['barrier_lowering_eV'], gamma['barrier_lowering_eV'], v, vg, *cases)]


def conditional_yield_factors(delta_lowering_eV):
    """Counterfactual rigid Eq. 7 shift, not a projection of the finite curve."""
    target = reference()
    # Keep the profile and sheet fixed so sensitivity isolates the shift.
    # This is neither retransport nor an energy-dependent data transfer.
    profile = target.profiles
    if not callable(profile):
        raise TypeError('Transported reference must supply a callable profile')
    snapshot = profile(1.)  # 2 keV nominal lab, 1 keV nominal DD CM
    def run(shift_eV):
        result = forward_yields(LAB_KEV/2, lambda _: snapshot, Ue_keV=shift_eV/1000,
                                measurement={'channel': 'p'}, stopping_fn=target.stopping_fn)
        return np.asarray(result['physical_yield_per_C'], dtype=float)
    baseline = 40.
    ref = run(baseline)
    changed = run(baseline+delta_lowering_eV)
    return [dict(E_lab_Dplus_keV=float(e), base_physical_DDp_events_per_true_C=float(a),
                 hypothetical_physical_DDp_events_per_true_C=float(b),
                 hypothetical_yield_ratio=float(b/a)) for e, a, b in zip(LAB_KEV, ref, changed)]


def project(paths):
    raw = [path.read_bytes() for path in paths]
    rows = collision_rows(*map(json.loads, raw))
    near = rows[0]
    # Work is supplied to the system. At fixed entrance kinetic energy, an
    # increase consumes more kinetic energy: equivalent lowering has minus sign.
    hypothetical_shift = -near['delta_encounter_work_eV']
    warm_max = .72  # observed warm single-site maximum, not a Gamma amplitude
    return dict(
        schema='coherencev1-pdd-collision-consumer-v1',
        sources={k: sha256(b).hexdigest() for k, b in
                 zip(('baseline_BO', 'displaced_BO', 'baseline_work', 'displaced_work'), raw)},
        snapshot=dict(projectile='D', incident_lab_keV=2., direction='[110]', impact_A=0.,
                      target='ordered PdD', displacement_A=[AMPLITUDE_A, 0., 0.], fixed_host='Pd'),
        contrast_by_R=rows,
        preparation=dict(material_contact_rate_s=None, reason='No prepared well or short-range continuation in BO exports.'),
        beam=dict(qualified_yield_ratio=None, near_R_A=near['R_A'],
                  largest_bare_entrance_turning_radius_A=float(C/(1000*LAB_KEV[0]/2)),
                  reason='One head-on path at 2 keV; inner barrier outside sampled support at these monoatomic entrance energies.',
                  encounter_excitation_convergence='Not checked in dt. Vacant drag is separate and not used as a contact shift.'),
        counterfactual=dict(kind='rigid Eq. 7 sensitivity, not an inferred yield or fitted amplitude',
                            fixed_reference_Ue_eV=40., fixed_profile_E_lab_keV=2.,
                            assumed_uniform_delta_Ue_eV=hypothetical_shift,
                            yield_factors=conditional_yield_factors(hypothetical_shift)),
        amplitude=dict(
            BO_secant_eV_per_A=near['delta_BO_lowering_eV']/AMPLITUDE_A,
            beam_lowering_secant_eV_per_A=hypothetical_shift/AMPLITUDE_A,
            targets=[dict(positive_shift_eV=e, needed_positive_slope_at_test_A_eV_per_A=e/AMPLITUDE_A,
                          needed_positive_slope_at_warm_max_A_eV_per_A=e/warm_max,
                          formal_A_at_work_secant_magnitude=e*AMPLITUDE_A/abs(hypothetical_shift)) for e in (500., 2000.)],
            interpretation='One negative secant, not a response law or bound; warm local displacement is not coherent Gamma amplitude.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('collision_dir', type=Path)
    parser.add_argument('--output', type=Path, default=HERE/'collision_projection.json')
    args = parser.parse_args()
    source = args.collision_dir
    result = project((source/'curves/baseline_BO_pair.json',
                      source/'curves/gamma_BO_pair.json',
                      source/'receipts/collision-base-periodic-pair-readout-001.json',
                      source/'receipts/collision-gamma-periodic-pair-readout-001.json'))
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
