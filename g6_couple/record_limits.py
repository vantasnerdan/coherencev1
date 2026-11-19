"""Price the point-charge scenario against actual native-step amplitude maxima.

A measured maximum bounds the SAVED discrete record only, never an unobserved
thermal tail or a future trajectory. This does not compute electronic response.
"""
from pathlib import Path
import json
from response import spectator_bound

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    manifest = json.loads((HERE/'trajectories.json').read_text())
    rows = []
    for file in manifest['cases']:
        record = json.loads((HERE/file).read_text())
        case = record['summary']['configuration']
        p = json.loads((ROOT/record['parameters_file']).read_text())['materials'][case['material']]
        stats = record['statistics']
        cap = stats['maximum_observed_site_displacement_A']
        envelope = spectator_bound(cap, a=p['a_A'])
        rows.append(dict(case=record['summary']['label'], ion_species=case['material'],
            max_observed_site_displacement_A=cap,
            max_observed_joint_NN_relative_displacement_A=stats['maximum_observed_pair_displacement_A'],
            prior_snapshot_control_A=.0432,
            prior_control_amplitude_exceeds_observed_max=bool(cap<=.0432),
            primary_collision_snapshot_A=.0768,
            screening_scenario_A=.541, reference_a_A=p['a_A'],
            screening_status='PdD .541A conditional scenario' if case['material']=='PdD' else 'NiH screening not identified; .541A is a shared sensitivity scenario, NOT a material bound',
            single_record_vs_undisplaced_spectator_absolute_envelope_eV=envelope,
            signed_spectator_change_eV=None,
            definition='Triangle envelope using the actual maximum over every site and observed native step; same .541A screened unit-charge scenario. It is not an electronic response, an observed signed shift, or a stochastic tail bound.'))
    contrasts = []
    for off, on in (('NiH_cold_h0','NiH_cold_h0.08'),
                    ('PdD_random_x1_h0','PdD_random_x1_h0.08'),
                    ('PdD_thermal_x0.875_h0','PdD_thermal_x0.875_h0.08')):
        a = next(r for r in rows if r['case']==off)
        b = next(r for r in rows if r['case']==on)
        contrasts.append(dict(unpumped=off,pumped=on,
            absolute_pumped_minus_unpumped_envelope_eV=
                a['single_record_vs_undisplaced_spectator_absolute_envelope_eV']+
                b['single_record_vs_undisplaced_spectator_absolute_envelope_eV'],
            equal_pole_linear_susceptibility_contrast_eV=0.,
            warning='Sum of two global triangle envelopes; cancellations and simultaneity discarded. Pure common Gamma displacement still cancels exactly. No uniform locked state exists in these PdD records.'))
    value=dict(records=rows,contrasts=contrasts,
        imposed_PdD_snapshots=[dict(amplitude_A=a, spectator_absolute_envelope_eV=spectator_bound(a),
            meaning='Imposed pointwise .541A-screened charge scenario, not a computed electronic response') for a in (.0432,.0768)],
        snapshot_interpretation='The primary imposed electronic snapshot is now .0768A (rounded cold pumped PdD record maximum); .0432A is a budget-conditional linearity check. Neither is an RMS-to-maximum conversion or a bound on thermal/future states.',
        finite_shell_note='The17^3-site envelope is a finite screened-shell sum; omitted far shells are negligible at .541A screening on these geometries, not a host-electronic error bound.')
    (HERE/'record_limits.json').write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    print(json.dumps(value,indent=2),flush=True)


if __name__=='__main__':
    main()
