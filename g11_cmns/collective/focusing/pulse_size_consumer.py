"""Settled extracted consumer: exact selected arrays and conditional pulse units."""
from pathlib import Path
import hashlib
import json
import math
import sys

from inward_consumer import consume

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'contact'))
from g12_readout import _PROTON_A3_PER_S, event_readout


def close(a, b):
    assert abs(a-b) < 5e-11, (a, b)


def pulses():
    result = json.loads((HERE/'pulse/results.json').read_text())
    bounds = json.loads((HERE/'pulse/class_bounds.json').read_text())
    for name, expected in result['input_sha256'].items():
        path = HERE/name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, name
    assert len(result['pulses']) == len(bounds['pulses']) == 6
    summaries = []
    for pulse, ceiling in zip(result['pulses'], bounds['pulses']):
        assert pulse['case'] == ceiling['case']
        assert pulse['closing_peak_A_fs'] > 0
        assert pulse['closing_peak_fs'] < pulse['interpolated_turns_fs'][1]
        assert pulse['sampled_min_R_A'] == ceiling['sampled_min_R_A']
        if pulse['first_original_exit_fs'] is not None:
            assert pulse['interpolated_turns_fs'][-1] < pulse['first_original_exit_fs']
        assert len(pulse['variants']) == 4
        assert {(v['join_A'], v['include_flagged']) for v in pulse['variants']} == {
            (1.28782409436586, True), (1.28782409436586, False), (1.6, True), (1.6, False)}
        for variant in pulse['variants']:
            cap = next(v for v in ceiling['caps'] if
                (v['join_A'], v['include_flagged']) == (variant['join_A'], variant['include_flagged']))
            for name, cap_name in (('lower_band', 'lower'), ('upper_band', 'upper')):
                band = variant[name]
                assert band['max_norm_error'] < 2e-9
                assert 0 <= band['closing_fraction_of_Xi'] <= 1
                assert band['log10_total_Xi_A_minus3_s'] <= band['finite_subspace_log10_Xi_cap_A_minus3_s']+1e-10
                assert band['log10_max_instantaneous_probability_R_lt_1fm'] <= band['finite_span_log10_snapshot_cap_R_lt_1fm']+1e-10
                close(band['finite_subspace_log10_Xi_cap_A_minus3_s'], cap[cap_name]['log10_Xi_cap_A_minus3_s'])
                close(band['log10_probability_proton_per_pair_pulse'],
                      band['log10_total_Xi_A_minus3_s']+math.log10(_PROTON_A3_PER_S))
                close(band['log10_probability_either_DD_branch_per_pair_pulse'],
                      band['log10_probability_proton_per_pair_pulse']+math.log10(1+53.701/55.576))
                for field, value in event_readout(band['log10_total_Xi_A_minus3_s']).items():
                    close(band['one_pulse_per_incident_D_total_G12'][field], value)
                excess = band['signed_excess_Xi']
                if excess['sign'] > 0:
                    assert band['one_pulse_per_incident_D_excess_G12'] == event_readout(excess['log10_magnitude_A_minus3_s'])
                else:
                    assert band['one_pulse_per_incident_D_excess_G12'] is None
                assert band['rms_R_mismatch_A'] >= 0
            close(variant['upper_minus_lower_log10_Xi'],
                variant['upper_band']['log10_total_Xi_A_minus3_s']-variant['lower_band']['log10_total_Xi_A_minus3_s'])
        summaries.append(dict(case=pulse['case'], launch=pulse['launch'],
            inward_eV=pulse['closing_peak_K_radial_eV'], sampled_min_R_A=pulse['sampled_min_R_A'],
            log10_proton_pulse_probability_range=[min(v[b]['log10_probability_proton_per_pair_pulse']
                for v in pulse['variants'] for b in ('lower_band', 'upper_band')),
                max(v[b]['log10_probability_proton_per_pair_pulse']
                for v in pulse['variants'] for b in ('lower_band', 'upper_band'))]))
    return summaries


if __name__ == '__main__':
    classical = consume()
    conditional = pulses()
    print(json.dumps(dict(verdict='settled extracted pulse/size consumer passed',
        classical_rows=len(classical['rows']), conditional_pulses=conditional,
        limitations='conditional finite-state radial examples and support-class ceilings; no inferred physical pulse chance, pair population or foil rate'), indent=2))
