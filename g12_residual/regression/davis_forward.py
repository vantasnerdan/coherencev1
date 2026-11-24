"""Replay the fixed G4/G8/G9 denominator for the published Davis Fig. 2 keys."""
import csv
import hashlib
import json
import math
import sys
from functools import cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / 'g8_why'), str(ROOT)]
import common as g8
from driven_target import KB, Motion
from reaction import resonance_parts

STUDY = 'karahadian2026'
DOI = 'https://doi.org/10.1038/s41467-026-74421-1'
OBS = 'g4_beam/forward/errbars.csv'
PROFILES = 'g4_beam/merit/profiles.json'
GATED = 'g9_edge/output/yields.csv'
FORMED = 'g9_edge/output/sensitivity.csv'
SMOOTH_CONTROL = 'g8_why/baseline.csv'
RESONANCE_CONTROL = 'g8_why/reaction_retained.csv'
THERMAL_CONTROL = 'g8_why/driven_target_floats.csv'
BASELINE = 'g12_davis_fixed_smooth40_ordinary_thermal_formed_sheet_no_join_k0'
FORMED_VARIANT = 'sheet_formed_bulk_wake'
REF = 6.5
JOIN_E = 4.5
GRID = tuple(float(e) for e in g8.ENERGIES)
VARIANTS = ('smooth40', 'no_edge', 'formed_sheet', 'gated_atom', 'threshold40meV', 'smooth40_join')
MOTION = Motion(48)


def csv_rows(relative):
    with (ROOT / relative).open(newline='') as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith('#')))


def row_key(row):
    return row['material'], float(row['E_cm_keV']), row['channel']


def channel(e):
    if e <= 4.0:
        return 'p'
    if e >= JOIN_E:
        return 'n'
    raise ValueError(f'No Fig. 2 channel assigned at {e} keV CM')


def reference_key(target):
    return target, REF, 'n'


def reference_group(target):
    return f'{STUDY}:fig2:{target}:n:{REF:.2f}'


def focus(e):
    return float(g8.DD.collected_fraction(e, 1., 0.))


def cup_current(e):
    return float(g8.DD.I_cup(e))


@cache
def smooth(target, e, ch):
    return g8.yield_C(target, e, ch)


@cache
def resonant(target, e, ch):
    # Sum the existing G8 smooth, interference and pole contractions at a=1.
    # The saved fitted "resonance_common" strength is NOT the literal source law.
    return float(sum(resonance_parts(target, e, ch))) if ch == 'p' else smooth(target, e, ch)


@cache
def settings(target, e):
    layers, sheet, state = g8.reference(target).profiles(e)
    return layers, float(sheet), float(state['T_K']), float(state['beam_time_s'])


@cache
def thermal(target, e, ch):
    _, _, temp, _ = settings(target, e)
    return MOTION.yield_C(target, e, KB * temp, ch)


@cache
def smooth40_join(target):
    return smooth(target, JOIN_E, 'n') / smooth(target, JOIN_E, 'p')


def readout_signal(yield_per_c, e, power):
    return yield_per_c / focus(e) / cup_current(e)**power


def normalize(yield_per_c, key, ref_yield_per_c, ref_key, power=0):
    target, e, _ = key
    if ref_key != reference_key(target):
        raise ValueError('Expected the same target’s 6.5-keV neutron reference')
    denom = readout_signal(ref_yield_per_c, REF, power)
    if not math.isfinite(denom) or denom <= 0.:
        raise ValueError(f'Nonpositive or nonfinite model reference: {ref_key}')
    return readout_signal(yield_per_c, e, power) / denom


def raw_variants(row):
    smooth_yield = float(row['model_smooth40_raw_reactions_per_true_C'])
    thermal_yield = smooth_yield * float(row['ordinary_thermal_raw_factor'])
    formed_yield = thermal_yield * float(row['formed_sheet_raw_factor'])
    return dict(smooth40=smooth_yield, no_edge=thermal_yield, formed_sheet=formed_yield,
                gated_atom=thermal_yield * float(row['gated_atom_raw_factor']),
                threshold40meV=formed_yield * float(row['threshold_reaction_raw_factor']),
                smooth40_join=formed_yield)


def add_clean_readouts(rows, anchors):
    references = {a['target']: a for a in anchors}
    assert set(references) == set(g8.MATERIALS) and len(anchors) == len(references)
    raw_reference = {}
    for target, anchor in references.items():
        assert anchor['channel'] == 'n' and float(anchor['E_cm_keV']) == REF
        variants = raw_variants(anchor)
        assert variants['formed_sheet'] == variants['threshold40meV']
        raw_reference[target] = variants
        anchor['normalizers'] = {name: dict(
            raw_neutron_reactions_per_true_C=y,
            k0_focus_corrected_signal=readout_signal(y, REF, 0),
            k1_focus_corrected_signal=readout_signal(y, REF, 1),
        ) for name, y in variants.items()}
        anchor.update(R_model=1., R_model_k1=1., log_residual=None, log_residual_k1=None)
    for row in rows:
        target, e, ch = key = row['target'], float(row['E_cm_keV']), row['channel']
        ref_key = (target, float(references[target]['E_cm_keV']), references[target]['channel'])
        raw = raw_variants(row)
        ref = raw_reference[target]
        modeled = {name: normalize(y, key, ref[name], ref_key) for name, y in raw.items()}
        if ch == 'p':
            modeled['smooth40_join'] *= float(row['sensitivity_smooth40_pn_join'])
        r_model = modeled['formed_sheet']
        row.update(model_corrected_raw_reactions_per_true_C=raw['formed_sheet'], R_model=r_model,
                   R_model_k1=normalize(raw['formed_sheet'], key, ref['formed_sheet'], ref_key, 1),
                   log_residual=math.log(float(row['R_observed']) / r_model))
        for name in VARIANTS:
            if name != 'formed_sheet':
                row[f'R_model_{name}'] = modeled[name]
                row[f'delta_log_residual_{name}_minus_primary'] = math.log(r_model / modeled[name])
        assert float(row['sigma_log_residual']) > 0. and math.isfinite(row['log_residual'])
        assert math.isclose(float(row['R_observed_k1']) / row['R_model_k1'],
                            float(row['R_observed']) / r_model, rel_tol=4e-15)
    key = ('Pd', .25, 'p')
    for wrong in (('PdD', REF, 'n'), ('Pd', REF, 'p'), ('Pd', 2.5, 'n')):
        try:
            normalize(1., key, 1., wrong)
        except ValueError as exc:
            assert 'same target' in str(exc)
        else:
            raise AssertionError(f'Accepted incompatible normalizer: {wrong}')
    return dict(incompatible_references_rejected=['wrong_target', 'wrong_channel', 'wrong_energy'],
                exact_digitized_key_count=len(rows), deterministic_anchor_count=len(anchors),
                matched_k1_common_current_factor_cancels=True,
                historical_k1_shift=add_current_reading(rows))


def edge_table(path, variant=None):
    rows = [row for row in csv_rows(path) if variant is None or row['variant'] == variant]
    grid = {(m, e, ch) for m in g8.MATERIALS for e in GRID for ch in ('p', 'n')}
    table = {row_key(row): row for row in rows}
    if len(rows) != len(grid) or set(table) != grid:
        raise ValueError(f'Unexpected G9 {path} grid or duplicate keys')
    if {row['response_proxy'] for row in rows} != {'pure_Pd_RPA_head_plus_Mermin'}:
        raise ValueError(f'Unexpected G9 proxy in {path}')
    return table


def edge_factor(table, key):
    target, e, ch = key
    if key in table:
        return float(table[key]['raw_yield_ratio']), 'saved_G9_node', e, e
    idx = int(np.searchsorted(GRID, e))
    if idx == 0 or idx == len(GRID):
        raise ValueError(f'Outside G9 support: {key}')
    lo, hi = GRID[idx - 1:idx + 1]
    a, b = (math.log(float(table[target, x, ch]['raw_yield_ratio'])) for x in (lo, hi))
    return math.exp(a + (e - lo) / (hi - lo) * (b - a)), 'interpolated_G9_log_ratio', lo, hi


def identity_fields(target, energy_text):
    e = float(energy_text)
    host = g8.DD.MATERIALS[target]['host']
    loaded = target.endswith('D')
    return dict(
        row_id=f'{STUDY}:fig2:{target}:{energy_text}', study_id=STUDY, lab='UC Davis', run_id=None,
        reaction='D+D', channel=channel(e), target=target, host=host, E_cm_keV=energy_text,
        target_preparation='beam_plus_electrochemical' if loaded else 'beam_only',
        host_Z=g8.DD.ELEM[host]['Z'], electrochemical_loading=int(loaded), foil_thickness_um=250,
        energy_frame_source='Fig. 2 nominal CM; published acceleration-potential convention',
        observable_kind='own_reference_relative_yield', observable_source=OBS, observable_unit='1',
        reference_E_cm_keV=REF, reference_group=reference_group(target),
        evidence_class='published_fig2_digitization', citation=DOI,
    )


def verify_raw():
    """Match raw contractions to saved controls without using historical fitted readouts."""
    def relative(actual, expected):
        return abs(actual / float(expected) - 1.)

    smooth_raw, reaction, motion, edge_raw, edge_ratio = ([] for _ in range(5))
    for row in csv_rows(SMOOTH_CONTROL):
        target, e = row['material'], float(row['E_cm_keV'])
        for ch in ('p', 'n'):
            smooth_raw.append(relative(smooth(target, e, ch), row[f'smooth40_{ch}_reactions_per_true_C']))
    for row in csv_rows(RESONANCE_CONTROL):
        target, e = row['material'], float(row['E_cm_keV'])
        p = target, e, 'p'
        n = target, e, 'n'
        reaction.extend((
            relative(resonant(*p) / smooth(*p), row['proton_same_E_factor']),
            relative((resonant(*p) + smooth(*n)) / (smooth(*p) + smooth(*n)), row['total_DD_same_E_factor']),
        ))
    for row in csv_rows(THERMAL_CONTROL):
        if row['state'] == 'thermal' and row['current_power'] == '0':
            target, e = row['material'], float(row['E_cm_keV'])
            for ch in ('p', 'n'):
                motion.append(relative(thermal(target, e, ch) / smooth(target, e, ch), row[f'physical_{ch}_factor']))
    for table in (edge_table(FORMED, FORMED_VARIANT), edge_table(GATED)):
        for key, row in table.items():
            edge_raw.append(relative(smooth(*key), row['smooth40_reactions_per_true_C']))
            edge_ratio.append(relative(float(row['total_reactions_per_true_C']) / smooth(*key), row['raw_yield_ratio']))
    checks = {}
    for name, errors in (('smooth_raw', smooth_raw), ('literal_resonance_raw', reaction), ('ordinary_thermal_raw', motion),
                         ('g9_smooth_raw', edge_raw), ('g9_raw_ratio', edge_ratio)):
        worst = max(errors)
        assert worst < 4e-13, (name, worst)
        checks[name] = dict(comparisons=len(errors), max_relative_error=worst)
    return checks


def add_current_reading(rows):
    """Historical G7/G8 k1: change the model readout, keep plotted R fixed."""
    plotted = {identity_fields(r['material'], r['E_cm'])['row_id']: float(r['R'])
               for r in csv_rows(OBS) if float(r['E_cm']) < REF}
    assert len(rows) == len(plotted) and {row['row_id'] for row in rows} == set(plotted)
    shifts = []
    for row in rows:
        r_obs = float(row['R_observed'])
        assert r_obs == plotted[row['row_id']]
        r_model_k1 = float(row['R_model_k1'])
        expected = math.log(cup_current(float(row['E_cm_keV'])) / cup_current(REF))
        y_k1 = math.log(r_obs / r_model_k1)
        shift_error = abs((y_k1 - float(row['log_residual'])) - expected)
        assert shift_error < 1e-14, (row['row_id'], shift_error)
        row['log_residual_k1'] = y_k1
        shifts.append(shift_error)
    return dict(rows=len(rows), unchanged_plot_observation=True, max_abs_shift_error=max(shifts))


def generate():
    observed = csv_rows(OBS)
    formed, gated = edge_table(FORMED, FORMED_VARIANT), edge_table(GATED)
    rows, anchors = [], []
    for obs in observed:
        target, e = obs['material'], float(obs['E_cm'])
        fields = identity_fields(target, obs['E_cm'])
        key, ref_key = (target, e, fields['channel']), reference_key(target)
        r_obs = float(obs['R'])
        y = smooth(*key)
        f_formed, support, lo, hi = edge_factor(formed, key)
        f_gated, _, _, _ = edge_factor(gated, key)
        factors = dict(model_smooth40_raw_reactions_per_true_C=y,
                       primary_pn_join=1., sensitivity_smooth40_pn_join=smooth40_join(target),
                       cup_focusing_fraction=focus(e), threshold_reaction_raw_factor=resonant(*key) / y,
                       ordinary_thermal_raw_factor=thermal(*key) / y,
                       formed_sheet_raw_factor=f_formed, gated_atom_raw_factor=f_gated)
        if key == ref_key:
            assert obs['method'] == 'ref' and r_obs == 1. and float(obs['half']) == float(obs['sigma_lnR']) == 0.
            anchors.append(dict(**fields, **factors, baseline_id=BASELINE,
                                status='deterministic_own_reference_not_fit_row',
                                R_observed=r_obs, sigma_log_residual=None, edge_energy_support=support))
            continue
        layers, sheet, temp, beam_time = settings(target, e)
        depth = float(g8.csda_nodes(2 * e, layers, g8.transport_stopping)[3]) * 1e9
        rows.append(dict(**fields, observable=obs['R'], sigma_minus=obs['half'], sigma_plus=obs['half'],
                         digitization_method=obs['method'], R_observed=r_obs, sigma_log_residual=float(obs['sigma_lnR']),
                         baseline_id=BASELINE, profile_id=f'{PROFILES}:{target}:{obs["E_cm"]}',
                         **factors, reaction_primary_case='smooth_Bosch_Hale',
                         threshold_sensitivity_case='external_literal_40meV_proton' if key[2] == 'p' else 'smooth_BH_neutron',
                         thermal_case='G8_ordinary_thermal_saved_T',
                         edge_case=FORMED_VARIANT, edge_response_proxy=formed[target, lo, key[2]]['response_proxy'],
                         edge_energy_support=support, edge_support_lo_E_cm_keV=lo, edge_support_hi_E_cm_keV=hi,
                         R_observed_k1=r_obs * cup_current(REF) / cup_current(e),
                         chart_cup_current_uA=cup_current(e), chart_current_status='published_SI_FigS3_setting_proxy',
                         model_beam_time_s=beam_time, model_target_temperature_K=temp,
                         model_csda_Dplus_depth_nm=depth, model_adsorbed_D_sheet_m2=sheet,
                         run_current_uA=None, run_time_s=None, run_dose_C=None,
                         run_covariate_status='unavailable_at_Fig2_row_level'))

    verification = add_clean_readouts(rows, anchors)
    verification['saved_raw_controls'] = verify_raw()
    emit(rows, anchors, verification)


def emit(rows, anchors, verification):
    provenance = metadata(rows, anchors, verification)
    with (OUT / 'davis_residuals.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    provenance['residual_csv_sha256'] = hashlib.sha256((OUT / 'davis_residuals.csv').read_bytes()).hexdigest()
    (OUT / 'davis_baseline.json').write_text(json.dumps(provenance, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(baseline_id=BASELINE, observed_rows=len(rows), anchors=len(anchors),
                          **provenance['summary'], verification=verification), indent=2))


def exposure_boundary(rows):
    """Observed-grid model change, not a time-controlled comparison."""
    by_key = {(r['target'], float(r['E_cm_keV'])): r for r in rows}
    results = []
    for target in g8.MATERIALS:
        left, right = (by_key[target, e] for e in (2., 2.5))
        t0, t1 = (float(r['model_beam_time_s']) for r in (left, right))
        assert (t0, t1) == (43200., 7200.)
        m0, m1 = (math.log(float(r['R_model'])) for r in (left, right))
        results.append(dict(target=target, E_cm_keV=[2., 2.5], nominal_beam_time_s=[t0, t1],
                            ln_R_model_k0=[m0, m1], delta_ln_R_model_k0=m1 - m0,
                            log_residual_k0=[float(r['log_residual']) for r in (left, right)]))
    return results


def metadata(rows, anchors, verification):
    sources = (OBS, PROFILES, 'g4_beam/forward/consumer.py', 'g4_beam/forward/model.py',
               'g4_beam/forward/profile_thermal.py', 'g4_beam/merit/api.py', 'g7_micro/analyze.py', 'g8_why/common.py',
               'g8_why/reaction.py', 'g8_why/reaction_source_notes.md',
               'g8_why/driven_target.py', 'g8_why/driven_target.md', 'g8_why/beam.md',
               'g10_quantum_site/replay.md', GATED, FORMED, 'g9_edge/replay.md',
               'g9_edge/first_principles/README.md', SMOOTH_CONTROL, RESONANCE_CONTROL, THERMAL_CONTROL)
    return dict(
        baseline_id=BASELINE, input_commit='bf114e0874f2a01dc62505e7c7520cadf46fd20f',
        source_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources},
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        observed=dict(
            citation=DOI, source=OBS, kind='dimensionless Fig. 2 relative yield as plotted (k0), not raw counts/C or S(E)',
            sigma_linear='source digitized mean symmetric half-width, copied to sigma_minus and sigma_plus',
            sigma_log_residual='source sigma_lnR: digitized relative bars, first-order log scale; observed only',
            shared_reference='one nuisance per foil; uncertainties of individual points and reference not separable from the published digitization; covariance/magnitude unknown, do not add a fabricated common variance',
            missing_keys=['PdD:6.00', 'TiD:3.50', 'TiD:4.00', 'TiD:5.50'],
            run_id='not published for these figure rows; row_id is a digitization ID, not a run join',
        ),
        transport=dict(
            source=PROFILES, profile_key=['material/condition', 'exact nominal CM energy'],
            reactions_unit='reactions per true charged coulomb',
            current_fractions={'D+': .85, 'D2+': .15, 'D3+': 0.}, species_atoms_per_charge=1.15,
            relative_energy_FWHM=.05, spread_hermite_nodes=5, csda_legendre_nodes_per_interval=12,
            lab_cutoff_keV=.1, static_screening_keV=.040,
            screening_radial_form='Eq. 7 constant energy-shift approximation, not a resolved radial V(r)',
            nominal_exposure='12 h for E_cm <= 2.0 keV; 2 h above; 1 h preload on loaded foils. Same saved profile throughout; discrete denominator change across 2.0→2.5 also includes energy effects, not an isolated time effect.',
            sheet='stored adsorbed-D column included once at each incident quadrature energy',
            high_energy_neutrals='not adopted: G8 tail sensitivities are not an established Davis beam component',
        ),
        normalization=dict(
            channel='proton at E<=4.0, neutron at E>=4.5; observed keys only',
            join='Primary p and n multipliers are both 1: no p/n amplitude imposed. Sensitivity only: multiply p by the same foil smooth40 raw n/p at 4.5 keV, n unchanged. No same-Davis fitted stitch.',
            focus='g4 collected_fraction(E,1,0), divide numerator and denominator each once; is a retained focusing convention, not measured detector acceptance',
            reference='same foil 6.5 keV n, with same baseline/variant as its numerator',
            formula='S_ch(E)=Y_ch(E)/c(E)/I_cup(E)^k; primary R_model=S_ch(E)/S_n(6.5)',
            primary_current_power=0,
            k0_paper_basis='Methods, Yield extraction and screening metrics: “All yields were reported per delivered Coulomb.”',
            k1='Historical G7/G8 sensitivity, not an established correction. Keep R_observed unchanged; R_model_k1 = R_model × I_cup(6.5)/I_cup(E). Thus log_residual_k1 = ln(R_observed/R_model_k1) = log_residual + ln(I_cup(E)/I_cup(6.5)). The observation log error and raw model physics are unchanged.',
            k1_paper_context='Dose was inferred from the Faraday-cup I–V characteristic, not per-run measured charge; the relative-yield equation is written in terms of a net-count ratio. This motivates a sensitivity, but does not establish a missing current division.',
            paired_axis_control='R_observed_k1 is the old paired-axis control (plotted R rescaled by I_cup(6.5)/I_cup(E)); comparing it to R_model_k1 cancels the factors and is NOT the historical k1 residual.',
            log_residual='natural log of R_observed/R_model, one ratio; observation bars do not include baseline sensitivities',
        ),
        layers=dict(
            smooth40='fixed G8 response, Bosch-Hale p and n',
            threshold_reaction='Not in primary. Named sensitivity uses the independent external Zr-derived, source-literal Eq. 17 proton amplitude: Gamma_p=0.040 eV, phase 115 degrees, radius 7 fm, ER→0; neutrons smooth. The Davis-fitted resonance_common and basis_common are not consumed.',
            ordinary_thermal='single G8 Motion(48) Maxwell target-rate case theta=k_B*saved conditional T, with relative flux and smooth40; ratios computed at every row and n reference. No thermal_ZPM, collision harmonic, g10 bound-state momentum, or hot100/hot1000 terms stacked. The G10 joint spectrum is a PdD case, not a four-foil replacement.',
            approach_edge=dict(
                adopted=FORMED_VARIANT, source=FORMED, alternate='gated_atom_at_t0', alternate_source=GATED,
                evidence='Native 0.25 keV CM nominal clean PdD(111) ray develops an electronic cloud before crossing the plane, motivating a nonzero sheet prehistory. It stops 0.5 Å short of the buried D center. It supplies no exact yield or evidence of a fully formed bulk wake in all foils. The adopted model changes the G9 surface-sheet prehistory only; it is not a solved vacuum approach.',
                proxy='All G9 foils use pure-Pd bulk RPA inverse-head response plus Mermin continuation; Ti/TiD transport profiles differ but no Ti dielectric calculation is used.',
                support='exact saved ratio only at G9 nodes; between nodes interpolate log(raw ratio) linearly in nominal E_cm within the same foil/channel; bracket recorded per row. No extrapolation, no channel averaging, not an interpolation error bound.',
                nodes_E_cm_keV=GRID,
            ),
            composition='Raw primary yield = smooth40 times ordinary thermal factor times formed-sheet G9 approach factor. Factors are each relative to the same smooth40 setting and channel. Factorizing yield ratios is not coupled reweighting; overlap of the thermal velocity, edge impulse and barrier response remains unresolved. G9 replaces its finite edge subtotal and subtracts settled projectile drag already in CSDA; no second drag, sheet or static barrier correction. Same product on n6.5 before normalization.',
        ),
        features=dict(
            E_cm_keV='nominal center-of-mass setting in keV; use as an independent design axis',
            host_Z='nominal host atomic number (not fitted); encodes the same contrast as host',
            electrochemical_loading='0 or 1, nominal preparation only, not a measured D/metal ratio',
            foil_thickness_um='source methods nominal common 250 μm membrane; no variation',
            chart_cup_current_uA='g4 literal SI Fig. S3 I–V interpolation at V_lab[kV]=2*E_cm[keV]; microamperes. Common across foils at a setting; not per-run current or an independent energy contrast.',
            model_beam_time_s='merit profiles exact nominal beam history in seconds; not actual beam-on time or detector live time',
            model_target_temperature_K='merit profiles conditional target temperature in kelvin; not a measured foil temperature',
            model_csda_Dplus_depth_nm='terminal depth in nm from fixed CSDA for monoenergetic nominal D+ at lab energy 2*E_cm to 0.1 keV lab cutoff or layer end; not SRIM range or a measured reaction-depth moment',
            model_adsorbed_D_sheet_m2='saved number of D atoms per square meter in modeled adsorbed sheet; not measured surface loading',
            missing='Blank run_current_uA, run_time_s, run_dose_C and run_id mean unavailable per Fig. 2 row, not zero. No imputed spot, flux or loading density.',
            eligibility='The chart and transport proxies describe settings or the denominator; not new independent measured variation. Do not regress on baseline terms, response-derived quantities, fitted U, or the plotted observation itself.',
        ),
        reference_anchors=anchors,
        nominal_exposure_boundary=exposure_boundary(rows),
        summary=dict(
            log_residual_min=min(float(r['log_residual']) for r in rows),
            log_residual_max=max(float(r['log_residual']) for r in rows),
            log_residual_k1_min=min(r['log_residual_k1'] for r in rows),
            log_residual_k1_max=max(r['log_residual_k1'] for r in rows),
            max_abs_log_residual_delta={name: max(abs(float(r[f'delta_log_residual_{name}_minus_primary']))
                                                 for r in rows)
                                        for name in VARIANTS if name != 'formed_sheet'},
        ),
        verification=verification,
    )


if __name__ == '__main__':
    generate()
