"""Compare the saved ten-valence speed contrast with the actual beam consumer.

Portable saved adapter; source identities are in ledger.json. Locators use cwd.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

from g4_beam.screen.dynamic import dynamic_results
from g4_beam.screen.dynamic.dynamic_results import load_models, screening_eV
from g4_beam.forward import consumer, model, profile_thermal
from g4_beam.merit import api
from g4_beam.screen.dynamic._upstream import electronic_speed


SCREENING_SCOPE = (
    'r->0 limit of the accepted low-q ten-valence C/lambda Yukawa proxy; '
    'one isotropic Huke Eq10 factor using finite-width band-RMS velocity; '
    'not a full large-q/contact dielectric response, core correction or U18'
)
DYNAMIC_SCOPE = (
    'frozen-incident-speed approximation: both microscopic terms use the same '
    'nominal incident E_lab=2*E_cm, at every depth including the adsorbed sheet; '
    'not component-resolved slowing-down dynamics; molecular, spread and neutral '
    'components receive the nominal-setting shift through the unchanged consumer'
)
RATE_MAPPING = (
    'Ue_shift_keV(E_cm,z)=[U10_Pd(2*E_cm)-U10_PdD(2*E_cm)]/1000; '
    'Pd transported profiles/stopping/adsorbate with common PdD phenomenological '
    'calibration; no absolute microscopic screening added to a fit'
)
UNCERTAINTY = (
    'Pd k12/k24 velocity sensitivity is available in the saved controls but only '
    'the primary Pd/PdD models enter this candidate; no PdD velocity-grid control '
    'or global convergence claim; static proxies retain prior conditional '
    'numerical allowance, not a global physical error bound'
)


def source_binding(path):
    path = Path(path).resolve()
    raw = path.read_bytes()
    return dict(path=os.path.relpath(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def comparison(models):
    pd = api.reference('Pd')
    # The correction is depth-independent: evaluate once per incident energy,
    # not once per CSDA depth node and beam component.
    differences = (screening_eV(models['Pd'], 2*api.ENERGIES)
                   - screening_eV(models['PdD'], 2*api.ENERGIES)) / 1000
    shifts = dict(zip(map(float, api.ENERGIES), map(float, differences)))

    def shift(E_cm_keV, z_m):
        return shifts[float(E_cm_keV)]

    candidate = api.Target(
        name='transported_Pd_with_tenvalence_speed_contrast_PdD_calibration',
        profiles=pd.profiles, material='PdD', stopping_fn=pd.stopping_fn,
        columns=pd.columns, measurement=dict(pd.measurement), Ue_shift_keV=shift,
    )
    rows = api.compare(candidate)
    for row in rows:
        energy = 2 * row['E_cm_keV']
        u_pd = float(screening_eV(models['Pd'], energy))
        u_pdd = float(screening_eV(models['PdD'], energy))
        row.update(
            E_lab_keV=energy, U10_Pd_eV=u_pd, U10_PdD_eV=u_pdd,
            U10_Pd_minus_PdD_eV=u_pd-u_pdd,
            U10_loading_PdD_minus_Pd_eV=u_pdd-u_pd,
            Ue_shift_keV=(u_pd-u_pdd)/1000,
            microscopic_screening_scope=SCREENING_SCOPE,
            microscopic_dynamic_scope=DYNAMIC_SCOPE,
            microscopic_rate_mapping=RATE_MAPPING,
            microscopic_uncertainty=UNCERTAINTY,
        )
    return candidate, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--moments', type=Path, required=True,
                        help='Actual saved Pd/PdD ten-valence moment JSON')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    models, controls, moments_binding = load_models(args.moments)
    candidate, rows = comparison(models)
    bindings = {
        'beam_screen_results': source_binding(__file__),
        'dynamic_results': source_binding(dynamic_results.__file__),
        'electronic_speed': source_binding(electronic_speed.__file__),
        'merit_api': source_binding(api.__file__),
        'forward_consumer': source_binding(consumer.__file__),
        'forward_model': source_binding(model.__file__),
        'profile_thermal': source_binding(profile_thermal.__file__),
        'merit_readme': source_binding(api.HERE/'README.md'),
        'transported_profiles': source_binding(api.HERE/'profiles.json'),
        'transported_cause_fits': source_binding(api.HERE/'fits_transported.json'),
        'paper_error_bars': source_binding(api.BEAM/'forward'/'errbars.csv'),
    }
    result = dict(
        schema='pass5-tenvalence-frozen-incident-beam-comparison-v1',
        source_moments=dict(path=os.path.relpath(args.moments.resolve()), **moments_binding),
        source_bindings=bindings,
        source_bindings_scope='Runtime-bound listed source/data files; not a complete transitive environment manifest',
        models=models, velocity_controls_not_applied=controls,
        candidate=dict(name=candidate.name, profile_reference='Pd',
                       stopping_reference='Pd', adsorbate_reference='Pd',
                       calibration_material=candidate.material,
                       unchanged_profiles_stopping_adsorbate_measurement=True),
        screening_scope=SCREENING_SCOPE, dynamic_scope=DYNAMIC_SCOPE,
        rate_mapping=RATE_MAPPING, uncertainty=UNCERTAINTY,
        fitted_values_scope=(
            'All original compare fields retained. matched_{Pd,PdD,Ti,TiD}_'
            'Ue_effective_{0,1,2}_keV are the actual corresponding paper Eq7 '
            'refits supplied by compare, not the quoted paper Ue or microscopic '
            'screening. Ue_effective_{0,1,2}_keV are the candidate diagnostic '
            'Eq7 fits. Powers 0/1/2 are additional inverse-current corrections '
            'to the native cause reading; cause_fit_current_power retains that '
            'reading. Existing cause parameters are not refitted.'
        ),
        no_new_SCF=True, no_new_solid=True, no_core_correction=True,
        native_paths_private=False, row_count=len(rows), rows=rows,
    )
    payload = json.dumps(result, indent=2, allow_nan=False) + '\n'
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/'comparison.json').write_text(payload)
    api.write_csv(args.output_dir/'comparison.csv', rows)


if __name__ == '__main__':
    main()
