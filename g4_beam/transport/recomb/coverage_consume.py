"""Gas equilibrium and bulk-fed frozen-run-end-v sheet coverage; no material defaults."""
from coverage_inputs import make_model
from coverage_boundary import gas_equilibrated, supplied_sheet
from coverage_solve import converge_loaded


def consume_cases(*, source, rows, thermodynamics, capture_probability, vacancies, provenance):
    if not provenance:
        raise ValueError('Qualified electronic/harmonic thermodynamic provenance required')
    for row in rows:
        supplied_sheet(thermodynamics, row['T_K'], capture_probability(row['T_K']))
    gasrow = next(row for row in rows if not row['loaded'])
    loadedrow = next(row for row in rows if row['loaded'])
    gas = gas_equilibrated(make_model(source, gasrow, s_effective=1.),
        sheet=thermodynamics(gasrow['T_K']), capture_s0=capture_probability(gasrow['T_K']))
    capture = capture_probability(loadedrow['T_K'])
    bulk = converge_loaded(source, loadedrow, thermodynamics=thermodynamics,
        capture_s0=capture, vacancies=vacancies, omit_overflow=True)
    compare = None
    if bulk['refined']['near_surface_reaches_095']:
        compare = converge_loaded(source, loadedrow, thermodynamics=thermodynamics,
            capture_s0=capture, vacancies=vacancies, omit_overflow=False)
    return dict(thermodynamic_provenance=provenance, gas=gas,
        bulk_fed_frozen_v=bulk, overflow_term_on_comparison=compare,
        vacancy_provenance=vacancies.provenance,
        beam_stimulated_desorption=False,
        source_current_conditions=dict(j_cell_A_m2=5000., eta_abs=.05, j_rel=1., n_j=1.),
        omitted_processes=['beam implantation', 'beam vacancy production', 'sputter erosion',
            'beam adsorbate stripping', 'overflow in the primary thermal readout',
            'vacancy equilibration beyond the recorded run end'],
        retained_processes=['original bulk isotherm', 'original Darken flux',
            'finite electrochemical current supply', 'trap equilibration at saved v(z)'])
