"""Join public-paper covariates to the adopted Fig. 2 digitization; no imputation."""
import argparse
import csv
from decimal import Decimal
from pathlib import Path


DOI = 'https://doi.org/10.1038/s41467-026-74421-1'
# Published supplement Table 1: E_cm, common calibration I, and per-foil S:B.
TABLE1 = {Decimal(e): tuple(map(int, values.split())) for e, values in {
    '0.25': '38 21 30 15 56', '0.5': '68 15 32 47 62', '1': '110 44 43 69 112',
    '1.5': '123 53 92 106 261', '2': '134 83 138 82 368', '2.5': '138 121 177 225 1082',
    '3': '145 443 465 1525 5427',
}.items()}
# SI Fig. 3 experimental marker centers rounded to 1 µA; calibrated from the axes.
FIG3_HIGH = {7: 156, 8: 165, 9: 176, 10: 182, 11: 207, 12: 222, 13: 239}
TARGETS = ('Pd', 'PdD', 'Ti', 'TiD')
SCREENING = dict(zip(TARGETS, ((.8, 3.28), (1.7, 3.68), (.3, 3.21), (.9, 2.75))))


def read_points(path):
    with path.open() as fd:
        return list(csv.DictReader(line for line in fd if not line.startswith('#')))


def point(row):
    target = row['material']
    host = target.removesuffix('D')
    loaded = target.endswith('D')
    energy = Decimal(row['E_cm'])
    reference = row['method'] == 'ref'
    voltage = 2 * energy
    low = energy <= 2
    screen_bulk, screen_surface = SCREENING[target]
    if energy in TABLE1:
        current, *signal_to_background = TABLE1[energy]
        current_source = 'SI Table 1 (cup I–V characterization)'
        signal_to_background = signal_to_background[TARGETS.index(target)]
    else:
        current = FIG3_HIGH.get(voltage, '')
        current_source = 'SI Fig. 3 marker (digitized, approx. 1 µA)' if current else ''
        signal_to_background = ''
    # The combined Fig. 2 series never gives a per-row run or channel ledger.
    channel = ('p+t' if energy < 3 else
               'n+3He (reported acquisition regime)' if energy > Decimal('4.5') else
               'p+t / n+3He (point mix unreported)')
    return dict(
        row_id=f'karahadian2026:fig2:{target}:{energy:.2f}', study_id='karahadian2026', lab='UC Davis', run_id='',
        reaction='D+D', channel=channel, target=target, host=host, target_geometry='thick',
        target_preparation=f'Vacuum annealed and same-host vacuum-side sputter renewed; {"Pd catalytic coating on electrolyte side" if host=="Ti" else "Pd membrane"}',
        loading_state='beam_plus_electrochemical' if loaded else 'beam_only',
        E_cm_keV=row['E_cm'], E_lab_nominal_Dplus_keV=str(voltage), energy_frame_source='Fig. 2 E_cm (stationary D)',
        observable_kind='relative_thick_target_yield', observable=row['R'], observable_unit='1',
        sigma_minus='' if reference else row['half'], sigma_plus='' if reference else row['half'],
        adopted_sigma_ln_R='' if reference else row['sigma_lnR'], digitization_method=row['method'],
        evidence_class='normalization_reference' if reference else 'figure_digitized_measurement',
        reference_E_cm_keV='6.50', reference_group=f'karahadian2026:fig2:{target}:6.50',
        citation=f'{DOI}; Fig. 2{"a" if host=="Pd" else "b"}; Methods; SI Figs. 2–3 and Tables 1–2',
        temperature_K='', actual_D_per_metal='', per_run_beam_current_uA='', current_density_uA_cm2='',
        run_duration_h='', dose_C='', accumulated_dose_C='',
        cup_IV_current_uA=current, cup_IV_current_source=current_source,
        cup_IV_repeat_spread_approx_uA=5 if current else '', cup_IV_nominal_voltage_kV=str(voltage) if current else '',
        cup_IV_bias_regime='internal_only' if voltage <= 5 else 'external_target_bias',
        proton_signal_background_ratio_author=signal_to_background,
        stated_run_duration_lo_h=12 if low else '', stated_run_duration_hi_h=24 if low else '',
        spot_area_typical_cm2=1, beam_current_typical_uA=150, D_flux_typical_per_cm2_s='1.1e15',
        Dplus_fraction_typical=.85, D2plus_fraction_typical=.15, retarding_spread_typical_percent=5,
        retarding_voltage_lo_kV=2.5, retarding_voltage_hi_kV=4.2,
        range_SRIM_reported_endpoint_nm={Decimal('.25'): 4.4, Decimal('6.5'): 67.5}.get(energy, ''),
        foil_thickness_um=250, cathode_diameter_cm=2.4, purity_stated_min_pct=99.95 if host=='Pd' else 98.5,
        anneal_temperature_C=400, anneal_duration_h=3, anneal_pressure_upper_Torr='1e-6',
        vacuum_side_renewal_nm=200, electrolyte_side_Pd_nm=250 if host=='Ti' else '',
        operating_pressure_lo_Torr='2.5e-4', operating_pressure_hi_Torr='3e-4', base_pressure_upper_Torr='1e-6',
        electrolyte_D2SO4_M=6 if loaded else '', galvanostatic_lo_mA_cm2=300 if loaded else '',
        galvanostatic_hi_mA_cm2=700 if loaded else '', alternative_constant_voltage_V=5 if loaded else '',
        host_atom_density_cm3='6.8e22' if host=='Pd' else '5.7e22', D_solubility_terminal_per_host=1 if host=='Pd' else 2,
        host_lattice_pristine='FCC' if host=='Pd' else 'HCP(alpha-Ti)',
        host_lattice_post_exposure='FCC' if host=='Pd' else 'HCP(alpha-Ti)+FCC(delta-TiD_x)',
        D_diffusivity_comparison='higher than Ti' if host=='Pd' else 'lower than Pd',
        electronic_conductivity_comparison='metallic, more itinerant carriers than Ti' if host=='Pd' else 'lower than Pd',
        D_diffusivity_m2_s='', electron_density_m3='', Debye_temperature_K='', conductivity_S_m='',
        fitted_bulk_Ue_keV=screen_bulk, fitted_bulk_Ue_sigma_keV=.2,
        fitted_surface_Ue_keV=screen_surface, fitted_crossover_E_cm_keV=.6,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('errbars', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    rows = [point(row) for row in read_points(args.errbars)]
    if len(rows) != 52 or len({r['row_id'] for r in rows}) != 52:
        raise ValueError('Expected the exact adopted 48 points plus 4 references')
    with args.output.open('x', newline='') as fd:
        writer = csv.DictWriter(fd, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
