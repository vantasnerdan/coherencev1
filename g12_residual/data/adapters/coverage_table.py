"""Count present point covariates by original study, reaction block, and evidence class."""
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

FIELDS = ('E_cm_keV', 'E_lab_keV', 'sigma_minus', 'temperature_K', 'actual_D_per_metal',
          'reported_inferred_D_per_metal', 'reported_series_current_uA',
          'per_run_beam_current_uA', 'current_density_uA_cm2', 'run_duration_h', 'dose_C',
          'accumulated_dose_C', 'ion_range_nm', 'operating_pressure_Torr', 'foil_thickness_um',
          'host_lattice', 'host_atom_density_cm3', 'electron_density_m3', 'Debye_temperature_K',
          'D_diffusivity_m2_s', 'D_solubility_terminal_per_host', 'conductivity_S_m',
          'cup_IV_current_uA', 'stated_run_duration_lo_h', 'range_SRIM_reported_endpoint_nm',
          'fitted_bulk_Ue_keV')
GROUP = ('measurement_block', 'study_id', 'observable_kind', 'evidence_class')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('observations', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    groups = defaultdict(Counter)
    with args.observations.open() as f:
        for r in csv.DictReader(f):
            group = groups[tuple(r[g] for g in GROUP)]
            group['rows'] += 1
            group.update(k for k in FIELDS if r.get(k, '').strip())
    with args.output.open('x', newline='') as f:
        w = csv.writer(f)
        w.writerow((*GROUP, 'rows', *(f'{k}_present' for k in FIELDS)))
        for group, c in sorted(groups.items()):
            w.writerow((*group, c['rows'], *(c[k] for k in FIELDS)))


if __name__ == '__main__':
    main()
