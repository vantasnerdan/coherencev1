"""Assemble checked observation blocks; preserve source-specific columns and blanks."""
import argparse
import csv
from decimal import Decimal
from pathlib import Path

from compile_external import FIELDS

BLOCKS = ('D-D', 'D on Li (comparison)', 'p on solids (comparison)',
          'Other reactions (foil-state controls)', 'Material characterization')
KASAGI_LAB = {'yuki_kasagi1998_jetp68', 'kasagi2002_jpsj71', 'kasagi1998_iccf7'}
PER_D_LAB = KASAGI_LAB | {'raiola2006_bochum_thesis'}


def tag(rows):
    for row in rows:
        if row['observable_kind'] in {'target_deuteron_saturation_density', 'proton_Knight_shift'}:
            block = BLOCKS[4]
        elif row['reaction'] == 'D+D':
            block = BLOCKS[0]
        elif row['reaction'].startswith('D+'):
            block = BLOCKS[1]
        elif row['reaction'].startswith('p+'):
            block = BLOCKS[2]
        elif row['study_id'] == 'cvetinovic2023' or row['reaction'] == '3He+D':
            block = BLOCKS[3]
        else:
            raise ValueError(f'Unexpected reaction for {row["row_id"]}: {row["reaction"]}')
        row['measurement_block'] = block
        study = row['study_id']
        if block == BLOCKS[0] and study in PER_D_LAB and row['E_lab_keV'] and not row['E_cm_keV']:
            row['E_cm_keV'] = Decimal(row['E_lab_keV']) / 2
            row['energy_frame_source'] += '; CM=lab/2 (stationary target D; lab energy per deuteron)'
        if study in KASAGI_LAB and row['observable_kind']=='relative_thick_target_yield' and row['reference_group'].startswith(study+':'):
            row['reference_E_cm_keV'] = '5'
        if row.get('actual_D_per_metal'):
            row['reported_inferred_D_per_metal'] = row['actual_D_per_metal']
            row['actual_D_per_metal'] = ''
        if row.get('per_run_beam_current_uA') and study == 'raiola2002_plb547':
            row['reported_series_current_uA'] = row['per_run_beam_current_uA']
            row['per_run_beam_current_uA'] = ''
        if study == 'yuki_kasagi1998_jetp68' and row['observable_kind'] == 'thick_target_yield_enhancement':
            row['evidence_class'] = 'derived_ratio_shared_numerator'
        if study == 'karahadian2026':
            row['E_lab_keV'] = row['E_lab_nominal_Dplus_keV']
            row['source_locator'] = f'Fig. 2 {row["target"]} at {row["E_cm_keV"]} keV CM, with Methods and SI context'
        yield row


def merge(paths):
    fields = ['measurement_block', *FIELDS, 'reported_inferred_D_per_metal', 'reported_series_current_uA']
    rows = []
    row_ids = set()
    for path in paths:
        with path.open() as f:
            reader = csv.DictReader(f)
            for k in reader.fieldnames:
                if k not in fields:
                    fields.append(k)
            for row in tag(reader):
                row_id = row['row_id']
                if row_id in row_ids:
                    raise ValueError(f'Duplicate source row {row_id} in {path}')
                if not row['citation'] or not row['observable_kind'] or not row['observable_unit']:
                    raise ValueError(f'Missing source, units or observable definition for {row_id}')
                if row.get('evidence_class') == 'reported_fit' and (row.get('E_cm_keV') or row.get('E_lab_keV')):
                    raise ValueError(f'A reported fit may not masquerade as a point: {row_id}')
                row_ids.add(row_id)
                rows.append(row)
    # Grouping is not a data-label transform; original relative order is stable within each block.
    rows.sort(key=lambda r: BLOCKS.index(r['measurement_block']))
    return fields, rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output', type=Path)
    p.add_argument('input', type=Path, nargs='+')
    args = p.parse_args()
    fields, rows = merge(args.input)
    with args.output.open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
