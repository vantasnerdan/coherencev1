"""Transcribe named public CSV/EXFOR tables without fitting or inventing settings."""
import argparse
import csv
import re
from decimal import Decimal
from pathlib import Path

FIELDS = 'row_id,study_id,lab,run_id,reaction,channel,target,host,target_preparation,E_cm_keV,E_lab_keV,energy_frame_source,observable_kind,observable,observable_unit,sigma_minus,sigma_plus,reference_E_cm_keV,reference_group,target_geometry,evidence_class,temperature_K,actual_D_per_metal,per_run_beam_current_uA,current_density_uA_cm2,run_duration_h,dose_C,accumulated_dose_C,ion_range_nm,operating_pressure_Torr,foil_thickness_um,host_lattice,host_atom_density_cm3,electron_density_m3,Debye_temperature_K,D_diffusivity_m2_s,D_solubility_terminal_per_host,conductivity_S_m,fitted_bulk_Ue_keV,citation,source_locator,conditions_note'.split(',')
X4_BASE = 'https://github.com/IAEA-NDS/exfor_master/blob/2193cdb7066421fb7931b696d6bab219f9cc8b03/exforall/'
EXFOR_FOLDERS = {'O1379': 'O13', 'O1282': 'O12', 'A0632': 'A06', 'O1084': 'O10'}
CZERSKI = 'https://arxiv.org/abs/2409.02112v3'


def x4_url(entry):
    return f'{X4_BASE}{EXFOR_FOLDERS[entry]}/{entry}.x4'


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def exfor_blocks(path):
    """Return subentry BIB/body and fixed-column DATA fields, preserving blanks."""
    blocks = {}
    for raw in re.split(r'^SUBENT +', path.read_text().replace('\r', ''), flags=re.M)[1:]:
        subid = raw.split()[0]
        lines = raw.splitlines()
        header = next((i for i, line in enumerate(lines) if line.startswith('DATA ')), None)
        if header is None:
            blocks[subid] = (raw, [])
            continue
        _, width, expected, *_ = lines[header].split()
        stride = 11
        n = int(width)
        rows = [tuple(line[stride*i:stride*(i+1)].strip() for i in range(n))
                for line in lines[header+3:header+3+int(expected)]]
        assert len(rows) == int(expected) and all(row[0] for row in rows), subid
        assert lines[header+3+int(expected)].startswith('ENDDATA'), subid
        blocks[subid] = ('\n'.join(lines[:header]), rows)
    return blocks


def temp_series(src):
    doi = 'https://doi.org/10.1088/0954-3899/31/11/002'
    blocks = exfor_blocks(src/'O1379.x4')
    outline = dict(study_id='raiola2005', lab='Bochum', reaction='D+D', channel='p+t',
                   target_preparation='in situ Kr sputter-cleaned', citation=f'{doi}; {x4_url("O1379")}')
    for last, host in [('2', 'Pt'), ('3', 'Hf')]:
        sid = 'O137900'+last
        _, data = blocks[sid]
        for i, (e, deg, s, err, sol) in enumerate(data, 1):
            note = f'Beam current <2 µA; EXFOR compiled E uncertainty 0.03%, energy digitization {"1.2" if host=="Pt" else "0.7"}%. '
            if sol:
                note += f'Attached EXFOR MISC solubility D/metal={sol}, not a point-by-point assay.'
            yield dict(outline, row_id=f'exfor:{sid}:{i:03}', target=host, host=host,
                       E_cm_keV=e, energy_frame_source='EXFOR EN-CM', temperature_K=Decimal(deg)+Decimal('273.15'),
                       observable_kind='S_factor', observable=s, observable_unit='keV*b', sigma_minus=err, sigma_plus=err,
                       reference_group=f'{sid}:absolute-S:{deg}C', evidence_class='figure_digitized_measurement',
                       source_locator=f'EXFOR {sid} DATA row {i}', conditions_note=note.strip())
    bib, _ = blocks['O1379001']
    host = ''
    summary = bib.split('material   Temp')[1].split('-by authors.')[0]
    row_pattern = re.compile(r'^\s*(?:(\w+)\s+)?(-?\d+)\s+(<\d+|\d+\+/-\d+)\s+(\d+\.\d+)(?:\s+(\S+))?\s*$')
    for line in summary.splitlines():
        match = row_pattern.match(line)
        if not match:
            continue
        next_host, deg, value, sol, electron = match.groups()
        host = next_host or host
        estimate, _, err = value.removeprefix('<').partition('+/-')
        upper_limit = value.startswith('<')
        note = f'Upper limit, not a detection. ' if upper_limit else 'Reported U fit at this temperature, not a measured energy point. '
        note += f'Authors list solubility D/metal={sol}; beam current <2 µA.'
        if electron:
            note += f' Authors report effective electron count {electron} per atom (not electron density).'
        yield dict(outline, row_id=f'raiola2005:O1379001:Ue:{host}:{deg}C', target=host, host=host,
                   observable_kind='fitted_screening_energy_upper_limit' if upper_limit else 'fitted_screening_energy',
                   observable=estimate, observable_unit='eV', sigma_minus=err, sigma_plus=err,
                   evidence_class='reported_fit', temperature_K=Decimal(deg)+Decimal('273.15'),
                   fitted_bulk_Ue_keV='' if upper_limit else Decimal(estimate)/1000,
                   source_locator='EXFOR O1379001 author COMMENT summary table', conditions_note=note)


def proton_lithium(src):
    blocks = exfor_blocks(src/'O1282.x4')
    targets = (('natural Li', 'Li', 'Natural Li metal (92.58% Li-7)', 'thick', ''),
               ('Li2WO4/steel', 'Li2WO4', 'Evaporated on steel; 360 µg/cm²', 'thin film', ''),
               ('PdLi', 'Pd', 'Plasma-produced alloy; vacuum anneal 850°C for 1 h; isotope/content per subentry', 'thick', 200))
    citation = f'https://doi.org/10.1016/j.physletb.2005.08.036; {x4_url("O1282")}'
    for j in range(2, 10):
        sid = f'O128200{j}'
        bib, rows = blocks[sid]
        isotope = 7 if j in (2, 3, 4, 8) else 6
        branch = 'alpha+alpha' if isotope == 7 else 'alpha+3He'
        # Experimental EN is the incident projectile frame (EXFOR EN, not EN-CM).
        refenergy = 83.3 if isotope == 7 else 81.6
        baseline = '54±4 µb at 83.3 keV lab (O1282008)' if isotope == 7 else '2.2±0.2 mb at 81.6 keV lab (O1282009)'
        part = dict(study_id='cruz2005', lab='Bochum/LUNA', reaction=f'p+{isotope}Li', channel=branch, citation=citation)
        for i, row in enumerate(rows, 1):
            e, value, err, *flag = row
            target, host, prep, geometry, thickness = targets[(j-2)%3 if not flag else int(Decimal(flag[0]))-1]
            cross = bool(flag)
            yield dict(part, row_id=f'exfor:{sid}:{i:03}', target=target, host=host, target_preparation=prep,
                       E_lab_keV=e, energy_frame_source='EXFOR EN (lab)', observable=value, sigma_minus=err, sigma_plus=err,
                       observable_kind='cross_section_reference_check' if cross else 'S_factor',
                       observable_unit=('microbarn' if isotope == 7 else 'mb') if cross else ('keV*b' if isotope == 7 else 'MeV*b'),
                       target_geometry=geometry, foil_thickness_um=thickness,
                       reference_group=f'cruz2005:p+{isotope}Li:{refenergy}-keV-lab-calibration',
                       evidence_class='table_measurement', source_locator=f'EXFOR {sid} DATA row {i}',
                       conditions_note=(f'Beam current <20 µA. S-factor linked to independent standard: {baseline}; '
                                        'keep reference dependence. X4 records thick-target differentiation for thick specimens. '
                                        'Cross-section check rows are not S-factor observations.'))
        if j <= 7:
            m = re.search(r'found as\s+(\d+)\+/-(\d+) eV', bib)
            u, err = m.groups()
            target, host, prep, geometry, thickness = targets[(j-2)%3]
            yield dict(part, row_id=f'exfor:{sid}:Ue', target=target, host=host, target_preparation=prep,
                       observable_kind='fitted_screening_energy', observable=u, observable_unit='eV', sigma_minus=err,
                       sigma_plus=err, evidence_class='reported_fit', fitted_bulk_Ue_keV=Decimal(u)/1000,
                       source_locator=f'EXFOR {sid} ADD-RES', target_geometry=geometry, foil_thickness_um=thickness,
                       conditions_note=f'Author fit to that subentry, not an energy-resolved point. Standard: {baseline}.')


def deuteron_lithium(src):
    sid = 'A0632003'
    _, rows = exfor_blocks(src/'A0632.x4')[sid]
    for i, (e, s, err, flag) in enumerate(rows, 1):
        if Decimal(flag) != 3:
            continue  # inverse kinematics target is D2 gas
        yield dict(row_id=f'exfor:{sid}:{i:03}', study_id='engstler1992', lab='Bochum/Münster', reaction='D+6Li',
                   channel='alpha+4He', target='LiF solid', host='LiF', E_cm_keV=e, energy_frame_source='EXFOR EN-CM',
                   observable_kind='S_factor', observable=s, observable_unit='MeV*b', sigma_minus=err, sigma_plus=err,
                   reference_group='engstler1992:A0632003:flag3', evidence_class='table_measurement',
                   citation=f'Engstler et al., Z. Phys. A 342, 471 (1992); {x4_url("A0632")}',
                   source_locator=f'EXFOR {sid} Table 1 DATA row {i}, FLAG=3',
                   conditions_note='D beam on solid LiF; no inverse-kinematics gas points. 10% common normalization uncertainty reported separately from quoted point errors.')


def czerski(src):
    points = read_csv(src/'czerski-figure_points.csv')
    memberships = read_csv(src/'czerski-membership.csv')
    assert len(points) == len(memberships)
    tables = {r['target']: r for r in read_csv(src/'czerski-table2.csv')}
    targets = {'Pd':'Pd plate (implanted, source D/Pd≈0.05)', 'Ti':'Ti0.993Pd0.007 plate (implanted, source D/Ti≈2)',
               'Zr':'Zr plate (implanted, source D/Zr≈2)'}
    regimes = {'Pd':'Constant power approx. 0.6 W, ~30–32°C, no deceleration; D-density fluctuation ~15%.',
               'Ti':'Air-cooled implanted Ti alloy under protocol similar to Zr; no point temperature stated.',
               'Zr':'Target kept below 45°C; air cooled and current reduced at higher energy; D-density fluctuations ~10%.'}
    head = dict(study_id='czerski2409v3', lab='Szczecin', reaction='D+D', channel='p+t', target_geometry='thick',
                foil_thickness_um=1000, citation=CZERSKI)
    for r, m in zip(points, memberships):
        assert (r['figure_panel'], r['target'], r['beam_species']) == (m['figure_panel'], m['target'], m['species'])
        if m['used_in_yield_objective'] != 'True':
            continue
        host = r['target']
        yield dict(head, row_id=f'czerski2409v3:{m["source_csv_line"]}', target=host, host=host,
                   target_preparation=targets[host], E_cm_keV=r['energy_cm_keV'], E_lab_keV=r['energy_lab_keV_per_deuteron'],
                   energy_frame_source='Source plotted lab energy per D; stationary equal-mass D, CM=lab/2',
                   observable_kind='detected_thick_target_yield', observable=r['value'], observable_unit=r['unit'],
                   evidence_class='figure_digitized_measurement', source_locator=f'Fig. {r["figure_panel"]} PDF ({r["x_pdf_pt"]},{r["y_pdf_pt"]}) pt',
                   conditions_note=f'{r["beam_species"]}. Counts are detector-normalized, not total emitted reactions. Published bars not digitized. {regimes[host]}')
    for host, r in tables.items():
        e = r['screening_energy_cold_eV']
        s = r['screening_energy_cold_uncertainty_eV']
        yield dict(head, row_id=f'czerski2409v3:table2:{host}:cold-Ue', target=host, host=host,
                   target_preparation=targets[host], observable_kind='fitted_screening_energy', observable=e,
                   observable_unit='eV', sigma_minus=s, sigma_plus=s, fitted_bulk_Ue_keV=Decimal(e)/1000,
                   evidence_class='reported_fit', source_locator='Table 2 (cold screening)',
                   conditions_note=f'{regimes[host]} Non-independent of the source yield set; inferred hot-track Ue excluded.')


def yuki(src):
    records = exfor_blocks(src/'O1084.x4')
    for sid, host, prep, fit, sig, sys, thickness in [
        ('O1084002','Yb','Beam implanted at 6.45 keV lab until near steady proton yield',81,10,
         'Stopping systematic +25/−35 eV; D inventory +13 eV',1000),
        ('O1084003','Ti','Gas-preloaded Ti, with subsequent beam implantation',19,12,
         'Stopping systematic ±12 eV; D inventory +4 eV',2000),
    ]:
        base = dict(study_id='yuki1997', lab='Tohoku/RIKEN', reaction='D+D', channel='p+t', target=host, host=host,
                    target_preparation=prep, target_geometry='thick', foil_thickness_um=thickness,
                    citation=f'https://doi.org/10.1088/0954-3899/23/10/036; {x4_url("O1084")}')
        for i, (e, value, percent) in enumerate(records[sid][1], 1):
            error = Decimal(value)*Decimal(percent)/100 if percent else ''
            yield dict(base, row_id=f'exfor:{sid}:{i:03}', E_cm_keV=Decimal(e)/2, E_lab_keV=e,
                       energy_frame_source='EXFOR lab; CM=lab/2', observable_kind='relative_thick_target_yield',
                       observable=value, observable_unit='1', sigma_minus=error, sigma_plus=error,
                       reference_E_cm_keV='3.225', reference_group=f'yuki1997:{host}:6.45-keV-lab',
                       evidence_class='figure_digitized_measurement', source_locator=f'EXFOR {sid} DATA row {i}',
                       conditions_note='Common digitization resolution 8% yield and 0.018 keV lab. Methods quote approximate target currents 100 µA at 2.55 and 500 µA at 6.45 keV lab; not assigned to nearby digitized energies. Point-specific charge not tabulated.')
        yield dict(base, row_id=f'exfor:{sid}:Ue', observable_kind='fitted_screening_energy', observable=fit,
                   observable_unit='eV', sigma_minus=sig, sigma_plus=sig, evidence_class='reported_fit',
                   fitted_bulk_Ue_keV=Decimal(fit)/1000, source_locator=f'EXFOR {sid} ADD-RES and paper §3',
                   conditions_note=f'Fitted to the same experimental yield series; errors in sigma fields are statistical. {sys}.')


def chen():
    common = dict(study_id='chen2025', lab='UBC/TRIUMF', reaction='D+D', channel='n+3He', host='Pd', target_geometry='thick',
                  foil_thickness_um=300, E_cm_keV=15, E_lab_keV=30, energy_frame_source='Nominal D+ sheath lab energy, CM=lab/2',
                  citation='https://doi.org/10.1038/s41586-025-09042-7', target_preparation='Annealed Pd after earlier beam-loading control (400°C for 1 h under vacuum)')
    for i, (a, ae, b, be, pct, pe) in enumerate([('135.5','.5','156.7','.4','15.1','.1'),
                                                  ('142.9','.5','159.2','.5','11.3','.1'),
                                                  ('138.6','.3','164.3','1.3','18.5','.2')], 1):
        base = dict(common, run_id=f'Fig 3 campaign {i}')
        for phase, value, e, condition in [('I', a, ae, 'beam only; source reports 1–2e-5 Torr during beam loading'),
                                           ('II', b, be, 'electrochemical loading; source reports final pressure usually 3.5–4e-5 Torr')]:
            yield dict(base, row_id=f'chen2025:fig3:{i}:{phase}', target=f'Pd phase {phase}',
                       observable_kind='neutron_rate_phase_saturation', observable=value, observable_unit='neutrons/s',
                       sigma_minus=e, sigma_plus=e, evidence_class='table_measurement', source_locator='Main text, electrochemically enhanced fusion',
                       conditions_note=f'{condition}. 0.5 mA reported setpoint, not actual point current; no direct in-run D/Pd. The ~0.7 XRD estimate was a separate test.')
        yield dict(base, row_id=f'chen2025:fig3:{i}:contrast', target='Pd phase II vs I',
                   observable_kind='neutron_rate_electrochemical_percent_increase', observable=pct, observable_unit='percent',
                   sigma_minus=pe, sigma_plus=pe, reference_E_cm_keV=15, reference_group=f'chen2025:fig3:{i}:I',
                   evidence_class='measured_intervention_ratio', source_locator='Fig. 3 caption, campaign contrasts',
                   conditions_note='Same phase rates underlie this contrast; do not treat all three rows as independent. At 15 keV CM, outside Davis range.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sources', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    rows = [row for adapter in (temp_series, proton_lithium, deuteron_lithium, czerski, yuki)
            for row in adapter(args.sources)]
    rows.extend(chen())
    if len({r['row_id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate source row ID')
    with args.output.open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
