"""Cross modeled screening inventories with explicit literature scenarios.

Only a necessary capacity test is possible for UC Davis (occupied D, not site
capacity). Other crossings compare modeled site capacity to stated scenarios,
not to measured inventories of the experimental foils. Uses saved fit results;
no fitting and no new trajectory integration. Run python -m g7_micro.site_targets.
"""
import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).parent
MURASE = 'https://doi.org/10.2320/matertrans.MT-M2020220'
FUKAI = 'https://doi.org/10.1103/PhysRevLett.73.1640'
HEUSER = 'https://www.osti.gov/biblio/5902107'
AMIN = 'https://doi.org/10.1016/j.actamat.2016.03.054'
XU = 'https://doi.org/10.1088/1757-899X/668/1/012006'
ZHANG = 'https://doi.org/10.2320/matertrans.M2018312'
HUANG = 'https://doi.org/10.1038/s41586-026-11238-4'


def band_gb(grain_nm):
    # 1 nm *total* band; exactly 1 site per metal in a cubic-grain band scenario.
    return -math.expm1(3*math.log1p(-1./grain_nm))


def scenarios():
    # Each numeric range conditional on its supplied conversion assumptions.
    # None is evidence absent, not a zero inventory. One independent D per site.
    return [
        ('Pd', 'vacancy cluster', 'ordinary annealed or cold rolled', 'unknown', None, None, FUKAI,
         'No ordinary foil concentration or available D per vacancy established; do not transfer the 5 GPa result.'),
        ('Pd', 'vacancy cluster', '5 GPa PdH high-temperature scenario', 'capacity', .15/.85, .21/.79, FUKAI,
         '15–21 at.% after H removal; assumes vacancy fraction per lattice site, one special D site per vacancy. Not a cold-rolled or annealed concentration.'),
        ('Ti', 'vacancy cluster', 'ambient or loaded Ti', 'unknown', None, None, 'https://doi.org/10.3390/met12040595',
         'PAS differs after annealing and deformation; no available D capacity per Ti established.'),
        ('Pd', 'dislocation core', 'annealed Pd', 'unknown', None, None, 'https://doi.org/10.1088/0953-8984/3/14/018',
         'NMR distinguishes normal interstitial and likely trapped signal; not a matched annealed inventory.'),
        ('Pd', 'dislocation core', 'cold-worked dilute PdD SANS', 'occupied_floor', 5e15*1e10/6.8e28, 5e15*3e10/6.8e28, HEUSER,
         '1–3 D per Å at line density ~5e15 m−2; retained fuel in those specimens, not an upper bound on available special sites.'),
        ('Pd', 'dislocation core', 'Pd film, as-sputtered to after beta hydrogen cycling', 'capacity', 5.4e16*1e9/6.8e28, 1.45e17*1e9/6.8e28, AMIN,
         'Illustrative one special site per nm of line. HRTEM densities 5.4e16–1.45e17 m−2 across ex situ film endpoints; not ordinary rolled or annealed foils.'),
        ('Pd', 'dislocation core', 'rolled PdD NMR extended neighborhood', 'capacity', 0., math.pi*(3e-9)**2*1e16, 'https://doi.org/10.1088/0953-8984/3/14/018',
         'Upper scenario with one candidate site per metal in 3 nm radius, no correction for overlap; NMR-inferred dilated occupied neighborhood is not measured high-U reach.'),
        ('Ti', 'dislocation core', 'vacuum annealed then H-permeated CP Ti', 'capacity', 8.537e13/2.951e-10/5.67e28, 3.5647e14/2.951e-10/5.67e28, XU,
         'XRD Table 1 line densities; assume one special site per 0.2951 nm of line. Neither uncharged Ti nor loaded TiD2 validated.'),
        ('Ti', 'dislocation core', 'vacuum annealed then H-permeated CP Ti cylinder', 'capacity', math.pi*.5e-9**2*8.537e13, math.pi*.5e-9**2*3.5647e14, XU,
         'Alternative geometric scenario: radius .5 nm, one site per metal inside; do not add to line-site case.'),
        ('Ti', 'dislocation core', 'cold-rolled CP Ti', 'unknown', None, None, 'https://doi.org/10.1007/s11666-025-02031-2',
         'PAS detects dislocations after rolling; no comparable quantitative capacity extracted.'),
        ('Pd', 'grain-boundary sites', 'unpolished as-annealed Pd, >200 µm grains', 'capacity', 0., band_gb(200e3), MURASE,
         'Cubic grains, total 1 nm band, one available special site per in-band metal. Surface EBSD; do not imply a bulk measurement.'),
        ('Pd', 'grain-boundary sites', 'rolled then polished surface; annealed then polished surface', 'capacity', band_gb(200.), band_gb(190.), MURASE,
         'Average high-angle grain diameters .19/.20 µm in mechanically polished surface, not bulk grain sizes; same 1 nm band scenario.'),
        ('Ti', 'grain-boundary sites', '50% cryo-rolled then 420/450 °C annealed', 'capacity', band_gb(4010.), band_gb(1740.), ZHANG,
         'Same cubic grain/band/site convention; nominal TEM sizes in metallic Ti, not TiD2 or unrolled annealed foil.'),
        ('Pd', 'subsurface layer', 'modeled 1–50 nm with distinct adsorbate', 'slab', None, None, MURASE,
         'Observed polishing damage in other foils ~1 µm, not the high-U thickness. Use modeled occupied slab inventory, no universal capacity.'),
        ('Ti', 'subsurface layer', 'modeled 1–50 nm with distinct adsorbate', 'slab', None, None, 'g4_beam/merit/profiles.json',
         'Use modeled occupied slab inventory, no measured site count or universal depth.'),
        ('Pd', 'H-affinitive second phase', 'annealed or cold-worked Pd', 'unknown', None, None, HUANG,
         'Inconel phase-partitioning observation provides no target phase abundance, local capacity or screening.'),
        ('Ti', 'H-affinitive second phase', 'Ti/TiD', 'unknown', None, None, HUANG,
         'Need phase identity, abundance, local capacity per metal, actual D partition and matched depth; do not apply Inconel concentrations.'),
    ]


def crossings(curve, lo, hi, mode):
    """Return all contiguous runs of matching sampled energies, no extrapolation."""
    spans = []
    opened = False
    for u, band_lo, band_hi in curve:
        match = band_lo is not None and band_lo <= hi and (mode == 'cap_only' or lo <= band_hi)
        if match:
            if opened:
                spans[-1][1] = u
            else:
                spans.append([u, u])
        opened = match
    return ';'.join(f'{a:.6g}–{b:.6g}' for a, b in spans)


def main():
    fit = json.loads((HERE/'occupancy_results.json').read_text())
    with (HERE/'depth_inventory_curves.csv').open() as f:
        depth = list(csv.DictReader(f))
    with (HERE/'inventories.csv').open() as f:
        slab = list(csv.DictReader(f))
    with (HERE/'depth_region.csv').open() as f:
        cells = [r for r in csv.DictReader(f) if r['power'] == '1']
    volume = {(r['material'], r['d_nm']): (float(r['top_fuel_per_metal']), float(r['fuel_per_metal'])) for r in slab}
    ti = [(r['U_high_keV'], r['N_s_per_Ti_band_low'], r['N_s_per_Ti_band_high'])
          for r in fit['Ti']['conditional_fits'][0]['curve']]
    pd = fit['Pd']['response_pivots'][0]['curve']
    for row in depth:
        if float(row['profile_nominal_DD_cm_key_keV']) != 1. or int(row['power']) != 1:
            raise ValueError('Site crossings require the 1-keV CM, current-power-1 depth projection')

    def depth_curve(material, fitted, *, slab_only=False):
        """Lower bound on required volume occupancy, never credited to the sheet."""
        result = []
        for row in depth:
            if row['material'] != material:
                continue
            u = float(row['U_high_keV'])
            if fitted and not slab_only:
                value = row['full_band_volume_occupied_min_per_metal']
                result.append((u, float(value) if value else None, None))
                continue
            cost_limit = float(row['full_min_chi2_one_foil'])+float(row['descriptive_delta_chi2'])
            required = []
            for cell in cells:
                if cell['material'] != material or float(cell['U_high_keV']) != u:
                    continue
                if fitted and float(cell['full_chi2_one_foil']) > cost_limit:
                    continue
                f_b = max(0., float(cell['feasible_f_b_min']))
                if slab_only and f_b > 1e-12:
                    continue
                f_s = max(0., float(cell['feasible_f_s_min']))
                top, bulk = volume[material, cell['d_nm']]
                required.append(f_s*top+f_b*(bulk-top))
            total_key = ('full_band' if fitted else 'ratio_only')+'_occupied_min_per_metal'
            all_sheet = max(float(r['sheet_fuel_per_metal']) for r in slab if r['material'] == material)
            low = max(min(required), float(row[total_key])-all_sheet, 0.) if required and row[total_key] else None
            result.append((u, low, None))
        return sorted(result)

    rows = []
    for mat, cls, state, kind, lo, hi, url, note in scenarios():
        matd = mat+'D'
        state_lo, state_hi = lo, hi
        if kind == 'slab':
            # Max allowed occupied volume fuel, if every D in the slab is special.
            values = [float(r['top_fuel_per_metal']) for r in slab if r['material'] == matd]
            state_lo, state_hi = 0., max(values)
        datasets = [('UC_Davis_'+matd+'_full', depth_curve(matd, True, slab_only=kind=='slab'), 'cap_only'),
                    ('UC_Davis_'+matd+'_ratio', depth_curve(matd, False, slab_only=kind=='slab'), 'cap_only')]
        if mat == 'Ti':
            datasets += [('Raiola_Ti_5keV_CM', ti, 'capacity')]
        else:
            datasets += [('Cvetinovic_FD_500keV_'+key, [(r['U_high_keV'], *r[key+'_band']) for r in pd], 'capacity')
                         for key in ('hard','soft','delta')]
        for name, curve, mode in datasets:
            allowed = kind == 'capacity' or (kind == 'slab' and mode == 'cap_only')
            if name.endswith('_delta'):
                # A generic absolute range cannot constrain a hard-minus-soft change.
                allowed = False
            target = crossings(curve, state_lo, state_hi, mode) if allowed else ''
            status = ('necessary_only' if mode == 'cap_only' else 'conditional_overlap') if target else (
                     'no_sampled_overlap' if allowed else 'not_identified')
            rows.append(dict(material=mat, site_class=cls, condition=state, inventory_kind=kind,
                             inventory_min_per_metal=state_lo, inventory_max_per_metal=state_hi,
                             dataset=name, result=status, sampled_U_high_keV=target,
                             source=url, assumptions_and_scope=note))
    # A matched geometric scenario can constrain an inventory difference.
    gdelta = band_gb(190.)-band_gb(200.)
    target = crossings([(r['U_high_keV'], *r['delta_band']) for r in pd], gdelta, gdelta, 'capacity')
    rows.append(dict(material='Pd', site_class='grain-boundary sites',
                     condition='mechanically polished rolled minus annealed Pd (geometric alternative)',
                     inventory_kind='capacity_difference', inventory_min_per_metal=gdelta,
                     inventory_max_per_metal=gdelta, dataset='Cvetinovic_FD_500keV_delta',
                     result='conditional_overlap' if target else 'no_sampled_overlap',
                     sampled_U_high_keV=target, source=MURASE,
                     assumptions_and_scope='Uses .19 and .20 µm as cubic diameters and a 1 nm band at one site per metal, not measured capacity change in Cvetinović foils.'))
    out = HERE/'site_targets.csv'
    with out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print(out)


if __name__ == '__main__':
    main()
