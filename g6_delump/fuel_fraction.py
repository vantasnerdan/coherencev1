"""Fuel fraction against screening: D beam into a metal holding tritium at T/M = x.

Diluting the fuel in a metal costs drag per fuel atom (the metal's stopping is
shared by fewer T nuclei) but measured screening rises as loading falls. This
script puts both on one footing:

  1. local break-even (sigma*Q = drag per T nucleus) against x at E_lab = 1, 2,
     5, 13 keV: the cross section needed, the screening that delivers it, and
     the ratio of the need to the screening ceiling. In the Eq 7 form the D-T
     cross section at fixed E peaks at U_e = BG^2 - E_cm (about 1.18 MeV), where
     the whole barrier is removed; no screening exceeds that ceiling;
  2. the thick-target energy return G(E0) on a grid of x and U_e;
  3. G for the measured loading series in Ti (Raiola et al. 2005, EXFOR O1379),
     with the measured U_e carried from D-D to D-T unchanged (both Z1*Z2 = 1),
     against TiT1.7, and the beam energy where dilution stops paying.

The ICRU velocity law is used below 10 keV/u (20 keV lab D), so beam energies
here stop at 13 keV.

Inputs are the same footing as g6_gap and g4_beam: ICRU-49 velocity-law
electronic stopping and ZBL nuclear stopping from g4_beam/forward/model.py,
Bragg additivity, and the Bosch-Hale D-T cross section in the Eq 7 screened
form from g6_gap/substrate_framework/lithium_fuel_return.py.

Run from the repository root:  python -m g6_delump.fuel_fraction
NumPy only.
"""
import csv
import json
from pathlib import Path

import numpy as np

from g4_beam.forward.model import ELEM, Se_atom, Sn_atom
from g6_gap.substrate_framework.lithium_fuel_return import (
    DT_BG, DT_Q_KEV, MASS_D, MASS_T, break_even_barns, dt_eq7_sigma_mb, thick_target_return)

OUT = Path(__file__).resolve().parent
E_FLOOR = 0.5*(MASS_D+MASS_T)/MASS_T      # lab keV at E_cm = 0.5 keV, Bosch-Hale fit floor
E_GRID = np.geomspace(E_FLOOR, 200., 4000)
X_GRID = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
UE_GRID = [0.04, 0.1, 0.3, 1.0, 2.0, 5.0]
E_LOCAL = [1., 2., 5., 13.]
E0_RETURN = [1., 2., 3., 5., 8., 13.]

# Measured D-D screening in Ti against D/Ti (Raiola et al., J. Phys. G 31, 1141 (2005); EXFOR O1379).
# The two highest-loading points are upper limits; they are used at the limit.
TI_SERIES = [
    dict(x=2.1, Ue_keV=0.030, note='upper limit, -10 C'),
    dict(x=1.1, Ue_keV=0.050, note='upper limit, 50 C'),
    dict(x=0.26, Ue_keV=0.250, note='250+-40 eV, 100 C'),
    dict(x=0.20, Ue_keV=0.290, note='290+-65 eV, 200 C'),
]


def drag_per_fuel(energy_lab, host, x):
    """Total stopping per T nucleus, eV/(1e15 T/cm^2), host M with T/M = x."""
    h = ELEM[host]
    per_host = (Se_atom(energy_lab, host)+Sn_atom(energy_lab, h['Z'], h['M'])
                + x*(Se_atom(energy_lab, 'H')+Sn_atom(energy_lab, 1, MASS_T)))
    return per_host/x


def required_ue(energy_lab, sigma_barns, hi=2000.):
    """Screening (keV) that lifts the Eq 7 D-T cross section to sigma_barns; bisection in log."""
    target = sigma_barns*1e3                        # mb
    lo_ue, hi_ue = 1e-6, hi
    if dt_eq7_sigma_mb(energy_lab, hi_ue) < target:
        return float('nan')
    if dt_eq7_sigma_mb(energy_lab, lo_ue) >= target:
        return 0.
    for _ in range(200):
        mid = np.sqrt(lo_ue*hi_ue)
        if dt_eq7_sigma_mb(energy_lab, mid) >= target:
            hi_ue = mid
        else:
            lo_ue = mid
    return float(hi_ue)


def sigma_ceiling_barns(energy_lab):
    """Largest Eq 7 D-T cross section at this energy over all U_e (at U_e = BG^2 - E_cm), barns."""
    e_cm = energy_lab*MASS_T/(MASS_D+MASS_T)
    return float(dt_eq7_sigma_mb(energy_lab, DT_BG**2-e_cm))/1e3


def g_return(host, x, ue):
    """Thick-target D-T energy return G at each E0 in E0_RETURN."""
    g = thick_target_return(E_GRID, drag_per_fuel(E_GRID, host, x), ue)
    return {e0: float(np.interp(e0, E_GRID, g)) for e0 in E0_RETURN}


def main():
    check = float(break_even_barns(drag_per_fuel(np.array([1.]), 'Ti', 1.7))[0])
    assert abs(check-206.)/206. < .01, f'TiT1.7 break-even at 1 keV is {check:.1f} b, expected 206 b'

    local = []
    for host in ('Pd', 'Ti'):
        for x in X_GRID:
            for e in E_LOCAL:
                barns = float(break_even_barns(drag_per_fuel(np.array([e]), host, x))[0])
                ceiling = sigma_ceiling_barns(e)
                local.append(dict(host=host, x=x, E_lab_keV=e, break_even_b=round(barns, 2),
                                  ceiling_b=round(ceiling, 2), drag_cut_needed=round(barns/ceiling, 3),
                                  Ue_needed_keV=round(required_ue(e, barns), 2)))

    grid = []
    for host in ('Pd', 'Ti'):
        for x in X_GRID:
            for ue in UE_GRID:
                g = g_return(host, x, ue)
                grid.append(dict(host=host, x=x, Ue_keV=ue, **{f'G_E0_{int(e0)}keV': g[e0] for e0 in E0_RETURN}))

    ti = []
    reference = g_return('Ti', 1.7, 0.030)
    for row in TI_SERIES:
        g = g_return('Ti', row['x'], row['Ue_keV'])
        ti.append(dict(**row, **{f'G_E0_{int(e0)}keV': g[e0] for e0 in E0_RETURN},
                       **{f'ratio_to_TiT1.7_E0_{int(e0)}keV': g[e0]/reference[e0] for e0 in E0_RETURN}))

    fine = np.geomspace(E_FLOOR*1.01, 13., 400)
    crossover = {}
    for row in TI_SERIES[2:]:
        dilute = thick_target_return(E_GRID, drag_per_fuel(E_GRID, 'Ti', row['x']), row['Ue_keV'])
        full = thick_target_return(E_GRID, drag_per_fuel(E_GRID, 'Ti', 1.7), .030)
        ratio = np.interp(fine, E_GRID, dilute)/np.interp(fine, E_GRID, full)
        above = np.nonzero(ratio < 1)[0]
        crossover[f"x{row['x']}"] = float(fine[above[0]]) if above.size else None

    for name, rows in (('break_even_local.csv', local), ('return_grid.csv', grid), ('ti_series.csv', ti)):
        with open(OUT/name, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    summary = dict(check_TiT1p7_break_even_1keV_b=round(check, 2), Q_keV=DT_Q_KEV, E_floor_lab_keV=E_FLOOR,
                   sigma_ceiling_b={str(e): round(sigma_ceiling_barns(e), 2) for e in E_LOCAL},
                   reference_TiT1p7_Ue30eV=reference,
                   dilution_pays_below_E0_keV=crossover)
    (OUT/'fuel_fraction.json').write_text(json.dumps(summary, indent=1)+'\n')


if __name__ == '__main__':
    main()
