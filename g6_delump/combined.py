"""Combined energy return: D beam into a tritium target, every pass 6 input on one footing.

For each target the script takes the drag per T nucleus, then gives the thick-target
D-T energy return G (fusion energy out over beam energy in) at beam energies 1, 2, 5
and 13 keV lab, and the two numbers a build has to hit:

  * at a given screening, the factor by which drag per T nucleus must fall for G = 1
    (G scales as 1/drag when the drag curve is scaled uniformly, so this is 1/G);
  * at the target's own drag, the screening U_e that gives G = 1, or none if even
    removing the whole barrier (the Eq 7 maximum) is not enough.

Targets (drag per T nucleus, eV per 1e15 T/cm^2):
  TiT1.7, PdT, Ti at T/Ti = 0.2   ICRU-49 velocity-law electronic + ZBL nuclear, Bragg sum,
                                  as g6_delump/fuel_fraction.py;
  LiT                             computed electronic drag of g6_gap (q4, LiD response at the
                                  same volume) + ZBL on Li6 and T, as g6_gap/tables;
  pure tritium                    ICRU-49 H electronic + ZBL on T;
  T as the only scatterer         ZBL on T alone: the floor for a host with no electronic drag.

Screening levels carried from D-D to D-T unchanged: 0.04 keV (computed static lowering,
g6_collision and g7_defect), 0.3 keV (measured D-D in metals at low loading), 1.4 keV
(the uniform UC Davis PdD fit). Channelling multiplies the PdT return by the g6_channel
gain at its computed point only (Pd(D/T) <110>, 13 keV, 2 deg: 1.777), a gain computed at
U_e = 0 and carried to the screened rows unchanged.

The slowing-down integral runs from E0 down to 0.1 keV lab, the g6_gap convention. Below
0.834 keV lab (0.5 keV centre of mass) the Bosch-Hale fit is an extrapolation. This floor
matters only where the screening needed is tens of keV or more; the U_e for G = 1 is
also given at floors 0.2 and 0.834 keV.

The screened cross section is the Eq 7 form used throughout g4_beam and g6_gap,
S(E)/sqrt(E(E+U)) exp(-B_G/sqrt(E+U)). The standard Assenbaum-Langanke-Rolfs
shift, S(E+U)/E exp(-B_G/sqrt(E+U)) ("form A"), keeps the true-energy flux factor
and evaluates S at E+U. Where U is much larger than E the two differ by several
times, so the U for G = 1 is also given in form A. Form A is capped at the s-wave
unitarity limit (604 b at 1 keV lab, 46.5 b at 13 keV); applying the same cap to
Eq 7 changes none of its U values.

Run from the repository root:  python -m g6_delump.combined
NumPy only.
"""
import csv
import json
from pathlib import Path

import numpy as np

from g4_beam.forward.model import Se_atom, Sn_atom
from g6_delump.fuel_fraction import E_FLOOR as FIT_FLOOR, drag_per_fuel
from g6_gap.substrate_framework.lithium_fuel_return import (
    DT_BG, DT_Q_KEV, MASS_D, MASS_T, break_even_barns, dt_eq7_sigma_mb, lit_zbl_nuclear, thick_target_return)

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
E_FLOOR = .1                                # lab keV, g6_gap inherited cutoff
SENSITIVITY_FLOORS = [.2, FIT_FLOOR]
E_GRID = np.geomspace(E_FLOOR, 13., 3000)
E0 = [1., 2., 5., 13.]
UE_LEVELS = {'computed_static': .04, 'measured_low_loading': .3, 'fitted_UC_Davis': 1.4}
CHANNEL_GAIN = dict(target='PdT', E0_keV=13., gain=1.7767, gain_se=.0728,
                    source='g6_channel density-dft-fit, Pd(D/T) <110>, 2 deg RMS, 300 K, U_e = 0')


def lit_computed_drag(energy):
    """LiT drag per formula (one T): g6_gap q4 electronic + ZBL on Li6 and T."""
    rows = list(csv.DictReader(open(ROOT/'g6_gap'/'tables'/'stopping-band.csv')))
    e = np.array([float(r['E_D_lab_keV']) for r in rows])
    se = np.array([float(r['Se_q4_eV_A']) for r in rows])
    sn = np.array([float(r['Sn_ZBL_D_on_Li6T_eV_A']) for r in rows])
    volume = 10*lit_zbl_nuclear(e)/sn                       # A^3 per formula, constant in the table
    assert np.ptp(volume)/volume.mean() < 1e-9
    return np.interp(energy, e, se)*volume[0]/10+lit_zbl_nuclear(energy)


def targets(grid=E_GRID):
    return {
        'TiT1.7': drag_per_fuel(grid, 'Ti', 1.7),
        'PdT': drag_per_fuel(grid, 'Pd', 1.),
        'Ti, T/Ti 0.2': drag_per_fuel(grid, 'Ti', .2),
        'LiT, computed drag': lit_computed_drag(grid),
        'pure tritium': Se_atom(grid, 'H')+Sn_atom(grid, 1, MASS_T),
        'T as only scatterer': Sn_atom(grid, 1, MASS_T),
    }


def dt_form_a_sigma_mb(energy_lab_keV, screening_keV=0.):
    """Assenbaum-Langanke-Rolfs shifted form, S(E+U)/E exp(-B_G/sqrt(E+U)), mb."""
    e = np.asarray(energy_lab_keV, dtype=float)*MASS_T/(MASS_D+MASS_T)
    x = e+screening_keV
    s = (6.927e4+x*(7.454e8+x*(2.050e6+x*5.2002e4)))/(1+x*(63.8+x*(-.995+x*(6.981e-5+x*1.728e-4))))
    return s/e*np.exp(-DT_BG/np.sqrt(x))


def dt_unitarity_mb(energy_lab_keV):
    """s-wave J = 3/2 D-T limit g pi lambda-bar^2, g = 4/6, in mb."""
    e_cm_mev = np.asarray(energy_lab_keV, dtype=float)*MASS_T/(MASS_D+MASS_T)*1e-3
    mu_mev = MASS_D*MASS_T/(MASS_D+MASS_T)*931.494
    return (2/3)*np.pi*197.327**2/(2*mu_mev*e_cm_mev)*10.          # fm^2 -> mb


def thick_target_return_form_a(energy, drag, screening_keV):
    """Same integration as g6_gap thick_target_return, form A capped at the unitarity limit."""
    integrand = np.minimum(dt_form_a_sigma_mb(energy, screening_keV), dt_unitarity_mb(energy))/drag
    probability = np.concatenate([[0.], np.cumsum(np.diff(energy)*(integrand[1:]+integrand[:-1])/2.)])*1e-9
    return DT_Q_KEV/energy*probability


def g_at(drag, ue):
    g = thick_target_return(E_GRID, drag, ue)
    return {e0: float(np.interp(e0, E_GRID, g)) for e0 in E0}


def ue_for_unity(drag, e0, gain=1., grid=E_GRID, ret=thick_target_return):
    """Smallest U_e (keV) with gain*G(e0) = 1, or None if no U_e reaches it."""
    def g(ue):
        return gain*float(np.interp(e0, grid, ret(grid, drag, ue)))
    scan = np.geomspace(1e-3, DT_BG**2, 400)
    values = np.array([g(u) for u in scan])
    best = int(np.argmax(values))
    if values[best] < 1:
        return None, float(values[best]), float(scan[best])
    first = int(np.nonzero(values >= 1)[0][0])
    lo, hi = (scan[first-1], scan[first]) if first else (scan[0], scan[0])
    for _ in range(100):
        mid = np.sqrt(lo*hi)
        lo, hi = (lo, mid) if g(mid) >= 1 else (mid, hi)
    return float(hi), float(values[best]), float(scan[best])


def main():
    drags = targets()
    check = float(np.interp(1., E_GRID, break_even_barns(drags['LiT, computed drag'])))
    assert abs(check-52.33)/52.33 < .005, f'LiT break-even at 1 keV is {check:.2f} b, expected 52.33 b'
    same = thick_target_return(E_GRID, drags['PdT'], 1.4)
    assert np.allclose(same, DT_Q_KEV/E_GRID*np.concatenate([[0.], np.cumsum(np.diff(E_GRID)*(
        (dt_eq7_sigma_mb(E_GRID, 1.4)/drags['PdT'])[1:]+(dt_eq7_sigma_mb(E_GRID, 1.4)/drags['PdT'])[:-1])/2.)])*1e-9, rtol=1e-12)
    check_ti = float(np.interp(1., E_GRID, break_even_barns(drags['TiT1.7'])))
    assert abs(check_ti-206.)/206. < .01, f'TiT1.7 break-even at 1 keV is {check_ti:.1f} b, expected 206 b'

    floors = {f: (np.geomspace(f, 13., 3000)) for f in SENSITIVITY_FLOORS}
    floor_drags = {f: targets(grid) for f, grid in floors.items()}

    rows = []
    for name, drag in drags.items():
        returns = {level: g_at(drag, ue) for level, ue in UE_LEVELS.items()}
        for e0 in E0:
            channel = name == CHANNEL_GAIN['target'] and e0 == CHANNEL_GAIN['E0_keV']
            for gain, label in ((1., ''), (CHANNEL_GAIN['gain'], ', channelled <110>')) if channel else ((1., ''),):
                ue, g_max, ue_at_max = ue_for_unity(drag, e0, gain)
                row = dict(target=name+label, E0_lab_keV=e0,
                           drag_per_T_at_E0=round(float(np.interp(e0, E_GRID, drag)), 4),
                           break_even_local_b=round(float(np.interp(e0, E_GRID, break_even_barns(drag))), 2))
                for level in UE_LEVELS:
                    g = gain*returns[level][e0]
                    row[f'G_{level}'] = g
                    row[f'drag_cut_for_G1_{level}'] = 1/g
                row['Ue_for_G1_keV'] = None if ue is None else round(ue, 2)
                row['G_max_over_Ue'] = g_max
                row['drag_cut_for_G1_at_G_max'] = None if ue is not None else 1/g_max
                ue_a, g_max_a, _ = ue_for_unity(drag, e0, gain, ret=thick_target_return_form_a)
                row['Ue_for_G1_keV_formA'] = None if ue_a is None else round(ue_a, 2)
                row['G_max_over_Ue_formA'] = g_max_a
                for f, grid in floors.items():
                    ue_f, g_max_f, _ = ue_for_unity(floor_drags[f][name], e0, gain, grid)
                    row[f'Ue_for_G1_keV_floor_{f:.3g}'] = None if ue_f is None else round(ue_f, 2)
                    row[f'G_max_over_Ue_floor_{f:.3g}'] = g_max_f
                rows.append(row)

    with open(OUT/'combined.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    summary = dict(check_LiT_break_even_1keV_b=round(check, 2), check_TiT1p7_break_even_1keV_b=round(check_ti, 2),
                   E_floor_lab_keV=E_FLOOR, Bosch_Hale_fit_floor_lab_keV=FIT_FLOOR, Ue_levels_keV=UE_LEVELS, channelling=CHANNEL_GAIN, rows=rows)
    (OUT/'combined.json').write_text(json.dumps(summary, indent=1)+'\n')
    for r in rows:
        print(f"{r['target']:32s} {r['E0_lab_keV']:5.0f}  G(.04)={r['G_computed_static']:.2e}  "
              f"G(.3)={r['G_measured_low_loading']:.2e}  G(1.4)={r['G_fitted_UC_Davis']:.2e}  "
              f"Ue_G1={r['Ue_for_G1_keV']}  Gmax={r['G_max_over_Ue']:.3g}  Ue_G1_A={r['Ue_for_G1_keV_formA']}  GmaxA={r['G_max_over_Ue_formA']:.3g}")


if __name__ == '__main__':
    main()
