"""Order-of-magnitude comparisons between Karahadian et al. 2026 and other low-energy d+d datasets.
All numbers printed here are my estimates built from the cited inputs; see report.md for sources.
Run: python3 -I compare.py <packet root>
"""
import sys
import numpy as np
from scipy.integrate import quad

S = sys.argv[1]
EG_BH = 31.3970 ** 2          # Bosch-Hale Gamow energy for DD, keV (986.0)


def S_p(E):                   # Bosch-Hale D(d,p)T S-factor, keV mb, E keV cm
    return 5.5576e4 + E * (2.1054e2 + E * (-3.2638e-2 + E * (1.4987e-6 + E * 1.8181e-10)))


def sigma(E, Ue):
    return S_p(E) / np.sqrt(E * (E + Ue)) * np.exp(-np.sqrt(EG_BH / (E + Ue)))


def Ytt(E, Ue):               # thick target, stopping ~ sqrt(E), E in keV cm
    return quad(lambda e: sigma(e, Ue) / np.sqrt(e), max(1e-3, 1e-3 * E), E, limit=200)[0]


# ---------- 1. Karahadian floors (Fig 2 digitized, relative to own 6.5 keV point)
d = np.genfromtxt(f"{S}/work/fig2_digitized.csv", delimiter=",", names=True, skip_header=2)
E = d["E_cm"]
print("1. Karahadian Fig 2 floors (mean of 0.5-2.0 keV points; relative to Y(6.5) and to Y(3.0))")
kfloor = {}
for m in ("Pd", "PdD", "Ti", "TiD"):
    y = d[m]
    sel = (E >= 0.5) & (E <= 2.0) & np.isfinite(y)
    fl = np.exp(np.mean(np.log(y[sel])))
    kfloor[m] = fl
    y3 = y[E == 3.0][0]
    print(f"   {m:4s} floor/Y(6.5)={fl:.2e}  floor/Y(3.0)={fl / y3:.2e}  Y(1.5)/Y(3.0)={y[E == 1.5][0] / y3:.3f}")

# ---------- 2. Czerski et al. plateau levels converted to a 'relative to own yield at E_cm=6.5 keV' scale
# Table 2 of arXiv:2409.02112v3 / 2605.27438: normalisation M, Ue (eV), plateau (counts per uA per ks).
# Their Eq. 23: Y(E) = M exp(-sqrt(EG/(E+Ue))), EG = 986 keV, E = cm energy.
EGc = 986.0
cz = {"Zr": (14.6e6, 0.340, 2.6e-4), "Ti": (12.69e6, 0.025, 1.4e-4), "Pd": (3.9e5, 0.400, 6.9e-3)}
print("\n2. Czerski et al. plateau relative to their own fitted screening curve (my estimate)")
for m, (M, Ue, pl) in cz.items():
    Y65 = M * np.exp(-np.sqrt(EGc / (6.5 + Ue)))
    Y30 = M * np.exp(-np.sqrt(EGc / (3.0 + Ue)))
    print(f"   {m}: Y(6.5)={Y65:.3g}  plateau/Y(6.5)={pl / Y65:.2e}  plateau/Y(3.0)={pl / Y30:.2e}")
print("   compare Karahadian Pd floor/Y(6.5) = %.1e, Ti = %.1e" % (kfloor["Pd"], kfloor["Ti"]))

# ---------- 3. Yuki et al. 1997 (Kasagi group) Ti thick target, EXFOR O1084 subentry 003
# E_d (lab keV) and yield relative to E_d = 6.45 keV
yuki_ti = np.array([[2.566, 4.1e-5], [2.765, 8.2e-5], [3.066, 3.3e-4], [3.272, 8.0e-4], [3.48, 2.21e-3],
                    [3.667, 2.87e-3], [3.858, 5.89e-3], [4.047, 9.2e-3], [4.263, 1.7e-2], [4.451, 2.46e-2],
                    [4.667, 4.4e-2], [4.855, 6.71e-2], [5.056, 0.1], [5.243, 0.14], [5.431, 0.203],
                    [5.631, 0.27], [5.845, 0.411], [6.057, 0.505], [6.244, 0.656]])
Ecm_y = yuki_ti[:, 0] / 2
r_y = np.interp(1.53, Ecm_y, yuki_ti[:, 1]) / np.interp(3.12, Ecm_y, yuki_ti[:, 1])
print("\n3. Ti: Y(E_cm~1.5)/Y(E_cm~3.1)")
print(f"   Yuki 1997 Ti (EXFOR O1084): {r_y:.2e}  (lowest point E_cm={Ecm_y[0]:.2f} keV, rel. yield {yuki_ti[0, 1]:.1e})")
for Ue in (0.0, 0.065, 0.3, 0.9):
    print(f"   thick-target model Ue={Ue:.3f} keV: {Ytt(1.53, Ue) / Ytt(3.12, Ue):.2e}")
print(f"   Karahadian Ti Fig 2: {d['Ti'][E == 1.5][0] / d['Ti'][E == 3.0][0]:.2e}")
# Yuki Ti relative yield at its lowest point vs Karahadian Ti floor on the same E_d=6.45 keV scale
k_ti_rel_323 = np.interp(3.225, E, d["Ti"])
print(f"   Karahadian Ti floor relative to its own Y(E_cm=3.225): {kfloor['Ti'] / k_ti_rel_323:.2e}"
      f"  vs Yuki Ti at E_cm=1.28 relative to E_cm=3.225: {yuki_ti[0, 1]:.1e}")

# ---------- 4. Lipson et al. 2005 glow discharge, Ti, E_d 0.8-2.45 keV, Ue = 610 eV
print("\n4. Lipson 2005: thick-target drop implied by Ue=0.61 keV between E_cm 1.225 and 0.4 keV")
print(f"   Y(0.4)/Y(1.225) = {Ytt(0.4, 0.61) / Ytt(1.225, 0.61):.1e}  (a flat floor would give ~1)")

# ---------- 5. Charge-exchange neutral fractions
kB = 1.380649e-23


def ncm3(p_pa, T=300.0):
    return p_pa / (kB * T) * 1e-6


print("\n5. Fast-neutral fraction f = 1-exp(-n sigma L)")
for label, p_pa, L in (("Karahadian chamber 3e-4 Torr, 4 cm", 3e-4 * 133.322, 4.0),
                       ("Karahadian, source-side 1e-3 Torr, 1 cm (assumed)", 1e-3 * 133.322, 1.0),
                       ("Czerski beamline 8e-9 mbar, 150 cm", 8e-9 * 100.0, 150.0)):
    n = ncm3(p_pa)
    fs = [1 - np.exp(-n * s * L) for s in (1e-16, 1e-15)]
    print(f"   {label}: n={n:.2e} cm^-3  f(sigma=1e-16)={fs[0]:.1e}  f(sigma=1e-15)={fs[1]:.1e}")
print("   Czerski measured neutral fraction ~3e-4 (D1+, D2+), arXiv:2409.02112v3")

# ---------- 6. Contaminant fraction needed if the floor were a fixed-energy component at E_x
# floor = f * Y(E_x) (+ nominal-energy yield, negligible below 2 keV). Optional correction for the
# possible extra 1/I factor in Fig 2 (brief item 5): divide floor by I(ref)/I(E) at the floor energies.
I_lab = {0.5: 38, 1: 68, 1.5: 94, 2: 110, 2.5: 118, 3: 123, 3.5: 130, 4: 134, 13: 238}  # SI Fig 3, uA vs kV
corr = np.mean([I_lab[13] / I_lab[2 * e] for e in (0.5, 1.0, 1.5, 2.0) if 2 * e in I_lab])
print(f"\n6. Fixed-energy contaminant fraction f needed (floor / Y(E_x)); mean I(ref)/I(E) over floor = {corr:.1f}")
for Ex in (2.5, 4.0, 5.0, 6.0):
    row = []
    for m in ("Pd", "PdD", "Ti", "TiD"):
        yx = np.interp(Ex, E[np.isfinite(d[m])], d[m][np.isfinite(d[m])])
        row.append(f"{m}:{kfloor[m] / yx:.2g}/{kfloor[m] / yx / corr:.2g}")
    print(f"   E_x={Ex} keV cm  (as plotted / with 1/I correction)  " + "  ".join(row))
