"""How much apparent Ue does a near-surface D-density excess mimic?
Thick-target yield with stopping ~ sqrt(E) (as in work/alr.py), paper Eq. 7
cross section, constant S.  A depth profile that raises the D density probed
at E relative to the 6.5 keV reference by factor F multiplies R(E) by ~F.
We find the Ue shift that multiplies R(E)=Y(E)/Y(6.5) by the same F.
NumPy/SciPy only.  python3 -I density_vs_ue.py
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

EG = 986.0  # keV, D+D Gamow energy

def sig(E, U):
    return np.exp(-np.sqrt(EG / (E + U))) / np.sqrt(E * (E + U))

def Y(E, U):
    return quad(lambda e: sig(e, U) / np.sqrt(e), 1e-4, E, limit=200)[0]

def R(E, U):
    return Y(E, U) / Y(6.5, U)

for U0 in (0.3, 1.1):
    for E in (2.5, 1.0, 0.5):
        for F in (1.5, 2.0, 2.2, 4.0):
            f = lambda d: R(E, U0 + d) / R(E, U0) - F
            try:
                d = brentq(f, 0.0, 5.0)
                print(f"U0={U0:.1f} keV  E={E:3.1f} keV  density factor F={F:3.1f} "
                      f"<-> dUe = {d*1e3:6.0f} eV")
            except ValueError:
                print(f"U0={U0} E={E} F={F}: no root below 5 keV")

# A flat yield floor read as a per-energy Ue: solve sigma(E,U)/sigma(6.5,U) = R
print("\nPer-energy Ue implied by a relative yield floor R at E (cross-section ratio, Eq. 7):")
for E in (0.25, 0.5, 1.0, 2.0):
    for Rf in (1e-4, 1e-3, 1e-2, 5e-2):
        g = lambda U: np.log(sig(E, U) / sig(6.5, U)) - np.log(Rf)
        try:
            U = brentq(g, 0.0, 50.0)
            print(f"  E={E:4.2f} keV  R={Rf:7.0e}  ->  Ue = {U*1e3:6.0f} eV")
        except ValueError:
            print(f"  E={E:4.2f} keV  R={Rf:7.0e}  ->  no solution < 50 keV")
