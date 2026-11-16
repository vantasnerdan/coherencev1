import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
BG=31.3970  # sqrt(keV), Bosch-Hale DD
EG=BG**2
def S_p(E):  # Bosch-Hale D(d,p)T S-factor keV*mb, E in keV cm
    return 5.5576e4+E*(2.1054e2+E*(-3.2638e-2+E*(1.4987e-6+E*1.8181e-10)))
def sigma(E,Ue):  # paper Eq 7
    return S_p(E)/np.sqrt(E*(E+Ue))*np.exp(-np.sqrt(EG/(E+Ue)))
def Y(E,Ue,stop=lambda e:np.sqrt(e)):   # thick target, stopping ~ sqrt(E)
    lo=max(1e-3,E*1e-3)
    return quad(lambda e: sigma(e,Ue)/stop(e), lo, E, limit=200, points=[E*0.5,E*0.9])[0]
Es=np.array([0.25,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5])
def R(Ue): 
    y=np.array([Y(e,Ue) for e in Es]); return y/y[-1]
print("E cm keV:", Es)
for Ue in [0,0.3,0.8,0.9,1.7,3.0]:
    print(f"R(E) Ue={Ue}:", " ".join(f"{v:.2e}" for v in R(Ue)))
# energy where each curve crosses 1e-4
for Ue in [0,0.3,0.8,0.9,1.7]:
    f=lambda e: np.log10(Y(e,Ue)/Y(6.5,Ue))+4
    print(f"Ue={Ue}: crosses 1e-4 at E={brentq(f,0.05,6.4):.2f} keV")
# cross-section enhancement at 0.25 keV (Eq 12)
for Ue in [2.14,2.25,2.46,3.0]:
    E=0.25; enh=np.sqrt(E/(E+Ue))*np.exp(np.sqrt(EG/E)-np.sqrt(EG/(E+Ue)))
    print(f"Eq12 enhancement at 0.25 keV, Ue={Ue}: {enh:.2e}")
np.save('work/Es.npy',Es)
