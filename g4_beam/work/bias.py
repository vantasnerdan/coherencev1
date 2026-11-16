import numpy as np
from scipy.optimize import minimize_scalar
exec(open('work/alr.py').read().split("Es=np.array")[0])
Ef=np.array([2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5])
def lsig(Ue): s=np.array([sigma(e,Ue) for e in Ef]); return np.log(s/s[-1])
def ltt(Ue):  y=np.array([Y(e,Ue) for e in Ef]); return np.log(y/y[-1])
def fit(data):
    return minimize_scalar(lambda u: np.sum((lsig(u)-data)**2),bounds=(0,6),method='bounded').x
print("A) thick-target truth, fitted with cross-section-ratio model (2.5-6.5 keV):")
for Ut in [0,0.1,0.3,0.5,0.8]:
    print(f"   Ue_true={Ut:.1f} -> Ue_fit={fit(ltt(Ut)):.2f} keV")
cup=np.where(Ef>=2.75,(6.5/Ef)**1.5,(6.5/2.75)**1.5)
print("B) plus Faraday-cup collected fraction ~V^1.5 used as charge (inflates R by (6.5/E)^1.5):")
for Ut in [0,0.1,0.3]:
    print(f"   Ue_true={Ut:.1f} -> Ue_fit={fit(ltt(Ut)+np.log(cup)):.2f} keV")
# digitized data fit with each model
d={'Pd':[7.04e-3,1.892e-2,4.183e-2,9.706e-2,0.21716,0.334,0.42129,0.68063,1],
   'PdD':[2.287e-2,4.389e-2,8.205e-2,0.15606,0.29025,0.40896,0.52445,np.nan,1],
   'Ti':[2.14e-3,1.062e-2,4.694e-2,9.132e-2,0.12954,0.27505,0.51398,0.62813,1],
   'TiD':[6.68e-3,2.429e-2,np.nan,np.nan,0.18624,0.35037,np.nan,0.76936,1]}
print("C) digitized Fig 2 points, 2.5-6.5 keV, unweighted log fit:")
for m,v in d.items():
    v=np.log(np.array(v)); ok=np.isfinite(v)
    us=minimize_scalar(lambda u: np.sum((lsig(u)-v)[ok]**2),bounds=(0,6),method='bounded').x
    ut=minimize_scalar(lambda u: np.sum((ltt(u)-v)[ok]**2),bounds=(0,6),method='bounded').x
    print(f"   {m}: cross-section model Ue={us:.2f}  thick-target model Ue={ut:.2f} keV")
