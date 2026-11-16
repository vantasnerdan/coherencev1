import numpy as np, sys
sys.path.insert(0,'work')
from PIL import Image
from scipy.optimize import brentq, minimize_scalar
exec(open('work/alr.py').read().split("Es=np.array")[0])
im=np.asarray(Image.open('paper/fig-001.png').convert('RGB')).astype(int)
XS=(1012-103.5)/6.5; Y0=51.5; YS=(736.5-51.5)/4
# find where solid curves reach the bottom axis: scan row just above bottom frame
row=im[730]
def runs(mask,x0):
    xs=np.where(mask)[0]; 
    if len(xs)==0: return []
    groups=np.split(xs,np.where(np.diff(xs)>3)[0]+1)
    return [round((g.mean()-x0)/XS,2) for g in groups]
r,g,b=row[:,0],row[:,1],row[:,2]
print("panel a row730 red-ish x(keV):",runs((r>110)&(g<60)&(b<60)&(np.arange(len(r))<1012),103.5))
print("panel a row730 black x(keV):",runs((r+g+b<120)&(np.arange(len(r))>110)&(np.arange(len(r))<1005),103.5))
print("panel a row730 gray x(keV):",runs((abs(r-g)<8)&(abs(g-b)<8)&(r>60)&(r<130)&(np.arange(len(r))>110)&(np.arange(len(r))<1005),103.5))
print("panel b row730 blue x(keV):",runs((b>100)&(r<60)&(np.arange(len(r))>1080),1073.5))
print("panel b row730 black x(keV):",runs((r+g+b<120)&(np.arange(len(r))>1080)&(np.arange(len(r))<1975),1073.5))
print("panel b row730 gray x(keV):",runs((abs(r-g)<8)&(abs(g-b)<8)&(r>60)&(r<150)&(np.arange(len(r))>1080)&(np.arange(len(r))<1975),1073.5))
yrow=-(730-Y0)/YS; print("row730 corresponds to log10 R =",round(yrow,3))
target=10**yrow
def Rsig(E,Ue): return sigma(E,Ue)/sigma(6.5,Ue)
def Rtt(E,Ue): return Y(E,Ue)/Y(6.5,Ue)
for Ue in [0,0.3,0.8,0.9,1.7]:
    es=brentq(lambda e: np.log10(Rsig(e,Ue)/target),0.02,6.4)
    et=brentq(lambda e: np.log10(Rtt(e,Ue)/target),0.02,6.4)
    print(f"Ue={Ue}: cross-section-ratio model reaches row at {es:.2f} keV; thick-target model at {et:.2f} keV")
