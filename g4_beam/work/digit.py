import numpy as np, sys
from PIL import Image
from scipy import ndimage as ndi
im=np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(int)
X0={'a':103.5,'b':1073.5}; XS=(1012-103.5)/6.5
Y0=51.5; YS=(736.5-51.5)/4.0
def classify(px):
    r,g,b=px[...,0],px[...,1],px[...,2]
    out={}
    out['red']=(r>120)&(r<190)&(g<40)&(b<40)
    out['blue']=(b>100)&(r<40)&(g<40)
    gray=(abs(r-g)<6)&(abs(g-b)<6)
    out['dgray']=gray&(r>85)&(r<120)
    out['lgray']=gray&(r>=125)&(r<165)
    return out
masks=classify(im)
for panel,(xlo,xhi) in {'a':(105,1011),'b':(1075,1981)}.items():
  for name,m in masks.items():
    mm=np.zeros_like(m); mm[53:735,xlo:xhi]=m[53:735,xlo:xhi]
    op=ndi.binary_opening(mm,structure=np.ones((9,9)))
    lab,n=ndi.label(op)
    res=[]
    for i in range(1,n+1):
        ys,xs=np.where(lab==i)
        if len(ys)<60: continue
        cx,cy=xs.mean(),ys.mean()
        E=(cx-X0[panel])/XS; logy=-(cy-Y0)/YS
        res.append((round(E,3),round(10**logy,5),len(ys)))
    if res: print(panel,name,sorted(res))
