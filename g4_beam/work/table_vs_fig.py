import numpy as np
E=np.array([3.0,2.5,2.0,1.5,1.0,0.5,0.25])
I=np.array([145,138,134,123,110,68,38.])
SB={'Pd':[443,121,83,53,44,15,21],'PdD':[465,177,138,92,43,32,30],
    'Ti':[1525,225,82,106,69,47,15],'TiD':[5427,1082,368,261,112,62,56]}
fig={'Pd':[1.89e-2,7.04e-3,5.14e-3,3.84e-3,4.12e-3,3.60e-3,1.656e-2],
     'PdD':[4.389e-2,2.287e-2,1.896e-2,1.466e-2,8.98e-3,1.668e-2,5.052e-2],
     'Ti':[1.062e-2,2.14e-3,8.3e-4,1.25e-3,1.09e-3,1.8e-3,1.95e-3],
     'TiD':[2.429e-2,6.68e-3,2.42e-3,2.01e-3,np.nan,1.51e-3,4.58e-3]}
bg=0.007
print("E    :",E)
print("I uA :",I)
for m in SB:
    sb=np.array(SB[m],float); f=np.array(fig[m])
    rate=sb*bg                      # cph
    ypc=rate/(I*1e-6*3600)          # counts per Coulomb
    rf=f/f[0]; r1=(sb/I)/(sb[0]/I[0]); r2=(sb/I**2)/(sb[0]/I[0]**2)
    print(f"\n{m}: signal rate cph   ", np.round(rate,3))
    print(f"{m}: counts/C (table)   ", np.round(ypc,2))
    print(f"{m}: counts in 12-24 h   ", np.round(rate*12,1), np.round(rate*24,1))
    print(f"{m}: fig  Y(E)/Y(3)      ", np.round(rf,3))
    print(f"{m}: tab  (SB/I) ratio   ", np.round(r1,3), " fig/tab:",np.round(rf/r1,2))
    print(f"{m}: tab  (SB/I^2) ratio ", np.round(r2,3), " fig/tab:",np.round(rf/r2,2))
