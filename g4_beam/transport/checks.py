"""Numerical checks: (1) grid two-point convergence, (2) D mass balance, (3) BCA ion-count convergence."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import transport as T
import bca

out = []
# (1) grid: default vs halved first cell and slower growth
for cond, kw in (("Pd", dict(material="Pd", E_cm=2.5)), ("PdD", dict(material="Pd", loaded=True, E_cm=2.5)),
                 ("Ti", dict(material="Ti", E_cm=0.25)), ("Pd", dict(material="Pd", E_cm=0.25))):
    vals = []
    for g in (None, dict(h0=0.125e-9, g1=1.03, hmid=2e-9, g2=1.08)):
        m = T.Model(dict(kw, grid=g)).run(n_out=4)
        d5, n5, v5 = m.window_DM(5e-9)
        d70 = m.window_DM(70e-9)[0]
        y = T.yield_per_charge(m, 0.0)
        vals.append((m.N, d5, d70, v5, y))
    (n1, a1, b1, c1, y1), (n2, a2, b2, c2, y2) = vals
    out.append(f"grid {cond} E={kw['E_cm']}: N {n1}->{n2}; DM5 {a1:.4g}->{a2:.4g} ({a2/a1-1:+.1%}); "
               f"DM70 {b1:.4g}->{b2:.4g} ({b2/b1-1:+.1%}); v5 {c1:.3g}->{c2:.3g}; Y(Ue=0) {y2/y1-1:+.1%}")
# (2) mass balance: Pd, no desorption, no gas, zero back flux; inventory vs retained beam dose
m = T.Model(dict(material="Pd", E_cm=2.5, s_front=1e-30, no_gas=True, t_run=600.0)).run(n_out=2)
x, xt, v, d = m.state()
inv = ((x - 1e-9 + xt) * m.h).sum() * m.nM
dose = m.F * (1 - m.RN) * 600.0
out.append(f"mass balance Pd 600 s, closed surfaces: inventory/retained dose = {inv/dose:.4f} "
           f"(erosion removes the rest)")
m = T.Model(dict(material="Ti", E_cm=2.5, no_gas=True, t_run=200.0, native_oxide=0.0)).run(n_out=2)
x, xt, v, d = m.state()
inv = ((x - 1e-6 + xt) * m.h).sum() * m.nM
dose = m.F * (1 - m.RN) * 200.0
out.append(f"mass balance Ti 200 s, no gas: inventory/retained dose = {inv/dose:.4f}")
# (3) BCA ion count
for E in (500, 5000):
    r1 = bca.run(E, [(46, 106.42, 1, 34.0)], 68.0, nions=4000, seed=3)
    r2 = bca.run(E, [(46, 106.42, 1, 34.0)], 68.0, nions=16000, seed=4)
    out.append(f"BCA Pd {E} eV: Rp {r1['Rp']:.2f}->{r2['Rp']:.2f} nm, R_N {r1['RN']:.3f}->{r2['RN']:.3f}, "
               f"vac/ion {r1['nv_per_ion']:.3f}->{r2['nv_per_ion']:.3f} (4000 -> 16000 ions)")
open(os.path.join(HERE, "checks.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
