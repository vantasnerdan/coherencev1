"""Q3 (Pd): can back-side loading reach the vacuum surface within the experiment? Beam off, loaded back side.
Writes reach.csv: time series of the vacuum-side solute x and the hydride front position."""
import os, sys, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import transport as T

CASES = [
    ("Pd clean front 300K", dict(material="Pd", s_front=0.5, T=300.0), 6 * 3600.0),
    ("Pd front s=1e-4 300K", dict(material="Pd", s_front=1e-4, T=300.0), 6 * 3600.0),
    ("Pd front s=1e-7 300K", dict(material="Pd", s_front=1e-7, T=300.0), 6 * 3600.0),
    ("Pd clean front eta=0.25", dict(material="Pd", s_front=0.5, T=300.0, eta_abs=0.25), 6 * 3600.0),
]   # Ti: see reach_ti.py (coarse two-sided grid, TEND env var; 330 K is stiff, use Arrhenius scaling)


def job(a):
    name, kw, tend = a
    c = dict(loaded=True, beam=False, t_pre=0.0, t_run=tend, E_cm=2.5, grid_both=True)
    c.update(kw)
    m = T.Model(c)
    te = np.unique(np.concatenate([np.geomspace(60, tend, 60), [600, 1800, 3600, 7200, 14400, 21600, 86400]]))
    te = te[te <= tend]
    t0 = time.time()
    m.run_te(te)
    rows = []
    N = m.N
    thr = 0.3 if c["material"] == "Pd" else 0.75
    for k, t in enumerate(m.t_beam):
        x = m.sol.y[:N, k]
        idx = np.nonzero(x > thr)[0]
        zf = m.zc[idx[0]] if len(idx) else m.e[-1]
        inv = (x * m.h).sum() * m.nM
        rows.append((name, t, x[0], m.window_DM(5e-9, k)[0], (L := T.L_FOIL) - zf, zf, inv))
    return rows, time.time() - t0


if __name__ == "__main__":
    with Pool(len(CASES)) as p:
        R = p.map(job, CASES)
    with open(os.path.join(HERE, "reach_pd.csv"), "w") as f:
        f.write("case,t_s,x_vacuum_surface_solute,DM_top5nm,hydride_front_depth_from_back_m,front_depth_from_vacuum_m,"
                "D_inventory_per_m2\n")
        for rows, wall in R:
            for r in rows:
                f.write(f"{r[0]},{r[1]:.4g},{r[2]:.4e},{r[3]:.4e},{r[4]:.4e},{r[5]:.4e},{r[6]:.4e}\n")
            print(rows[0][0], "wall", round(wall, 1))
