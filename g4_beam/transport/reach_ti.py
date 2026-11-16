"""Q3 for Ti on a coarser two-sided grid (50 nm first cells): delta-hydride front from the loaded back side."""
import os, sys, time
import numpy as np
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import transport as T
TEND = float(os.environ.get("TEND", 8.64e5))


def job(Tk):
    c = dict(material="Ti", loaded=True, beam=False, t_pre=0.0, t_run=TEND, E_cm=2.5, grid_both=True,
             grid=dict(h0=50e-9, g1=1.05, hmid=0.5e-6, g2=1.12), T=Tk)
    m = T.Model(c)
    te = np.geomspace(60, TEND, 40)
    t0 = time.time()
    m.run_te(te)
    rows = []
    for k, t in enumerate(m.t_beam):
        x = m.sol.y[:m.N, k]
        low = np.nonzero(x < 0.75)[0]           # hydride front from the back face: last non-hydride cell
        front_from_back = T.L_FOIL - m.e[low[-1] + 1] if len(low) else T.L_FOIL
        # vacuum-side hydride grown from chamber gas (beam off), for reference
        hi = np.nonzero(x > 0.75)[0]
        vac_side = m.e[hi[hi < m.N // 2][-1] + 1] if np.any(hi < m.N // 2) else 0.0
        rows.append((f"Ti {Tk:.0f}K", t, x[0], front_from_back, vac_side, x[m.N // 2]))
    return rows, time.time() - t0


if __name__ == "__main__":
    with Pool(2) as p:
        R = p.map(job, [300.0, 330.0])
    with open(os.path.join(HERE, "reach_ti.csv"), "w") as f:
        f.write("case,t_s,x_vacuum_side_surface,delta_front_depth_from_back_m,delta_depth_from_vacuum_side_gas_m,x_midplane\n")
        for rows, wall in R:
            for r in rows:
                f.write(f"{r[0]},{r[1]:.4g},{r[2]:.4e},{r[3]:.4e},{r[4]:.4e},{r[5]:.4e}\n")
            print(rows[0][0], "wall", round(wall, 1))
