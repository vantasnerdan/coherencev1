"""Export full end-of-run depth profiles (Model.profile) for the baseline cases as JSON, for the forward model.

Run: python3 export_profiles.py   (16 cases, about 2 CPU-minutes on 12 cores; uses bca_tables.npz if present)
Writes baseline/profiles_full.json. Each record holds, at the last stored time of the run:
  case (condition, E_cm keV, beam time s, foil temperature K),
  layers: list of [z_top_m, z_bot_m, n_host_m3, n_D_m3, Z_host, M_host_amu, n_metal_m3] from the surface down
          (first layer is the overlayer when one is present; none in the baseline),
  nads_m2: chemisorbed D areal density (m^-2), counted at the surface,
  state: solute x, trapped xt, vacancy fraction v per host cell, overlayer thickness d (m), cell centres zc (m),
  end: the same D/M window summaries as summary.csv, for a consistency check against base runs.
"""
import json
import os
import sys
import time
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import transport as T  # noqa: E402
from run_all import COND, summarize  # noqa: E402

CASES = [(c, e) for c in ("Pd", "PdD", "Ti", "TiD") for e in (0.25, 0.5, 2.5, 6.5)]


def one(args):
    cond, E = args
    t0 = time.time()
    mat, loaded = COND[cond]
    m = T.Model(dict(material=mat, loaded=loaded, E_cm=E))
    tr = m.c["t_run"]
    te = np.unique(np.concatenate([[600.0, 1800.0, 3600.0], np.linspace(0, tr, 13)[1:]]))
    m.run_te(te[te <= tr])
    layers, nads, st = m.profile(-1)
    return dict(
        case=dict(condition=cond, E_cm_keV=E, beam_time_s=float(m.t_beam[-1]), T_K=float(m.T)),
        layers=[[float(v) for v in L] for L in layers],
        nads_m2=float(nads),
        state={k: (np.asarray(v, float).tolist() if np.ndim(v) else float(v)) for k, v in st.items()},
        end={k: float(v) for k, v in summarize(m, cond, -1).items()},
        wall_s=time.time() - t0,
    )


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "baseline"), exist_ok=True)
    with Pool(min(12, os.cpu_count() or 1)) as p:
        out = p.map(one, CASES, chunksize=1)
    with open(os.path.join(HERE, "baseline", "profiles_full.json"), "w") as f:
        json.dump(dict(units="SI (m, m^-3, m^-2, s, K); E_cm in keV", records=out), f)
    for r in out:
        print(r["case"]["condition"], r["case"]["E_cm_keV"], "DM5 %.4g DM70 %.4g layers %d wall %.1fs"
              % (r["end"]["DM5"], r["end"]["DM70"], len(r["layers"]), r["wall_s"]))
