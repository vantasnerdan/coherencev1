"""Driver: base profiles (profiles.csv, summary.csv), knob table (knobs.csv), back-side reach (reach.csv), figure.
Run: python3 run_all.py [base|knobs|reach|fig|all]
"""
import os
import sys
import time
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import transport as T  # noqa: E402

COND = {"Pd": ("Pd", False), "PdD": ("Pd", True), "Ti": ("Ti", False), "TiD": ("Ti", True)}
UE_PAPER = {"Pd": 0.8, "PdD": 1.7, "Ti": 0.3, "TiD": 0.9}        # paper's fitted constant Ue (keV)
ECM = [0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5]
ZOUT = np.concatenate([np.arange(0.125, 20, 0.25), np.arange(20.5, 100, 1.0), np.arange(102.5, 300, 5.0)])
PER_C = 1e-13 / T.QE          # (mb * nm^-2) -> probability, per coulomb of charge


def nominal_layers(cond):
    mat = "Pd" if cond.startswith("Pd") else "Ti"
    return T.uniform_layers(mat, 1.0 if mat == "Pd" else 2.0), 0.0


def sample_profile(m, k):
    layers, nads, st = m.profile(k)
    out = np.zeros((len(ZOUT), 4))       # D/M, nD (nm^-3), solute x, vacancy fraction
    zt = np.array([L[0] for L in layers]) * 1e9
    zb = np.array([L[1] for L in layers]) * 1e9
    for j, z in enumerate(ZOUT):
        i = np.searchsorted(zb, z)
        if i >= len(layers):
            continue
        a, b, nh, nd, Z, M, nmet = layers[i]
        out[j, 0] = nd / nmet
        out[j, 1] = nd * 1e-27
        if m.ol is not None and st["d"] > 1e-12 and i == 0:
            out[j, 2] = np.nan
            out[j, 3] = np.nan
        else:
            ii = np.clip(np.searchsorted(m.e + st["d"], z * 1e-9) - 1, 0, m.N - 1)
            out[j, 2] = st["x"][ii]
            out[j, 3] = st["v"][ii]
    return out


def yields(m, k, cond):
    prof = m.profile(k)
    layers, nads, st = prof
    res = {}
    for tag, Ue in (("Ue0", 0.0), ("Ue1", 1.0), ("UeP", UE_PAPER[cond])):
        y = T.yield_per_charge(m, Ue, k, profile=prof) * PER_C
        yn = T.yield_per_charge(m, Ue, k, profile=(layers, 0.0, st)) * PER_C
        res[tag] = y
        res[tag + "_noads"] = yn
        res[tag + "_nom"] = T.yield_per_charge(m, Ue, k, profile=(*nominal_layers(cond), st)) * PER_C
        # same mean D content over the top 70 nm, spread uniformly (shape effect only)
        dm, nd, v = m.window_DM(70e-9, k)
        uni = T.uniform_layers(m.mat["name"], dm)
        res[tag + "_uni70"] = T.yield_per_charge(m, Ue, k, profile=(uni, 0.0, st)) * PER_C
    return res


def summarize(m, cond, k):
    out = {}
    for dep in (5, 20, 70):
        dm, nd, v = m.window_DM(dep * 1e-9, k)
        out[f"DM{dep}"], out[f"nD{dep}"], out[f"v{dep}"] = dm, nd, v
    x, xt, v, d = m.state(k)
    out["d_ol_nm"] = d * 1e9
    nads, mus = m.surface_ads(x[0], d)
    out["ads_ML"] = nads / m.mat["ns"]
    out["x_surf_solute"] = x[0]
    return out


def run_case(args):
    tag, cond, case = args
    t0 = time.time()
    mat, loaded = COND[cond]
    c = dict(material=mat, loaded=loaded)
    c.update(case)
    m = T.Model(c)
    tr = m.c["t_run"]
    te = np.unique(np.concatenate([[600.0, 1800.0, 3600.0], np.linspace(0, tr, 13)[1:]]))
    te = te[te <= tr]
    try:
        m.run_te(te)
    except Exception as e:  # noqa: BLE001
        return dict(tag=tag, cond=cond, case=case, error=str(e))
    k30 = int(np.argmin(np.abs(m.t_beam - 1800.0)))
    res = dict(tag=tag, cond=cond, case=case, T=m.T, F=m.F, RN=m.RN, Rp_nm=m.Rp * 1.0, nv_ion=m.nv_ion,
               Ysp=m.Ysp, t_run=tr, wall=time.time() - t0)
    res["end"] = summarize(m, cond, -1)
    res["m30"] = summarize(m, cond, k30)
    res["yield"] = yields(m, -1, cond)
    res["prof"] = sample_profile(m, -1)
    layers, nads, state = m.profile(-1)
    res["profile_SI"] = dict(
        case=dict(condition=cond, E_cm_keV=float(m.c["E_cm"])),
        layers=[[float(value) for value in layer] for layer in layers],
        nads_m2=float(nads), state={key: np.asarray(value).tolist() for key, value in state.items()},
        end=res["end"])
    # overlayer lifetime: first time thickness < 0.1 nm
    if m.ol is not None:
        d = m.sol.y[3 * m.N]
        gone = np.nonzero(d < 0.1e-9)[0]
        res["ol_life_s"] = float(m.t_beam[gone[0]]) if len(gone) else np.inf
        res["Ysp_ol"] = m.Ysp_ol
    else:
        res["ol_life_s"] = np.nan
    return res


def base_jobs():
    return [(f"base_{c}_{E}", c, dict(E_cm=E)) for c in COND for E in ECM]


KNOBS = [
    ("ALD Al2O3 2nm", dict(overlayer="Al2O3", d_ol=2e-9)),
    ("ALD Al2O3 5nm", dict(overlayer="Al2O3", d_ol=5e-9)),
    ("ALD Al2O3 1nm", dict(overlayer="Al2O3", d_ol=1e-9)),
    ("ALD TiO2 2nm", dict(overlayer="TiO2", d_ol=2e-9)),
    ("Pt 2nm", dict(overlayer="Pt", d_ol=2e-9)),
    ("Pd cap 2nm", dict(overlayer="Pd", d_ol=2e-9)),
    ("current x2", dict(j_rel=2.0)),
    ("T +30 K", dict(dT=30.0)),
    ("T -30 K", dict(dT=-30.0)),
    ("strain +0.5%", dict(strain=0.005)),
    ("strain -0.5%", dict(strain=-0.005)),
    ("PdAg23 alloy", dict(alloy="PdAg23")),
    ("PdPt10 alloy", dict(alloy="PdPt10")),
    # model uncertainties, not knobs
    ("[unc] vsat 0.003", dict(vsat=0.003)),
    ("[unc] vsat 0.03", dict(vsat=0.03)),
    ("[unc] front sticking 1e-4", dict(s_front=1e-4)),
    ("[unc] front sticking 1e-7", dict(s_front=1e-7)),
    ("[unc] eta_abs 0.01", dict(eta_abs=0.01)),
    ("[unc] eta_abs 0.25", dict(eta_abs=0.25)),
    ("[unc] no chamber gas", dict(no_gas=True)),
]


def knob_jobs():
    jobs = []
    for c in COND:
        mat, loaded = COND[c]
        for E in (0.25, 2.5):
            jobs.append((f"knob_{c}_{E}_base", c, dict(E_cm=E)))
            for name, kw in KNOBS:
                if name == "Pd cap 2nm" and mat == "Pd":
                    continue
                if "alloy" in name and mat == "Ti":
                    continue
                if (name == "current x2" or "eta_abs" in name) and not loaded:
                    continue
                if "front sticking" in name and mat == "Ti":
                    continue
                kk = dict(E_cm=E)
                kk.update(kw)
                jobs.append((f"knob_{c}_{E}_{name}", c, kk))
    return jobs


def pool_run(jobs, nproc=12):
    with Pool(nproc) as p:
        return p.map(run_case, jobs, chunksize=1)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    import pickle
    if what in ("base", "all"):
        t = time.time()
        R = pool_run(base_jobs())
        pickle.dump(R, open(os.path.join(HERE, "base.pkl"), "wb"))
        print("base done", time.time() - t, [r["tag"] for r in R if "error" in r])
    if what in ("knobs", "all"):
        t = time.time()
        R = pool_run(knob_jobs())
        pickle.dump(R, open(os.path.join(HERE, "knobs.pkl"), "wb"))
        print("knobs done", time.time() - t, [r["tag"] for r in R if "error" in r])
