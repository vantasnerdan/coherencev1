"""Post-processing: profiles.csv, summary.csv, ratios.csv, knobs.csv, profiles.png from base.pkl / knobs.pkl."""
import os
import sys
import pickle
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from run_all import ZOUT, ECM  # noqa: E402

B = pickle.load(open(os.path.join(HERE, "base.pkl"), "rb"))
K = pickle.load(open(os.path.join(HERE, "knobs.pkl"), "rb")) if os.path.exists(os.path.join(HERE, "knobs.pkl")) else []
base = {(r["cond"], r["case"]["E_cm"]): r for r in B if "error" not in r}
errs = [r["tag"] + ": " + r["error"] for r in B + K if "error" in r]
if errs:
    print("ERRORS:", *errs, sep="\n  ")
CONDS = ["Pd", "PdD", "Ti", "TiD"]

# ---- profiles.csv (long format)
with open(os.path.join(HERE, "profiles.csv"), "w") as f:
    f.write("condition,E_cm_keV,T_K,beam_time_h,z_nm,D_per_M,nD_per_nm3,x_solute,vacancy_fraction\n")
    for c in CONDS:
        for E in ECM:
            r = base.get((c, E))
            if r is None:
                continue
            P = r["prof"]
            for j, z in enumerate(ZOUT):
                f.write(f"{c},{E},{r['T']:.1f},{r['t_run'] / 3600:.0f},{z:.3f},{P[j, 0]:.4e},{P[j, 1]:.4e},"
                        f"{P[j, 2]:.4e},{P[j, 3]:.4e}\n")

# ---- summary.csv
cols = ["DM5", "DM20", "DM70", "nD5", "nD20", "nD70", "v5", "v20", "v70", "d_ol_nm", "ads_ML", "x_surf_solute"]
ycols = ["Ue0", "Ue0_noads", "Ue0_nom", "Ue0_uni70", "UeP", "UeP_noads", "UeP_nom", "UeP_uni70", "Ue1"]
with open(os.path.join(HERE, "summary.csv"), "w") as f:
    f.write("condition,E_cm_keV,T_K,flux_D_m2s,R_N,Rp_nm,vac_per_D,Y_sputter,beam_time_h," + ",".join(cols) +
            "," + ",".join("Ypc_" + y for y in ycols) + "\n")
    for c in CONDS:
        for E in ECM:
            r = base.get((c, E))
            if r is None:
                continue
            f.write(f"{c},{E},{r['T']:.1f},{r['F']:.3e},{r['RN']:.3f},{r['Rp_nm']:.2f},{r['nv_ion']:.3f},"
                    f"{r['Ysp']:.4f},{r['t_run'] / 3600:.0f}," + ",".join(f"{r['end'][k]:.4e}" for k in cols) + "," +
                    ",".join(f"{r['yield'][y]:.4e}" for y in ycols) + "\n")

# ---- ratios.csv: hydride/metal and the profile factors
SB = {3.0: (145, 443, 465, 1525, 5427), 2.5: (138, 121, 177, 225, 1082), 2.0: (134, 83, 138, 82, 368),
      1.5: (123, 53, 92, 106, 261), 1.0: (110, 44, 43, 69, 112), 0.5: (68, 15, 32, 47, 62),
      0.25: (38, 21, 30, 15, 56)}


def meas_per_C(E, c):
    I, *sb = SB[E]
    return sb[CONDS.index(c)] * 0.007 / (I * 1e-6 * 3600)


with open(os.path.join(HERE, "ratios.csv"), "w") as f:
    f.write("E_cm_keV," + ",".join(f"shape_{c},content_{c},ads_share_{c}" for c in CONDS) +
            ",PdD/Pd_Ue0,PdD/Pd_Ue1,PdD/Pd_UeP,PdD/Pd_meas,TiD/Ti_Ue0,TiD/Ti_Ue1,TiD/Ti_UeP,TiD/Ti_meas,"
            "TiD/PdD_UeP,TiD/PdD_meas," + ",".join(f"eff_{c}" for c in CONDS) + "\n")
    for E in ECM:
        if any((c, E) not in base for c in CONDS):
            continue
        y = {c: base[(c, E)]["yield"] for c in CONDS}
        row = [E]
        for c in CONDS:
            row += [y[c]["UeP"] / y[c]["UeP_uni70"], y[c]["UeP"] / y[c]["UeP_nom"], 1 - y[c]["UeP_noads"] / y[c]["UeP"]]
        row += [y["PdD"]["Ue0"] / y["Pd"]["Ue0"], y["PdD"]["Ue1"] / y["Pd"]["Ue1"], y["PdD"]["UeP"] / y["Pd"]["UeP"]]
        row += [meas_per_C(E, "PdD") / meas_per_C(E, "Pd") if E in SB else np.nan]
        row += [y["TiD"]["Ue0"] / y["Ti"]["Ue0"], y["TiD"]["Ue1"] / y["Ti"]["Ue1"], y["TiD"]["UeP"] / y["Ti"]["UeP"]]
        row += [meas_per_C(E, "TiD") / meas_per_C(E, "Ti") if E in SB else np.nan]
        row += [y["TiD"]["UeP"] / y["PdD"]["UeP"], meas_per_C(E, "TiD") / meas_per_C(E, "PdD") if E in SB else np.nan]
        # implied detector efficiency: measured counts/C over model fusions/C (p branch), paper Ue per condition
        row += [meas_per_C(E, c) / (y[c]["UeP"]) if E in SB else np.nan for c in CONDS]
        f.write(",".join(f"{v:.4g}" for v in row) + "\n")

# ---- knobs.csv
if K:
    kb = {}
    for r in K:
        if "error" in r:
            continue
        name = r["tag"].split("_", 3)[3]
        kb[(r["cond"], r["case"]["E_cm"], name)] = r
    with open(os.path.join(HERE, "knobs.csv"), "w") as f:
        f.write("condition,E_cm_keV,knob,T_K,DM5_30min,DM5_end,DM20_end,DM70_end,nD5_end,v5_end,v70_end,overlayer_nm_30min,"
                "overlayer_life_h,Y_UeP_per_C_end,ratio_DM5_30min,ratio_DM5_end,ratio_DM70_end,ratio_v5_end,"
                "ratio_v70_end,ratio_Y_end,ratio_Y_30min_proxy\n")
        for (c, E, name), r in sorted(kb.items(), key=lambda t: (CONDS.index(t[0][0]), t[0][1], t[0][2] != "base")):
            b = kb.get((c, E, "base"))
            if b is None:
                continue

            def rt(a, bb):
                return a / bb if bb > 0 else np.nan
            life = r.get("ol_life_s", np.nan)
            f.write(f"{c},{E},{name},{r['T']:.1f},{r['m30']['DM5']:.4e},{r['end']['DM5']:.4e},{r['end']['DM20']:.4e},"
                    f"{r['end']['DM70']:.4e},{r['end']['nD5']:.4e},{r['end']['v5']:.4e},{r['end']['v70']:.4e},"
                    f"{r['m30']['d_ol_nm']:.3f},{life / 3600 if np.isfinite(life) else np.inf:.3g},"
                    f"{r['yield']['UeP']:.4e},"
                    f"{rt(r['m30']['DM5'], b['m30']['DM5']):.4g},{rt(r['end']['DM5'], b['end']['DM5']):.4g},"
                    f"{rt(r['end']['DM70'], b['end']['DM70']):.4g},{rt(r['end']['v5'], b['end']['v5']):.4g},"
                    f"{rt(r['end']['v70'], b['end']['v70']):.4g},{rt(r['yield']['UeP'], b['yield']['UeP']):.4g},"
                    f"{rt(r['m30']['nD5'], b['m30']['nD5']):.4g}\n")

# ---- figure
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
COL = {"Pd": "#2a78d6", "PdD": "#eb6834", "Ti": "#1baf7a", "TiD": "#eda100"}
LS = {"Pd": "-", "PdD": (0, (5, 2)), "Ti": "-", "TiD": (0, (5, 2))}
fig, axs = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True, facecolor=SURF)
for ax, E in zip(axs, (0.25, 2.5)):
    ax.set_facecolor(SURF)
    for c in CONDS:
        r = base.get((c, E))
        if r is None:
            continue
        P = r["prof"]
        sel = ZOUT <= 100
        ax.plot(ZOUT[sel], P[sel, 0], color=COL[c], ls=LS[c], lw=2, label=f"{c} ({r['T']:.0f} K)")
    # direct labels (Ti and TiD coincide; Pd and PdD offset to avoid collision)
    yl = {c: base[(c, E)]["prof"][ZOUT <= 100, 0][-1] for c in CONDS}
    ax.text(101, yl["TiD"], "Ti, TiD", color=INK2, fontsize=9, va="center")
    lo, hi = sorted([yl["Pd"], yl["PdD"]])
    ax.text(101, max(hi, 2e-5) * 1.35, "PdD" if yl["PdD"] >= yl["Pd"] else "Pd", color=INK2, fontsize=9, va="center")
    ax.text(101, max(lo, 1.2e-5) / 1.35, "Pd" if yl["PdD"] >= yl["Pd"] else "PdD", color=INK2, fontsize=9,
            va="center")
    ax.set_yscale("log")
    ax.set_xlim(0, 116)
    ax.set_ylim(1e-5, 4)
    ax.set_title(f"E_cm = {E} keV (E_lab {2 * E:g} keV), end of {base[('Pd', E)]['t_run'] / 3600:.0f} h beam",
                 color=INK, fontsize=11, loc="left")
    ax.set_xlabel("depth below vacuum surface (nm)", color=INK2)
    ax.grid(True, color=GRID, lw=0.8)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=INK2)
axs[0].set_ylabel("D per metal atom (solute + trapped; oxide per Ti)", color=INK2)
axs[1].legend(frameon=False, fontsize=9, loc="lower left", labelcolor=INK)
fig.suptitle("Model D/M depth profiles under the beam (clean Pd surface, native TiO2 on Ti, chamber D2 0.037 Pa;"
             " plus ~1 monolayer of chemisorbed D on every surface, not plotted)",
             color=INK, fontsize=10, x=0.01, ha="left")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "profiles.png"), dpi=150, facecolor=SURF)
print("written")

# ---- Ue bias: fit a constant Ue with a uniform stoichiometric target to the model's own R(E) over 2.5-6.5 keV
import transport as T  # noqa: E402
from scipy.optimize import minimize_scalar  # noqa: E402
from run_all import UE_PAPER, PER_C  # noqa: E402
EF = [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5]


def nominal_Y(mat, E, Ue):
    lay = T.uniform_layers(mat, 1.0 if mat == "Pd" else 2.0)
    El = 2e3 * E
    return 0.85 * T.yield_per_ion(lay, 0.0, El, Ue) + 0.30 * T.yield_per_ion(lay, 0.0, El / 2, Ue)


with open(os.path.join(HERE, "ue_bias.csv"), "w") as f:
    f.write("condition,Ue_true_keV,Ue_fit_uniform_keV,bias_keV\n")
    for c in CONDS:
        mat = "Pd" if c.startswith("Pd") else "Ti"
        Rm = np.array([base[(c, E)]["yield"]["UeP"] for E in EF])
        Rm = Rm / Rm[-1]

        def cost(Ue):
            Yn = np.array([nominal_Y(mat, E, Ue) for E in EF])
            return np.sum(np.log((Yn / Yn[-1]) / Rm) ** 2)
        r = minimize_scalar(cost, bounds=(0.0, 4.0), method="bounded", options=dict(xatol=0.01))
        f.write(f"{c},{UE_PAPER[c]},{r.x:.3f},{r.x - UE_PAPER[c]:+.3f}\n")
        print(c, "Ue true", UE_PAPER[c], "fit", round(r.x, 3))
