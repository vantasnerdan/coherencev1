"""Screening-potential analysis: correlations of measured d+d Ue with material
properties, and model numbers (Debye, Thomas-Fermi, dielectric scaling,
screening length / lattice spacing) for Pd, PdD, Ti, TiD2.

Run:  python3 -I analysis.py      (reads ue_table.csv, writes correlations.csv,
                                    model_numbers.txt)
NumPy/SciPy only.  Property tables below carry their source in a comment.
"""
import csv
import os
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(1)

# ---------------------------------------------------------------- constants
E2 = 14.399645          # e^2/(4 pi eps0) in eV*Angstrom
A0 = 0.529177           # Bohr radius, Angstrom
HA = 27.211386          # Hartree, eV
KB = 8.617333e-5        # eV/K
E2_EPS0 = 4 * np.pi * E2  # e^2/eps0 in eV*Angstrom

# ------------------------------------------------------------ property tables
# group number (IUPAC); lanthanides given 3
GROUP = dict(Be=2, Mg=2, Al=13, Sc=3, Ti=4, V=5, Cr=6, Mn=7, Fe=8, Co=9, Ni=10,
             Cu=11, Zn=12, Sr=2, Y=3, Zr=4, Nb=5, Mo=6, Ru=8, Rh=9, Pd=10,
             Ag=11, Cd=12, In=13, Sn=14, Sb=15, Ba=2, La=3, Hf=4, Ta=5, W=6,
             Re=7, Ir=9, Pt=10, Au=11, Tl=13, Pb=14, Bi=15, Lu=3, Ce=3, Pr=3,
             Nd=3, Sm=3, Eu=3, Gd=3, Tb=3, Dy=3, Ho=3, Er=3, Tm=3, Yb=3)
LANTH = {"La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
         "Yb", "Lu"}

# density g/cm3 and molar mass g/mol (standard handbook values, CRC)
DENS_M = dict(Be=(1.848, 9.012), Mg=(1.738, 24.305), Al=(2.70, 26.982),
              Sc=(2.985, 44.956), Ti=(4.506, 47.867), V=(6.11, 50.942),
              Cr=(7.19, 51.996), Mn=(7.21, 54.938), Fe=(7.874, 55.845),
              Co=(8.90, 58.933), Ni=(8.908, 58.693), Cu=(8.96, 63.546),
              Zn=(7.14, 65.38), Sr=(2.64, 87.62), Y=(4.472, 88.906),
              Zr=(6.52, 91.224), Nb=(8.57, 92.906), Mo=(10.28, 95.95),
              Ru=(12.45, 101.07), Rh=(12.41, 102.91), Pd=(12.02, 106.42),
              Ag=(10.49, 107.87), Cd=(8.65, 112.41), In=(7.31, 114.82),
              Sn=(7.287, 118.71), Sb=(6.697, 121.76), Ba=(3.51, 137.33),
              La=(6.162, 138.91), Hf=(13.31, 178.49), Ta=(16.69, 180.95),
              W=(19.25, 183.84), Re=(21.02, 186.21), Ir=(22.56, 192.22),
              Pt=(21.45, 195.08), Au=(19.30, 196.97), Tl=(11.85, 204.38),
              Pb=(11.34, 207.2), Bi=(9.78, 208.98), Lu=(9.841, 174.97),
              Ce=(6.77, 140.12), Pr=(6.77, 140.91), Nd=(7.01, 144.24),
              Sm=(7.52, 150.36), Eu=(5.264, 151.96), Gd=(7.90, 157.25),
              Tb=(8.23, 158.93), Dy=(8.54, 162.50), Ho=(8.79, 164.93),
              Er=(9.066, 167.26), Tm=(9.32, 168.93), Yb=(6.90, 173.05))

# Sommerfeld coefficient gamma, mJ mol^-1 K^-2 (Kittel, Introduction to Solid
# State Physics, Table 6.2; transcribed).  gamma/2.357 = N(EF)(1+lambda) in
# states/eV/atom.  Lanthanides omitted (4f contributions).
GAMMA = dict(Be=0.17, Mg=1.3, Al=1.35, Sc=10.7, Ti=3.35, V=9.26, Cr=1.40,
             Mn=9.20, Fe=4.98, Co=4.73, Ni=7.02, Cu=0.695, Zn=0.64, Sr=3.6,
             Y=10.2, Zr=2.80, Nb=7.79, Mo=2.0, Ru=3.3, Rh=4.9, Pd=9.42,
             Ag=0.646, Cd=0.688, In=1.69, Sn=1.78, Sb=0.11, Ba=2.7, La=10.0,
             Hf=2.16, Ta=5.9, W=1.3, Re=2.3, Ir=3.1, Pt=6.8, Au=0.729,
             Tl=1.47, Pb=2.98, Bi=0.008)

# work function, eV: midpoint of the range listed on Wikipedia "Work function"
# (CRC Handbook values), fetched 2026-11-16.
WF = dict(Ag=4.50, Al=4.16, Au=5.285, Ba=2.61, Be=4.98, Bi=4.31, Cd=4.08,
          Ce=2.9, Co=5.0, Cr=4.5, Cu=4.815, Eu=2.5, Fe=4.74, Gd=2.90, Hf=3.90,
          In=4.09, Ir=5.335, La=3.5, Lu=3.3, Mg=3.66, Mn=4.1, Mo=4.655,
          Nb=4.41, Nd=3.2, Ni=5.195, Pb=4.25, Pd=5.41, Pt=5.525, Re=4.72,
          Rh=4.98, Ru=4.71, Sb=4.625, Sc=3.5, Sm=2.7, Sn=4.42, Sr=2.59,
          Ta=4.40, Tb=3.00, Ti=4.33, Tl=3.84, V=4.3, W=4.435, Y=3.1, Yb=2.60,
          Zn=4.265, Zr=4.05)

# d-band centre relative to EF, eV, close-packed surface (Ruban, Hammer,
# Stoltze, Skriver, Norskov, J. Mol. Catal. A 115 (1997) 421; values as
# reproduced by Takigawa et al., RSC Adv. 2016).  Fe..Cu checked against
# Takigawa; the rest transcribed.
EPS_D = dict(Fe=-0.92, Co=-1.17, Ni=-1.29, Cu=-2.67, Ru=-1.41, Rh=-1.73,
             Pd=-1.83, Ag=-4.30, Ir=-2.11, Pt=-2.25, Au=-3.56)

# exothermic hydrogen solution at infinite dilution ("exothermic occluders",
# Fukai, The Metal-Hydrogen System): groups 3, 4, 5, lanthanides, Pd.
EXO = {"Sc", "Y", "La", "Ti", "Zr", "Hf", "V", "Nb", "Ta", "Pd"} | LANTH
IONIC = {"Mg", "Sr", "Ba"}      # stable ionic/saline hydrides


def omega(el):
    """atomic volume, Angstrom^3"""
    rho, m = DENS_M[el]
    return 1.66054 * m / rho


def z_total(el):
    g = GROUP[el]
    return g if g <= 12 else g - 10


def z_sp(el):
    g = GROUP[el]
    if g <= 2:
        return g
    if g <= 11:
        return 1          # one s-like electron for transition metals
    if g == 12:
        return 2
    return g - 10


def rs(el, z):
    n = z / omega(el)                       # electrons per A^3
    return (3.0 / (4 * np.pi * n)) ** (1 / 3) / A0


def dcount(el):
    g = GROUP[el]
    if 3 <= g <= 11 and el not in LANTH:
        return g - 1
    return None


# ----------------------------------------------------------------- load data
def load():
    rows = list(csv.DictReader(open(os.path.join(HERE, "ue_table.csv"))))
    return {r["id"]: r for r in rows}


METALS_EXCL = {"C", "Si", "Ge", "B", "BeO", "Al2O3", "CaO2"}


def raiola(rows, prefix):
    out = {}
    for k, r in rows.items():
        if not k.startswith(prefix + "_"):
            continue
        el = r["material"]
        if el in METALS_EXCL:
            continue
        out[el] = (float(r["Ue_eV"]), r["kind"])
    return out


def corr(x, y, nboot=2000):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    rp, pp = stats.pearsonr(x, y)
    rs_, ps = stats.spearmanr(x, y)
    # vectorised bootstrap of Spearman rho (ranks with ties averaged)
    idx = rng.integers(0, n, (nboot, n))
    rx = stats.rankdata(x[idx], axis=1)
    ry = stats.rankdata(y[idx], axis=1)
    rx -= rx.mean(axis=1, keepdims=True)
    ry -= ry.mean(axis=1, keepdims=True)
    den = np.sqrt((rx ** 2).sum(1) * (ry ** 2).sum(1))
    with np.errstate(invalid="ignore", divide="ignore"):
        bs = (rx * ry).sum(1) / den
    lo, hi = np.nanpercentile(bs, [2.5, 97.5])
    return n, rp, pp, rs_, ps, lo, hi


def main():
    rows = load()
    r04 = raiola(rows, "r04")
    r02 = raiola(rows, "r02")
    out = []
    lines = []

    def P(s=""):
        lines.append(s)
        print(s)

    num04 = {el: v for el, (v, k) in r04.items() if k == "value"}
    all04 = {el: v for el, (v, k) in r04.items()}   # upper limits at limit
    P(f"Raiola 2004 metals: {len(num04)} numeric, {len(all04)} incl. upper limits")

    tests = [
        ("gamma (DOS proxy, mJ/molK2)", lambda e: GAMMA.get(e)),
        ("log10 gamma", lambda e: np.log10(GAMMA[e]) if e in GAMMA else None),
        ("work function (eV)", lambda e: WF.get(e)),
        ("d-band centre (eV)", lambda e: EPS_D.get(e)),
        ("d count (group-1, TM only)", dcount),
        ("exothermic H solution (0/1)", lambda e: float(e in EXO)),
        ("H-absorber incl. ionic hydrides (0/1)",
         lambda e: float(e in EXO or e in IONIC)),
        ("r_s, all valence e", lambda e: rs(e, z_total(e))),
        ("r_s, s-p electrons only", lambda e: rs(e, z_sp(e))),
        ("atomic density (1/A3)", lambda e: 1.0 / omega(e)),
        ("group number", lambda e: float(GROUP[e])),
    ]
    P("\nCorrelation of Ue (Raiola 2004, RT) with material properties")
    P(f"{'property':42s} {'set':8s} {'n':>3s} {'Pearson':>8s} {'p':>8s} "
      f"{'Spearman':>8s} {'p':>8s} {'95% CI (boot)':>16s}")
    for name, f in tests:
        for label, data in (("numeric", num04), ("all", all04)):
            xs, ys = [], []
            for el, u in data.items():
                v = f(el)
                if v is None:
                    continue
                xs.append(v)
                ys.append(u)
            if len(xs) < 5 or np.ptp(xs) == 0:
                continue
            n, rp, pp, rsp, ps, lo, hi = corr(xs, ys)
            P(f"{name:42s} {label:8s} {n:3d} {rp:8.2f} {pp:8.3f} {rsp:8.2f} "
              f"{ps:8.3f}   [{lo:5.2f},{hi:5.2f}]")
            out.append(dict(dataset="Raiola2004_" + label, property=name, n=n,
                            pearson=round(rp, 3), p_pearson=round(pp, 4),
                            spearman=round(rsp, 3), p_spearman=round(ps, 4),
                            spearman_ci_lo=round(lo, 2),
                            spearman_ci_hi=round(hi, 2)))

    # Mann-Whitney: exothermic absorbers vs others (all metals incl. limits)
    a = [u for el, u in all04.items() if el in EXO]
    b = [u for el, u in all04.items() if el not in EXO]
    mw = stats.mannwhitneyu(a, b, alternative="less")
    P(f"\nExothermic absorbers (n={len(a)}): median Ue {np.median(a):.0f} eV; "
      f"others (n={len(b)}): median {np.median(b):.0f} eV; "
      f"Mann-Whitney one-sided p = {mw.pvalue:.2e}")
    a2 = [u for el, u in num04.items() if el in EXO]
    P(f"  exothermic absorbers with a numeric Ue at RT: "
      f"{sorted((el, u) for el, u in num04.items() if el in EXO)}")

    # reproducibility within one laboratory: 2002 vs 2004
    common = sorted(set(r02) & set(r04))
    x = [r02[e][0] for e in common]
    y = [r04[e][0] for e in common]
    n, rp, pp, rsp, ps, lo, hi = corr(x, y)
    P(f"\nRaiola 2002 vs 2004, same metals (n={n}): Pearson {rp:.2f} (p={pp:.3f}),"
      f" Spearman {rsp:.2f} (p={ps:.3f}) CI [{lo:.2f},{hi:.2f}]")
    ratio = [max(a_, b_) / max(min(a_, b_), 1) for a_, b_ in zip(x, y)]
    big = [(e, r02[e][0], r04[e][0]) for e, rr in zip(common, ratio) if rr > 2]
    P(f"  metals changing by >2x between the two surveys: {len(big)}/{n}: {big}")
    out.append(dict(dataset="Raiola2002_vs_2004", property="same metal, repeat",
                    n=n, pearson=round(rp, 3), p_pearson=round(pp, 4),
                    spearman=round(rsp, 3), p_spearman=round(ps, 4),
                    spearman_ci_lo=round(lo, 2), spearman_ci_hi=round(hi, 2)))

    # cross-laboratory spread for the same metal (d+d, room temperature,
    # metal or deuteride; excluding Karahadian and glow discharges)
    P("\nSame metal, different laboratories (d+d; RT; accelerator):")
    groups = {
        "Pd": ["r02_Pd", "r04_Pd", "y98_Pd", "k02_Pd", "c04_Pd", "h08_Pd",
               "cz26_Pd"],
        "Ti": ["r04_Ti", "y97_Ti", "y98_Ti", "k02_Ti", "by12_Ti", "cz26_Ti"],
        "Zr": ["r02_Zr", "r04_Zr", "c01_Zr", "c08_Zr", "k17_Zr", "cz26_Zr",
               "by12_Zr"],
        "Ta": ["r02e_Ta", "r02_Ta", "r04_Ta", "c01_Ta", "c04_Ta", "by12_Ta"],
        "Al": ["r02_Al", "r04_Al", "c01_Al"],
        "Au": ["r02_Au", "r04_Au", "y98_Au", "k02_Au"],
        "Fe": ["r02_Fe", "r04_Fe", "k02_Fe"],
        "Pt": ["r02_Pt", "r04_Pt"],
    }
    for el, ids in groups.items():
        vals = [(i, float(rows[i]["Ue_eV"]), rows[i]["kind"]) for i in ids]
        v = [max(u, 1.0) for _, u, _ in vals]
        P(f"  {el}: n={len(v)}  min {min(v):.0f}  max {max(v):.0f} eV  "
          f"max/min {max(v)/min(v):.1f}   "
          + ", ".join(f"{i}={u:.0f}{'(<=)' if k=='upper' else ''}"
                      for i, u, k in vals))

    with open(os.path.join(HERE, "correlations.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    # ------------------------------------------------------------ models
    P("\n================ MODEL NUMBERS ================")
    # Debye (Raiola/Bonomo): U = 2.09e-11 sqrt(n_eff rho_a / T) eV, SI rho_a
    rhoPd = 1e30 / omega("Pd")          # atoms per m^3
    T = 293.0
    for neff in (1.0, 6.3):
        U = 2.09e-11 * np.sqrt(neff * rhoPd / T)
        RD = 69 * np.sqrt(T / (neff * rhoPd)) * 1e10   # Angstrom
        nD = 4 / 3 * np.pi * RD ** 3 * neff / omega("Pd")
        P(f"Debye Pd 293 K n_eff={neff}: U={U:.0f} eV, R_D={RD:.4f} A, "
          f"electrons in Debye sphere={nD:.1e}")
    for Ut in (300, 800, 1700, 3000):
        neff = (Ut / 2.09e-11) ** 2 * T / rhoPd
        Tn = T * (2.09e-11 * np.sqrt(6.3 * rhoPd / T) / Ut) ** 2
        P(f"  Debye: Ue={Ut} eV needs n_eff={neff:.1f} at 293 K, or T={Tn:.0f} K "
          f"with n_eff=6.3")
    P("  Debye sensitivities: dlnU = 0.5 dln n_eff + 0.5 dln rho_a - 0.5 dln T;"
      " +1% linear strain -> -1.5%; 293->600 K -> x0.70")

    # Thomas-Fermi with DFT-like DOS (states/eV/atom or per formula unit),
    # and free-electron TF; Czerski fitted free-electron polarization
    mats = {
        # name: (Omega A^3 per metal atom, N(EF) per metal atom [est], d_MM, d_DD,
        #        z_sp, z_tot)
        "Pd": (14.72, 2.3, 2.751, None, 1, 10),
        "PdD0.7": ((4.025 ** 3) / 4, 0.4, 4.025 / np.sqrt(2), 4.025 / np.sqrt(2),
                   1.7, 10.7),
        "Ti": (17.65, 0.9, 2.92, None, 1, 4),
        "TiD2": ((4.44 ** 3) / 4, 2.0, 4.44 / np.sqrt(2), 4.44 / 2, 3, 6),
    }
    P("\nStatic screening lengths and implied Ue = e^2/lambda "
      "(N(EF) values are estimates; see models.md)")
    P(f"{'mat':7s} {'lam_TF(DOS)':>11s} {'Ue':>6s} {'lam/dMM':>8s} "
      f"{'lam/dDD':>8s} | {'lam_TF free sp':>14s} {'Ue':>5s} {'ratio':>6s} | "
      f"{'lam free all':>12s} {'Ue':>5s} | {'Czerski 250/sqrt(rs)':>20s}")
    ratios = {}
    for m, (Om, N, dmm, ddd, zsp, zt) in mats.items():
        g = N / Om                                  # states/eV/A^3
        lam = 1 / np.sqrt(E2_EPS0 * g)
        def tf_free(z):
            n = z / Om
            kF = (3 * np.pi ** 2 * n) ** (1 / 3)
            k = np.sqrt(4 * kF / (np.pi * A0))
            return 1 / k
        lsp, lall = tf_free(zsp), tf_free(zt)
        rs_all = (3 / (4 * np.pi * zt / Om)) ** (1 / 3) / A0
        ucz = 250 / np.sqrt(rs_all)
        ddd_ = ddd if ddd else dmm
        ratios[m] = (lam / dmm, lam / ddd_, lsp / dmm, lall / dmm)
        P(f"{m:7s} {lam:9.3f} A {E2/lam:6.1f} {lam/dmm:8.3f} {lam/ddd_:8.3f} | "
          f"{lsp:12.3f} A {E2/lsp:5.1f} {lsp/dmm:6.3f} | {lall:10.3f} A "
          f"{E2/lall:5.1f} | rs={rs_all:.2f} U={ucz:.0f} eV")
    for Ut in (300, 800, 1700, 3000):
        lam = E2 / Ut
        P(f"  Ue={Ut} eV requires static lambda={lam:.4f} A = {lam/2.751:.4f} d_nn(Pd);"
          f" turning point at E=Ue is the same length")
    # TF local density needed
    for Ut in (100, 1700):
        k = Ut / E2                                   # 1/A
        kF = np.pi * A0 * k ** 2 / 4
        n = kF ** 3 / (3 * np.pi ** 2)
        P(f"  TF: Ue={Ut} eV needs local electron density {n:.3g} e/A^3 "
          f"({n*A0**3:.3g} e/a0^3); Pd valence-electron density is "
          f"{10/14.72:.2f} e/A^3")
    # constant-shift validity: lambda/r_c = E/Ue
    P("  Constant-shift validity needs lambda >> r_c(E), i.e. Ue << E:")
    for E, U in ((6.5, 1.7), (2.5, 1.7), (0.25, 3.0), (0.25, 0.8)):
        P(f"    E={E} keV, Ue={U} keV: lambda/r_c = E/Ue = {E/U:.2f}")

    # dynamic (Lifschitz-Arista) factor for the beam deuteron vs v_F(sp) of Pd
    vF = 1.157e6 * (3 * np.pi ** 2 * 1 / 14.72) ** (1 / 3)   # m/s, z=1 per Pd
    P(f"\nDynamic screening factor U_dyn/U_ad (Lifschitz-Arista), v_F={vF:.2e} m/s")
    for Ecm in (0.25, 1.0, 2.5, 6.5):
        Elab = 2 * Ecm * 1e3                       # eV, deuteron on resting d
        v = np.sqrt(2 * Elab / 1875.613e6) * 2.998e8
        r = v / vF
        br = 0.5 + (1 - r ** 2) / (4 * r) * np.log(abs((1 + r) / (1 - r)))
        P(f"  E_cm={Ecm} keV  v/vF={r:.2f}  U_dyn/U_ad={np.sqrt(br):.3f}")

    # sensitivity of tunnelling to Ue at the paper energies, and the Ue shift
    # equivalent to a yield factor
    EG = 986.0
    P("\nYield sensitivity: dln(sigma)/dUe = 0.5*sqrt(EG)*E^-1.5 (keV^-1)")
    for E in (0.25, 1.0, 2.5, 6.5, 15.0):
        s = 0.5 * np.sqrt(EG) * E ** -1.5
        P(f"  E={E:5.2f} keV: {s:7.3f} /keV -> factor 2 in yield = dUe "
          f"{np.log(2)/s*1e3:6.0f} eV")

    with open(os.path.join(HERE, "model_numbers.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
