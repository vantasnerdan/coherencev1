"""Material layer: D depth profile and beam damage in the top ~100 nm of a Pd or Ti foil under a keV D beam,
with optional electrochemical loading from the back. 1D finite volume over the full 250 um foil.

SCOPE (written before the code)
Keeps
  - Source: per-energy flux from the Faraday-cup IV (1.15 D per charge over 1 cm2), D+ at E_lab and D2+ fragments
    at E_lab/2; implantation pdf, reflection and NRT vacancy production per ion from bca.py (ZBL + LS stopping).
  - Transport: flux driven by the gradient of the D chemical potential mu(x) (Darken form, lattice-gas mobility),
    with a non-ideal isotherm per host: Pd-D alpha / plateau / beta, Ti-D alpha / plateau / delta-TiD_x.
    Reduces to Fick in the dilute limit; two-phase coexistence is a near-flat mu, so hydride fronts move.
  - Traps: beam vacancies (state variable), m D per vacancy with binding E_b, relaxation to local equilibrium with
    the solute mu. Vacancy saturation from recombination with interstitials (v_sat), vacancy diffusion to the
    surface (Arrhenius), thermal/superabundant vacancy supply from the surface (equilibrium with H-enhanced
    formation), removal by erosion.
  - Vacuum surface: net desorption 2 s Gamma1 (p_eq(mu) - p_gas) (= k_r c^2 - absorption in the dilute limit),
    chamber D2 gas at 0.04 Pa, sputter erosion (moving frame, advection), chemisorbed D monolayer in equilibrium
    with the surface mu. Overlayer (native TiO2 or ALD) as a lumped layer: thickness ODE (sputter vs regrowth),
    series permeance, its own sticking, its own D content (Langmuir in mu) for the yield.
  - Back surface: reservoir at the electrochemical fugacity (fixed mu) or zero flux or sink.
  - Temperature: one foil temperature from beam power against electrolyte or edge-conduction cooling.
  - Yield: straight-ahead CSDA integral of n_D(z) sigma(E_cm) over the path with Bragg stopping, Ue fixed.
  - Knobs: overlayer, back fugacity (current), temperature, biaxial strain (V_H sigma_h shift of mu), alloying.
Drops
  - 3D structure (grain boundaries, blisters, hydride platelets, beam edges); ballistic mixing and recoil
    implantation of D; ion-induced desorption beyond sputtering and adsorbate stripping.
  - Vacancy clustering and loops; hydride coherency stress and hysteresis (one plateau per T).
  - Run-to-run history (each run starts from an annealed foil unless v_init is passed).
  - Electrochemical kinetics (the back is a fugacity, not a Volmer/Tafel/Heyrovsky model).
  - Surface chemistry beyond one sticking coefficient per surface (carbon build-up, S from the electrolyte).
  - Angular scattering in the yield integral (depth = path length).
Settling observable: near-surface (top 5 nm) D/M and its hydride/metal ratio vs BRIEF fact 6.

Units SI unless noted (energies in eV). x = D per host metal atom.
"""
import os
import sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.sparse import lil_matrix, csc_matrix
from scipy.optimize._numdiff import group_columns

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import bca  # noqa: E402
KB = 8.617333e-5           # eV/K
QE = 1.602177e-19
NA = 6.02214e23
AMU = 1.660539e-27
M_D2 = 4.028 * AMU

# ----------------------------------------------------------------------------------------------------------
# Parameters. Citations: see PARAM_SOURCES at the bottom and report. "est." = this model's estimate.
# ----------------------------------------------------------------------------------------------------------
PD = dict(
    name="Pd", nM=6.80e28, Z=46, M=106.42, Ed=34.0,   # Ed: Phys. Rev. 153, 735 (1967), 34 +- 2 eV
    # D diffusivity in alpha-Pd, Voelkl & Alefeld 1978: D0 1.7e-7 m2/s, Ea 0.206 eV (D faster than H at RT)
    D0=1.7e-7, Ea=0.206, D0_hyd=1.7e-7, Ea_hyd=0.206,
    # Pd-D isotherm: plateau from van 't Hoff (PdD), alpha_max / beta_min at 300 K, critical point (est. shape)
    dH=-35.0e3, dS=-92.0, a300=0.012, b300=0.60, xc=0.27, Tc=555.0, tha=1.0, ths=1.0,
    w=0.465,          # beta-branch interaction, eV per unit x (est.; fits x~0.7 at 1 bar, ~0.9 at ~1e5 bar)
    tilt=0.002,       # plateau tilt, eV across the gap (numerical regularisation, est.)
    # traps: H-vacancy binding 0.23 eV, up to 6 D per vacancy (Myers et al. 1992 RMP 64, 559; Besenbacher 1987)
    Eb=0.23, mtrap=6, Ef_v=1.5, Sf_v=2.0, Em_v=1.0, vsat=0.01,
    x_lmm=0.95,       # local-mixing re-emission threshold (est.)
    Eads=0.45,        # chemisorbed D relative to 1/2 D2, eV per D (Christmann 1988 Surf. Sci. Rep. 9, 1)
    ns=1.53e19,       # surface sites per m2 (Pd(111))
    s_clean=0.5,      # D2 sticking on clean Pd (Conrad, Ertl & Latta 1974 Surf. Sci. 41, 435: s0 ~ 0.6 for H2)
    sp="Pd",
    kth=71.8,         # W/m/K
    VH=1.73e-6, Emod=121e9, nu=0.39,   # partial molar volume of H (Peisl 1978), Young's modulus, Poisson ratio
    bca="Pd",
)
TI = dict(
    name="Ti", nM=5.67e28, Z=22, M=47.87, Ed=30.0,
    # alpha-Ti: Wasilewski & Kehl 1954 (H): D0 1.8e-6 m2/s, Ea 0.537 eV; isotope factor ~1 at 300 K
    # (Lu & Zhang 2013 J. Appl. Phys. 113, 193502)
    D0=1.8e-6, Ea=0.537,
    # delta-TiD_x: Ea(D) 0.60 eV (Kappesser & Wipf 1996; Wipf et al. 2000 J. Alloys Compd. 310, 190;
    # Kaess et al. 1997); D0 3e-7 m2/s is an estimate within the 1e-7..1e-6 range -> 2.5e-17 m2/s at 300 K
    D0_hyd=3.0e-7, Ea_hyd=0.60,
    # Ti-D: alpha + delta plateau from van 't Hoff, TSS in alpha, delta lower boundary (San-Martin & Manchester 1987)
    dH=-125.0e3, dS=-120.0, a300=0.001, b300=1.50, xc=None, Tc=None, tha=1.0, ths=2.0,
    # TSS ~0.1 at% at RT (Lenning, Craighead & Jaffee 1954); delta TiH1.5-1.99 (San-Martin & Manchester 1987)
    w=0.07, tilt=0.002,
    Eb=0.20, mtrap=4, Ef_v=2.0, Sf_v=2.0, Em_v=0.5, vsat=0.01,   # Eb, Ef_v est.; Em_v ~0.5 eV (low confidence)
    Eads=0.5, ns=1.4e19, s_clean=0.2,
    x_lmm=1.85,       # implanted Ti saturates as a TiD~1.8 layer at RT (6.6 keV D+, Radiat. Eff. 1980,
                      # doi 10.1080/00337578008209259); excess re-emitted (local mixing model, est. rate)
    sp="Ti",
    kth=20.0, VH=1.6e-6, Emod=116e9, nu=0.32,
    bca="TiD2",
)

# Overlayers (lumped). n: host formula units per m3; xm: max D per formula unit; mu0: Langmuir mid-point (eV vs
# 1/2 D2 at 1 Pa); perm: D permeability, D atoms /(m s Pa^0.5); s: D2 sticking on its outer surface;
# sp_Q, sp_Eth: sputter fit; Z, M: mean atom for stopping; natoms: atoms per formula unit; G0: regrowth (m/s)
OVERLAYERS = dict(
    # perm in D atoms/(m s Pa^0.5): Al2O3 chosen so a 2 nm film cuts the Pd outflux ~100x, TiO2 ~10x (est.)
    TiO2=dict(n=3.19e28, natoms=3, nmet=1, Z=10.7, M=26.6, xm=0.5, mu0=-0.3, perm=3e12, s=1e-4, sp="TiO2"),
    Al2O3=dict(n=1.80e28, natoms=5, nmet=2, Z=10.0, M=20.4, xm=0.05, mu0=0.3, perm=3e11, s=1e-5, sp="Al2O3"),
    Pt=dict(n=6.62e28, natoms=1, nmet=1, Z=78, M=195.1, xm=1.0, mu0=0.55, perm=1e17, s=0.3, sp="Pt"),
    Pd=dict(n=6.80e28, natoms=1, nmet=1, Z=46, M=106.42, xm=1.0, mu0=None, perm=1e17, s=0.5, sp="Pd"),
)

# Faraday-cup IV (SI Fig 3, read from plot): V_lab kV -> I uA
IV_V = np.array([0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5,
                 12, 12.5, 13])
IV_I = np.array([38, 68, 94, 110, 118, 123, 130, 134, 136, 138, 140, 145, 152, 156, 160, 165, 171, 176, 180, 182,
                 192, 207, 217, 222, 228, 238])
P_GAS = 2.75e-4 * 133.322   # chamber D2, Pa (BRIEF fact 1)
L_FOIL = 250e-6


def beam_current(E_lab_keV):
    return np.interp(E_lab_keV, IV_V, IV_I) * 1e-6


def beam_flux(E_lab_keV, area_m2=1e-4):
    """Deuterons per m2 per s (Eq. 1 of the paper with the IV current)."""
    return beam_current(E_lab_keV) * (0.85 + 2 * 0.15) / QE / area_m2


def gamma1(T):
    """D2 impingement rate per Pa (molecules / m2 / s / Pa)."""
    return 1.0 / np.sqrt(2 * np.pi * M_D2 * KB * QE * T)


# Sputter yields (atoms per incident D, normal incidence) at E_lab = 0.25, 0.5, 1, 2.5, 6.5, 13 keV.
# Ti: Yamamura formula with Yamamura-Tawara parameters (Q 0.54, Us 4.84 eV). Pd: RustBCA scaled by the Ti
# fit/BCA ratio. TiO2: RustBCA O+Ti halved. Al2O3: 0.6 x TiO2 (est.). Pt: est. (threshold ~150 eV).
SP_E = np.array([0.25, 0.5, 1.0, 2.5, 6.5, 13.0])
SP_TAB = dict(Ti=[0.0087, 0.013, 0.015, 0.0125, 0.0071, 0.0040],
              Pd=[0.008, 0.017, 0.023, 0.023, 0.015, 0.010],
              TiO2=[0.0185, 0.0185, 0.0157, 0.0124, 0.0065, 0.0038],
              Al2O3=[0.011, 0.011, 0.0094, 0.0074, 0.0039, 0.0023],
              Pt=[0.001, 0.004, 0.008, 0.010, 0.008, 0.006])
SP_ETH = dict(Ti=37.0, Pd=57.0, TiO2=30.0, Al2O3=35.0, Pt=150.0)


def sputter_yield(E_eV, key):
    E = np.asarray(E_eV, float) / 1e3
    y = np.exp(np.interp(np.log(np.maximum(E, 1e-6)), np.log(SP_E), np.log(SP_TAB[key])))
    Eth = SP_ETH[key] / 1e3
    low = E < SP_E[0]
    ramp = np.clip((E - Eth) / (SP_E[0] - Eth), 0, 1) ** 2
    return np.where(low, SP_TAB[key][0] * ramp, y)


# ----------------------------------------------------------------------------------------------------------
# Isotherms: mu(x) in eV relative to 1/2 D2 gas at 1 Pa, so p_eq = exp(2 mu / kT) Pa.
# ----------------------------------------------------------------------------------------------------------
class Isotherm:
    def __init__(self, mat, T, dmu=0.0, alloy=None):
        self.T, self.kT = T, KB * T
        m = dict(mat)
        if alloy:
            m.update(alloy)
        self.tha, self.ths, self.w, self.tilt = m["tha"], m["ths"], m["w"], m["tilt"]
        R = 8.314462
        p_pl = np.exp(m["dH"] / (R * T) - m["dS"] / R) * 101325.0       # Pa
        if m.get("Tc"):
            s0 = np.sqrt(max(1 - 300.0 / m["Tc"], 1e-9))
            s = np.sqrt(max(1 - T / m["Tc"], 0.0)) / s0
            # alpha_max: Arrhenius, ~2x from 298 to 353 K (0.11 eV, est. from PdH data); beta_min: mean-field
            a = min(m["a300"] * np.exp(-0.11 / KB * (1 / T - 1 / 300.0)), m["xc"] - 1e-3)
            b = max(m["xc"] + (m["b300"] - m["xc"]) * s, a + 2e-3)
        else:   # Ti: TSS Arrhenius-like growth with T (est.), delta boundary fixed
            a = m["a300"] * np.exp(-0.3 / KB * (1 / T - 1 / 300.0))
            b = m["b300"]
        self.a, self.b = a, b
        self.mupl = 0.5 * self.kT * np.log(p_pl) - dmu
        self.p_pl = p_pl

    def _la(self, x):
        return np.log(x / (self.tha - x))

    def _ls(self, x):
        return np.log(x / (self.ths - x))

    def mu(self, x):
        a, b, kT = self.a, self.b, self.kT
        x = np.clip(x, 1e-30, self.ths * (1 - 1e-12))
        xa = np.minimum(x, a)
        xb = np.maximum(x, b)
        mu_a = self.mupl + kT * (self._la(xa) - self._la(a))
        mu_p = self.mupl + self.tilt * (np.clip(x, a, b) - a) / (b - a)
        mu_b = self.mupl + self.tilt + kT * (self._ls(xb) - self._ls(b)) + self.w * (xb - b)
        return np.where(x < a, mu_a, np.where(x <= b, mu_p, mu_b))

    def mu0(self):
        """Dilute reference: mu = mu0 + kT ln(x/(1-x)) in the alpha branch."""
        return self.mupl - self.kT * self._la(self.a)

    def x_of_mu(self, mu):
        f = lambda x: float(self.mu(np.array([x]))[0]) - mu
        lo, hi = 1e-28, self.ths * (1 - 1e-10)
        if f(lo) > 0:
            return lo
        if f(hi) < 0:
            return hi
        return float(np.exp(brentq(lambda lx: f(np.exp(lx)), np.log(lo), np.log(hi), xtol=1e-12)))


def langmuir_x(mu, mu0, xm, kT):
    return xm / (1 + np.exp(np.clip((mu0 - mu) / kT, -200, 200)))


# ----------------------------------------------------------------------------------------------------------
# Implantation and damage tables from the BCA (cached)
# ----------------------------------------------------------------------------------------------------------
BCA_TARGETS = {"Pd": ([(46, 106.42, 1, 34.0)], 68.0),
               "PdD0.7": ([(46, 106.42, 1, 34.0), (1, 2.014, 0.7, 0)], 68.0 * 1.7 / 1.035),
               "Ti": ([(22, 47.87, 1, 30.0)], 56.7),
               "TiD2": ([(22, 47.87, 1, 30.0), (1, 2.014, 2, 0)], 45.3 * 3)}
BCA_E = np.array([250, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000,
                  8000, 9000, 10000, 11000, 12000, 13000], float)
BCA_FILE = os.path.join(HERE, "bca_tables.npz")


def _bca_job(args):
    import bca
    name, E = args
    comps, N = BCA_TARGETS[name]
    r = bca.run(E, comps, N, nions=4000, seed=int(E) + len(name))
    return name, E, r


def build_bca_tables(nproc=12):
    from multiprocessing import Pool
    jobs = [(n, E) for n in BCA_TARGETS for E in BCA_E]
    with Pool(nproc) as p:
        res = p.map(_bca_job, jobs)
    out = {}
    for name, E, r in res:
        k = f"{name}_{int(E)}"
        out[k + "_imp"] = r["imp"]
        out[k + "_vac"] = r["vac"]
        out[k + "_meta"] = np.array([r["RN"], r["RE"], r["Rp"], r["dRp"], r["nv_per_ion"]])
    out["z"] = res[0][2]["z"]
    np.savez_compressed(BCA_FILE, **out)


_BCA = None


def bca_profile(target, E_lab_eV):
    """Returns z (nm bin centres), implantation pdf per incident ion (per bin), vacancy count per ion per bin, R_N.
    Linear interpolation between tabulated energies after scaling depth with the range."""
    global _BCA
    if _BCA is None:
        if not os.path.exists(BCA_FILE):
            build_bca_tables()
        _BCA = dict(np.load(BCA_FILE))
    z = _BCA["z"]
    E = float(np.clip(E_lab_eV, BCA_E[0], BCA_E[-1]))
    j = np.searchsorted(BCA_E, E)
    j = min(max(j, 1), len(BCA_E) - 1)
    E0, E1 = BCA_E[j - 1], BCA_E[j]
    f = (np.log(E) - np.log(E0)) / (np.log(E1) - np.log(E0))
    out = []
    for key in ("imp", "vac"):
        a0 = _BCA[f"{target}_{int(E0)}_{key}"]
        a1 = _BCA[f"{target}_{int(E1)}_{key}"]
        out.append((1 - f) * a0 + f * a1)
    m0 = _BCA[f"{target}_{int(E0)}_meta"]
    m1 = _BCA[f"{target}_{int(E1)}_meta"]
    meta = (1 - f) * m0 + f * m1
    return z, out[0], out[1], meta


# ----------------------------------------------------------------------------------------------------------
# Thermal estimate
# ----------------------------------------------------------------------------------------------------------
def foil_temperature(E_lab_keV, mat, cooling, T_bath=300.0, h_el=5000.0, r_spot=5.64e-3, r_clamp=1.0e-2):
    """Spot-average foil temperature. cooling='electrolyte': conduction through the foil plus a convective film
    (h_el ~ 5e3 W/m2/K for a gas-evolving electrode, est.). cooling='dry': radial conduction in the foil to a
    clamp at r_clamp (uniform disc source; radiation neglected, < 1 W/m2/K)."""
    P = beam_current(E_lab_keV) * E_lab_keV * 1e3
    k, t = mat["kth"], L_FOIL
    if cooling == "electrolyte":
        q = P / (np.pi * r_spot ** 2)
        return T_bath + q * t / k + q / h_el
    return T_bath + P / (8 * np.pi * k * t) + P * np.log(r_clamp / r_spot) / (2 * np.pi * k * t)


# ----------------------------------------------------------------------------------------------------------
# Grid
# ----------------------------------------------------------------------------------------------------------
def make_grid(h0=0.25e-9, g1=1.05, hmid=3e-9, g2=1.13, L=L_FOIL, both=False):
    if both:   # refined at both faces (for back-side loading fronts)
        half = make_grid(h0, g1, hmid, g2, L / 2)
        return np.concatenate([half, L - half[::-1][1:]])
    e = [0.0]
    h = h0
    while e[-1] < L:
        e.append(e[-1] + h)
        h = h * (g1 if h < hmid else g2)
    e = np.array(e)
    e[-1] = L
    if e[-1] - e[-2] < 0.3 * (e[-2] - e[-3]):
        e = np.delete(e, -2)
    return e


# ----------------------------------------------------------------------------------------------------------
# The model
# ----------------------------------------------------------------------------------------------------------
def default_case(material="Pd", loaded=False, E_cm=2.5):
    return dict(
        material=material, loaded=loaded, E_cm=E_cm,
        t_pre=3600.0 if loaded else 0.0,        # electrochemical pre-loading before beam on (s)
        t_run=None,                               # beam time (s); default 12 h for E_cm <= 2, 2 h above
        T=None, cooling=None, T_bath=300.0, dT=0.0,
        s_front=None,                             # D2 sticking on the vacuum surface (None: material clean value)
        p_back=None, j_rel=1.0, n_j=1.0,          # back fugacity (Pa); scaled by j_rel**n_j
        j_cell=5000.0, eta_abs=0.05,              # cell current density (A/m2) and fraction entering the metal
        x_back=0.80,                              # target x at the back for the reference current (Pd)
        back_unloaded="zero",                     # 'zero' or 'sink'
        strain=0.0,                               # biaxial in-plane strain
        alloy=None,                               # dict of overrides for the host, see ALLOYS
        overlayer=None, d_ol=0.0,                 # name in OVERLAYERS and initial thickness (m)
        native_oxide=None,                        # Ti: initial native TiO2 thickness (m); None -> 3 nm for Ti
        pH2O_torr=3e-7,                           # oxidising residual gas for native-oxide regrowth
        sigma_des=3e-20,                          # beam stripping cross-section of chemisorbed D, m2 (est.)
        v_init=None,                              # initial vacancy fraction profile (callable z->v) for history
        vsat=None, Em_v=None, Eb=None,
        beam=True, flux_scale=1.0,
        no_gas=False,                             # remove the chamber D2 (uncertainty check)
        grid_both=False,                          # refine the grid at the back face as well
        grid=None,                                # dict of make_grid overrides (h0, g1, hmid, g2) for convergence
    )


ALLOYS = {
    # Pd77Ag23: no miscibility gap at 300 K, higher dilute solubility, lower x at 1 bar, D ~ 1/3 (est. from
    # Holleck 1970 J. Phys. Chem. 74, 503; Serra et al. 1998; see report)
    "PdAg23": dict(Tc=250.0, a300=0.30, b300=0.30, xc=0.30, dH=-50.0e3, dS=-92.0, w=0.75, D0=1.7e-7 / 3,
                   D0_hyd=1.7e-7 / 3),
    # Pd90Pt10: smaller gap, plateau ~10x higher, lower beta x (est. from Flanagan & Oates reviews)
    "PdPt10": dict(dH=-31.0e3, dS=-96.0, a300=0.008, b300=0.45, Tc=480.0, D0=1.7e-7 * 0.8, D0_hyd=1.7e-7 * 0.8),
}


class Model:
    def __init__(self, case):
        c = dict(default_case(case.get("material", "Pd"), case.get("loaded", False), case.get("E_cm", 2.5)))
        c.update(case)
        self.c = c
        mat = dict(PD if c["material"] == "Pd" else TI)
        if c["alloy"]:
            mat.update(ALLOYS[c["alloy"]] if isinstance(c["alloy"], str) else c["alloy"])
        for k in ("vsat", "Em_v", "Eb"):
            if c[k] is not None:
                mat[k] = c[k]
        self.mat = mat
        E_lab = 2 * c["E_cm"] * 1e3
        self.E_lab = E_lab
        if c["t_run"] is None:
            c["t_run"] = 12 * 3600.0 if c["E_cm"] <= 2.0 else 2 * 3600.0
        if c["cooling"] is None:
            c["cooling"] = "electrolyte" if c["loaded"] else "dry"
        if c["T"] is None:
            c["T"] = foil_temperature(E_lab / 1e3, mat, c["cooling"], c["T_bath"]) + c["dT"]
        T = c["T"]
        self.T, self.kT = T, KB * T
        # strain: biaxial in-plane strain -> hydrostatic stress, mu shift -V_H sigma_h (tensile raises solubility)
        sig_h = 2.0 / 3.0 * mat["Emod"] / (1 - mat["nu"]) * c["strain"]
        self.dmu_strain = mat["VH"] / NA * sig_h / QE
        self.iso = Isotherm(mat, T, dmu=self.dmu_strain)
        self.ths = mat["ths"]
        self.e = make_grid(both=c["grid_both"], **(c["grid"] or {}))
        self.h = np.diff(self.e)
        self.zc = 0.5 * (self.e[1:] + self.e[:-1])
        self.N = len(self.h)
        self.dzf = np.diff(self.zc)
        self.nM = mat["nM"]
        self.Dl = mat["D0"] * np.exp(-mat["Ea"] / self.kT)
        self.Dh = mat["D0_hyd"] * np.exp(-mat["Ea_hyd"] / self.kT)
        self.Dv = 1.5e-6 * np.exp(-mat["Em_v"] / self.kT)
        self.G1 = gamma1(T)
        self.p_gas = 0.0 if c["no_gas"] else P_GAS
        self.s_front = c["s_front"] if c["s_front"] is not None else mat["s_clean"]
        # beam
        self.F = beam_flux(E_lab / 1e3) * c["flux_scale"] if c["beam"] else 0.0
        zb, imp1, vac1, m1 = bca_profile(mat["bca"], E_lab)
        _, imp2, vac2, m2 = bca_profile(mat["bca"], E_lab / 2)
        f1, f2 = 0.85 / 1.15, 0.30 / 1.15
        imp = f1 * imp1 + f2 * imp2          # per incident deuteron, sums to 1 - R_N
        vac = f1 * vac1 + f2 * vac2
        self.RN = f1 * m1[0] + f2 * m2[0]
        self.nv_ion = vac.sum()
        self.Rp = (zb * imp).sum() / imp.sum()
        ze = np.concatenate([[0.0], zb + 0.5 * (zb[1] - zb[0])]) * 1e-9
        self.cimp = (ze, np.concatenate([[0.0], np.cumsum(imp)]))
        self.cvac = (ze, np.concatenate([[0.0], np.cumsum(vac)]))
        self.Ysp = float(sputter_yield(E_lab, mat["sp"])) * f1 + float(sputter_yield(E_lab / 2, mat["sp"])) * f2
        # overlayer
        ol = c["overlayer"]
        d0 = c["d_ol"]
        if ol is None and c["material"] == "Ti":
            ol = "TiO2"
            d0 = 3e-9 if c["native_oxide"] is None else c["native_oxide"]
        self.ol = OVERLAYERS[ol] if ol else None
        self.ol_name = ol
        self.d0 = d0 if ol else 0.0
        if self.ol:
            o = self.ol
            self.Ysp_ol = float(sputter_yield(E_lab, o["sp"])) * f1 + float(sputter_yield(E_lab / 2, o["sp"])) * f2
            native = (ol == "TiO2" and c["material"] == "Ti" and c["overlayer"] is None)
            # regrowth from residual H2O: impingement x sticking 0.5 / O density, attenuated exp(-d/0.4 nm) (est.)
            G_h2o = 3.513e22 * c["pH2O_torr"] / np.sqrt(18 * T) * 1e4 * 0.5 / (2 * o["n"])
            self.G0 = G_h2o if native else 0.0
            if ol == "Pd":
                self.iso_ol = Isotherm(PD, T)
        else:
            self.Ysp_ol, self.G0 = 0.0, 0.0
        # back side
        if c["loaded"]:
            if c["p_back"] is None:
                # reference fugacity: x_back in Pd at this T (Ti back sees the same Pd-overlayer fugacity)
                isoPd = Isotherm(PD, T)
                mu_ref = float(isoPd.mu(np.array([c["x_back"]]))[0])
                p_ref = np.exp(2 * mu_ref / self.kT)
            else:
                p_ref = c["p_back"]
            self.p_back = p_ref * c["j_rel"] ** c["n_j"]
            self.mu_back = 0.5 * self.kT * np.log(self.p_back)
            # electrochemical supply limit: J_in <= eta_abs * j / e
            self.Jmax_back = c["eta_abs"] * c["j_cell"] * c["j_rel"] / QE
        else:
            self.p_back = None
        # thermal / superabundant vacancy surface source (equilibrium with H-enhanced formation)
        self._build_sparsity()

    # ---------------- helpers
    def mu(self, x):
        return self.iso.mu(x)

    def D_of_x(self, x):
        if self.mat["name"] == "Pd":
            return np.full_like(x, self.Dl)
        a, b = self.iso.a, self.iso.b
        f = np.clip((x - a) / (b - a), 0, 1)
        return np.exp((1 - f) * np.log(self.Dl) + f * np.log(self.Dh))

    def ol_state(self, d):
        """Overlayer: thickness factor, permeance (atoms/m2/s/Pa^0.5), sticking."""
        if self.ol is None or d <= 1e-12:
            return 0.0, np.inf, self.s_front
        o = self.ol
        w = d / (d + 0.1e-9)
        return w, o["perm"] / d, o["s"] * w + self.s_front * (1 - w)

    def front_flux(self, mu0, d):
        """Net D atoms per m2 per s leaving the metal at the vacuum side (positive = out)."""
        w, Pi, s = self.ol_state(d)
        k = 2 * s * self.G1
        a = np.exp(np.clip(mu0 / self.kT, -300, 300))   # sqrt(p_eq)
        g = self.p_gas
        if not np.isfinite(Pi):
            return k * (a * a - g), a
        # k b^2 + Pi b - (Pi a + k g) = 0
        cq = Pi * a + k * g
        b = 2 * cq / (Pi + np.sqrt(Pi * Pi + 4 * k * cq))
        return Pi * (a - b), b

    def v_eq(self, mu):
        m = self.mat
        occ = np.log1p(np.exp(np.clip((mu - self.iso.mu0() + m["Eb"]) / self.kT, -200, 200)))
        return np.minimum(np.exp(m["Sf_v"] + np.clip(-m["Ef_v"] / self.kT + m["mtrap"] * occ, -300, 0)), 0.25)

    def _build_sparsity(self):
        N = self.N
        n = 3 * N + 1
        S = lil_matrix((n, n), dtype=np.int8)
        for i in range(N):
            for j in (i - 1, i, i + 1):
                if 0 <= j < N:
                    for a in range(3):
                        for b in range(3):
                            S[a * N + i, b * N + j] = 1
            S[i, 3 * N] = 1
            S[N + i, 3 * N] = 1
            S[2 * N + i, 3 * N] = 1
        for j in range(3 * N + 1):
            S[3 * N, j] = 1 if j in (0, N, 2 * N, 3 * N) else 0
        self.sparsity = S.tocsr()
        coo = S.tocoo()
        self._rows, self._cols = coo.row, coo.col
        self._groups = group_columns(self.sparsity)
        self._ng = self._groups.max() + 1
        self._floor = np.concatenate([np.full(3 * N, 1e-10), [1e-13]])

    def jac(self, t, y, beam_on=True):
        """Grouped sparse forward-difference Jacobian with fixed relative steps (avoids the adaptive-step
        overflow of scipy's num_jac on near-inert bulk cells)."""
        f0 = self.rhs(t, y, beam_on)
        hs = 1e-7 * np.maximum(np.abs(y), self._floor)
        vals = np.empty(len(self._rows))
        gc = self._groups[self._cols]
        for g in range(self._ng):
            cols = self._groups == g
            yp = y.copy()
            yp[cols] += hs[cols]
            df = self.rhs(t, yp, beam_on) - f0
            sel = gc == g
            vals[sel] = df[self._rows[sel]] / hs[self._cols[sel]]
        n = len(y)
        return csc_matrix((vals, (self._rows, self._cols)), shape=(n, n))

    # ---------------- right-hand side
    def rhs(self, t, y, beam_on=True):
        N, h, nM, kT = self.N, self.h, self.nM, self.kT
        x = np.maximum(y[:N], 1e-30)
        xt = y[N:2 * N]
        v = np.maximum(y[2 * N:3 * N], 0.0)
        d = max(y[3 * N], 0.0)
        mu = self.mu(x)
        F = self.F if beam_on else 0.0
        # erosion velocity of the metal (only once the overlayer is gone)
        w_ol = d / (d + 0.1e-9) if self.ol is not None else 0.0
        u = self.Ysp * F / nM * (1 - w_ol)
        # --- solute flux between cells (Darken, lattice-gas mobility, log-type mean)
        ths = self.ths
        lg = np.log(x / (ths - np.minimum(x, ths * (1 - 1e-12))))
        dl = np.diff(lg)
        dx = np.diff(x)
        small = np.abs(dl) < 1e-10
        xm = np.where(small, 0.5 * (x[1:] + x[:-1]) * (1 - 0.5 * (x[1:] + x[:-1]) / ths),
                      dx / np.where(small, 1.0, dl))
        Df = np.sqrt(self.D_of_x(x[1:]) * self.D_of_x(x[:-1]))
        J = -Df * nM * xm * np.diff(mu) / kT / self.dzf          # atoms/m2/s, positive = into depth
        # --- front
        Jf, _ = self.front_flux(mu[0], d)
        # --- back
        if self.p_back is not None:
            # reservoir: large conductance (fixed mu); conductance of a 0.5*h_last diffusion layer
            hb = 0.5 * h[-1]
            Db = self.D_of_x(x[-1:])[0]
            JD = Db * nM * x[-1] * (self.mu_back - mu[-1]) / kT / hb   # into the foil, fixed-mu limit
            Jb = JD * self.Jmax_back / (abs(JD) + self.Jmax_back)       # supply-limited (harmonic)
        elif self.c["back_unloaded"] == "sink":
            Jb = -self.D_of_x(x[-1:])[0] * nM * x[-1] / (0.5 * h[-1])
        else:
            Jb = 0.0
        dx_dt = np.zeros(N)
        dx_dt[:-1] -= J / (nM * h[:-1])
        dx_dt[1:] += J / (nM * h[1:])
        dx_dt[0] -= Jf / (nM * h[0])
        dx_dt[-1] += Jb / (nM * h[-1])
        # --- beam source, shifted by the overlayer thickness; overlayer share routed to cell 0
        if F > 0:
            ze, C = self.cimp
            Ce = np.interp(self.e + d, ze, C, right=C[-1])
            src = F * np.diff(Ce)
            src[0] += F * np.interp(d, ze, C)
            dx_dt += src / (nM * h)
            zev, Cv = self.cvac
            Cve = np.interp(self.e + d, zev, Cv, right=Cv[-1])
            gv = F * np.diff(Cve) / (nM * h)                    # vacancies per metal atom per s
        else:
            gv = 0.0
        # --- traps (relaxation to local equilibrium with the solute mu)
        m = self.mat
        theta_eq = 1.0 / (1.0 + np.exp(np.clip(-(mu - self.iso.mu0() + m["Eb"]) / kT, -200, 200)))
        cap = m["mtrap"] * v
        ktr = 1e-2 + 4 * np.pi * 3e-10 * self.D_of_x(x) * nM * (x + 1e-6)   # capture-limited rate (s^-1)
        ktr = np.minimum(ktr, 1e3)
        R = ktr * (cap * theta_eq - xt)
        dx_dt -= R
        dxt_dt = R.copy()
        # --- local-mixing re-emission above saturation (released D leaves the sample)
        lmm = 1.0 * np.maximum(x + xt - m["x_lmm"], 0.0)
        dx_dt -= lmm
        # --- vacancies: production with saturation, diffusion with surface at v_eq, erosion advection
        dv_dt = gv * np.clip(1 - v / m["vsat"], 0, None) if F > 0 else np.zeros(N)
        Jv = -self.Dv * np.diff(v) / self.dzf
        dv_dt[:-1] -= Jv / h[:-1]
        dv_dt[1:] += Jv / h[1:]
        veq_s = self.v_eq(mu[0])
        dv_dt[0] += self.Dv * (veq_s - v[0]) / (0.5 * h[0]) / h[0]
        # --- erosion advection (frame moves with the surface): q_i gains u q_{i+1}/h_i, loses u q_i/h_i
        if u > 0:
            for q, dq in ((x, dx_dt), (xt, dxt_dt), (v, dv_dt)):
                dq[:-1] += u * (q[1:] - q[:-1]) / h[:-1]      # upwind; cell 0 content leaves with the metal
        # --- overlayer thickness
        if self.ol is not None:
            Rsp = self.Ysp_ol * F / (self.ol["n"] * self.ol["natoms"]) * (d / (d + 0.02e-9))
            dd_dt = self.G0 * np.exp(-d / 0.4e-9) - Rsp
        else:
            dd_dt = -y[3 * N]          # no overlayer: keep the slot at zero with a nonzero Jacobian entry
        return np.concatenate([dx_dt, dxt_dt, dv_dt, [dd_dt]])

    # ---------------- run
    def initial(self):
        N = self.N
        x0 = np.full(N, 1e-9 if self.mat["name"] == "Pd" else 1e-6)
        v0 = np.zeros(N) if self.c["v_init"] is None else self.c["v_init"](self.zc)
        return np.concatenate([x0, np.zeros(N), v0, [self.d0]])

    def run(self, n_out=60):
        return self.run_te(np.linspace(0, self.c["t_run"], n_out + 1)[1:])

    def run_te(self, te):
        y0 = self.initial()
        sols = []
        t0 = 0.0
        kw = dict(method="BDF", rtol=1e-5, atol=1e-13)
        if self.c["t_pre"] > 0:
            s = solve_ivp(lambda t, y: self.rhs(t, y, beam_on=False), (0, self.c["t_pre"]), y0,
                          jac=lambda t, y: self.jac(t, y, False), **kw)
            y0 = s.y[:, -1]
            t0 = self.c["t_pre"]
            sols.append(s)
        s = solve_ivp(lambda t, y: self.rhs(t, y, beam_on=True), (0, self.c["t_run"]), y0, t_eval=te,
                          jac=lambda t, y: self.jac(t, y, True), **kw)
        if not s.success:   # retry: Radau with the same sparsity, then BDF with looser tolerances
            kw2 = dict(kw, method="Radau")
            s = solve_ivp(lambda t, y: self.rhs(t, y, beam_on=True), (0, self.c["t_run"]), y0, t_eval=te,
                          jac=lambda t, y: self.jac(t, y, True), **kw2)
        if not s.success:
            kw3 = dict(kw, rtol=3e-5, atol=1e-12)
            s = solve_ivp(lambda t, y: self.rhs(t, y, beam_on=True), (0, self.c["t_run"]), y0, t_eval=te,
                          jac=lambda t, y: self.jac(t, y, True), **kw3)
        if not s.success:
            raise RuntimeError(s.message)
        self.sol = s
        self.t_beam = s.t
        return self

    # ---------------- results
    def state(self, k=-1):
        N = self.N
        y = self.sol.y[:, k]
        return y[:N], y[N:2 * N], y[2 * N:3 * N], max(y[3 * N], 0.0)

    def surface_ads(self, x0, d):
        """Chemisorbed D per m2 on the outer surface, in equilibrium with the outer-surface mu, reduced by beam
        stripping against re-adsorption."""
        mu0 = float(self.mu(np.array([x0]))[0])
        Jf, b = self.front_flux(mu0, d)
        mu_s = self.kT * np.log(max(b, 1e-300))
        w, Pi, s = self.ol_state(d)
        Eads = self.mat["Eads"] if (self.ol is None or w < 0.5) else 0.2
        th = 1 / (1 + np.exp(np.clip(-(mu_s + Eads) / self.kT, -200, 200)))
        k_ex = 2 * s * self.G1 * (self.p_gas + b * b) / self.mat["ns"] + 1e-6
        th = th / (1 + self.c["sigma_des"] * self.F / k_ex)
        return th * self.mat["ns"], mu_s

    def profile(self, k=-1):
        """Depth profile including overlayer and adsorbate: list of layers (z_top, z_bot, n_host, n_D, Z, M)."""
        x, xt, v, d = self.state(k)
        mu = self.mu(x)
        layers = []
        nads, mu_s = self.surface_ads(x[0], d)
        if self.ol is not None and d > 1e-12:
            o = self.ol
            mu_mid = 0.5 * (mu[0] + mu_s)
            if self.ol_name == "Pd":
                xo = float(self.iso_ol.x_of_mu(mu_mid))
            else:
                xo = float(langmuir_x(mu_mid, o["mu0"], o["xm"], self.kT))
            layers.append((0.0, d, o["n"] * o["natoms"], o["n"] * xo, o["Z"], o["M"], o["n"] * o["nmet"]))
        top = d
        # host density falls with D content (lattice expansion): Pd ~ +11% volume at x 0.7; delta-TiD2 Ti density
        # 4.53e28 vs 5.67e28 m^-3 in alpha-Ti
        kexp = 0.16 if self.mat["name"] == "Pd" else 0.125
        for i in range(self.N):
            nh = self.nM / (1 + kexp * (x[i] + xt[i]))
            layers.append((top + self.e[i], top + self.e[i + 1], nh, nh * (x[i] + xt[i]),
                           self.mat["Z"], self.mat["M"], nh))
        return layers, nads, dict(x=x, xt=xt, v=v, d=d, zc=self.zc + d)

    def window_DM(self, depth, k=-1, with_ads=False):
        """Mean D per metal atom in [0, depth] (overlayer counted per its metal atom), mean D density (nm^-3),
        and mean vacancy fraction in the metal part of the window."""
        layers, nads, st = self.profile(k)
        nD = nads if with_ads else 0.0
        nH = 0.0
        for zt, zb, nh, nd, Z, M, nmet in layers:
            if zt >= depth:
                break
            w = min(zb, depth) - zt
            nD += nd * w
            nH += nmet * w
        zc, v = st["zc"], st["v"]
        e = self.e + st["d"]
        wv = np.clip(np.minimum(e[1:], depth) - e[:-1], 0, None)
        vbar = (v * wv).sum() / max(wv.sum(), 1e-30)
        return nD / nH, nD / depth * 1e-27, vbar


# ----------------------------------------------------------------------------------------------------------
# Yield: straight-ahead CSDA through the layered profile
# ----------------------------------------------------------------------------------------------------------
BG = 31.3970


def S_bh(E):   # Bosch-Hale D(d,p)T S-factor, keV mb (as in work/alr.py)
    return 5.5576e4 + E * (2.1054e2 + E * (-3.2638e-2 + E * (1.4987e-6 + E * 1.8181e-10)))


def sigma(Ecm_keV, Ue):  # paper Eq 7, mb
    E = np.maximum(Ecm_keV, 1e-6)
    return S_bh(E) / np.sqrt(E * (E + Ue)) * np.exp(-BG / np.sqrt(E + Ue))


def stopping(E_eV, nh, Z, M, nD):
    """eV per nm for D at lab energy E in a layer with host density nh (m^-3, mean Z, M) and D density nD."""
    se = bca.se_ls(E_eV, Z) + bca.sn_zbl(E_eV, Z, M)
    sd = bca.se_ls(E_eV, 1) + bca.sn_zbl(E_eV, 1, 2.014)
    return (nh * se + nD * sd) * 1e-27        # (eV nm^2)(nm^-3)


def yield_per_ion(layers, nads, E_lab_eV, Ue_keV, ds=0.05):
    """Fusion probability per incident deuteron (arbitrary consistent units: mb * D/nm^2)."""
    E = float(E_lab_eV)
    Y = nads * 1e-18 * sigma(E / 2e3, Ue_keV)          # adsorbed monolayer, D per nm^2
    z = 0.0
    for zt, zb, nh, nd, Z, M, _ in layers:
        zt9, zb9 = zt * 1e9, zb * 1e9
        nd9 = nd * 1e-27
        z = zt9
        while z < zb9 - 1e-9 and E > 0.03 * E_lab_eV:
            S = stopping(E, nh, Z, M, nd)
            st = min(zb9 - z, max(ds, 0.003 * E / S))
            Em = max(E - 0.5 * S * st, 1.0)
            Y += nd9 * sigma(Em / 2e3, Ue_keV) * st
            E -= stopping(Em, nh, Z, M, nd) * st
            z += st
        if E <= 0.03 * E_lab_eV:
            break
    return Y


def yield_per_charge(model, Ue_keV, k=-1, profile=None):
    layers, nads, _ = model.profile(k) if profile is None else profile
    E = model.E_lab
    return 0.85 * yield_per_ion(layers, nads, E, Ue_keV) + 0.30 * yield_per_ion(layers, nads, E / 2, Ue_keV)


def uniform_layers(material, x, depth=400e-9, n=200):
    mat = PD if material == "Pd" else TI
    nM = mat["nM"] / (1 + (0.16 if material == "Pd" else 0.125) * x)
    e = np.linspace(0, depth, n + 1)
    return [(e[i], e[i + 1], nM, nM * x, mat["Z"], mat["M"], nM) for i in range(n)]


PARAM_SOURCES = """see report"""
