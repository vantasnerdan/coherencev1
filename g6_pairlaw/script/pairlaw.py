"""Normalized s-wave pair-at-rest model; run from any directory (NumPy/SciPy).

Length: Angstrom. Energy: eV. Rates: per normalized pair per second.
This is a radial ground-state model, not a muonic three-body solver.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import os
import platform
import time

import numpy as np
import scipy
from scipy.constants import alpha, c, elementary_charge, hbar, physical_constants
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq, minimize_scalar
from numpy.polynomial import chebyshev as cheb

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
A0 = physical_constants["Bohr radius"][0] / 1e-10
EH = physical_constants["Hartree energy"][0] / elementary_charge
MD = physical_constants["deuteron mass"][0]
MT = physical_constants["triton mass"][0]
MUON_RATIO = 206.768262  # match the contact-source unit, not a fitted parameter
C = alpha * hbar * c / elementary_charge / 1e-10
MUS = {"DD": MD / 2, "DT": MD * MT / (MD + MT)}
T = {key: hbar**2 / (2 * mu) / elementary_charge / 1e-20
     for key, mu in MUS.items()}
LN10 = math.log(10)
TARGETS = [3e-13, 1e-9, 4e-9, 1.0]
# Ratio 7 denominator is an independently supplied physical lattice-site spacing.
LATTICE_A = 2.78


def s0_kev_b(isotope):
    if isotope == "DD":
        return (55.576, 53.701)  # Bosch-Hale branch intercepts, keV barn
    # Kamimura et al. 2112.08399v3 Fig.4 EXP fit: E in MeV, S in MeV barn.
    return (1000 * 26 / (1 + (0.0479 / 0.0392)**2),)


def coefficient(isotope, spin="allowed_statistical"):
    """A in A^3/s. Beam S already includes entrance spin/exchange weights.

    DD normalized L=0 allowed-spin mixture S=0:2 = 1:5 needs 3/4 of
    beam sigma*v/C_eta^2. DT default is unpolarized, quartet uses 3/2.
    Individual polarized DD rates are not identifiable from summed S.
    """
    s_si = sum(s0_kev_b(isotope)) * 1000 * elementary_charge * 1e-28
    spin_factor = 0.75 if isotope == "DD" else 1.0
    if isotope == "DT" and spin == "quartet":
        spin_factor = 1.5
    return spin_factor * s_si / (math.pi * alpha * MUS[isotope] * c) / 1e-30


class MolecularCurve:
    """Pachucki two-electron BO curve, with an explicit united-atom core."""
    def __init__(self, stride=1, core="quadratic"):
        data = np.genfromtxt(ROOT / "source-data-if-needed/pachucki-table4.csv",
                             delimiter=",", names=True)
        r, e = data["r_bohr"], data["E_hartree"]
        keep = np.unique(np.r_[np.arange(0, len(r), stride), len(r)-1])
        # Interpolate the smooth electronic energy, not the singular 1/r.
        self.spline = CubicSpline(r[keep], e[keep] - 1 / r[keep])
        self.core = core
        self.he = -2.9037243770341196  # infinite-mass He nonrelativistic limit
        self.core_a = ((e[0] - 1/r[0]) - self.he) / r[0]**2
        self.core_b = ((e[1] - 1/r[1]) - self.he) / r[1]**2

    def __call__(self, r_a, scale=1.0):
        x = scale * np.asarray(r_a) / A0
        # Shift dissociation to zero; this changes E, not psi/contact.
        el = self.spline(np.clip(x, .1, 20))
        a = self.core_a if self.core == "quadratic" else self.core_b
        el = np.where(x < .1, self.he + a*x*x, el)
        # Never use this long-range continuation in a consumed bound state.
        total = 1/x + el + 1
        total = np.where(x > 20, -6.499026705/x**6, total)
        return scale * EH * total


def solve(potential, isotope, rmax, n=6500, match_height=4.0):
    """Log-grid radial eigenproblem plus stable regular log-tail integration.

    u=sqrt(r)*phi, x=ln(r); -t phi''+t phi/4+r^2 V phi=E r^2 phi.
    Whiten by r, diagonalize tridiagonal matrix. The tiny inward eigenvector
    tail is NOT trusted: propagate the regular solution outward in logarithms
    and match to the normalized, numerically resolved inner-well amplitude.
    """
    t = T[isotope]
    # Restrict the well eigensolve only when a remote origin would undersample
    # a narrow well. A 14-width forbidden tail is excluded, not renormalized
    # into a different physical state; the full origin is propagated below.
    minimum = minimize_scalar(lambda z: float(potential(z)), bounds=(1e-8,rmax),
                              method="bounded", options={"xatol": 1e-12}).x
    delta = max(minimum*1e-3, 1e-8)
    curvature = (potential(minimum+delta)-2*potential(minimum)+potential(minimum-delta))/delta**2
    well_width = (2*t/float(curvature))**.25
    rlow = max(1e-9, minimum-14*well_width)
    x = np.linspace(math.log(rlow), math.log(rmax), n+2)
    r = np.exp(x)
    h = x[1]-x[0]
    v = np.asarray(potential(r), dtype=float)
    ri = r[1:-1]
    diagonal = t * (2/h**2 + .25) / ri**2 + v[1:-1]
    off = -t / h**2 / (ri[:-1]*ri[1:])
    energy, vec = eigh_tridiagonal(diagonal, off, select="i", select_range=(0, 1),
                                    tol=1e-10, lapack_driver="stebz")
    e, gap = float(energy[0]), float(energy[1]-energy[0])
    # vec has Euclidean norm 1; u=vec/sqrt(r*h) has integral norm 1.
    u = np.zeros_like(r)
    u[1:-1] = np.abs(vec[:, 0]) / np.sqrt(ri*h)
    norm = simpson(u*u*r, x=x)
    u /= math.sqrt(norm)
    mean = float(simpson(u*u*r*r, x=x))
    rms = float(math.sqrt(simpson(u*u*r**3, x=x)))
    imax = int(np.argmax(u))
    candidates = np.where((np.arange(len(r)) < imax) & (v-e > match_height*gap))[0]
    if not len(candidates):
        raise ArithmeticError("No inner forbidden-region matching point")
    j = int(candidates[-1])
    if u[j] < np.max(u)*1e-10:
        raise ArithmeticError("Matching point lies below eigenvector resolution")
    origin = math.log(1e-9)
    xmatch = float(x[j])
    def qvalue(xx):
        rr = np.exp(xx)
        return .25+rr*rr*(potential(rr)-e)/t
    # y=phi'/phi, z=ln(phi); the regular solution stays positive in this
    # forbidden region. Automatic stiff integration avoids a fixed Numerov
    # step's denominator pole when confinement creates a very high barrier.
    def ode(xx, yz):
        return [qvalue(xx)-yz[0]**2, yz[0]]
    def jac(xx, yz):
        return np.array([[-2*yz[0],0.],[1.,0.]])
    small = C*math.exp(origin)/(2*t)
    initial = [.5+small/(1+small), .5*origin+math.log1p(small)]
    tail = solve_ivp(ode,(origin,xmatch),initial,method="LSODA",jac=jac,
                     rtol=2e-10,atol=1e-10)
    if not tail.success:
        raise ArithmeticError(tail.message)
    logphi = float(tail.y[1,-1])
    logslope = math.log(u[j]) - .5*xmatch - logphi
    logcontact = (2*logslope - math.log(4*math.pi))/LN10
    # Langer/log-coordinate WKB at the SAME matched well amplitude.
    wx = np.linspace(origin,xmatch,n)
    q = qvalue(wx)
    action = float(simpson(np.sqrt(q), x=wx))
    lw = math.log(u[j]) - .5*xmatch + .25*math.log(q[-1]/.25) - action - .5*origin
    wkbcontact = (2*lw-math.log(4*math.pi))/LN10
    lograte = math.log10(coefficient(isotope)) + logcontact
    return dict(E_eV=e, gap_eV=gap, mean_A=mean, rms_A=rms,
                log10_contact_A_minus3=logcontact, log10_rate_s=lograte,
                wkb_minus_direct_dex=wkbcontact-logcontact,
                match_A=float(r[j]), outer_tail_probability=float(simpson(u[-100:]**2*r[-100:], x=x[-100:])),
                n=n, rmax_A=rmax, well_rmin_A=rlow, tail_evaluations=int(tail.nfev))


def molecule(scale, isotope, n=6500, curve=None):
    if scale <= 0:
        raise ValueError("Positive molecular scale required")
    if curve is None:
        curve = BO
    result = solve(lambda r: curve(r, scale), isotope, 15*A0/scale, n=n)
    if result["E_eV"] >= 0:
        raise ArithmeticError("Molecular state is not bound below dissociation")
    return result


def harmonic(screen_A, center_A, stiffness, isotope, n=6500, match_height=4.0):
    if screen_A <= 0 or center_A < 0 or stiffness <= 0:
        raise ValueError("Require screening>0, harmonic centre>=0, stiffness>0")
    t = T[isotope]
    width = (2*t/stiffness)**.25
    rmax = max(center_A+12*width, 3*(C/stiffness)**(1/3), .1)
    def potential(r):
        return C*np.exp(-r/screen_A)/r + .5*stiffness*(r-center_A)**2
    out = solve(potential, isotope, rmax, n, match_height)
    out.update(screen_A=screen_A, center_A=center_A, stiffness_eV_A2=stiffness,
               ratio7=screen_A/LATTICE_A, lattice_spacing_A=LATTICE_A)
    return out


def csvout(name, rows):
    with (OUT/name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def jsonout(name, value):
    (OUT/name).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def fit_scaled(samples, isotope):
    rows = [s for s in samples if s["isotope"] == isotope]
    # Smooth dimensionless log scale, 6 coefficients, interleaved holdouts.
    x = np.array([2*math.log(s["scale"])/math.log(MUON_RATIO)-1 for s in rows])
    y = np.array([s["log10_rate_s"] for s in rows])
    coefficients = cheb.chebfit(x[::2], y[::2], 5)
    error = np.abs(cheb.chebval(x[1::2], coefficients)-y[1::2])
    return dict(coordinates="x=2*ln(m)/ln(206.768262)-1; log10(rate/s^-1)=sum c_j T_j(x)",
                coefficients=coefficients.tolist(), scale_domain=[1,MUON_RATIO],
                heldout_max_dex=float(max(error)), heldout_rms_dex=float(np.sqrt(np.mean(error**2))),
                isotope=isotope, interpretation="two-electron scaled-BO model only; never a muonic law")


def on_shell_s(energy_kev):
    # Exact literal g4_beam/forward/model.py reaction_S, proton branch.
    e = energy_kev
    sb = 5.5576e4 + e*(2.1054e2+e*(-3.2638e-2+e*(1.4987e-6+e*1.8181e-10)))
    sr = math.sqrt(2*986)*197326.9804**3*40e-6/(7*937806.5**1.5*e**2)*10
    return (sb+sr+2*math.sqrt(sb*sr/3)*math.cos(math.radians(115)))/sb


def main():
    start = time.perf_counter()
    OUT.mkdir(exist_ok=True)
    scaled = []
    for isotope in MUS:
        for scale in np.geomspace(1,MUON_RATIO,33):
            scaled.append(dict(isotope=isotope, scale=float(scale),
                               electronic_length_A=A0/scale, ratio7_electronic_scale=A0/scale/LATTICE_A,
                               **molecule(float(scale),isotope)))
    fits = {isotope: fit_scaled(scaled,isotope) for isotope in MUS}
    csvout("scaled_molecule.csv", scaled)
    jsonout("scaled_fit.json", fits)

    # Numerical uncertainty decisions: dominant log-grid, origin matching,
    # outer boundary and source interpolation, not repeated anchor fitting.
    checks = []
    for isotope, scale in [("DD",1),("DD",MUON_RATIO),("DT",1),("DT",MUON_RATIO)]:
        baseline = molecule(scale,isotope)
        fine = molecule(scale,isotope,n=11000)
        coarse_curve = molecule(scale,isotope,curve=MolecularCurve(stride=2))
        core = molecule(scale,isotope,curve=MolecularCurve(core="second_point"))
        checks.append(dict(case=f"{isotope}_m{scale:g}",
                           fine_minus_base_dex=fine["log10_rate_s"]-baseline["log10_rate_s"],
                           source_thinning_shift_dex=coarse_curve["log10_rate_s"]-baseline["log10_rate_s"],
                           core_change_dex=core["log10_rate_s"]-baseline["log10_rate_s"],
                           wkb_minus_direct_dex=baseline["wkb_minus_direct_dex"]))

    # Material-readable response map, and closed-form slices on actual inputs.
    # Centre is NOT mean separation. The grid includes stiff compact controls.
    screens = [.01,.03,.1,.4,.54]
    stiffnesses = [10.,100.,1e4,1e6]
    centers = [.0,.02,.05,.1,.25,.75,2.78]
    mapping, slices, thresholds = [], [], []
    for isotope in MUS:
        for screen in screens:
            for stiffness in stiffnesses:
                cache = {}
                def get(center):
                    key = float(center)
                    if key not in cache:
                        cache[key] = harmonic(screen,key,stiffness,isotope)
                    return cache[key]
                for center in centers:
                    mapping.append(dict(isotope=isotope, **get(center)))
                # Chebyshev in centre, fitted on 13 Lobatto nodes; 12 interleaved
                # nodes independently decide whether the closed form is usable.
                nodes = np.sort(np.cos(np.arange(13)*math.pi/12))
                domain = .75
                values = np.array([get((x+1)*domain/2)["log10_rate_s"] for x in nodes])
                cf = cheb.chebfit(nodes,values,12)
                hold = np.cos((np.arange(12)+.5)*math.pi/12)
                err = [abs(cheb.chebval(x,cf)-get((x+1)*domain/2)["log10_rate_s"]) for x in hold]
                slices.append(dict(isotope=isotope, screen_A=screen, stiffness_eV_A2=stiffness,
                                   center_domain_A=[0,domain], coefficients=cf.tolist(),
                                   coordinates="x=2*r_center/0.75-1; log10(rate/s^-1)=sum c_j T_j(x)",
                                   heldout_max_dex=float(max(err)), usable_max_0p15dex=bool(max(err)<=.15)))
                for gain in [1.,1e3,1e6]:
                    for target in TARGETS:
                        goal = math.log10(target)-math.log10(gain)
                        lo,hi=get(0),get(3.0)
                        row = dict(isotope=isotope, screen_A=screen, stiffness_eV_A2=stiffness,
                                   target_s=target, contact_gain_control=gain,
                                   interpretation="background" if gain==1 else "unidentified constant-contact-gain sensitivity, not resonance prediction",
                                   ratio7=screen/LATTICE_A, lattice_spacing_A=LATTICE_A,
                                   minimum_mean_A=lo["mean_A"], center_A=None, mean_A=None,
                                   achieved_log10_rate_s=None)
                        if lo["log10_rate_s"] < goal:
                            row["status"]="no_root_nonnegative_center_rate_too_low"
                        elif hi["log10_rate_s"] > goal:
                            row["status"]="root_beyond_3A_domain"
                        else:
                            root=brentq(lambda z:get(z)["log10_rate_s"]-goal,0,3,xtol=2e-7)
                            point=get(root)
                            row.update(status="root",center_A=root,mean_A=point["mean_A"],
                                       achieved_log10_rate_s=point["log10_rate_s"]+math.log10(gain))
                        thresholds.append(row)
    csvout("harmonic_map.csv",mapping)
    jsonout("harmonic_fit.json",slices)
    csvout("harmonic_thresholds.csv",thresholds)

    molecule_thresholds=[]
    for isotope in MUS:
        for gain in [1.,1e3,1e6]:
            for target in TARGETS:
                goal=math.log10(target)-math.log10(gain)
                root=brentq(lambda lm:molecule(math.exp(lm),isotope)["log10_rate_s"]-goal,0,math.log(MUON_RATIO),xtol=1e-8)
                scale=math.exp(root)
                point=molecule(scale,isotope)
                molecule_thresholds.append(dict(isotope=isotope,target_s=target,contact_gain_control=gain,
                                                scale=scale,mean_A=point["mean_A"],electronic_length_A=A0/scale,
                                                ratio7_electronic_scale=A0/scale/LATTICE_A,
                                                lattice_spacing_A=LATTICE_A))
    csvout("scaled_thresholds.csv",molecule_thresholds)

    for isotope,screen,center,stiffness in [("DD",.54,2.78,10.),("DD",.01,0.,1e6),
                                           ("DT",.54,2.78,10.),("DT",.03,.05,1e4),
                                           ("DD",.54,2.78,1e6)]:
        base=harmonic(screen,center,stiffness,isotope)
        fine=harmonic(screen,center,stiffness,isotope,n=11000)
        match=harmonic(screen,center,stiffness,isotope,match_height=3.)
        width=(2*T[isotope]/stiffness)**.25
        p=lambda r:C*np.exp(-r/screen)/r+.5*stiffness*(r-center)**2
        outer=solve(p,isotope,base["rmax_A"]+4*width)
        checks.append(dict(case=f"{isotope}_lambda{screen:g}_centre{center:g}_k{stiffness:g}",
                           fine_minus_base_dex=fine["log10_rate_s"]-base["log10_rate_s"],
                           source_thinning_shift_dex=None,core_change_dex=None,
                           wkb_minus_direct_dex=base["wkb_minus_direct_dex"],
                           matching_shift_dex=match["log10_rate_s"]-base["log10_rate_s"],
                           outer_extension_shift_dex=outer["log10_rate_s"]-base["log10_rate_s"]))
    jsonout("numerical_checks.json",checks)

    amu=A0/MUON_RATIO
    ddcontact=2.43871205e-6/amu**3
    dtcontact=8.870847e-7/amu**3
    end={isotope: next(s for s in scaled if s["isotope"]==isotope and s["scale"]==MUON_RATIO) for isotope in MUS}
    source_comparison=dict(
        D2_reference_rate_s=3e-64,
        D2_model_over_reference=10**(scaled[0]["log10_rate_s"]-math.log10(3e-64)),
        dd_mu_ground_L0_contact_A_minus3=ddcontact,
        dd_mu_ground_L0_contact_rate_s=coefficient("DD")*ddcontact,
        dd_scaled_BO_contact_over_true_L0=10**end["DD"]["log10_contact_A_minus3"]/ddcontact,
        dd_observed_J1v1_effective_s=3.81e8,dd_deexcitation_s=2e7,
        dd_observed_p_wave_fusion_s=3.61e8,
        dd_measured_rate_is_not_s_wave_comparison=True,
        dt_mu_ground_L0_contact_A_minus3=dtcontact,
        dt_contact_quartet_rate_s=coefficient("DT","quartet")*dtcontact,
        dt_CC_J0v0_quartet_rate_s=1.15e12,
        dt_contact_quartet_over_CC=coefficient("DT","quartet")*dtcontact/1.15e12,
        dt_scaled_BO_contact_over_true_L0=10**end["DT"]["log10_contact_A_minus3"]/dtcontact,
        dd_scaled_BO_rate_over_J1_observation_nonlike_state=10**end["DD"]["log10_rate_s"]/3.61e8,
        dt_scaled_BO_quartet_rate_over_CC=1.5*10**end["DT"]["log10_rate_s"]/1.15e12)
    jsonout("source_comparison.json",source_comparison)
    csvout("on_shell_resonance_control.csv",[dict(E_cm_keV=e,proton_S_over_BH=on_shell_s(e),
                                                    interpretation="positive-energy scattering formula ONLY")
                                              for e in [1e-6,1e-4,.001,.01,.1,1.,6.5]])
    # The changed consumer: use actual mean, retain lattice spacing, no fabricated
    # short-pair claim from a minimum of the external harmonic well.
    physical_rows=[r for r in mapping if r["screen_A"]==.54 and r["center_A"]==2.78 and r["stiffness_eV_A2"]==10]
    checkmax=max(abs(r["fine_minus_base_dex"]) for r in checks)
    evidence=dict(runtime_s=time.perf_counter()-start,python=platform.python_version(),numpy=np.__version__,
                  scipy=scipy.__version__,threads={k:os.environ.get(k) for k in ["OPENBLAS_NUM_THREADS","OMP_NUM_THREADS"]},
                  source_sha256=hashlib.sha256((ROOT/"source-data-if-needed/pachucki-table4.csv").read_bytes()).hexdigest(),
                  physical_reference=physical_rows,max_grid_shift_dex=checkmax,
                  scaled_fits=fits,material_fit_usable=sum(r["usable_max_0p15dex"] for r in slices),
                  material_fit_total=len(slices),threshold_roots=sum(r["status"]=="root" for r in thresholds),
                  threshold_no_roots=sum(r["status"]!="root" for r in thresholds),
                  resonance_bound_matrix_element="not identifiable from reaction_S positive-energy fit",
                  part_c="not run: no actual material total pair potential delivered")
    jsonout("summary.json",evidence)
    # Settled numerical/consumer acceptance, expressed in log rates.
    assert checkmax < .12, ("dominant-grid shift exceeds declared 0.12 dex",checkmax)
    assert all(abs(r["achieved_log10_rate_s"]-math.log10(r["target_s"]))<.005
               for r in thresholds if r["status"]=="root")
    assert all(r["mean_A"] > r["center_A"] for r in physical_rows)
    assert all(r["ratio7"] == .54/LATTICE_A for r in physical_rows)
    assert all(r["log10_rate_s"] < -13 for r in physical_rows)
    print(json.dumps(dict(decisive_checks="passed",summary=evidence,source_comparison=source_comparison),indent=2))


BO = MolecularCurve()
if __name__ == "__main__":
    main()
