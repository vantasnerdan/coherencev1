"""Conditional one-path l=0 D–D contact; NumPy/SciPy, eV/Å/fs.

No electronic jobs, physical foil occupancy or angular population are inferred.
Run python contact.py --pair pair.json --output spectrum.json [--refine].
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import os

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np
from scipy.constants import alpha, c, elementary_charge, hbar, physical_constants
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator
from scipy.linalg import eigh_tridiagonal
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
BULK_PAIR = HERE / "bulk_pair.json"
C = alpha * hbar * c / elementary_charge / 1e-10
MU = physical_constants["deuteron mass"][0] / 2
T = hbar**2 / (2 * MU) / elementary_charge / 1e-20
HBAR = hbar / elementary_charge / 1e-15
KB = physical_constants["Boltzmann constant in eV/K"][0]
LN10 = math.log(10)


class WellRay:
    """Single-site energy at a positive displacement along a declared ray."""
    def __init__(self, data):
        self.k, self.k4 = data["K_eV_A2"], data["K4_eV_A4"]
        self.cut = data["origin_cut_A"]
        self.slope0 = data.get("dE_dq_at_zero_eV_A", 0.)
        x = np.r_[self.cut, data["q_A"]]
        e = np.r_[self.slope0*self.cut+.5*self.k*self.cut**2+self.k4*self.cut**4/4, data["E_eV"]]
        s = np.r_[self.slope0+self.k*self.cut+self.k4*self.cut**3, data["dE_dq_eV_A"]]
        if self.cut < 0 or self.k <= 0 or len(x) < 2:
            raise ValueError("A well ray needs a nonnegative origin cut, positive K, and a nonzero knot")
        self.body = CubicHermiteSpline(x, e, s, extrapolate=False)
        self.end = (float(x[-1]), float(e[-1]), float(s[-1]))

    def __call__(self, q, derivative=0):
        x = np.clip(q, self.cut, self.end[0])
        if derivative:
            near = self.slope0+self.k*q+self.k4*q**3
        else:
            near = self.slope0*q+.5*self.k*q*q+self.k4*q**4/4
        end, value, slope = self.end
        beyond = np.maximum(q-end, 0.)
        outer = slope+self.k*beyond if derivative else value+slope*beyond+.5*self.k*beyond*beyond
        return np.where(q<self.cut, near, np.where(q>end, outer, self.body(x, derivative)))


class SitePair:
    """Two D sites on one registered equal-and-opposite displacement path."""
    def __init__(self, path):
        raw = path.read_bytes()
        data = json.loads(raw)
        if data["schema"] != "g11-site-pair-v1":
            raise ValueError("Expected g11-site-pair-v1 input")
        self.sha256 = hashlib.sha256(raw).hexdigest()
        self.id = data["pair_id"]
        self.r0 = float(data["R0_A"])
        self.axis = tuple(data["axis"])
        sites = data["sites"]
        if len(sites) != 2 or sites[0]["id"] == sites[1]["id"] or self.r0 <= 0:
            raise ValueError("Supply distinct sites and a positive rest separation")
        def rays(model):
            inward = WellRay(model["inward"])
            outward = inward if model["outward"] == "reflect" else WellRay(model["outward"])
            if inward.slope0 != -outward.slope0:
                raise ValueError("Opposite rays must have opposite slopes at zero displacement")
            return inward, outward
        models = {name:rays(model) for name, model in data["wells"].items()}
        a, b = (models[site["well"]] for site in sites)
        self.inward = (a[0], b[0])
        self.outward = (a[1], b[1])
        self.outward_curvature = (a[1].k+b[1].k)/4
        self.connected = data["connected"]
        self.inner = data["inner"]
        self.joins = tuple(data["joins_A"])
        if not self.joins or len(set(self.joins)) != len(self.joins) or any(
            not self.inner["R_A"][0] <= r <= self.inner["R_A"][-1] or r >= self.r0 for r in self.joins
        ):
            raise ValueError("Joins must be distinct, covered by the inner curve, and inside R0")
        self.preparation = data.get("preparation", {})
        self.pulse_drive = data.get("pulse_drive")
        self.g12 = data.get("readout") == "g12-pdd-0.25"

    def potentials(self):
        flag_modes = (True, False) if any(self.connected["flagged"]) else (True,)
        return [Potential(self, join, flagged) for flagged in flag_modes for join in self.joins]


class Potential:
    def __init__(self, pair, join, include_flagged):
        self.pair, self.join, self.include_flagged = pair, join, include_flagged
        w = pair.connected
        keep = np.ones(len(w["u_A"]), dtype=bool) if include_flagged else ~np.array(w["flagged"], dtype=bool)
        if len(keep) != len(w["u_A"]) or not np.any(keep):
            raise ValueError("Each requested connected-curve treatment needs at least one retained knot")
        x = np.r_[0., np.array(w["u_A"])[keep]]
        e = np.r_[0., np.array(w["W_eV"])[keep]]
        s = np.r_[0., np.array(w["dW_du_eV_A"])[keep]]
        self.connected = CubicHermiteSpline(x, e, s, extrapolate=False)
        self.w_end = (float(x[-1]), float(e[-1]), float(s[-1]))
        core = pair.inner
        self.core_r = np.array(core["R_A"])
        self.electronic = PchipInterpolator(self.core_r, np.array(core["V_total_eV"])-C/self.core_r,
                                           extrapolate=False)
        self.shift = float(self.lattice(join)-self.inner(join))


    def w(self, u, derivative=0):
        end, value, slope = self.w_end
        x = np.clip(u, 0., end)
        bridge = np.full_like(x, slope) if derivative else value+slope*(u-end)
        return np.where(u<=0, 0., np.where(u>end, bridge, self.connected(x, derivative)))

    def lattice(self, r, derivative=0):
        u = (self.pair.r0-np.asarray(r))/2
        q = np.abs(u)
        i, o = self.pair.inward, self.pair.outward
        def sum_rays(rays):
            a, b = rays
            return 2*a(q, derivative) if a is b else a(q, derivative)+b(q, derivative)
        inward = sum_rays(i)
        outward = inward if i == o else sum_rays(o)
        if derivative:
            return -.5*(np.where(u<0, -outward, inward)+self.w(u, 1))
        return np.where(u<0, outward, inward)+self.w(u)

    def inner(self, r, derivative=0):
        r = np.asarray(r)
        a = self.core_r[0]
        er, slope = float(self.electronic(a)), float(self.electronic(a, 1))
        el = self.electronic(np.clip(r, a, self.core_r[-1]), derivative)
        continuation = slope*r/a if derivative else er+slope*(r*r-a*a)/(2*a)
        el = np.where(r<a, continuation, el)
        return (-C/r**2 if derivative else C/r)+el

    def __call__(self, r):
        r = np.asarray(r)
        return np.where(r>=self.join, self.lattice(r), self.inner(r)+self.shift)

    def diagnostics(self):
        grid = np.linspace(1e-5, self.join, 12001)
        v = self(grid)
        # A nonincreasing inward barrier has no artificial inner local minimum.
        return dict(join_A=self.join, include_flagged=self.include_flagged,
                    energy_shift_eV=self.shift,
                    continuity_eV=float(self.lattice(self.join)-self.inner(self.join)-self.shift),
                    slope_lattice_eV_A=float(self.lattice(self.join, 1)),
                    slope_inner_eV_A=float(self.inner(self.join, 1)),
                    slope_mismatch_eV_A=float(self.lattice(self.join, 1)-self.inner(self.join, 1)),
                    inner_monotone_outward=bool(np.all(np.diff(v)<0)),
                    inner_minimum_eV=float(np.min(v)),
                    connected_last_anchor_A=self.w_end[0])


def contact_log(potential, energy, r, amplitude, match_height, rtol):
    """Contact from the regular Coulomb solution, not an eigensolver noise tail.

    u=sqrt(r)phi, x=ln r, y=phi'/phi, z=ln phi. Regular u/r -> 1
    fixes the origin slope; match z to normalized resolved eigen-amplitude.
    """
    candidates = np.flatnonzero((r<potential.join+1.) & (potential(r)>energy+match_height))
    if not len(candidates):
        raise ArithmeticError("No radial barrier sample for the regular contact match")
    j = int(candidates[-1])
    # The outermost barrier point must lie inside the first classical turning point.
    if abs(amplitude[j]) < np.max(np.abs(amplitude))*1e-8:
        raise ArithmeticError("Contact match below eigenvector resolution")
    x0, xm = math.log(1e-9), math.log(float(r[j]))
    small = C*math.exp(x0)/(2*T)
    def ode(x, yz):
        rr = math.exp(x)
        return [.25+rr*rr*(float(potential(rr))-energy)/T-yz[0]**2, yz[0]]
    solution = solve_ivp(ode, (x0, xm), [.5+small/(1+small), .5*x0+math.log1p(small)],
                         method="DOP853", rtol=rtol, atol=rtol*.2)
    if not solution.success:
        raise ArithmeticError(solution.message)
    logslope = math.log(abs(float(amplitude[j])))-.5*xm-float(solution.y[1, -1])
    return 2*logslope-math.log(4*math.pi), math.copysign(1., amplitude[j]), float(r[j])


def solve(potential, h=.001, rmax=9., ecut=8., match_height=.15, rtol=2e-9):
    n = round(rmax/h)-1
    h = rmax/(n+1)
    r = np.arange(1, n+1)*h
    v = potential(r)
    diagonal = 2*T/h**2+v
    off = np.full(n-1, -T/h**2)
    energy, vectors = eigh_tridiagonal(diagonal, off, select="v", select_range=(-1., ecut),
                                      tol=1e-11, lapack_driver="stebz")
    amplitudes = vectors/math.sqrt(h)
    lc, signs, matches = [], [], []
    for i, en in enumerate(energy):
        logc, sign, match = contact_log(potential, float(en), r, amplitudes[:, i], match_height, rtol)
        lc.append(logc)
        signs.append(sign)
        matches.append(match)
    lc = np.array(lc)
    # Fix every eigenfunction to positive regular origin slope, preserving coherent signs.
    amplitudes *= np.array(signs)[None, :]
    rows = []
    pinput = potential.pair.preparation
    for temperature in pinput.get("thermal_K", ()):
        lw = -(energy-energy[0])/(KB*temperature)
        logz = float(logsumexp(lw))
        terms = lw+lc
        lognumerator = float(logsumexp(terms))
        contact = (lognumerator-logz)/LN10
        contribution = np.exp(terms-lognumerator)
        peak = int(np.argmax(terms))
        rows.append(dict(preparation=f"thermal_{temperature:g}K", log10_contact_A_minus3=contact,
                         contact_peak_energy_eV=float(energy[peak]),
                         contact_last_1eV_fraction=float(np.sum(contribution[energy>ecut-1.])),
                         partition_last_1eV_fraction=float(np.exp(logsumexp(lw[energy>ecut-1.])-logz)),
                         contact_weighted_mean_energy_eV=float(contribution@energy)))
    rows.append(dict(preparation="cold_ground", log10_contact_A_minus3=float(lc[0]/LN10)))
    details = dict(ground_energy_eV=float(energy[0]))
    if pinput.get("projected_approach_A") is not None:
        concentrated, details = concentrated_pulse(potential, pinput, h, r, energy, amplitudes, lc)
        rows.append(concentrated)
    diagnostics = potential.diagnostics()
    diagnostics.update(n_grid=n, h_A=h, rmax_A=rmax, energy_cutoff_eV=ecut, levels=len(energy),
                       match_height_eV=match_height, ode_rtol=rtol,
                       ground_match_A=matches[0],
                       maximum_last_100_grid_probability=float(np.max(h*np.sum(amplitudes[-100:]**2, axis=0))))
    return dict(model=diagnostics, preparation=details, observables=rows,
                spectrum=dict(energy_eV=energy.tolist(), log10_contact_A_minus3=(lc/LN10).tolist()))


def concentrated_pulse(potential, pinput, h, r, energy, amplitudes, lc):
    delta = pinput["projected_approach_A"]
    r0 = potential.pair.r0
    work = float(potential.lattice(r0-delta)-potential.lattice(r0))
    gap = float(energy[1]-energy[0])
    pop = work/gap
    if not 0<pop<1:
        raise ArithmeticError("Two-level energy matching cannot represent this preparation")
    r01 = float(h*np.dot(amplitudes[:, 0]*r, amplitudes[:, 1]))
    csign = -math.copysign(1., r01)
    duration = pinput["window_fs"]
    sinc = float(np.sinc(gap*duration/HBAR/math.pi))
    scale = max(lc[0], lc[1])
    c0, c1 = math.exp(lc[0]-scale), math.exp(lc[1]-scale)
    cross = 2*csign*math.sqrt(pop*(1-pop)*c0*c1)*sinc
    averaged = (1-pop)*c0+pop*c1+cross
    excess = pop*(c1-c0)+cross
    if min(averaged, excess)<=0:
        raise ArithmeticError("Selected window does not have positive excess contact")
    mean0 = float(h*np.dot(amplitudes[:, 0]**2, r))
    mean1 = float(h*np.dot(amplitudes[:, 1]**2, r))
    initial_mean = (1-pop)*mean0+pop*mean1+2*csign*math.sqrt(pop*(1-pop))*r01
    row = dict(preparation=f"concentrated_2level_0_{duration/1000:g}ps",
               log10_contact_A_minus3=float((scale+math.log(averaged))/LN10))
    details = dict(ground_energy_eV=float(energy[0]), gap_eV=gap, level1_log10_contact=float(lc[1]/LN10),
                   projected_turning_work_eV=work, excited_population=pop,
                   initial_mean_R_A=initial_mean, ground_mean_R_A=mean0,
                   actual_mean_inward_shift_A=mean0-initial_mean, cross_sinc=sinc,
                   log10_excess_Xi_A_minus3_s=float((scale+math.log(excess)+math.log(duration*1e-15))/LN10),
                   concentrated_over_ground=float(averaged/c0),
                   highest_preparation_energy_eV=float(energy[1]))
    return row, details


def attach_readout(result):
    from g12_readout import event_readout, rate_readout
    for row in result["observables"]:
        row["steady_one_pair"] = rate_readout(row["log10_contact_A_minus3"])
    prep = result["preparation"]
    if "log10_excess_Xi_A_minus3_s" in prep:
        prep["excess_event"] = event_readout(prep["log10_excess_Xi_A_minus3_s"])


def consumer_proof(pair, output, output_path):
    """One settled regression at final proton/G12 columns, using frozen source rows."""
    import csv
    from scipy.optimize import brentq
    source = json.loads((HERE/"sources/pair-wells.json").read_text())
    axis = np.array(source["closing_axis"])
    ia, ib = source["target_indices"]
    projection_errors = []
    for i, row in enumerate(source["rows"]):
        fa = np.array(row["one_A"]["forces_eV_A"])
        fb = np.array(row["one_B"]["forces_eV_A"])
        fab = np.array(row["pair_AB"]["forces_eV_A"])
        slope = -float(axis@(fab[ia]-fa[ia]))+float(axis@(fab[ib]-fb[ib]))
        projection_errors.extend([
            abs(slope-pair.connected["dW_du_eV_A"][i]),
            abs(row["pair_nonadditive_energy_eV"]-pair.connected["W_eV"][i])])
        assert row["widest_tabulated_Pd_H_partial_wave_overlap"] == pair.connected["flagged"][i]
    arc = json.loads((HERE/"sources/baseline_BO_pair.json").read_text())
    assert tuple(axis) == pair.axis
    assert arc["R_A"] == pair.inner["R_A"]
    assert arc["V_pair_total_eV"] == pair.inner["V_total_eV"]
    with (HERE/"sources/part1-rates.csv").open() as stream:
        rates = {row["reading"]: row for row in csv.DictReader(stream)
                 if row["target"] == "PdD" and float(row["E_cm_keV"]) == .25}
    readout_errors, shelves = [], []
    for variant in output["variants"]:
        for row in variant["observables"]:
            for key, column in (
                ("log10_emitted_p_per_s", "conditional_stationary_emitted_per_s_per_Cdomain"),
                ("log10_emitted_p_per_true_C", "conditional_stationary_emitted_per_C_per_Cdomain")):
                expected = row["log10_contact_A_minus3"]+math.log10(float(rates["k0"][column]))
                readout_errors.append(abs(row["steady_one_pair"][key]-expected))
            for reading in ("k0", "k1"):
                expected = row["log10_contact_A_minus3"]+math.log10(float(rates[reading]["conditional_stationary_deltaR_per_Cdomain"]))
                readout_errors.append(abs(row["steady_one_pair"][f"log10_delta_R_{reading}"]-expected))
        prep = variant["preparation"]
        for key, column in (
            ("log10_emitted_p_per_s", "conditional_event_emitted_per_s_per_Xi"),
            ("log10_emitted_p_per_true_C", "conditional_event_emitted_per_C_per_Xi")):
            expected = prep["log10_excess_Xi_A_minus3_s"]+math.log10(float(rates["k0"][column]))
            readout_errors.append(abs(prep["excess_event"][key]-expected))
        for reading in ("k0", "k1"):
            expected = prep["log10_excess_Xi_A_minus3_s"]+math.log10(float(rates[reading]["conditional_event_deltaR_per_Xi"]))
            readout_errors.append(abs(prep["excess_event"][f"log10_delta_R_{reading}"]-expected))
        p = Potential(pair, variant["model"]["join_A"], variant["model"]["include_flagged"])
        radius = np.linspace(.1, p.join, 12001)
        slope = p.inner(radius, 1)
        roots = [brentq(lambda z: float(p.inner(z, 1)), radius[i], radius[i+1])
                 for i in np.flatnonzero(slope[:-1]*slope[1:]<0)]
        values = [float(p(z)) for z in roots]
        energy = np.array(variant["spectrum"]["energy_eV"])
        shelf_levels = energy[(energy>=min(values)) & (energy<=max(values))]
        assert len(shelf_levels) == 0
        hot_weights = -(energy-energy[0])/(KB*300)
        above = energy>variant["model"]["inner_minimum_eV"]
        shelves.append(dict(join_A=p.join, include_flagged=p.include_flagged,
                            extrema_R_A=roots, extrema_U_eV=values,
                            shelf_depth_eV=max(values)-min(values),
                            localization_kinetic_scale_eV=T/(roots[1]-roots[0])**2,
                            eigenvalues_in_shelf_interval=shelf_levels.tolist(),
                            thermal_300K_Z_weight_above_inner_floor=float(np.exp(logsumexp(hot_weights[above])-logsumexp(hot_weights)))))
    maximum_change = max(abs(change) for v in output["variants"]
                         for change in v.get("convergence", {}).get("final_log10_contact_change", {}).values()) if all("convergence" in v for v in output["variants"]) else None
    event_change = max(abs(v["convergence"]["final_log10_excess_Xi_change"]) for v in output["variants"]) if maximum_change is not None else None
    assert max(projection_errors)<1e-12 and max(readout_errors)<1e-12
    if maximum_change is not None:
        assert max(maximum_change, event_change)<1e-4
    return dict(schema="g11-settled-consumer-proof-v1",
                source_projection_max_error=max(projection_errors),
                adopted_one_body_source_check=json.loads((HERE/"source_proof.json").read_text()),
                final_log10_readout_vs_actual_part1_max_error=max(readout_errors),
                final_readout_assertions=len(readout_errors),
                final_log10_contact_refinement_max_change=maximum_change,
                final_log10_event_refinement_max_change=event_change,
                interpolation_shelves=shelves,
                results_sha256=hashlib.sha256(output_path.read_bytes()).hexdigest(),
                independent_scope="PairRateReadout independently bound Part1 preparation/proton/G12 semantics; radial numerics and integration are delivery-owner checks.",
                limitation="Finite conditional l=0 model comparison, not physical foil rate or rigorous continuum/electronic enclosure.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair", type=Path, default=BULK_PAIR, help="Site-pair JSON description")
    parser.add_argument("--output", type=Path, required=True, help="Write a new spectrum; archived results are not defaults")
    parser.add_argument("--refine", action="store_true", help="Compare contact with a finer grid and moved cusp match")
    parser.add_argument("--no-readout", action="store_true", help="Omit the optional G12 normalization")
    parser.add_argument("--bulk-proof-output", type=Path, help="Run the historical bulk consumer check")
    args = parser.parse_args()
    pair = SitePair(args.pair)
    g12 = pair.g12 and not args.no_readout
    if args.bulk_proof_output and (args.pair.resolve() != BULK_PAIR or not g12):
        parser.error("--bulk-proof-output requires the shipped bulk pair and its readout")
    variants = []
    for potential in pair.potentials():
        result = solve(potential)
        if g12:
            attach_readout(result)
        if args.refine:
            refined = solve(potential, h=.0005, rmax=10., ecut=10., match_height=.25, rtol=3e-10)
            if g12:
                attach_readout(refined)
            result["convergence"] = dict(
                refined_model=refined["model"],
                final_log10_contact_change={row["preparation"]: row["log10_contact_A_minus3"]-base["log10_contact_A_minus3"]
                                            for row, base in zip(refined["observables"], result["observables"])},
                refined_thermal_support=[row for row in refined["observables"] if row["preparation"].startswith("thermal_")])
            if "log10_excess_Xi_A_minus3_s" in result["preparation"]:
                k = "log10_excess_Xi_A_minus3_s"
                result["convergence"]["final_log10_excess_Xi_change"] = refined["preparation"][k]-result["preparation"][k]
        variants.append(result)
    output = dict(schema="g11-conditional-contact-v1", pair_id=pair.id, inputs_sha256=pair.sha256,
                  constants=dict(C_eV_A=C, kinetic_eV_A2=T, relative_mass_kg=MU, hbar_eV_fs=HBAR), variants=variants)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    if args.bulk_proof_output:
        proof = consumer_proof(pair, output, args.output)
        args.bulk_proof_output.write_text(json.dumps(proof, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
