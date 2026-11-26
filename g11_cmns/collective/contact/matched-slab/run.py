"""Matched frozen D26/D28 slab pair, using the unchanged G11 radial consumer."""
from copy import deepcopy
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
G11 = HERE.parent
SLAB = ROOT / "g10_quantum_site/pdd-covered-slab-pair"
SITES = ROOT / "g10_quantum_site/pdd-surface-vacancy/site-wells.json"
TI = ROOT / "g10_quantum_site/tid2-site/pair"
BOUND = {
    "g10_quantum_site/pdd-covered-slab-pair/pair.json": "e3760bfca4a7d8d78627ad248589d2252813250dd7083d2f394873c6f0ac9036",
    "g10_quantum_site/pdd-surface-vacancy/site-wells.json": "5b1d19e247bcaf546d3e38fe387b65ef7da78e4db26379688d69331075ea94f7",
    "g10_quantum_site/tid2-site/pair/result.json": "5150f3451159318d3e7ae7a46350e3561d572d6f6075c747f11aac5601ddbc46",
    "g10_quantum_site/tid2-site/pair/bulk-reference.json": "c8f0efa347d5c503bfa5c5b2d989ebd5b1982eb5e439887466f29acc0acc0b19",
    "g11_cmns/collective/contact/bulk_pair.json": "7bcd74342a1ff22599eb0f98e0f8ed55aa5a722d450c40331de4a3c7c0163966",
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_bytes())


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def main():
    spec = importlib.util.spec_from_file_location("g11_contact", G11 / "contact.py")
    contact = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(contact)
    sources = {}

    def source(path):
        key = str(path.relative_to(ROOT))
        sources[key] = digest(path)
        require(key not in BOUND or sources[key] == BOUND[key], "Changed accepted input: " + key)
        return load(path)

    slab = source(SLAB / "pair.json")
    surface = source(SITES)["surface"]
    bulk = source(G11 / "bulk_pair.json")
    core = source(G11 / "sources/baseline_BO_pair.json")
    ti = source(TI / "result.json")
    ti_pair = source(TI / "pair.json")
    pd = source(TI / "bulk-reference.json")
    old = source(G11 / "surface_steady_results.json")
    require(core["R_A"] == bulk["inner"]["R_A"] and core["V_pair_total_eV"] == bulk["inner"]["V_total_eV"], "G6 samples differ")
    require(ti_pair["inner"]["R_A"] == bulk["inner"]["R_A"] and ti_pair["inner"]["V_total_eV"] == bulk["inner"]["V_total_eV"], "Ti inner differs")
    require(ti["descriptor_sha256"] == digest(TI / "pair.json") and pd["inputs_sha256"] == digest(G11 / "bulk_pair.json"), "Unbound comparator")
    require(ti["solver_sha256"] == digest(G11 / "contact.py"), "Comparator used a different consumer")
    sources[str((G11 / "contact.py").relative_to(ROOT))] = digest(G11 / "contact.py")
    sources[str(Path(__file__).relative_to(ROOT))] = digest(Path(__file__))

    native = {}
    def point(mesh, role, amp):
        name = f"m{mesh}-{role}-{amp}.json"
        if name not in native:
            native[name] = source(SLAB / "native" / name)
        return native[name]

    zero3 = point(3, "zero", 0)
    xyz = np.array(zero3["positions_A"])
    axis = np.array(zero3["point"]["axis_A_to_B"])
    r0 = float(np.linalg.norm(xyz[28]-xyz[26]))
    require(abs(r0-slab["rest_separation_A"]) < 2e-12 and abs(r0-2.4370398315432267) < 2e-12, "Wrong rest separation")
    require(np.max(np.abs((xyz[28]-xyz[26])/r0-axis)) < 2e-12, "Wrong actual axis")
    require(zero3["point"]["original_D_indices"] == [26, 28] and zero3["point"]["B_image_cell_translation"] == [0, 0, 0], "Wrong sites/image")
    nz = float(axis[2])
    normal = {}
    for label in ("subsurface", "adsorbed"):
        normal[label] = {p["point"]["displacement_A"]:p for p in surface["directions"][label] if p["point"]["axis"] == "normal"}
    require(surface["zero"]["reference"]["sha256"] == zero3["reference"]["sha256"], "Different relaxed reference")
    for label, index in (("subsurface", 26), ("adsorbed", 28)):
        require(all(p["moving_index"] == index for p in normal[label].values()), "Normal site labels changed")
    near_k = {label:(rows[-.12]["force_parallel_eV_A"]-rows[.12]["force_parallel_eV_A"])/.24*nz*nz for label, rows in normal.items()}

    measured = {}
    descriptors = {}
    for mesh, exported in ((3, slab["triplets"]), (2, slab["coarse"]["triplets"])):
        zero = point(mesh, "zero", 0)
        require(zero["reference"] == zero3["reference"] and zero["positions_A"] == zero3["positions_A"], "Mesh reference/geometry changed")
        e0 = zero["energy_free_eV"]
        slopes0 = [-float(np.dot(zero["forces_eV_A"][26], axis)), float(np.dot(zero["forces_eV_A"][28], axis))]
        rows = []
        for row in exported:
            u = row["displacement_per_D_A"]
            amp = round(u*100)
            a, b, ab = [point(mesh, role, amp) for role in ("a", "b", "ab")]
            for p, role in ((a, "A"), (b, "B"), (ab, "AB")):
                require(p["reference"] == zero["reference"] and p["calculator"] == zero["calculator"] and
                        p["symbols"] == zero["symbols"] and p["cell_A"] == zero["cell_A"], "Unmatched control")
                require(p["point"]["role"] == role and p["point"]["original_D_indices"] == [26, 28], "Wrong native role")
                expected = xyz.copy()
                if role in ("A", "AB"):
                    expected[26] += u*axis
                if role in ("B", "AB"):
                    expected[28] -= u*axis
                require(np.max(np.abs(expected-np.array(p["positions_A"]))) < 3e-12, "Nonregistered frozen path")
            for ref in row["sources"]:
                path = SLAB / "native" / ref["file"]
                require(path.stat().st_size == ref["bytes"] and digest(path) == ref["sha256"], "Native control hash mismatch")
            va, vb, total = [p["energy_free_eV"]-e0 for p in (a, b, ab)]
            w = math.fsum((ab["energy_free_eV"], -a["energy_free_eV"], -b["energy_free_eV"], e0))
            da = -float(np.dot(a["forces_eV_A"][26], axis))
            db = float(np.dot(b["forces_eV_A"][28], axis))
            dt = -float(np.dot(ab["forces_eV_A"][26], axis))+float(np.dot(ab["forces_eV_A"][28], axis))
            dw = dt-da-db
            require(max(abs(va-row["A_eV"]), abs(vb-row["B_eV"]), abs(total-row["AB_eV"]),
                        abs(w-row["connected_eV"]), abs(dt-row["AB_dE_du_eV_A"]), abs(dw-row["connected_dE_du_eV_A"])) < 2e-12, "Native/export mismatch")
            rows.append(dict(u_A=u, R_A=r0-2*u, VA_eV=va, VB_eV=vb, U_eV=total, W_eV=w,
                             dVA_du_eV_A=da, dVB_du_eV_A=db, dU_du_eV_A=dt, dW_du_eV_A=dw,
                             PAW_overlap=row["any_tagged_Pd_H_overlap"]))
        wells = {}
        for label, index, name, sign, slope0 in (("subsurface", 26, "A", 1, slopes0[0]), ("adsorbed", 28, "B", -1, slopes0[1])):
            inward = dict(K_eV_A2=near_k[label], K4_eV_A4=0., origin_cut_A=0., dE_dq_at_zero_eV_A=slope0,
                q_A=[r["u_A"] for r in rows], E_eV=[r[f"V{name}_eV"] for r in rows],
                dE_dq_eV_A=[r[f"dV{name}_du_eV_A"] for r in rows])
            out = [normal[label][-sign*q] for q in (.12, .36)]
            outward = dict(K_eV_A2=near_k[label], K4_eV_A4=0., origin_cut_A=0., dE_dq_at_zero_eV_A=-slope0,
                q_A=[abs(p["point"]["displacement_A"])/nz for p in out],
                E_eV=[p["energy_free_eV"]-surface["zero"]["energy_free_eV"] for p in out],
                dE_dq_eV_A=[sign*nz*p["force_parallel_eV_A"] for p in out])
            wells[name] = dict(inward=inward, outward=outward)
        desc = dict(schema="g11-site-pair-v1", pair_id=f"covered-pdd-D26-D28-m{mesh}-matched",
            R0_A=r0, axis=axis.tolist(), sites=[dict(id="D26", well="A"), dict(id="D28", well="B")], wells=wells,
            connected=dict(u_A=[r["u_A"] for r in rows], W_eV=[r["W_eV"] for r in rows],
                dW_du_eV_A=[r["dW_du_eV_A"] for r in rows], flagged=[r["PAW_overlap"] for r in rows],
                flag_meaning="Native tagged D-Pd PAW overlap flags retained; all five knots are primary, no error bound."),
            inner=deepcopy(bulk["inner"]), joins_A=bulk["joins_A"], preparation=dict(thermal_K=[77., 300.]),
            assumptions=dict(origin="Own matched mesh signed actual-axis zero force, not a stationary minimum. Outward origin slope is its negative for C1.",
                outward="Same original mesh2 signed normal samples relative to their own accepted zero, q=|z|/|n_z| and slope=-F_z*dz/dq. Ignore transverse response. First segment bridges the own-mesh actual-axis zero slope to mesh2 endpoints, preserving endpoints, not native mesh3 outward data.",
                exterior="Signed endpoint quadratic one-body rays, K from original +/-0.12 normal slope secant projected by n_z^2. W=0 outward, tangent after last inward knot.",
                inner="Common bulk PdD G6 total incl Coulomb and unchanged smooth electronic core; additive energy match, no slab electronic measurement or derivative fit."))
        path = HERE / f"pair-m{mesh}.json"
        save(path, desc)
        descriptors[f"m{mesh}"] = (desc, path)
        measured[f"m{mesh}"] = dict(rows=rows, native_zero_slopes_eV_A=slopes0)

    alternate = deepcopy(descriptors["m3"][0])
    alternate["pair_id"] += "-terminal-inward-curvature"
    terminal = {}
    for name, model in alternate["wells"].items():
        ray = model["inward"]
        h = CubicHermiteSpline([0., *ray["q_A"]], [0., *ray["E_eV"]], [ray["dE_dq_at_zero_eV_A"], *ray["dE_dq_eV_A"]])
        k = float(h(ray["q_A"][-1], 2))
        require(k > 0, "Terminal inward curvature not confining")
        terminal[name] = dict(default_K_eV_A2=ray["K_eV_A2"], terminal_K_eV_A2=k)
        ray["K_eV_A2"] = k
    alternate["joins_A"] = [bulk["joins_A"][0]]
    alternate["assumptions"]["exterior"] += " Sensitivity changes ONLY inward quadratic curvature beyond u=.48 to final native Hermite curvature; W tangent and all nominal knots unchanged."
    altpath = HERE / "pair-short-terminal.json"
    save(altpath, alternate)
    descriptors["short_terminal"] = (alternate, altpath)

    variants = []
    errors = []
    for case, (desc, path) in descriptors.items():
        pair = contact.SitePair(path)
        # False modes are not run: overlap flags remain, all native rows stay primary.
        for potential in (p for p in pair.potentials() if p.include_flagged):
            for row in measured["m2" if case == "m2" else "m3"]["rows"]:
                errors.extend([abs(float(potential.lattice(row["R_A"]))-row["U_eV"]),
                    abs(-2*float(potential.lattice(row["R_A"], 1))-row["dU_du_eV_A"]),
                    abs(float(potential.w(row["u_A"]))-row["W_eV"])])
            base = contact.solve(potential, ecut=3.)
            refined = contact.solve(potential, h=.0005, rmax=10., ecut=4., match_height=.25, rtol=3e-10)
            changes = {r["preparation"]:r["log10_contact_A_minus3"]-b["log10_contact_A_minus3"]
                       for r, b in zip(refined["observables"], base["observables"])}
            refined["case"] = case
            refined["descriptor_sha256"] = pair.sha256
            refined["convergence"] = dict(base_model=base["model"], base_observables=base["observables"], final_log10_contact_change=changes)
            energy = np.array(refined["spectrum"]["energy_eV"])
            logc = np.array(refined["spectrum"]["log10_contact_A_minus3"])*math.log(10)
            support = []
            for temp in (77., 300.):
                lw = -(energy-energy[0])/(contact.KB*temp)
                terms = lw+logc
                weights = np.exp(terms-logsumexp(terms))
                below = energy < 3.
                truncated = (logsumexp(terms[below])-logsumexp(lw[below]))/math.log(10)
                full = (logsumexp(terms)-logsumexp(lw))/math.log(10)
                support.append(dict(temperature_K=temp,
                    contact_99_9_percent_below_energy_eV=float(energy[np.searchsorted(np.cumsum(weights), .999)]),
                    contact_fraction_above_2eV=float(weights[energy>2].sum()),
                    log10_contact_fraction_top_1eV=float((logsumexp(terms[energy>3.])-logsumexp(terms))/math.log(10)),
                    truncating_4_to_3eV_change_dex=float(full-truncated)))
            refined["contact_support"] = support
            refined["join_support"] = dict(u_per_D_A=(r0-potential.join)/2,
                per_D_extension_A=max(0., (r0-potential.join)/2-.48), last_native_R_A=r0-.96,
                status="interpolated A/B/W" if potential.join >= r0-.96 else "conditional continuation: quadratic A/B plus tangent W")
            variants.append(refined)
            print(case, potential.join, {r["preparation"]:r["log10_contact_A_minus3"] for r in refined["observables"]}, flush=True)
    require(max(errors) < 3e-12, "Actual G11 outer model fails native energies/derivatives")
    comparisons = []
    for v in variants:
        if v["case"] != "m3":
            continue
        join = v["model"]["join_A"]
        prv = next(p for p in pd["variants"] if p["model"]["include_flagged"] and p["model"]["join_A"] == join)
        trv = next(t for t in ti["variants"] if t["model"]["include_flagged"] and t["model"]["join_A"] == join)
        for row in v["observables"]:
            prep = row["preparation"]
            p = next(r for r in prv["observables"] if r["preparation"] == prep)["log10_contact_A_minus3"]
            t = next(r for r in trv["convergence"]["refined_observables"] if r["preparation"] == prep)["log10_contact_A_minus3"]
            s = row["log10_contact_A_minus3"]
            comparisons.append(dict(join_A=join, preparation=prep, slab_log10_rho_A_minus3=s,
                Pd_log10_rho_A_minus3=p, Ti_log10_rho_A_minus3=t, log10_slab_over_Pd=s-p, log10_slab_over_Ti=s-t))
    oldlong = next(v for v in old["cases"]["near"]["variants"] if v["join_A"] == 1.6 and v["include_flagged"])
    finer = deepcopy(slab["connected_mesh_4_control"])
    for record in [finer["mesh_4_zero"], *finer["mesh_4_triplet"]["sources"]]:
        record["file"] = str((SLAB / "native" / record["file"]).relative_to(ROOT))
    save(HERE / "result.json", dict(schema="covered-slab-matched-contact-v1",
        base_commit="6d795cc90fbf62cfca46d286447d1e46262831e4", source_sha256=sources,
        normalization="One normalized pure spatial l=0 radial pair; integral |u_R|^2 dR=1; rho=|u_R'(0)|^2/(4*pi), A^-3; canonical radial equilibrium, not lattice/foil rate.",
        geometry=dict(A="D26 subsurface +axis inward", B="D28 adsorbed -axis inward", axis=axis.tolist(), R0_A=r0),
        measured=measured, single_finer_mesh_control=finer,
        outward_proxy=dict(mesh=[2,2,1], source_zero_energy_eV=surface["zero"]["energy_free_eV"],
            near_signed_normal_K_eV_A2=near_k,
            source_normal_zero_slopes_eV_A=dict(A=float(nz*surface["zero"]["forces_eV_A"][26][2]), B=float(-nz*surface["zero"]["forces_eV_A"][28][2])),
            maximum_ignored_transverse_displacement_A=.36/nz*math.sqrt(1-nz*nz)),
        terminal_inward_sensitivity=terminal, actual_consumer_native_max_error_eV_or_eV_A=max(errors),
        constants=dict(C_eV_A=contact.C, kinetic_eV_A2=contact.T, reduced_mass_kg=contact.MU, KB_eV_K=contact.KB),
        variants=variants, comparisons=comparisons,
        old_borrowed_bulk_W_long_join=dict(cold_log10_rho=oldlong["ground"]["log10_contact_A_minus3"],
            thermal_300K_log10_rho=oldlong["thermal_300K"]["log10_contact_A_minus3"],
            meaning="Historical near-exterior provisional row reused unchanged; replacement changes inward A/B and W, not solely W."),
        limitations=["No slab short-range electronic measurement; common borrowed PdD G6 total core including Coulomb, energy shift only, derivative mismatch retained.",
            "Frozen one-path pure l=0 normalized equilibrium model, not angular population, pair occupancy, dynamics, foil emission or experimental material contact.",
            "Four shortest rows retain PAW-overlap flags; coarse/finer and continuation contrasts are not error bounds. Outward normal proxy is mixed-mesh with modeled C1 origin bridge and quadratic confinement."]))


if __name__ == "__main__":
    main()
