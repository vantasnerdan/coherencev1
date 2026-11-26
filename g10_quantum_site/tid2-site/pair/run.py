"""Registered TiD2 closing work and conditional G11 radial outer-wall sensitivity.

Run from the repository root: python g10_quantum_site/tid2-site/pair/run.py
The only solver is the existing G11 SitePair + solve consumer. No DFT is run.
"""
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import math
import os
import sys

sys.dont_write_bytecode = True

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"

import numpy as np
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NATIVE = ROOT / "g10_quantum_site/tid2-fluorite/native"
G11 = ROOT / "g11_cmns/collective/contact"
BULK_INPUT_SHA = "7bcd74342a1ff22599eb0f98e0f8ed55aa5a722d450c40331de4a3c7c0163966"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    spec = importlib.util.spec_from_file_location("g11_contact", G11 / "contact.py")
    contact = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(contact)
    sources = {}

    def source(path):
        sources[str(path.relative_to(ROOT))] = sha(path)
        return read(path)

    pair_source = source(NATIVE / "pair.json")
    cuts = source(NATIVE / "cuts.json")
    setting = source(NATIVE / "setting.json")
    reference = source(NATIVE / "relax/reference.json")
    zero = source(NATIVE / "zero/point.json")
    geometry = np.asarray(zero["positions_A"])
    cell = np.asarray(zero["cell_A"])
    invcell = np.linalg.inv(cell)
    register = pair_source["pair"]
    r0 = register["d0_A"]
    require(register["A_index"] == 8 and register["B_index"] == 19 and
            register["B_image_cell_translation"] == [0, 0, 0] and
            r0 == 2.211551631258434 and register["A_to_B_A"] == [r0, 0., 0.],
            "Unexpected native A/B/image/axis register")
    require(pair_source["reference"]["sha256"] == sha(NATIVE / "relax/reference.json") ==
            cuts["reference"]["sha256"] == zero["reference"]["sha256"], "Reference mismatch")
    require(reference["pair"] == register, "Reference pair mismatch")
    require(np.allclose(geometry[19]-geometry[8], [r0, 0., 0.], atol=1e-12, rtol=0),
            "Actual zero geometry differs from registered image")
    points = {"zero": zero}

    def point(name):
        if name not in points:
            p = source(NATIVE / name / "point.json")
            require(p["reference"] == zero["reference"] and p["settings"] == zero["settings"] and
                    p["atom_index"] == zero["atom_index"] and p["physical_species"] == zero["physical_species"] and
                    p["cell_A"] == zero["cell_A"], "Unmatched control " + name)
            expected = geometry.copy()
            for index, displacement in p["spec"]["moves_A"]:
                expected[index, 0] += displacement
            require(p["spec"]["direction"] == "100" and
                    np.allclose(p["positions_A"], expected, atol=2e-12, rtol=0), "Nonregistered motion " + name)
            points[name] = p
        return points[name]

    def symmetry_permutation(transformed):
        delta = (transformed[:, None, :] - geometry[None, :, :]) @ invcell
        species = np.array(zero["physical_species"])
        matches = (np.abs(delta-np.rint(delta)) < 1e-10).all(axis=2)
        matches &= species[:, None] == species[None, :]
        require((matches.sum(0) == 1).all() and (matches.sum(1) == 1).all(), "Symmetry not a species-preserving bijection")
        return matches.argmax(axis=1)

    # Full periodic geometry establishes exchange of these sites, not a label mirror.
    midpoint_perm = symmetry_permutation(geometry[8]+geometry[19]-geometry)
    require(midpoint_perm[8] == 19 and midpoint_perm[19] == 8, "Midpoint inversion does not exchange A/B")
    rotation = np.diag([-1., 1., -1.])
    site_perm = symmetry_permutation((geometry-geometry[8]) @ rotation.T + geometry[8])
    e0 = zero["energy_free_eV"]
    native_rows, exchange, signed = [], [], []
    register_error = []
    for row in pair_source["rows"]:
        amp = round(1000*row["displacement_per_site_A"])
        a, b, ab = (point(name) for name in (f"100p{amp:03d}", f"pair-b-{amp:03d}", f"pair-ab-{amp:03d}"))
        for label, name in (("0", "zero"), ("A", f"100p{amp:03d}"),
                            ("B", f"pair-b-{amp:03d}"), ("AB", f"pair-ab-{amp:03d}")):
            require(row["controls"][label]["sha256"] == sha(NATIVE / name / "point.json"), "Control hash mismatch")
        fa, fb, fab = (np.asarray(p["forces_eV_A"]) for p in (a, b, ab))
        ea, eb = a["energy_free_eV"]-e0, b["energy_free_eV"]-e0
        total = ab["energy_free_eV"]-e0
        connected = math.fsum((ab["energy_free_eV"], -a["energy_free_eV"], -b["energy_free_eV"], e0))
        da, db = -fa[8, 0], fb[19, 0]
        dtotal = -fab[8, 0]+fab[19, 0]
        dw = dtotal-da-db
        separation = np.asarray(ab["positions_A"])[19, 0]-np.asarray(ab["positions_A"])[8, 0]
        require(abs(separation-(r0-2*amp/1000)) < 2e-12, "Actual separation mismatch")
        register_error.extend([abs(connected-row["delta_energy_free_eV"]),
                               abs(dw+row["delta_compression_path_force_eV_A"])])
        native_rows.append(dict(u_per_site_A=amp/1000, R_A=separation,
            A_inward_energy_eV=ea, B_inward_energy_eV=eb, one_body_sum_eV=ea+eb,
            total_closing_work_eV=total, connected_W_eV=connected,
            dVA_du_eV_A=float(da), dVB_du_eV_A=float(db), dU_du_eV_A=float(dtotal),
            dW_du_eV_A=float(dw), dU_dR_eV_A=float(-dtotal/2), dW_dR_eV_A=float(-dw/2),
            restoring_generalized_force_eV_A=float(-dtotal),
            pair_mesh_check=False, PAW_overlap_checked=None, PAW_overlap_present=None,
            shortest_pair_PAW_caveat=amp >= 480))
        exchange.append(dict(u_A=amp/1000, B_minus_A_energy_eV=eb-ea,
            maximum_inverted_force_mismatch_eV_A=float(np.max(np.linalg.norm(fb[midpoint_perm]+fa, axis=1)))))
    require(max(register_error) < 1e-12, "Recomputed connected register mismatch")
    for amp in (120, 240, 480, 800):
        minus, plus = point(f"100m{amp:03d}"), point(f"100p{amp:03d}")
        fm, fp = (np.asarray(p["forces_eV_A"]) for p in (minus, plus))
        signed.append(dict(q_A=amp/1000,
            E_plus_eV=plus["energy_free_eV"]-e0, E_minus_eV=minus["energy_free_eV"]-e0,
            dE_plus_dq_eV_A=float(-fp[8, 0]), dE_minus_dq_eV_A=float(fm[8, 0]),
            maximum_rotated_force_mismatch_eV_A=float(np.max(np.linalg.norm(fm[site_perm]-fp@rotation.T, axis=1)))))

    def ray(names, index, direction):
        rows = [point(name) for name in names]
        q = [abs(p["spec"]["moves_A"][0][1]) for p in rows]
        en = [p["energy_free_eV"]-e0 for p in rows]
        slope = [-direction*p["forces_eV_A"][index][0] for p in rows]
        # Local curvature is inferred from the first Hermite segment, not measured at zero.
        k = 6*en[0]/q[0]**2-2*slope[0]/q[0]
        return dict(K_eV_A2=k, K4_eV_A4=0., origin_cut_A=0.,
                    dE_dq_at_zero_eV_A=0., q_A=q, E_eV=en, dE_dq_eV_A=slope)

    inward_a = ray([f"100p{amp:03d}" for amp in (120, 240, 480, 600, 800)], 8, 1)
    outward_a = ray([f"100m{amp:03d}" for amp in (120, 240, 480, 800)], 8, -1)
    inward_b = ray([f"pair-b-{amp:03d}" for amp in (120, 240, 480, 600)], 19, -1)
    for key in ("q_A", "E_eV", "dE_dq_eV_A"):
        inward_b[key].append(inward_a[key][-1])
    bulk_descriptor = source(G11 / "bulk_pair.json")
    require(sha(G11 / "bulk_pair.json") == BULK_INPUT_SHA, "Bulk descriptor binding changed")
    inner_source = source(G11 / "sources/baseline_BO_pair.json")
    require(inner_source["R_A"] == bulk_descriptor["inner"]["R_A"] and
            inner_source["V_pair_total_eV"] == bulk_descriptor["inner"]["V_total_eV"],
            "Borrowed G6 samples do not match the accepted bulk descriptor")
    descriptor = dict(schema="g11-site-pair-v1", pair_id="cubic-fluorite-tid2-D8-D19-100",
        R0_A=r0, axis=[1., 0., 0.], sites=[dict(id="D8", well="A"), dict(id="D19", well="B")],
        wells={"A":dict(inward=inward_a, outward=outward_a),
               "B":dict(inward=inward_b, outward=outward_a)},
        connected=dict(u_A=[r["u_per_site_A"] for r in native_rows],
            W_eV=[r["connected_W_eV"] for r in native_rows],
            dW_du_eV_A=[r["dW_du_eV_A"] for r in native_rows],
            flagged=[False, False, True, True], source_sha256=sha(NATIVE / "pair.json"),
            flag_meaning="Qualification marker ONLY: shortest two pair rows have no PAW-overlap check. Overlap is unknown, not asserted. All rows are primary; false treatment omits only these connected knots."),
        inner={**bulk_descriptor["inner"], "source":"Borrowed bulk PdD G6 occupied/vacant total interaction; NOT a Ti electronic cut."},
        joins_A=bulk_descriptor["joins_A"], preparation=dict(thermal_K=[77., 300.]),
        assumptions=dict(energy_zero="Native material-specific frozen reference: Ulat(R0)=0; no cross-material absolute electronic energy comparison.",
            origin="Zero slopes imposed by ideal crystal symmetry; native residual forces retained in result.",
            one_body_B="Independent B(-[100]) .12/.24/.48/.60 knots used. B inward .80 and every B outward knot use the explicitly geometry-checked midpoint inversion of A, not independent B SCFs.",
            exterior="W=0 for R>R0; distinct signed A rays and symmetry-derived B outward ray. Endpoint quadratic one-body confinement is modeled.",
            nearest_pair="W Hermite on one path only, W(0)=W'(0)=0; tangent continuation after last retained knot. No additional J, coupling matrix or lattice graph.",
            contact="Common Coulomb plus borrowed PdD electronic remainder, energy-continuous joins, no force matching or Ti G12 readout."))
    save(HERE / "pair.json", descriptor)
    pair = contact.SitePair(HERE / "pair.json")
    potentials = pair.potentials()
    bulk_pair = contact.SitePair(G11 / "bulk_pair.json")
    bulk_primary = bulk_pair.potentials()[0]
    pdd_native = []
    csvpath = ROOT / "g10_quantum_site/pdd-nearest-pair/pair-cut.csv"
    sources[str(csvpath.relative_to(ROOT))] = sha(csvpath)
    with csvpath.open() as stream:
        for row in csv.DictReader(stream):
            pdd_native.append({key:(value.lower() == "true" if value.lower() in ("true", "false") else float(value)) for key, value in row.items()})
    common = []
    for r in (1.971551631258434, 1.92782409436586, 1.73, 1.6, 1.28782409436586):
        common.append(dict(R_A=r, Ti_model_closing_work_eV=float(potentials[0].lattice(r)),
            Pd_model_closing_work_eV=float(bulk_primary.lattice(r)),
            Ti_connected_eV=float(potentials[0].w((r0-r)/2)),
            Pd_connected_eV=float(bulk_primary.w((bulk_pair.r0-r)/2))))
    reused = source(HERE / "bulk-reference.json")
    require(reused["inputs_sha256"] == bulk_pair.sha256 and reused["pair_id"] == bulk_pair.id,
            "Reused bulk spectrum not bound to the current pair")
    variants = []
    actual_consumer_error = []
    for potential in potentials:
        for row in native_rows:
            u, r = row["u_per_site_A"], row["R_A"]
            if potential.include_flagged:
                actual_consumer_error.extend([abs(float(potential.lattice(r))-row["total_closing_work_eV"]),
                    abs(float(potential.lattice(r, 1))-row["dU_dR_eV_A"]),
                    abs(float(potential.w(u))-row["connected_W_eV"])])
        result = contact.solve(potential)
        refined = contact.solve(potential, h=.0005, rmax=10., ecut=10., match_height=.25, rtol=3e-10)
        result.pop("spectrum")
        result["row_treatment"] = "all_four_primary" if potential.include_flagged else "omit_shortest_two_connected_only_tangent_from_u0.24"
        result["convergence"] = dict(refined_model=refined["model"],
            final_log10_contact_change={row["preparation"]:row["log10_contact_A_minus3"]-base["log10_contact_A_minus3"]
                for row, base in zip(refined["observables"], result["observables"])},
            refined_observables=refined["observables"])
        # Explicit contact-numerator support, not a partition-only cutoff check.
        energy = np.array(refined["spectrum"]["energy_eV"])
        logc = np.array(refined["spectrum"]["log10_contact_A_minus3"])*math.log(10)
        support = []
        for temperature in (77., 300.):
            lw = -(energy-energy[0])/(contact.KB*temperature)
            weights = np.exp(lw+logc-logsumexp(lw+logc))
            cumulative = np.cumsum(weights)
            support.append(dict(temperature_K=temperature,
                contact_99_9_percent_below_energy_eV=float(energy[np.searchsorted(cumulative, .999)]),
                contact_fraction_above_2eV=float(weights[energy>2].sum())))
        result["refined_contact_support"] = support
        variants.append(result)
        print(result["row_treatment"], potential.join,
              {row["preparation"]:round(row["log10_contact_A_minus3"], 8) for row in result["observables"]}, flush=True)
    require(max(actual_consumer_error) < 2e-12, "Actual SitePair lattice does not reproduce co-registered native rows")
    comparisons = []
    for ti in variants:
        if not ti["model"]["include_flagged"]:
            continue
        pd = next(v for v in reused["variants"] if v["model"]["include_flagged"] and v["model"]["join_A"] == ti["model"]["join_A"])
        for tirow in ti["observables"]:
            pdrow = next(r for r in pd["observables"] if r["preparation"] == tirow["preparation"])
            comparisons.append(dict(join_A=ti["model"]["join_A"], preparation=tirow["preparation"],
                Ti_log10_rho_A_minus3=tirow["log10_contact_A_minus3"],
                Pd_log10_rho_A_minus3=pdrow["log10_contact_A_minus3"],
                log10_Ti_over_Pd=tirow["log10_contact_A_minus3"]-pdrow["log10_contact_A_minus3"]))
    output = dict(schema="tid2-registered-pair-result-v1", descriptor_sha256=pair.sha256,
        source_sha256=sources, solver_path="g11_cmns/collective/contact/contact.py",
        solver_sha256=sha(G11 / "contact.py"), constants=dict(C_eV_A=contact.C, kinetic_eV_A2=contact.T, relative_mass_kg=contact.MU),
        normalization="Pure single-path spatial l=0; integral |u_R|^2 dR=1; rho(0)=|u_R'(0)|^2/(4*pi), A^-3. Canonical radial equilibrium, not driven; no point probability, occupancy or G12/foil normalization.",
        measured=dict(native_Ti_rows=native_rows, signed_A_100=signed, midpoint_B_controls=exchange,
            midpoint_species_permutation=midpoint_perm.tolist(), site_rotation_species_permutation=site_perm.tolist(),
            native_zero_A_force_eV_A=zero["forces_eV_A"][8], native_zero_B_force_eV_A=zero["forces_eV_A"][19],
            connected_register_max_error_eV_or_eV_A=max(register_error),
            finer_mesh=dict(point="111m480 and its own zero only", delta_change_eV=setting["change_delta_eV"],
                maximum_force_response_change_eV_A=setting["max_change_force_response_eV_A"], pair_checked=False),
            Pd_native_rows=pdd_native),
        modeled=dict(common_separations=common,
            common_separation_meaning="All-four Ti and all-five Pd G11 outer models evaluated at the same R. Hermite interpolation, not new native calculations; both energies referenced to their own R0.",
            near_origin_connected_curvature_d2W_du2_eV_A2=float(potentials[0].connected(0., 2)),
            curvature_meaning="Near-zero curvature of the declared origin-to-.12 Hermite interpolant only, not a measured Hessian, coupling matrix or full graph.",
            actual_consumer_native_max_error_eV_or_eV_A=max(actual_consumer_error)),
        variants=variants, primary_comparison=comparisons,
        reused_bulk=dict(path="g10_quantum_site/tid2-site/pair/bulk-reference.json", sha256=sha(HERE / "bulk-reference.json"),
            spectrum_sha256=reused["spectrum_sha256"], inputs_sha256=reused["inputs_sha256"],
            historical_refinement=reused["historical_refinement"]),
        largest_omissions=["No Ti electronic inner curve: the contact comparison borrows the same PdD G6 electronic shape and Coulomb core. It measures Ti outer-wall sensitivity, not a physical material reaction-rate ratio.",
            "Frozen constrained fully occupied cubic fluorite, one path and pure l=0 normalization; no relaxation, angular population, lattice bath, driven preparation or established foil phase.",
            "No dedicated Ti pair mesh comparison; shortest two pair rows have no PAW-overlap check. Their overlap/error is unknown. The all-row primary and connected-only omission are sensitivities, not error bounds."])
    save(HERE / "result.json", output)
    print("actual-consumer native register maximum error:", max(actual_consumer_error), flush=True)


if __name__ == "__main__":
    main()
