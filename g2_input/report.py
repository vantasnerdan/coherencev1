"""Assemble the accepted pump/strain inputs and retained quartic saturation estimates."""
import hashlib
import json
from pathlib import Path
import math

from mesh_readout import load_case


def number(x):
    return "none" if x is None else f"{x:.6g}"




def strain_response(metal, mesh):
    base, source = load_case(metal, mesh)
    K0 = base["K_Gamma_eV_A2"]
    response, sources = {}, [source]
    for direction in ("hydro", "parallel", "perpendicular"):
        minus, sm = load_case(metal, mesh, case=direction+"_minus")
        plus, sp = load_case(metal, mesh, case=direction+"_plus")
        Km, Kp = minus["K_Gamma_eV_A2"], plus["K_Gamma_eV_A2"]
        slope = (Kp-Km)/0.02
        response[direction] = dict(K_minus_eV_A2=Km, K_plus_eV_A2=Kp, dK_depsilon_eV_A2=slope,
                                   dlnK_depsilon=slope/K0 if K0 > 0 else None,
                                   finite_interval_log_secant=math.log(Kp/Km)/0.02
                                   if min(Kp, Km) > 0 else None)
        sources.extend((sm, sp))
    response["hydro"]["dlnK_dlnV"] = response["hydro"]["dlnK_depsilon"]/3 if K0 > 0 else None
    return response, sources


def table(headers, rows):
    return "| "+" | ".join(headers)+" |\n|"+"|".join(["---"]+["---:"]*(len(headers)-1))+"|\n" + "\n".join(
        "| "+" | ".join(str(v) for v in row)+" |" for row in rows)


def main():
    convergence = json.loads(Path("convergence.json").read_text())
    wells = json.loads(Path("wells.json").read_text())
    cost_ledger = json.loads(Path("cost-ledger.json").read_text())["materials"]
    collective = json.loads(Path("collective-wells.json").read_text())
    sign_check, sign_source = load_case("Pd", 16, case="hydro_plus")
    inputs, sources, convergence_rows, strain_rows, quartic_rows = [], [], [], [], []
    for material in convergence["materials"]:
        metal, mesh = material["metal"], material["accepted_kmesh"]
        if mesh is None:
            raise ValueError(f"{metal}: no mesh meets the accepted convergence gate")
        base, source = load_case(metal, mesh)
        K = base["K_Gamma_eV_A2"]
        response, strain_sources = strain_response(metal, mesh)
        sources.extend(strain_sources)
        cutoff_check = None
        if material["PW600_check_required"]:
            high, source = load_case(metal, mesh, 600)
            sources.append(source)
            change = high["K_Gamma_eV_A2"]/K-1
            cutoff_check = dict(K_PW600_eV_A2=high["K_Gamma_eV_A2"],
                                quantum_PW600_meV=high["harmonic_quantum_meV"], relative_K_change=change,
                                within_5percent=abs(change) <= 0.05)
        original = next(w for w in wells["materials"] if w["material"] == ("NiH" if metal == "Ni" else "PdD"))
        quartics = {direction: original["directions"][direction] for direction in ("100", "111")}
        row = dict(material=original["material"], accepted_kmesh=mesh, cutoff_eV=500,
                   K_Gamma_eV_A2=K, harmonic_quantum_meV=base["harmonic_quantum_meV"],
                   measured_peak_meV=base["measured_peak_meV"],
                   relative_frequency_difference=base["harmonic_quantum_meV"]/base["measured_peak_meV"]-1,
                   cutoff_check=cutoff_check, strain=response,
                   saturation_estimates=quartics,
                   saturation_method="LOCAL one-H/D displacement in a64-atom PW400 Gamma-only cell; all other nuclei fixed. Rough saturation estimate, distinct from the primitive collective potential.",
                   local_reference=dict(onsite_K2_eV_A2=original["coupling"]["on_site_k2_eV_A2"],
                                        sum_neighbour_Jxx_eV_A2=original["coupling"]["sum_neighbour_pulling_J_xx_eV_A2"],
                                        same_matrix_H_K_Gamma_eV_A2=original["coupling"]["H_only_q0_stiffness_eV_A2"],
                                        method=wells["method"]),
                   cost=cost_ledger[metal])
        inputs.append(row)
        comparisons = {p["mesh"]: p["relative_K_difference"] for p in material["adjacent_comparisons"]}
        for case in material["cases"]:
            difference = comparisons.get(case["kmesh"])
            convergence_rows.append([row["material"], str(case["kmesh"])+("*" if case["kmesh"] == mesh else ""),
                                      number(case["K_Gamma_eV_A2"]), number(case["harmonic_quantum_meV"]),
                                      number(100*difference) if difference is not None else "—"])
        strain_rows.append([row["material"], *[number(response[d]["dlnK_depsilon"])
                                             for d in ("hydro", "parallel", "perpendicular")],
                            number(response["hydro"]["dlnK_dlnV"])])
        quartic_rows.append([row["material"], *[number(quartics[d]["k4_eV_A4"]) for d in ("100", "111")]])
    paths = [Path(p) for p in ("convergence.json", "wells.json", "primitive.py", "gpaw_model.py", "converge.py", "strain.py", "mesh_readout.py", "report.py", "receipts/primitive-mesh-source.py")]
    paths += [Path(p) for p in ("cost-ledger.json", "mpi-prep-cost.json", "accounting.py")]
    paths += [Path(p) for p in ("collective-wells.json", "collective.py", "collective_readout.py")]
    sources.append(sign_source)
    manifests = {s["path"]: s for s in convergence["source_files"]+sources+collective["source_files"]}
    for path in paths:
        raw = path.read_bytes()
        manifests[str(path)] = dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
    summary = dict(model_inputs=inputs, convergence=convergence, collective_potentials=collective,
                   strained_PdD_k16_sign_control_eV_A2=sign_check["K_Gamma_eV_A2"],
                   strain_definition="Centre-normalized central differences [K(+1%)-K(-1%)]/(0.02*K(0)). Hydro epsilon=delta a/a; uniaxial epsilon=epsilon_xx, H displaced parallel x or perpendicular y. Metals fixed, electrons self-consistent.",
                   source_files=list(manifests.values()),
                   total_cpu_core_minutes=sum(row["cost"]["cpu_core_minutes"] for row in inputs))
    Path("model-inputs.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
    report = "# GPAW pump and strain inputs\n\n"
    report += "PW500, PBE PAW10, Fermi0.1eV, Gamma-centred meshes, H displacement+/-0.01A; frozen metal. * marks the accepted mesh. Gate:5% K agreement with the next denser mesh, with no later contradictory pair. k20 is omitted when k12/k16 already agree. Reference INS peaks: NiH89meV; PdD37meV.\n\n"
    report += table(["Material", "k mesh", "K Gamma (eV/A2)", "harmonic quantum (meV)", "K difference to next (%)"], convergence_rows)
    report += "\n\n## Central strain response\n\n"+table(
        ["Material", "hydro dlnK/depsilon_a", "uniax parallel", "uniax perpendicular", "hydro dlnK/dlnV"], strain_rows)
    report += "\n\n"+summary["strain_definition"]+" The finite-interval log secants are separate JSON fields, not the reported local response.\n\n"
    report += "## Quartic saturation estimates\n\n"+table(["Material", "K4 100 (eV/A4)", "K4 111 (eV/A4)"], quartic_rows)
    report += "\n\nV(u)-V(0)=K2*u²/2+K4*u⁴/4. These LOCAL one-H/D PW400 Gamma-only64-atom coefficients use0.12/0.24A samples with all other nuclei fixed. They are not primitive collective coefficients.\n\n"
    matched_rows = []
    for state in collective["states"]:
        label = ("PdD" if state["metal"] == "Pd" else "NiH")+(" +1% hydro" if state["strain_case"] == "hydro_plus" else " zero")
        for direction, fit in state["directions"].items():
            matched_rows.append([label, direction, number(fit["K4_eV_A4"]),
                                 number(fit["energy_RMS_meV"]), number(fit["force_prediction_RMS_eV_A"])])
    report += "## Matched primitive collective potentials\n\n"+table(
        ["State", "Direction", "K4 (eV/A4)", "energy-fit RMS (meV)", "force-prediction RMS (eV/A)"], matched_rows)
    report += "\n\n"+collective["definition"]+" K2 is fixed by the matching small-step force response. Raw energies and projected forces through0.30A, shared zero energies and fit definitions are in collective-wells.json. "+collective["normalization"]+"\n\n"
    strained = next(s for s in collective["states"] if s["metal"] == "Pd" and s["strain_case"] == "hydro_plus")
    report += "PdD at+1% hydrostatic strain: the negative bare harmonic curvature persists at the next mesh.\n\n"+table(
        ["k mesh", "K Gamma (eV/A2)"], [[12, number(strained["K2_force_eV_A2"])], [16, number(sign_check["K_Gamma_eV_A2"])]])+"\n\n"
    for row in inputs:
        report += f"{row['material']}: accepted harmonic frequency differs from INS by{100*row['relative_frequency_difference']:.2f}%. "
        if row["cutoff_check"]:
            report += f"PW600 changes K by{100*row['cutoff_check']['relative_K_change']:.3f}%."
        else:
            report += "PW600 omitted under the20% frequency-mismatch rule."
        report += "\n\n"
    report += "Charged DFT CPU: "+", ".join(f"{r['material']} {r['cost']['cpu_core_minutes']:.3f} core-min" for r in inputs)+f"; total{summary['total_cpu_core_minutes']:.3f}. Includes all failed/stopped attempts.\n\n"
    report += "Scripts: converge.py runs one mesh/cutoff case; mesh_readout.py selects the mesh; strain.py runs one of hydro_minus/hydro_plus/parallel_minus/parallel_plus/perpendicular_minus/perpendicular_plus at that mesh; report.py consumes the completed saved outputs. All SCFs use primitive.py. GPAW_SETUP_PATH must name the matching Ni.10/Pd.10/H PAW data. Raw profile identities and source-byte hashes are retained in the JSON.\n"
    Path("README.md").write_text(report)
    print(report)


if __name__ == "__main__":
    main()
