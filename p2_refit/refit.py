"""Replace two inputs of the accepted stationary chain; do not retune eta/Ec.

Run: python p2_refit/refit.py . NumPy/SciPy only; no DFT, PySR or SINDy fit.
The accepted p1_rank/ranking.json supplies the old powers and fixed near eta.
"""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from p0_chain.pass0 import Chain, domain
from p1_rank.sensitivity import prepared, stationary, varied, grids, summarize
from p2_refit.coarse import features


def stage(name, Cv, gap, near_eta, knobs, baseline_cases=None):
    rows = []
    for eta in (0.0, near_eta, 1.0):
        model = replace(prepared(eta), gap_eV=gap)
        if baseline_cases is not None and eta in baseline_cases:
            state = baseline_cases[eta]
            origin = "saved accepted p1_rank central result, not replayed"
        else:
            state = stationary(model, knobs, Cv)
            origin = "accepted stationary closure re-evaluated"
        final_model = replace(model, density_m3=model.density_m3 * state["availability"], T_K=state["T_K"])
        rows.append(dict(eta=eta, central=state, power_origin=origin,
                         initial_domain=domain(model), steady_domain=domain(final_model),
                         initial_coarse=features(model)))
    return dict(name=name, Cv_J_m3_K=Cv, gap_eV=gap, cases=rows)


def ranking(eta, knobs, Cv, Cv_range, gap, near_eta):
    model = replace(prepared(eta), gap_eV=gap)
    axes = grids()
    axes["Cv"] = (np.linspace(*Cv_range, 9).tolist() + [Cv], "a/b", "accepted PartA working range, not the old 100-fold placeholder")
    axes["gap_eV"] = ([gap], "a/b (pinned)", "one harmonic Gamma approximation fixed; physical uncertainty NOT established")
    rows = []
    for parameter, (values, kind, basis) in axes.items():
        if parameter == "eta":
            values = values + [0.0, near_eta]
        samples = [dict(value=v, **stationary(*varied(parameter, v, model, knobs, Cv))) for v in sorted(set(values))]
        summary = summarize(samples)
        summary.pop("sampled_values")
        rows.append(dict(parameter=parameter, kind=kind, basis=basis,
                         sampled_values=[x["value"] for x in samples],
                         positive_count=sum(x["heat_W_m3"] > 0 for x in samples),
                         zero_count=sum(x["heat_W_m3"] == 0 for x in samples),
                         **summary))
    rows.sort(key=lambda x: (-(x["positive_decades"] if x["positive_decades"] is not None else -1), x["parameter"]))
    return dict(eta=eta, ranking=rows)


def main():
    started = time.perf_counter()
    baseline = json.loads((ROOT / "p1_rank/ranking.json").read_text())
    inputs = json.loads((ROOT / "p2_refit/inputs.json").read_text())
    near_eta = baseline["near_eta"]
    knobs = Chain(**baseline["chain_knobs"])
    old_gap = baseline["cases"][0]["initial_model"]["gap_eV"]
    old_cases = {row["eta"]: row["central"] for row in baseline["cases"]}
    Cv, gap = inputs["Cv_J_m3_K"], inputs["gap_eV"]
    stages = [stage("accepted_original", baseline["Cv"], old_gap, near_eta, knobs, old_cases),
              stage("Cv_only", Cv, old_gap, near_eta, knobs),
              stage("Cv_and_harmonic_gap", Cv, gap, near_eta, knobs)]
    movement = []
    for i, old in enumerate(stages[0]["cases"]):
        p0 = old["central"]["heat_W_m3"]
        p_cv = stages[1]["cases"][i]["central"]["heat_W_m3"]
        p_both = stages[2]["cases"][i]["central"]["heat_W_m3"]
        movement.append(dict(eta=old["eta"], old_W_m3=p0, Cv_only_W_m3=p_cv,
                             Cv_and_gap_W_m3=p_both, difference_W_m3=p_both-p0,
                             relative_change=(p_both/p0-1) if p0 else None))
    ranked = [ranking(eta, knobs, Cv, inputs["Cv_range_J_m3_K"], gap, near_eta) for eta in (0.0, near_eta)]
    source_paths = ("p0_chain/pass0.py", "p1_rank/sensitivity.py", "p1_rank/ranking.json",
                    "p2_refit/inputs.json", "p2_refit/coarse.py", "p2_refit/refit.py")
    output = dict(inputs=inputs, fixed_near_eta=near_eta, fixed_Ec_eV=prepared(0).Ec_eV,
                  chain_knobs=asdict(knobs), stages=stages, movement=movement, ranking_cases=ranked,
                  source_sha256={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in source_paths},
                  limits=["Long-time passive closure, not a new transient simulation or a nuclear mechanism.",
                          "Gap row is pinned by instruction, not established physically exact; 27.071 versus measured 37meV remains a material-model discrepancy.",
                          "All-zero axes have undefined positive span, not demonstrated irrelevance. A zero span with one positive point is not insensitivity.",
                          "Exact PySR threshold is re-evaluated. Approximate C/energy/heat fits are not promoted over the existing exact equilibrium.",
                          "SINDy material points exceed its declared local ranges. No new trajectories exist to identify a correction fit; none is fabricated.",
                          "Hot fixed-cell outputs are extrapolations, not demonstrated phase-stable PdD predictions.",
                          "Domain geometry and finite deterministic lifetime remain undefined; static order has zero locked frequency, not a finite-frequency domain."],
                  elapsed_s=time.perf_counter()-started)
    target = ROOT / "p2_refit/refit.json"
    target.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(json.dumps(dict(movement=movement,
                          new_initial=[dict(eta=x["eta"], ordered=x["initial_domain"]["ordered"],
                                            rho_threshold=x["initial_domain"]["rho_threshold"],
                                            alpha=x["initial_coarse"]["alpha"],
                                            outside_sindy=x["initial_coarse"]["outside_declared_sindy_parameter_ranges"])
                                       for x in stages[-1]["cases"]],
                          ranking=[dict(eta=x["eta"], rows=[dict(parameter=r["parameter"], decades=r["positive_decades"],
                                                               switch=r["zero_nonzero_switch"], positive_count=r["positive_count"],
                                                               zero_count=r["zero_count"]) for r in x["ranking"]]) for x in ranked],
                          elapsed_s=output["elapsed_s"]), indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
