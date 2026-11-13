#!/usr/bin/env python3
"""Actual bounded PySR fits to sweep.csv; no hand-written regression stand-in.

Requires PySR (with working Julia runtime), NumPy and the delivered sweep.csv.
Run with JULIA_NUM_THREADS=OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1.
Search budgets are 20s each, serial, seed20261113, not a convergence claim.
The supplied static symmetry reduction is explicit, NOT an inferred law:
G=f*rho^2/(1+eta*rho^2), w=tanh(1/(2theta)); C and dU/(Ec*r) depend only
on G,theta. The thermal feature w is provided, not claimed discovered by PySR.
Threshold fit uses w,f,eta, NOT the known target denominator or rho_c itself.
C and energy fit the ordered branch; normal-branch values are exactly zero.
Heat: the thermal-resistance candidate is
Rhat=(tau_th/Cv)/(1000s/[3e6 J/(m^3 K)]), dimensionless.
Steady heat depends on tau_th,Cv via that ratio ONLY with all other inputs
fixed, not independently of bridging, nuclear Q, refill or domain groups.
The FULL complexity/error front is written, not one selected winner.
No holdout, cross-validation, precision ladder or model-existence inference.
"""
import argparse
import csv
import json
import time
import numpy as np
import pysr
from pysr import PySRRegressor

SEED = 20261113


def read_table():
    with open("sweep.csv") as stream:
        rows = list(csv.DictReader(stream))
    values = {key: np.array([float(row[key]) for row in rows])
              for key in ("rho", "r", "theta", "eta", "f", "C", "dU_Ec", "ordered",
                          "k_site_s", "Q_MeV", "heat_fraction", "Cv_J_m3_K",
                          "tau_th_s", "tau_refill_s", "T0_K", "n_ref_m3", "heat_long_W_m3")}
    values["rho_c"] = np.array([float(row["rho_c"]) if row["rho_c"] else np.nan for row in rows])
    return values


def fit_target(name, X, y, names, original_scale):
    estimator = PySRRegressor(
        niterations=40, timeout_in_seconds=20, parallelism="serial", deterministic=True,
        random_state=SEED, populations=4, population_size=30, ncycles_per_iteration=150,
        maxsize=18, precision=64, binary_operators=["+", "-", "*", "/"],
        unary_operators=["sqrt", "exp", "log", "tanh"],
        batching=True, batch_size=min(256, len(y)), model_selection="accuracy",
        output_directory="pysr-work", run_id=name, progress=False, verbosity=0,
        turbo=False, bumper=False)
    start = time.perf_counter()
    estimator.fit(X, y, variable_names=names)
    elapsed = time.perf_counter()-start
    # Hall of fame: retain every nondominated complexity/loss entry.
    front = []
    previous_loss = np.inf
    for index, row in estimator.equations_.sort_values("complexity").iterrows():
        if float(row["loss"]) >= previous_loss:
            continue
        previous_loss = float(row["loss"])
        prediction = estimator.predict(X, index=index)
        front.append({"complexity": int(row["complexity"]), "search_loss": float(row["loss"]),
                      "equation": str(row["equation"]), "sympy": str(row["sympy_format"]),
                      "training_mse_fit_target": float(np.mean((prediction-y)**2)),
                      "training_mse_original_target": float(np.mean(((prediction-y)*original_scale)**2)),
                      "training_max_abs_original_target": float(np.max(abs((prediction-y)*original_scale)))})
    return {"name": name, "rows": len(y), "variables": names, "elapsed_s": elapsed, "front": front}


def heat_target(data):
    ordered = data["ordered"] == 1
    Rref = 1000/3e6  # m^3 K/W
    resistance = data["tau_th_s"]/data["Cv_J_m3_K"]
    power_scale = data["n_ref_m3"]*data["Q_MeV"]*1e6*1.602176634e-19*data["heat_fraction"]*data["k_site_s"]
    load = resistance*power_scale/data["T0_K"]
    X = np.c_[data["f"]*data["rho"]**2, data["eta"]*data["rho"]**2, data["theta"],
              data["k_site_s"]*data["tau_refill_s"], load, resistance/Rref]
    result = fit_target("heat_fraction_of_available_power", X[ordered],
                        data["heat_long_W_m3"][ordered]/power_scale[ordered],
                        ["coupling_num", "diamag_num", "theta", "refill_load", "thermal_load", "Rhat"],
                        power_scale[ordered])
    result["reconstruction"] = "P_heat=n_ref*Qe*h*k_site*F(coupling_num=f*rho^2,diamag_num=eta*rho^2,theta,refill_load=k_site*tau_refill,thermal_load,Rhat)"
    result["thermal_ratio_scope"] = "Rhat=(tau_th/Cv)/(1000/3e6); thermal_load=(tau_th/Cv)*n_ref*Qe*h*k_site/T0; other inputs not held irrelevant"
    return result


def inverse_threshold_target(data):
    finite = np.isfinite(data["rho_c"])
    w = np.tanh(1/(2*data["theta"]))
    return fit_target("inverse_rho_c_squared",
                      np.c_[w[finite], data["f"][finite], data["eta"][finite]],
                      1/data["rho_c"][finite]**2, ["w", "f", "eta"], np.ones(np.sum(finite)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--heat-only", action="store_true",
                       help="Append newly requested heat target without replaying completed searches")
    group.add_argument("--inverse-threshold-only", action="store_true",
                       help="One authorized well-posed target diagnostic; no other fit replay")
    args = parser.parse_args()
    data = read_table()
    if args.heat_only or args.inverse_threshold_only:
        with open("pysr-results.json") as stream:
            output = json.load(stream)
    else:
        w = np.tanh(1/(2*data["theta"]))
        G = data["f"]*data["rho"]**2/(1+data["eta"]*data["rho"]**2)
        finite = np.isfinite(data["rho_c"])
        ordered = data["ordered"] == 1
        threshold = fit_target("rho_c", np.c_[w[finite], data["f"][finite], data["eta"][finite]],
                               data["rho_c"][finite], ["w", "f", "eta"], np.ones(np.sum(finite)))
        X = np.c_[G[ordered], data["theta"][ordered], w[ordered]]
        C = fit_target("C", X, data["C"][ordered], ["G", "theta", "w"], np.ones(np.sum(ordered)))
        energy = fit_target("dU_over_gap", X, data["dU_Ec"][ordered]/data["r"][ordered],
                            ["G", "theta", "w"], data["r"][ordered])
        output = {"pysr_version": pysr.__version__, "seed": SEED,
                  "static_coordinates": {"G": "f*rho^2/(1+eta*rho^2)", "w": "tanh(1/(2theta))"},
                  "known_threshold_for_comparison_NOT_seeded": "rho_c=1/sqrt(f*w-eta), only f*w>eta",
                  "targets": [threshold, C, energy],
                  "energy_reconstruction": "dU/Ec = r*F_energy(G,theta,w) on ordered branch, else0",
                  "normal_branch": "C=dU/Ec=0 when (f*w-eta)*rho^2<=1",
                  "limitations": "training errors only; bounded search, physically supplied features; no guarantee of known-law recovery or physical effect"}
    if not args.inverse_threshold_only:
        output["targets"].append(heat_target(data))
    if not args.heat_only:
        output["targets"].append(inverse_threshold_target(data))
    with open("pysr-results.json", "w") as stream:
        json.dump(output, stream, indent=2, allow_nan=False)
    with open("pysr-front.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["target", "complexity", "search_loss", "equation", "sympy", "training_mse_fit_target", "training_mse_original_target", "training_max_abs_original_target"])
        writer.writeheader()
        for target in output["targets"]:
            writer.writerows({"target": target["name"], **row} for row in target["front"])
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
