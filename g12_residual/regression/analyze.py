#!/usr/bin/env python3
"""G12 fixed-forward residual consumer. No density fitting or mechanism inputs.

Run from any directory: python g12_residual/regression/analyze.py [--search]
Searches are actual serial weighted PySR; without --search, reuse saved fronts.
Only public figure/table transcriptions and fixed-model producer outputs enter.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
from importlib.metadata import version

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE
FOILS = ("Pd", "PdD", "Ti", "TiD")
COLORS = dict(zip(FOILS, ("#3366aa", "#aa6633", "#339966", "#aa3399")))
SEEDS = (20261124, 20261125)
FEATURES = {
    "E_cm_keV": ("keV", "nominal setting", "search"),
    "host_Z": ("1", "host identity", "search"),
    "electrochemical_loading": ("1", "preparation indicator, not measured D loading", "search"),
    "foil_thickness_um": ("um", "source foil property, constant", "diagnostic"),
    "chart_cup_current_uA": ("uA", "published setting proxy, not run current", "diagnostic"),
    "model_beam_time_s": ("s", "fixed-model assumed time, not measured run", "diagnostic"),
    "model_target_temperature_K": ("K", "fixed-model assumed temperature, not measured run", "diagnostic"),
    "model_csda_Dplus_depth_nm": ("nm", "fixed-model derived depth, not independent measurement", "diagnostic only"),
    "run_current_uA": ("uA", "run measurement, unavailable unless published", "missing"),
    "run_time_s": ("s", "run measurement, unavailable unless published", "missing"),
    "run_dose_C": ("C", "run measurement, unavailable unless published", "missing"),
}
SOURCE_NUMERIC = {
    "temperature_K": ("K", "unpublished run temperature"),
    "actual_D_per_metal": ("1", "unpublished actual loading"),
    "current_density_uA_cm2": ("uA/cm2", "unpublished run current density"),
    "accumulated_dose_C": ("C", "unpublished accumulated dose"),
    "cup_IV_current_uA": ("uA", "published I–V setting proxy, not run current"),
    "stated_run_duration_lo_h": ("h", "stated lower duration bound, not assigned exposure"),
    "stated_run_duration_hi_h": ("h", "stated upper duration bound, not assigned exposure"),
    "spot_area_typical_cm2": ("cm2", "typical apparatus area, not point measurement"),
    "beam_current_typical_uA": ("uA", "typical apparatus current"),
    "D_flux_typical_per_cm2_s": ("1/cm2/s", "typical apparatus flux"),
    "Dplus_fraction_typical": ("1", "typical ion composition"),
    "D2plus_fraction_typical": ("1", "typical ion composition"),
    "retarding_spread_typical_percent": ("percent", "typical energy spread, convention unspecified"),
    "range_SRIM_reported_endpoint_nm": ("nm", "published simulated endpoint, intermediate values absent"),
    "purity_stated_min_pct": ("percent", "host purity bound"),
    "anneal_temperature_C": ("degC", "preparation temperature, not point temperature"),
    "anneal_duration_h": ("h", "preparation duration, not exposure"),
    "vacuum_side_renewal_nm": ("nm", "surface preparation"),
    "electrolyte_side_Pd_nm": ("nm", "opposite-face coating"),
    "operating_pressure_lo_Torr": ("Torr", "apparatus pressure bound"),
    "operating_pressure_hi_Torr": ("Torr", "apparatus pressure bound"),
    "electrolyte_D2SO4_M": ("mol/L", "preparation electrolyte"),
    "galvanostatic_lo_mA_cm2": ("mA/cm2", "preparation current bound, not beam current"),
    "galvanostatic_hi_mA_cm2": ("mA/cm2", "preparation current bound, not beam current"),
    "host_atom_density_cm3": ("1/cm3", "source host context, not actual D density"),
    "D_solubility_terminal_per_host": ("1", "source terminal solubility, not point loading"),
    "D_diffusivity_m2_s": ("m2/s", "unpublished numeric diffusivity"),
    "electron_density_m3": ("1/m3", "unpublished electron density"),
    "Debye_temperature_K": ("K", "unpublished Debye temperature"),
    "conductivity_S_m": ("S/m", "unpublished conductivity"),
}
FEATURES.update({f"source_{key}": (unit, role, "diagnostic; not independent feature")
                 for key, (unit, role) in SOURCE_NUMERIC.items()})


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith("#")))


def save_csv(name, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with (OUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def save_json(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def number(row, key):
    value = row.get(key, "")
    return float(value) if value not in ("", None) else np.nan


def column(rows, key):
    return np.array([number(row, key) for row in rows])


def metrics(actual, predicted, sigma, mask=None):
    mask = np.ones(len(actual), dtype=bool) if mask is None else np.asarray(mask)
    finite = mask & np.isfinite(actual) & np.isfinite(predicted) & np.isfinite(sigma)
    n = int(finite.sum())
    if not n:
        return {"n": 0, "weighted_rmse_log": None, "rmse_log": None, "chi2": None, "max_abs_log": None}
    error = predicted[finite] - actual[finite]
    weights = sigma[finite] ** -2
    return {"n": n, "weighted_rmse_log": float(np.sqrt(weights @ error**2 / weights.sum())),
            "rmse_log": float(np.sqrt(np.mean(error**2))),
            "chi2": float(np.sum((error / sigma[finite])**2)), "max_abs_log": float(np.max(abs(error)))}


def checked_davis(joined_input=None, catalog=None):
    rows = read_csv(HERE / "davis_residuals.csv")
    baseline = json.loads((HERE / "davis_baseline.json").read_text())
    assert baseline["baseline_id"] == "g12_davis_fixed_smooth40_ordinary_thermal_formed_sheet_no_join_k0", "use accepted clean reference, never provisional fitted-stitch front"
    assert hashlib.sha256((HERE / "davis_residuals.csv").read_bytes()).hexdigest() == baseline["residual_csv_sha256"]
    assert len(rows) == 48, "producer must supply the 48 non-anchor plotted rows"
    assert len({r["row_id"] for r in rows}) == 48
    energy = column(rows, "E_cm_keV")
    assert not np.any(energy == 6.5), "reference anchors are deterministic, not residual observations"
    observed, model = column(rows, "R_observed"), column(rows, "R_model")
    target, sigma = column(rows, "log_residual"), column(rows, "sigma_log_residual")
    assert np.all(observed > 0) and np.all(model > 0) and np.all(sigma > 0)
    mismatch = float(np.max(abs(target - (np.log(observed) - np.log(model)))))
    assert mismatch < 1e-10, "target must be ln(Robs/Rmodel), not log of that log"
    assert len({r["lab"] for r in rows}) == 1
    assert all(r["reaction"] for r in rows)
    source_path = joined_input or catalog or (HERE.parent / "data" / "observations.csv")
    public = read_csv(source_path)
    join_input_sha256 = hashlib.sha256(source_path.read_bytes()).hexdigest()
    if joined_input:
        source_by_id = {r["row_id"]: r for r in public}
        assert len(source_by_id) == 48
        for row in rows:
            source = source_by_id[row["row_id"]]
            assert float(source["R_observed"]) == float(row["R_observed"])
            assert float(source["sigma_log_residual"]) == float(row["sigma_log_residual"])
            row["modeled_channel"] = row["channel"]
            row.update({k: v for k, v in source.items() if k.startswith("source_")})
    else:
        source_by_key = {(r["target"], float(r["E_cm_keV"])): r for r in public if r["study_id"] == "karahadian2026"}
        assert len(source_by_key) == 52
        for row in rows:
            source = source_by_key[(row["target"], float(row["E_cm_keV"]))]
            assert float(source["observable"]) == float(row["R_observed"])
            assert float(source["adopted_sigma_ln_R"]) == float(row["sigma_log_residual"])
            row["modeled_channel"] = row["channel"]
            for key in ("channel", "evidence_class", "target_geometry", "target_preparation",
                        "reference_group", "citation", "loading_state", "cup_IV_current_source",
                        "cup_IV_bias_regime", "proton_signal_background_ratio_author", *SOURCE_NUMERIC):
                row[f"source_{key}"] = source[key]
            row["source_reference_current_uA"] = source_by_key[(row["target"], 6.5)]["cup_IV_current_uA"]
    save_csv("analysis_rows.csv", rows)
    save_json("consumer_proof.json", {"rows": 48, "unique_ids": 48, "deterministic_anchors_excluded": True,
        "max_log_identity_error": mismatch, "observed_uncertainty_only": True,
        "baseline_id": baseline["baseline_id"],
        "input_sha256": hashlib.sha256((HERE / "davis_residuals.csv").read_bytes()).hexdigest(),
        "public_join_input_sha256": join_input_sha256,
        "public_join_input_kind": "explicit curated analysis snapshot" if joined_input else "source-owner public intake",
        "joined_public_rows": 48, "public_anchors_excluded": 4,
        "observed_channel_kept_separate_from_modeled_splice": True,
        "target": "ln(R_observed)-ln(R_model), exactly once", "anchor_normalizers": "davis_baseline.json",
        "group_splits": "entire foil held out; one Davis lab cannot supply laboratory holdout"})
    return rows, target, sigma


def characterize(rows, target, sigma):
    inventory, aliases, support = [], [], []
    for key, (unit, role, use) in FEATURES.items():
        x = column(rows, key)
        finite = np.isfinite(x)
        same_at_energy = all(len(set(x[(column(rows, "E_cm_keV") == e) & finite])) <= 1
                             for e in np.unique(column(rows, "E_cm_keV")))
        inventory.append(dict(feature=key, unit=unit, meaning=role, use=use,
                              available=int(finite.sum()), missing=int((~finite).sum()),
                              distinct=len(set(x[finite])), common_by_energy=same_at_energy if finite.any() else ""))
        for foil in FOILS:
            ix = np.array([r["target"] == foil for r in rows]) & finite
            varied = ix.sum() >= 3 and len(set(x[ix])) > 1
            rho = float(spearmanr(x[ix], column(rows, "E_cm_keV")[ix]).statistic) if varied else None
            aliases.append(dict(feature=key, target=foil, rows=int(ix.sum()), spearman_vs_energy=rho,
                                independent_run_measurement=role.startswith("run measurement") and bool(ix.any())))
    save_csv("feature_inventory.csv", inventory)
    save_csv("energy_alias_audit.csv", aliases)
    for foil in FOILS:
        energies = set(column([r for r in rows if r["target"] == foil], "E_cm_keV"))
        for other in FOILS:
            other_e = set(column([r for r in rows if r["target"] == other], "E_cm_keV"))
            support.append(dict(lab=rows[0]["lab"], target=foil, compared_target=other, rows=len(energies),
                                overlap_energy_count=len(energies & other_e), E_min_keV=min(energies),
                                E_max_keV=max(energies), unknown_run_id=all(not r.get("run_id") for r in rows)))
    save_csv("foil_support.csv", support)
    ncols = 4
    nrows = (len(FEATURES) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(18, 3.5 * nrows), constrained_layout=True)
    for ax, (key, (unit, role, _)) in zip(axes.flat, FEATURES.items()):
        for foil in FOILS:
            ix = np.array([r["target"] == foil for r in rows])
            x = column(rows, key)
            finite = ix & np.isfinite(x)
            if finite.any():
                ax.errorbar(x[finite], target[finite], yerr=sigma[finite], fmt="o", ms=3,
                            color=COLORS[foil], label=foil, alpha=.8, capsize=2)
        ax.set_xlabel(f"{key} [{unit}]")
        ax.set_ylabel("ln(Robs/Rmodel)")
        ax.set_title(role, fontsize=8)
        if not np.isfinite(column(rows, key)).any():
            ax.text(.5, .5, "NOT PUBLISHED — no zero substitution", ha="center", transform=ax.transAxes, fontsize=8)
    for ax in list(axes.flat)[len(FEATURES):]:
        ax.axis("off")
    axes.flat[0].legend()
    fig.savefig(OUT / "raw_covariates.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for ax, foil in zip(axes.flat, FOILS):
        ix = np.array([r["target"] == foil for r in rows])
        ax.errorbar(column(rows, "E_cm_keV")[ix], target[ix], yerr=sigma[ix], fmt="o", color=COLORS[foil])
        ax.set(title=foil, xlabel="E_cm [keV]", ylabel="ln(Robs/Rmodel)")
    fig.savefig(OUT / "residual_by_foil.png", dpi=150)
    plt.close(fig)
    # Same-energy contrasts cancel shared setting proxies but not preparation/host/surface confounding.
    contrasts = []
    for a, b in (("PdD", "Pd"), ("TiD", "Ti"), ("Pd", "Ti"), ("PdD", "TiD")):
        left = {r["E_cm_keV"]: i for i, r in enumerate(rows) if r["target"] == a}
        right = {r["E_cm_keV"]: i for i, r in enumerate(rows) if r["target"] == b}
        for e in sorted(left.keys() & right.keys(), key=float):
            i, j = left[e], right[e]
            contrasts.append(dict(target_a=a, target_b=b, E_cm_keV=float(e),
                log_residual_contrast=float(target[i]-target[j]),
                sigma_independent_plotbars_only=float(np.hypot(sigma[i], sigma[j])),
                shared_reference_covariance="unknown; not included"))
    save_csv("same_energy_contrasts.csv", contrasts)


def features(rows, mode):
    e, z, load = (column(rows, k) for k in ("E_cm_keV", "host_Z", "electrochemical_loading"))
    if mode == "raw":
        return np.column_stack((e, z, load)), ["E_keV", "Z", "prepared"]
    if mode == "dimensionless":
        return np.column_stack((e / 6.5, z / 46., load)), ["e", "z", "prepared"]
    # Diagnostic setting-proxy challenge; not an independent measured intervention.
    current = column(rows, "source_cup_IV_current_uA")
    return np.column_stack((e / 6.5, z / 46., load, current / 100.)), ["e", "z", "prepared", "Isetting"]


def baselines(rows, target, sigma, train, test):
    energy = column(rows, "E_cm_keV")
    model = column(rows, "R_model")
    weights = sigma[train] ** -2
    q = np.sqrt(6.5 / energy)
    predictions = {"fixed_forward_reference_only": np.zeros(len(rows)),
                   "training_constant": np.full(len(rows), np.average(target[train], weights=weights))}
    for degree in (1, 2):
        design = np.column_stack([q**j for j in range(degree + 1)])
        beta = np.linalg.lstsq(design[train] * np.sqrt(weights)[:, None], target[train] * np.sqrt(weights), rcond=None)[0]
        predictions[f"energy_only_q_degree{degree}"] = design @ beta
        raw_beta = np.linalg.lstsq(design[train] * np.sqrt(weights)[:, None],
            np.log(column(rows, "R_observed"))[train] * np.sqrt(weights), rcond=None)[0]
        predictions[f"raw_logR_energy_q_degree{degree}"] = design @ raw_beta - np.log(model)
    # Rhat=(Rmodel+c)/(1+c) preserves the deterministic 6.5-keV anchor.
    def floor_prediction(logc):
        return np.logaddexp(np.log(model), logc) - np.logaddexp(0., logc) - np.log(model)
    fit = minimize_scalar(lambda lc: np.sum(weights * (floor_prediction(lc)[train]-target[train])**2),
                          bounds=(-35., 2.), method="bounded")
    predictions["baseline_plus_global_floor"] = floor_prediction(fit.x)
    result = []
    for name, pred in predictions.items():
        result.append(dict(model=name, training=metrics(target[train], pred[train], sigma[train]),
                           held=metrics(target[test], pred[test], sigma[test]),
                           floor_R=float(np.exp(fit.x)) if name == "baseline_plus_global_floor" else None))
    return result, predictions


def symbolic_estimator(run_id, seed, maxsize=12):
    from pysr import PySRRegressor
    return PySRRegressor(niterations=50, timeout_in_seconds=15,
        parallelism="serial", deterministic=True, random_state=seed,
        populations=4, population_size=30, ncycles_per_iteration=150,
        maxsize=maxsize, precision=64, binary_operators=["+", "-", "*", "/"],
        unary_operators=["sqrt", "log"], model_selection="accuracy",
        output_directory=str(OUT / "pysr-work"), run_id=run_id,
        progress=False, verbosity=0, turbo=False, bumper=False)


def search(rows, target, sigma):
    splits = [("all", np.ones(len(rows), bool), np.zeros(len(rows), bool))]
    splits += [(foil, np.array([r["target"] != foil for r in rows]), np.array([r["target"] == foil for r in rows])) for foil in FOILS]
    runs, predictions, baseline_rows, split_rows = [], [], [], []
    for held, train, test in splits:
        base, _ = baselines(rows, target, sigma, train, test)
        baseline_rows += [dict(split=held, **b) for b in base]
        for mode in ("raw", "dimensionless", "alias_challenge"):
            X, names = features(rows, mode)
            if not np.all(np.isfinite(X)):
                runs.append(dict(split=held, mode=mode, status="not evaluable: missing feature; no imputation"))
                continue
            seeds = SEEDS if mode != "alias_challenge" else SEEDS[:1]
            inside = np.all((X >= X[train].min(axis=0)) & (X <= X[train].max(axis=0)), axis=1)
            energy_supported = (X[:, 0] >= X[train, 0].min()) & (X[:, 0] <= X[train, 0].max())
            combination_seen = np.array([any((r["host"], r["electrochemical_loading"]) ==
                (s["host"], s["electrochemical_loading"]) for s, keep in zip(rows, train) if keep) for r in rows])
            for i in np.flatnonzero(test):
                split_rows.append(dict(split=held, mode=mode, row_id=rows[i]["row_id"],
                    energy_supported=bool(energy_supported[i]), feature_box_supported=bool(inside[i]),
                    host_preparation_combination_seen=bool(combination_seen[i]),
                    group_unseen=True, support_meaning="box is marginal, NOT joint physical support"))
            for seed in seeds:
                run_id = f"{held}-{mode}-{seed}"
                print("PySR", run_id, flush=True)
                estimator = symbolic_estimator(run_id, seed)
                weights = sigma[train] ** -2
                weights /= weights.mean()
                estimator.fit(X[train], target[train], weights=weights, variable_names=names)
                front, best_loss = [], np.inf
                for index, eq in estimator.equations_.sort_values("complexity").iterrows():
                    with np.errstate(all="ignore"):
                        pred = np.asarray(estimator.predict(X, index=index), dtype=float)
                    if pred.ndim == 0:
                        pred = np.full(len(rows), float(pred))
                    entry = dict(index=int(index), complexity=int(eq["complexity"]), search_loss=float(eq["loss"]),
                        equation=str(eq["equation"]), sympy=str(eq["sympy_format"]),
                        pareto=bool(eq["loss"] < best_loss),
                        training=metrics(target[train], pred[train], sigma[train]),
                        held_all=metrics(target[test], pred[test], sigma[test]),
                        held_box_supported=metrics(target[test], pred[test], sigma[test], inside[test]),
                        held_joint_host_preparation_supported=metrics(target[test], pred[test], sigma[test], combination_seen[test]),
                        held_nonfinite=int(np.sum(~np.isfinite(pred[test]))))
                    best_loss = min(best_loss, float(eq["loss"]))
                    front.append(entry)
                    for i in np.flatnonzero(test):
                        predictions.append(dict(run=run_id, split=held, mode=mode, seed=seed,
                            complexity=entry["complexity"], row_id=rows[i]["row_id"],
                            log_residual_prediction=float(pred[i]) if np.isfinite(pred[i]) else "",
                            log_R_observed=float(np.log(number(rows[i], "R_observed"))),
                            log_R_reconstructed=float(pred[i]+np.log(number(rows[i], "R_model"))) if np.isfinite(pred[i]) else "",
                            energy_supported=bool(energy_supported[i]), feature_box_supported=bool(inside[i]),
                            combination_seen=bool(combination_seen[i])))
                eligible = [f for f in front if f["complexity"] <= 9]
                selected = min(eligible, key=lambda f: (f["search_loss"], f["complexity"]))
                runs.append(dict(run=run_id, split=held, mode=mode, seed=seed, status="completed",
                    variables=names, train_ids=[r["row_id"] for r, keep in zip(rows, train) if keep],
                    held_ids=[r["row_id"] for r, keep in zip(rows, test) if keep],
                    selection="minimum TRAINING loss among complexity <=9; no held-group selection",
                    selected_complexity=selected["complexity"], selected_equation=selected["equation"], front=front))
                save_json("pysr_fronts.json", {"pysr_version": version("pysr"), "runs": runs})
    save_csv("held_predictions.csv", predictions)
    save_csv("split_support.csv", split_rows)
    save_json("simple_baselines.json", baseline_rows)
    return runs


def denominator_diagnostic(rows, target, sigma):
    energy = column(rows, "E_cm_keV")
    logobs, logmodel = np.log(column(rows, "R_observed")), np.log(column(rows, "R_model"))
    q = np.sqrt(6.5 / energy)
    stats, detail = [], []
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    for foil in FOILS:
        ix = np.array([r["target"] == foil for r in rows])
        for label, values, ax in (("log measured R", logobs, axes[0]), ("minus log model R", -logmodel, axes[1]), ("log residual", target, axes[2])):
            ax.errorbar(q[ix], values[ix], yerr=sigma[ix] if label != "minus log model R" else None,
                        fmt="o", ms=4, label=foil, color=COLORS[foil])
            ax.set(title=label, xlabel="sqrt(6.5 keV / E_cm)", ylabel="natural log, dimensionless")
        for low, name in ((energy <= 2.5, "0.25–2.5 keV"), (np.ones(len(rows), bool), "all observed energies")):
            mask = ix & low
            stats.append(dict(target=foil, band=name, rows=int(mask.sum()),
                log_measured_range=float(np.ptp(logobs[mask])), minus_log_model_range=float(np.ptp(logmodel[mask])),
                log_residual_range=float(np.ptp(target[mask])),
                corr_residual_minus_log_model=float(np.corrcoef(target[mask], -logmodel[mask])[0, 1])))
        for i in np.flatnonzero(ix):
            detail.append(dict(row_id=rows[i]["row_id"], target=foil, E_cm_keV=energy[i], q=q[i],
                ln_R_observed=logobs[i], minus_ln_R_model=-logmodel[i], log_residual=target[i],
                identity_error=target[i]-(logobs[i]-logmodel[i])))
    axes[0].legend()
    fig.savefig(OUT / "denominator_shape.png", dpi=150)
    plt.close(fig)
    save_csv("denominator_shape.csv", detail)
    save_csv("denominator_diagnostic.csv", stats)


def summarize_search():
    path = OUT / "pysr_fronts.json"
    if not path.exists():
        return
    runs = json.loads(path.read_text())["runs"]
    scores, recurrence = [], []
    for run in runs:
        if run["status"] != "completed":
            continue
        for entry in run["front"]:
            scores.append(dict(run=run["run"], split=run["split"], mode=run["mode"], seed=run["seed"],
                complexity=entry["complexity"], equation=entry["equation"], search_loss=entry["search_loss"],
                selected=entry["complexity"] == run["selected_complexity"],
                train_wrmse=entry["training"]["weighted_rmse_log"],
                held_wrmse=entry["held_all"]["weighted_rmse_log"], held_n=entry["held_all"]["n"],
                held_box_n=entry["held_box_supported"]["n"], held_joint_n=entry["held_joint_host_preparation_supported"]["n"],
                held_nonfinite=entry["held_nonfinite"]))
        eq = run["selected_equation"]
        recurrence.append(dict(run=run["run"], split=run["split"], mode=run["mode"], seed=run["seed"],
            equation=eq, complexity=run["selected_complexity"], uses_sqrt="sqrt" in eq, uses_log="log" in eq,
            uses_host=bool(re.search(r"\b(Z|z)\b", eq)), uses_preparation="prepared" in eq,
            uses_current_proxy="Isetting" in eq, uses_model_depth="depth_model" in eq,
            inference="syntax audit only; recurrence is not mechanism identification"))
    save_csv("complexity_scores.csv", scores)
    save_csv("expression_recurrence.csv", recurrence)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for run in runs:
        if run["status"] != "completed" or run["mode"] != "dimensionless":
            continue
        f = run["front"]
        color = COLORS.get(run["split"], "black")
        axes[0].plot([a["complexity"] for a in f], [a["training"]["weighted_rmse_log"] for a in f], "o-", alpha=.65, color=color, label=run["split"] if run["seed"] == SEEDS[0] else None)
        if run["split"] != "all":
            axes[1].plot([a["complexity"] for a in f], [a["held_all"]["weighted_rmse_log"] for a in f], "o-", alpha=.65, color=color)
    axes[0].set(title="Training weighted error: all retained complexities", xlabel="complexity", ylabel="weighted RMSE [ln ratio]")
    axes[1].set(title="Held foil — training-only model selection", xlabel="complexity", ylabel="weighted RMSE [ln ratio]")
    axes[0].legend()
    fig.savefig(OUT / "pareto_fronts.png", dpi=150)
    plt.close(fig)


def readout_sensitivity(rows, target, sigma):
    primary = target.copy()
    k1 = np.log(column(rows, "R_observed")) - np.log(column(rows, "R_model_k1"))
    factor = column(rows, "R_model_k1") / column(rows, "R_model")
    error = float(np.max(abs(k1 - primary + np.log(factor))))
    assert error < 1e-12
    detail, summary = [], []
    for i, row in enumerate(rows):
        detail.append(dict(row_id=row["row_id"], target=row["target"], E_cm_keV=row["E_cm_keV"],
            log_residual_k0=primary[i], log_residual_k1_fixed_observation=k1[i],
            delta_k1_minus_k0=k1[i]-primary[i], model_current_factor=factor[i],
            sigma_observation_only=sigma[i]))
    for key in (k for k in rows[0] if k.startswith("delta_log_residual_") and k.endswith("_minus_primary")):
        values = column(rows, key)
        if np.isfinite(values).all():
            summary.append(dict(case=key, max_abs_delta_log=float(np.max(abs(values))),
                meaning="model sensitivity, not observation uncertainty; no common covariance invented"))
    summary.append(dict(case="k1_fixed_observations_minus_k0", max_abs_delta_log=float(np.max(abs(k1-primary))),
        analytic_identity_error=error, meaning="readout ambiguity sensitivity; not a second measurement"))
    save_csv("readout_sensitivity.csv", detail)
    save_csv("model_sensitivity_summary.csv", summary)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for foil in FOILS:
        ix = np.array([r["target"] == foil for r in rows])
        energy = column(rows, "E_cm_keV")
        for ax, values, label in ((axes[0], primary, "k0 as plotted"), (axes[1], k1, "k1 model-only extra current power")):
            ax.errorbar(energy[ix], values[ix], yerr=sigma[ix], fmt="o", ms=3, color=COLORS[foil], label=foil)
            ax.set(title=label, xlabel="E_cm [keV]", ylabel="ln(Robs/Rmodel)")
    axes[0].legend()
    fig.savefig(OUT / "current_readings.png", dpi=150)
    plt.close(fig)
    return k1


def reconstruct_held_observation(rows, target, sigma):
    detail, score_rows = [], []
    logmodel = np.log(column(rows, "R_model"))
    logobs = np.log(column(rows, "R_observed"))
    pred_path = OUT / "held_predictions.csv"
    all_predictions = read_csv(pred_path) if pred_path.exists() else []
    fronts = json.loads((OUT / "pysr_fronts.json").read_text())["runs"] if pred_path.exists() else []
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, foil in zip(axes.flat, FOILS):
        test = np.array([r["target"] == foil for r in rows])
        train = ~test
        scores, predictions = baselines(rows, target, sigma, train, test)
        score_rows += [dict(split=foil, **s) for s in scores]
        energy = column(rows, "E_cm_keV")
        ix = np.flatnonzero(test)
        ax.errorbar(energy[test], logobs[test], yerr=sigma[test], fmt="ko", ms=4, label="observed")
        for name, pred in predictions.items():
            for i in ix:
                detail.append(dict(split=foil, model=name, row_id=rows[i]["row_id"],
                    E_cm_keV=energy[i], log_R_observed=logobs[i], log_R_reconstructed=logmodel[i]+pred[i],
                    log_residual_prediction=pred[i], training_only=True))
            if name in ("baseline_plus_global_floor", "energy_only_q_degree2", "raw_logR_energy_q_degree2"):
                ax.plot(energy[test], (pred+logmodel)[test], "o-", ms=2, label=name)
        if all_predictions:
            selected = next(r for r in fronts if r.get("split") == foil and r.get("mode") == "dimensionless" and r.get("seed") == SEEDS[0])
            by_id = {r["row_id"]: number(r, "log_R_reconstructed") for r in all_predictions
                     if r["run"] == selected["run"] and int(r["complexity"]) == selected["selected_complexity"]}
            ax.plot(energy[test], [by_id[rows[i]["row_id"]] for i in ix], "o--", ms=2,
                    label="PySR dim; fixed seed; train-selected")
        ax.set(title=f"Unseen {foil}: raw measured log R", xlabel="E_cm [keV]", ylabel="ln R")
    axes.flat[0].legend(fontsize=6)
    fig.savefig(OUT / "held_observed_reconstruction.png", dpi=150)
    plt.close(fig)
    save_csv("baseline_held_predictions.csv", detail)
    save_json("simple_baselines.json", score_rows)


def external_analysis(path):
    """Keep native response classes separate; no inferred Ue is made into energy rows."""
    from collections import defaultdict
    import sys
    from scipy.integrate import quad
    sys.path.insert(0, str(HERE.parents[1]))
    from g4_beam.forward import model as nuclear
    source = read_csv(path)
    assert len({r["row_id"] for r in source}) == len(source)
    blocks, output = defaultdict(list), []
    def source_integral(e, screening=.040):
        return quad(lambda x: float(nuclear.sigma(x, screening, "p"))/np.sqrt(x),
                    .05, e, epsabs=1e-18, epsrel=1e-8)[0]
    def thick_shape_ti(e):
        # Source homogeneous D/Ti=2; changed approximation: existing ICRU electronic
        # coefficient, not the unavailable exact Anderson–Ziegler source curve.
        host = nuclear.ELEM["Ti"]
        def integrand(x):
            lab = 2*x
            stop = (nuclear.Se_atom(lab, "Ti") + nuclear.Sn_atom(lab, host["Z"], host["M"])
                    + 2*nuclear.Sn_atom(lab, 1, nuclear.M_D))
            return float(nuclear.sigma(x, .040, "p") / stop)
        return quad(integrand, .05, e, epsabs=1e-18, epsrel=1e-8)[0]
    for row in source:
        kind, reaction = row["observable_kind"], row["reaction"]
        e, observed = number(row, "E_cm_keV"), number(row, "observable")
        native, model, reason = "", np.nan, "source comparison only; no compatible pointwise denominator"
        if reaction != "D+D":
            reason = "different reaction; never mixed into DD targets"
        elif row["evidence_class"] == "reported_fit":
            reason = "author fitted/inferred parameter or limit, not a yield; unknown E does not define unique residual"
        elif kind == "S_factor" and row["observable_unit"] == "keV*b" and np.isfinite(e):
            factor40 = float(np.sqrt(e/(e+.040)) * np.exp(nuclear.BG*(1/np.sqrt(e)-1/np.sqrt(e+.040))))
            ta = row["study_id"] == "raiola2002_epja13_ta"
            model = ((53.0+.48*e) if ta else (43.0+.54*e))*factor40
            native = "DD_extracted_S_over_source_linear_Sbare_times_f40"
            exact = row["study_id"] in ("raiola2002_epja13_ta", "raiola2002_plb547", "raiola2006_bochum_thesis")
            reason = ("source bare polynomial " if exact else "conditional source-family bare polynomial; exact article normalization not independently matched; ") + ("Ta53+.48E" if ta else "43+.54E") + " keV*b times Eq7 f40; standard S=E*sigma*exp(2pi eta), BG31.397; extracted S not per-charge yield; BH40 sensitivity separate"
        elif row["study_id"] == "huke-berlin" and kind == "normalized_thick_target_enhancement":
            ref = number(row, "reference_E_cm_keV")
            model = (source_integral(e)/source_integral(e, 0.))/(source_integral(ref)/source_integral(ref, 0.))
            native = "DD_reduced_enhancement_over_conditional_BH_gas_f40"
            reason = "source reduced F=y(E)/y(12.5CM); conditional BH gas integral weighted1/sqrt(E), normalized at same preparation; exact source gas evaluation not reproduced, no thin-yield substitution"
        elif row["study_id"] in ("yuki_kasagi1998_jetp68", "kasagi2002_jpsj71", "kasagi1998_iccf7") and kind == "relative_thick_target_yield":
            ref = number(row, "reference_E_cm_keV")
            model = source_integral(e)/source_integral(ref)
            native = "DD_source_relative_proton_yield_over_homogeneous_sqrt_stopping40"
            reason = "BH p40eV homogeneous velocity-proportional stopping comparator, .05CM cutoff, own reported5CM p reference; source-concentration/layered-target profile unresolved; conditional, not Davis profile"
        elif row["study_id"] == "karahadian2026":
            reason = "handled by clean exact fixed-profile Davis producer; four deterministic reference rows not independent fits"
        elif row["study_id"] == "yuki1997" and kind == "relative_thick_target_yield" and row["host"] == "Ti" and np.isfinite(e):
            ref = number(row, "reference_E_cm_keV")
            assert ref > .05 and e > .05
            model = thick_shape_ti(e) / thick_shape_ti(ref)
            native = "DD_Ti_relative_yield40_changed_stopping_approximation"
            reason = "conditional homogeneous D/Ti=2; BH p40eV; host-only ICRU electronic plus host/D ZBL nuclear stopping, not exact source Anderson–Ziegler; own3.225CM reference; cutoff .05CM; no Davis focusing/join"
        elif row["study_id"] == "yuki1997" and kind == "relative_thick_target_yield":
            reason = "Yb source relative yields retained; exact source stopping coefficient unavailable, no unsupported element extrapolation or invented profile"
        elif kind == "detected_thick_target_yield":
            reason = "detected counts/(uA*ks)=1000 counts/C, not emitted reactions/C; acceptance and unfitted source normalizer unavailable; distinct source beam components"
        elif kind.startswith("neutron_rate"):
            reason = "native rate/intervention contrast at15keV; no point current/loading/beam spectrum normalization, not DD yield per charge"
        residual = float(np.log(observed / model)) if np.isfinite(model) and model > 0 and observed > 0 else ""
        sigma_minus, sigma_plus = number(row, "sigma_minus"), number(row, "sigma_plus")
        output.append(dict(row_id=row["row_id"], study_id=row["study_id"], lab=row["lab"], reaction=reaction,
            measurement_block=row.get("measurement_block", ""), evidence_class=row["evidence_class"],
            channel=row["channel"], target=row["target"], host=row["host"], E_cm_keV=row["E_cm_keV"],
            target_preparation=row["target_preparation"], conditions_note=row["conditions_note"],
            source_sigma_minus=row["sigma_minus"], source_sigma_plus=row["sigma_plus"],
            per_run_beam_current_uA=row.get("per_run_beam_current_uA", ""),
            actual_D_per_metal=row.get("actual_D_per_metal", ""),
            reported_inferred_D_per_metal=row.get("reported_inferred_D_per_metal", ""),
            temperature_K=row.get("temperature_K", ""), observable_kind=kind, observable=row["observable"],
            observable_unit=row["observable_unit"], reference_group=row["reference_group"],
            reference_E_cm_keV=row["reference_E_cm_keV"],
            native_residual_kind=native, model_same_observable=float(model) if np.isfinite(model) else "",
            log_residual_BH40_S_sensitivity=(float(np.log(observed/(nuclear.sigma(e, .040, "p")*e*np.exp(nuclear.BG/np.sqrt(e))/1000.)))
                if native.startswith("DD_extracted_S") else ""),
            log_residual=residual,
            sigma_log_first_order_minus=float(sigma_minus/observed) if native and observed > 0 and np.isfinite(sigma_minus) else "",
            sigma_log_first_order_plus=float(sigma_plus/observed) if native and observed > 0 and np.isfinite(sigma_plus) else "",
            status=reason, citation=row["citation"], source_locator=row["source_locator"]))
        blocks[(row["study_id"], row["lab"], reaction, kind, row["channel"], row["target_geometry"], row["evidence_class"])].append(output[-1])
    grid = []
    for key, block in blocks.items():
        energies = np.array([number(r, "E_cm_keV") for r in block])
        finite = np.isfinite(energies)
        grid.append(dict(study_id=key[0], lab=key[1], reaction=key[2], observable_kind=key[3],
            channel=key[4], target_geometry=key[5], rows=len(block), energies_known=int(finite.sum()),
            evidence_class=key[6], hosts=";".join(sorted({r["host"] for r in block})),
            unique_targets=len({r["target"] for r in block}),
            actual_run_current_rows=sum(bool(r["per_run_beam_current_uA"]) for r in block),
            actual_loading_rows=sum(bool(r["actual_D_per_metal"]) for r in block),
            inferred_loading_rows=sum(bool(r["reported_inferred_D_per_metal"]) for r in block),
            temperature_rows=sum(bool(r["temperature_K"]) for r in block),
            sigma_rows=sum(bool(r["source_sigma_plus"]) for r in block),
            E_min_keV=float(energies[finite].min()) if finite.any() else "",
            E_max_keV=float(energies[finite].max()) if finite.any() else "",
            native_residual_rows=sum(r["log_residual"] != "" for r in block),
            unique_reference_groups=len({r["reference_group"] for r in block if r["reference_group"]}),
            pooled_with_Davis=False, limitation=block[0]["status"]))
    save_csv("external_residuals.csv", output)
    save_csv("external_support_grid.csv", grid)
    save_json("external_input.json", {"public_repository_path": "g12_residual/data/observations.csv", "data_commit": "e89133b",
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "rows": len(source),
        "all_source_rows_preserved_in_support_accounting": True,
        "study_keys": len({r["study_id"] for r in source}), "source_rows_not_independent_samples": True,
        "native_S_convention": "standard S=E*sigma*exp(2pi eta); source-family thesis Sbare43+.54E, not exact inaccessible IOP coefficient; BG=31.397 from existing kernel",
        "pooled_Davis_fixed_profile_LOLO": "not evaluable: no second laboratory has the same published profile/readout target",
        "source_only_is_not_zero_residual": True})
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    for host in ("Pt", "Hf"):
        subset = [r for r in output if r["study_id"] == "raiola2005" and r["target"] == host and r["native_residual_kind"].startswith("DD_extracted_S")]
        for temp in sorted({r["temperature_K"] for r in subset}):
            rr = [r for r in subset if r["temperature_K"] == temp]
            axes[0, 0].errorbar(column(rr, "E_cm_keV"), column(rr, "log_residual"),
                yerr=column(rr, "sigma_log_first_order_plus"), fmt="o", ms=3, label=f"{host} {temp} K")
    for host in ("Zr", "Ti", "Pd"):
        rr = [r for r in output if r["study_id"] == "czerski2409v3" and r["target"] == host
              and r["observable_kind"] == "detected_thick_target_yield"]
        if rr:
            axes[0, 1].plot(column(rr, "E_cm_keV"), np.log(column(rr, "observable")), "o", ms=3, label=host)
    for host in ("Ti", "Yb"):
        rr = [r for r in output if r["study_id"] == "yuki1997" and r["target"] == host
              and r["observable_kind"] == "relative_thick_target_yield"]
        if rr:
            axes[1, 0].plot(column(rr, "E_cm_keV"), np.log(column(rr, "observable")), "o", ms=3, label=host)
    for host in sorted({r["target"] for r in source if r["study_id"] == "raiola2005"}):
        rr = [r for r in source if r["study_id"] == "raiola2005" and r["target"] == host
              and r["observable_kind"] == "fitted_screening_energy" and r["temperature_K"]]
        if rr:
            axes[1, 1].errorbar(column(rr, "temperature_K"), column(rr, "observable"),
                               yerr=column(rr, "sigma_plus"), fmt="o-", ms=3, label=host)
    for ax in axes.flat:
        ax.legend(fontsize=6)
    axes[0, 0].set(title="O1379 conditional native S / thesis-bare×40eV", xlabel="E_cm [keV]", ylabel="ln S ratio")
    axes[0, 1].set(title="Czerski: detected yields; bars unavailable", xlabel="E_cm [keV]", ylabel="ln(counts/(uA ks))")
    axes[1, 0].set(title="Yuki: same-host per-dose relative yields", xlabel="E_cm [keV]", ylabel="ln R (3.225-keV own reference)")
    axes[1, 1].set(title="Reported Ue estimates — not regression targets", xlabel="temperature [K]", ylabel="reported Ue [eV]")
    fig.savefig(OUT / "independent_observables.png", dpi=150)
    plt.close(fig)
    interventions = []
    for row in source:
        if row["observable_kind"] == "neutron_rate_electrochemical_percent_increase":
            fraction = number(row, "observable") / 100.
            interventions.append(dict(row_id=row["row_id"], E_cm_keV=row["E_cm_keV"],
                log_rate_ratio=float(np.log1p(fraction)),
                sigma_log_first_order=number(row, "sigma_plus")/(100.*(1+fraction)),
                kind="measured intervention rate contrast, not per-charge fixed-profile residual",
                citation=row["citation"]))
    save_csv("independent_interventions.csv", interventions)
    source_native_added(source, output, nuclear)


def source_native_added(source, native, nuclear):
    """Source-response decompositions, never mixing S, rates and DD yields."""
    from collections import defaultdict
    from scipy.integrate import quad
    result, grouped = [], defaultdict(list)
    for row in native:
        if row["log_residual"] == "":
            continue
        observed, model = float(row["observable"]), float(row["model_same_observable"])
        group = f'{row["study_id"]}:{row["target"]}:{row["temperature_K"]}:{row["reference_group"]}:{row["native_residual_kind"]}'
        unit = row["observable_unit"]
        if row["native_residual_kind"].startswith("DD_reduced_enhancement"):
            e, ref = float(row["E_cm_keV"]), float(row["reference_E_cm_keV"])
            def bare_shape(x):
                return quad(lambda q: float(nuclear.sigma(q, 0., "p"))/np.sqrt(q),
                            .05, x, epsabs=1e-18, epsrel=1e-8)[0]
            gas_ratio = bare_shape(e)/bare_shape(ref)
            observed, model = observed*gas_ratio, model*gas_ratio
            unit = "1, reconstructed same-preparation relative yield under conditional BH gas integral"
        result.append(dict(row_id=row["row_id"], study_id=row["study_id"], series=group,
            E_cm_keV=row["E_cm_keV"], observed_native=observed, model_native=model,
            added_native=observed-model, native_unit=unit,
            log_ratio=float(np.log(observed/model)), primary_Davis_compatible=False,
            qualification=row["status"], citation=row["citation"]))
    for row in source:
        if row["study_id"] == "czerski2409v3" and row["observable_kind"] == "detected_thick_target_yield":
            panel = row["source_locator"].split(" PDF")[0]
            beam = row["conditions_note"].split(".")[0]
            grouped[(row["target"], panel, beam)].append(row)
    def cold(e):
        return quad(lambda x: float(nuclear.sigma(x, .040, "p"))/np.sqrt(x),
                    .05, e, epsabs=1e-18, epsrel=1e-8)[0]
    for (target, panel, beam), rr in grouped.items():
        anchor = max(rr, key=lambda r: number(r, "E_cm_keV"))
        ref, anchor_y = number(anchor, "E_cm_keV"), number(anchor, "observable")
        for row in rr:
            e = number(row, "E_cm_keV")
            observed = number(row, "observable")/anchor_y
            model = cold(e)/cold(ref)
            result.append(dict(row_id=row["row_id"], study_id=row["study_id"],
                series=f'czerski:{target}:{panel}:{beam}', E_cm_keV=e,
                observed_native=observed, model_native=model, added_native=observed-model,
                native_unit="1, self-rebased detected yield",
                log_ratio=float(np.log(observed/model)), primary_Davis_compatible=False,
                anchor_row_id=anchor["row_id"], anchor_E_cm_keV=ref,
                qualification="secondary shape diagnostic: highest-energy measured point is shared reference; covariance unknown; homogeneous sqrt-stopping BH40 cold kernel, source neutral components unmodeled, no fitted M/Ue or acceptance borrowed",
                citation=row["citation"]))
    save_csv("external_added_terms.csv", result)
    groups = defaultdict(list)
    for row in result:
        groups[row["series"]].append(row)
    summaries = []
    for group, rr in groups.items():
        rr = sorted(rr, key=lambda r: float(r["E_cm_keV"]))
        low = [r for r in rr if float(r["E_cm_keV"]) <= 2.5]
        low = low if len(low) >= 3 else rr[:3]
        obs = column(low, "observed_native")
        model = column(low, "model_native")
        added = column(low, "added_native")
        summaries.append(dict(series=group, study_id=rr[0]["study_id"], low_rows=len(low),
            low_E_min_keV=float(low[0]["E_cm_keV"]), low_E_max_keV=float(low[-1]["E_cm_keV"]),
            native_unit=rr[0]["native_unit"], observed_log_span=float(np.ptp(np.log(obs))),
            model_log_span=float(np.ptp(np.log(model))), added_min=float(added.min()), added_max=float(added.max()),
            all_low_added_positive=bool(np.all(added > 0)), primary_Davis_compatible=False,
            qualification=rr[0]["qualification"]))
    save_csv("external_added_summary.csv", summaries)


def profile_policy_view(rows, target, sigma):
    energy = column(rows, "E_cm_keV")
    logobs = np.log(column(rows, "R_observed"))
    logmodel = np.log(column(rows, "R_model"))
    fig, axes = plt.subplots(3, 4, figsize=(16, 10), sharex=True, constrained_layout=True)
    records = []
    for j, foil in enumerate(FOILS):
        ix = np.array([r["target"] == foil for r in rows])
        for ax, values, label in zip(axes[:, j], (logobs, logmodel, target),
                                     ("ln observed R", "ln fixed model R", "ln observed/model")):
            ax.errorbar(energy[ix], values[ix], yerr=sigma[ix] if label != "ln fixed model R" else None,
                        fmt="o-", ms=3, color=COLORS[foil])
            ax.axvspan(2., 2.5, color="grey", alpha=.18)
            ax.set_ylabel(label)
        axes[0, j].set_title(foil + ": fixed exposure policy changes")
        axes[2, j].set_xlabel("E_cm [keV]; shaded 2.0→2.5")
        weights = sigma[ix] ** -2
        centered = target[ix] - np.average(target[ix], weights=weights)
        for i, value in zip(np.flatnonzero(ix), centered):
            records.append(dict(row_id=rows[i]["row_id"], target=foil, E_cm_keV=energy[i],
                model_assumed_time_h=number(rows[i], "model_beam_time_s")/3600.,
                log_R_observed=logobs[i], log_R_model=logmodel[i], log_residual=target[i],
                weighted_within_foil_centered_residual=value,
                reference_nuisance="unknown shared log offset; centering is diagnostic, not a correction/variance estimate",
                boundary="energy and exposure change simultaneously; no isolated time effect"))
    fig.savefig(OUT / "profile_policy_step.png", dpi=150)
    plt.close(fig)
    save_csv("profile_policy_step.csv", records)


def laboratory_shape_stress(rows, target, sigma):
    """True held-lab training; secondary shape transfer, NOT pooled-profile validation."""
    native = read_csv(HERE / "external_residuals.csv")
    yuki = [r for r in native if r["study_id"] == "yuki1997" and r["target"] == "Ti"
            and r["native_residual_kind"] == "DD_Ti_relative_yield40_changed_stopping_approximation"]
    davis_ix = np.array([r["target"] == "Ti" and r["modeled_channel"] == "p" for r in rows])
    if len(yuki) < 3:
        save_json("laboratory_shape_stress.json", {"status": "not evaluable", "reason": "missing conditional Yuki Ti native block"})
        return
    groups = {
        "UC Davis": (column(rows, "E_cm_keV")[davis_ix], target[davis_ix], sigma[davis_ix],
                     [r["row_id"] for r, keep in zip(rows, davis_ix) if keep]),
        "Tohoku": (column(yuki, "E_cm_keV"), column(yuki, "log_residual"),
                   (column(yuki, "sigma_log_first_order_minus")+column(yuki, "sigma_log_first_order_plus"))/2,
                   [r["row_id"] for r in yuki]),
    }
    runs = []
    for training_lab, held_lab in (("UC Davis", "Tohoku"), ("Tohoku", "UC Davis")):
        e, y, s, ids = groups[training_lab]
        he, hy, hs, held_ids = groups[held_lab]
        support = (he >= e.min()) & (he <= e.max())
        order = np.argsort(he[support])
        actual_delta = np.diff(hy[support][order])
        def shape_metric(pred):
            delta_error = np.diff(pred[support][order]) - actual_delta
            return {"supported_rows": int(support.sum()), "adjacent_differences": len(delta_error),
                    "adjacent_difference_rmse_log": float(np.sqrt(np.mean(delta_error**2))) if len(delta_error) and np.isfinite(delta_error).all() else None,
                    "nonfinite_predictions": int(np.sum(~np.isfinite(pred)))}
        estimator = symbolic_estimator("shape-transfer-" + training_lab.replace(" ", "_"), SEEDS[0])
        errors_available = bool(np.all(np.isfinite(s) & (s > 0)))
        weights = s**-2 if errors_available else np.ones(len(y))
        weights /= weights.mean()
        estimator.fit((e/6.5)[:, None], y, weights=weights, variable_names=["e"])
        weight_status = "published plot-uncertainty weights" if errors_available else "unit optimization weights only; most source errors unpublished, no uncertainty imputed"
        front = []
        for index, eq in estimator.equations_.sort_values("complexity").iterrows():
            pred = np.asarray(estimator.predict((he/6.5)[:, None], index=index))
            if pred.ndim == 0:
                pred = np.full(len(he), float(pred))
            score = shape_metric(pred)
            front.append(dict(complexity=int(eq["complexity"]), equation=str(eq["equation"]),
                search_loss=float(eq["loss"]), held_shape=score))
        selected = min((f for f in front if f["complexity"] <= 9), key=lambda f: (f["search_loss"], f["complexity"]))
        simple = []
        for degree in (0, 1, 2):
            q, hq = np.sqrt(6.5/e), np.sqrt(6.5/he)
            design = np.column_stack([q**j for j in range(degree+1)])
            held_design = np.column_stack([hq**j for j in range(degree+1)])
            beta = np.linalg.lstsq(design*np.sqrt(weights)[:, None], y*np.sqrt(weights), rcond=None)[0]
            simple.append(dict(model=f"energy_q_degree{degree}", held_shape=shape_metric(held_design@beta)))
        # Missing source errors stay blank; unit optimization weights are not measurement errors.
        runs.append(dict(training_lab=training_lab, held_lab=held_lab, train_ids=ids, held_ids=held_ids,
            training_weight_status=weight_status,
            feature="e=E_cm/(6.5keV); no yield/denominator input", front=front,
            selected_complexity=selected["complexity"], selected_equation=selected["equation"],
            selection="training loss only at complexity<=9", simple_baselines=simple,
            supported_held_ids=[held_ids[i] for i in np.flatnonzero(support)],
            excluded_held_ids=[held_ids[i] for i in np.flatnonzero(~support)]))
    save_json("laboratory_shape_stress.json", {"status": "conditional domain-transfer stress only",
        "differences": "adjacent-energy log-residual differences within overlap, cancel own-reference level; no held-lab intercept fitted",
        "limitations": "Ti only; proton numerator; Davis n6.5 vs Yuki p3.225 reference; fixed transported vs homogeneous x2/changed stopping; adjacent differences dependent, no significance",
        "runs": runs})


def package_capsule():
    import gzip
    import io
    import tarfile
    generated = (
        "analysis_rows.csv", "consumer_proof.json", "feature_inventory.csv", "energy_alias_audit.csv",
        "foil_support.csv", "raw_covariates.png", "residual_by_foil.png", "same_energy_contrasts.csv",
        "denominator_shape.png", "denominator_shape.csv", "denominator_diagnostic.csv",
        "pysr_fronts.json", "held_predictions.csv", "split_support.csv", "simple_baselines.json",
        "complexity_scores.csv", "expression_recurrence.csv", "pareto_fronts.png",
        "readout_sensitivity.csv", "model_sensitivity_summary.csv", "current_readings.png",
        "baseline_held_predictions.csv", "held_observed_reconstruction.png",
        "profile_policy_step.png", "profile_policy_step.csv", "laboratory_shape_stress.json",
        "additive_rows.csv", "additive_summary.csv", "additive_loaded_ratios.csv",
        "SI1_normalizations.csv", "SI1_normalization_summary.csv", "measured_model_overlay.png",
        "additive_scatter.png", "SI1_normalization_shapes.png", "additive_fronts.json",
        "additive_held_predictions.csv", "additive_simple_baselines.json", "additive_fold_scores.csv",
        "additive_recurrence.csv", "additive_fronts.png")
    names = ["analyze.py", "davis_forward.py", "davis_residuals.csv", "davis_baseline.json",
             "README.md", "inputs_ledger.md", "external_residuals.csv", "external_support_grid.csv",
             "external_input.json", "independent_observables.png", "independent_interventions.csv",
             "external_added_terms.csv", "external_added_summary.csv", "deposition.csv", "deposition_provenance.json"]
    names += list(generated) + ["k1/" + name for name in generated]
    payload = {name: (HERE / name).read_bytes() for name in sorted(names)}
    manifest = {"format": "g12-regression-whitelist-v1", "base_commit": "bf114e0874f2a01dc62505e7c7520cadf46fd20f",
        "excludes": ["provisional_preledger", "pysr-work", "logs", "PDFs", "underlying raw datasets", "source-owner data partition"],
        "files": [{"path": "g12_residual/regression/" + name, "bytes": len(data),
                   "sha256": hashlib.sha256(data).hexdigest()} for name, data in payload.items()]}
    manifest_bytes = (json.dumps(manifest, indent=2) + "\n").encode()
    (HERE / "capsule_manifest.json").write_bytes(manifest_bytes)
    payload["capsule_manifest.json"] = manifest_bytes
    target = HERE / "g12-regression.tar.gz"
    with target.open("wb") as stream, gzip.GzipFile(filename="", mode="wb", fileobj=stream, mtime=0) as zipped:
        with tarfile.open(fileobj=zipped, mode="w") as archive:
            for name, data in sorted(payload.items()):
                info = tarfile.TarInfo("g12_residual/regression/" + name)
                info.size, info.mtime, info.mode = len(data), 0, 0o644
                archive.addfile(info, io.BytesIO(data))
    with tarfile.open(target, "r:gz") as archive:
        members = archive.getmembers()
        assert len(members) == len(payload) and all(m.isfile() and m.name.startswith("g12_residual/regression/") and ".." not in Path(m.name).parts for m in members)
        for member in members:
            name = member.name.removeprefix("g12_residual/regression/")
            assert archive.extractfile(member).read() == payload[name]
    print(json.dumps({"path": "g12_residual/regression/" + target.name, "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
                      "bytes": target.stat().st_size, "media_type": "application/gzip", "regular_members_verified": len(payload)}))


def prepare_deposition():
    """Extend only the two missing low-energy G8 deposition settings, no transport replay."""
    import sys
    sys.path[:0] = [str(HERE.parents[1] / "g8_why"), str(HERE.parents[1])]
    import common as g8
    from g4_beam.forward.consumer import _csda_segments
    from g4_beam.forward.profile_thermal import transport_stopping
    saved_path = HERE.parents[1] / "g8_why" / "driven_target_budget.csv"
    saved = {(r["material"], float(r["E_cm_keV"])): r for r in read_csv(saved_path)}
    output = []
    for foil in FOILS:
        for e in (.25, .5, 1., 1.5, 2., 2.5, 6.5):
            key = foil, e
            final_energies = []
            if key in saved:
                deposited = float(saved[key]["deposited_eV_per_incident_charge"])
                origin = "G8 saved matched frozen-profile CSDA deposition"
            else:
                layers, _, _ = g8.reference(foil).profiles(e)
                frozen = tuple(tuple(float(v) for v in layer) for layer in layers)
                deposited = 0.
                for n, fraction in ((1, .85), (2, .15)):
                    for g, weight in zip(g8.GX, g8.GW):
                        lab = 2*e/n*(1+.05/2.3548*g)
                        parts = _csda_segments(float(lab), frozen, transport_stopping, .1)
                        final = float(parts[-1][3](parts[-1][1])[0]) if parts else lab
                        final_energies.append(final)
                        deposited += n*fraction*weight*(lab-final)*1000.
                origin = "G12 extension of identical G8 frozen-profile CSDA accumulator"
                assert max(abs(np.array(final_energies)-.1)) < 1e-6
            cutoff_accounting = 2000*e - 115.
            assert abs(deposited-cutoff_accounting) < .001
            output.append(dict(target=foil, E_cm_keV=e, deposited_eV_per_incident_charge=deposited,
                deposited_J_per_true_C=deposited,
                feature_deposited_J_per_true_C=round(deposited, 3),
                cutoff_only_accounting_J_per_C=cutoff_accounting,
                energy_balance_error_J_per_C=deposited-cutoff_accounting,
                max_terminal_lab_keV=max(final_energies) if final_energies else "",
                min_terminal_lab_keV=min(final_energies) if final_energies else "",
                source=origin, meaning="model deposited stopping energy to0.1keV/D cutoff, not energy in a selected mode"))
    save_csv("deposition.csv", output)
    save_json("deposition_provenance.json", {"source": "g8_why/driven_target_budget.csv",
        "source_sha256": hashlib.sha256(saved_path.read_bytes()).hexdigest(),
        "profile_sha256": hashlib.sha256((HERE.parents[1]/"g4_beam/merit/profiles.json").read_bytes()).hexdigest(),
        "saved_rows_reused": 20, "new_matched_settings": 8, "no_profile_evolution": True,
        "charge_conversion": "one eV per elementary charge equals one J/C numerically; 85/15 electrical mix gives1.15D/e",
        "cutoff_accounting": "2000*E_cm_keV−1.15*100 J/C; all incident paths stop at0.1keV/D; no coupling-to-mode assertion",
        "feature_resolution": "rounded to0.001J/C, removing sub-micro-eV ODE noise, not a measurement uncertainty",
        "max_abs_saved_or_extended_balance_error_J_per_C": max(abs(r["energy_balance_error_J_per_C"]) for r in output)})


def additive_analysis(rows, reading, do_search=False):
    """Direct signed added term; unit transforms invert before every held score."""
    energy = column(rows, "E_cm_keV")
    obs, model = column(rows, "R_observed"), column(rows, "R_model")
    added = obs-model
    sigma = (column(rows, "sigma_minus")+column(rows, "sigma_plus"))/2
    low = energy <= 2.5
    assert low.sum() == 24 and np.all(added[low] > 0)
    dep_by_key = {(r["target"], float(r["E_cm_keV"])): r for r in read_csv(HERE / "deposition.csv")}
    dep = np.array([number(dep_by_key[(r["target"], float(r["E_cm_keV"]))], "feature_deposited_J_per_true_C")
                    if float(r["E_cm_keV"]) <= 2.5 else np.nan for r in rows])
    dep_ref = np.array([number(dep_by_key[(r["target"], 6.5)], "feature_deposited_J_per_true_C") for r in rows])
    anchors = json.loads((HERE / "davis_baseline.json").read_text())["reference_anchors"]
    cref = {r["target"]: r["cup_focusing_fraction"] for r in anchors}
    focus_relative = column(rows, "cup_focusing_fraction") / np.array([cref[r["target"]] for r in rows])
    current_relative = column(rows, "source_cup_IV_current_uA") / column(rows, "source_reference_current_uA")
    primary_model = column(rows, "R_model_k0") if reading == "k1" else model
    model_extra_factor = column(rows, "R_model_k1") / primary_model
    cup_factor = 1/model_extra_factor if reading == "k1" else np.ones(len(rows))
    factors = {"plotted": np.ones(len(rows)), "cup_C": cup_factor,
        "true_C_or_D": cup_factor*focus_relative, "time": cup_factor*current_relative,
        "deposited_J": cup_factor*focus_relative*dep_ref/dep}
    output, summaries, ratios, absolute = [], [], [], []
    for i, row in enumerate(rows):
        item = dict(row_id=row["row_id"], target=row["target"], E_cm_keV=energy[i], reading=reading,
            R_observed=obs[i], R_model=model[i], delta_R=added[i], sigma_delta_R_observation_only=sigma[i],
            low_energy=bool(low[i]), measured_own_reference_absolute_yield="unknown",
            source_I_setting_uA=number(row, "source_cup_IV_current_uA"),
            model_deposited_J_per_true_C=dep[i] if low[i] else "",
            normalizer_status="dimensionless own6.5-relative proxies; chart I not run current; c/deposition modeled")
        for name, factor in factors.items():
            item["delta_"+name] = added[i]*factor[i] if np.isfinite(factor[i]) else ""
        output.append(item)
        sb = number(row, "source_proton_signal_background_ratio_author")
        current_A = number(row, "source_cup_IV_current_uA")*1e-6
        if low[i] and np.isfinite(sb):
            rate = sb*.007/3600.
            qcup = rate/current_A
            qtrue = qcup*number(row, "cup_focusing_fraction")
            lo_h, hi_h = (number(row, "source_stated_run_duration_"+side+"_h") for side in ("lo", "hi"))
            absolute.append(dict(row_id=row["row_id"], target=row["target"], E_cm_keV=energy[i],
                detected_per_s_proxy=rate, detected_per_cup_C_proxy=qcup,
                detected_per_true_C_model_collection_proxy=qtrue,
                detected_per_incident_D_model_mix_proxy=qtrue*1.602176634e-19/1.15,
                detected_per_model_deposited_J_proxy=qtrue/dep[i],
                model_deposited_J_per_true_C=dep[i],
                conditional_counts_over_stated_low_duration=rate*lo_h*3600 if np.isfinite(lo_h) else "",
                conditional_counts_over_stated_high_duration=rate*hi_h*3600 if np.isfinite(hi_h) else "",
                uncertainty="not published; no fabricated error", actual_run_current_A="",
                status="SI1 detected signal proxy; approximate added counts only where physical baseline negligible; no exact acceptance/6.5 cross-calibration"))
    save_csv("additive_rows.csv", output)
    save_csv("SI1_normalizations.csv", absolute)
    for foil in FOILS:
        mask = low & np.array([r["target"] == foil for r in rows])
        for name, factor in factors.items():
            value, err = added[mask]*factor[mask], sigma[mask]*factor[mask]
            weights = err**-2
            design = np.column_stack((np.ones(mask.sum()), energy[mask]))
            beta = np.linalg.lstsq(design*np.sqrt(weights)[:, None], value*np.sqrt(weights), rcond=None)[0]
            keep = energy[mask] > .25
            dd, ww = design[keep], weights[keep]
            beta_without = np.linalg.lstsq(dd*np.sqrt(ww)[:, None], value[keep]*np.sqrt(ww), rcond=None)[0]
            formal_slope_sigma = float(np.sqrt(np.linalg.inv(design.T@(weights[:, None]*design))[1, 1]))
            summaries.append(dict(target=foil, reading=reading, normalization=name, rows=int(mask.sum()),
                weighted_level=float(np.average(value, weights=weights)), minimum=float(value.min()), maximum=float(value.max()),
                log_span=float(np.ptp(np.log(value))), population_SD_log=float(np.std(np.log(value))),
                weighted_slope_per_keV=float(beta[1]), weighted_slope_without_025_per_keV=float(beta_without[1]),
                endpoint_net_slope_per_keV=float((value[-1]-value[0])/(energy[mask][-1]-energy[mask][0])),
                formal_diagonal_sigma_slope=formal_slope_sigma,
                formal_diagonal_one_sigma_includes_zero=bool(abs(beta[1]) <= formal_slope_sigma),
                formal_diagonal_one_sigma_low=float(beta[1]-formal_slope_sigma),
                formal_diagonal_one_sigma_high=float(beta[1]+formal_slope_sigma),
                formal_diagonal_sigma_level=float(1/np.sqrt(weights.sum())),
                covariance="shared reference/point correlations unknown; formal slope scale not significance",
                units="relative proxy in stated normalization, not absolute common rate"))
    for a, b in (("PdD", "Pd"), ("TiD", "Ti")):
        left = {float(r["E_cm_keV"]): r for r in output if r["target"] == a and r["low_energy"]}
        right = {float(r["E_cm_keV"]): r for r in output if r["target"] == b and r["low_energy"]}
        for e in sorted(left.keys() & right.keys()):
            la, lb = left[e], right[e]
            ratio = la["delta_R"]/lb["delta_R"]
            firstorder = abs(ratio)*np.hypot(la["sigma_delta_R_observation_only"]/la["delta_R"], lb["sigma_delta_R_observation_only"]/lb["delta_R"])
            ratios.append(dict(numerator=a, denominator=b, E_cm_keV=e, reading=reading,
                ratio_of_separately_normalized_added_curves=ratio, first_order_plotbars_uncorrelated_only=firstorder,
                max_normalizer_ratio_difference=max(abs(la["delta_"+n]/lb["delta_"+n]-ratio) for n in factors),
                meaning="not absolute loaded/unloaded yield enhancement; own reference constants unknown"))
    save_csv("additive_summary.csv", summaries)
    save_csv("additive_loaded_ratios.csv", ratios)
    si_summary = []
    for foil in FOILS:
        rr = [r for r in absolute if r["target"] == foil]
        for key in ("detected_per_s_proxy", "detected_per_cup_C_proxy",
                    "detected_per_true_C_model_collection_proxy", "detected_per_incident_D_model_mix_proxy",
                    "detected_per_model_deposited_J_proxy"):
            v, ee = column(rr, key), column(rr, "E_cm_keV")
            si_summary.append(dict(target=foil, normalization=key, rows=len(rr), minimum=float(v.min()), maximum=float(v.max()),
                population_SD_log=float(np.std(np.log(v))), endpoint_ratio=float(v[-1]/v[0]),
                unweighted_signed_slope_per_keV=float(np.polyfit(ee, v, 1)[0]),
                slope_without_025=float(np.polyfit(ee[1:], v[1:], 1)[0]),
                uncertainty="unknown; no formal confidence; actual time/charge not reconstructed"))
    save_csv("SI1_normalization_summary.csv", si_summary)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, foil in zip(axes.flat, FOILS):
        ix = np.array([r["target"] == foil for r in rows])
        ax.errorbar(energy[ix], np.log(obs[ix]), yerr=column(rows, "sigma_log_residual")[ix], fmt="o-", label="observed ln R")
        ax.plot(energy[ix], np.log(model[ix]), "o-", label="clean model ln R")
        ax.axvspan(2., 2.5, color="grey", alpha=.2)
        ax.set(title=foil, xlabel="E_cm [keV]", ylabel="natural log, identical axis")
    axes.flat[0].legend()
    fig.savefig(OUT / "measured_model_overlay.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), constrained_layout=True)
    for j, foil in enumerate(FOILS):
        ix = np.array([r["target"] == foil for r in rows])
        for ax, mask, label in ((axes[0, j], ix, "signed full range"), (axes[1, j], ix & low, "low-energy added term")):
            ax.errorbar(energy[mask], added[mask], yerr=sigma[mask], fmt="o-", ms=3, color=COLORS[foil])
            ax.axhline(0., color="grey", lw=.7)
            ax.axvspan(2., 2.5, color="grey", alpha=.15)
            ax.set(title=foil+": "+label, xlabel="E_cm [keV]", ylabel="Robs−Rmodel, own6.5 units")
    fig.savefig(OUT / "additive_scatter.png", dpi=150)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for ax, foil in zip(axes.flat, FOILS):
        rr = [r for r in absolute if r["target"] == foil]
        for key, label in (("detected_per_s_proxy", "per second"), ("detected_per_cup_C_proxy", "per cup C"),
                           ("detected_per_model_deposited_J_proxy", "per modeled deposited J")):
            values = column(rr, key)
            ax.plot(column(rr, "E_cm_keV"), values/np.exp(np.mean(np.log(values))), "o-", label=label)
        ax.set(title=foil+": SI1 scales / geometric mean", xlabel="E_cm [keV]", ylabel="dimensionless shape")
    axes.flat[0].legend()
    fig.savefig(OUT / "SI1_normalization_shapes.png", dpi=150)
    plt.close(fig)
    if do_search:
        additive_search(rows, added, sigma, low, factors, current_relative, dep/dep_ref, reading)
    additive_readout()


def additive_search(rows, added, sigma, low, factors, current, deposition, reading):
    """Short direct-added-term search scored in original ΔR after inversion."""
    energy = column(rows, "E_cm_keV")
    X = np.column_stack((energy/6.5, current, column(rows, "host_Z")/46.,
                         column(rows, "electrochemical_loading"), deposition))
    names = ["e", "Isetting", "z", "prepared", "dep_model"]
    splits = [("all", low.copy(), np.zeros(len(rows), bool))]
    splits += [(foil, low & np.array([r["target"] != foil for r in rows]),
                low & np.array([r["target"] == foil for r in rows])) for foil in FOILS]
    runs, predicted, simple_rows = [], [], []
    for norm, factor in factors.items():
        transformed, transformed_sigma = added*factor/.01, sigma*factor/.01
        for held, train, test in splits:
            cup = factors["cup_C"]
            cup_values, cup_sigma = added[train]/cup[train], sigma[train]/cup[train]
            cup_mean = np.average(cup_values, weights=cup_sigma**-2)
            candidate_mean = np.average(transformed[train], weights=transformed_sigma[train]**-2)
            design = np.column_stack((np.ones(len(rows)), energy))
            weights_raw = sigma[train]**-2
            beta = np.linalg.lstsq(design[train]*np.sqrt(weights_raw)[:, None], added[train]*np.sqrt(weights_raw), rcond=None)[0]
            nulls = {"constant_per_cup_C": cup*cup_mean, "signed_linear_deltaR": design@beta,
                     "constant_in_candidate_units": np.full(len(rows), candidate_mean)*.01/factor}
            for key, pred in nulls.items():
                simple_rows.append(dict(normalization=norm, split=held, model=key,
                    training_original_deltaR=metrics(added[train], pred[train], sigma[train]),
                    held_original_deltaR=metrics(added[test], pred[test], sigma[test])))
            for seed in (SEEDS if held == "all" else SEEDS[:1]):
                estimator = symbolic_estimator(f"add-{norm}-{held}-{seed}", seed, maxsize=7)
                weights = transformed_sigma[train]**-2
                weights /= weights.mean()
                estimator.fit(X[train], transformed[train], weights=weights, variable_names=names)
                front = []
                for index, eq in estimator.equations_.sort_values("complexity").iterrows():
                    with np.errstate(all="ignore"):
                        pred_scaled = np.asarray(estimator.predict(X[low], index=index), dtype=float)
                    if pred_scaled.ndim == 0:
                        pred_scaled = np.full(low.sum(), float(pred_scaled))
                    pred = np.full(len(rows), np.nan)
                    pred[low] = pred_scaled*.01/factor[low]
                    front.append(dict(complexity=int(eq["complexity"]), equation=str(eq["equation"]),
                        search_loss_candidate_scaled=float(eq["loss"]),
                        training_original_deltaR=metrics(added[train], pred[train], sigma[train]),
                        held_original_deltaR=metrics(added[test], pred[test], sigma[test]),
                        held_nonfinite=int(np.sum(~np.isfinite(pred[test])))))
                    for i in np.flatnonzero(test):
                        predicted.append(dict(reading=reading, normalization=norm, held_foil=held, seed=seed,
                            complexity=int(eq["complexity"]), row_id=rows[i]["row_id"],
                            delta_R=added[i], predicted_delta_R=float(pred[i]) if np.isfinite(pred[i]) else "",
                            reconstructed_R=float(number(rows[i], "R_model")+pred[i]) if np.isfinite(pred[i]) else "",
                            sigma_observation=sigma[i], energy_supported=True, host_preparation_combination_seen=False))
                selected = min((f for f in front if f["complexity"] <= 5),
                               key=lambda f: (f["search_loss_candidate_scaled"], f["complexity"]))
                runs.append(dict(reading=reading, normalization=norm, split=held, seed=seed,
                    train_ids=[r["row_id"] for r, keep in zip(rows, train) if keep],
                    held_ids=[r["row_id"] for r, keep in zip(rows, test) if keep],
                    selected_complexity=selected["complexity"], selected_equation=selected["equation"], front=front,
                    selection="minimum TRAINING candidate-unit loss at complexity<=5, no held selection",
                    features="e=E/6.5; Isetting=source chart I/Iref; z=Z/46; prepared0/1; dep_model=deposition/ref; energy-linked, not independent controls",
                    scale="direct transformed ΔR divided by fixed0.01; held scores invert to original ΔR"))
                save_json("additive_fronts.json", {"runs": runs})
    save_csv("additive_held_predictions.csv", predicted)
    save_json("additive_simple_baselines.json", simple_rows)


def additive_readout():
    path = OUT / "additive_fronts.json"
    if not path.exists():
        return
    data = json.loads(path.read_text())
    scores, recurrence = [], []
    for run in data["runs"]:
        for front in run["front"]:
            for key in ("training_original_deltaR", "held_original_deltaR"):
                front[key] = {name.replace("_log", "_deltaR"): value for name, value in front[key].items()}
            scores.append(dict(reading=run["reading"], normalization=run["normalization"], split=run["split"],
                seed=run["seed"], complexity=front["complexity"], equation=front["equation"],
                selected=front["complexity"] == run["selected_complexity"],
                training_candidate_scaled_loss=front["search_loss_candidate_scaled"],
                training_wrmse_deltaR=front["training_original_deltaR"]["weighted_rmse_deltaR"],
                held_wrmse_deltaR=front["held_original_deltaR"]["weighted_rmse_deltaR"],
                held_max_abs_deltaR=front["held_original_deltaR"]["max_abs_deltaR"],
                held_rows=front["held_original_deltaR"]["n"], held_nonfinite=front["held_nonfinite"],
                units="original own6.5-relative ΔR; each candidate transform inverted"))
        eq = run["selected_equation"]
        recurrence.append(dict(reading=run["reading"], normalization=run["normalization"], split=run["split"],
            seed=run["seed"], equation=eq, complexity=run["selected_complexity"],
            uses_energy=bool(re.search(r"\be\b", eq)), uses_host=bool(re.search(r"\bz\b", eq)),
            uses_preparation="prepared" in eq, uses_setting_current="Isetting" in eq,
            uses_model_deposition="dep_model" in eq,
            meaning="syntax recurrence; host has only two levels, current/deposition alias energy; not a mechanism"))
    data["score_unit"] = "original own6.5-relative ΔR, not logarithmic"
    data["pysr_version"] = version("pysr")
    save_json("additive_fronts.json", data)
    simple_path = OUT / "additive_simple_baselines.json"
    simple = json.loads(simple_path.read_text())
    for row in simple:
        for key in ("training_original_deltaR", "held_original_deltaR"):
            row[key] = {name.replace("_log", "_deltaR"): value for name, value in row[key].items()}
    save_json("additive_simple_baselines.json", simple)
    save_csv("additive_fold_scores.csv", scores)
    save_csv("additive_recurrence.csv", recurrence)
    norms = ("plotted", "cup_C", "true_C_or_D", "time", "deposited_J")
    fig, axes = plt.subplots(2, 5, figsize=(18, 7), constrained_layout=True)
    for j, norm in enumerate(norms):
        for run in data["runs"]:
            if run["normalization"] != norm:
                continue
            f = run["front"]
            if run["split"] == "all":
                axes[0, j].plot([a["complexity"] for a in f],
                    [a["training_original_deltaR"]["weighted_rmse_deltaR"] for a in f], "o-", ms=3,
                    label=str(run["seed"]))
            else:
                axes[1, j].plot([a["complexity"] for a in f],
                    [a["held_original_deltaR"]["weighted_rmse_deltaR"] for a in f], "o-", ms=3,
                    color=COLORS[run["split"]], label=run["split"])
        axes[0, j].set_title(norm+": training")
        axes[1, j].set_title(norm+": held foil")
        for ax in axes[:, j]:
            ax.set(xlabel="complexity", ylabel="weighted RMSE, original ΔR")
    axes[0, 0].legend(fontsize=6)
    axes[1, 0].legend(fontsize=6)
    fig.savefig(OUT / "additive_fronts.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search", action="store_true")
    parser.add_argument("--reading", choices=("k0", "k1"), default="k0")
    parser.add_argument("--external", type=Path)
    parser.add_argument("--external-only", action="store_true")
    parser.add_argument("--transfer", action="store_true")
    parser.add_argument("--joined-input", type=Path, help="explicit curated per-row source snapshot for portable Davis replay")
    parser.add_argument("--catalog", type=Path, help="joined source catalog; staging path may differ from public repository path")
    parser.add_argument("--package", action="store_true")
    parser.add_argument("--prepare-deposition", action="store_true")
    parser.add_argument("--additive-search", action="store_true")
    parser.add_argument("--additive-only", action="store_true")
    args = parser.parse_args()
    if args.package:
        package_capsule()
        return
    assert all(os.environ.get(k) == "1" for k in ("JULIA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")), "set numerical threads explicitly to one"
    global OUT
    if args.reading == "k1":
        OUT = HERE / "k1"
        OUT.mkdir(exist_ok=True)
    if args.prepare_deposition:
        prepare_deposition()
        return
    if args.external:
        external_analysis(args.external.resolve())
    if args.external_only:
        return
    rows, target, sigma = checked_davis(args.joined_input, args.catalog or args.external)
    k1_target = readout_sensitivity(rows, target, sigma)
    if args.reading == "k1":
        for row, value in zip(rows, k1_target):
            row["R_model_k0"] = row["R_model"]
            row["log_residual_k0"] = row["log_residual"]
            row["R_model"] = row["R_model_k1"]
            row["log_residual"] = value
        target = k1_target
        save_csv("analysis_rows.csv", rows)
    additive_analysis(rows, args.reading, do_search=args.additive_search)
    if args.additive_only:
        return
    characterize(rows, target, sigma)
    denominator_diagnostic(rows, target, sigma)
    profile_policy_view(rows, target, sigma)
    if args.search:
        search(rows, target, sigma)
    if args.transfer:
        laboratory_shape_stress(rows, target, sigma)
    summarize_search()
    reconstruct_held_observation(rows, target, sigma)
    print("Davis consumer finished", flush=True)


if __name__ == "__main__":
    main()
