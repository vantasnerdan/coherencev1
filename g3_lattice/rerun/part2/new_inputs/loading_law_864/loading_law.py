"""Anchored directional-policy tangents; not a relaxed-well or full tensor law.

Each observable is y_x1 + c_epsilon*epsilon + c_v*n_vac_NN, optionally +c_pair*pairs.
The spacing coordinate follows the approved composition/reference-spacing path.
It is NOT an independent local elastic strain. All y and coefficients use eV/A2.
"""
import numpy as np

SCHEMA = "loading-anchored-affine-law-v2"
A_REFERENCE = 4.084


def features(n_vac_NN, a_reference_A, *, vacant_NN_pair_count=None):
    n, a = float(n_vac_NN), float(a_reference_A)
    if not np.isfinite([n, a]).all() or not 0 <= n <= 12 or a <= 0:
        raise ValueError("Expected finite 0<=n_vac_NN<=12 and positive reference spacing")
    value = [a / A_REFERENCE - 1, n]
    if vacant_NN_pair_count is not None:
        pairs = float(vacant_NN_pair_count)
        if not np.isfinite(pairs) or not 0 <= pairs <= min(24, max(0, n * (n - 1) / 2)):
            raise ValueError("Invalid adjacent-pair count within the vacant nearest-neighbour shell")
        value.append(pairs)
    return np.array(value)


def _row_features(row, with_motif):
    return features(row["n_vac_NN"], row["a_reference_A"],
                    vacant_NN_pair_count=row["vacant_NN_pair_count"] if with_motif else None)


def _coefficients(x, delta):
    scale = np.max(np.abs(x), axis=0) if len(x) else np.zeros(x.shape[1])
    if np.any(scale == 0):
        return None
    coefficient, _, rank, singular = np.linalg.lstsq(x / scale, delta, rcond=None)
    if rank != x.shape[1]:
        return None
    return coefficient / scale, scale, singular


def _fit_observable(rows, anchor, name, descriptors):
    selected = [r for r in rows if r["observables"][name] is not None]
    x = np.array([_row_features(r, len(descriptors) == 3) for r in selected]).reshape(-1, len(descriptors))
    y0 = anchor["observables"][name]
    y = np.array([r["observables"][name] for r in selected])
    solution = _coefficients(x, y - y0)
    if solution is None:
        return dict(identified=False, reason=f"{len(descriptors)} slope coefficients are not identifiable",
                    anchor_eV_A2=y0, included_rows=[r["id"] for r in selected])
    coefficient, scale, singular = solution
    residual = y0 + x @ coefficient - y
    from scipy.spatial import ConvexHull
    domain = np.vstack([np.zeros(len(descriptors)), x]) / scale
    hull = ConvexHull(domain)
    folds = []
    for group in sorted({r["config_id"] for r in selected}):
        keep = np.array([r["config_id"] != group for r in selected])
        train = _coefficients(x[keep], y[keep] - y0)
        fold = dict(held_out_config=group, identified=train is not None)
        if train is not None:
            fold["predictions"] = [dict(row=selected[i]["id"], actual_eV_A2=float(y[i]),
                predicted_eV_A2=float(y0 + x[i] @ train[0]),
                error_eV_A2=float(y0 + x[i] @ train[0] - y[i])) for i in np.flatnonzero(~keep)]
        folds.append(fold)
    return dict(identified=True, anchor_eV_A2=y0,
                coefficients_eV_A2=dict(zip(descriptors, map(float, coefficient))),
                scaled_design_condition_number=float(singular[0] / singular[-1]),
                residual_convention="prediction minus actual; anchor is fixed, not a fitted observation",
                residuals=[dict(row=r["id"], predicted_eV_A2=float(y0 + xx @ coefficient),
                    actual_eV_A2=float(yy), residual_eV_A2=float(rr))
                    for r, xx, yy, rr in zip(selected, x, y, residual)],
                rms_residual_eV_A2=float(np.sqrt(np.mean(residual**2))),
                max_abs_residual_eV_A2=float(np.abs(residual).max()),
                residual_degrees_of_freedom=len(selected) - len(descriptors),
                leave_one_configuration_out=folds,
                domain=dict(feature_order=descriptors,
                    scale=scale.tolist(), normalized_halfspaces=hull.equations.tolist(),
                    normalized_vertices=domain[hull.vertices].tolist(),
                    interpretation="Convex hull for THIS observable; no unsampled configuration guarantee"))


def fit(rows, anchor, *, with_motif=False):
    descriptors = ["epsilon", "n_vac_NN"] + (["vacant_NN_pair_count"] if with_motif else [])
    if not rows or _row_features(anchor, with_motif).any():
        raise ValueError("Need loading rows and the fixed zero-descriptor x1 anchor")
    names = sorted(anchor["observables"])
    if any(set(r["observables"]) != set(names) for r in rows):
        raise ValueError("Observable/shell inventories differ")
    return dict(schema=SCHEMA, reference_a_A=A_REFERENCE,
                formula="y=y_x1+c_epsilon*(a_reference_A/4.084-1)+c_v*n_vac_NN"
                    + ("+c_pair*vacant_NN_pair_count" if with_motif else ""),
                descriptor_names=descriptors,
                model_id="single_pair_motif_alternative" if with_motif else "baseline_two_descriptor",
                descriptor_selection="Baseline preregistered776; single motif alternative authorized after the motif contrast by axiom852",
                anchor=anchor, training_rows=rows,
                models={name: _fit_observable(rows, anchor, name, descriptors) for name in names},
                uncertainty="Equal-weight descriptive fit, not a statistical confidence interval. "
                    "SCF/mesh/finite-cell/reference-spacing and descriptor errors remain. "
                    "Configuration holdouts reuse the same x1 anchor.",
                scope="Kii_dd and mean J_dd per existing periodic bond under the registered probe-axis policy. "
                    "No loading quartic, relaxed minimum, full Hessian, independent hydrostatic derivative "
                    "or global heterogeneous Gamma is inferred. Apply lattice occupancy once.")


def evaluate(n_vac_NN, a_reference_A, law, *, vacant_NN_pair_count=None, extrapolate=False):
    if law["schema"] != SCHEMA or law["reference_a_A"] != A_REFERENCE:
        raise ValueError("Unsupported loading law/reference")
    with_motif = len(law["descriptor_names"]) == 3
    if with_motif and vacant_NN_pair_count is None:
        raise ValueError("Motif law requires the actual vacant_NN_pair_count")
    x = features(n_vac_NN, a_reference_A,
                 vacant_NN_pair_count=vacant_NN_pair_count if with_motif else None)
    values, outside = {}, []
    for name, model in law["models"].items():
        if not model["identified"]:
            raise ValueError(f"Unidentified observable: {name}")
        domain = model["domain"]
        equations = np.array(domain["normalized_halfspaces"])
        is_outside = bool(np.any(equations[:, :-1] @ (x / domain["scale"]) + equations[:, -1] > 1e-10))
        if is_outside:
            outside.append(name)
        coefficient = np.array([model["coefficients_eV_A2"][key] for key in law["descriptor_names"]])
        values[name] = float(model["anchor_eV_A2"] + x @ coefficient)
    if outside and not extrapolate:
        raise ValueError("Outside sampled convex hull: " + ", ".join(outside))
    return dict(values_eV_A2=values, extrapolated_observables=outside,
                direction_policy="Clean +100; vacancy-adjacent ray toward original vacant H5 (110 family). "
                    "Isotropic/scalar application at other directions is an unmeasured model closure.")
