"""Fit a low-order cubic material-frame closure to the complete saved tensor table."""
import hashlib
import json
from pathlib import Path

import numpy as np
from strain_tensor_law import NAMES, feature_tensor, stiffness_tensor

UPPER = np.triu_indices(3)


def solve(rows, K0, size):
    X = np.concatenate([feature_tensor(row["strain"])[:size, UPPER[0], UPPER[1]].T for row in rows])
    y = np.concatenate([(np.asarray(row["K_material_eV_A2"])-K0*np.eye(3))[UPPER] for row in rows])
    scale = np.linalg.norm(X, axis=0)
    if np.any(scale == 0):
        raise ValueError("Unidentified strain coefficient")
    values, residual, rank, singular = np.linalg.lstsq(X/scale, y, rcond=1e-10)
    if rank != size:
        raise ValueError(f"Unidentified strain basis: rank{rank}/{size}")
    return values/scale, dict(rank=int(rank), scaled_condition=float(singular[0]/singular[-1]),
                              scaled_singular_values=singular.tolist())


def errors(observed, predicted):
    delta = predicted-observed
    return dict(max_abs_entry_eV_A2=float(np.max(np.abs(delta))),
                operator_norm_eV_A2=float(np.linalg.norm(delta,ord=2)),
                max_eigenvalue_error_eV_A2=float(np.max(np.abs(np.linalg.eigvalsh(predicted)-np.linalg.eigvalsh(observed)))))


def candidate(table, size, box):
    rows, K0 = table["rows"], table["K0_eV_A2"]
    coefficients, conditioning = solve(rows,K0,size)
    model = dict(K0_eV_A2=K0,coefficients_eV_A2=coefficients.tolist(),calibration_box=box)
    residuals = []
    squared = []
    for row in rows:
        observed = np.asarray(row["K_material_eV_A2"])
        predicted = stiffness_tensor(row["strain"],model)
        residuals.append(dict(name=row["name"],**errors(observed,predicted)))
        squared.extend((predicted-observed)[UPPER]**2)
    loo = []
    for index,row in enumerate(rows):
        if row["family"] == "cubic_reference":
            continue  # K0 is an external accepted reference, not a fitted intercept.
        try:
            c, condition = solve(rows[:index]+rows[index+1:],K0,size)
        except ValueError as exc:
            loo.append(dict(name=row["name"],identified=False,reason=str(exc)))
            continue
        held_model = dict(model,coefficients_eV_A2=c.tolist())
        predicted = stiffness_tensor(row["strain"],held_model,extrapolate=True)
        loo.append(dict(name=row["name"],identified=True,**errors(np.asarray(row["K_material_eV_A2"]),predicted)))
    return dict(model,degree=1 if size==3 else 2,basis_names=list(NAMES[:size]),conditioning=conditioning,
                rms_independent_entry_error_eV_A2=float(np.sqrt(np.mean(squared))),
                max_abs_entry_error_eV_A2=max(r["max_abs_entry_eV_A2"] for r in residuals),
                residuals=residuals,leave_one_case_out=loo)


def main():
    raw = Path("strain-table.json").read_bytes()
    table = json.loads(raw)
    E = np.asarray([row["strain"] for row in table["rows"]])
    diagonal = np.diagonal(E,axis1=-2,axis2=-1)
    off = E[:,[0,0,1],[1,2,2]]
    box = dict(normal_min=float(diagonal.min()),normal_max=float(diagonal.max()),
               shear_abs_max=float(np.max(np.abs(off))))
    linear, quadratic = candidate(table,3,box), candidate(table,12,box)
    # The 0.01 target selects complexity; actual achieved errors remain explicit.
    chosen = linear if linear["max_abs_entry_error_eV_A2"] <= .01 else quadratic
    model = dict(chosen,reference_a_A=4.084,material="PdD",cutoff_eV=500,kmesh=12,
                 frame="Material-frame Biot strain U-I; lab tensor is R K R-transpose",
                 operator="Collective Gamma, not local Kii; local contrast requires frozen sum(J)",
                 calibration_families="Hydro/uniaxial +/-1%; shear100 0..8% and one110 8%; biaxial -4..4% with c_z zero-stress relaxation",
                 applicability="Bounding box is necessary, not certification of arbitrary mixed strains. Extrapolation must be requested explicitly. Bare eigenvalues are not clipped.",
                 selection="Linear if maximum entry error <=0.01eV/A2, otherwise quadratic; no universal accuracy claim",
                 alternatives=[dict(degree=c["degree"],max_abs_entry_error_eV_A2=c["max_abs_entry_error_eV_A2"],
                                    rms_independent_entry_error_eV_A2=c["rms_independent_entry_error_eV_A2"]) for c in (linear,quadratic)],
                 source_table=dict(path="strain-table.json",bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    Path("strain-law.json").write_text(json.dumps(model,indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(degree=model["degree"],coefficients=dict(zip(model["basis_names"],model["coefficients_eV_A2"])),
                         max_abs_entry_error_eV_A2=model["max_abs_entry_error_eV_A2"],
                         rms_independent_entry_error_eV_A2=model["rms_independent_entry_error_eV_A2"],
                         conditioning=model["conditioning"],leave_one_case_out=model["leave_one_case_out"]),indent=2))


if __name__ == "__main__":
    main()
