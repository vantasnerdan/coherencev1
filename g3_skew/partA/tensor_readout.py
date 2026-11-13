"""Assemble the actual force columns, retaining finite-step and frame residuals."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def column(spec, raw, axis):
    folder = Path("primitive_Pd")/(spec["name"]+f"_axis{axis}")
    meta = json.loads((folder/"tensor-meta.json").read_text())
    result = json.loads((folder/"result.json").read_text())
    expected = dict(metal="Pd", kmesh=12, cutoff_eV=500, step_A=0.01,
                    displacement_axis=axis, smearing_eV=0.1)
    if any(result[k] != v for k, v in expected.items()):
        raise ValueError(f"Wrong material/numerical profile in {folder}")
    if meta["input_sha256"] != hashlib.sha256(raw).hexdigest() or meta["input_bytes"] != len(raw):
        raise ValueError(f"Wrong deformation input binding in {folder}")
    if not np.array_equal(result["deformation"], spec["deformation"]):
        raise ValueError(f"Wrong physical deformation in {folder}")
    moving = meta["moving_index"]
    if meta["atomic_numbers"][moving] != 1:
        raise ValueError("Force column does not select hydrogen")
    minus, plus = result["points"]
    if [minus["displacement_A"], plus["displacement_A"]] != [-0.01, 0.01]:
        raise ValueError("Wrong central-displacement ordering")
    K = -(np.asarray(plus["forces_eV_A"])[moving]-np.asarray(minus["forces_eV_A"])[moving])/0.02
    return K, dict(path=str(folder/"result.json"), bytes=(folder/"result.json").stat().st_size,
                   sha256=hashlib.sha256((folder/"result.json").read_bytes()).hexdigest())


def read_tensor(spec_path):
    raw = spec_path.read_bytes()
    spec = json.loads(raw)
    columns = [column(spec, raw, axis) for axis in range(3)]
    Kraw = np.column_stack([value for value, source in columns])
    if not np.isfinite(Kraw).all():
        raise ValueError("Nonfinite force Jacobian")
    antisym = float(np.max(np.abs(Kraw-Kraw.T)))
    if antisym > 0.005:
        raise ValueError(f"Unresolved finite-step Hessian antisymmetry: {antisym} eV/A2")
    K = (Kraw+Kraw.T)/2
    F = np.asarray(spec["deformation"])
    left, singular, rightT = np.linalg.svd(F)
    R = left@rightT
    U = rightT.T@np.diag(singular)@rightT
    material = R.T@K@R
    residuals = dict(polar=float(np.max(np.abs(R@U-F))),
                     frame_roundtrip=float(np.max(np.abs(R@material@R.T-K))))
    if max(residuals.values()) > 1e-12 or np.linalg.det(R) <= 0:
        raise ValueError(f"Unresolved polar/frame convention: {residuals}")
    output = dict(spec=spec, K_raw_lab_eV_A2=Kraw.tolist(), K_lab_eV_A2=K.tolist(),
                  K_material_eV_A2=material.tolist(), eigenvalues_eV_A2=np.linalg.eigvalsh(K).tolist(),
                  R=R.tolist(), U=U.tolist(), Biot_strain=(U-np.eye(3)).tolist(),
                  max_antisymmetry_eV_A2=antisym, frame_residuals=residuals,
                  source_files=[source for value, source in columns],
                  convention="Collective Gamma central force Jacobian, fixed metal, self-consistent electrons; not local Kii.")
    destination = Path("tensor-results")/(spec["name"]+".json")
    destination.parent.mkdir(exist_ok=True)
    destination.write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spec", type=Path)
    args = parser.parse_args()
    print(json.dumps(read_tensor(args.spec), indent=2, allow_nan=False))
