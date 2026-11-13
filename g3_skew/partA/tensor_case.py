"""Run one PdD central-force column at an immutable homogeneous deformation."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from ase.build import bulk
from primitive import A, calculate


def run(spec_path, axis):
    raw = spec_path.read_bytes()
    spec = json.loads(raw)
    F = np.asarray(spec["deformation"], dtype=float)
    if F.shape != (3, 3) or not np.isfinite(F).all() or np.linalg.det(F) <= 0:
        raise ValueError("Expected a finite orientation-preserving deformation")
    name = spec["name"]
    if not name or not all(c.isalnum() or c == "_" for c in name):
        raise ValueError("Case name must be a single identifier")
    name += f"_axis{axis}"
    folder = Path("primitive_Pd")/name
    folder.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).read_bytes()
    atoms = bulk("PdH", "rocksalt", a=A["Pd"])
    meta = dict(input_bytes=len(raw), input_sha256=hashlib.sha256(raw).hexdigest(),
                driver_bytes=len(source), driver_sha256=hashlib.sha256(source).hexdigest(),
                displacement_axis=axis, moving_index=int(np.flatnonzero(atoms.numbers == 1)[0]),
                atomic_numbers=atoms.numbers.tolist())
    (folder/"tensor-meta.json").write_text(json.dumps(meta, indent=2)+"\n")
    result = calculate("Pd", name, F, kmesh=12, cutoff_eV=500, axis=axis, step_A=0.01)
    print("TENSOR_COLUMN_DONE", name, result["K_Gamma_eV_A2"], flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("spec", type=Path)
    parser.add_argument("axis", type=int, choices=(0, 1, 2))
    args = parser.parse_args()
    run(args.spec, args.axis)
