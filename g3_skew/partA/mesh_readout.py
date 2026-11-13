"""Consume every requested PW500 mesh point; accept a stable 5% adjacent-mesh tail."""
import hashlib
import json
from pathlib import Path

import numpy as np

MESHES = (6, 8, 12, 16)


def require_profile(data, metal, mesh, cutoff, name, axis=0, deformation=None):
    expected = dict(metal=metal, kmesh=mesh, cutoff_eV=cutoff, name=name,
                    displacement_axis=axis, step_A=0.01, smearing_eV=0.1)
    for key, value in expected.items():
        if data[key] != value:
            raise ValueError(f"Wrong response profile {key}: {data[key]!r}; expected {value!r}")
    if deformation is None:
        deformation = np.eye(3)
    if not np.array_equal(data["deformation"], deformation):
        raise ValueError("Wrong imposed deformation for this response")


def load_case(metal, mesh, cutoff=500, case=None):
    name = f"k{mesh:02d}_pw{cutoff}"+("_"+case if case else "")
    path = Path("primitive_"+metal)/name/"result.json"
    raw = path.read_bytes()
    data = json.loads(raw)
    axis = 1 if case and case.startswith("perpendicular") else 0
    if case:
        eps = -0.01 if case.endswith("minus") else 0.01
        deformation = np.eye(3)*(1+eps) if case.startswith("hydro") else np.diag([1+eps, 1., 1.])
    else:
        deformation = np.eye(3)
    require_profile(data, metal, mesh, cutoff, name, axis, deformation)
    source = dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
    return data, source


def main():
    materials, sources = [], []
    for metal in ("Ni", "Pd"):
        cases = []
        for mesh in MESHES:
            case, source = load_case(metal, mesh)
            cases.append(case)
            sources.append(source)
        K12, K16 = cases[-2]["K_Gamma_eV_A2"], cases[-1]["K_Gamma_eV_A2"]
        if K16 == 0 or abs(K12-K16)/abs(K16) > 0.05:
            case, source = load_case(metal, 20)
            cases.append(case)
            sources.append(source)
        adjacent = []
        for left, right in zip(cases[:-1], cases[1:]):
            K = right["K_Gamma_eV_A2"]
            difference = abs(left["K_Gamma_eV_A2"]-K)/abs(K) if K != 0 else None
            adjacent.append(dict(mesh=left["kmesh"], next_mesh=right["kmesh"],
                                 relative_K_difference=difference,
                                 passes=difference is not None and difference <= 0.05))
        index = next((i for i in range(len(adjacent)) if all(p["passes"] for p in adjacent[i:])), None)
        accepted = cases[index] if index is not None else None
        frequency = accepted["harmonic_quantum_meV"] if accepted else None
        measured = cases[0]["measured_peak_meV"]
        mismatch = abs(frequency/measured-1) if frequency is not None else None
        materials.append(dict(metal=metal, status="accepted" if accepted else "unconverged",
                              accepted_kmesh=accepted["kmesh"] if accepted else None,
                              accepted_K_Gamma_eV_A2=accepted["K_Gamma_eV_A2"] if accepted else None,
                              accepted_quantum_meV=frequency, measured_peak_meV=measured,
                              relative_frequency_mismatch=mismatch,
                              PW600_check_required=mismatch is not None and mismatch > 0.20,
                              adjacent_comparisons=adjacent, cases=cases))
    result = dict(criterion="5% K agreement with next denser mesh; reject an apparent plateau contradicted by later requested meshes.",
                  method="PW500, PBE PAW10, Fermi0.1eV, Gamma-centred meshes, +/-0.01A H displacement; frozen metal.",
                  materials=materials, source_files=sources)
    Path("convergence.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    for material in materials:
        print(material["metal"], material["status"], "mesh", material["accepted_kmesh"],
              "K", material["accepted_K_Gamma_eV_A2"], "meV", material["accepted_quantum_meV"],
              "PW600?", material["PW600_check_required"])
        for case in material["cases"]:
            print(case["kmesh"], case["K_Gamma_eV_A2"], case["harmonic_quantum_meV"])


if __name__ == "__main__":
    main()
