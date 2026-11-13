"""Collect the entire requested tensor table plus accepted cubic normal stencils."""
import hashlib
import json
from pathlib import Path

import numpy as np
from mesh_readout import load_case

TENSOR_NAMES = (
    "g3_shear100_g002", "g3_shear100_g004", "g3_shear100_g006", "g3_shear100_g008",
    "g3_shear110_g110008", "g3_biax_m0p0400", "g3_biax_m0p0200",
    "g3_biax_0p0000", "g3_biax_0p0200", "g3_biax_0p0400",
)


def saved(path):
    raw = path.read_bytes()
    return json.loads(raw), dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def normal_rows():
    base, source = load_case("Pd", 12)
    K0 = base["K_Gamma_eV_A2"]
    sources = [source]
    rows = [dict(name="cubic_zero", family="cubic_reference", strain=np.zeros((3,3)).tolist(),
                 K_material_eV_A2=(K0*np.eye(3)).tolist(), construction="Cubic equivalence of the measured x column")]
    for side in ("minus", "plus"):
        hydro, source = load_case("Pd", 12, case="hydro_"+side)
        sources.append(source)
        rows.append(dict(name="hydro_"+side, family="hydrostatic",
                         strain=(np.asarray(hydro["deformation"])-np.eye(3)).tolist(),
                         K_material_eV_A2=(hydro["K_Gamma_eV_A2"]*np.eye(3)).tolist(),
                         construction="Cubic equivalence of the measured x column"))
        parallel, sp = load_case("Pd", 12, case="parallel_"+side)
        perpendicular, st = load_case("Pd", 12, case="perpendicular_"+side)
        sources.extend((sp,st))
        rows.append(dict(name="uniaxial_"+side, family="uniaxial_normal",
                         strain=(np.asarray(parallel["deformation"])-np.eye(3)).tolist(),
                         K_material_eV_A2=np.diag([parallel["K_Gamma_eV_A2"],
                             perpendicular["K_Gamma_eV_A2"],perpendicular["K_Gamma_eV_A2"]]).tolist(),
                         construction="Measured parallel x/perpendicular y; z=y by residual tetragonal symmetry"))
    return K0, rows, sources


def collect():
    K0, rows, sources = normal_rows()
    for name in TENSOR_NAMES:
        data, source = saved(Path("tensor-results")/(name+".json"))
        if data["spec"]["name"] != name:
            raise ValueError("Wrong saved tensor case")
        for binding in data["source_files"]:
            raw = Path(binding["path"]).read_bytes()
            if len(raw) != binding["bytes"] or hashlib.sha256(raw).hexdigest() != binding["sha256"]:
                raise ValueError("Changed native source column: "+binding["path"])
        sources.extend([source]+data["source_files"])
        E, K = np.asarray(data["Biot_strain"]), np.asarray(data["K_material_eV_A2"])
        if E.shape != (3,3) or K.shape != (3,3) or not np.isfinite(E).all() or not np.isfinite(K).all():
            raise ValueError("Unsupported tensor layout")
        rows.append(dict(name=name, family=data["spec"]["family"], strain=E.tolist(),
                         K_material_eV_A2=K.tolist(), spec=data["spec"],
                         max_antisymmetry_eV_A2=data["max_antisymmetry_eV_A2"],
                         construction="All three actual central displacement columns"))
    result = dict(K0_eV_A2=K0, rows=rows, source_files=list({s["path"]:s for s in sources}.values()),
                  strain="Material-frame Biot U-I", operator="Fixed-metal collective Gamma central force Jacobian")
    Path("strain-table.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    return result


if __name__ == "__main__":
    table = collect()
    for row in table["rows"]:
        print(row["name"], np.linalg.eigvalsh(row["K_material_eV_A2"]).tolist())
