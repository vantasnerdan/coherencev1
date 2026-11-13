"""Zero-pressure centred cubic PdD reference and its matched collective stiffness."""
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from ase.build import bulk
from scipy.optimize import brentq

from gpaw_model import calculator
from primitive import A, calculate
from mesh_readout import load_case


def binding(path):
    raw = Path(path).read_bytes()
    return dict(path=str(path), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def main():
    folder = Path("cubic_Pd")
    folder.mkdir(exist_ok=False)
    source_paths = ["cubic_spacing.py", "gpaw_model.py", "primitive.py",
                    "film_Pd/g3_biax_0p0000/result.json"]
    sources = [binding(path) for path in source_paths]
    film = json.loads(Path(source_paths[-1]).read_text())
    initial = film["points"][0]
    if film["profile"] != dict(metal="Pd", kmesh=12, cutoff_eV=500, smearing_eV=0.1):
        raise ValueError("Wrong reused zero-strain stress profile")
    if not np.array_equal(initial["deformation"], np.eye(3)):
        raise ValueError("Reused stress is not the cubic4.084A point")
    records = [dict(initial, scale=1.0, reused_from=source_paths[-1])]
    cache = {1.0: float(np.trace(initial["stress_eV_A3"])/3)}

    def stress(scale):
        if scale in cache:
            return cache[scale]
        trial = folder/f"trial{len(records):02d}"
        trial.mkdir()
        atoms = bulk("PdH", "rocksalt", a=A["Pd"]*scale)
        atoms.calc = calculator("Pd", 12, trial, 500)
        start = time.process_time()
        print("CUBIC_SCF_START", scale, flush=True)
        energy = atoms.get_potential_energy(force_consistent=True)
        sigma = atoms.get_stress(voigt=False)
        mean = float(np.trace(sigma)/3)
        row = dict(scale=scale, a_A=A["Pd"]*scale, energy_free_eV=float(energy),
                   stress_eV_A3=sigma.tolist(), mean_stress_eV_A3=mean,
                   scf_cpu_s=time.process_time()-start)
        atoms.calc.write(str(trial/"state.gpw"), mode="all")
        records.append(row)
        cache[scale] = mean
        with (folder/"points.jsonl").open("a") as stream:
            stream.write(json.dumps(row, allow_nan=False)+"\n")
        print("CUBIC_SCF_DONE", scale, mean, flush=True)
        return mean

    brentq(stress, 1.0, 1.04, xtol=5e-5, rtol=1e-8, maxiter=12)
    scale = min(cache, key=lambda value: abs(cache[value]))
    if abs(cache[scale]) > 0.0005:
        raise ValueError("Cubic stress did not reach the declared0.0005eV/A3 tolerance")
    gamma = calculate("Pd", "g3_cubic_pbe_equilibrium", np.eye(3)*scale,
                      kmesh=12, cutoff_eV=500, axis=0, step_A=0.01)
    base, base_source = load_case("Pd",12)
    output = dict(reference="Centred cubic zero-pressure PBE free-energy reference; nuclei fixed fractionally, no nuclear ZPE",
                  reference_a_A=A["Pd"], PBE_cubic_a_A=A["Pd"]*scale,
                  PBE_scale=scale, mean_stress_eV_A3=cache[scale],
                  K_Gamma_at_PBE_a_eV_A2=gamma["K_Gamma_eV_A2"],
                  K_Gamma_at_4p084_eV_A2=base["K_Gamma_eV_A2"],
                  points=records, equilibrium_gamma=gamma,
                  source_files=sources+[base_source,binding("primitive_Pd/g3_cubic_pbe_equilibrium/result.json")])
    (folder/"result.json").write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    print("CUBIC_REFERENCE_DONE", output["PBE_cubic_a_A"], output["K_Gamma_at_PBE_a_eV_A2"], flush=True)


if __name__ == "__main__":
    main()
