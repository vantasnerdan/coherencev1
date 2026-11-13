"""Relax only PdD c_z at fixed biaxial in-plane strain and fractional sites."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from ase.build import bulk
from scipy.optimize import brentq

from gpaw_model import calculator
from primitive import A

STRESS_TOL = 0.0005  # eV/A3, approximately 0.08 GPa


def relax(inplane):
    token = f"{inplane:.4f}".replace("-", "m").replace(".", "p")
    name = "g3_biax_"+token
    folder = Path("film_Pd")/name
    folder.mkdir(parents=True, exist_ok=False)
    sources = {}
    for filename in ("film_relax.py", "gpaw_model.py", "primitive.py"):
        raw = Path(filename).read_bytes()
        sources[filename] = dict(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
    records, cache = [], {}

    def stress(cz):
        if cz in cache:
            return cache[cz]["stress_zz_eV_A3"]
        trial = folder/f"trial{len(records):02d}"
        trial.mkdir()
        F = np.diag([1+inplane, 1+inplane, cz])
        atoms = bulk("PdH", "rocksalt", a=A["Pd"])
        atoms.set_cell(atoms.cell.array@F.T, scale_atoms=True)
        atoms.calc = calculator("Pd", 12, trial, 500)
        start, cpu = time.perf_counter(), time.process_time()
        print("FILM_SCF_START", name, cz, flush=True)
        energy = atoms.get_potential_energy(force_consistent=True)
        sigma = atoms.get_stress(voigt=False)
        row = dict(cz_ratio=cz, deformation=F.tolist(), energy_free_eV=float(energy),
                   stress_eV_A3=sigma.tolist(), stress_zz_eV_A3=float(sigma[2, 2]),
                   scf_cpu_s=time.process_time()-cpu, scf_elapsed_s=time.perf_counter()-start)
        atoms.calc.write(str(trial/"state.gpw"), mode="all")
        records.append(row)
        cache[cz] = row
        with (folder/"points.jsonl").open("a") as stream:
            stream.write(json.dumps(row, allow_nan=False)+"\n")
        print("FILM_SCF_DONE", name, cz, row["stress_zz_eV_A3"], flush=True)
        return row["stress_zz_eV_A3"]

    initial = stress(1.0)
    if abs(initial) > STRESS_TOL:
        for width in (0.04, 0.08, 0.16, 0.24):
            lo, hi = (1-width, 1.0) if initial > 0 else (1.0, 1+width)
            slo, shi = stress(lo), stress(hi)
            if slo <= 0 <= shi:
                brentq(stress, lo, hi, xtol=1e-4, rtol=1e-8, maxiter=12)
                break
        else:
            raise ValueError("No stable zero-stress c_z bracket inside the bounded search")
    best = min(records, key=lambda row: abs(row["stress_zz_eV_A3"]))
    if abs(best["stress_zz_eV_A3"]) > STRESS_TOL:
        raise ValueError("Out-of-plane stress has not reached the declared tolerance")
    spec = dict(name=name, family="biaxial_relaxed", inplane_strain=inplane,
                deformation=best["deformation"], reference_a_A=A["Pd"],
                cz_ratio=best["cz_ratio"], stress_zz_eV_A3=best["stress_zz_eV_A3"],
                constraint="Fixed x/y lengths and reference fractional sites; c_z zero-stress stationary state, not a global phase search.")
    output = dict(spec=spec, points=records, source_files=sources,
                  profile=dict(metal="Pd", kmesh=12, cutoff_eV=500, smearing_eV=0.1),
                  stress_tolerance_eV_A3=STRESS_TOL)
    (folder/"result.json").write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    destination = Path("tensor-specs")/(name+".json")
    destination.parent.mkdir(exist_ok=True)
    destination.write_text(json.dumps(spec, indent=2, allow_nan=False)+"\n")
    print("FILM_RELAXED", name, best["cz_ratio"], best["stress_zz_eV_A3"], flush=True)
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inplane", type=float, choices=(-0.04, -0.02, 0.0, 0.02, 0.04))
    args = parser.parse_args()
    relax(args.inplane)
