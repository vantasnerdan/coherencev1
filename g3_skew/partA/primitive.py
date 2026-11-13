"""Finite-displacement primitive hydride response at an explicit electronic mesh."""
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from ase.build import bulk
from gpaw_model import calculator

A = {"Ni": 3.732, "Pd": 4.084}


def calculate(metal, name, deformation, *, kmesh, cutoff_eV=500, axis=0, step_A=0.01):
    folder = Path("primitive_"+metal)/name
    folder.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).read_bytes()
    calculator_source = Path("gpaw_model.py").read_bytes()
    atoms = bulk(metal+"H", "rocksalt", a=A[metal])
    atoms.set_cell(atoms.cell.array @ deformation.T, scale_atoms=True)
    reference = atoms.positions.copy()
    moving = int(np.flatnonzero(atoms.numbers == 1)[0])
    atoms.calc = calculator(metal, kmesh, folder, cutoff_eV)
    records = []
    for displacement in (-step_A, step_A):
        xyz = reference.copy()
        xyz[moving, axis] += displacement
        atoms.set_positions(xyz)
        print(f"PRIMITIVE_SCF_START {metal} {name} u_A={displacement}", flush=True)
        start, cpu0 = time.perf_counter(), time.process_time()
        energy = atoms.get_potential_energy(force_consistent=True)
        forces = atoms.get_forces()
        row = dict(displacement_A=displacement, energy_free_eV=float(energy),
                   forces_eV_A=forces.tolist(), force_parallel_eV_A=float(forces[moving, axis]),
                   elapsed_s=time.perf_counter()-start, cpu_s=time.process_time()-cpu0)
        with (folder/"points.jsonl").open("a") as stream:
            stream.write(json.dumps(row, allow_nan=False)+"\n")
        records.append(row)
        atoms.calc.write(str(folder/("minus.gpw" if displacement < 0 else "plus.gpw")), mode="all")
        print(f"PRIMITIVE_SCF_DONE {metal} {name} u_A={displacement} Faxis={row['force_parallel_eV_A']}", flush=True)
    K = -(records[1]["force_parallel_eV_A"]-records[0]["force_parallel_eV_A"])/(2*step_A)
    result = dict(metal=metal, name=name, K_Gamma_eV_A2=float(K), deformation=deformation.tolist(),
                  reference_positions_A=reference.tolist(), cell_A=atoms.cell.array.tolist(),
                  displacement_axis=axis, cutoff_eV=cutoff_eV, kmesh=kmesh, step_A=step_A,
                  smearing_eV=0.1, point_group="Actual displaced-cell subgroup; not undisplaced crystal symmetry.",
                  n_ibz_kpoints=len(atoms.calc.get_ibz_k_points()),
                  producer_bytes=len(source), producer_sha256=hashlib.sha256(source).hexdigest(), points=records,
                  calculator_bytes=len(calculator_source), calculator_sha256=hashlib.sha256(calculator_source).hexdigest())
    (folder/"result.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    return result
