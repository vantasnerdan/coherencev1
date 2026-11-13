"""Collective primitive H/D energy/force grid at the accepted electronic sampling.

100 includes the shared zero. 111 uses the same zero from 100 in the readout.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import traceback

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers
from gpaw import GPAW
from gpaw_model import calculator
from mesh_readout import load_case


def restart_plus(base, path, folder):
    calc = GPAW(str(path), txt=str(folder/"gpaw.out"))
    atoms = calc.get_atoms()
    atoms.calc = calc
    expected = np.array(base["reference_positions_A"])
    expected[1, 0] += 0.01
    cell = np.array(base["cell_A"])
    delta = (atoms.positions-expected) @ np.linalg.inv(cell)
    delta -= np.rint(delta)
    if not np.array_equal(atoms.numbers, [atomic_numbers[base["metal"]], 1]) or not np.allclose(atoms.cell.array, cell, rtol=0, atol=1e-10) or np.max(np.abs(delta @ cell)) > 1e-10:
        raise ValueError("Checkpoint is not the accepted plus-displacement state")
    point = base["points"][1]
    if abs(atoms.get_potential_energy(force_consistent=True)-point["energy_free_eV"]) > 1e-6 or not np.allclose(atoms.get_forces(), point["forces_eV_A"], rtol=0, atol=1e-6):
        raise ValueError("Checkpoint energy/forces disagree with the bound response record")
    return atoms


def sample(metal, mesh, direction, strain_case):
    folder = Path("collective_"+metal)/f"k{mesh}_{strain_case}_{direction}"
    folder.mkdir(parents=True, exist_ok=True)
    base, source = load_case(metal, mesh, case="hydro_plus" if strain_case == "hydro_plus" else None)
    code_sources = []
    for name in ("collective.py", "gpaw_model.py"):
        raw = Path(name).read_bytes()
        code_sources.append(dict(path=name, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    reference = np.array(base["reference_positions_A"])
    vector = np.array([1., 0., 0.]) if direction == "100" else np.ones(3)/np.sqrt(3)
    if direction == "100":
        state = f"k{mesh:02d}_pw500"+("_hydro_plus" if strain_case == "hydro_plus" else "")
        checkpoint = Path("primitive_"+metal)/state/"plus.gpw"
        atoms = restart_plus(base, checkpoint, folder)
        with checkpoint.open("rb") as stream:
            checkpoint_source = dict(path=str(checkpoint), bytes=checkpoint.stat().st_size,
                                     sha256=hashlib.file_digest(stream, "sha256").hexdigest())
        displacements = (0., 0.06, 0.12, 0.18, 0.24, 0.30)
    else:
        atoms = Atoms(numbers=[atomic_numbers[metal], 1], positions=reference, cell=base["cell_A"], pbc=True)
        atoms.calc = calculator(metal, mesh, folder)
        checkpoint_source = None
        displacements = (0.06, 0.12, 0.18, 0.24, 0.30)
    points = []
    for u in displacements:
        positions = reference.copy()
        positions[1] += u*vector
        atoms.set_positions(positions)
        print(f"COLLECTIVE_SCF_START {metal} k{mesh} {direction} u_A={u}", flush=True)
        start, cpu0 = time.perf_counter(), time.process_time()
        energy = atoms.get_potential_energy(force_consistent=True)
        forces = atoms.get_forces()
        point = dict(u_A=u, energy_free_eV=float(energy), force_parallel_eV_A=float(forces[1] @ vector),
                     forces_eV_A=forces.tolist(), elapsed_s=time.perf_counter()-start, cpu_s=time.process_time()-cpu0)
        with (folder/"points.jsonl").open("a") as stream:
            stream.write(json.dumps(point, allow_nan=False)+"\n")
        points.append(point)
        atoms.calc.write(str(folder/f"u{u:.2f}.gpw"), mode="all")
        print(f"COLLECTIVE_SCF_DONE {metal} k{mesh} {direction} u_A={u} E={energy}", flush=True)
    result = dict(metal=metal, kmesh=mesh, cutoff_eV=500, smearing_eV=0.1, direction=direction,
                  strain_case=strain_case, deformation=base["deformation"], K2_force_eV_A2=base["K_Gamma_eV_A2"],
                  displacement_unit_vector=vector.tolist(), base_source=source, checkpoint_source=checkpoint_source,
                  definition="Energy per primitive cell/per H or D for uniform collective displacement; fixed metal, self-consistent electrons.",
                  zero_source="Shared zero is the u=0 point in the 100 grid.", points=points, code_sources=code_sources)
    (folder/"result.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps(result, indent=2, allow_nan=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metal", choices=("Ni", "Pd"))
    parser.add_argument("kmesh", type=int, choices=(8, 12))
    parser.add_argument("direction", choices=("100", "111"))
    parser.add_argument("strain_case", choices=("zero", "hydro_plus"))
    args = parser.parse_args()
    try:
        sample(args.metal, args.kmesh, args.direction, args.strain_case)
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
