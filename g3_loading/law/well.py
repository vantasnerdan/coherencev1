"""Single-nucleus NiH/PdD potentials; PBE PAW/PW400, frozen host.

python well.py Ni
python well.py Pd
"""
import argparse
import json
from pathlib import Path
import resource
import sys
import time
import traceback

import numpy as np
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac, Mixer
from scipy.constants import e, epsilon_0

MATERIALS = {"Ni": 3.732, "Pd": 4.084}
DIRECTIONS = {"100": np.array([1., 0., 0.]), "111": np.ones(3) / np.sqrt(3)}
SMEAR_EV = 0.1


def cpu_s():
    u = resource.getrusage(resource.RUSAGE_SELF)
    return u.ru_utime + u.ru_stime


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def electronic_inputs(calc, volume_A3):
    ef = calc.get_fermi_level()
    weights = calc.get_k_point_weights()
    dos = 0.0
    for kpt, weight in enumerate(weights):
        bands = calc.get_eigenvalues(kpt=kpt, spin=0)
        f = 1 / (1 + np.exp((bands - ef) / SMEAR_EV))
        dos += 2 * weight * np.sum(f * (1 - f)) / SMEAR_EV
    q2 = e * dos / (epsilon_0 * volume_A3 * 1e-30)
    return dict(valence_electrons=float(calc.get_number_of_electrons()),
                valence_electron_density_m3=float(calc.get_number_of_electrons() / (volume_A3 * 1e-30)),
                Fermi_energy_eV=float(ef), thermal_DOS_states_per_eV_cell=float(dos),
                Thomas_Fermi_length_A=float(1e10 / np.sqrt(q2)) if q2 > 0 else None,
                screening_definition="Independent-electron finite-T compressibility, dn/dmu from0.1eV Fermi occupations; not a measured dielectric length.",
                conduction_carrier_density_m3=None, dynamical_effective_charge=None)


def baseline(metal):
    folder = Path(metal)
    folder.mkdir(exist_ok=True)
    atoms = bulk(metal + "H", "rocksalt", a=MATERIALS[metal], cubic=True).repeat((2, 2, 2))
    atoms.calc = GPAW(mode=PW(400), setups={metal: "10"}, xc="PBE",
                      kpts={"size": (1, 1, 1), "gamma": True}, nbands=-20,
                      occupations=FermiDirac(SMEAR_EV), spinpol=False,
                      symmetry={"point_group": False}, mixer=Mixer(0.05, 5, 50),
                      maxiter=80, convergence={"energy": 1e-7, "density": 1e-6},
                      txt=str(folder / "base.out"))
    hydrogen = np.flatnonzero(atoms.numbers == 1)
    centre = atoms.cell.sum(axis=0) / 2
    moving = int(hydrogen[np.argmin(np.linalg.norm(atoms.positions[hydrogen] - centre, axis=1))])
    started, cpu0 = time.perf_counter(), cpu_s()
    print(f"BASE_START {metal} atoms={len(atoms)} moving={moving}", flush=True)
    energy = atoms.get_potential_energy(force_consistent=True)
    forces = atoms.get_forces()
    record = dict(metal=metal, a_A=MATERIALS[metal], atoms=len(atoms), moving_index=moving,
                  cell_A=atoms.cell.tolist(), reference_positions_A=atoms.positions.tolist(),
                  energy_free_eV=float(energy), force_H_eV_A=forces[moving].tolist(),
                  atomic_numbers=atoms.numbers.tolist(),
                  elapsed_s=time.perf_counter() - started, cpu_s=cpu_s() - cpu0,
                  electronic=electronic_inputs(atoms.calc, atoms.get_volume()))
    save(folder / "base.json", record)
    atoms.calc.write(str(folder / "base.gpw"), mode="all")
    print(f"BASE_DONE {metal} cpu_s={record['cpu_s']} wall_s={record['elapsed_s']}", flush=True)
    return record


def points(metal, kind):
    folder = Path(metal)
    base = json.loads((folder / "base.json").read_text())
    calc = GPAW(str(folder / "base.gpw"), txt=str(folder / f"{kind}.out"))
    atoms = calc.get_atoms()
    atoms.calc = calc
    reference = np.array(base["reference_positions_A"])
    moving = base["moving_index"]
    if kind == "coupling":
        cases = [("100", -0.03), ("100", 0.03)]
    else:
        cases = [(name, u) for name in DIRECTIONS for u in (0.12, 0.24)]
    for number, (name, u) in enumerate(cases):
        direction = DIRECTIONS[name]
        xyz = reference.copy()
        xyz[moving] += u * direction
        atoms.set_positions(xyz)
        started, cpu0 = time.perf_counter(), cpu_s()
        print(f"POINT_START {metal} {kind} {name} u_A={u}", flush=True)
        energy = atoms.get_potential_energy(force_consistent=True)
        force = atoms.get_forces()
        row = dict(metal=metal, kind=kind, direction=name, displacement_A=u,
                   energy_free_eV=float(energy), delta_energy_eV=float(energy-base["energy_free_eV"]),
                   force_H_eV_A=force[moving].tolist(), force_parallel_eV_A=float(force[moving] @ direction),
                   net_force_eV_A=force.sum(axis=0).tolist(),
                   elapsed_s=time.perf_counter()-started, cpu_s=cpu_s()-cpu0)
        if kind == "coupling":
            hydrogen = np.flatnonzero(atoms.numbers == 1)
            row["H_indices"] = hydrogen.tolist()
            row["forces_H_eV_A"] = force[hydrogen].tolist()
        with (folder / "points.jsonl").open("a") as stream:
            stream.write(json.dumps(row, allow_nan=False) + "\n")
        print(f"POINT_DONE {metal} {kind} {name} u_A={u} dE_eV={row['delta_energy_eV']} cpu_s={row['cpu_s']}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metal", choices=tuple(MATERIALS))
    metal = parser.parse_args().metal
    baseline(metal)
    points(metal, "coupling")
    print(f"G2_{metal}_COUPLING_COMPLETED", flush=True)
    points(metal, "well")
    print(f"G2_{metal}_SAMPLING_COMPLETED", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        save(f"failure-{sys.argv[1]}.json", dict(type=type(error).__name__, message=str(error)))
        traceback.print_exc()
        raise SystemExit(1)
