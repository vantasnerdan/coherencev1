"""Single authorized GPAW route: at most three PdD SCFs, no retry.

Two-atom rocksalt cell, a=4.084 A, PBE PAW/PW500 eV, 6x6x6 k mesh.
D is represented electronically by H with atomic mass2.014101778 u.
Displace every periodic D by +/-0.02 A along Cartesian cubic x, Pd fixed.
The Gamma optic uses the Pd/D reduced mass, not the D mass alone.
"""
import json
import math
from pathlib import Path
import resource
import time
import traceback

import numpy as np
from ase import Atoms
from gpaw import GPAW, PW, FermiDirac

A = 4.084
DELTA_A = 0.02
M_PD = 106.42
M_D = 2.014101778
AMU_KG = 1.66053906892e-27
EV_J = 1.602176634e-19
HBAR = 6.62607015e-34/(2*math.pi)


def cpu_seconds():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_utime+usage.ru_stime


def main():
    start = time.perf_counter()
    cpu_start = cpu_seconds()
    atoms = Atoms("PdH",scaled_positions=[(0,0,0),(0.5,0.5,0.5)],
                  cell=[(0,A/2,A/2),(A/2,0,A/2),(A/2,A/2,0)],pbc=True)
    atoms.set_masses([M_PD,M_D])
    atoms.calc = GPAW(mode=PW(500),xc="PBE",kpts={"size":(6,6,6),"gamma":True},
                      nbands=16,occupations=FermiDirac(0.1),spinpol=False,
                      symmetry={"point_group":False},maxiter=100,
                      convergence={"energy":1e-6,"density":1e-5,"eigenstates":1e-8},
                      txt="gpaw.out")
    reference = atoms.positions.copy()
    records = []
    for name, displacement in (("base",0.0),("plus",DELTA_A),("minus",-DELTA_A)):
        positions = reference.copy()
        positions[1,0] += displacement
        atoms.set_positions(positions)
        print(f"SCF_START {name} displacement_A={displacement}",flush=True)
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        row = dict(name=name,displacement_A=displacement,energy_eV=float(energy),
                   forces_eV_A=forces.tolist(),elapsed_s=time.perf_counter()-start,
                   cpu_s=cpu_seconds()-cpu_start)
        records.append(row)
        with Path("scf-results.jsonl").open("a") as stream:
            stream.write(json.dumps(row)+"\n")
        print(f"SCF_DONE {name} energy_eV={energy} force_D_x={forces[1,0]}",flush=True)
    curvature = -(records[1]["forces_eV_A"][1][0]-records[2]["forces_eV_A"][1][0])/(2*DELTA_A)
    reduced_u = M_PD*M_D/(M_PD+M_D)
    omega2 = curvature*EV_J/1e-20/(reduced_u*AMU_KG)
    magnitude_eV = HBAR*math.sqrt(abs(omega2))/EV_J
    output = dict(functional="PBE",mode="PW500eV",kpoints=[6,6,6],
                  electronic_smearing_eV=0.1,lattice_constant_A=A,
                  displacement_A=DELTA_A,masses_u=[M_PD,M_D],reduced_mass_u=reduced_u,
                  curvature_eV_A2=curvature,omega2_s2=omega2,
                  harmonic_optic_energy_eV=magnitude_eV if omega2>0 else None,
                  imaginary_energy_magnitude_eV=magnitude_eV if omega2<0 else None,
                  measured_peak_eV=0.037,
                  difference_from_peak_eV=magnitude_eV-0.037 if omega2>0 else None,
                  scf_count=len(records),scf_results=records,
                  elapsed_s=time.perf_counter()-start,cpu_s=cpu_seconds()-cpu_start)
    Path("result.json").write_text(json.dumps(output,indent=2)+"\n")
    print("GPAW_ROUTE_COMPLETED",flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        Path("failure.json").write_text(json.dumps({"type":type(error).__name__,
                                                    "message":str(error)},indent=2)+"\n")
        traceback.print_exc()
        raise SystemExit(1)
