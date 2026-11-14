"""Read completed single-nucleus GPAW points; fit V=k2*u²/2+k4*u⁴/4.

The harmonic coefficient and H-H Hessian come from the same small-step force
response. Larger energy/force rows constrain the directional quartic term.
"""
import hashlib
from importlib.metadata import version
import json
from pathlib import Path

import numpy as np
from scipy.constants import e, hbar, m_p, physical_constants


def fit_direction(rows, mass, k2):
    u = np.array([x["displacement_A"] for x in rows])
    energy = np.array([x["delta_energy_eV"] for x in rows])
    force = np.array([x["force_parallel_eV_A"] for x in rows])
    scale = max(abs(u))
    design_energy = np.c_[u*u/2, u**4/4]
    design_force = np.c_[u, u**3]
    matrix = np.vstack([design_energy, scale*design_force])
    target = np.r_[energy, -scale*force]
    quartic_column = matrix[:, 1]
    k4 = np.dot(quartic_column, target-matrix[:, 0]*k2) / np.dot(quartic_column, quartic_column)
    efit, ffit = design_energy @ [k2, k4], -design_force @ [k2, k4]
    quantum = hbar * np.sqrt(abs(k2)*e/1e-20/mass) / e
    return dict(k2_eV_A2=float(k2), k4_eV_A4=float(k4),
                coefficient_of_u2_eV_A2=float(k2/2), coefficient_of_u4_eV_A4=float(k4/4),
                harmonic_quantum_eV=float(quantum) if k2 > 0 else None,
                imaginary_harmonic_magnitude_eV=float(quantum) if k2 < 0 else None,
                energy_fit_RMS_eV=float(np.sqrt(np.mean((efit-energy)**2))),
                force_fit_RMS_eV_A=float(np.sqrt(np.mean((ffit-force)**2))),
                sampled_positive_range_A=[float(min(u)), float(max(u))],
                inversion_assumption="Ideal centrosymmetric rocksalt O site; V(-u)=V(u).")


def coupling(base, rows):
    pair = sorted((r for r in rows if r["kind"] == "coupling"), key=lambda r: r["displacement_A"])
    if [r["displacement_A"] for r in pair] != [-0.03, 0.03]:
        raise ValueError("Incomplete or duplicated small-step force-response pair")
    minus, plus = pair
    if minus["H_indices"] != plus["H_indices"]:
        raise ValueError("Hydrogen index order changed between force responses")
    indices = plus["H_indices"]
    step = plus["displacement_A"]-minus["displacement_A"]
    phi = -(np.array(plus["forces_H_eV_A"])-np.array(minus["forces_H_eV_A"]))/step
    positions = np.array(base["reference_positions_A"])[indices]
    reference = np.array(base["reference_positions_A"])[base["moving_index"]]
    lengths = np.diag(base["cell_A"])
    offsets = positions-reference
    offsets -= np.rint(offsets/lengths)*lengths
    distances = np.linalg.norm(offsets, axis=1)
    source = indices.index(base["moving_index"])
    shells = {}
    atom_rows = []
    for i, atom in enumerate(indices):
        atom_rows.append(dict(atom_index=atom, R_A=offsets[i].tolist(), distance_A=float(distances[i]),
                              Phi_ix_jx_eV_A2=float(phi[i, 0]), Phi_iy_jx_eV_A2=float(phi[i, 1]),
                              Phi_iz_jx_eV_A2=float(phi[i, 2]), pulling_J_xx_eV_A2=float(-phi[i, 0])))
        if i != source:
            shells.setdefault(round(float(distances[i]), 6), []).append(float(phi[i, 0]))
    shell_rows = [dict(distance_A=r, distinct_supercell_H_sites=len(values),
                       mean_Phi_xx_eV_A2=float(np.mean(values)),
                       sum_Phi_xx_eV_A2=float(sum(values)), sum_pulling_J_xx_eV_A2=float(-sum(values)),
                       min_Phi_xx_eV_A2=min(values), max_Phi_xx_eV_A2=max(values))
                  for r, values in sorted(shells.items())]
    total_phi = -(np.array(plus["net_force_eV_A"])-np.array(minus["net_force_eV_A"]))/step
    return dict(definition="Phi_i(alpha),j(x) = -[F_i(alpha)(+0.03A)-F_i(alpha)(-0.03A)]/0.06A; energy Hessian.",
                displaced_direction="100", fixed_metal=True, self_consistent_electrons=True,
                on_site_k2_eV_A2=float(phi[source, 0]),
                sum_neighbour_pulling_J_xx_eV_A2=float(-(phi[:, 0].sum()-phi[source, 0])),
                H_only_q0_stiffness_eV_A2=float(phi[:, 0].sum()),
                q0_definition="k2_onsite-sum(J_neighbours), J=-Phi_xx. H/D subspace only, metals fixed; not the all-ion acoustic sum.",
                H_column_sum_vector_eV_A2=phi.sum(axis=0).tolist(),
                metal_reaction_column_sum_eV_A2=(total_phi-phi.sum(axis=0)).tolist(),
                all_ion_acoustic_residual_eV_A2=total_phi.tolist(),
                shells=shell_rows, atoms=atom_rows)


def main():
    output = dict(convention="V(u)-V(0)=k2*u²/2+k4*u⁴/4; k2 from the small-step force Hessian.",
                  method="GPAW PBE PAW,10-valence Ni/Pd setups,PW400eV,Gamma-only,Fermi0.1eV,nonmagnetic,frozen host.",
                  isotope_rule="Local frozen-host vibration uses the proton/deuteron mass. Force constants already contain self-consistent electronic screening.",
                  software={name: version(name) for name in ("gpaw", "ase", "numpy", "scipy")}, materials=[])
    for metal, isotope in (("Ni", "H"), ("Pd", "D")):
        base = json.loads((Path(metal)/"base.json").read_text())
        rows = [json.loads(x) for x in (Path(metal)/"points.jsonl").read_text().splitlines()]
        response = coupling(base, rows)
        mass = m_p if isotope == "H" else physical_constants["deuteron mass"][0]
        directions = {}
        for direction in ("100", "111"):
            selected = [r for r in rows if r["kind"] == "well" and r["direction"] == direction]
            if sorted(r["displacement_A"] for r in selected) != [0.12, 0.24]:
                raise ValueError(f"Incomplete or duplicated {metal}/{direction} well data")
            directions[direction] = fit_direction(selected, mass, response["on_site_k2_eV_A2"])
        output["materials"].append(dict(material=metal+isotope, a_A=base["a_A"],
                                        supercell_atoms=base["atoms"], image_repeat_A=2*base["a_A"],
                                        directions=directions, coupling=response,
                                        electronic=base["electronic"], raw_points=rows))
    output["source_sha256"] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [Path("well.py"), Path("analyse.py"), *sorted(Path("setups").glob("*.gz"))]}
    Path("wells.json").write_text(json.dumps(output, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"materials": [{k:v for k,v in m.items() if k != "raw_points"} for m in output["materials"]]},
                     indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
