"""Read actual loading SCF records without importing GPAW or reading GPW/PAW bytes.

--config CONFIG --mesh 3 --probe clean|near [--output UNIQUE.json].
Reads one actual central pair along that probe's registered axis and shared zero.
Outputs registered-site tangents, not relaxed minima or quartic/well parameters.
"""
import argparse
from pathlib import Path
import uuid

import numpy as np

from dense_pd import ROOT, check_profile, commit_json, manifest, read_json, require
from loading_pd import (SCHEMA, cases, config, direction_vector, folder_for,
                        load_point, probe_for, profile)


def central_column(p, probe, direction, zero, minus, plus, bindings):
    g = p["geometry"]
    source = probe_for(p, probe)["moving_index"]
    numbers = np.array(g["atomic_numbers"])
    original = np.array(g["original_atom_indices"])
    vector = direction_vector(p, probe)
    h = 0.03
    require(minus["displacement_A"] == -h and plus["displacement_A"] == h
            and minus["direction"] == plus["direction"] == direction
            and minus["probe_id"] == plus["probe_id"] == probe, "Not a matching actual central pair")
    fm = np.array(minus["forces_eV_A"])
    fp = np.array(plus["forces_eV_A"])
    f0 = np.array(zero["forces_eV_A"])
    phi = -(fp - fm) / (2 * h)
    even = (fp + fm) / 2 - f0
    kii = float(phi[source] @ vector)
    data = dict(record_type="measured_central_force_column", direction=direction,
                probe=probe_for(p, probe), measured_points=bindings, independent_negative_scf=True,
                source_displacement_vector=vector.tolist(), h_A=h,
                definition="Phi_i(alpha),j(d)=-(F_i(alpha)(+h*d)-F_i(alpha)(-h*d))/(2h)",
                scope="One actual all-atom column; not a full Hessian or reciprocal bond law",
                K2_eV_A2=kii, F0_source_eV_A=f0[source].tolist(),
                F0_parallel_eV_A=float(f0[source] @ vector),
                F0_source_norm_eV_A=float(np.linalg.norm(f0[source])),
                K2_energy_central_eV_A2=float((plus["energy_free_eV"] + minus["energy_free_eV"]
                                               - 2 * zero["energy_free_eV"]) / h**2),
                energy_odd_slope_plus_F0_eV_A=float((plus["energy_free_eV"] - minus["energy_free_eV"])
                                                   / (2 * h) + f0[source] @ vector),
                force_evenness=dict(definition="(F(+h)+F(-h))/2-F(0); no inversion assumption",
                    source_vector_eV_A=even[source].tolist(), source_parallel_eV_A=float(even[source] @ vector),
                    max_all_atom_norm_eV_A=float(np.linalg.norm(even, axis=1).max()),
                    all_atom_vectors_eV_A=even.tolist()),
                all_atom_columns=[dict(atom_index=i, original_atom_index=int(original[i]),
                                      atomic_number=int(numbers[i]), Phi_vector_eV_A2=phi[i].tolist())
                                  for i in range(len(numbers))],
                all_ion_column_sum_eV_A2=phi.sum(axis=0).tolist(),
                all_ion_sum_scope="Acoustic diagnostic for the full-ion source column, not an H-only stiffness",
                stationarity="Registered-site finite-step tangent only; F0 is retained, no relaxed-minimum claim")
    data["coupling"] = hydrogen_column(p, source, phi, vector)
    return data


def hydrogen_column(p, source, phi, vector):
    g = p["geometry"]
    numbers = np.array(g["atomic_numbers"])
    indices = np.flatnonzero(numbers == 1)
    original = np.array(g["original_atom_indices"])
    positions = np.array(g["reference_positions_A"])
    cell = np.array(g["cell_A"])
    unwrapped = positions[indices] - positions[source]
    shifts = -np.rint(unwrapped @ np.linalg.inv(cell)).astype(int)
    offsets = unwrapped + shifts @ cell
    distances = np.linalg.norm(offsets, axis=1)
    hphi = phi[indices]
    atoms, shells = [], {}
    # Geometry was source-index checked by the producer's shared config loader.
    base = read_json(ROOT / "Pd/base.json")
    all_h = np.flatnonzero(np.array(base["atomic_numbers"]) == 1)
    base_xyz = np.array(base["reference_positions_A"]) * (p["a_A"] / base["a_A"])
    reference = base_xyz[all_h] - positions[source]
    reference -= np.rint(reference @ np.linalg.inv(cell)) @ cell
    reference_radii = np.linalg.norm(reference, axis=1)
    for i, atom in enumerate(indices):
        atoms.append(dict(atom_index=int(atom), original_atom_index=int(original[atom]),
                          is_source=int(atom) == source, R_reference_unwrapped_A=unwrapped[i].tolist(),
                          minimum_image_translation=shifts[i].tolist(), R_A=offsets[i].tolist(),
                          distance_A=float(distances[i]), radius_over_a=float(distances[i] / p["a_A"]),
                          Phi_i_alpha_jd_vector_eV_A2=hphi[i].tolist(),
                          Phi_dd_eV_A2=float(hphi[i] @ vector),
                          pulling_Jdd_eV_A2=None if atom == source else float(-hphi[i] @ vector)))
        if atom != source:
            shells.setdefault(round(float(distances[i] / p["a_A"]), 10), []).append(i)
    shell_rows = []
    # Include vacancy-only shells as empty inventories, never as invented bonds.
    radii = sorted(set(shells) | {round(float(r / p["a_A"]), 10) for r in reference_radii if r > 1e-8})
    vacancies = set(p["vacancy_original_indices"])
    for radius in radii:
        members = shells.get(radius, [])
        values = hphi[members]
        projected = values @ vector
        reference_members = [int(atom) for atom, r in zip(all_h, reference_radii)
                             if abs(r / p["a_A"] - radius) < 1e-9]
        occupied = indices[members]
        vacant = sorted(vacancies & set(reference_members))
        require(sorted(original[occupied].tolist()) == sorted(set(reference_members) - vacancies),
                "Shell source-index inventory mismatch")
        shell_rows.append(dict(radius_over_a=radius, distance_A=radius * p["a_A"],
            occupied_multiplicity=len(members), original_reference_multiplicity=len(reference_members),
            occupied_atom_indices=occupied.tolist(), occupied_original_indices=original[occupied].tolist(),
            vacant_original_indices=vacant, vacant_multiplicity=len(vacant),
            sum_Phi_vector_eV_A2=values.sum(axis=0).tolist(),
            mean_Phi_vector_eV_A2=values.mean(axis=0).tolist() if members else None,
            sum_Phi_dd_eV_A2=float(projected.sum()),
            mean_Phi_dd_eV_A2=float(projected.mean()) if members else None,
            sum_pulling_Jdd_eV_A2=float(-projected.sum()),
            mean_pulling_Jdd_eV_A2=float(-projected.mean()) if members else None,
            full_reference_equivalent_Jdd_sum_eV_A2=float(-projected.mean() * len(reference_members)) if members else None))
    kii = float(phi[source] @ vector)
    sum_j = float(-((hphi @ vector).sum() - kii))
    return dict(source_displacement_direction=vector.tolist(), H_atoms=atoms, shells=shell_rows,
                normalization="Fit mean_pulling_Jdd per existing bond; lattice occupancy is applied once. Full-reference equivalent is diagnostic, not an unfolded shell.",
                onsite_Kii_dd_eV_A2=kii, sum_neighbour_Jdd_eV_A2=sum_j, Kii_minus_sumJ_eV_A2=kii - sum_j,
                H_only_column_sum_vector_eV_A2=hphi.sum(axis=0).tolist(),
                metal_reaction_column_sum_eV_A2=phi[numbers == 46].sum(axis=0).tolist(),
                interpretation="H-only heterogeneous column sum, metals fixed; NOT a global Gamma eigenmode, "
                               "not an all-ion acoustic sum, not an unscreened reciprocal bond law",
                folding="Each occupied periodic H site once. Source-image terms remain folded in onsite. "
                        "Nearest-image offsets label finite-cell shells; half-cell ties use numpy.rint. "
                        "No boundary multiplicity replication, missing neighbours absent.",
                supercell_image_repeat_A=2 * p["a_A"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--mesh", type=int, choices=(3,), required=True)
    parser.add_argument("--probe", choices=("clean", "near"), required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    c = config(args.config)
    folder = ROOT / "loading_Pd" / c["config_id"] / f"k{args.mesh}"
    saved = read_json(folder / "profile.json")
    p = profile(args.config, args.mesh, frozen_execution=saved)
    check_profile(saved, p)
    require(folder == folder_for(p), "Folder/profile differs")
    zero = load_point(folder, p, "zero", None, None, 0.0)
    zero_binding = manifest(folder / "zero.json")
    f0 = np.array(zero["forces_eV_A"])
    selected = cases("probe", args.probe)
    minus, plus = [load_point(folder, p, *case) for case in selected]
    bindings = dict(zero=zero_binding,
                    minus=manifest(folder / f"{selected[0][0]}.json"),
                    plus=manifest(folder / f"{selected[1][0]}.json"))
    response = central_column(p, args.probe, "registered_probe_axis", zero, minus, plus, bindings)
    coupling = response.pop("coupling")
    value = dict(schema=SCHEMA, record_type="frozen_registered_site_loading_readout", profile=p,
                 probe=probe_for(p, args.probe), zero_point=zero_binding,
                 zero_energy_free_eV=zero["energy_free_eV"], zero_forces_eV_A=f0.tolist(),
                 zero_max_all_atom_force_eV_A=float(np.linalg.norm(f0, axis=1).max()),
                 zero_net_force_eV_A=f0.sum(axis=0).tolist(), response=response, coupling=coupling,
                 limitations="No GPW/PAW bytes required for this frozen readout. No spatial inversion reconstruction, "
                             "no quartic/well fit, no relaxation, no full Hessian or global Gamma mode. "
                             "Finite-step SCF, mesh, finite-cell and lattice-source errors are not bounded here.")
    path = args.output or folder / f"readout-{args.probe}-{uuid.uuid4().hex}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    commit_json(path, value)
    import json
    print(json.dumps(manifest(path), allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
