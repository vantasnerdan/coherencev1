"""Frozen-output consumer for dense_pd.py; performs no electronic calculation.

--mesh 3|4 --a 4.084 [--probe-only] [--output UNIQUE.json]. Probe-only reads the
actual central force pair; the default also requires zero and both complete
local well directions. Point physical profiles and source-byte provenance must
match; GPW bytes are needed for Native restarts, not frozen force/energy readout.
The consumer never imports or substitutes primitive collective stiffness records.
"""
import argparse
import json
from pathlib import Path
import uuid

import numpy as np
from scipy.constants import e, hbar, physical_constants

from dense_pd import (ROOT, SCHEMA, CASES, check_profile, commit_json, digest,
                      inversion_permutation, load_point, manifest, profile,
                      read_json, require)


def central_column(p, minus, plus, minus_manifest, plus_manifest):
    geometry = p["geometry"]
    permutation = np.array(inversion_permutation(geometry))
    measured = np.array(plus["forces_eV_A"])
    measured_minus = np.array(minus["forces_eV_A"])
    h = plus["displacement_A"]
    require(plus["point_id"] == "probe" and minus["point_id"] == "probe_minus"
            and plus["direction"] == minus["direction"] == "100"
            and h == 0.03 and minus["displacement_A"] == -h,
            "Central column requires actual matching -0.03/+0.03 A x probes")
    phi = -(measured - measured_minus) / (2 * h)
    indices = np.array(geometry["H_indices"])
    source = geometry["moving_index"]
    source_h = list(indices).index(source)
    hphi = phi[indices]
    positions = np.array(geometry["reference_positions_A"])
    cell = np.array(geometry["cell_A"])
    unwrapped = positions[indices] - positions[source]
    fractional = unwrapped @ np.linalg.inv(cell)
    image_shifts = -np.rint(fractional).astype(int)
    offsets = unwrapped + image_shifts @ cell
    distances = np.linalg.norm(offsets, axis=1)
    atom_rows, shells = [], {}
    for i, atom in enumerate(indices):
        atom_rows.append(dict(atom_index=int(atom), inversion_partner=int(permutation[atom]),
                              is_source=int(atom) == source, R_reference_unwrapped_A=unwrapped[i].tolist(),
                              minimum_image_translation=image_shifts[i].tolist(), R_A=offsets[i].tolist(),
                              distance_A=float(distances[i]), Phi_vector_eV_A2=hphi[i].tolist(),
                              Phi_ix_jx_eV_A2=float(hphi[i, 0]), Phi_iy_jx_eV_A2=float(hphi[i, 1]),
                              Phi_iz_jx_eV_A2=float(hphi[i, 2]), pulling_J_xx_eV_A2=float(-hphi[i, 0])))
        if i != source_h:
            shells.setdefault(round(float(distances[i]), 6), []).append(i)
    shell_rows = []
    for radius, members in sorted(shells.items()):
        values = hphi[members]
        shell_rows.append(dict(distance_A=radius, distinct_supercell_H_sites=len(members),
                               atom_indices=indices[members].tolist(),
                               sum_Phi_vector_eV_A2=values.sum(axis=0).tolist(),
                               mean_Phi_xx_eV_A2=float(values[:, 0].mean()),
                               sum_Phi_xx_eV_A2=float(values[:, 0].sum()),
                               sum_pulling_J_xx_eV_A2=float(-values[:, 0].sum()),
                               min_Phi_xx_eV_A2=float(values[:, 0].min()),
                               max_Phi_xx_eV_A2=float(values[:, 0].max())))
    kii = float(hphi[source_h, 0])
    sum_j = float(-(hphi[:, 0].sum() - hphi[source_h, 0]))
    inversion_error = measured_minus + measured[permutation]
    return dict(record_type="measured_central_force_column",
                measured_positive_point=plus_manifest, measured_negative_point=minus_manifest,
                independent_negative_scf=True,
                definition="Phi_i(alpha),j(x)=-(F_i(alpha)(+h)-F_i(alpha)(-h))/(2h); h=0.03 A",
                column_scope="Actual all-atom responses to source x; not a full Hessian, "
                             "source-y/z columns or extracted reciprocal bond kernel.",
                inversion_diagnostic=dict(used_for_reconstruction=False,
                    max_H_force_error_eV_A=float(np.linalg.norm(inversion_error[indices], axis=1).max()),
                    max_H_column_component_difference_eV_A2=float(np.abs(inversion_error[indices]).max() / (2 * h))),
                on_site_k2_eV_A2=kii, on_site_definition="-(F_jx(+h)-F_jx(-h))/(2h); finite-step periodic-site local curvature",
                sum_neighbour_pulling_J_xx_eV_A2=sum_j,
                H_only_q0_stiffness_eV_A2=float(hphi[:, 0].sum()),
                q0_definition="Kii-sum_neighbour(Jxx); H-only subspace, metals fixed, not all-ion acoustic zero",
                H_column_sum_vector_eV_A2=hphi.sum(axis=0).tolist(),
                metal_reaction_column_sum_eV_A2=phi[np.array(geometry["atomic_numbers"]) == 46].sum(axis=0).tolist(),
                all_ion_acoustic_residual_eV_A2=phi.sum(axis=0).tolist(),
                all_ion_Phi_vectors_eV_A2=phi.tolist(),
                image_folding="Each distinct periodic H site counted once; response sums all images of the periodic source. "
                              "Offsets label nearest-image distance shells, not unfolded infinite-lattice shells. "
                              "Half-cell ties use numpy.rint (ties to even). Source diagonal also contains folded "
                              "source-image terms; boundary images are not multiplied or divided.",
                supercell_image_repeat_A=2 * p["a_A"], H_atoms=atom_rows, shells=shell_rows)


def fit_direction(rows, zero, k2):
    u = np.array([row["displacement_A"] for row in rows])
    energy = np.array([row["energy_free_eV"] - zero["energy_free_eV"] for row in rows])
    force = np.array([row["force_parallel_eV_A"] for row in rows])
    require(sorted(u.tolist()) == [0.12, 0.24], "Incomplete or duplicated local well direction")
    scale = float(max(abs(u)))
    design_energy = np.c_[u * u / 2, u ** 4 / 4]
    design_force = np.c_[u, u ** 3]
    matrix = np.vstack([design_energy, scale * design_force])
    target = np.r_[energy, -scale * force]
    column = matrix[:, 1]
    k4 = float(np.dot(column, target - matrix[:, 0] * k2) / np.dot(column, column))
    efit = design_energy @ [k2, k4]
    ffit = -design_force @ [k2, k4]
    mass = physical_constants["deuteron mass"][0]
    quantum = float(hbar * np.sqrt(abs(k2) * e / 1e-20 / mass) / e)
    return dict(k2_eV_A2=k2, k4_eV_A4=k4, coefficient_of_u2_eV_A2=k2 / 2,
                coefficient_of_u4_eV_A4=k4 / 4, energy_fit_RMS_eV=float(np.sqrt(np.mean((efit - energy) ** 2))),
                force_fit_RMS_eV_A=float(np.sqrt(np.mean((ffit - force) ** 2))),
                energy_fit_max_abs_eV=float(max(abs(efit - energy))),
                force_fit_max_abs_eV_A=float(max(abs(ffit - force))),
                least_squares="Fixed probe K2, joint energy and force; force rows scaled by max(u)=0.24 A, as analyse.py",
                harmonic_quantum_eV=quantum if k2 > 0 else None,
                imaginary_harmonic_magnitude_eV=quantum if k2 < 0 else None,
                vibrational_mass_kg=mass, sampled_positive_range_A=[0.12, 0.24],
                k2_assumption="Cubic x=1 reference gives isotropic quadratic coefficient; finite +100 probe supplies both directions",
                inversion_assumption="Even local potential about the ideal centrosymmetric site; no negative SCF measurements added",
                residual_convention="fitted minus measured; forces projected along the displacement unit vector",
                points=[dict(point_id=row["point_id"], displacement_A=float(u[i]),
                             delta_energy_free_eV=float(energy[i]), force_parallel_eV_A=float(force[i]),
                             fitted_energy_eV=float(efit[i]), fitted_force_eV_A=float(ffit[i]),
                             energy_residual_eV=float(efit[i] - energy[i]),
                             force_residual_eV_A=float(ffit[i] - force[i])) for i, row in enumerate(rows)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=int, choices=(3, 4), required=True)
    parser.add_argument("--a", type=float, required=True)
    parser.add_argument("--probe-only", action="store_true", help="Read actual local central column without requiring well stages")
    parser.add_argument("--output", type=Path, help="New immutable JSON path; default unique readout in the mesh folder")
    args = parser.parse_args()
    folder = ROOT / "dense_Pd" / f"k{args.mesh}"
    frozen = read_json(folder / "profile.json")
    p = profile(args.mesh, args.a, frozen_execution=frozen)
    check_profile(frozen, p)
    stages = ("probe",) if args.probe_only else tuple(CASES)
    rows, sources = {}, {}
    for stage in stages:
        for name, direction, u in CASES[stage]:
            rows[name] = load_point(folder, p, name, direction, u)
            sources[name] = manifest(folder / f"{name}.json")
    column = central_column(p, rows["probe_minus"], rows["probe"], sources["probe_minus"], sources["probe"])
    output = dict(schema=SCHEMA, record_type="derived_local_readout", profile=p,
                  profile_id=p["profile_id"], probe_only=args.probe_only,
                  measured_point_sources=sources, consumer_source=manifest(Path(__file__)),
                  convention="V(u)-V(0)=K2*u^2/2+K4*u^4/4; one D displaced, other D and Pd fixed",
                  operator_separation="LOCAL periodic-site column and local wells only; no primitive collective data consumed",
                  coupling=column, directions={})
    if not args.probe_only:
        for direction in ("100", "111"):
            selected = [rows[name] for name, _, _ in CASES[f"well{direction}"]]
            output["directions"][direction] = fit_direction(selected, rows["zero"], column["on_site_k2_eV_A2"])
        output["zero_point"] = dict(measured_source=sources["zero"], energy_free_eV=rows["zero"]["energy_free_eV"],
                                    force_source_eV_A=rows["zero"]["force_source_eV_A"])
    output["content_id"] = digest(output)
    path = args.output or folder / f"readout-{uuid.uuid4().hex}.json"
    commit_json(path, output)
    print(json.dumps(manifest(path), allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
