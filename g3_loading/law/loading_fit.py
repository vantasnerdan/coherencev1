"""Fit the preregistered anchored law from complete actual frozen readouts.

--dense X1.json --readout CASE-clean.json --readout CASE-near.json ... --output UNIQUE.json
Require each declared probe set; use --near-only-config CASE for a commissioned near-only case.
"""
import argparse
import copy
from pathlib import Path
from itertools import combinations

import numpy as np

from dense_pd import ROOT, commit_json, manifest, read_json, require
from loading_lattice import a as lattice_a
from loading_law import A_REFERENCE, fit
from loading_readout import hydrogen_column

PHYSICAL_FIELDS = ("kmesh", "gamma_centered", "cutoff_eV", "xc", "setups", "spinpol",
                   "smearing_eV", "nbands", "mixer", "maxiter", "convergence",
                   "eigensolver", "symmetry_settings", "software")


def observables(coupling):
    values = {"Kii_dd": coupling["onsite_Kii_dd_eV_A2"]}
    inventory = {}
    total = 0.0
    for shell in coupling["shells"]:
        key = f"Jbar_shell_{shell['radius_over_a']:.10f}"
        count, mean = shell["occupied_multiplicity"], shell["mean_pulling_Jdd_eV_A2"]
        summed = shell["sum_pulling_Jdd_eV_A2"]
        require((count == 0 and mean is None and summed == 0)
                or (count > 0 and mean is not None and abs(count * mean - summed) < 1e-10),
                "Per-existing-bond normalization differs from occupied shell sum")
        require(key not in values, "Duplicate shell")
        values[key] = mean
        inventory[key] = dict(radius_over_a=shell["radius_over_a"],
            occupied_multiplicity=count, reference_multiplicity=shell["original_reference_multiplicity"],
            sum_Jdd_eV_A2=summed, mean_Jdd_eV_A2=mean,
            occupied_original_indices=shell["occupied_original_indices"],
            vacant_original_indices=shell["vacant_original_indices"])
        total += summed
    require(abs(total - coupling["sum_neighbour_Jdd_eV_A2"]) < 1e-10, "Shell/neighbour sums differ")
    return values, inventory


def vacant_nn_motif(probe):
    vacant = [item for item in probe["NN_images"] if not item["occupied"]]
    require(sorted(item["original_atom_index"] for item in vacant) == probe["vacant_NN_original_indices"]
            and len(vacant) == probe["n_vac_NN"], "Vacant nearest-neighbour descriptor inventory differs")
    pairs = []
    for left, right in combinations(vacant, 2):
        delta = np.array(left["offset_over_a"]) - right["offset_over_a"]
        delta -= 2 * np.rint(delta / 2)  # The registered periodic cell is exactly 2a.
        if abs(delta @ delta - .5) < 1e-10:
            pairs.append(sorted([left["original_atom_index"], right["original_atom_index"]]))
    return len(pairs), sorted(pairs)


def dense_anchor(path):
    dense = read_json(path)
    p = dense["profile"]
    require(p["a_A"] == A_REFERENCE and p["loading_x"] == 1
            and dense["coupling"]["independent_negative_scf"], "Not the measured x1 anchor")
    # Reuse the one shell normalizer on the ACTUAL dense all-ion force column.
    # The canonical fully occupied source has the identity original-index map.
    adapted = copy.deepcopy(p)
    adapted["geometry"]["original_atom_indices"] = list(range(len(p["geometry"]["atomic_numbers"])))
    adapted["vacancy_original_indices"] = []
    phi = np.array(dense["coupling"]["all_ion_Phi_vectors_eV_A2"])
    require(phi.shape == (64, 3) and np.isfinite(phi).all(), "Dense column shape differs")
    coupling = hydrogen_column(adapted, p["geometry"]["moving_index"], phi, np.array([1., 0., 0.]))
    require(abs(coupling["onsite_Kii_dd_eV_A2"] - dense["coupling"]["on_site_k2_eV_A2"]) < 1e-12
            and abs(coupling["sum_neighbour_Jdd_eV_A2"]
                    - dense["coupling"]["sum_neighbour_pulling_J_xx_eV_A2"]) < 1e-12,
            "Dense-to-loading normalizer changed adopted values")
    values, inventory = observables(coupling)
    return dict(id="x1:clean", config_id="x1", probe_id="clean", loading_x=1,
                a_reference_A=A_REFERENCE, n_vac_NN=0, vacant_NN_pair_count=0,
                vacant_NN_pairs_original_indices=[], direction_vector=[1., 0., 0.],
                observables=values, shell_inventory=inventory, source=manifest(path)), p


def loading_row(path, dense_profile):
    r = read_json(path)
    p, probe, response = r["profile"], r["probe"], r["response"]
    require(r["record_type"] == "frozen_registered_site_loading_readout"
            and response["independent_negative_scf"], "Not an actual frozen loading pair")
    require(all(p[key] == dense_profile[key] for key in PHYSICAL_FIELDS)
            and p["charge"] == 0, "Loading/anchor electronic profiles differ")
    require(p["a_A"] == lattice_a(p["loading_x"])["a_angstrom"], "Loading reference spacing differs")
    require(np.allclose(response["source_displacement_vector"], probe["direction_vector"], atol=1e-12, rtol=0),
            "Probe axis differs")
    values, inventory = observables(r["coupling"])
    pair_count, pairs = vacant_nn_motif(probe)
    return dict(id=f"{p['config_id']}:{probe['probe_id']}", config_id=p["config_id"],
                probe_id=probe["probe_id"], loading_x=p["loading_x"], a_reference_A=p["a_A"],
                arrangement=p["arrangement"], n_vac_NN=probe["n_vac_NN"],
                vacant_NN_pair_count=pair_count, vacant_NN_pairs_original_indices=pairs,
                direction_vector=probe["direction_vector"],
                F0_source_eV_A=response["F0_source_eV_A"],
                F0_parallel_eV_A=response["F0_parallel_eV_A"],
                Kii_energy_eV_A2=response["K2_energy_central_eV_A2"],
                raw_all_ion_column_sum_eV_A2=response["all_ion_column_sum_eV_A2"],
                observables=values, shell_inventory=inventory, source=manifest(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dense", type=Path, required=True)
    parser.add_argument("--readout", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--motif", action="store_true",
                        help="Add only the count of adjacent pairs among vacant nearest neighbours")
    parser.add_argument("--near-only-config", action="append", default=[],
                        help="Explicitly commissioned config with near pair only, not missing clean data")
    args = parser.parse_args()
    anchor, profile = dense_anchor(args.dense)
    rows = [loading_row(path, profile) for path in args.readout]
    require(len({r["id"] for r in rows}) == len(rows), "Duplicate probe rows")
    groups = {r["config_id"] for r in rows}
    near_only = set(args.near_only_config)
    require(near_only <= groups, "Near-only declaration names an absent config")
    sampling = {group: ["near"] if group in near_only else ["clean", "near"] for group in groups}
    for group, expected in sampling.items():
        require({r["probe_id"] for r in rows if r["config_id"] == group} == set(expected),
                "Readout set differs from the explicitly declared sampling plan")
    law = fit(rows, anchor, with_motif=args.motif)
    law["commissioned_probe_sets"] = dict(sorted(sampling.items()))
    law["vacancy_motif_definition"] = ("Number of unordered pairs of vacant sites in the source's "
        "12-site nearest-neighbour shell that are themselves nearest neighbours, r/a=1/sqrt(2), "
        "with minimum images in the registered 2a periodic cell. One extra descriptor only.")
    law["source_files"] = [manifest(ROOT / name) for name in
        ("loading_fit.py", "loading_law.py", "loading_readout.py", "loading_lattice.py")]
    law["direction_closure"] = ("Dense100 is the fixed x1 anchor; loading clean100 and near110-family "
        "are measured under that policy. A scalar/isotropic transfer is a closure, not measured tensors.")
    commit_json(args.output, law)
    import json
    print(json.dumps(manifest(args.output), allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
