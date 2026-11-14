"""Produce registered-site PdD vacancy configs and an exhaustive geometry inventory.

Usage: python loading_geometry.py [--run-id UNIQUE_NAME]
Outputs: loading-configs/x0875-{spread4,bunched4}.json and x0750-bunched8.json;
uniquely named inventory/full-findings JSON under evidence/loading-geometry/.
Publication is atomic and no-clobber. Reuse requires reading the existing configs,
not rerunning this producer over them. No relaxation, GPAW, MPI or fit is run.
The vacancy objective is a graph edge count, not an energy or softness prediction.
"""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import combinations
from math import comb, sqrt
from pathlib import Path
import re
import resource
import sys
import time

import numpy as np
from ase.build import bulk

from dense_pd import ROOT, commit_json, digest, manifest, read_json, require
from loading_lattice import a

BASE_SHA256 = "466191aede61708287ae5433147748a4bf113c3c50afc6251656a5fa277c9a7c"
CLEAN = 25
SPECS = (("x0875-spread4", 0.875, "spread4", 4, 2, "min"),
         ("x0875-bunched4", 0.875, "bunched4", 4, 2, "max"),
         ("x0750-bunched8", 0.75, "bunched8", 8, 4, "max"))
SOURCE_NAMES = ("loading_geometry.py", "loading_lattice.py",
                "loading-lattice-source.json", "Pd/base.json", "dense_pd.py")
LIMITATIONS = [
    "Registered ideal sites with all remaining nuclei fixed; not relaxed minima.",
    "Vacancy graph extrema are not energy, softness or stability extrema.",
    "Lattice reference is beta-PdD at77 K; .875 is model interpolation and .75 is authorized extrapolation; source pointwise uncertainty is unknown.",
    "No electronic, Native, force, stiffness or quartic calculations; actual consumer verification belongs to Main.",
    "Exact finite enumeration covers this19-site admissible set and the fixed occupied near probe, not arbitrary supercells or vacancy counts.",
]


def canonical_sites():
    """Fence original bytes and ASE ordering before deriving any vacancy geometry."""
    require(manifest(ROOT / "Pd/base.json")["sha256"] == BASE_SHA256,
            "Canonical Pd/base.json bytes differ")
    base = read_json(ROOT / "Pd/base.json")
    canonical = bulk("PdH", "rocksalt", a=base["a_A"], cubic=True).repeat((2, 2, 2))
    require(base["atoms"] == 64 and base["moving_index"] == CLEAN,
            "Canonical count/source index differs")
    require(base["atomic_numbers"] == canonical.numbers.tolist(), "Species/order differ")
    position_error = float(np.max(np.abs(np.asarray(base["reference_positions_A"]) - canonical.positions)))
    cell_error = float(np.max(np.abs(np.asarray(base["cell_A"]) - canonical.cell)))
    require(position_error < 1e-10 and cell_error < 1e-10,
            "Source is not ASE cubic rocksalt2x2x2")
    half_a = np.asarray(base["reference_positions_A"]) / (base["a_A"] / 2)
    integer_sites = np.rint(half_a).astype(int)
    require(np.max(np.abs(half_a - integer_sites)) < 1e-10, "Nonregistered source positions")
    hydrogen = [i for i, z in enumerate(base["atomic_numbers"]) if z == 1]
    require(len(hydrogen) == 32 and base["atomic_numbers"].count(46) == 32,
            "Expected Pd32H32 source")
    return base, hydrogen, integer_sites, dict(
        canonical_bytes_sha256=BASE_SHA256, ASE_builder="bulk('PdH','rocksalt',a=4.084,cubic=True).repeat((2,2,2))",
        species_and_order_match=True, maximum_position_error_A=position_error,
        maximum_cell_error_A=cell_error, integer_half_a_coordinates=integer_sites.tolist())


def fcc_neighbors(hydrogen, sites):
    """12 unique FCC octahedral neighbours using all three periodic axes."""
    graph = {}
    for i in hydrogen:
        adjacent = []
        for j in hydrogen:
            delta = (sites[j] - sites[i] + 2) % 4 - 2
            if int(delta @ delta) == 2:
                adjacent.append(j)
        require(len(adjacent) == 12 and i not in adjacent, "Not a12-NN FCC graph")
        graph[i] = adjacent
    require(all(i in graph[j] for i in hydrogen for j in graph[i]), "Asymmetric NN graph")
    return graph


def choose_near(hydrogen, allowed, graph):
    """Fix one shared occupied near site before optimizing any objective."""
    candidates = []
    for near in hydrogen:
        if near == CLEAN:
            continue
        eligible = [i for i in allowed if i != near]
        adjacent = len(set(eligible).intersection(graph[near]))
        counts = {}
        for k, target in ((4, 2), (8, 4)):
            remote = len(eligible) - adjacent
            counts[str(k)] = (comb(adjacent, target) * comb(remote, k - target)
                              if adjacent >= target and remote >= k - target else 0)
        candidates.append(dict(original_atom_index=near, available_vacancy_sites=len(eligible),
                               available_NN_sites=adjacent, target_feasible_counts=counts))
    feasible = [row for row in candidates if all(row["target_feasible_counts"].values())]
    return (feasible[0]["original_atom_index"] if feasible else None), candidates


def exhaustive_inventory(allowed, graph, near, k):
    """Enumerate all19-choose-k; separately expose all and occupied-near constraints."""
    adjacency = [sum(1 << j for j, original in enumerate(allowed) if original in graph[i])
                 for i in allowed]
    near_mask = sum(1 << j for j, original in enumerate(allowed) if original in graph[near])
    near_bit = (1 << allowed.index(near)) if near in allowed else 0
    all_edges, by_count, extrema = Counter(), defaultdict(Counter), {}
    occupied_count = 0
    for subset in combinations(range(len(allowed)), k):
        mask = sum(1 << j for j in subset)
        edges = sum((adjacency[j] & mask).bit_count() for j in subset) // 2
        all_edges[edges] += 1
        if mask & near_bit:
            continue
        occupied_count += 1
        count = (near_mask & mask).bit_count()
        by_count[count][edges] += 1
        vacancies = tuple(allowed[j] for j in subset)
        if count not in extrema:
            extrema[count] = {mode: dict(edge_count=edges, tie_count=0,
                                       vacancy_original_indices=list(vacancies)) for mode in ("min", "max")}
        for mode in ("min", "max"):
            row = extrema[count][mode]
            better = edges < row["edge_count"] if mode == "min" else edges > row["edge_count"]
            if better:
                row.update(edge_count=edges, tie_count=1, vacancy_original_indices=list(vacancies))
            elif edges == row["edge_count"]:
                row["tie_count"] += 1
    total = sum(all_edges.values())
    require(total == comb(len(allowed), k), "Enumeration count differs from combinatorial total")
    return dict(vacancy_count=k, total_19_choose_k=total,
                all_subsets_edge_count_histogram=dict(sorted(all_edges.items())),
                occupied_near_feasible_count=occupied_count,
                occupied_near_by_n_vac_NN={str(n): dict(feasible_count=sum(hist.values()),
                    edge_count_histogram=dict(sorted(hist.items())), extrema=extrema[n])
                    for n, hist in sorted(by_count.items())})


def probe_metadata(probe_id, original, retained, vacancies, graph, sites):
    vacant = sorted(set(graph[original]).intersection(vacancies))
    images = []
    for neighbor in graph[original]:
        raw = sites[neighbor] - sites[original]
        shift = -np.floor_divide(raw + 2, 4)
        images.append(dict(original_atom_index=neighbor, image_shift=shift.tolist(),
                           offset_over_a=((raw + 4 * shift) / 2).tolist(),
                           occupied=neighbor not in vacancies))
    facing = next((item for item in images if item["original_atom_index"] == min(vacant)), None) if vacant else None
    if probe_id == "near":
        require(facing is not None, "Near probe has no vacant nearest neighbour")
        ray = np.asarray(facing["offset_over_a"], dtype=float)
        direction = (ray / np.linalg.norm(ray)).tolist()
        label = "vacancy_facing_110_family"
    else:
        direction, label = [1.0, 0.0, 0.0], "100"
    return dict(probe_id=probe_id, moving_index=retained.index(original),
                original_atom_index=original, n_vac_NN=len(vacant),
                NN_original_indices=graph[original], vacant_NN_original_indices=vacant,
                NN_images=images, direction_vector=direction, direction_label=label,
                vacancy_facing_original_index=facing["original_atom_index"] if probe_id == "near" else None,
                vacancy_facing_image_shift=facing["image_shift"] if probe_id == "near" else None,
                vacancy_facing_offset_over_a=facing["offset_over_a"] if probe_id == "near" else None)


def vacancy_metrics(vacancies, graph, sites, lattice_a):
    pairs, distance_counts = [], Counter()
    for i, j in combinations(vacancies, 2):
        delta = (sites[j] - sites[i] + 2) % 4 - 2
        squared = int(delta @ delta)
        distance_counts[squared] += 1
        if j in graph[i]:
            pairs.append([i, j])
    return dict(NN_edge_count=len(pairs), NN_edges_original_indices=pairs,
                all_pair_count=comb(len(vacancies), 2),
                minimum_image_pair_distance_histogram=[dict(squared_half_a_distance=squared,
                    r_over_a=sqrt(squared) / 2, r_A=lattice_a * sqrt(squared) / 2, count=count)
                    for squared, count in sorted(distance_counts.items())],
                edge_convention="Undirected original-site pairs at r=a/sqrt(2), full3D PBC; no farther shells or image multiplicity")


def make_config(spec, optimum, base, sites, graph, near, sources, allowed, feasible_count):
    config_id, x, arrangement, k, target, mode = spec
    vacancies = optimum["vacancy_original_indices"]
    retained = [i for i in range(64) if i not in vacancies]
    lattice = a(x)
    factor = lattice["a_angstrom"] / base["a_A"]
    geometry = dict(atomic_numbers=[base["atomic_numbers"][i] for i in retained],
                    reference_positions_A=(np.asarray(base["reference_positions_A"])[retained] * factor).tolist(),
                    cell_A=(np.asarray(base["cell_A"]) * factor).tolist(), pbc=[True] * 3,
                    original_atom_indices=retained)
    probes = [probe_metadata(name, original, retained, vacancies, graph, sites)
              for name, original in (("clean", CLEAN), ("near", near))]
    require(probes[0]["n_vac_NN"] == 0 and probes[1]["n_vac_NN"] == target,
            "Selected geometry does not meet probe constraints")
    require(len(vacancies) == k and geometry["atomic_numbers"].count(1) / 32 == x,
            "Selected loading differs")
    metrics = vacancy_metrics(vacancies, graph, sites, lattice["a_angstrom"])
    require(metrics["NN_edge_count"] == optimum["edge_count"], "Objective/config mismatch")
    return dict(schema="loading-pd-config-v2", config_id=config_id, loading_x=x,
                arrangement=arrangement, a_A=lattice["a_angstrom"], lattice_source=lattice,
                source_files=sources, vacancy_original_indices=vacancies, geometry=geometry,
                probes=probes, vacancy_pair_metrics=metrics,
                original_to_retained_atom_indices=[retained.index(i) if i in retained else None for i in range(64)],
                vacancy_design=dict(allowed_original_indices=allowed, fixed_near_original_index=near,
                    target_near_vac_NN=target, objective="vacancy_vacancy_NN_edge_count", extremum=mode,
                    feasible_count=feasible_count, extremal_tie_count=optimum["tie_count"],
                    tie_break="lexicographically smallest sorted original vacancy tuple"),
                reference_status="scaled_registered_sites; unrelaxed; no energy_or_softness_optimization")


def main():
    start_wall, start_cpu = time.monotonic(), time.process_time()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    args = parser.parse_args()
    require(re.fullmatch(r"[A-Za-z0-9_.-]+", args.run_id) is not None, "Unsafe run id")
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))
    config_dir, evidence_dir = ROOT / "loading-configs", ROOT / "evidence/loading-geometry"
    inventory_path = evidence_dir / f"inventory-{args.run_id}.json"
    findings_path = evidence_dir / f"full-findings-{args.run_id}.json"
    targets = [config_dir / f"{spec[0]}.json" for spec in SPECS] + [inventory_path, findings_path]
    require(not any(path.exists() for path in targets), "Immutable publication target already exists")
    sources = [manifest(ROOT / name) for name in SOURCE_NAMES]
    base, hydrogen, sites, canonical_validation = canonical_sites()
    graph = fcc_neighbors(hydrogen, sites)
    allowed = [i for i in hydrogen if i != CLEAN and i not in graph[CLEAN]]
    require(len(allowed) == 19, "Expected19 allowed vacancy sites")
    near, candidates = choose_near(hydrogen, allowed, graph)
    inventories = ({str(k): exhaustive_inventory(allowed, graph, near, k) for k in (4, 8)}
                   if near is not None else {})
    selected, configs = [], []
    for spec in SPECS:
        config_id, x, arrangement, k, target, mode = spec
        constrained = inventories.get(str(k), {}).get("occupied_near_by_n_vac_NN", {}).get(str(target))
        if constrained is None:
            selected.append(dict(config_id=config_id, disposition="infeasible_under_exact_constraints",
                                 target_near_vac_NN=target, config_written=False))
            continue
        optimum = constrained["extrema"][mode]
        configs.append(make_config(spec, optimum, base, sites, graph, near, sources,
                                   allowed, constrained["feasible_count"]))
        selected.append(dict(config_id=config_id, disposition="exact_target_feasible", config_written=True,
                             feasible_count=constrained["feasible_count"], target_near_vac_NN=target,
                             objective_extremum=mode, **optimum))
    require(sources == [manifest(ROOT / name) for name in SOURCE_NAMES], "Sources changed during inventory")
    config_dir.mkdir(parents=True, exist_ok=True)
    evidence_dir.mkdir(parents=True, exist_ok=True)
    for config in configs:
        commit_json(config_dir / f"{config['config_id']}.json", config)
    inventory = dict(schema="loading-geometry-inventory-v1", run_id=args.run_id, source_files=sources,
        canonical_validation=canonical_validation, clean_original_H=CLEAN, fixed_near_original_H=near,
        near_selection_rule="lowest feasible original H index distinct from clean; fixed before objective optimization",
        near_candidate_exact_combinatorial_counts=candidates, allowed_vacancy_original_indices=allowed,
        original_H_indices=hydrogen, full_PBC_FCC_NN_original_indices={str(i): graph[i] for i in hydrogen},
        graph_undirected_edge_count=sum(map(len, graph.values())) // 2,
        objective="undirected vacancy-vacancy FCC NN edge count", enumeration=inventories,
        selections=selected, limitations=LIMITATIONS,
        cost=dict(cpu_s=time.process_time() - start_cpu, wall_s=time.monotonic() - start_wall))
    commit_json(inventory_path, inventory)
    outputs = [manifest(config_dir / f"{config['config_id']}.json") for config in configs]
    outputs.append(manifest(inventory_path))
    findings = dict(schema="loading-geometry-full-findings-v1", run_id=args.run_id,
        disposition="complete" if len(configs) == 3 else "actual_constraint_infeasibility",
        source_files=sources, output_files=outputs,
        actual_command=[sys.executable, "loading_geometry.py", "--run-id", args.run_id],
        findings=selected, inventory=manifest(inventory_path),
        source_profile_id=digest(sources), common_clean_original_index=CLEAN, common_near_original_index=near,
        retained_atom_counts={config["config_id"]: len(config["geometry"]["atomic_numbers"]) for config in configs},
        lattice_values={config["config_id"]: config["lattice_source"] for config in configs},
        checks_run=["One actual bounded exhaustive geometry inventory; canonical ASE geometry comparison and graph/config invariants inside producer"],
        checks_not_run=["No synthetic tests/suites, build, lint, formatter, GPAW, MPI, Native or electronic calculation"],
        limitations=LIMITATIONS, cost=dict(cpu_s=time.process_time() - start_cpu,
                                         wall_s=time.monotonic() - start_wall),
        documentation="Module docstring documents CLI and JSON implementation plan records boundaries and exact objective/constraints.")
    commit_json(findings_path, findings)
    print(__import__("json").dumps(dict(disposition=findings["disposition"],
          full_findings=manifest(findings_path), output_files=outputs, selections=selected), indent=2))


if __name__ == "__main__":
    main()
