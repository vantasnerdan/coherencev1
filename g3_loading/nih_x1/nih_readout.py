"""Frozen NiH k3 readout and source-bound comparison to the original coarse data."""
import argparse
import json
from pathlib import Path

import numpy as np

from analyse import coupling
from dense_pd import ROOT, commit_json, manifest, read_json, require
from nih_profile import FOLDER, load_point, profile


KEYS = ("on_site_k2_eV_A2", "sum_neighbour_pulling_J_xx_eV_A2", "H_only_q0_stiffness_eV_A2")


def coupling_rows(p, minus, plus):
    indices = p["geometry"]["H_indices"]
    return [dict(kind="coupling", direction="100", displacement_A=r["displacement_A"],
                 H_indices=indices, forces_H_eV_A=np.array(r["forces_eV_A"])[indices].tolist(),
                 net_force_eV_A=np.array(r["forces_eV_A"]).sum(axis=0).tolist()) for r in (minus, plus)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    saved = read_json(FOLDER / "profile.json")
    p = profile(frozen=saved)
    require(saved == p, "Frozen NiH profile differs")
    zero = load_point(p, "zero", 0.)
    minus = load_point(p, "probe_minus", -.03)
    plus = load_point(p, "probe_plus", .03)
    current = coupling(p["geometry"], coupling_rows(p, minus, plus))
    original_base = read_json(ROOT / "Ni/base.json")
    original_rows = [json.loads(line) for line in (ROOT / "Ni/points.jsonl").read_text().splitlines()]
    old = coupling(original_base, original_rows)
    changes = {key: (current[key] - old[key]) / abs(old[key]) for key in KEYS}
    energy_k = (plus["energy_free_eV"] + minus["energy_free_eV"] - 2 * zero["energy_free_eV"]) / .03**2
    force_k = current[KEYS[0]]
    phi = -(np.array(plus["forces_eV_A"]) - minus["forces_eV_A"]) / .06
    raw_acoustic = np.array(current["all_ion_acoustic_residual_eV_A2"])
    value = dict(schema="nih-k3-control-readout-v1", record_type="actual_signed_NiH_control",
        profile=p, point_sources={name: manifest(FOLDER / f"{name}.json")
            for name in ("zero", "probe_minus", "probe_plus")},
        coupling=current, Kii_energy_eV_A2=float(energy_k),
        relative_energy_force_difference=float(abs(energy_k-force_k)/max(abs(energy_k), abs(force_k))),
        raw_acoustic_over_abs_Kii=float(np.linalg.norm(raw_acoustic)/abs(force_k)),
        all_ion_Phi_vectors_eV_A2=phi.tolist(),
        zero_source_force_eV_A=zero["forces_eV_A"][25],
        original_coarse_coupling=old, relative_changes=changes,
        update_decision="update" if any(abs(v) >= .1 for v in changes.values()) else "retain",
        decision_rule="Axiom's approximately10percent retain/update control; report every actual relative change",
        original_sources=[manifest(ROOT / name) for name in ("Ni/base.json", "Ni/points.jsonl", "analyse.py")],
        limitations="Actual one-H source-x central column, finite periodic image-folded H interactions. "
            "H-onlyq0 is the same-matrix fullyoccupied cubic H-subspace value with metals fixed, "
            "not a coupled all-ion acoustic mode. PW400Gamma/density1e-6 to PW500k3/density1e-4/force1e-4 "
            "is a profile correction, not an isolated mesh derivative. No NiH quartic, relaxation, "
            "principal-axis tensor or new mesh-convergence certificate.")
    commit_json(args.output, value)
    print(json.dumps(manifest(args.output), allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
