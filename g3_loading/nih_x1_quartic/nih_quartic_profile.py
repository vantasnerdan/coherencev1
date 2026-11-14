"""One-point quartic extension of the unchanged, completed NiH k3 source epoch."""

import numpy as np

from dense_pd import ROOT, digest, manifest, read_json, require
from nih_profile import FOLDER as BASE_FOLDER, profile as base_profile

FOLDER = BASE_FOLDER / "quartic"
SCHEMA = "nih-one-point-quartic-v1"
CODE = ("nih_quartic_profile.py", "nih_quartic.py", "nih_quartic_readout.py", "run-nih-quartic.sh")
HARMONIC = BASE_FOLDER / "readout-20261114T0917Z.json"
HARMONIC_SHA = "ca59c4098b718b67647d493cf98e5cde872d0ba1bdbafab5e6c0b496456e00c3"
DIRECTIONS = {"100": [1., 0., 0.], "111": (np.ones(3) / np.sqrt(3)).tolist()}


def profile(*, frozen=None):
    p = base_profile(frozen=frozen["base_profile"] if frozen else None)
    require(read_json(BASE_FOLDER / "profile.json") == p, "Base NiH profile differs")
    harmonic = manifest(HARMONIC)
    require(harmonic["sha256"] == HARMONIC_SHA, "Fixed NiH harmonic source differs")
    r = read_json(HARMONIC)
    require(r["profile"] == p, "Harmonic NiH profile differs")
    value = dict(schema=SCHEMA, base_profile=p, base_profile_source=manifest(BASE_FOLDER / "profile.json"),
        harmonic_readout=harmonic, K2_fixed_eV_A2=r["coupling"]["on_site_k2_eV_A2"],
        amplitude_A=.24, directions=DIRECTIONS, extension_sources=[manifest(ROOT / name) for name in CODE],
        scope="One actual positive amplitude per direction; fixed source-x K2 transferred by cubic symmetry. "
              "No additional harmonic, negative-amplitude, higher-order or convergence calculation.")
    value["profile_id"] = digest(value)
    return value


def point_path(direction):
    require(direction in DIRECTIONS, "Unknown NiH quartic direction")
    return FOLDER / f"q{direction}_u024.json"


def positions(p, direction):
    value = np.array(p["base_profile"]["geometry"]["reference_positions_A"])
    value[25] += p["amplitude_A"] * np.array(p["directions"][direction])
    return value


def load_point(p, direction, *, verify_checkpoint=False):
    path = point_path(direction)
    r = read_json(path)
    require(r["schema"] == SCHEMA and r["record_type"] == "measured_scf" and r["profile"] == p,
            "NiH quartic point/profile differs")
    require(r["direction"] == direction and r["displacement_A"] == .24
            and np.allclose(r["direction_vector"], p["directions"][direction], atol=1e-12, rtol=0)
            and np.allclose(r["positions_A"], positions(p, direction), atol=1e-10, rtol=0),
            "NiH quartic signed direction/positions differ")
    force = np.asarray(r["forces_eV_A"])
    require(force.shape == (64, 3) and np.isfinite(force).all() and np.isfinite(r["energy_free_eV"])
            and abs(r["force_parallel_eV_A"] - force[25] @ p["directions"][direction]) < 1e-12,
            "Incomplete NiH quartic values")
    layout = r["parallel"]
    require(r["electron_count"] == 352 and layout["world_size"] == layout["kpt"] > 0
            and layout["domain"] == layout["band"] == 1 and r["point_group_operations"] == 1,
            "NiH quartic electronic/layout identity differs")
    affinity = r["rank_affinity"]
    require(len(affinity) == layout["world_size"]
            and all(a["rank"] == i and len(a["actual_cpus"]) == 1
                    and set(a["actual_cpus"]) <= set(a["inherited_cpus"]) for i, a in enumerate(affinity))
            and len({a["actual_cpus"][0] for a in affinity}) == len(affinity), "NiH quartic affinity differs")
    require(r["restart_point"] == manifest(BASE_FOLDER / "zero.json")
            and r["initialization"] == "full_self_consistent_zero_checkpoint_restart",
            "NiH quartic does not use the common zero")
    if verify_checkpoint:
        require(r["checkpoint"] == manifest(path.with_suffix(".gpw")), "NiH quartic checkpoint differs")
    return r
