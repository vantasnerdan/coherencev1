"""NiH x=1 k3 control identity; original geometry, no new quartic or relaxation."""
from importlib.metadata import version
from pathlib import Path

import numpy as np
from ase.build import bulk

from dense_pd import ROOT, digest, manifest, read_json, require

SCHEMA = "nih-local-k3-v1"
BASE_SHA = "d6b76bdbcb4434aab337326469c70e1aca5a3eb7f63173f9442470830b7d70fd"
FOLDER = ROOT / "dense_Ni" / "k3"
CODE = ("nih_profile.py", "nih_dense.py", "nih_readout.py", "run-nih.sh",
        "dense_pd.py", "loading_pd.py", "loading_lattice.py", "loading-lattice-source.json",
        "analyse.py", "Ni/base.json", "Ni/points.jsonl")
CASES = {"zero": [("zero", 0.)], "probe": [("probe_minus", -.03), ("probe_plus", .03)]}


def profile(*, frozen=None):
    base = read_json(ROOT / "Ni/base.json")
    require(manifest(ROOT / "Ni/base.json")["sha256"] == BASE_SHA, "Original NiH geometry source differs")
    atoms = bulk("NiH", "rocksalt", a=base["a_A"], cubic=True).repeat((2, 2, 2))
    require(base["metal"] == "Ni" and base["a_A"] == 3.732 and base["moving_index"] == 25
            and base["atomic_numbers"] == atoms.numbers.tolist()
            and np.allclose(base["reference_positions_A"], atoms.positions, atol=1e-10, rtol=0)
            and np.allclose(base["cell_A"], atoms.cell, atol=1e-10, rtol=0), "Not the original Ni32H32 registered source")
    sources = [manifest(ROOT / name) for name in CODE]
    if frozen is None:
        sources += [manifest(ROOT / "setups" / name) for name in ("Ni.10.PBE.gz", "H.PBE.gz")]
        software = {name: version(name) for name in ("gpaw", "ase", "numpy", "scipy")}
    else:
        require(frozen["sources"][:len(sources)] == sources and len(frozen["sources"]) == len(sources) + 2,
                "Frozen NiH source provenance differs")
        sources, software = frozen["sources"], frozen["software"]
    value = dict(schema=SCHEMA, material="NiH", loading_x=1, a_A=base["a_A"],
        operator="LOCAL_one_H_other_nuclei_fixed", electronic_species="H", vibrational_isotope="H",
        kmesh=[3, 3, 3], cutoff_eV=500, xc="PBE", setups={"Ni": "10"},
        gamma_centered=True, charge=0, spinpol=False, smearing_eV=.1, nbands=-20,
        mixer=[.05, 5, 50], maxiter=80, eigensolver={"name": "davidson", "niter": 2},
        symmetry_settings={"point_group": False},
        convergence={"energy": 1e-7, "density": 1e-4, "forces": 1e-4},
        expected_electrons=352, geometry=dict(atomic_numbers=atoms.numbers.tolist(),
            reference_positions_A=atoms.positions.tolist(), cell_A=atoms.cell.array.tolist(),
            moving_index=25, H_indices=np.flatnonzero(atoms.numbers == 1).tolist(), pbc=[True]*3),
        sources=sources, software=software, scope="Numerical-profile correction, not an isolated k-point study. No new NiH quartic.")
    value["profile_id"] = digest(value)
    return value


def positions(p, u):
    value = np.array(p["geometry"]["reference_positions_A"])
    value[p["geometry"]["moving_index"], 0] += u
    return value


def load_point(p, name, u, *, verify_checkpoint=False):
    path = FOLDER / f"{name}.json"
    r = read_json(path)
    require(r["schema"] == SCHEMA and r["record_type"] == "measured_scf" and r["profile"] == p,
            "NiH point/profile identity differs")
    require(r["point_id"] == name and r["displacement_A"] == u
            and np.allclose(r["positions_A"], positions(p, u), atol=1e-10, rtol=0), "NiH displacement differs")
    force = np.asarray(r["forces_eV_A"])
    require(force.shape == (64, 3) and np.isfinite(force).all() and np.isfinite(r["energy_free_eV"])
            and abs(r["force_parallel_eV_A"] - force[25, 0]) < 1e-12, "Incomplete NiH SCF values")
    layout = r["parallel"]
    require(r["electron_count"] == 352 and layout["world_size"] == layout["kpt"] > 0
            and layout["domain"] == layout["band"] == 1 and r["point_group_operations"] == 1,
            "NiH electronic/layout provenance differs")
    affinity = r["rank_affinity"]
    require(len(affinity) == layout["world_size"]
            and all(a["rank"] == i and len(a["actual_cpus"]) == 1
                    and set(a["actual_cpus"]) <= set(a["inherited_cpus"]) for i, a in enumerate(affinity))
            and len({a["actual_cpus"][0] for a in affinity}) == len(affinity), "NiH affinity differs")
    if u:
        require(r["restart_point"] == manifest(FOLDER / "zero.json")
                and r["initialization"] == "full_self_consistent_zero_checkpoint_restart", "Not common-zero NiH restart")
    else:
        require(r["restart_point"] is None and r["initialization"] == "fresh_zero_atomic_density", "NiH zero was not fresh")
    if verify_checkpoint:
        require(r["checkpoint"] == manifest(path.with_suffix(".gpw")), "NiH checkpoint bytes differ")
    return r
