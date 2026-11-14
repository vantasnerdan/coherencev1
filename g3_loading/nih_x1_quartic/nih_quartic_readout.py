"""Read actual +0.24 A points; fixed-K2, explicitly one-point directional quartics."""
import argparse
import json
from pathlib import Path

from dense_pd import commit_json, manifest, read_json, require
from nih_profile import FOLDER as BASE_FOLDER, load_point as load_base
from nih_quartic_profile import DIRECTIONS, FOLDER, load_point, point_path, profile


def one_point(k2, u, energy, force):
    """Same energy + u-scaled-force least-squares convention as analyse.fit_direction."""
    energy_column = u**4 / 4
    force_column = u**4
    energy_target = energy - k2*u*u/2
    force_target = -u*force - k2*u*u
    joint = (energy_column*energy_target + force_column*force_target) / (energy_column**2 + force_column**2)
    energy_only = energy_target / energy_column
    force_only = force_target / force_column
    return dict(K4_joint_eV_A4=float(joint), K4_energy_eV_A4=float(energy_only),
        K4_force_eV_A4=float(force_only), coefficient_of_u4_eV_A4=float(joint/4),
        energy_residual_eV=float(k2*u*u/2 + joint*u**4/4-energy),
        force_residual_eV_A=float(-k2*u-joint*u**3-force),
        relative_energy_force_K4_difference=float(abs(energy_only-force_only)/max(abs(energy_only),abs(force_only)))
        if energy_only or force_only else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--direction", choices=DIRECTIONS, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(len(args.direction) == len(set(args.direction)), "Duplicate quartic direction")
    saved = read_json(FOLDER / "profile.json")
    p = profile(frozen=saved)
    require(saved == p, "Frozen NiH quartic profile differs")
    zero = load_base(p["base_profile"], "zero", 0.)
    results = {}
    for direction in args.direction:
        r = load_point(p, direction)
        delta = r["energy_free_eV"]-zero["energy_free_eV"]
        fit = one_point(p["K2_fixed_eV_A2"], p["amplitude_A"], delta, r["force_parallel_eV_A"])
        results[direction] = dict(point_source=manifest(point_path(direction)), delta_energy_eV=delta,
            measured_force_parallel_eV_A=r["force_parallel_eV_A"], amplitude_A=p["amplitude_A"], **fit)
    value = dict(schema="nih-one-point-quartic-readout-v1", profile=p, directions=results,
        K2_fixed_eV_A2=p["K2_fixed_eV_A2"], zero_source=manifest(BASE_FOLDER / "zero.json"),
        model="V(u)-V(0)=K2*u^2/2+K4*u^4/4; joint residuals in energy and u-scaled force",
        limitations="One positive amplitude per reported direction, not a multi-amplitude Taylor convergence "
            "or spectrum. Negative amplitudes use ideal cubic inversion symmetry; fixed source-x K2 uses "
            "cubic isotropy. No new harmonic/negative point/relaxation/higher-order or mesh study.")
    commit_json(args.output, value)
    print(json.dumps(manifest(args.output), allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
