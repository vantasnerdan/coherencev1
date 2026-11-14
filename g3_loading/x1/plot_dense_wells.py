"""Plot actual local energy/force samples and the frozen quartic fit; no DFT."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    if data["probe_only"] or set(data["directions"]) != {"100", "111"}:
        raise ValueError("Figure requires both complete actual local-well directions")
    if data["profile"]["operator"] != "LOCAL_one_D_other_nuclei_fixed":
        raise ValueError("Figure input is not the local frozen-host operator")
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 5.4), sharex=True)
    zero_force = np.array(data["zero_point"]["force_source_eV_A"])
    for col, direction in enumerate(("100", "111")):
        fit = data["directions"][direction]
        points = fit["points"]
        if [p["displacement_A"] for p in points] != [0.12, 0.24]:
            raise ValueError("Figure requires actual +0.12/+0.24 A points")
        k2, k4 = fit["k2_eV_A2"], fit["k4_eV_A4"]
        u = np.linspace(0, 0.24, 160)
        measured_u = np.array([0.] + [p["displacement_A"] for p in points])
        measured_e = np.array([0.] + [p["delta_energy_free_eV"] for p in points])
        d = np.array([1., 0., 0.]) if direction == "100" else np.ones(3) / np.sqrt(3)
        measured_f = np.array([float(zero_force @ d)] + [p["force_parallel_eV_A"] for p in points])
        energy, force = axes[0, col], axes[1, col]
        energy.plot(u, 1000 * (k2 * u**2 / 2 + k4 * u**4 / 4), label="joint quartic fit")
        energy.scatter(measured_u, 1000 * measured_e, marker="o", color="C1", label="actual SCF")
        force.plot(u, k2 * u + k4 * u**3)
        force.scatter(measured_u, -measured_f, marker="o", color="C1")
        energy.set_title(f"[{direction}]  K4 = {k4:.4g} eV/A⁴")
        force.set_xlabel("D displacement u (A)")
        for ax in (energy, force):
            ax.grid(alpha=0.25)
            ax.set_xlim(-0.005, 0.245)
            ax.tick_params(labelsize=9)
    axes[0, 0].set_ylabel("V(u) - V(0) (meV)")
    axes[1, 0].set_ylabel("Restoring force -F (eV/A)")
    axes[0, 0].legend(fontsize=8, loc="upper left")
    fig.suptitle(f"PdD local x=1, a={data['profile']['a_A']:.3f} A; K2={k2:.5f} eV/A²", fontsize=12)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=130)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
