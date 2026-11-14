"""Plot actual ±.03 probes, explicitly removing the measured origin force tilt.

python plot_loading_case.py CLEAN-readout.json NEAR-readout.json OUTPUT.png
Guides use the central finite-step coefficients, not a measured minimum/well fit.
"""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent


def main():
    clean, near, output = map(Path, sys.argv[1:])
    records = [json.loads(p.read_text()) for p in (clean, near)]
    assert records[0]["profile"]["config_id"] == records[1]["profile"]["config_id"]
    fig, axes = plt.subplots(2, 2, figsize=(9, 6))
    for row, r in enumerate(records):
        response = r["response"]
        d = np.array(response["source_displacement_vector"])
        i = r["probe"]["moving_index"]
        points = [json.loads((ROOT / response["measured_points"][name]["path"]).read_text())
                  for name in ("minus", "zero", "plus")]
        u = np.array([p["displacement_A"] for p in points])
        force = np.array([np.array(p["forces_eV_A"])[i] @ d for p in points])
        energy = np.array([p["energy_free_eV"] for p in points])
        tilt = response["F0_parallel_eV_A"]
        x = np.linspace(-.03, .03, 150)
        axes[row, 0].scatter(u, 1000 * (energy - energy[1] + tilt * u), color="tab:orange", zorder=3)
        axes[row, 0].plot(x, 1000 * .5 * response["K2_energy_central_eV_A2"] * x**2, color="tab:blue")
        axes[row, 1].scatter(u, force - tilt, color="tab:orange", zorder=3)
        axes[row, 1].plot(x, -response["K2_eV_A2"] * x, color="tab:blue")
        label = f"{r['probe']['probe_id']}: {r['probe']['n_vac_NN']} vacant NN, {r['probe']['direction_label']}"
        axes[row, 0].set_title(label, fontsize=10)
        axes[row, 1].set_title(f"Kii={response['K2_eV_A2']:.5f} eV/A²; F0={tilt:.5f} eV/A", fontsize=10)
        axes[row, 0].set_ylabel("E(u)-E(0)+F0 u (meV)")
        axes[row, 1].set_ylabel("Fparallel(u)-F0 (eV/A)")
        for ax in axes[row]:
            ax.set_xlabel("Registered displacement u (A)")
            ax.axhline(0, color="0.7", linewidth=.5)
            ax.grid(alpha=.15)
    p = records[0]["profile"]
    fig.suptitle(f"PdD x={p['loading_x']}, {p['arrangement']}, reference a={p['a_A']:.6f} A\n"
                 "Origin force tilt REMOVED for display; registered sites are not relaxed minima", fontsize=11)
    fig.text(.5, .015, "Dots: actual SCFs. Lines: central harmonic guides within sampled window only. No loading quartic measured.",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .035, 1, .91))
    fig.savefig(output, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
