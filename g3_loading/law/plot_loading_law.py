"""Actual loading tangents and fitted in-domain curves; no invented error bars."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from loading_law import features


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("law", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    law = json.loads(args.law.read_text())
    rows = [law["anchor"], *law["training_rows"]]
    names = ["Kii_dd", *sorted(n for n in law["models"] if n != "Kii_dd")]
    with_motif = len(law["descriptor_names"]) == 3
    keys = sorted({(r["n_vac_NN"], r["vacant_NN_pair_count"] if with_motif else None) for r in rows})
    colors = dict(zip(keys, plt.get_cmap("tab10").colors))
    groups = list(dict.fromkeys(r["config_id"] for r in rows))
    markers = dict(zip(groups, ("*", "o", "^", "s")))
    spacings = np.linspace(min(r["a_reference_A"] for r in rows), 4.084, 240)
    fig, axes = plt.subplots(3, 2, figsize=(10, 10), sharex=True)
    for ax, name in zip(axes.flat, names):
        model = law["models"][name]
        for row in rows:
            value = row["observables"][name]
            if value is not None:
                key = (row["n_vac_NN"], row["vacant_NN_pair_count"] if with_motif else None)
                ax.scatter(row["a_reference_A"], value, color=colors[key],
                           marker=markers[row["config_id"]], s=65 if row["config_id"] == "x1" else 35,
                           zorder=3, edgecolors="black", linewidths=.4)
        if model["identified"]:
            domain = model["domain"]
            equations = np.array(domain["normalized_halfspaces"])
            coefficient = np.array([model["coefficients_eV_A2"][key] for key in law["descriptor_names"]])
            for n, pairs in keys:
                x = np.array([features(n, a, vacant_NN_pair_count=pairs) for a in spacings])
                inside = np.all((x / domain["scale"]) @ equations[:, :-1].T + equations[:, -1] <= 1e-10, axis=1)
                values = model["anchor_eV_A2"] + x @ coefficient
                ax.plot(spacings, np.where(inside, values, np.nan), color=colors[(n, pairs)], linewidth=1.3)
            ax.text(.03, .96, f"fit RMS {model['rms_residual_eV_A2']:.4g} eV/A²",
                    transform=ax.transAxes, va="top", fontsize=9)
        else:
            ax.text(.03, .96, "slopes not identifiable", transform=ax.transAxes, va="top", fontsize=9)
        title = "Onsite Kii along probe axis" if name == "Kii_dd" else "Mean Jdd per existing bond: r/a=" + name.removeprefix("Jbar_shell_").rstrip("0").rstrip(".")
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("eV/A²")
        ax.axhline(0, color="0.7", linewidth=.6)
        ax.grid(alpha=.15)
    for ax in axes[-1]:
        ax.set_xlabel("Loading reference a (A), not extra local strain")
    legend = [Line2D([], [], color=colors[(n, pairs)],
                    label=f"{n} vacant NN" + (f", {pairs} adjacent pairs" if with_motif else ""))
              for n, pairs in keys]
    legend += [Line2D([], [], color="0.3", marker=markers[g], linestyle="none", label=g) for g in groups]
    fig.legend(handles=legend, loc="lower center", ncol=4, fontsize=8, bbox_to_anchor=(.5, .018))
    fig.suptitle("PdD registered-site tangents: " + law["model_id"] + "\n"
                 "Clean100 / vacancy-facing110; curves only inside each observable's sampled hull", fontsize=12)
    fig.tight_layout(rect=(0, .085, 1, .94))
    fig.savefig(args.output, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
