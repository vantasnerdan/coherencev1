"""Plot the actual Born-stiffness samples, not an extrapolated spatial map."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    table = json.loads(Path("strain-table.json").read_text())
    rows = table["rows"]
    fig, axes = plt.subplots(1,2,figsize=(8.8,3.5),layout="constrained")
    shear = sorted([r for r in rows if r["family"] == "shear100"],key=lambda r:r["spec"]["engineering_gamma"])
    gamma = [0]+[100*r["spec"]["engineering_gamma"] for r in shear]
    lowest = [table["K0_eV_A2"]]+[np.linalg.eigvalsh(r["K_material_eV_A2"])[0] for r in shear]
    axes[0].plot(gamma,lowest,"o-",label="shear along [100]")
    for row in rows:
        if row["family"] == "shear110":
            axes[0].scatter([100*row["spec"]["engineering_gamma"]],[np.linalg.eigvalsh(row["K_material_eV_A2"])[0]],marker="s",label="[110] sample")
    axes[0].set(xlabel="engineering shear (%)",ylabel="lowest Gamma stiffness (eV/A²)",title="Shear splits the collective modes")
    film = sorted([r for r in rows if r["family"] == "biaxial_relaxed"],key=lambda r:r["spec"]["inplane_strain"])
    eps = [100*r["spec"]["inplane_strain"] for r in film]
    inplane = [np.trace(np.asarray(r["K_material_eV_A2"])[:2,:2])/2 for r in film]
    normal = [r["K_material_eV_A2"][2][2] for r in film]
    axes[1].plot(eps,inplane,"o-",label="in-plane (x/y mean)")
    axes[1].plot(eps,normal,"s-",label="across layer (z)")
    axes[1].set(xlabel="biaxial in-plane strain (%)",ylabel="Gamma stiffness (eV/A²)",title="Out-of-plane length relaxed")
    for ax in axes:
        ax.axhline(0,color="0.5",lw=.7)
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
    fig.suptitle("PdD, a reference 4.084 Å; PW500 / Gamma12",fontsize=11)
    fig.supxlabel("Fixed-metal collective Born response at registered D sites; not a local well or quantum-stability test",fontsize=8)
    fig.savefig("strain-family.png",dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
