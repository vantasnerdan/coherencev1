"""Electronic calculation profile shared by the primitive response and well grids."""
from pathlib import Path
from gpaw import GPAW, PW, FermiDirac, Mixer


def calculator(metal, kmesh, folder, cutoff_eV=500):
    return GPAW(mode=PW(cutoff_eV), setups={metal: "10"}, xc="PBE",
                kpts={"size": (kmesh,)*3, "gamma": True}, nbands=-8,
                occupations=FermiDirac(0.1), spinpol=False,
                symmetry={"point_group": True}, mixer=Mixer(0.05, 5, 50), maxiter=80,
                convergence={"energy": 1e-7, "density": 1e-6}, txt=str(Path(folder)/"gpaw.out"))
