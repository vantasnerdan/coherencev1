# g3_skew: Ni on Pd, the natural mismatch pattern and the deuterium wells on it

Deposition cannot dial in an arbitrary skew or twist between layers. What it gives is the fixed lattice mismatch: Ni (3.52 Angstrom) on Pd (3.89 Angstrom), about 9.5 percent. This directory models that pattern and makes a first estimate of what it does to the deuterium wells.

    python3 mismatch.py relax --case pd389
    python3 mismatch.py relax --case pd4084
    python3 mismatch.py maps

Needs ASE (for its built-in EMT potential), NumPy and Matplotlib. About two CPU-minutes.

## What is built

A (001) slab of six Pd layers under four Ni layers, bottom two Pd layers fixed, relaxed with EMT.

| Case | Pd spacing | Coincidence | Atoms | In-plane repeat |
|---|---:|---|---:|---:|
| `pd389` | 3.89 Angstrom (bare Pd) | 11 Ni on 10 Pd | 2168 | 38.9 Angstrom |
| `pd4084` | 4.084 Angstrom (PdD spacing) | 7 Ni on 6 Pd | 824 | 24.5 Angstrom |

Loading the Pd swells it, which raises the mismatch to about 16 percent and shortens the pattern. EMT is used for the positions of the metal atoms only; it knows nothing about hydrogen.

## What comes out

- **Pattern.** A smooth square-symmetric corrugation whose strongest modulation runs along the cell diagonals, with wavelength 27.5 Angstrom (bare Pd) and 17.3 Angstrom (PdD spacing). No dislocation network forms in this relaxation.
- **Depth.** The strain modulation is strongest in the Pd layer at the interface and falls to roughly a third of that by the fourth layer. The slab is too thin to give a decay length.
- **Wells (first estimate).** The strain at each site is turned into a collective deuterium stiffness with the strain derivatives from `g2_input` (`adopted-strains.json`). At the bare-Pd spacing every site stays stiff (lowest 0.37 eV/A^2). At the PdD spacing 6 percent of deuterium sites come out with negative stiffness and 15 percent of metal sites below half the bulk value, in the two Pd layers nearest the interface, arranged as patches on the pattern. `stiffness_maps.png` shows them.

## How far to trust the well map

It is an illustration of where soft sites would sit, not a prediction that they exist. The strain derivatives were computed for strains of plus or minus 1 percent; the strains here run from -18 to +8 percent, so every site is an extrapolation. The response to shear is not known and is left out, while shear reaches 10 percent. Hydrogen is absent from the relaxation. The proper strain and shear laws (part A) and one GPAW check on a cut-out (part D) are still to come and will replace this estimate; the `.npz` arrays let the law be swapped without redoing the geometry.

## Files

`pd389/`, `pd4084/`: reference and relaxed structures, per-site strain and stiffness maps for metal sites and for deuterium sites, and the same on the grid the lattice simulation uses. `summary.json`: all numbers. `strain_maps.png`, `stiffness_maps.png`, `penetration.png`.
