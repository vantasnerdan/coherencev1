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

The arrays and figures here are the originals from the machine that ran the lattice simulation, so the two stay consistent. A rerun on a second machine reproduced the PdD-spacing case exactly and the bare-Pd case to within 1 percent (the relaxation path differs slightly).

## Strain and shear laws from GPAW (`partA/`)

Two-atom PdD cell at 4.084 Angstrom, 500 eV, k = 12, metal held fixed. These are stiffnesses for all deuterium moving together. 200 core-minutes.

**Shear.** Lowest stiffness as the cell is sheared along [100]:

| Shear | 0 | 2% | 4% | 6% | 8% |
|---|---:|---:|---:|---:|---:|
| Lowest stiffness (eV/A^2) | 0.181 | 0.155 | 0.119 | 0.078 | 0.035 |

Shear along [110] at 8 percent gives 0.030. Shear softens steadily and would reach zero near 9 to 10 percent.

**Thin film.** In-plane strain imposed, layer free to relax its height:

| In-plane strain | Height / reference | In-plane stiffness | Across-layer stiffness |
|---:|---:|---:|---:|
| -4% | 1.103 | 0.868 | -0.341 |
| -2% | 1.066 | 0.395 | -0.290 |
| 0 | 1.032 | 0.084 | -0.181 |
| +2% | 1.000 | -0.127 | 0.057 |
| +4% | 0.971 | -0.281 | 0.421 |

A compressed film softens the motion across the layer; a stretched film softens the motion within it. Between them, at every strain in this range, one direction has a negative bare stiffness. These are small-vibration values at the site centre: zero-point motion and the anharmonic part of the well add stiffness, so a negative entry marks a direction where the well is flat or double-bottomed, not a demonstrated instability. The film with no imposed strain still grows 3.2 percent taller, because this method's own preferred spacing is larger than 4.084 Angstrom.

`strain_tensor_law.py` with `strain-law.json` returns the stiffness tensor for a given strain inside the computed range and refuses outside it unless asked to extrapolate. None of the sites in the EMT mismatch maps above fall inside that range, so the well map there remains an extrapolation. `strain-family.png` shows both families.


## The spacing the method itself prefers

All PdD numbers here use a lattice constant of 4.084 Angstrom (an extrapolation of measurements at loadings 0.80 to 0.98 to full loading; Schirber and Morosin, Phys. Rev. B 12, 117 (1975)). With the cell free to relax, the same GPAW setup settles at 4.134 Angstrom, 1.2 percent larger. There the collective small-vibration stiffness is -0.055 eV/A^2, against +0.181 at 4.084. So the PdD inputs describe a crystal this method regards as slightly compressed, and at the method's own spacing the bottom of the well is flat to slightly humped: the deuterium is held by the anharmonic walls and its zero-point motion. Strain thresholds quoted from 4.084 Angstrom would be about 1.2 percentage points smaller if counted from 4.134.

This point was not used in fitting the strain law above. The law predicts -0.058 for it, an error of 0.003 eV/A^2, which is a fair independent check of the law for uniform expansion.
