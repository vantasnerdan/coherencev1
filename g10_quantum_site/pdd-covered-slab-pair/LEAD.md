# Pair cut at the covered PdD(111) face (lead's note)

Giuliano's capsule. On the relaxed covered slab of `../pdd-surface-vacancy/`, the adsorbed deuteron and the first deuteron below the surface, which rest 2.44 A apart, are pushed straight toward each other in equal steps with every other atom held. `README.md` is the technical page.

## What it shows

| Separation | Total cost | Subsurface deuteron alone | Adsorbed deuteron alone | Mutual part |
|---:|---:|---:|---:|---:|
| 2.32 A | 0.008 eV | 0.000 eV | 0.007 eV | 0.002 eV |
| 2.20 A | 0.037 | 0.005 | 0.025 | 0.007 |
| 1.96 A | 0.150 | 0.024 | 0.092 | 0.033 |
| 1.72 A | 0.336 | 0.055 | 0.191 | 0.090 |
| 1.48 A | 0.589 | 0.095 | 0.310 | 0.185 |

This pair closes far more cheaply than a pair in the bulk. Bringing it to 1.72 A costs 0.34 eV, against 0.54 eV for the bulk PdD pair at 1.73 A (`../tid2-site/pair/README.md`, interpolated) and that bulk pair starts from 2.89 A. The subsurface deuteron moves toward the surface almost freely: 0.095 eV for 0.48 A.

## Assumptions and limits

- Frozen neighbours, one 2x2 slab, PBE. These are the 3x3x1 mesh values. The first mesh (2x2x1) is kept in the files as the coarse curve; it gave totals 0.004 to 0.039 eV higher.
- One point, 1.72 A, was repeated at 4x4x1: the mutual part moved by -0.002 eV (2.5 percent) and the total by +0.005 eV. The other points have no finer check.
- At the four shortest separations a moving deuteron comes closer to a palladium atom than the sum of the pseudo-atom sphere radii. Those rows are flagged in `pair.json`; the size of the resulting error is not known.

## Lead's check

Archive SHA-256 2734ee7808f3a3ce6f452292fc017dfdd804d4e63adeb1e4d0e3ceb0fb321049, all 49 file sums checked. The delivered reader (NumPy only) runs on the lead's machine, and the lead also recomputed every total, single-deuteron and mutual energy at all three meshes from the per-point records with a separate short script: they match `pair.json`. The first-principles runs themselves were not rerun.
