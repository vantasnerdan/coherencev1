# TiD2 wells and pair cut (lead's note)

Giuliano's first-principles wells for a deuteron in cubic fluorite TiD2, the second hydride after PdD to get signed one-deuteron cuts and a pair cut. `README.md` is the technical page. The cell is the constrained cubic phase with every tetrahedral site filled; it is a fingerprint of that phase and does not say which titanium deuteride phase a foil holds.

## What it shows

**The lattice.** Lattice parameter 4.423 A. Nearest deuterons rest 2.21 A apart, against 2.89 A in PdD.

**One deuteron displaced, all else fixed** (energy rise in eV):

| Displacement | Along the cube axis (toward the neighbouring deuteron) | Along the face diagonal | Along the body diagonal, toward the titanium atom | Along the body diagonal, away from it |
|---:|---:|---:|---:|---:|
| 0.12 A | 0.040 | 0.041 | 0.047 | 0.034 |
| 0.24 A | 0.157 | 0.165 | 0.219 | 0.117 |
| 0.48 A | 0.582 | 0.703 | 1.203 | 0.334 |
| 0.80 A | 1.361 | 2.131 | 5.294 | 0.585 |

The well is stiff and strongly lopsided along the body diagonal: at 0.8 A it costs nine times more to move toward the titanium atom than away from it, into the empty octahedral site. The stiffness at the bottom (0.040 eV at 0.12 A) corresponds to a harmonic vibration quantum of about 107 meV for a deuteron (lead's arithmetic), against a first transition of 42.6 meV in PdD.

**Two neighbouring deuterons pushed toward each other along the cube axis:**

| Separation | Total cost | Of which mutual (beyond the two single-site costs) |
|---:|---:|---:|
| 1.97 A | 0.086 eV | 0.006 eV |
| 1.73 A | 0.342 | 0.028 |
| 1.25 A | 1.321 | 0.157 |
| 1.01 A | 2.017 | 0.297 |

For comparison, the PdD pair at 1.29 A costs 1.39 eV, of which 0.38 eV is mutual (`../pdd-nearest-pair/`). In TiD2 the pair starts closer, and its mutual repulsion at about the same separation is less than half as large; the cost is mostly each deuteron climbing its own well.

## Assumptions and limits

Frozen neighbours and cell, PBE, one 24-atom cell. The mesh check (4 to 6 points per axis) was done at one geometry, the 0.48 A point toward titanium: it moved the rise by 0.002 eV of 1.203 eV. The pair points have no mesh check of their own, and the two shortest have no separate check of the pseudo-atom spheres at close approach.

## Lead's check

Archive SHA-256 53734575158a4e9645775f380fc600e3eb25ff6080cefd3e10a30d229b3a9aea, all 53 file sums checked. The delivered replay script needs ASE, which is not on the lead's machine. The lead recomputed every cut energy and all four mutual pair energies directly from the per-point records with a separate short script: they match the delivered tables exactly. The first-principles runs themselves were not rerun.
