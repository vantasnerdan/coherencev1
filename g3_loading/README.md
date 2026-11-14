# g3_loading: how the well and the coupling change with loading

Aim: compute a few loadings once in a 64-atom PdD cell (32 Pd, 32 sites), fit a law, and let the lattice simulation use that law for any loading and arrangement.

## x = 1 (`x1/`)

One D is moved while the other 31 D and 32 Pd stay in place, a = 4.084 A, GPAW PBE, plane waves at 500 eV, 3x3x3 k-points. Details and the exact arrays are in `x1/dense-README.md` and `x1/dense-x1-summary.json`.

| Quantity | Value |
|---|---:|
| On-site stiffness, from energy | 0.599 eV/A^2 |
| On-site stiffness, from force | 0.591 eV/A^2 |
| Sum of neighbour couplings | 0.389 eV/A^2 |
| Collective (all D together) stiffness | 0.202 eV/A^2 |
| Quartic along 100 | 11.6 eV/A^4 |
| Quartic along 111 | 2.23 eV/A^4 |
| Harmonic quantum of one D | 35.0 meV (measured peak: 37 meV) |

The well is `V(u) = K2 u^2/2 + K4 u^4/4`, sampled to 0.24 A. These replace the earlier Gamma-point-only 64-atom values (on-site 1.09, coupling 0.064), which were not converged in k-points.

## x = 0.875, vacancies spread out (`x0875_spread/`)

Four of the 32 sites are empty, as far apart as the cell allows, a = 4.0616 A (the spacing that goes with this loading). Two D atoms are probed: one with no empty nearest-neighbour site, moved along 100, and one with two empty nearest-neighbour sites, moved along the 110 direction that faces a vacancy. Details in `x0875_spread/loading-README.md`.

| D atom | On-site stiffness | Sum of neighbour couplings | Collective | Force at the site |
|---|---:|---:|---:|---:|
| No vacant neighbour | 0.703 eV/A^2 | 0.485 | 0.218 | 0.002 eV/A |
| Two vacant neighbours | 0.529 eV/A^2 | 0.399 | 0.131 | 0.070 eV/A |

A D atom with no vacancy beside it is stiffer than at full loading (0.703 against 0.591), because the lattice is 0.55 percent smaller. A D atom beside vacancies is softer (0.529), and it is not at rest on its lattice site: it feels a steady force of 0.07 eV/A along the probed direction, which a linear estimate turns into an off-centre shift of roughly 0.13 A. The atoms were not relaxed, so that shift is an estimate only.

## x = 0.75, vacancies bunched (`x0750_bunched/`)

Eight of the 32 sites are empty, grouped together, a = 4.0393 A (extrapolated along the same loading line; not a measured point). Same two probes.

| D atom | On-site stiffness | Sum of neighbour couplings | Collective | Force at the site |
|---|---:|---:|---:|---:|
| No vacant neighbour | 0.732 eV/A^2 | 0.508 | 0.223 | 0.002 eV/A |
| Four vacant neighbours | 0.477 eV/A^2 | 0.246 | 0.231 | 0.136 eV/A |

## Across the cases so far

| Loading | a (A) | Vacant neighbours | On-site stiffness | Sum of couplings | Force toward vacancy |
|---|---:|---:|---:|---:|---:|
| 1.0 | 4.084 | 0 | 0.591 | 0.389 | 0 |
| 0.875 | 4.062 | 0 | 0.703 | 0.485 | 0.002 |
| 0.875 | 4.062 | 2 | 0.529 | 0.399 | 0.070 |
| 0.75 | 4.039 | 0 | 0.732 | 0.508 | 0.002 |
| 0.75 | 4.039 | 4 | 0.477 | 0.246 | 0.136 |

Two trends. Lower loading shrinks the lattice and stiffens a D that has all its neighbours. Each vacant neighbour softens the D beside it and adds a push toward the vacancy of roughly 0.035 eV/A per vacant neighbour (0.070 for two, 0.136 for four).

The last case (0.875 bunched, near-vacancy probe only) and the fitted law follow.
