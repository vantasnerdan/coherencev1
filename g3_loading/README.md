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

The well is `V(u) = K2 u^2/2 + K4 u^4/4`, sampled to 0.24 A. These replace the earlier Gamma-point-only 64-atom values (on-site 1.09, coupling 0.064), which used coarser settings throughout (400 eV cutoff, one k-point, looser force convergence).

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

## Across the cases

| Loading | a (A) | Vacant neighbours | On-site stiffness | Sum of couplings | Force toward vacancy |
|---|---:|---:|---:|---:|---:|
| 1.0 | 4.084 | 0 | 0.591 | 0.389 | 0 |
| 0.875 | 4.062 | 0 | 0.703 | 0.485 | 0.002 |
| 0.875 | 4.062 | 2 | 0.529 | 0.399 | 0.070 |
| 0.75 | 4.039 | 0 | 0.732 | 0.508 | 0.002 |
| 0.75 | 4.039 | 4 | 0.477 | 0.246 | 0.136 |
| 0.875 bunched | 4.062 | 2 | 0.431 | 0.128 | 0.109 |

(The unlabelled 0.875 rows are the spread arrangement; the 0.75 rows are bunched.)

Three results. Lower loading shrinks the lattice and stiffens a D that has all its neighbours. A vacant neighbour softens the D beside it and pushes it toward the vacancy. And the arrangement matters, not only the count: at the same loading, lattice spacing and number of vacant neighbours (two), a D in the bunched cell is softer (0.431 against 0.529), far more weakly coupled to its neighbours (0.128 against 0.399) and pushed harder (0.109 against 0.070 eV/A) than in the spread cell. A law in only (lattice spacing, number of vacant neighbours) cannot reproduce both; the fit keeps that misfit visible.

## x = 0.875, vacancies bunched (`x0875_bunched/`)

Four vacancies grouped together, a = 4.0616 A. Only the D beside the vacancies was probed (three calculations), to save budget.

## The fitted law (`law/`)

`law/loading_law.py` with `law/loading-law-base.json`: each quantity is its x = 1 value plus a term linear in the lattice spacing and a term linear in the number of vacant nearest neighbours,

    y = y(x=1) + c_eps (a/4.084 - 1) + c_n n_vac

For the on-site stiffness, c_eps = -14.4 eV/A^2 and c_n = -0.077 eV/A^2 per vacant neighbour. RMS residual 0.045 eV/A^2, largest 0.085, largest error predicting a left-out cell 0.108. The neighbour couplings are fitted shell by shell in the same form.

A second version with one extra term (number of adjacent vacancy pairs among the nearest neighbours, `loading-law-motif.json`) fits the computed points slightly better (RMS 0.043) but predicts a left-out cell much worse (error up to 0.40), so the two-term law is the one to use. It carries the trends with loading; it does not distinguish bunched from spread arrangements, and the three computed cells are the reference for that.

Cost of the PdD loading set, including the x = 1 well: about 2,510 core-minutes.

## NiH x = 1 check (`nih_x1/`)

The NiH single-hydrogen well in `g2_input` was computed with a 400 eV cutoff and one k-point. Repeating the zero and +/-0.03 A probe with the settings adopted above (500 eV, 3x3x3 k-points, tight forces), a = 3.732 A:

| eV/A^2 | g2_input | Here | Change |
|---|---:|---:|---:|
| On-site stiffness | 1.468 | 1.712 | +17% |
| Sum of neighbour couplings | 0.234 | 0.676 | +189% |
| Collective (on-site minus couplings) | 1.234 | 1.036 | -16% |

The collective value now agrees with the two-atom-cell result (1.050) to 1.3 percent, which the old values did not. These replace the `g2_input` NiH local inputs. One further point at the same settings, +0.24 A along 100 (`nih_x1_quartic/`), gives a 100 quartic of 16.2 eV/A^4 (from one displaced point with the harmonic term fixed; `g2_input` had 18.3). The 111 quartic was not recomputed and remains the `g2_input` value, 1.13.

Total cost of this directory: about 2,870 core-minutes.
