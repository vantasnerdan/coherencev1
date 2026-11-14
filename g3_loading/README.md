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

The reduced loadings (0.875 spread, 0.875 bunched, 0.75) and the fitted law follow.
