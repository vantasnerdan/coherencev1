# g0_stability: does the published ordering survive without the rotating-wave approximation?

Gamberale and Modanese (arXiv:2205.14015) derive a coherent ordered state of charged oscillators in a lattice, using the rotating-wave approximation. Their Hamiltonian is quadratic, so the question can be settled exactly by finding the normal modes of the collective ("bright") matter and field modes. Run:

    python3 stability.py > stability.json

NumPy only.

## Two versions of the same Hamiltonian

**As printed.** The paper removes the diamagnetic `A^2` term mode by mode, which turns the field frequency into `omega' = sqrt(omega^2 + omega_p^2)`. It then sums the `p.A` coupling over all propagation directions, which multiplies the coupling by `sqrt(s)` with `s = 8 pi/3`. With every counter-rotating term kept, the bright block is

    H = omega' (x^2 + p^2 + q^2 + r^2)/2 + sqrt(s) omega_p p q

and its lowest mode has `nu^2 = omega'^2 (1 - sqrt(s) eps)`, with `eps = omega_p/omega'`.

**Same field in both terms.** If the field that enters `p.A` (uniform over the sites, summed over directions) is also the field that enters `A^2`, the direction sum multiplies both. The bright block is then

    H = (P + g Q)^2/2 + omega'^2 X^2/2 + Pi^2/2 + omega^2 Q^2/2,   g = sqrt(s) omega_p

which is a sum of squares and cannot have an unstable mode.

## Result

Lowest normal-mode frequency squared, in units of `omega'^2`:

| `eps` | As printed | Same field in both terms |
|---:|---:|---:|
| 0.25 | 0.276 | 0.471 |
| 0.3455 | 0 | 0.348 |
| 0.52 | -0.505 | 0.192 |
| 0.56 | -0.621 | 0.165 |
| 0.69 | -0.997 | 0.097 |
| 0.99 | -1.865 | 0.002 |

- The printed reduction is unstable above `eps = sqrt(3/(8 pi)) = 0.3455`. This is exact, not an artefact of the rotating-wave approximation.
- The consistent version is stable for every `eps` below 1.
- The instability therefore comes from applying the direction-sum factor to the coupling and not to the diamagnetic term. The starting Hamiltonian, `(p + eA)^2/2m` plus a non-negative potential, is non-negative, so no consistent reduction of it can lower the energy below the uncoupled ground state.

`stability.json` also holds the authors' fixed-axis comparison (`s = 4`).

## Limits

This is a statement about the quadratic model in the dipole approximation, as the paper sets it up. It does not cover a finite domain treated beyond the dipole approximation, anharmonic wells, or any mechanism outside minimal coupling. The paper's abstract refers to a gauge argument about what counts as radiation; that argument was not evaluated here. Both matrices were also diagonalised directly as a cross-check and agree with the closed forms.
