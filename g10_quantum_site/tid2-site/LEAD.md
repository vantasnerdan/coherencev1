# The deuteron in the TiD2 well, and the TiD2 pair beside PdD (lead's note)

Marvin's capsule, built on the first-principles wells in `../tid2-fluorite/`. `README.md` is the technical page. Everything here is for the constrained cubic phase with every site filled and neighbours frozen; it does not say which phase a titanium foil holds.

## What it shows

**A tightly held deuteron.** The first vibration step is 105 meV for D (147 meV for H, 86 meV for T). At 300 K the mean-square displacement per axis is 0.0101 A^2, against 0.082 A^2 in PdD, and the density falls to a millionth of its peak at 0.56 A, against 0.79 A in PdD. The chance of finding the deuteron beyond 0.8 A is 2e-9. The two sides of the body diagonal are kept apart in the tables: the deuteron leans away from the near titanium atom.

**The vibration softens with amplitude.** Along the neighbour direction the first two steps are 106.33 and 104.92 meV for D, a shift of -1.41 meV per quantum. In PdD the shift is +5.34 meV (hardening). The sign decides where a ripple can feed on the vibration: a softening vibration can be unstable where its band curves upward, a hardening one only where it curves downward (near the band top in PdD, see `g13_rogue/`). The TiD2 band itself is not computed, so this gives the sign and size of the non-linearity and no threshold.

**The pair.** At the same separations, bringing two neighbours together costs less in TiD2 than in PdD:

| Separation | Cost in TiD2 | Cost in PdD | Mutual part in TiD2 | Mutual part in PdD |
|---:|---:|---:|---:|---:|
| 1.97 A | 0.086 eV | 0.278 eV | 0.006 eV | 0.065 eV |
| 1.60 A | 0.552 | 0.733 | 0.051 | 0.177 |
| 1.29 A | 1.228 | 1.390 | 0.142 | 0.376 |

(TiD2 values between its computed points are interpolated.)

**Contact.** Passed through the `g11_cmns` contact calculation, with the inner part of the pair interaction borrowed from PdD because none exists for titanium, the density of the pair at contact is 3.5 to 3.7 powers of ten above the PdD bulk pair at 300 K and 7 powers of ten above it for the cold ground state. Dropping the two shortest pair points raises the TiD2 figure by a further 0.06 to 0.5 powers of ten. Most of the difference comes from the pair starting 2.21 A apart instead of 2.89 A. This is a sensitivity to the outer wall under a borrowed interior, not a rate for a titanium foil.

## Assumptions and limits

- The three-dimensional well is filled in between a few computed directions. Changing that fill moves the H tail radius by 0.035 A; the tail digits are not unique.
- The pair points have no mesh check of their own, and the two shortest have no check of the pseudo-atom spheres.
- No neighbour coupling for TiD2 exists yet beyond this one frozen pair, so there is no band, no coupled lattice and no envelope run.

## Lead's check

Archive SHA-256 179927f23f57c09976d67f28b364c90de98e69339e6473927233cea5d92115d0, all 90 file sums checked. The lead reran the solver, its numerical comparisons, the reader exercise, the pair calculation and the report from the inputs alone, on a clean copy of the repository at commit 463910d (NumPy, SciPy, Matplotlib). Every table reproduces: distributions agree to better than 1 part in 1e6 and the contact figures to all printed digits. The stored eigenvectors differ by the free choice of sign and mixing within degenerate levels, which does not affect any table. The three eigenvector files (50 MB each) stay in the release archive and are not in the repository.
