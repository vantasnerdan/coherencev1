# Lead's note on the bulk refit through 0.8 A

The PdD and NiH wells were refitted with the computed points out to 0.8 A from `bulk-wide-wells/`.

- **PdD.** The deuteron's mean-square displacement is 0.0816 A^2 at 300 K and 0.0638 A^2 at 77 K, which is 1.3 and 0.5 percent above the original quartic well. With points to 0.8 A, the choice of how to continue the well beyond the last point no longer changes these numbers (it moved them by 2 percent when the points stopped at 0.48 A). The chance of finding the deuteron more than 0.8 A from its site is about 7e-6 at 300 K and 2e-7 at 77 K.
- **The old well was too soft toward the palladium atom.** At 0.8 A along [100] the computed energy is 2.03 eV; the original quartic gave 1.38 eV, and the range allowed by the 0.48 A refit (0.88 to 1.50 eV) did not contain the computed value. Toward the neighbouring deuteron and toward the tetrahedral site the computed values are softer than the quartic.
- **NiH.** With the new [111] points the hydrogen spreads 0.4 to 0.9 percent more. Only one direction was computed, so how the softening is shared between directions is open and is now the larger uncertainty.
- One changed site in the coupled 32-site sample moves the same way by a similar amount. No sampled path reaches 0.8 A, so the coupled sample says nothing about the far tail.

So the single-site picture of `g10_quantum_site` stands, with the well now computed over the range the deuteron occupies. This is still one deuteron in a frozen lattice: not the pair cut, a surface site or a vacancy.

Archive SHA-256 aa45d70094240eb39245ad18aaa18e277ebb1b46ad2c2920f342f1ffdb3b6601; all file sums check. The lead reran `solve.py --variant wide_quadratic` in a separate tree built from the repository plus this archive (NumPy and SciPy only): the mean-square displacements and kinetic energies agree with the page to the eight digits printed. The other variant, the NiH fit and the coupled reweighting were not rerun. The archive's embedded copy of the `bulk-wide-wells` input is identical to the committed one and is not repeated here. One interpreter path quoted in `replay.md`, `prerequisites.json` and `package.py` was replaced by a description, so those three files do not match `SHA256SUMS`.
