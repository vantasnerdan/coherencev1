# g7_defect: screening at a vacancy and at low loading in Pd

The question: published data show screening rising with cold work and falling loading, while our perfect-crystal PdD gives about 40 eV. Which of the two, defects or low loading, raises the static screening, and by how much? The target is about 1 keV for D-D, 25 times the perfect-PdD value.

## Result

Static Born-Oppenheimer lowering at R = 0.1 A (same method as `g6_collision/`, all cells at the PdD lattice constant, positions relaxed):

| Case | Lowering at 0.1 A (eV) |
|---|---:|
| A: Pd with one D (low loading) | 37.2 |
| B: perfect PdD (from `g6_collision/`) | 40.5 |
| C: one D trapped at a Pd vacancy, low loading | 38.3 |
| D: Pd vacancy in PdD (high loading) | 37.9 |

Neither variable moves the screening by more than 10 percent. Loading changes it by 1.09 times in the perfect cell. A vacancy changes it by 1.03 times at low loading and 0.93 times at high loading. The full curves from 0.1 to 1.6 A are in `capsule/README.md` and are exported in the `g6_pairlaw` pair schema under `capsule/curves/`.

In D the target D did not stay in the vacancy: it relaxed 2.4 A away, so D is a high-loading vacancy cell, not an occupied vacancy-hydrogen complex.

What it means: static electron screening of a deuteron in Pd stays at 37-41 eV across loading and a vacancy. The 5.7-times cold-work effect and the Ti loading trend in the published data are not reproduced by static screening at a vacancy or at low loading. The dislocation-core case was not run, because nothing here came near the 3-times trigger. The keV must come from something these static cells leave out: the near-surface region, motion during the collision beyond the 2 keV point already computed, or the beam-side causes in `g6_plateau/`.

## Checks

- Grid spacing (h 0.20 to 0.16 A) on C: +0.43 eV, within max(5 eV, 10 percent).
- Cell size: a converged check was not run. A 27-site cell with the C geometry embedded and the far field unrelaxed gives 38.4 eV (+0.10 eV). This is a sensitivity estimate only.
- Forces: below 0.1 eV/A in every case (A 0.046, C 0.079, D 0.064). Geometries were relaxed with the standard H setup; all energies were then recomputed with the all-electron H setup.

## Three omissions

1. Converged cell size and k sampling; vacancy density is high in the 8-site cell.
2. The smaller lattice of low-loading Pd (a = 3.89 A); all cells use the PdD lattice.
3. Static BO only: no thermal ensemble, zero-point motion or collision dynamics, and one local basin for each case.

## Lead's note

`capsule/` is the delivered archive `defect-screening-delivery.tar.gz` (SHA-256 878f69b11cb238d1b763a96e0b0f2f4bfdc9d49b7602c62c2aff0a490bff5228), unchanged. Replay with `python g7_defect/capsule/reader.py g7_defect/capsule` (NumPy only). The reader compares the replayed readout for exact equality. Under NumPy 1.26.4 it refuses on floating-point roundoff; the lead's field-by-field comparison found no non-numeric differences and a largest numeric difference of 7.1e-15 eV. The producer's replay under NumPy 2.5.3 is exact. Compute used: 287 core-min of an 800 core-min budget.
