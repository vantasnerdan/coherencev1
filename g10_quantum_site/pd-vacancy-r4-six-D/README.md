# Pd monovacancy with six D: relaxed 64-site cell

The fixed 4×4×4 primitive fcc Pd repeat (63 Pd plus six D) relaxed with neutral PBE, Pd10/H PAW, 500 eV plane waves and a 2×2×2 k-mesh to a largest force of 0.04690563 eV/Å, within the 0.05 eV/Å target. All six D remained near the vacancy. Their 15 distances are in `native/reference.json`; the shortest is 1.96607463 Å (D0–D3, effectively tied with D0–D1) and the longest is 3.28974836 Å. The nearest D pair across periodic vacancy images is 8.36357812 Å.

Against the corresponding relaxed 27-site reference in `native/small-cell/reference.json`, the largest labelled D–D change is 0.01220219 Å (D1–D3), within the specified 0.05 Å comparison target. The minimum increases from 1.95933232 to 1.96607463 Å. This checks this six-D basin at two cell sizes; it is not a convergence claim about other occupancies. The result does not require a 64-site two-D repeat under the agreed rule, and it does trigger a frozen pair cut because the larger-cell minimum is below 2.20 Å.

`sha256sum -c SHA256SUMS` verifies the two native records and the producer files. External literature and PAW files are identified by size and SHA-256 in the native records and are not bundled.
