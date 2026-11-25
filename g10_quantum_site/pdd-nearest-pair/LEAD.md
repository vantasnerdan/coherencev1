# Lead's note on the nearest-neighbour pair cut

First-principles energy of two neighbouring deuterons in fully loaded PdD pushed toward each other along the line joining them, each by the same distance, with every other atom held in place. Their distance at rest is 2.888 A. The "extra" column is what the pair costs beyond the cost of moving each one alone with the other at rest.

| Distance between the two (A) | Total energy rise (eV) | Sum of the two moved alone (eV) | Extra for the pair (eV) | Closest palladium (A) | Inside the nominal overlap range of the setups |
|---:|---:|---:|---:|---:|---|
| 2.648 | 0.012 | 0.009 | 0.003 | 1.959 | no |
| 2.408 | 0.054 | 0.041 | 0.013 | 1.880 | no |
| 2.168 | 0.146 | 0.112 | 0.034 | 1.805 | no |
| 1.928 | 0.316 | 0.242 | 0.074 | 1.736 | yes |
| 1.288 | 1.390 | 1.014 | 0.376 | 1.581 | yes |

- **Two deuterons 1.29 A apart cost 1.39 eV**, of which 1.01 eV is the cost of each leaving its own site and 0.38 eV is the two acting on each other. The pair term is repulsive at every distance computed and grows steadily; nothing in this range lowers the cost of approach.
- This joins the range of the single-site wells (down to about 2.1 A between neighbours) to the outermost point of the earlier fast-ion pair calculations (1.6 A).
- **Toward the empty tetrahedral site** ([111], one deuteron, others fixed): 0.301, 0.530, 0.678 and 0.813 eV at 0.80, 1.10, 1.40 and 1.768 A from the site, the last being the tetrahedral position itself. The force still points back at every point, so with the neighbours frozen there is no barrier and no second minimum on this line.

**Not checked yet.** The points use the `g3_loading` settings (3x3x3 k-point mesh). The finer-mesh run at the closest point and the repeat with an all-electron treatment of the hydrogen did not finish; their reference rows are in the archive and the comparisons are blank. At the two shortest distances, and at all four [111] points, a deuteron sits closer to a palladium atom than the sum of the setup radii (1.799 A), so those rows carry a flag. The supplement with both comparisons is requested.

Archive SHA-256 716978d9797e12683a937a5f1fe79c715b4d3d000e8ac4da5614cc2f25ef7777; all file sums check. The lead recomputed the table and the four [111] energies directly from the energies in `pair-native/`, `bulk-native/` and `references/Pd_zero.json`: all agree with `pair-cut.csv` and the README to 5e-12 eV. The delivered replay scripts were not run.
