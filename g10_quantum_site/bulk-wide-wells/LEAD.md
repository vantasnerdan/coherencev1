# Lead's note on the wide bulk wells

First-principles energy and force for one deuteron (or hydrogen) pushed away from its site in fully loaded PdD and NiH, out to 0.8 A, with all other atoms held in place. The earlier wells stopped at 0.24 A.

| Material and direction | What lies that way | Energy rise at 0.48 A (eV) | at 0.8 A (eV) | Restoring force at 0.8 A (eV/A) |
|---|---|---:|---:|---:|
| PdD [100] | a palladium atom | 0.244 | 2.035 | 11.8 |
| PdD [110] | a neighbouring deuteron | 0.121 | 0.507 | 1.72 |
| PdD [111] | the empty tetrahedral site | 0.093 | 0.301 | 0.79 |
| NiH [111] | the empty tetrahedral site | 0.197 | 0.501 | 0.99 |

- The well is strongly direction-dependent once the deuteron is far out: a wall toward the metal atom, and a soft way out toward the empty tetrahedral site. A simple spring with the site stiffness would give 0.19 eV at 0.8 A in every direction.
- Along [111] the force grows slowly (0.46 to 0.79 eV/A between 0.48 and 0.8 A in PdD, 0.81 to 0.99 in NiH), which is how a path toward a barrier and a second site looks. The barrier itself is beyond these points.
- Along [110] at 0.8 A the deuteron is 2.09 A from its neighbour, and that costs 0.51 eV with the neighbour held still.

Positive and negative displacements were computed separately along [110] at three distances and agree to 1e-5 eV and 3e-4 eV/A. The other negative rows are the positive ones mirrored and are marked so.

Archive SHA-256 8085c5623dfbbf137a097249e6637c42182565dcb3760ee544e36df1f90553f4; all file sums check. There is no script to rerun; the table above is read directly from `bulk-wells.json`, and the spring comparison and the 2.09 A distance are the lead's arithmetic. The points use the `g3_loading` settings (3x3x3 k-point mesh); no second mesh was run for these wide points.
