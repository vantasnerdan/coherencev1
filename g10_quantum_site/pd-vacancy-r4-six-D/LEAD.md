# Six deuterons in a palladium vacancy, larger cell (lead's note)

Giuliano's capsule. The same arrangement as `../pd-vacancy-r3-six-D/` (one palladium atom removed, six deuterons in the empty cage), repeated in a 64-site cell and relaxed until the largest force is 0.047 eV/A. `README.md` is the technical page.

## What it shows

The short distances survive the larger cell. All six deuterons stay in the vacancy.

| Distance (64-site cell) | Number of pairs | Same pairs in the 27-site cell |
|---:|---:|---:|
| 1.97 A | 4 | 1.96 A |
| 1.97 A | 2 | 1.98 A |
| 2.49 A | 1 | 2.50 A |
| 2.66 to 2.78 A | 6 | 2.64 to 2.77 A |
| 3.29 A | 2 | 3.29 A |

No distance moves by more than 0.012 A between the two cells. Six pairs rest closer than 2.0 A, against 2.89 A for neighbours in bulk PdD, 2.44 A at the covered PdD(111) face and 2.21 A in cubic TiD2. Each deuteron sits 1.73 to 1.75 A from its nearest palladium atom. The nearest deuteron of the next vacancy image is 8.4 A away.

## What follows

The two-deuteron vacancy is not repeated in the larger cell (the six-deuteron distances moved by less than 0.05 A). The energy cut for the closest pair in this cell is in `../pd-vacancy-r4-six-D-pair/LEAD.md`.

## Assumptions and limits

One relaxed arrangement from one starting arrangement, at two cell sizes; not shown to be the lowest one, and nothing is said about other fillings.

## Lead's check

Archive SHA-256 df29d40aa7c7c18e2750cf355c5ba231b22c86fa3dccfd2ed2e2fd254ee4cb53, all file sums checked. The lead recomputed all fifteen deuteron distances and the deuteron to palladium distances in both cells from the stored positions with a separate script (they match), and confirmed that the small-cell record here is byte-identical to the one committed in `../pd-vacancy-r3-six-D/`. The first-principles relaxation was not rerun.
