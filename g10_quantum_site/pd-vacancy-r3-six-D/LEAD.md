# Six deuterons in a palladium vacancy, small cell (lead's note)

Giuliano's capsule. One palladium atom is removed from a 27-site cell and six deuterons are placed in the empty cage, then every atom is let go until the largest force is 0.045 eV/A. `README.md` is the technical page.

## What it shows

All six deuterons stay in the vacancy. Their fifteen mutual distances at rest:

| Distance | Number of pairs |
|---:|---:|
| 1.96 A | 4 |
| 1.98 A | 2 |
| 2.49 A | 1 |
| 2.64 to 2.77 A | 6 |
| 3.29 A | 2 |

Six pairs rest closer than 2.0 A. For comparison, neighbours rest 2.89 A apart in bulk PdD, 2.44 A for the closest pair at the covered PdD(111) face and 2.21 A in cubic TiD2. Each deuteron sits 1.71 to 1.74 A from its nearest palladium atom.

## Assumptions and limits

One relaxed arrangement from one starting arrangement in the small cell; not shown to be the lowest one. The larger-cell repeat is in `../pd-vacancy-r4-six-D/`: the distances survive to within 0.012 A.

## Lead's check

Archive SHA-256 ce6bc94c40711b410643a3844d484865684dbe611fbf59abbdb4fcbc21f1fb0f, all file sums checked as delivered. The lead recomputed all fifteen deuteron distances and the deuteron to palladium distances from the stored positions with a separate script (they match). The first-principles relaxation was not rerun.

One line of `replay/surface_vacancy.py` was changed by the lead: a default input folder that named a location on the machine that ran it is now the current folder. Its two lines in the sum lists were regenerated; `SHA256SUMS.delivered` is the delivered list. The native records still cite the delivered file's hash.
