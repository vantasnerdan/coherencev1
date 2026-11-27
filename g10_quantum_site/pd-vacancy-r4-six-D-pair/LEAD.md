# The closest pair in the palladium vacancy pushed together (lead's note)

Giuliano's capsule. In the relaxed 64-site cell of `../pd-vacancy-r4-six-D/` (one palladium atom removed, six deuterons in the empty cage) the two closest deuterons rest 1.966 A apart. Here they are pushed toward each other along the line joining them, each by the same distance, with every other atom held in place. The "extra" column is what the pair costs beyond the cost of moving each one alone with the other at rest. `README.md` is the technical page.

| Distance between the two (A) | Total energy rise (eV) | Sum of the two moved alone (eV) | Extra for the pair (eV) | Closest palladium (A) |
|---:|---:|---:|---:|---:|
| 1.966 (rest) | 0 | 0 | 0 | 1.727 |
| 1.846 | 0.016 | 0.014 | 0.001 | 1.719 |
| 1.726 | 0.061 | 0.055 | 0.005 | 1.688 |
| 1.486 | 0.231 | 0.207 | 0.024 | 1.632 |

- **Two deuterons 1.49 A apart cost 0.23 eV in the vacancy.** In bulk PdD (`../pdd-nearest-pair/LEAD.md`) the pair rests at 2.89 A, and reaching 1.93 A already costs 0.32 eV and 1.29 A costs 1.39 eV. The vacancy pair rests, at no cost, at a distance the bulk pair reaches only by paying 0.3 eV, and it then reaches 1.49 A for 0.23 eV more, where the bulk pair pays roughly 0.7 eV for the same step (read between its 1.93 and 1.29 A points).
- **The two barely act on each other.** The pair term is 0.024 eV at 1.49 A, a tenth of the total; in the bulk it was 0.074 eV at 1.93 A, about a quarter of the total. Nearly all the cost in the vacancy is each deuteron leaving its own resting place.
- **The pair term is repulsive at every distance computed** and grows steadily, as in the bulk. Nothing in this range lowers the cost of approach.

## Check on the mesh

At 1.486 A a finer k-point mesh (3x3x3 against 2x2x2, each with its own rest energy) moves the pair term from 23.6 to 24.5 meV (3.8 percent) and the total from 0.231 to 0.243 eV.

## Limits

- At rest and at every point a deuteron sits closer to a palladium atom than the sum of the setup radii, so every row carries the flag the two shortest bulk rows carried. In the bulk an all-electron hydrogen moved the pair term by 1 percent and the single-site cost by 0.08 eV at a larger displacement; that test was not repeated here.
- The other four deuterons and the palladium atoms are frozen. One line of approach, three distances, shortest 1.49 A.
- The two-deuteron vacancy was not repeated in the 64-site cell, as agreed.

## Lead's check

Archive SHA-256 7d141c91c70230943d744a9d12e07936051677a39582edf7ac99559456b3b21f; all file sums in both `SHA256SUMS` lists check. The lead ran the delivered NumPy readout on the native records: its output is byte-identical to `pair.json`. The first-principles runs were not rerun.
