# Zirconium calibration: does oxygen next to the target make a strong site?

Part (b) of `g7_micro/`, first stage. Kowalska et al. (2025) report enhanced screening in zirconium with oxygen present. This run asks whether a static calculation can see any of that: the same pair calculation as `g6_collision/` and `g7_hydride/`, for a deuteron in ordered ZrD2, with and without an interstitial oxygen atom beside it. Bare incoming deuteron.

## Result

Oxygen next to the target changes the static lowering by less than 1 eV at 0.1 A and by about 2 eV at contact. It lowers the screening slightly; it does not raise it.

| Settings | Clean, 0.1 A (eV) | With O, 0.1 A (eV) | Clean, contact (eV) | With O, contact (eV) |
|---|---:|---:|---:|---:|
| h 0.20 A, Gamma | 39.1 | 38.5 | 41.7 | 39.2 |
| h 0.16 A, Gamma | 39.9 | 39.4 | 41.9 | 41.3 |
| h 0.20 A, 2x2x2 k-points | 37.4 | 37.4 | 39.9 | 38.1 |

Contact values are the lowering at R = 0 (the delivered signed W0 with its sign reversed). The full curves from 0.1 to 1.6 A are in `curves.csv`; the oxygen curve lies 0.2-1.1 eV below the clean curve at every separation.

What it means:
- ZrD2 joins PdD, Pd at low loading, Pd vacancies, LiD, the D atom and D2: all give 37-51 eV at 0.1 A. An oxygen neighbour does not change that.
- The measured oxygen effect in Zr is therefore not a static property of a deuteron site with an oxygen atom beside it. If oxygen matters, it acts through something this cell does not hold: an oxide layer or interface, a different deuterium population near the surface, or the measurement itself (see `g7_gap/`).

## Checks

- Grid spacing and k-points move the oxygen-minus-clean difference by 0.05 and 0.58 eV at 0.1 A, and by 1.8 and 0.6 eV at contact. All are within the 5 eV tolerance.
- The lead recomputed the 0.1 A lowerings from the delivered pair energies and confirmed the differences in `paired-comparisons.json`.

## Three omissions

1. No cell-size check. The larger oxygen cell did not finish relaxing (forces reached 0.29 eV/A) and the larger clean cell failed at its contact step. The lead dropped the check: a sub-eV difference cannot grow to hundreds of eV with cell size.
2. One ordered ZrD2 proxy with one O per 8 Zr. It is not the mixed-phase experimental target, and not an oxide layer.
3. Static calculation with a bare probe only; the slope at contact was not computed.

## Lead's note

`delivered_README.md`, `paired-comparisons.json` and `Zr-delivery-manifest-001.json` are the delivered files, unchanged, and their SHA-256 values match the release assets. The delivered curve and contact files hold machine paths in their provenance blocks, so they are not in this public repository; `curves.csv` carries their pair energies. The 3.3 GB raw archive (SHA-256 3686cd12f43c3867dc2ab3623fb593f9d90009231c82cb6226e43d73b478829d) was not downloaded or replayed by the lead. Compute: 563 core-min of the 1500 core-min budget for part (b).
