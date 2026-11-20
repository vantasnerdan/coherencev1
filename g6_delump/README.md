# g6_delump: energy return with the pass 6 inputs

The pass 6 rerun of the energy return (fusion energy out over beam energy in) for a deuteron beam into a target that holds tritium. It has two parts: the trade-off between fuel fraction and screening in a metal host (`fuel_fraction.py`), and one table that puts every target, screening level and the channelling gain on the same footing (`combined.py`).

## Fuel fraction against screening

`fuel_fraction.py` (run `python -m g6_delump.fuel_fraction` from the repository root) treats a metal M holding tritium at T/M = x. It writes `break_even_local.csv`, `return_grid.csv`, `ti_series.csv` and `fuel_fraction.json`.

Diluting the tritium has two effects. The metal's drag is shared by fewer T nuclei. But measured screening rises as loading falls: in Ti, D-D screening is under 30 eV at D/Ti 2.1 and 250-295 eV at 0.2-0.26 (Raiola et al. 2005, EXFOR O1379).

Inputs are the same as in `g6_gap/`:
- ICRU-49 velocity-law electronic stopping and ZBL nuclear stopping (`g4_beam/forward/model.py`);
- the Bosch-Hale D-T cross section in the Eq 7 screened form;
- measured screening carried from D-D to D-T unchanged.

The velocity law holds below 20 keV lab, so beam energies stop at 13 keV. The script checks that it reproduces the 206 b break-even for TiT1.7 at 1 keV.

**What break-even needs.** Local break-even means the fusion energy gained per T nucleus passed equals the drag per T nucleus.

| Target | Cross section needed at 1 keV lab | Screening that delivers it |
|---|---|---|
| PdT (x = 1) | 314 b | none |
| TiT1.7 | 206 b | none |
| Metal at x = 2 | 184-188 b | none |
| Pure tritium, metal drag law | 62 b | 126 keV |
| LiT, computed drag (`g6_gap/`) | 52 b | 107 keV |

In the Eq 7 form, screening raises the D-T cross section only up to the point where the whole barrier is removed, at U_e of about 1.18 MeV. That ceiling is 159 b at 1 keV, 115 b at 2 keV, 76 b at 5 keV and 52 b at 13 keV. For a metal host to reach break-even, the drag per T nucleus must therefore fall by:
- 1.2 to 2 times at 1 keV (x = 2 to 1);
- 11 to 19 times at 13 keV.
Diluting the fuel raises the cut needed (8 times at 1 keV for x = 0.2). Hosts already under the ceiling are the low-drag fuels: LiT, and pure tritium. Pure tritium's drag is 81 percent electronic at 1 keV, so a gapped solid T2 would sit lower than the metal law gives.

**What dilution buys.** The thick-target return for Ti at the measured low loadings (x = 0.2, U_e 290 eV; x = 0.26, 250 eV), relative to TiT1.7 at 30 eV:

| Beam energy (keV) | 1 | 2 | 3 | 5 | 13 |
|---|---|---|---|---|---|
| x = 0.2 | 204 | 3.3 | 0.94 | 0.39 | 0.20 |
| x = 0.26 | 108 | 2.8 | 0.94 | 0.44 | 0.25 |

Dilution pays below 2.9 keV, where the screening gain beats the extra drag per T nucleus. Above that the drag wins. This table integrates down to the fit floor (0.83 keV lab). Integrating down to 0.1 keV instead, as `combined.py` does, raises the 1 keV ratios to 225 and 117 and leaves the 2.9 keV crossover unchanged.

**Levers, with the number each must reach:**
- Drag per fuel nucleus, the binding spec. It must fall under the ceiling: below 2.8 eV per 10^15 T/cm^2 at 1 keV. Insulating, fuel-dense hosts (LiT, solid T2) are already there and then need U_e of about 107-126 keV.
- Screening, in a host under the ceiling. Measured D-D values in metals at low loading are 0.3-0.8 keV (`g4_beam/priordata/`), and the UC Davis fits give 0.25-2 keV, against the 107 keV target. `g7_defect/` finds that neither a vacancy nor low loading does: static screening stays at 37-41 eV in all four cells.
- Fuel fraction. A metal host gains from dilution only below about 3 keV beam energy.

## Combined energy return

`combined.py` (run `python -m g6_delump.combined` from the repository root) gives the thick-target energy return G for a deuteron beam stopping in each target. It writes `combined.csv` and `combined.json`. Drag is per T nucleus in every row:
- TiT1.7, PdT and Ti at T/Ti = 0.2 use the same stopping as above.
- LiT uses the computed electronic drag from `g6_gap/` (q4) plus ZBL nuclear drag on Li and T. The script checks that it reproduces the 52.3 b break-even at 1 keV.
- Pure tritium uses ICRU-49 hydrogen electronic stopping plus ZBL on T.
- "T as only scatterer" keeps only ZBL nuclear drag on T. It is the floor for a host with no electronic drag.

Three screening levels are carried from D-D to D-T unchanged: 0.04 keV (computed static lowering, `g6_collision/` and `g7_defect/`), 0.3 keV (measured D-D in metals at low loading) and 1.4 keV (the uniform UC Davis PdD fit). Channelling multiplies the PdT return by 1.78 at the one point where `g6_channel/` computed it: Pd(D/T) `<110>`, 13 keV, 2 degrees. That gain was computed without screening and is applied to the screened rows unchanged. The slowing-down integral runs from the beam energy down to 0.1 keV lab, the `g6_gap/` convention. Below 0.83 keV lab the Bosch-Hale fit is an extrapolation.

**Screening needed for G = 1**, keV. "None" means that even removing the whole barrier leaves G below 1; the largest G reachable is given in brackets.

| Target | 1 keV beam | 2 keV | 5 keV | 13 keV |
|---|---|---|---|---|
| TiT1.7 | 297 | none (0.96) | none (0.55) | none (0.29) |
| PdT | none (0.94) | none (0.66) | none (0.37) | none (0.19) |
| PdT, channelled `<110>` | | | | none (0.34) |
| Ti, T/Ti 0.2 | none (0.23) | none (0.16) | none (0.09) | none (0.05) |
| LiT, computed drag | 102 | 116 | 202 | none (0.84) |
| Pure tritium | 87 | 111 | 196 | none (0.92) |
| T as only scatterer | 45 | 42 | 38 | 31 |

**Return at the fitted screening.** At U_e = 1.4 keV, G is 8e-11 to 2e-5 for every target with electronic drag. It is 1.7e-3 at 13 keV for T as only scatterer, 580 times short of G = 1. At the computed static 0.04 keV, G is lower again: 3 times lower with a 13 keV beam and up to 5e8 times lower with a 1 keV beam. Channelling raises the PdT return at 13 keV from 4.1e-6 to 7.3e-6 at 1.4 keV screening.

**Integration floor.** The floor matters only where tens of keV of screening are needed. Stopping at 0.2 keV instead of 0.1 keV raises the screening needed for LiT and pure tritium at 1 to 2 keV by 8 to 19 percent, and for T as only scatterer by under 7 percent. Stopping at the fit floor (0.83 keV) removes most of a 1 keV beam's path and is not used. Both floors are in `combined.csv`.

**Levers, with the number each must reach:**
- Host. No metal tritide reaches G = 1 at any screening, except TiT1.7 with a 1 keV beam at 297 keV of screening, which depends on the integration floor. Channelling does not change this. The fuel must sit in a host with no metal electrons on the beam's path.
- Electronic drag. With it, as in LiT or pure tritium, a 1-2 keV beam needs 87-116 keV of screening, and no screening gives G = 1 at 13 keV. Without it, the need falls to 31-45 keV at every beam energy from 1 to 13 keV, and the largest reachable G is 12-21. A band gap cuts electronic stopping for slow ions, so a gapped, fuel-dense solid such as solid T2 lies between these two rows at 1-2 keV.
- Screening in that host: 31-45 keV, against 0.3-2 keV measured in metals and fitted to the UC Davis data, and 0.04 keV computed statically.

