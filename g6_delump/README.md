# g6_delump: energy return with the pass 6 inputs

The pass 6 rerun of the energy return (fusion energy out over beam energy in) for a deuteron beam into a target that holds tritium. Its first piece is here: the trade-off between fuel fraction and screening in a metal host. The remaining pieces will be added to this directory as they are rerun: channelling, LiT with computed drag, and the defect screening from `g7_defect/`.

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

Dilution pays below 2.9 keV, where the screening gain beats the extra drag per T nucleus. Above that the drag wins.

**Levers, with the number each must reach:**
- Drag per fuel nucleus, the binding spec. It must fall under the ceiling: below 2.8 eV per 10^15 T/cm^2 at 1 keV. Insulating, fuel-dense hosts (LiT, solid T2) are already there and then need U_e of about 107-126 keV.
- Screening, in a host under the ceiling. Measured D-D values in metals at low loading are 0.3-0.8 keV (`g4_beam/priordata/`), and the UC Davis fits give 0.25-2 keV, against the 107 keV target. `g7_defect/` finds that neither a vacancy nor low loading does: static screening stays at 37-41 eV in all four cells.
- Fuel fraction. A metal host gains from dilution only below about 3 keV beam energy.
