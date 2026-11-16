# screening: measured screening potentials and screening models

- `ue_table.csv`: 184 measured d+d (and converted) screening potentials, with material, D content where known, method and source.
- `models.md`: formulas and numbers for Debye, Thomas-Fermi, Lindhard, the Czerski dielectric model, the dynamic correction, threshold resonance, thermal spike, the Kim-Zubarev momentum tail and defect pictures.
- `analysis.py`: Spearman tests of Ue against N(E_F), work function, d-band centre, d count, r_s, hydride formation, group number and Hall density (`correlations.csv`, `model_numbers.txt`). `density_vs_ue.py` gives how much a change in near-surface D density mimics a change in Ue.

Results:
- Static screening models give 25-135 eV for these metals. Reaching 0.3-3 keV needs a static screening length of 0.005-0.05 A.
- Across 32 metals with a finite Ue, neither d-band centre (rho 0.09, n 11) nor N(E_F) (rho 0.04) tracks Ue.
- D content does track it, negatively (rho -0.80; hydride formers sit low).
- Ratio 7 from static Thomas-Fermi: Pd 0.07-0.21, PdD0.7 0.17-0.19, Ti 0.11-0.20, TiD2 0.07-0.16.

Three biggest omissions: the primary per-metal D/M and Hall data (taken from Huke 2008); no computed electronic structure; no survey that pairs Ue with a measured near-surface D profile.
