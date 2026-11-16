# priordata: other low-energy d+d data in metals

- `experiments.csv`: one row per dataset (group, target, energy range, beam, vacuum, charge measurement, detection, lowest energy, plateau or not, source and how it was read). `scripts/make_experiments_csv.py` writes it.
- `scripts/compare.py` (run from `g4_beam/` as `python3 priordata/scripts/compare.py .`): plateau levels relative to each dataset's own 3.0 and 6.5 keV yields, the Tohoku Ti thick-target ratio (EXFOR O1084) against the bare and screened model, and the beam fraction a fixed-energy component would need (`compare_out.txt`).
- `poster_neutron_rates.csv`: neutron rates digitised from the same group's APS March Meeting 2025 poster (OSTI 2569797).

Main comparison: Czerski et al. (arXiv:2409.02112v3, 2605.27438) see a low-energy plateau in Pd at a similar level relative to the 3 keV yield (0.44 against 0.22 here). Relative to each target's own 6.5 keV yield, their Ti and Zr plateaus are about 1e-3 of their Pd value; here Ti's is about 0.3 of Pd's. Tohoku Ti (Yuki et al. 1997), Kasagi et al. 2002 and the Greife et al. 1995 gas target follow the bare or weakly screened cross section down to 1.25-1.6 keV.

Data still to digitise: Czerski arXiv:2409.02112v3 Figs 2, 3, 4, 7 and Table 2; Schenkel et al. 2019 Fig 9; Kasagi et al. 2002 excitation function.
