# work: digitised data and first checks

- `fig2_digitized.csv`: Fig 2 marker centres for Pd, PdD, Ti and TiD, R(E) = Y(E)/Y(6.5 keV). Log-axis calibration on the published Fig 2 image: panel a x=0 at px 103.5, x=6.5 at px 1012; panel b x=0 at px 1073.5, x=6.5 at px 1982.5; y=1e0 at px 51.5, y=1e-4 at px 736.5.
- `si_tables.csv`: SI Table 1 (signal-to-background and cup current per energy), SI Table 2 (sigmoid parameters), SI Fig 3 IV points.
- `alr.py`: Bosch-Hale D(d,p)T S-factor, the paper's Eq 7 screened cross-section, thick-target integral with stopping proportional to sqrt(E).
- `table_vs_fig.py`: Fig 2 points against SI Table 1. Fig 2 tracks S:B/I^2 within 1.22-1.31 at every energy, S:B/I does not; the forward model carries both readings.
- `curvecheck.py`, `bias.py`: the solid Fig 2 curves are the cross-section ratio sigma(E,Ue)/sigma(6.5,Ue); a thick-target fit gives Ue about 0.4 keV higher. `curvecheck.py`, `axes.py` and `digit.py` need the Fig 2 image, which is not in the repo.

Run from `g4_beam/`: `python3 work/table_vs_fig.py`.
