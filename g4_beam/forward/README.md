# forward: thick-target yield model with a measurement layer

`model.py` gives the yield per incident deuteron for each low-energy candidate; `fits.py` fits every candidate under six charge-normalisation hypotheses (M0-M5) to the 48 digitised points with error bars (`errbars.csv`), writing `fits.csv` and `fits.json` (about 6 CPU-minutes). `extra_checks.py` and `analysis.py` give the per-energy Ue each candidate implies (`ue_implied.csv`), absolute counts per coulomb (`abs_yield.csv`) and ranges (`ranges.csv`).

Keeps: Bosch-Hale D(d,p)T and D(d,n)3He S-factors; the paper's Eq 7; thick-target integral from 0.1 keV lab; ICRU-49 electronic stopping (Pd A1 5.238, Ti 4.858, H 1.254 eV/(1e15 at/cm2)/sqrt(keV/amu)) plus ZBL nuclear stopping, Bragg additivity for hydrides; 85% D+ / 15% D2+, 5% FWHM energy spread. Measurement layer: charge per coulomb or per coulomb x current, Faraday-cup collected-fraction exponents above 5.5 kV and below 5 kV from the SI Fig 3 IV curve, one proton-to-neutron stitch factor applied at E <= 4.0 keV, fast neutrals at a fixed energy that the cup does not count.

Low-energy candidates: the paper's sigmoid (as a cross-section ratio, and inside the thick-target integral), a step Ue(z) surface layer, a Czerski-type threshold resonance, a free floor.

Drops: target evolution, straggling, backscatter, depth profile of D (uniform), neutral energy spread, correlated error of the 6.5 keV reference.

`ranges.csv` uses a first-moment projected-range estimate; the binary-collision ranges in `../transport/` (Pd 7.3 nm at 0.5 keV lab, 83 nm at 13 keV) replace it. Only the surface-layer thickness numbers depend on it.

Three biggest omissions: target evolution during 12-24 h runs; the charge actually delivered to the foil in each run and the neutral energy distribution; the D depth profile and screening below about 50 eV.
