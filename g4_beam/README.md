# g4_beam: low-energy D-D beam fusion in Pd and Ti

A forward model of the relative D(d,p)T yields that Karahadian et al. report for Pd, PdD, Ti and TiD at centre-of-mass energies 0.25 to 6.5 keV (Nat. Commun. 17:8845 (2026), doi:10.1038/s41467-026-74421-1; source data on Zenodo, doi:10.5281/zenodo.20435774, open from 2026-12-01). Alongside it, a near-surface deuterium transport model turns material knobs (alloying, strain, temperature, overlayers, loading current) into the deuterium profile the beam sees.

## Place in the map

Outside the map. The beam does not change the lattice. This chain reads three quantities from the material layer:
- the static dielectric response eps(q), which sets ratio 7 (screening length over lattice spacing);
- the deuterium loading in the top 100 nm;
- the metal vacancies the beam makes, which trap deuterium. These are not the empty deuterium sites of ratio 6.

Ratio 7 pulls in opposite directions in the two chains:
- **Lattice map:** the static ordering branch in `g1_rank/` needs a screening length of at least 1.09 A (PdD) or 1.55 A (NiH).
- **Beam yield:** a static screening length that gave a 1 keV screening energy would be about 0.014 A (U_e = e^2/lambda).

One static number cannot serve both, so the shared deliverable is eps(q) over the full q range from first principles. The ordering branch reads it at q of order 1 per A, and the fusion screening reads it at much larger q.

## Layout

| Directory | What it is |
|---|---|
| `work/` | Digitised Fig 2, SI Tables 1 and 2, SI Fig 3 IV points; first checks |
| `forward/` | Thick-target yield model with a measurement layer, fitted to Fig 2; independent second route; callable consumer |
| `czerski/` | Czerski et al. (Pd, Ti, Zr) and Tohoku (Ti, EXFOR O1084) fitted with the same candidate set; shared thermal amplitude |
| `transport/` | 1D deuterium transport with beam damage, binary-collision ranges, knob sweeps |
| `screening/` | 184 measured screening potentials with sources, screening-model formulas, correlation tests |
| `screen/` | First-principles screening (GPAW, PBE): density of states and d-band centres for 14 hosts, static RPA screening length for 9, the neutral adiabatic screening cycle for Pd and Pd3Rh, and ratio 7 |
| `priordata/` | Other low-energy d+d datasets in one table, plateau comparison |

Run every script from its own directory, except `work/` scripts and `priordata/scripts/compare.py .`, which run from `g4_beam/`.

## Current best numbers

Every fit now carries a fixed threshold-resonance term (Czerski et al., arXiv:2409.02112v3: Gamma_p 40 meV, phase 115 degrees, 7 fm) as common reaction physics; see `forward/README.md` and `forward/bulk_screening.csv` for the Eq 7-only values beside these.

Bulk screening potential U_e (keV) from a thick-target fit to the 2.5-6.5 keV points (28 points). Fig 2 and SI Table 1 differ by one factor of cup current, so the model carries the charge normalisation as a hypothesis:

| Charge normalisation | chi2 | Pd | PdD | Ti | TiD |
|---|---|---|---|---|---|
| As plotted (M0) | 37.85 | 1.25 | 2.04 | 0.98 | 1.32 |
| One power of cup current (M1, M2, M4) | 19.9-20.0 | 0.70 | 1.40 | 0.45 | 0.84 |
| Two powers (M3, M5) | 25.58 | 0.25 | 0.81 | 0.03 | 0.45 |

Each power of current moves U_e by 0.39-0.64 keV. The step loading adds stays at PdD-Pd 0.55-0.79 keV and TiD-Ti 0.34-0.42 keV under all six hypotheses.

Low-energy plateau, best fit of each candidate on all 48 points (`forward/fits.json`):

| Candidate | Hypothesis | k | chi2 | Key parameters |
|---|---|---|---|---|
| Sigmoid screening inside the thick target | M3 | 10 | 36.92 | surface U_e 1.29-1.69 keV, E0 1.63 keV |
| Sub-monolayer surface layer | M2 | 10 | 40.29 | surface U_e 3.3-5.9 keV, depth 0.052 nm |
| Fast neutrals, fitted energy | M2 | 7 | 46.31 | 3.8 keV lab, flux 1.59 x extracted current |
| Fast neutrals at 5 keV lab | M4 | 6 | 49.65 | flux 0.267 |
| Thermal spike, one amplitude | M4 | 6 | 132.23 | |
| No additional plateau | M0 | 5 | 181.45 | |

The first three lie within 3.4 in AIC. With the fixed resonance in place the no-plateau control misses the lowest points by factors of 4 to 5; without it (Eq 7 only) it had chi2 459.

Electronic screening from first principles (`screen/README.md`, valence-PAW PBE, static low-q RPA). The screening length lambda gives a Yukawa static U_e = e^2/lambda:

| Host | N(E_F) per metal atom (/eV) | Pd d-band centre (eV) | lambda (A) | static U_e (eV) | ratio 7, lambda/D-D |
|---|---|---|---|---|---|
| Pd -2% | 1.937 | -2.684 | 0.391 | 36.85 | 0.145 |
| Pd | 2.055 | -2.424 | 0.397 | 36.23 | 0.145 |
| Pd +2% | 2.169 | -2.193 | 0.404 | 35.68 | 0.144 |
| Pd3Rh | 1.848 | -2.547 | 0.404 | 35.62 | 0.146 |
| Pd3Ag | 1.357 | -2.074 | 0.478 | 30.10 | 0.169 |
| PdD | 0.497 | -2.684 | 0.541 | 26.61 | 0.185 |
| Ni | 1.276 | -2.226 | 0.376 | 38.32 | 0.155 |
| NiH | 0.727 | -2.214 | 0.451 | 31.90 | 0.175 |
| TiD2 | 2.452 | -2.588 | 0.435 | 33.08 | 0.197 |

(The alloy rows give the Pd-projected d-band centre; the README also lists all-metal centres.) Across the Pd strain series and the two alloys, static U_e changes by -8.35 eV per eV of Pd d-band shift (interval -22.3 to +5.4 under a +/-10 % allowance on lambda). The neutral adiabatic cycle gives 42.25 eV for Pd and 42.26 eV for Pd3Rh, a slope of -0.066 eV/eV. Loading Pd to PdD lowers static U_e from 36.2 to 26.6 eV.

Across laboratories (`czerski/README.md`): the palladium plateau can share layer or sigmoid parameters between Szczecin and UC Davis (mismatch factors 1.2-1.4); titanium cannot once the Tohoku data are included. One thermal-spike amplitude for both laboratories fails: A_UC Davis / A_Szczecin = 622.5 (as plotted) and 283.8 (one cup-current power).

Near-surface deuterium (`transport/README.md`, `transport/knobs.csv`). In Pd the top 5 nm hold about 0.01 D per metal atom, almost all in beam-made vacancies; Ti and TiD both reach TiD ~1.85 under the beam. Factors on D/M in the top 5 nm:

| Knob | Pd / PdD | Ti / TiD |
|---|---|---|
| Pd77Ag23 / Pd90Pt10 | x6.5-12 / x0.27-0.31 | - |
| Pd80Rh20 (estimate from H-electrosorption data) | x0.31 at 0.25 keV, yield x0.96 | - |
| -30 K / +30 K | x2.3-2.5 / x0.40-0.63 | ~x1 |
| Strain +0.5% / -0.5% | x1.4 / x0.68 | x1.00 |
| Back-side current x2 | x1.13-1.14 (end of run) | x1.00 |
| 2 nm Al2O3 while it lasts | x0.04 (0.25 keV) to x3.3 (PdD, 2.5 keV) | Al2O3 5 nm x0.19-0.30 |

A beam-side surface that almost stops recombining deuterium (sticking 1e-7) lets the hydride reach the beam side (0.68-0.72 D/Pd). First-principles inputs: the converged 1 ML adsorption of D on clean Pd(111) is 0.466 eV per D below half a D2 molecule, beside the model's 0.45 eV; the Pd vacancy binding stays at the 0.23 eV estimate (8-Pd pilot 0.131 eV, unrelaxed 27-Pd trial 0.242 eV); the D2 recombination barrier was not converged. Pd80Rh20 lowers near-surface D but barely moves the yield, because the adsorbed monolayer stays nearly full.

## What the models keep, drop, and leave out

`forward/` keeps the Bosch-Hale S-factors, ICRU-49 electronic and ZBL nuclear stopping with Bragg additivity, the 85/15 D+/D2+ beam, and a measurement layer: charge normalisation, Faraday-cup collected fraction in each beam-delivery regime, proton-to-neutron stitch factor, and fast neutrals. It drops target evolution, straggling and backscatter, and the deuterium depth profile.

`transport/` keeps the isotherm-driven diffusion over the full foil, vacancy trapping, surface desorption against chamber gas, sputter erosion, overlayers and beam heating. It drops microstructure, blisters and run history.

The three biggest omissions are:
1. Beam-side surface chemistry, which sets the adsorbed deuterium that carries 34-94% of the modelled Pd yield.
2. The real charge on the foil in each run.
3. Screening below about 50 eV.

## Open work

Open work is tracked on Fleet nodes, one directory per node when delivered:
- `g4_zenodo`: refit on the per-run data once it opens.
