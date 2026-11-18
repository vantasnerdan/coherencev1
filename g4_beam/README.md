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
| `merit/` | Pass 5: yield of any layered target relative to the transported UC Davis PdD foil, side by side with the UC Davis values (`compare.csv`); causes refitted on transported profiles |
| `causes/` | Pass 5: how each plateau cause responds to beam and target changes, surface-density and sheet-occupancy scans, Zenodo column list |
| `transport/recomb/` | Pass 5: first-principles D sheet on clean Pd(111), no-beam coverage and D2 release, connected clean recombination route |
| `ald/` | Pass 5: 245 ALD designs in six families (sealing oxides, deuteride caps, deuterated oxides, donor promoters, BZY proton conductor, Pd/oxide nanolaminates), with run-integrated merit, erosion and cap sputter thresholds |
| `pairs/` | Pass 5: D clusters in a Pd vacancy (8-Pd cells, 1 to 6 D) and screening by mobile deuterons, through the merit consumer |

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

The first three lie within 3.4 in AIC. With the fixed resonance in place the no-plateau control misses the lowest points by factors of 4 to 5; without it (Eq 7 only) it had chi2 459. These fits use uniform targets. Pass 5 refits them with the transported deuterium profiles; see below.

Pass 5 (`merit/`): the three causes refitted with the transported, run-end deuterium profiles of all four foils plus the adsorbed sheet, the reference for every pass 5 design. Best chi2 over the 0, 1 and 2 extra current powers: rising screening 51.48, fast neutrals 69.89, surface layer 262.24. The surface layer no longer fits once the real profile is used. The refitted model gives PdD/Pd 4.9-7.2 over 0.25-3 keV against the SI Table 1 counts per charge of 0.98-2.13.

Merit is the D(d,p)T yield of a design over the yield of the transported UC Davis PdD foil at the same energy:

| | 0.25 keV | 1 keV | 2.5 keV | 6.5 keV |
|---|---|---|---|---|
| Extra uniform U_e for 10x, rising screening (keV) | 0.51 | 0.53 | 0.78 | 2.82 |
| Extra uniform U_e for 10x, fast neutrals (keV) | 0.99 | 0.86 | 0.95 | 2.29 |
| Ideal beam side at 0.7 D/Pd with the sheet kept, rising screening | 136 | 41 | 20 | 29 |
| The same, fast neutrals | 35 | 21 | 21 | 30 |

Against the uniform 0.7 D/Pd target instead, the same ideal beam side gives 1.02-1.05.

Pass 5 (`causes/`): the adsorbed D sheet decides which cause fits.
- **Without the sheet:** refitting the surface cause gives chi2 44.21, the best shape fit, while rising screening falls to 110 and neutrals to 159.
- **Target-side scans:** no setting reconciles the curve shapes with the SI Table 1 counts. Scaling near-surface bulk D (0.1-100, shared and per foil) and varying sheet occupancy (0-1 ML per foil) were tried for all three causes. The canonical reference gives a model/measured PdD/Pd of 2.55-6.83, and the closest variants 1.5-3. The charge per run is the next hypothesis, for `g4_zenodo`.
- **Ideal beam side at 0.7 D/Pd with no sheet** (rising screening, best joint fit): 137/41/22/32 at 0.25/1/2.5/6.5 keV. Against the zero-sheet surface reference it gives 53/47/45/47.

Measurements that separate the causes, as changed/base counts per charge:

| Change | Rising screening | Surface layer | Fast neutrals |
|---|---|---|---|
| Near-surface D x50, 0.5 keV | 14.5 | 1.01 | 13.8 |
| D2+ current fraction 0.85, 0.25 keV | 3.08 | 3.31 | 1.00 |
| Half current, 0.5 keV | 1.00 | 1.00 | 2.00 |

Pass 5 (`ald/`): the merit is the charge-integrated D(d,p)T yield over the paper's run at each energy, divided by the same integral for the transported PdD reference.
- **Erosion decides it.** Sealing oxides and promoters sputter away within about 1 to 2 h at the UC Davis current. The 2 nm Al2O3 seal integrates to 9.8/1.2/1.8/3.8 at 0.25/1/2.5/6.5 keV (rising screening), against a peak of 128.
- **Best design: the 20 nm TiD2 cap**, with retained D and sputter x0.1, under rising screening:

| | 0.25 keV | 1 keV | 2.5 keV | 6.5 keV |
|---|---|---|---|---|
| Model | 210 | 18 | 37 | 38 |
| Scaled to the SI Table 1 TiD/PdD counts | 41 | 8.5 | 9.4 | - |

- **The cap's other endpoints.** At nominal sputter the scaled values are 12.5/1.6/7.0. A cap that loses its D gives 0.1 to 1.
- **Measured anchor.** The SI Table 1 TiD/PdD counts are 1.9 to 2.8 at 0.25-2 keV, 6.1 at 2.5 keV and 11.7 at 3 keV. Model TiD/PdD runs 1.2 to 8.4 times above them.
- **Sputter thresholds** for 10x are in `cap_sputter_thresholds.csv`. Scaled to the measured counts, 1 keV needs sputter about 30 to 40 times below nominal. Rising screening at 2.5 keV and neutrals at 0.25 keV stay below 10 even with no erosion.

Pass 5 (`transport/recomb/`): the D sheet on clean Pd(111) from first principles (PBE, relaxed, local harmonic corrections calibrated at 1 ML).
- **Binding.** The electronic binding is 0.489 eV/D at 1 ML and 0.518 eV/D at 0.5 ML. Against half a D2 at 1 Pa and 305 K, the binding is 0.117 eV/D.
- **No-beam sheet**, with capture assumed at 0.5:

| State | Sheet (ML) | D2 release (m^-2 s^-1) |
|---|---|---|
| Pd, equilibrated with 0.037 Pa D2 | 0.78 | 6.6e19 |
| PdD, fed from the bulk | 1.00 | 7.8e20 |

- **Clean recombination route.** The connected bridge-to-fcc route peaks 98 meV above clean surface plus gas, at the highest sampled frame. This is a restricted route, and the peak is not k- or slab-qualified. As a sensitivity, it lowers D2 release by a factor of 42 to 44 and leaves the sheet unchanged.
- **The transport model keeps its surface inputs.** The sheet the data prefer depends on beam stripping, which none of these calculations includes.

Pass 5 (`pairs/`): two candidates for keV screening that the static calculation leaves out.
- **Mobile deuterons.** A drift-diffusion estimate using the transported profiles gives at most 4e-5 eV in the top 100 nm of PdD and 0.11 eV over the whole foil. In BZY with the measured 3 mol% H, the static upper bound is 24 eV. All are far below the 100 eV that would matter.
- **D in a Pd vacancy (8-Pd cells, PBE, lowest enumerated basins).**
  - Binding per added D is 0.179, 0.230, -0.014 and -0.275 eV for 1, 2, 2 to 4 and 4 to 6 D. The last two are block means.
  - The closest D-D distance is 2.78 A, and the valence density at its midpoint is 0.028 e/A^3, against 0.295 in bulk PdD.
  - Through the merit consumer, the largest merit is 1.0003.
- **Ledger.** A 27-Pd size check of the binding was not reached. transport.py keeps Eb 0.23 eV and 6 D per vacancy. The capacity of 6 is untested; in the axis basin the binding levels off near 4 D.

Open transport item: `transport.py` lets superabundant vacancies reach 0.25 per Pd with 6 D each, so an infinite-time bulk-fed steady state holds 2.3 D/Pd. The 2 to 12 h run-end profiles used here are unaffected (vacancies move about 1 nm in 12 h). Superabundant-vacancy formation needs a kinetic limit before any long-time state is used.

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
| Pd80Rh20 (estimate from H-electrosorption data) | x0.31-0.33; yield x0.95-0.96 at 0.25 keV, x0.68-0.75 at 2.5 keV | - |
| -30 K / +30 K | x2.3-2.5 / x0.40-0.63 | ~x1 |
| Strain +0.5% / -0.5% | x1.4 / x0.68 | x1.00 |
| Back-side current x2 | x1.13-1.14 (end of run) | x1.00 |
| 2 nm Al2O3 while it lasts | x0.04 (0.25 keV) to x3.3 (PdD, 2.5 keV) | Al2O3 5 nm x0.19-0.30 |

A beam-side surface that almost stops recombining deuterium (sticking 1e-7) lets the hydride reach the beam side (0.68-0.72 D/Pd). First-principles inputs: the converged 1 ML adsorption of D on clean Pd(111) is 0.466 eV per D below half a D2 molecule, beside the model's 0.45 eV; the Pd vacancy binding stays at the 0.23 eV estimate (8-Pd pilot 0.131 eV, unrelaxed 27-Pd trial 0.242 eV); the D2 recombination barrier was not converged. Pd80Rh20 lowers near-surface D to about a third; at 0.25 keV the yield barely moves because the adsorbed monolayer stays nearly full, while at 2.5 keV it falls to x0.68-0.75. In the screening cells Rh moves the Pd d band away from the Fermi level (Pd3Rh -2.547 eV against Pd -2.424 eV); only tensile strain moves it toward it.

## What the models keep, drop, and leave out

`forward/` keeps the Bosch-Hale S-factors, ICRU-49 electronic and ZBL nuclear stopping with Bragg additivity, the 85/15 D+/D2+ beam, and a measurement layer: charge normalisation, Faraday-cup collected fraction in each beam-delivery regime, proton-to-neutron stitch factor, and fast neutrals. It drops target evolution, straggling and backscatter. It reads the deuterium depth profile when given one; `merit/` supplies the transported profiles.

`transport/` keeps the isotherm-driven diffusion over the full foil, vacancy trapping, surface desorption against chamber gas, sputter erosion, overlayers and beam heating. It drops microstructure, blisters and run history.

The three biggest omissions are:
1. Beam-side surface chemistry, which sets the adsorbed deuterium that carries 34-94% of the modelled Pd yield.
2. The real charge on the foil in each run.
3. Screening below about 50 eV.

## Open work

Open work is tracked on Fleet nodes, one directory per node when delivered:
- `g4_zenodo`: refit on the per-run data once it opens.
