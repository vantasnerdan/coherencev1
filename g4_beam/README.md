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
| `forward/` | Thick-target yield model with a measurement layer, fitted to Fig 2 (about 6 CPU-minutes) |
| `transport/` | 1D deuterium transport with beam damage, binary-collision ranges, knob sweeps |
| `screening/` | 184 measured screening potentials with sources, screening-model formulas, correlation tests |
| `priordata/` | Other low-energy d+d datasets in one table, plateau comparison |

Run every script from its own directory, except `work/` scripts and `priordata/scripts/compare.py .`, which run from `g4_beam/`.

## Current best numbers

Bulk screening potential U_e (keV) from a thick-target fit to the 2.5-6.5 keV points. The value depends on how the beam charge is normalised. Fig 2 and SI Table 1 differ by one factor of cup current, so the model carries the normalisation as a hypothesis:

| Charge normalisation | chi2 (28 pts) | Pd | PdD | Ti | TiD |
|---|---|---|---|---|---|
| As plotted | 33.9 | 1.27 | 2.20 | 0.99 | 1.34 |
| One power of cup current | 19.9 | 0.70 | 1.41 | 0.45 | 0.84 |
| Two powers | 25.5 | 0.25 | 0.80 | 0.03 | 0.45 |

The hydride-minus-metal step is stable across these: PdD-Pd 0.55-0.93 keV, TiD-Ti 0.35-0.42 keV.

Low-energy plateau. These three candidates fit all 48 points at chi2/dof ~ 1 and lie within delta-BIC 2.2 of each other:
- a sub-monolayer surface deuterium population with U_e 6-12 keV;
- the paper's sigmoid U_e placed inside the thick-target integral, with run-to-run values 1.6-2.2 keV;
- a fast-neutral beam component of at least 0.18 of the cup current at about 5 keV lab.

These do not fit: no plateau (chi2 459), a threshold resonance with one width (chi2 >= 136), and any deuterium layer thicker than 0.2 nm.

Near-surface deuterium in the top 5 nm, against the base case, from `transport/`:

| Knob | Pd / PdD | Ti / TiD |
|---|---|---|
| Pd77Ag23 / Pd90Pt10 | x6.5-12 / x0.27-0.31 | - |
| -30 K / +30 K | x2.3-2.5 / x0.40-0.63 | ~x1 |
| Strain +0.5% / -0.5% | x1.4 / x0.68 | x1.00 |
| Back-side current x2 | x1.13-1.25 | x1.00 |
| 2 nm Al2O3, while it lasts (0.9-1.7 h) | x0.04 (0.25 keV) to x3.3 (PdD, 2.5 keV) | Al2O3 5 nm x0.19-0.30 |

Beam-side surface recombination is the largest lever in the transport model. At sticking <= 1e-6, the beta phase reaches the beam surface (D/Pd 0.6-0.7, against ~0.01 for a clean surface).

## What the models keep, drop, and leave out

`forward/` keeps the Bosch-Hale S-factors, ICRU-49 electronic and ZBL nuclear stopping with Bragg additivity, the 85/15 D+/D2+ beam, and a measurement layer: charge normalisation, Faraday-cup collected fraction in each beam-delivery regime, proton-to-neutron stitch factor, and fast neutrals. It drops target evolution, straggling and backscatter, and the deuterium depth profile.

`transport/` keeps the isotherm-driven diffusion over the full foil, vacancy trapping, surface desorption against chamber gas, sputter erosion, overlayers and beam heating. It drops microstructure, blisters and run history.

The three biggest omissions are:
1. Beam-side surface chemistry, which sets the adsorbed deuterium that carries 34-94% of the modelled Pd yield.
2. The real charge on the foil in each run.
3. Screening below about 50 eV.

## Open work

Open work is tracked on Fleet nodes, one directory per node when delivered:
- `g4_screen`: eps(q), N(E_F), d-band centre and the adiabatic screening limit from GPAW. This gives ratio 7 for the map and the U_e response to d-band shifts.
- `g4_forward`: second-route check of the forward model, and a fit to the Czerski et al. data (arXiv:2409.02112).
- `g4_transport`: second-route check of the transport model, and DFT inputs for the Pd surface and vacancy trapping.
- `g4_zenodo`: refit on the per-run data once it opens.
