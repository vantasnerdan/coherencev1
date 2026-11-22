# g8_why: why and how a plain foil gives more D-D fusion than the model

The results to explain are plain metal foils under a low-energy deuteron beam: the D-D yields for Pd, PdD, Ti and TiD of Karahadian et al. (Nat. Commun. 17, 8845 (2026)). Earlier nodes kept explaining the excess by where the deuterium sits, because the deuterium profile was the only free quantity left in the forward model. This node holds the profile fixed at the transported one from `g4_beam/merit` and works the other side: the reaction, the collision, the target while the beam is on, and what the beam delivers.

Rules: the profile is not a variable; no site fractions, layers or site counts; each branch ends in a mechanism, its effect on the yield as a ratio against the 40 eV model, and one observable other than total yield; no target design.

## How large the excess is

With the profile fixed and 40 eV of screening, the measured yield, relative to each foil's own value at 6.5 keV centre of mass, exceeds the model by:

| Foil | 0.25 keV | 0.5 keV | 1 keV | 2.5 keV |
|---|---:|---:|---:|---:|
| Pd | 5e16 | 4e9 | 5e4 | 3.0 |
| PdD | 2e17 | 2e10 | 2e5 | 12.7 |
| Ti | 2e17 | 2e10 | 8e4 | 2.1 |
| TiD | 5e17 | 2e10 | 8e4 | 6.8 |

The measured proton yield at 0.25 keV is 0.2 to 0.65 of its value at 2.5 keV. Any collision at the nominal energy through a Coulomb barrier falls by thirteen orders of magnitude or more over that range. The measured curve is close to flat at the bottom, so the excess is not a correction to barrier penetration at the nominal energy. Whatever reacts at the lowest beam settings is not meeting at the nominal collision energy.

## Branch: the reaction (`reaction.md`, Marvin)

A threshold resonance in 4He (Czerski, Phys. Rev. C 106, L011601 (2022)) changes the reaction strength itself. With the published strength it multiplies the proton yield by 5.9-6.4 at 0.25 keV, 1.9-2.1 at 0.5 keV and 1.1 at 1 keV, and lowers it by 4-6 percent at 2.5 and 6.5 keV. That is against a required 1e16 or more at 0.25 keV. For any strength and any phase of that term, the proton yield at 0.25 keV over 2.5 keV stays below 5e-14; the data give 0.2-0.65. Changing the resonance per foil does not close it. Host effects on the exit channels are of order 1e-5.

The reaction does not explain the excess. Observable that would show a reaction change: the proton, neutron and electron-pair branching against energy and foil.

## Branch: the target while the beam is on (`driven_target.md`, Marvin)

If the barrier stays ordinary, the flat bottom of the curve needs deuteron pairs meeting at keV-scale relative energy. One beam ion at 0.5 keV cannot give more than 0.5 keV to any pair, so a keV tail would have to pool energy from several beam ions or from another reservoir.

- Ordinary thermal and zero-point motion: 1.04 times at 0.25 keV, less above.
- Every watt the beam deposits given to the deuterons as a 100 eV population that lives 100 ps: 1.25 times at 0.25 keV in Pd and PdD, 1.015 in Ti and TiD, nothing above.
- A 1 keV population large enough to give the measured floor would have to hold its energy for 40 to 1600 seconds at the full beam power. Measured and computed vibration lifetimes in PdH are under a picosecond; electron-to-lattice transfer takes picoseconds.
- Two beam ions overlap within 10 nm and 100 ps with probability about 1e-7.

Heating and beam-driven excitation of the target do not explain the excess. Observable that would show pooling: a floor that pooled energy from two ions would grow with beam current per unit charge; a single-ion process gives a constant yield per coulomb.

`method.md` is Marvin's note on the shared baseline and on what the published numbers mean. It records the measured low-energy signal as detected counts per coulomb of cup current: between 0.25 and 2.5 keV these stay within a factor of 2 to 8 of each other in every foil (Pd 1.1, 0.4, 0.8, 1.7; TiD 2.9, 1.8, 2.0, 15.2).

## Branches to come

`collision` (giuliano) and `beam` (Marvin).

## Lead's note

`reaction.*`, `common.py` and `baseline.csv` are the delivered archive (SHA-256 c3653ee63b327c24eae7ec456251dc781ee8a28bb2264e659adf0cf3f0219962), unchanged. `driven_target.*`, `readout_floats.*` and `method.md` (delivered as its README) come from the second archive (SHA-256 ae4dd59696931dd33a0782fc2d97287599afa8777d611168a6edca188d187928); the lead read the page and checked its tables against the delivered CSVs but did not rerun it. Run `python -m g8_why.reaction` from the repository root with `g8_why` on the path (NumPy 2). The lead's rerun reproduced the tables to rounding. The excess table uses the yields as plotted; under the other reading of the current normalisation (`g7_gap/`) the 0.25 keV values are about six times smaller, which changes nothing above.
