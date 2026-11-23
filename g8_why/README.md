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

## Branch: what the beam delivers (`beam.md`, Marvin)

The beam is not only D+ at the set energy. Three ways it could differ were sized.

- Composition. The reported mix is about 85 percent D+ and 15 percent D2+ by current. A molecular ion splits its energy between its deuterons, so more molecules lower the yield. Pure D+ gives at most 1.18 times.
- Two deuterons from one D2+ arriving together. They land within a femtosecond and an angstrom of each other, but D2+ supplies only a few parts in 1e9 of the 0.25 keV yield, and the first deuteron's wake changes the second's by parts in 1e10.
- Fast neutrals. Deuterons that were neutralised while still at high energy would reach the foil above the set energy. A share of 3e-4 neutral atoms per ion at 13 keV multiplies the proton yield by 1e15 to 7e16 at 0.25 keV, 7e8 to 9e9 at 0.5 keV and 1e4 to 5e4 at 1 keV. Matching the measured 0.25 keV points takes 0.1 to 5 percent of 13 keV neutrals per ion (Pd 1.2, PdD 5.1, Ti 0.09, TiD 0.22 percent; about six times less under the other current reading). At 5 keV it takes 1.7 to 33 neutrals per ion.

Fast neutrals are the only mechanism in the four branches that reaches the size of the excess. Two things stand against them as described. The gun as documented accelerates directly to the set voltage up to 5 kV, with no stage where ions are faster and then slowed, so no source of 13 keV atoms at the 0.5 kV setting is known. And one fixed neutral share at 13 keV does not follow the four curves across energy (chi-squared 96 to 122 on 16 points). Chamber gas over the 4 cm flight could neutralise 0.3 to 3.8 percent of the beam, but those atoms keep the set energy and add nothing.

Observable: the energy spectrum of ions and neutrals at the foil plane against gun voltage, with the actual grid potentials. A constant upstream neutral rate would also give a yield per coulomb that falls as the current rises.

## Branch: the collision with a bound target (`collision.md`, giuliano)

The target deuteron is not at rest. It sits in a well, with zero-point and thermal motion, and its momentum could add to the beam's.

- Ordinary motion in the well (harmonic, 300 K): 1.035 times at 0.25 keV, 1.011 at 0.5 keV, 1.003 at 1 keV, nothing above.
- Any momentum distribution at all. Treat the target as a bound spectator that keeps its energy accounts: whatever momentum it has, taking it out of the well costs energy, so the collision energy cannot exceed the beam energy plus the excitation available. With up to 100 eV of excitation, the PdD proton yield at 0.25 keV relative to the reference cannot exceed 8e-9. The measured value is 0.05, six million times higher. This holds for every distribution, not only a Gaussian one.
- At 0.5 keV the same bound is 0.006 to 0.020 against a measured 0.017, so that point is not excluded under one of the two current readings.
- Measured momentum distributions of hydrogen in ZrH2 (neutron Compton scattering) are Gaussian with about 0.11 eV of kinetic energy, with no heavy tail.

A bound target does not explain the lowest point. The statement is made inside the spectator model; it is not a full many-body result. The computed well is strongly anharmonic, but the high-momentum tail of the real wavefunction and converged anharmonic yields were not computed, because the bound above makes them unnecessary.

Observable: the Doppler width of the fusion proton line. A target population at 100, 300 or 1000 eV would broaden the 3.02 MeV line to 3.4, 4.9 or 6.5 keV (one sigma, at 90 degrees).

The work also found that the older hydrogen setup at 0.20 A grid spacing changes the energy by 0.32 eV when the whole crystal is shifted rigidly by 0.1 A. No number in this node uses that well. The pair curves in `g6_collision`, `g7_defect` and `g7_hydride` are four-term differences at matched geometry, which cancel this term; single-well energies from that setup would not.

## The four branches together

With the deuterium profile fixed and 40 eV of screening, the measured curve is 1e16 to 1e17 above the model at 0.25 keV and nearly flat below 2.5 keV.

| Branch | Most it gives at 0.25 keV | Reaches the excess |
|---|---:|---|
| Reaction (threshold resonance) | 6 times | no |
| Bound target (any momentum distribution) | 8e-9 of the reference against 0.05 measured | no |
| Driven target (beam heating, 100 eV for 100 ps) | 1.25 times | no |
| Beam composition and paired arrival | 1.18 times | no |
| Fast neutrals at 13 keV, 0.1 to 5 percent per ion | 1e15 to 7e16 | yes, if present |

Three of the four branches change the yield by order one. None of them makes the bottom of the curve flat. The only route that reaches the size is deuterium arriving well above the set energy, and the documented gun has no source for it at low settings. The flat floor says the same thing from the data side: what reacts at the lowest settings is not meeting at the set collision energy.

The measurements that would separate these are an energy spectrum of ions and neutrals at the foil plane against gun voltage, the yield per coulomb against beam current at fixed energy, and the width of the proton line.

## Lead's note

`reaction.*`, `common.py` and `baseline.csv` are the delivered archive (SHA-256 c3653ee63b327c24eae7ec456251dc781ee8a28bb2264e659adf0cf3f0219962), unchanged. `driven_target.*`, `readout_floats.*` and `method.md` (delivered as its README) come from the second archive (SHA-256 ae4dd59696931dd33a0782fc2d97287599afa8777d611168a6edca188d187928); the lead read the page and checked its tables against the delivered CSVs but did not rerun it. `beam.*` is the third archive (SHA-256 5b6460beadbc919ab3c7799b39f700a2b540d586009000840485c692663f88c5), unchanged; the lead's rerun of `python g8_why/beam.py` reproduced the tables to 4 parts in 1e15 and the fitted tail shares to 1 part in 1e8. `method.md` is now the README of that third archive. `collision.md` and `collision/` are the fourth archive (SHA-256 2fc68b604d43b812243fbd901e4d01e4c193f34c6a3da1fde8b56e527d37ee66), unchanged; from `g8_why`, `python collision/reader.py .` replays the 32 reported bounds (the lead's run: passed, largest difference 2 parts in 1e15). Run `python -m g8_why.reaction` from the repository root with `g8_why` on the path (NumPy 2). The lead's rerun reproduced the tables to rounding. The excess table uses the yields as plotted; under the other reading of the current normalisation (`g7_gap/`) the 0.25 keV values are about six times smaller, which changes nothing above.
