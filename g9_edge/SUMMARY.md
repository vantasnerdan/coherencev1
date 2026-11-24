# g9_edge: the instant a deuteron arrives (lead's summary)

Every earlier model in this repository treats the arrival of a beam deuteron as settled at each instant: one static screening energy, a tabulated drag, a target at rest. This node builds the part that was lumped away, the change from no ion to ion at the face of the foil. It reports peaks and durations, not lifetimes. The deuterium profile is fixed at the transported one in `g4_beam/merit`, as in `g8_why`.

Two parts. The model part (Marvin) is `README.md` and `replay.md` in this directory, with its script and tables. The first-principles part (giuliano, a deuteron entering a PdD slab from vacuum, electrons followed in time) is in `first_principles/` and is summarised below.

## Model part: what it is

A charged deuteron is absent, then appears at the surface plane and moves in at the beam speed. The electrons respond through the dielectric response of palladium computed in `g6_response`, continued with a free-electron form where that calculation runs out, plus a simple image term for the surface. The same calculation is run for an ion that appears inside the bulk, so the part that belongs to the surface can be read off. The settled value used for comparison is the same model's own static limit (25 eV at 0.1 A, where the first-principles cells give about 40 eV).

It keeps the time dependence of the electron response and the force on the neighbouring deuterons. It drops the ion's approach through the vacuum, its pick-up of an electron, any non-linear response, and motion of the nuclei during the edge. Palladium's response stands in for all four foils.

## What it gives

**The screening is absent at first, then overshoots.** At 0.25 keV centre of mass and 0.1 A ahead of the ion:

| | Departure from the settled value | When | For how long |
|---|---:|---:|---:|
| At switch-on | 25 eV less screening (none has formed) | 0 | 0.07 fs |
| Overshoot, entry at the surface | 17 eV more screening | 0.11 fs | 0.07 fs |
| Overshoot, same ion in the bulk | 2.9 eV more screening | 0.05 fs | 0.04 fs |

The overshoot at the surface is about six times the bulk one. The transient is the same size as the settled number it replaces, and it is over in about 0.15 fs, during which the ion has moved 0.2 to 0.5 A at the four lower energies.

**Most of the low-energy yield of the model sits inside that edge.** In the fixed profile, the share of the model's reactions that happen within the edge (the adsorbed sheet at the face and the first 0.2 to 0.5 A behind it) is:

| Foil | 0.25 keV | 0.5 keV | 1 keV | 2.5 keV | 6.5 keV |
|---|---:|---:|---:|---:|---:|
| Pd | 97% | 92% | 86% | 73% | 64% |
| PdD | 96% | 91% | 83% | 66% | 44% |
| Ti, TiD | 34% | 13% | 5% | 2% | 1% |

So for palladium the settled screening number has been applied almost entirely in the one place where it has not had time to form. This is the main thing the node has shown so far, and it does not depend on the details of the response model.

**What that does to the yield depends on the ion's history before the surface.** With the ion appearing unscreened at the plane, the proton yield against the 40 eV model is 0.28 (Pd), 0.28 (PdD), 0.76 (Ti, TiD) at 0.25 keV, rising to 0.99 or more at 6.5 keV: the sheet reacts before any screening exists. If instead the screening at the sheet is taken as already formed, as for an ion that has been travelling, the Pd ratios are 1.017, 1.011, 1.005, 1.002, 1.001. The real case lies in the approach through the vacuum, where the surface starts to respond before the ion arrives. The model does not contain that approach; the first-principles part does.

**Speed.** One angstrom of travel takes 3.6 electron response periods at 0.25 keV and 0.7 at 6.5 keV. On this measure the electrons keep up better at the low energies than at the high ones.

**The kick to the neighbouring deuterons** from the unsettled field is small: 31 micro-eV to the target deuteron and 0.37 meV summed over it and its twelve neighbours at 0.25 keV (0.16 to 0.47 meV across the five energies).

**Stores the arrival could interrupt** (`stored_energy.csv`). The one with a number from this repository is a region of hydrogens moving in step: 0.6 meV per deuteron above the ground state in the pumped PdD lattice run, 19.7 eV in the whole 32,768-deuteron box. The strain of implanted deuterium is about 1.8 meV per host atom at 1 percent strain. For the screening clouds, the surface dipole and the bias field, the page gives no excited energy: at rest these are ground-state quantities, and what an arrival can release from them is not known. How much of any store reaches how many pairs is left open in every row.

Lead's arithmetic on the first row: to put 1 keV on one pair from a region pumped at that level, the motion of about 1.7 million deuterons would have to arrive on that pair.

**Observable other than yield:** the recoil of the target shifts the 3.02 MeV proton line by -15 to +17 eV along the beam axis during the edge at 0.25 keV.

## First-principles part: one ion entering a PdD slab

A bare deuteron at 0.5 keV starts 5 A above a clean, frozen PdD(111) slab and moves along the surface normal toward the first deuteron below the surface. The electrons are followed in time, and the same geometries are computed in the settled (static) way for comparison. The run stops 0.5 A short of the target.

**The electron cloud forms just before the surface.** The ion stays bare until about 1.5 A out. It reaches half of the settled charge at 1.4 to 1.6 A above the outer plane, picking it up at 1.5 to 2.5 electrons per femtosecond at the fastest. It is still 0.19 to 0.28 of an electron short where the run ends. So the real arrival lies between the two cases of the model part: not switched on unscreened at the plane, and not fully formed.

**What the target deuteron sees, as the departure from the settled picture:**

| | Peak | Width | Steepest rate |
|---|---:|---:|---:|
| Potential energy | 5.2 eV | 0.45 fs | 19 eV per fs |
| Force, outward lobe | 2.7 eV/A | 0.21 fs | 66 eV/A per fs |
| Force, inward lobe | 4.7 eV/A | cut off by the end of the run | |

**What it hands to the target nucleus is small.** Integrating the force on a target held in place gives 0.014 eV of recoil, of which 0.002 eV is the excess over the settled picture.

The transient has the shape the model part found, a departure followed by a reversal within a fraction of a femtosecond, and is a few eV in size where the model gave 17 to 25 eV. The time step (4 attoseconds) was not checked against a finer one, by decision: the result was judged enough to place the effect. One ray, a clean face with no adsorbed deuterium, and nuclei held fixed.

## Three largest omissions

1. The approach through the vacuum, electron pick-up and non-linear response at the surface. The yield result turns on this.
2. How much excited energy the foil holds where the ion lands, and how it would be passed to a pair.
3. Nuclei that move, and tunnelling that changes, during the edge.

## Lead's note

The archive is SHA-256 d4da0afd626eaee6e941a073ffc4eeae9045473e42def02e504ca0577468eebd, unchanged; `SHA256SUMS` checks. The script needs the palladium response matrices from the `g6_response` delivery (SHA-256 1be3de744451f7938d44bd3dced26854849643e80e8c84f1f35e77991a384c10), which are not in the repository, and NumPy 2. The lead's rerun, `python g9_edge/edge.py --matrix-dir <g6_response/data/Pd> --output <new dir>` from the repository root, reproduced every table (largest difference 1 part in 1e10). The image-term cutoff changes one near-ion curve by 14 eV and the yield ratio by 0.006; the surface treatment is the least settled part of the model.

First-principles archive: SHA-256 a84c480f784a2edd5279175fd3dff3e05e0c8617db6d5818f52dd2b054c426fe, unchanged; its file sums check. The lead reran `readout/edge_readout.py` on the raw journals in `first_principles/native/`: the potential, force and charge tables are byte-identical and the recoil table agrees to 3 parts in 1e10.
