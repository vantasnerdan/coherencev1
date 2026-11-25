# g13_rogue: rogue waves in the beam-driven lattice (lead's summary)

The idea: a foil under a continuous beam is never quiet. Each impact launches waves that travel through the lattice, waves from different impacts overlap, and in a medium whose vibration frequency depends on amplitude one vibration can grow at the expense of its neighbours. That is modulational instability, the mechanism behind ocean rogue waves, and its result is a small hot spot with far more energy per site than its surroundings. `g11_cmns/collective/focusing` found the tallest peak that can be made on purpose from a quiet start. This node asks what a lattice that is driven all the time does without being arranged.

Pass 1 (Marvin) is the coarse model: one envelope equation for the deuterium vibration in PdD, with every coefficient taken from wells and couplings already in the repository. `README.md` in this directory is the full technical page, with equations, inputs and every table. An addendum is commissioned (see the last section) and will supersede parts of this.

## What the model is

One complex amplitude per deuteron site for the vibration along the nearest-neighbour direction. It keeps the real well's first two transition energies, the signed couplings to neighbours, a loss time, random impacts at the beam current, waves arriving from other impact positions, and heat at 300 K. It is run twice, once with the amplitude dependence of the frequency and once without, on the same random numbers, so the difference is what the non-linearity adds.

It drops any new electronic calculation, the electron sea beyond what is already inside the wells, three-dimensional hot spots (the grown one is uniform across the other two directions), and a measured figure for how much impact energy ends in the deuterium vibration. That last figure is carried as two assumed values, 0.1 percent and 10 percent.

## What it gives

**The instability exists, and its threshold is low.** The deuterium vibration hardens with amplitude: the first step is 42.6 meV and the second 48.0 meV, a shift of 5.3 meV per quantum. A hardening vibration can only feed a ripple near the top of its band, where the band curves downward. There a ripple grows faster than it is lost once the vibration carries:

| Loss time of the amplitude | Threshold energy per deuteron | Ripple wavelength that grows first |
|---|---:|---:|
| 8.2 ps | 0.64 meV | 11.3 nm |
| 4.1 ps | 1.28 meV | 7.8 nm |
| 2.05 ps | 2.57 meV | 5.4 nm |

0.64 meV is 0.015 of a vibration quantum per site, a displacement of 0.038 A. This is for motion in step (one coherent wave). Partial loading moves the threshold to 0.69 to 0.84 meV, but those rows rest on wells probed only to 0.03 A and are comparisons, not predictions.

**What a ripple grows into.** A wave prepared at 5.3 meV per deuteron with a 5 percent ripple grows, in 7.5 ps, into a spot holding 69.7 meV per deuteron at its peak: thirteen times the starting level, 2.6 nm wide along the wave, lasting 0.97 ps at half height. A neighbouring pair in that spot closes with 0.086 eV and reaches 2.35 A at its nearest (the bulk spacing is 2.89 A). The wave was put there and held by an explicit pump; the model does not say the beam prepares it.

**What the beam maintains is far below the threshold.** Standing level of the deuterium vibration kept up by the beam, with 10 percent of each recoil channel converted at each step:

| Beam energy | Beam current | From deuterium recoils | From lattice waves arriving |
|---|---:|---:|---:|
| 0.25 keV | 38 uA | 5e-17 eV per deuteron | 6e-15 eV per deuteron |
| 2.5 keV | 138 uA | 1e-15 | 3e-14 |
| 6.5 keV (scaled from 2.5) | 140 uA | 3e-15 | 8e-14 |

The threshold is 6e-4 eV per deuteron, ten to eleven powers of ten higher. Lead's arithmetic on the same numbers: if every electronvolt of the 0.25 keV beam went into the deuterium vibration of the layer the beam stops in, the level would be 3.5e-12 eV per deuteron, still short by a factor of about 2e8. Put the other way, the beam supplies about 0.9 eV per deuteron per second to that layer, so reaching the threshold needs the energy to be held for about a millisecond, where the vibration holds it for picoseconds.

**Heat is already there at a higher level.** At 300 K the same vibration holds 10.2 meV per deuteron, sixteen times the threshold, but as incoherent motion, to which the threshold does not apply directly.

**How many impacts overlap** at one point, counting an impact only once its wave has had time to arrive, at 150 uA over 1 cm2:

| Distance a wave survives | Impacts present at once |
|---|---:|
| 3 nm | 7e-11 |
| 100 nm | 2e-6 |
| 5 um | 0.3 |
| 10 um | 2.5 |

So long-range waves do overlap continuously, as the idea requires. The count says how many impacts contribute, not how much energy they bring; the energy is the table above.

**Peak statistics.** In windows around arriving wave packets, with heat present, the tallest energies seen are 0.12 to 0.15 eV on a site and 0.11 to 0.23 eV closing a pair. These peaks are thermal. The packets from the beam add nothing that can be told from zero in twelve random samples (at 2.5 keV the excess rate of closings above 0.093 eV is -44 +/- 29 per deuteron per second in the delivered run and -0.02 +/- 0.02 in the lead's rerun). In both runs the non-linear model's tallest pair closing is about 1.6 times the linear one (0.18 against 0.11 eV delivered, 0.23 against 0.14 eV rerun); with twelve samples that is an observation, not a measured tail.

**As a reaction chance.** The grown hot spot's closing episode, passed unchanged to the contact calculation of `g11_cmns`, gives a proton chance of 1e-114 to 1e-110 per pair per pulse.

**Scaling with beam current** at fixed total charge: the contribution of single impacts does not depend on current; a contribution from two overlapping impacts rises in proportion to current. The model's supplied level rises with beam energy, so as built it is not flat in energy.

**The electron sea term could not be filled in.** It needs two numbers the repository does not hold: how much the deuterium vibration frequency changes per change in local electron density, and how large an electron density ripple a lattice strain wave carries. The palladium response tables do not give either. The term is left empty, not set to zero.

## Where the beam stops there is almost no deuterium (lead)

The stopping record this pass reads (`g11_cmns/collective/part1/stopping.json`) gives 0.06 eV of a 435 eV deuteron going directly to lattice deuterons and 79 eV to palladium. The reason is in the fixed profile (`g4_beam/merit/profiles.json`), not in the stopping law. In the PdD foil of that profile the ratio of deuterium to palladium is:

| Depth | D per Pd |
|---|---:|
| first 0.25 nm | 0.012 |
| 2 nm to 150 nm | about 1e-4 |
| 20 um | 0.008 |
| 31 um | 0.07 |
| 51 um | 0.43 |
| 200 um and deeper | 0.61 |

The beam stops within 24 nm at 0.25 keV and 129 nm at 2.5 keV. So in this profile the impacts happen in nearly deuterium-free palladium, and the loaded lattice that the envelope describes begins tens of micrometres deeper. This pass placed the fully loaded lattice at the stopping depth, which the profile does not contain. The only link between the two regions is a lattice wave that travels 30 um or more, and the transport table here stops at 5 um.

## Commissioned next (addendum to this pass)

1. State where each part of the model sits in the profile.
2. The standing level in the loaded region from waves launched where the beam stops, over 30, 50 and 200 um, with the attenuation length as a named sensitivity up to the foil thickness and reflection from the foil faces as its own row.
3. What the instability means among isolated deuterons at 1e-4 to 1e-2 per palladium atom, where there is no deuterium band.
4. Random realisations that reproduce across machines, and enough of them that the non-linear against linear comparison has a usable error.

## Three largest omissions

1. A measured conversion from each recoil channel into lattice waves and deuterium vibration, with real attenuation and foil boundaries.
2. A non-linear model beyond the first two transitions and one polarisation, calibrated at partial loading.
3. The electron density term, and the quantum state of the pair at contact.

## Corrections this pass made to earlier pages

- The loss time of 4.1 ps used in `g11_cmns/collective/focusing` is described there as an amplitude decay time. The engine applies it to the velocity, so 4.1 ps is the lifetime of the energy and the amplitude lives 8.2 ps. No focusing number changes; only the label.
- In `g11_cmns/collective/part1/stopping.json` the field `terminal_depth_A` holds metres, not angstroms. The `depth_A` values in its timeline are in angstroms as named.

## Lead's note

Archive SHA-256 c82650415c8aa2188251cc918c89aec0960cb3b805e0eb9eebc78b95404a2ea8; all 32 file sums check. The lead reran `run.py` from scratch and then `readout.py` (NumPy, SciPy, Matplotlib) on a clean copy of the repository at commit 02fe7cc plus this directory. The transitions, thresholds, dispersion, loading, drive, transport and overlap tables, the grown hot spot and the contact numbers reproduce (largest difference 6 parts in 1e6, in quantities that are near zero; the contact range agrees to 3 parts in 1e11). The random peak statistics do not reproduce number for number: the generator is seeded, but the thermal noise is built on eigenvectors whose choice differs between machines, so the rerun is a different draw from the same ensemble. Both draws are quoted above. The committed outputs are the delivered ones. The count of impacts in the overlap table is Marvin's causal count; the lead's first estimate in the node used a different convention and was three times larger.
