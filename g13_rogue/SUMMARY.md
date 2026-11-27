# g13_rogue: rogue waves in the beam-driven lattice (lead's summary)

The idea: a foil under a continuous beam is never quiet. Each impact launches waves that travel through the lattice, waves from different impacts overlap, and in a medium whose vibration frequency depends on amplitude one vibration can grow at the expense of its neighbours. That is modulational instability, the mechanism behind ocean rogue waves, and its result is a small hot spot with far more energy per site than its surroundings. `g11_cmns/collective/focusing` found the tallest peak that can be made on purpose from a quiet start. This node asks what a lattice that is driven all the time does without being arranged.

Pass 2 is in `pass2/` (lead's summary: `pass2/LEAD.md`). It models the whole foil as a vibrating body under a live, uneven beam and counts impacts into moving places, converging waves, and rogue peaks and holes. Pass 4 (`pass4/SUMMARY.md`) carries the same whole foil up to the fast, short waves where an impact's energy starts. Pass 3 (`pass3/`, first step `pass3/rung1/LEAD.md`) refines that same whole-foil model toward shorter waves in steps. It replaces this page's picture of identical small packets launched into a quiet lattice; the instability and deuterium-vibration results below stand.

Pass 1 (Marvin) is the coarse model: one envelope equation for the deuterium vibration in PdD, with every coefficient taken from wells and couplings already in the repository. `README.md` in this directory is the full technical page, with equations, inputs and every table.

The pass has two layers, and they must be read apart:

- **The loaded-lattice model** (first delivery): a fully loaded PdD lattice, its instability, a grown hot spot, and peak statistics. It was placed at the depth where the beam stops.
- **The foil as the fixed profile has it** (addendum): where the beam stops there is almost no deuterium, and the loaded region begins tens of micrometres deeper. The addendum carries the beam's lattice waves down to that region.

## What the model is

One complex amplitude per deuteron site for the vibration along the nearest-neighbour direction. It keeps the real well's first two transition energies, the signed couplings to neighbours, a loss time, random impacts at the beam current, waves arriving from other impact positions, and heat at 300 K. It is run twice, once with the amplitude dependence of the frequency and once without, on the same random numbers, so the difference is what the non-linearity adds.

It drops any new electronic calculation, the electron sea beyond what is already inside the wells, three-dimensional hot spots (the grown one is uniform across the other two directions), and a measured figure for how much impact energy ends in the deuterium vibration. That last figure is carried as assumed fractions (0.1 percent and 10 percent at each step).

## The loaded lattice: the instability and what it grows into

**The instability exists, and its threshold is low.** The deuterium vibration hardens with amplitude: the first step is 42.6 meV and the second 48.0 meV, a shift of 5.3 meV per quantum. A hardening vibration can only feed a ripple near the top of its band, where the band curves downward. There a ripple grows faster than it is lost once the vibration carries:

| Loss time of the amplitude | Threshold energy per deuteron | Ripple wavelength that grows first |
|---|---:|---:|
| 8.2 ps | 0.64 meV | 11.3 nm |
| 4.1 ps | 1.28 meV | 7.8 nm |
| 2.05 ps | 2.57 meV | 5.4 nm |

0.64 meV is 0.015 of a vibration quantum per site, a displacement of 0.038 A. This is for motion in step (one coherent wave) at full loading. Partial loading moves the threshold to 0.69 to 0.84 meV, but those rows rest on wells probed only to 0.03 A and are comparisons, not predictions.

**What a ripple grows into.** A wave prepared at 5.3 meV per deuteron with a 5 percent ripple grows, in 7.5 ps, into a spot holding 69.7 meV per deuteron at its peak: thirteen times the starting level, 2.6 nm wide along the wave, lasting 0.97 ps at half height. A neighbouring pair in that spot closes with 0.086 eV and reaches 2.35 A at its nearest (the bulk spacing is 2.89 A). The wave was put there and held by an explicit pump; the model does not say the beam prepares it.

**As a reaction chance.** That closing episode, passed unchanged to the contact calculation of `g11_cmns`, gives a proton chance of 1e-114 to 1e-110 per pair per pulse.

**Heat is already there at a higher level.** At 300 K the same vibration holds 10.2 meV per deuteron, sixteen times the threshold, but as incoherent motion, to which the threshold does not apply directly.

## The foil as the fixed profile has it

**Where the deuterium is.** In the PdD foil of `g4_beam/merit/profiles.json` (250 um thick), deuterium per palladium atom is:

| Depth | 0.25 keV record | 2.5 keV record |
|---|---:|---:|
| Along the beam's track (averaged) | 0.0008 | 0.012 |
| 30 um | 0.012 | 0.011 |
| 50 um | 0.35 | 0.32 |
| 200 um | 0.61 | 0.61 |

The beam stops within 24 nm at 0.25 keV and 129 nm at 2.5 keV. Along that track the deuterium fraction ranges from 4e-5 to 0.02. Nowhere does the foil reach the full loading the instability table above was computed for; the deepest layers hold 0.61. In the Pd foil the fraction is 3e-5 at all these depths.

This is why the stopping record gives 0.06 eV of a 435 eV deuteron going directly to lattice deuterons and 79 eV to palladium: the stopping law treats deuterium and palladium alike (per atom, deuterium takes 0.6 of what palladium takes at 500 eV), and there is simply little deuterium on the track.

**Among isolated deuterons there is no band to go unstable.** On the track each deuteron vibrates in its own palladium cage with, on average, far less than one deuterium neighbour. The instability above needs a connected deuterium vibration band carrying a wave over many sites. Rare neighbouring pairs exist, but they do not form that wave, so the 0.64 meV threshold has no meaning there. A coupling between distant deuterons through the palladium would be a different mechanism, not built here.

**Lattice waves reaching the loaded region.** The palladium recoils on the track (91 and 631 eV per unit charge at 0.25 and 2.5 keV) launch waves downward at 4 km/s. With 10 percent of that energy launched as travelling waves and an attenuation length of 250 um:

| Beam energy | Depth | Impacts whose waves are present at once | Wave energy per local deuteron | Of which converted to deuterium vibration |
|---|---:|---:|---:|---:|
| 0.25 keV | 30 um | 66,000 | 1e-11 eV | 7e-17 eV |
| 0.25 keV | 50 um | 66,000 | 3e-13 | 2e-18 |
| 0.25 keV | 200 um | 54,000 | 4e-14 | 3e-19 |
| 2.5 keV | 30 um | 240,000 | 3e-10 | 2e-15 |
| 2.5 keV | 50 um | 240,000 | 8e-12 | 5e-17 |
| 2.5 keV | 200 um | 200,000 | 1e-12 | 7e-18 |

So the overlap the idea asks for is there in full: with waves that cross the foil, every point of the loaded region is under the waves of tens to hundreds of thousands of impacts at any instant. With a 30 um attenuation length the count at 50 um falls to 58 (0.25 keV) and 210 (2.5 keV), and with 5 um to below one.

What they bring is small. The wave energy per deuteron is at most 3e-10 eV, and the part handed to the deuterium vibration (a further 10 percent of what is absorbed locally, held for the 4.1 ps energy lifetime) is 2e-15 eV or less. Against the full-loading threshold of 6e-4 eV per deuteron that is eleven to fifteen powers of ten short. One reflection from the back face of the foil (90 percent reflected) adds 3 to 46 percent to the wave energy, most at 200 um.

The same arithmetic for the first delivery's scenario, where a loaded lattice sits on the track: 5e-17 to 8e-14 eV per deuteron maintained, and 3.5e-12 eV if every electronvolt of the 0.25 keV beam went into the vibration (lead's figure). The beam supplies about 0.9 eV per deuteron per second to that layer, so reaching the threshold would need the energy held for about a millisecond, where the vibration holds it for picoseconds.

## Peak statistics

These are for the loaded lattice with heat at 300 K and beam wave packets arriving, in windows around each arrival. The addendum redrew the random numbers in a form that does not depend on the machine, with 64 samples for every case and 256 for the 2.5 keV case.

- **The peaks are thermal.** The tallest closing energies on a pair are 0.15 eV (linear model) and 0.24 eV (non-linear) in 256 samples; the typical ones are 0.09 and 0.12 eV.
- **The beam's packets add nothing that can be told from zero.** The excess rate of pair closings above 0.093 eV caused by an arriving packet, non-linear minus linear, is +76 +/- 140 per deuteron per second, with a 95 percent interval of -900 to +1050. For scale, heat alone gives of order 1e8 to 1e9 such window maxima per deuteron per second.
- **The non-linear model's taller peaks are not concentration.** The non-linear tallest closing is 1.67 times the linear one, and 1.31 times in the typical sample. Taking the peak apart shows where it comes from: the hardened vibration runs at a higher frequency, so the same displacement means a higher velocity. Removing only that frequency term leaves a difference of -2.3 +/- 1.3 meV, that is, none. The energy gathered on a site is, if anything, slightly lower in the non-linear model (0.091 against 0.102 eV on a common measure). So in this ensemble the non-linearity does not gather energy into rare tall spots.

## Scaling

At fixed total charge, the contribution of single impacts does not depend on beam current, and a contribution from overlapping impacts rises in proportion to current. Every level in the depth table is proportional to current. The supplied level rises with beam energy, so as built the model is not flat in energy.

## The electron sea term could not be filled in

It needs two numbers the repository does not hold: how much the deuterium vibration frequency changes per change in local electron density, and how large an electron density ripple a lattice strain wave carries. The palladium response tables do not give either. The term is left empty, not set to zero.

## What this pass shows

- The deuterium vibration of loaded PdD is modulationally unstable at a low coherent level, and the instability does build a hot spot when the wave is prepared.
- The beam's waves do reach and overlap in the loaded region, in very large numbers if they cross the foil.
- The energy those waves carry to each deuteron is ten or more powers of ten below the level at which the instability starts, and the random thermal motion that is far above that level does not concentrate in the model.
- What would change this is something that holds vibration energy for milliseconds instead of picoseconds, a coherent source much stronger than the beam's recoils, or a term the model lacks. The electron density term is the one missing piece named by the idea itself.

## Three largest omissions

1. A measured conversion from each recoil channel into lattice waves and deuterium vibration, with real attenuation and all foil boundaries (only one back-face reflection is included).
2. A non-linear model of the partly loaded, disordered lattice (deuterium fraction 0.3 to 0.6) that the foil actually contains at depth; everything non-linear here is calibrated at full loading, two transitions and one polarisation.
3. The electron density term, and the quantum state of the pair at contact.

## Corrections this pass made to earlier pages

- The loss time of 4.1 ps used in `g11_cmns/collective/focusing` is described there as an amplitude decay time. The engine applies it to the velocity, so 4.1 ps is the lifetime of the energy and the amplitude lives 8.2 ps. No focusing number changes; only the label.
- In `g11_cmns/collective/part1/stopping.json` the field `terminal_depth_A` holds metres, not angstroms. The `depth_A` values in its timeline are in angstroms as named.

## Lead's note

First delivery: archive SHA-256 c82650415c8aa2188251cc918c89aec0960cb3b805e0eb9eebc78b95404a2ea8, all 32 file sums checked. The lead reran `run.py` from scratch and then `readout.py` (NumPy, SciPy, Matplotlib) on a clean copy of the repository at commit 02fe7cc plus this directory. The transitions, thresholds, dispersion, loading, drive, transport and overlap tables, the grown hot spot and the contact numbers reproduce (largest difference 6 parts in 1e6, in quantities that are near zero; the contact range agrees to 3 parts in 1e11). Its twelve-sample random statistics did not reproduce number for number: the thermal noise was built on eigenvectors whose choice differs between machines. Those files are kept unchanged in `output/` and are superseded by `output/statistics_addendum/`.

Addendum: archive SHA-256 f1e5b67286bff3e65a157e7a7f35b9333f50fb4b0b244ef121d344ffa6fa079c, all 50 file sums checked as delivered. The lead reran `depth.py`: the profile table and summary are byte-identical and the wave table agrees to 9 parts in 1e16. The new sampling needs Numba, which is not on the lead's machine, so the samples were not redrawn; the lead ran `statistics_addendum.py --readout-only` on the stored arrays (with a pass-through stand-in for the Numba import) and every readout table reproduces (largest difference 1 part in 1e10). The profile fractions in the table above were also read directly from `g4_beam/merit/profiles.json` by the lead.

The delivered archive carried its own edit of this page. This page is the lead's and replaces it, so its line in `SHA256SUMS` was regenerated and `manifest.json` still records the delivered one; `SHA256SUMS.delivered` is the delivered list.
