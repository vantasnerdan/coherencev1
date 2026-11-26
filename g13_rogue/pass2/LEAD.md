# g13_rogue pass 2: the whole foil under a live, structured beam (lead's summary)

Pass 1 treated every impact as the same small packet launched into a quiet lattice and asked whether the deuterium vibration then concentrates. Pass 2 (Marvin) turns that around, on direction from the person leading the work, whose picture was hail on a pool: there is a mean splash and a mean ripple, and then there are splashes far above the mean, from double impacts, from waves converging from several directions, and from an impact landing in a rogue hole (a rare deep trough of the waves already there). With impacts arriving continuously the surface is never at rest, and the beam itself is not an even rain.

`README.md` and `readout.md` are the technical pages, `live_method.md` the equations, `inputs_ledger.md` the inputs and their sources.

## What the model is

The whole palladium foil (250 um thick, the area and mass of the 24 mm disk, beam on the central square centimetre) as a vibrating elastic body, with palladium's direction-dependent stiffness, grains of two sizes, and reflection at both faces and at all four edges. It keeps the foil's own vibrations between 5 and 250 kHz: 53 modes, with wavelengths of millimetres. Impacts arrive at the charted beam current in real time, none left out.

Three things are new against pass 1:

- **The launch depends on where the impact lands.** Each impact's push on the foil is drawn from a spread of sizes, and a larger push is more likely where the foil is already displaced or strained. The strength of that dependence is an assumed number, not a measured one.
- **The beam can be uneven.** Three shapes across the spot (even, peaked, peaked with two hot spots) and three in time (steady, rippling at 100 kHz, rippling with random bursts), all at the same total charge.
- **Three kinds of event are built and counted separately**, along with peaks and holes that occur whether or not an impact lands in them:
  - A: an impact into a place still disturbed by earlier impacts;
  - B: waves from different parts of the beam spot converging, from directions at least a right angle apart, on a place outside the spot;
  - C: an impact into an existing peak or hole more than three times the usual size of the motion there.

It drops the short waves (anything above 250 kHz, so nothing of the splash at atomic scale), lasting damage, the foil's slowest drum-like motions below 5 kHz, and any calibrated link from foil vibration to deuterium vibration.

## The time frames

| Wave | Speed | Centre of the foil to the edge and back | Impacts on the spot in that time, at 38 / 138 uA |
|---|---:|---:|---:|
| Fast bulk wave | 4.6 km/s | 5.2 us | 1.2 / 4.5 billion |
| Slow bulk wave | 2.0 km/s | 12 us | 2.9 / 10 billion |
| Bending wave, 10 mm wavelength | 0.33 km/s | 74 us | 17 / 64 billion |

A bulk wave crosses the thickness in 54 ns. Whether a wave survives to the edge depends on how far it travels before it dies, which is not known for this foil. With an energy range of 2 mm almost nothing returns (6 parts in a million); with 2 cm, 29 percent of the fast wave's energy comes back from the edge; with 20 cm, 86 percent.

## What it shows

**The foil is never at rest.** Once the vibrations live 50 us or more and travel 2 cm or more, 90 to 99 percent of impacts land on a place that is already moving because of earlier impacts. With a 2 mm range the figure stays near 40 percent. At the lower beam energy and current (0.25 keV), it is 28 to 44 percent.

**Convergence from several directions is ordinary.** A place outside the beam spot is under converging, reinforcing waves from different parts of the spot 16 to 56 percent of the time.

**Impacts into rogue peaks and holes happen, at a percent or so.** Between 0 and 5 percent of impacts land in an existing peak or hole beyond three times the usual size. Holes and peaks occur in both signs. A place is strained by about 2 parts in a million and displaced by a few hundred nanometres at such moments.

**In these slow, foil-wide vibrations the beam is far above heat.** At the same place and in the same vibrations, the beam-driven level is 1e10 to 1e11 times what heat at 300 K puts there. Heat's energy sits almost entirely in fast, short vibrations that this model does not hold; in the slow ones the beam has the field to itself. This rests on the assumed 1 percent of each impact's recoil energy going into these vibrations.

**How tall the tall events are.** The tallest peak or hole in any region of the 48 cases is 4 to 5 times the usual size of the motion; none reached 10 times. In energy, the tallest moment in a small region is about 4 to 6 times that region's mean; larger regions are steadier (2.2 times for the whole foil). Each case is 0.2 ms of foil time, so events rarer than a few thousand per second per region could not have been seen.

**Whether an uneven beam changes any of this is not settled by this pass.** The delivered table shows the three kinds of event moving up and down between beam shapes with no consistent direction. The lead's rerun, which drew different random numbers on a different machine, moved each entry by as much as the beam shapes differ from one another (for example, convergence with a peaked steady beam: 56 percent delivered, 42 percent in the rerun; impacts into peaks and holes with an even bursting beam: 1.1 and 4.7 percent). Two short records per case are too few to tell the beam shapes apart.

**The stiffness law adds little.** In the fixed-launch comparison, a foil whose stiffness depends on strain, at assumed strengths far beyond anything measured, changed the tallest peak by at most 2 percent.

## Working backwards

The reference from pass 1 is 0.64 meV per deuteron in coherent deuterium vibration, computed for a fully loaded lattice (the foil at depth holds 0.3 to 0.6, so it is a reference, not a known threshold).

- The tallest foil-wave moment holds about 1e-10 eV per local deuteron as elastic energy.
- If that energy drained into the deuterium vibration as fast as the foil waves decay, and the deuterium vibration held it for its 4.1 ps, the deuterium vibration would hold about 5e-16 eV per deuteron: a factor of 1e12 below the reference.
- If instead everything the best observed event (7.65 us long) dissipates were handed to the deuterium vibration and held for ever, it would reach 3e-10 eV per deuteron: a factor of 2 million below.

So within this model the open quantities, in order of leverage, are: how much of an impact goes into waves and how that depends on the state it lands in (assumed here); the real shape of the beam in space and time (assumed here); how long the deuterium vibration, or anything coupled to it, can hold energy (4.1 ps here); and the part of the wave field this model does not hold at all, the short waves near the impact.

## Three largest omissions

1. The measured launch: what fraction of an impact's energy becomes waves, at which wavelengths, and how that changes on a moving surface. The 1 percent and the state dependence are assumptions.
2. Everything above 250 kHz, including the micrometre and nanometre waves that carry most of an impact's energy and reach the deuterium vibration most directly; and everything below 5 kHz.
3. The real foil: its grains, how it is clamped, how far its waves travel, and its elastic constants when loaded (pure palladium values are used throughout).

## Inputs

Beam current, ion species, recoil energies and the fixed deuterium profile are those already in the repository. New in this pass: the elastic constants of palladium (227.1, 176.0 and 71.7 GPa, Rayne 1960, taken from a later paper's transcription) and its density, both from published measurements cited in `inputs_ledger.md`.

## Lead's note

Archive SHA-256 cf46ddea5eb224f393ae9d3a2380e82ffdc06b5150464e1e06c0165a1f2010b8; all 46 file hashes in `manifest.json` checked; no machine paths.

The lead reran both models with NumPy and SciPy on a clean copy of the repository at commit 9c26f8c plus `coarse.py` and `live.py`.

- **Fixed-launch model** (`coarse.py`): every table reproduces, largest difference 3 parts in 1e9.
- **Live model** (`coarse.py --live --exposure 0.0001 --realizations 2`): the time-frame table, edge-survival table and beam parameters are byte-identical, and every case received the same number of impacts. The event tables did not reproduce number for number: the random marks and signs come out differently on the lead's machine. The rerun is therefore a second independent draw. Across the eleven summary rows the two draws give A 89.6 to 98.7 percent at 2.5 keV, C 0 to 4.7 percent, B 16 to 56 percent; the largest C over the whole scan is 3.5 percent delivered and 5.3 percent in the rerun; the tallest peak or hole is 4.2 and 5.3 times the usual size; the best conversion figure is 4.7e-16 and 5.9e-16 eV per deuteron. The ranges in this page are taken across both draws. The delivered files are kept unchanged.

Asked of Marvin after this delivery: draws that do not depend on the machine, and enough records per beam shape to tell the shapes apart.
