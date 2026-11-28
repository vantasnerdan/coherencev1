# g15_track pass 5: sound with wave shape over the whole foil (lead's summary)

Pass 4 (`../pass4/LEAD.md`) built travelling, coupled fields over the whole foil, but carried sound above 250 kHz as energy without crests and troughs, resolved nothing between 0.4 µm and 1.25 mm, and counted the meetings of deuterium waves without building them. This pass (Marvin) works on those three, on the same model and the same two seconds of beam. It also takes the computed split of the hot-electron stiffening from `../../g16_hot_electrons/followup/` as a third case beside pass 4's two extremes.

`README.md` is Marvin's technical page, `run.py` the script, `consumer_comparison.png` the picture.

## What is built

- **Signed sound in every band.** Each impact launches sound in all eight bands, 250 kHz to 2.5 THz, with the sign and timing its track gives it. At a receiver the waves of all impacts over the whole beam, direct and reflected from faces and edges, are added as amplitudes before anything is squared.
- **Receivers from the foil down to the atom.** The same signed sum is evaluated at nested points from the 1.25 mm grid down to 0.125 nm, at actual depths, with every impact near and far contributing. There is no separate small model.
- **Heat to the steady state**, resolved through the first 250 nm, for three mounts.
- **Meetings of deuterium waves from different impacts**, drawn as events and followed through the meeting (162 cases and 54 matched controls).

## What came out

1. **Signed sound is a few hundred times larger than pass 4 could show.** The largest strain sampled is 4.7e-6, against 1.8e-8 in pass 4. It shifts the deuterium vibration frequency by 1.0e-5 (pass 4: 2.5e-8). That is still 0.2 percent of the shift made by slow heat.
2. **The tall events are not yet pinned down.** Each receiver has 256 sampled moments over the two seconds. When the sampling of nearby impacts was made eight times finer, the measure of the typical size fell by a third and the measure of the rare tall events rose eightfold. The largest values above are therefore sampled values, not the height or frequency of the tallest events.
3. **Heat.** The burst-beam hot spot reaches 308.5 K after 2 s and 311.8 K in the steady state (nominal mount). Across the first 250 nm the temperature differs by 0.00008 K. No deuterium wave was turned back or held by the heat-made frequency landscape in the paths followed.
4. **Meetings.** Deuterium waves from different impacts overlap 1e7 to 1e12 times in two seconds over the foil, depending on beam energy, beam shape and the set of rates (2.5e10 for the nominal burst case at 2.5 keV). A wave travels at most 29 nm in these cases. For any one incoming deuteron the chance of passing through such a meeting is about 3e-13.
5. **The computed split of the stiffening.** With each impact's launch held as in pass 4, changing only how the waves travel alters the yield by 4e-20 of the unexplained signal. Applied to the saved electron energies around one 2.5 keV track, it changes the local deuterium frequency by up to 0.8 percent. The launch itself was not rebuilt with it.
6. **Yield.** The largest addition is 7.6e-7 of the unexplained signal of `g12_residual/` (2.5 keV, burst beam), with a sampling error of 6.5e-7, nearly all from the foil being warm. Signed sound contributes 1.5e-12 at most and the meetings 4e-18.

## What the result rests on

1. **The inside of the launch region and the rare tall events.** Sound is summed in its far-field form; within a few wavelengths of a track, where strain is largest, the near field is left out, and the tall-event tail has not settled under finer sampling (item 2).
2. **What counts as yield.** The yield here is the beam deuteron, or a struck one, meeting lattice deuterons that move with the built fields. Nothing in it lets two lattice deuterons react with each other, because the pair energy inside the closest computed distance is not joined to contact. A signal that hardly depends on beam energy is not well matched to a yield that needs the beam deuteron to arrive.
3. **The launch with the computed stiffening.** Electron energy around a track is saved in space and time only at 2.5 keV, so the launch was not rebuilt.

Marvin's own three: the source interior, phases and rare recoils; late returns, slow thermoelastic motion and phase correlations between far impacts; exchange between deuterium waves of close wavelength and the step from neighbour distance to contact.

## Lead's note

Core archive SHA-256 f151f98294cc21cb593bd03a0b0fb84231f1876bb7e281111186576da6c084f6; evidence archive SHA-256 9a0c9d17af70e5d535aa762fc368881d8c7d8def52ab85b0f2cdf66598f5848a (not in the repository). No machine paths. The core archive carries no manifest of its own; the output hashes listed in `results.json` were checked for the ten outputs it contains.

The lead reran `run.py --consume` with NumPy and SciPy on a clean copy of the two archives alone. It verified the inherited inputs and rebuilt `yield_outcomes.csv`, `frozen_launch_controls.csv`, both `probe_projection_*.csv` and `consumer_checks.json` byte for byte. The receiver, heat and meeting runs themselves were not rerun.
