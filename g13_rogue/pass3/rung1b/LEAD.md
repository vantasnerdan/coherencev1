# g13_rogue pass 3, step 1b: every impact pushing the same way (lead's summary)

In pass 2 and in step 1 each impact was given a push of random direction, so the pushes cancel on average and only the leftover drives the foil. In the foil every deuteron arrives from the same side and heats the same face. This step (Marvin) adds that same-direction part to the step 1 model, on the whole foil, and reads it apart from the random launch and together with it. It also repeats the step 1 comparison with sixteen records in place of two.

`README.md` is the technical page, `inputs_ledger.md` the inputs and their sources, `signed.py` the single script.

Step 2 (`../rung2/LEAD.md`) runs the built sources alone, without the assumed random launch. There an uneven beam does matter: bursts make the tallest events 5.9 times the typical level against 1.7 for a steady beam.

## What was added

- **Momentum.** Each deuteron stops in the foil and hands over its momentum, inward. At 138 uA that is a steady force of 2.1 millionths of a newton on the beam spot.
- **Heat.** The beam's 0.69 W is deposited within about 100 nm of the face. The heated layer expands, and the heat spreads into the foil by palladium's own conduction. Expansion, heat capacity and conductivity are published values for pure palladium (Arblaster; Furukawa; Ho, Powell and Liley).
- **Beam cases.** Even and steady; hot spots; hot spots with bursts; and an even beam with a ripple at 100 kHz, at 98.63 kHz (on one of the foil's own vibrations) and at 100.48 kHz (between two).

The random launch of step 1 (one percent of the recoil energy into the 5 to 250 kHz range, assumed) is kept as a separate term. Both rules for how fast waves die are run. Every case has sixteen independent records, each run with and without the added terms on the same impacts.

## What the step shows

**A steady beam bends the foil and warms it.** With the edge held at fixed temperature the centre is 7 to 9 K warmer and the foil bows by 16 to 23 nm, nearly all of it from heat. This is a standing shape, not a vibration.

**As a source of vibration in this range the same-direction push is small beside the assumed random launch.** In one followed record (hot spots with bursts) the largest displacements are:

| Source | Short-lived waves | Long-lived waves |
|---|---:|---:|
| Random launch (assumed one percent) | 349 nm | 1071 nm |
| Heat | 0.51 nm | 0.92 nm |
| Momentum | 0.003 nm | 0.006 nm |

At the tallest events the added terms account for at most 6 parts in a million of the elastic energy. Adding them changes the height of the tallest event against the typical level by between -1 and +1.4 percent, in all 24 comparisons, including the ripple tuned onto a foil vibration.

**This puts the one percent in a new light.** The momentum and heat terms are built from the beam and from measured properties of palladium. The random launch is an assumed fraction. In this range the built terms move the foil several hundred times less than the assumed one. So the levels in pass 2 and in step 1 (energy per deuteron, the comparison with heat) stand or fall with that one percent, and a single short, small impulse would not supply it either (`../rung1/pulse_sensitivities.csv`). What an impact really launches goes mostly into much shorter waves than this model holds yet.

**With sixteen records, shorter waves do make the tall events taller against the typical level.** Comparing the step 1 model (waves down to 4.6 mm) with the pass 2 model (13.7 mm) on the same impacts, both source terms included:

| | Range over 24 cases |
|---|---:|
| Tallest over typical, step 1 against pass 2 model | 1.09 to 1.33 times (22 of 24 clear of no change) |
| Typical energy density | 4.6 to 11 times |
| Elastic energy per deuteron at the tallest event | 5.9 to 12 times |
| Tallest over typical, in the step 1 model itself | 2.9 to 3.6 |

Two records in step 1 could not show this; sixteen do, under both loss rules.

**The counts.** Impacts onto a place already moving: 94 to 95 percent with the push independent of the surface state and 98 percent with the push dependent on it (short-lived waves). Converging waves: 32 to 36 percent of the time. Impacts into a rogue peak or hole: 0.2 to 0.5 percent, and 1.3 to 1.7 percent with the state-dependent push. With long-lived waves the last figure is 20 to 60 percent, but those runs had not settled and the figure is not a steady rate.

**Widths and energies** are as in step 1: the tallest events are 1.7 to 7 mm across, most at the limit of what the step can hold. At the tallest events the elastic energy per deuteron at 30 um is 3e-8 eV (short-lived waves) and 1.4e-7 eV (long-lived), and the part passed to the deuterium vibration under the earlier assumptions is 3e-15 and 8e-16 eV.

## Three largest omissions

1. What an impact launches, by wavelength. The built terms are small here and the assumed term carries the result.
2. How fast these waves die in a palladium foil holding deuterium.
3. The mount. The foil is free in the model; a clamp changes its vibrations and how heat leaves.

## Lead's note

Archive SHA-256 ccd791792998b61b5bb13be5406402a5c92265c1b0a183b285ffd59f10109a97; all 711 file hashes in `manifest.json` checked; no machine paths.

The archive is 1.4 GB. This directory holds everything except the 96 sampled energy histories (`*_histories.npz`) and the two full animations (`clip_r1_hot_bursts_b035_retained.png`, `clip_r1_hot_bursts_b035_Cu_proxy.png`); their hashes are in `manifest.json`, the two `*_preview.png` stills are here, and `signed.py run` regenerates the histories.

The lead reran one case with and without the added terms (step 1 model, hot spots with bursts, state-dependent push, short-lived waves, all sixteen records) with NumPy, SciPy and a C compiler on a clean copy of the repository at commit fe514a1 plus this directory. Records, signed-event and energy tables are byte-identical; peak, tail and same-event tables agree to 1 part in 1e10 or better. The other 190 cases were not rerun; the ranges above were read by the lead from `summary_r1.csv` and `contrasts_r1.csv`.

The lead asked for this step expecting the same-direction push to be where an uneven beam enters. In this range of slow vibrations it is not.
