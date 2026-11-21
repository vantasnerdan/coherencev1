# g7_oxide_beam: what an oxide cap does to a beam measurement besides static screening

The static oxide sites are all ordinary (`g7_micro/pd_oxide/`, `g7_micro/zr_calibration/`): 38-39 eV at 0.1 A. Kasagi et al. (J. Phys. Soc. Jpn. 71, 2881 (2002)) fit about 310 eV for Pd and 600 eV for PdO on Pd. This node sizes two other ways a cap could make a 40 eV site read high in such a fit: deuterium piled up near the surface, and charge on the cap. Marvin's calculation, with the `g4_beam` forward model; NumPy and SciPy only.

## Result

**Deuterium piled up near the surface.** The published fits assume a uniform deuterium density. A deuterium-rich surface layer raises the low-energy yield relative to the 10 keV reference and so reads as extra screening.

| Profile (true screening 40 eV everywhere) | Apparent screening in a uniform fit (eV) | Worst mismatch in curve shape |
|---|---:|---:|
| Uniform | 40 | 0 |
| 1 nm layer, 100 times bulk density, exposed Pd | 206 | 25 percent |
| 1 nm layer inside the oxide, 30 times bulk density | 209 | 24 percent |
| Infinitely thin, infinitely rich sheet (the ceiling) | 240 | 29 percent |
| Rich layer buried under 30 nm of oxide | below 40 | fit fails |

A surface-rich profile can turn 40 eV into about 200 eV. It cannot reach 310 or 600 eV at any thickness or density: at 2.5 keV the ceiling is 4.1 times the uniform yield, and the 310 and 600 eV curves need 9.9 and 62 times. The curve shape it gives is also wrong by a quarter.

**Charge on the cap.** A constant energy gain of the deuteron before the collision imitates extra screening almost exactly (shape within 1 percent). Reading 300 eV more needs about -650 V at the collision.

- The beam drives the wrong sign. Arriving D+ and escaping electrons charge a surface positive, which slows the beam.
- Holding 650 V across a 30 nm cap is a field of 22 GV/m, about a hundred times the breakdown field of a good engineered oxide.
- PdO conducts. Measured film resistivity is about 0.001-0.0025 ohm m; holding the voltage at the beam current needs about 1e9 ohm m.
- A field across a cap on grounded metal gives no net energy gain to a deuteron that reacts in the metal.

Charge on the oxide does not explain the PdO result.

## What it means

- Part of any beam-fit screening number can be depth profile. A nanometre-scale deuterium-rich skin is worth up to about 160 eV of apparent screening on top of the true value. This is a lever the experiments do not control and the fits do not include.
- Neither mechanism, alone or together, reaches the measured oxide contrast. With the static sites, that is three candidates sized and found too small: site chemistry, depth profile, cap charging.
- What remains: how the yield is normalised and which bare cross-section and stopping law the fit assumes (`g7_gap/`; Kasagi's own reference curve and stopping differ from the ones used here), sites with screening above 40 eV that no static cell has yet produced (`g7_micro/` part a), and the collision dynamics.
- One measurement separates profile from charging: on one foil, at fixed and independently measured beam energy, step the beam current up and down while measuring the shallow deuterium depth profile and the surface potential, with a grounded bare-Pd control.

## Three omissions

1. Straight-line slowing with tabulated stopping. No scattering, no roughness, and oxygen stopping estimated (halving or doubling it moves the in-oxide case by 5 eV).
2. Fixed profiles. No beam-driven change of the profile during a run.
3. The comparison curves are this model's 310 and 600 eV curves, not a refit of the original counts.

## Lead's note

Delivered archive SHA-256 1bf0b4d6c6fe8d07f0f4686d120f2e132acc5df8ca8c0463013ebe4197e4c535, unchanged; the producer's full record is `delivered_README.md` and the cited sources are in `material_source_notes.md`. Run `python g7_oxide_beam/oxide_beam.py` from the repository root (needs NumPy 2 for `np.trapezoid`; the lead ran it under NumPy 1.26 with that name aliased). The lead's rerun reproduced every number in `results.json` to better than 1e-6 and three of the four tables byte for byte. Zirconium with oxygen was not run, because no Pd mechanism reached the measured contrast.
