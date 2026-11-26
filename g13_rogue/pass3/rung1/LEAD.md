# g13_rogue pass 3, step 1: the same whole foil, holding waves three times shorter (lead's summary)

Pass 2 modelled the whole palladium foil under a live beam and counted impacts into moving places, converging waves, and rogue peaks and holes. Pass 3 (Marvin) refines that same model in steps. Each step lets it hold shorter waves, still across the whole foil and under the whole beam. No small patch replaces the foil. The question across the steps is whether the tall events get taller against the typical level, narrower, and richer in energy per deuteron.

`README.md` is the technical page, `inputs_ledger.md` the inputs and their sources, `rungs.py` the single script.

## How short the waves of pass 2 really were

Pass 2 kept the foil's vibrations between 5 and 250 kHz, and its pages speak of millimetre waves. This step measured what the pass 2 model can carry as a travelling wave. A wave counts as carried when the model captures at least 95 percent of it, for every direction along the foil.

| Wavelength | Captured by pass 2 | Captured by this step |
|---:|---:|---:|
| 20 mm | 99.9 % | 100 % |
| 13.7 mm | 95 % | 100 % |
| 10 mm | 57 % | 100 % |
| 7 mm | 16 % | 100 % |
| 4.6 mm | 1.5 % | 96 % |
| 3.5 mm | 0.6 % | 36 % |
| 1.5 mm | 0.007 % | 2 % |

So pass 2 carried waves of 13.7 mm and longer, about two thirds of the foil's width. Its "converging waves" are the whole foil sloshing. This step carries waves down to 4.6 mm. It does so with 190 vibration shapes in place of 53, inside the same 5 to 250 kHz range.

## What was kept the same

The foil, its grains, six free faces, the beam (138 uA, 2.5 keV), the impacts drawn, and the launch: one percent of the recoil energy goes into this frequency range, as in pass 2, with the push either independent of the surface state or dependent on it. The step adds no energy; it only lets the same energy take shorter shapes. Three beam shapes are run (even and steady; hot spots; hot spots with bursts), each twice.

Two rules for how fast the waves die are run side by side. The first is the pass 2 rule (50 us life, 0.2 m range). The second is taken from a published measurement on copper (Ono 2020), stretched down to these frequencies, which lets the waves live for 0.1 to 6 ms. Neither is a measurement on palladium.

## What the step shows

Figures are for the layer 0 to 50 um deep under the centre of the beam, with the pass 2 loss rule. Each range covers the three beam shapes, both launch rules and both records.

| | Pass 2 model | This step |
|---|---:|---:|
| Tallest energy density | 0.08 to 0.59 J/m3 | 0.66 to 4.8 J/m3 |
| Tallest over typical | 2.3 to 3.6 | 2.2 to 4.7 |
| Width of the tallest event | 7 mm to the whole foil | 1.8 to 4.5 mm |
| How long it lasts | 6 to 18 us | 4 to 9 us |
| Elastic energy per deuteron at 30 um | 4e-10 to 5e-9 eV | 7e-10 to 6e-8 eV |

- **The tall events get narrower and richer.** With shorter waves allowed, the same launched energy gathers into events a few millimetres across, with five to fifteen times the peak energy.
- **Their height against the typical level does not move by an amount two records can show.** The typical level rises by about the same factor as the peaks. With the pass 2 loss rule the ratio goes up in some cases; with the copper rule it goes down in some of the same cases (10 of 18 comparisons change direction).
- **The widths are at the limit of what the step can hold.** 14 of the 24 widths are below half the shortest carried wavelength, so the true events may be narrower still.
- **With long-lived waves the levels are higher again.** The copper rule gives peak energy densities of 4 to 17 J/m3 and up to 1.7e-7 eV per deuteron at 30 um. These runs had not settled: the waves live longer than the record.
- **Impacts onto a moving place** stay at 91 to 98 percent (pass 2 rule). Impacts into a rogue peak or hole stay at 0 to 1 percent, higher with the state-dependent push, as in pass 2.
- **What made the tallest events** is not separated: they carry signs of several causes at once.

The part of that energy passed to the deuterium vibration, under the same assumptions as before (all locally lost energy transferred, held 4.1 ps), is 4e-16 to 6e-15 eV per deuteron at 30 um.

## Three largest omissions

1. How fast these waves die in a palladium foil holding deuterium, measured by wave type and frequency. Two assumed rules already disagree on the direction of one result.
2. What an impact launches. One percent of recoil energy into this range, with a random sign for each impact, is assumed. A single short, small impulse would by itself put far less than one percent here (`pulse_sensitivities.csv`), so the launch into these slow vibrations has to come from somewhere else.
3. Waves shorter than 4.6 mm, and frequencies above 250 kHz.

## Lead's note

Archive SHA-256 b8d24dcdd0f1e3cd3fdd2ef85a5f23831efaf807bbeee39c71a1bf1e55c5be2c; all 27 file hashes in `manifest.json` checked; no machine paths.

The lead reran `rungs.py compare` and `rungs.py readout` (NumPy, SciPy and a C compiler) on a clean copy of the repository at commit b0a67c0 plus this directory. All 24 cases reproduce: the record table, the signed-event table and the energy budget are byte-identical, and every derived table agrees to 2 parts in 1e8 or better. The frozen vibration shapes (`model.npz`) were used as delivered and not rebuilt.

The delivery stopped before the next step on the ground that the loss rule is not measured. The lead's decision is to continue with both loss rules carried side by side, since their disagreement is itself a result, and to take up omission 2 first. In this model each impact pushes with a random sign, so pushes cancel on average. In the foil every deuteron arrives from the same side and heats the same face, so that part of the push has one sign for all impacts. A steady, even beam then only bends the foil, but any change of the beam in time or across the spot drives the foil's vibrations in step. That part is not in pass 2 or in this step.
