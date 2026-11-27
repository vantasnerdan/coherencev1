# g13_rogue pass 3, step 2: the same whole foil at three times finer still (lead's summary)

Step 1 let the whole-foil model hold waves down to 4.6 mm, and step 1b added the push every impact shares. This step (Marvin) refines the model again, to hold waves down to 1.5 mm, and runs each case three ways: with the assumed random launch alone, with the built sources (deuteron momentum and heating of the face) added, and with the built sources alone.

`README.md` is the technical page, `inputs_ledger.md` the inputs and their sources, `evaluate.py` and `spatial.py` the scripts (the driver of step 1b is reused unchanged).

## The slow vibrations have no shorter waves to give

The model can now represent a 1.5 mm wave (it captures 99.8 percent of one). But the foil's vibrations between 5 and 250 kHz do not contain one: the shortest bending wave at 250 kHz is about 2.2 to 2.4 mm long. So this step adds resolution and finds nothing new to resolve. Comparing it with step 1 on the same impacts, the height of the tallest event against the typical level changes by a factor of 0.88 to 1.09, with no case clear of no change.

Shorter waves mean higher frequencies. Extending this model to 475 to 625 kHz was priced and not run: under the one launch shape we have for an impact, that window would receive 2e-16 W.

## The built sources alone, and an uneven beam

With the assumed one percent switched off, what moves the foil is the beam's own momentum and heat, including the natural unevenness in when the deuterons arrive.

| Beam | Tallest over typical | Typical energy density, 0 to 50 um under the beam | Against heat at 300 K in the same vibrations |
|---|---:|---:|---:|
| Even, steady | 1.7 | 2.6e-8 J/m3 | 2,800 times |
| Hot spots, steady | 1.8 | 4.0e-8 J/m3 | 4,300 times |
| Hot spots with bursts | 5.9 | 7.2e-7 J/m3 | 78,000 times |

(Short-lived waves; sixteen records each. With long-lived waves the levels are five to eight times higher and the tallest events 2.8 to 3.5 times the typical level.)

- **Here an uneven beam matters.** Bursts raise the typical level about thirty times and make the tallest events 5.9 times the typical level, where a steady beam gives 1.7. Under the assumed random launch, beam shape made no difference; under the sources built from the beam itself, it does.
- **The levels are small.** The largest elastic energy per deuteron seen with the built sources alone is 4e-13 eV, and the part passed to the deuterium vibration under the earlier assumptions is 2e-20 eV.
- **They are still far above heat in these particular vibrations**, by thousands to hundreds of thousands. That comparison covers only the foil's vibrations between 5 and 250 kHz, not the heat in the lattice as a whole.

## With the assumed random launch

Adding the built sources to the assumed one percent changes the height of the tallest events by -1.6 to +1.8 percent, as in step 1b. The counts are those of step 1b: 95 percent of impacts land on a moving place (98 percent with the state-dependent push), waves converge 31 to 35 percent of the time, and 0.4 to 0.7 percent of impacts land in a rogue peak or hole (1.3 to 2.0 percent with the state-dependent push). The tallest events are 2.9 to 3.7 times the typical level and 0.5 to 7 mm across.

## Where the pass stands

- Inside the slow vibrations the refinement is finished: 13.7 mm, 4.6 mm, then no shorter wave to find.
- Shorter waves made the tall events narrower (from the size of the foil to a few millimetres), richer (6 to 12 times the energy per deuteron) and somewhat taller against the typical level (1.1 to 1.3 times).
- Every energy level under the assumed launch rests on the one percent. The sources built from the beam give 1e-13 eV per deuteron where the assumed launch gives 1e-9 to 1e-8.
- An impact's energy starts in waves far shorter and faster than these. Pass 4 (`../../pass4/`, in progress) carries those across the same whole foil, band by band.

## Three largest omissions

1. What an impact launches at each wavelength.
2. Everything above 250 kHz.
3. How fast the waves die in a palladium foil holding deuterium, and the mount.

## Lead's note

Main archive SHA-256 a932a8884677f61ad40bc76392febabcbbd3bdf325de43631abca27c55f15980, all 306 file hashes in `manifest.json` checked, no machine paths. Evidence archive (frozen vibration shapes, sampled histories, animations; 1.0 GB) SHA-256 e61ab2af46343f68cf682cab4ace136e6fdc8fbddbf5f37148605d87e5845638; it is not in the repository, and `evidence_manifest.json` lists its contents.

The lead reran two cases (hot spots with bursts, short-lived waves: built sources alone, and both sources with the state-dependent push; sixteen records each) with NumPy, SciPy and a C compiler on a clean copy of the repository at commit 650a2eb plus this directory and the frozen shapes from the evidence archive. Records, signed-event and energy tables are byte-identical; peak, tail and same-event tables agree to 2 parts in 1e11 or better. The other 28 cases were not rerun; the ranges above were read by the lead from `summary_sources_r2.csv` and `README.md`.
