# g15_track pass 3: hot electrons changing the stiffness of every site at once (lead's summary)

Pass 2 (`../pass2/LEAD.md`) showed that an impact hardly pushes its patch of deuterons in step. It left one handle: while the electrons around a site are hot, the stiffness of the site itself may change. That acts on every deuteron of the patch at nearly the same moment and needs no lopsided site. A common, sudden stiffness change makes the spread of each deuteron's existing zero-point and thermal motion breathe in step, at twice the vibration frequency. This pass (Marvin) builds that on the same whole model: pass 1's energy budgets, tracks, fixed profiles, three sets of rates and three beam shapes.

`README.md` is Marvin's technical page, `inputs_ledger.md` the inputs, `run.py` the script.

## The one open number

How much the stiffness changes per eV of electron energy held in a site has not been measured or computed for PdD (node `g16_hot_electrons` is computing it). It is carried here as 0.1, 1 and 10 percent per eV, stiffer or softer, beside no change.

## What came out

1. **The electron energy per site is small.** One impact's 5 keV is shared among roughly a hundred thousand sites, a few hundredths of an eV each, more close to the track. So even at 10 percent per eV the stiffness of a typical site changes by a fraction of a percent, and it stays positive everywhere.
2. **The breathing is slight.** The spread of a deuteron's position changes by at most 0.11 percent (0.22 percent in its square), in the vacancy cage, and 0.10 percent in bulk. The largest extra energy in a site is 0.005 eV.
3. **It is in step, and brief.** The breathing of the widths across a patch is almost perfectly in step (0.9998 of full alignment) and falls to about a third of its peak in 0.3 ps. In the four patch-wide patterns of motion followed separately it falls that far in 0.011 ps in bulk and 0.3 ps in the cage. A rise of 15 fs gives about a third of what an instant step would; a rise of 500 fs gives thirty times less again, not zero.
4. **Pair distances barely move.** The chance of finding the closest bulk pair at 2.41 Å or the cage pair at 1.73 Å rises by about 5 and 13 parts in 100,000 at the best moment, for a site next to the track.
5. **Whole foil and yield.** With all six beam settings and three beam shapes, the breathing adds at most 1.1e-14 of the unexplained signal of `g12_residual/`. Two impacts' patches overlap in at most 4 cases in 10,000, bursts included.
6. **The same event's sound.** The breathing presses on the palladium and launches a strain of at most 3e-10, carrying 5e-12 eV per impact, beside the launches already built in `g13_rogue/`.

Energy is followed throughout: the work the stiffness change does on the deuterons is taken from, or returned to, the electrons, to 2.4e-7 eV per impact.

## What the result rests on

1. **The stiffness change per eV of electron energy in PdD.** The result scales with it. To matter it would have to be far above the 10 percent per eV carried here.
2. **The state of the six deuterons in a vacancy cage.** The cage is described by its closest-pair cut only, not by a full set of cage stiffnesses.
3. **From the closest computed pair distance to contact, and how many vacancy cages a foil holds.**

## Lead's note

Archive SHA-256 78c6b5a5a06a243cc20c8b4b9ca56f7823e598a9fb6c7071d1373c3e7a747856, all 31 file hashes in `manifest_capsule.json` checked, no machine paths. Four tables are left out of the repository for size: `phase_evolution.csv` (29 MB), `yield_consumers.csv` (10 MB), `acoustic_transmission.csv` (8 MB) and `whole_foil.csv` (5 MB). The three evidence archives (3.9 GB together; SHA-256 0d7121cb…3bc95e, 3c4664ad…aa82b0, 497216d6…e26718 as reported by Marvin) were not downloaded.

The lead has not rerun this pass. The numbers above are read from the delivered tables and `results.json`.
