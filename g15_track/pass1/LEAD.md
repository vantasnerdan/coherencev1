# g15_track pass 1: where an impact's energy goes (lead's summary)

`g13_rogue` followed only the orderly push an impact gives the palladium and filed the rest as heat. This pass (Marvin) follows all 5 keV of an impact through every carrier it passes, from the first quarter femtosecond until the last of it has left the foil about 45 seconds later, and then does the same for the whole foil under the continuous beam. Nothing is binned as heat: the run accounts for the energy to 1 part in 1e10.

`README.md` is the technical page, `inputs_ledger.md` the inputs with their sources, `run.py` the script, `single_impact_channels.png` the picture of one impact.

## The route of one impact

| Stage | Carrier | Energy of the 5,000 eV | Holds it for |
|---|---|---:|---|
| Arrival (0.3 to 0.4 ps) | electrons set moving by the ion | 4,290 eV | 15 to 500 fs |
| | recoiling palladium atoms | 606 eV | 0.05 to 1 ps |
| | the slowed ion itself | 100 eV | 0.1 to 2 ps |
| | deuterons struck directly | 1.8 eV | |
| | orderly palladium sound (all of g13) | 1.8 eV | microseconds |
| Handover | **the deuterons' vibration in their sites** | 1,140 to 4,030 eV pass through it | 0.08 to 2 ps |
| | palladium vibrations, top of the band | 3,600 to 3,900 eV pass through | 0.15 to 5 ps |
| | palladium vibrations, lower band | 320 to 1,180 eV pass through | 0.3 to 10 ps |
| Spreading | warm lattice around the track | all 5,000 eV | about 1 ns |
| Leaving | heat flow through the foil to the mount (95 percent), the faces and radiation | all 5,000 eV | seconds |

Energy passes through several carriers in turn, so the middle rows add to more than 5,000 eV.

**Between a quarter and nine tenths of the beam's power passes through the deuterium vibration**: 0.16, 0.30 or 0.55 W of the 0.69 W, in the three cases carried. No measurement fixes how hot electrons in PdD divide their energy between the deuterium and the palladium, or how long the deuterium vibration holds it. The three cases span published calculations: the share of the electron coupling that sits in the hydrogen vibration (Errea 2013: 75 to 79 percent), and a calculated lifetime of that vibration in PdH (Paulatto 2015, 0.08 ps) against the 2 ps this project assumed before.

## The whole foil under the beam

- **A large flow, a small standing amount.** Because the vibration gives its energy up within a picosecond or two, what is stored in it across the whole foil at any moment is 1e-14 to 1e-12 J, against 0.16 to 0.55 W flowing through.
- **It is concentrated, not spread.** Each impact pumps a patch 8 to 60 nm across (assumed sizes) to 0.002 to 0.019 eV per deuteron, a twentieth to a half of one quantum of the vibration, in a lattice whose thermal energy is 0.026 eV. A given patch is then quiet until the next impact: of 221,184 moments looked at, 5 found a patch still ringing, and two tracks overlapping while both ring has a chance of 3e-10 to 8e-6.
- **On the lattice layer's map** (pump depth over pump threshold): if all the deuterons in a patch are pumped in step, one impact puts the patch at 1.2 to 5 times the threshold at which the vibration feeds itself. If one percent of the motion is in step, 0.012 to 0.05 of it. Whether an impact pumps its patch in step is not in this model; it is assumed in neither direction.
- **Temperature.** The foil sits 11.5 K above its surroundings with an even beam and 12.2 K with hot spots, for the middle guess at how well the foil is held; 7 K for a well-held foil and about 170 K for a poorly held one. How the foil is mounted is not given in the paper.
- **The paired arrivals** (15 percent of impacts): the second deuteron of a D2+ travels through the electrons the first has just set moving; its field is in `second_constituent_prior.csv`.

## Yield

Three ways the pumped lattice could add D-D reactions were computed with the cross section and profile of `g12_residual/inputs_ledger.md`, with the vibration's energy shared out at random among the deuterons:

- a beam deuteron meeting a target that is moving: at most 2e-10 of the unexplained signal;
- deuterons struck directly, hitting others: at most 6e-8 of it;
- neighbouring deuterons in a pumped patch moving toward each other: common out to the distances for which the pair energy has been computed (1.3 to 1.5 A), with no figure beyond that, because the pair energy at shorter distances and the number of vacancies are not known. The table (in the release archive, 8 MB) gives what those two would have to be.

With the energy shared at random the pumped lattice does not account for the signal. The case not yet built is the one the threshold figures point to: the patch moving in step.

## What the result rests on

1. **How hot electrons in PdD divide their energy** between deuterium and palladium vibrations, and how the palladium recoils do.
2. **How long the deuterium vibration holds energy, and whether an impact starts it in step.** This is the loss rate that has been open since the first campaign (AGENTS.md, open item 3).
3. **The pair energy closer than 1.3 A, and how many vacancies there are.**

## Lead's note

Archive SHA-256 9d6418bf5e6440bda1501056994b680d1f429b52b3232f93aeed6cb352139fa7, all 31 file hashes in `manifest.json` checked, no machine paths. Two files are left out of the repository for size: `single_response.npz` (18 MB) and `pair_consumers.csv` (8 MB); `run.py --all` regenerates both.

The lead reran `run.py --all` with NumPy, SciPy and Matplotlib on a clean copy of the repository at commit d8b935d plus this directory. Every table and array agrees with the delivery to 1e-12 or better, apart from entries that are zero to rounding.
