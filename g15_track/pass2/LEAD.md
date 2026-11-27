# g15_track pass 2: does an impact set its patch of deuterons moving in step? (lead's summary)

Pass 1 (`../pass1/LEAD.md`) showed that a quarter to nine tenths of the beam's power passes through the deuterium vibration, and that the patch an impact pumps would be above the lattice layer's pump threshold if its deuterons moved in step. This pass (Marvin) builds the in-step part. It keeps pass 1's energy budgets, tracks, profiles, three sets of rates and three beam shapes, and works out for each deuteron in a patch the push it actually receives, then lets the patch move in the lattice layer's own equations (`g3_lattice/coherence/`) with the computed pair energies of `g10_quantum_site/` for the closest pair, in bulk PdD and in the six-deuteron vacancy cage.

`README.md` is Marvin's technical page, `inputs_ledger.md` the inputs and the force laws, `run.py` the script, `phase_map.png` and `phase_evolution.png` the pictures.

## The two pushes

- **The passing ion.** Its electric force gives each nearby deuteron a directed kick. This is computed from the same screened interaction as the stopping.
- **The hot electrons.** They carry most of the impact's energy, but heating the electrons evenly around a deuteron in a symmetric site pushes it nowhere. A directed push needs either a slope in the electron heating across the site or a lopsided site such as a vacancy cage. How strong that push is in PdD has not been measured or computed. It is carried as five settings: none, one natural unit in either direction (the deuteron feels the full slope of the electron energy held in its site), and a million units in either direction, limited to the energy actually available.

## What came out

Highest beam setting, middle set of rates, single deuteron, a patch of 128,000 deuterons, unless stated.

| | Energy |
|---|---:|
| Passing through the deuterium vibration of the patch (pass 1) | 2,136 eV |
| Directed work by the ion's kick | 0.024 eV |
| Directed work by the electrons at one unit, bulk | 1e-7 to 5e-7 eV |
| Directed work by the electrons at one unit, vacancy cage everywhere | 0.05 eV |

1. **Almost none of the pump arrives in step.** About one part in 100,000 of the energy passing through the patch's deuterium vibration is directed motion; the rest reaches the deuterons with no common timing. In the patch-wide patterns of motion, the electron push at one unit puts in at most 5e-8 eV.
2. **The stiffness is barely modulated.** The lattice layer locks when the stiffness of the sites is modulated at twice the vibration frequency above a threshold. The largest modulation realised in any case is 0.2 percent of that threshold, and that includes the million-unit settings, where most of the vibration's energy is directed. Pushing deuterons in step moves them; it does not by itself modulate their stiffness.
3. **The motion dies.** With the ion's kick alone or with one unit of electron push, the directed motion peaks at 0.023 eV for the patch and is down to 0.02 percent of that after 4.4 ps. Nothing grows.
4. **Closest approach.** The nearest pair in bulk comes from 2.89 Å at rest to 2.58 Å (2.57 Å over all settings); the vacancy pair from 1.97 Å to 1.84 Å (1.73 Å over all settings). Largest relative kinetic energy 0.05 eV in bulk and 0.04 eV in the cage.
5. **Whole foil.** Every impact brings some deuteron in its patch above 0.01 eV (8.6e14 patches a second). A molecular ion brings one above 0.1 eV (largest 0.106 eV, against 0.040 eV for a single deuteron), 1.3e14 patches a second. None reaches 1 eV. Two impacts overlap in one patch within its lifetime in fewer than 3 in ten million cases, bursts included.
6. **Yield.** With the built motion, the beam meeting moving targets adds at most 2.5e-8 of the unexplained signal of `g12_residual/`, and struck deuterons 4.3e-8, both largest at the highest beam setting.

How strong the electron push would have to be for the patch to reach threshold, read from the linear part of the model: 4,000 to 2.5 million units, depending on the set of rates and the beam setting.

## What the result rests on

1. **How hot electrons act on a deuteron in PdD.** Both the directed force and, separately, the direct change of site stiffness when the electrons around a site are hot. The second acts on every site at once, symmetric or not, and is the quantity the pump threshold is written in; here it is carried only as a declared setting. Neither has been measured or computed for PdD.
2. **Passage of motion between neighbouring parts of the patch and inside a vacancy cage.** The patch is built from eight-site blocks joined by the linear lattice response; the cage uses the bulk couplings.
3. **From the closest computed pair distance to contact, and how many vacancy cages a foil holds.** The pair energies end at 2.65 Å in bulk and 1.85 Å in the cage.

The lattice numbers are those of pass 1 (the earlier `g3_lattice/coherence/` inputs, not the corrected rerun).

## Lead's note

Archive SHA-256 f960b918ddb322f86915a8beb89c7d66662d464bbcc96ba655a83830a1ef066d, all 34 file hashes in `manifest.json` checked, no machine paths. Three tables are left out of the repository for size: `pair_consumers.csv` (48 MB), `yield_consumers.csv` and `whole_patch_modes.csv` (7 MB each). The evidence archive (site-level arrays and trajectories, 79 files) has SHA-256 21ed2f9807d40d54f873f53bd20a305d9b85b159961682c1c614980193867535 and is not in the repository.

The lead reran `run.py --consume` with NumPy, SciPy and Matplotlib on a clean copy of the needed directories at commit c037c8c plus both archives. This rebuilds every table from the saved trajectories; all agree with the delivery to 1e-9 or better. The 72 trajectory runs themselves (`--source-only`) were not rerun.
