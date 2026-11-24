# g12_residual: the gap taken apart from the data (lead's summary)

Earlier nodes started from a mechanism and built toward the measured curves. This node starts from the difference between the measurements and the model and asks what that difference depends on. `data/` (giuliano) is the table of published measurements with their reported conditions. `regression/` (Marvin) is the analysis. `inputs_ledger.md` lists every input of the model side, where it comes from, and whether a later node replaced it.

## The model side was cleaned first

The model the measurements are compared with now contains nothing fitted to the 48 UC Davis points: the transported deuterium profiles, standard stopping, the standard smooth reaction rate, a common 40 eV of screening for all four foils, and the paper's own normalisation at 6.5 keV. Several earlier results carried pieces fitted to the same data (a factor joining the proton and neutron points, a knee chosen from the data, a threshold-resonance width of about 28 eV where the externally measured value is 40 meV). They are listed at the end of `inputs_ledger.md` and are not used here.

## What the difference is

**It is an added signal, not a multiplier.** Between 0.25 and 2.5 keV the measured yield changes by a factor of 2.6 to 6, depending on the foil. Over the same range the model falls by a factor of about 1e16 to 1e17. Dividing one by the other gives a steep curve whose shape is the model's own collapse. Subtracting gives the informative quantity: a low-energy signal that barely depends on beam energy.

| Foil | Added signal below 2.5 keV, as a fraction of that foil's own yield at 6.5 keV |
|---|---:|
| Pd | 0.0044 |
| PdD | 0.0113 |
| Ti | 0.0010 |
| TiD | 0.0019 |

- **No slope with energy.** No foil's slope differs from zero by two of its formal error widths, and the foils do not share a sign.
- **Loading raises it.** Point by point, PdD is 2.2 to 4.6 times Pd and TiD is 0.8 to 4.9 times Ti. Each foil is in its own units, because the paper gives no common absolute scale, so these ratios are indicative.

**The counts behind it are few.** The paper gives a background of 0.007 counts per hour in the proton window and signal-to-background ratios down to 15 at the lowest energies, with runs of 12 to 24 hours. That implies about 1 to 4 detected protons at the lowest palladium and titanium points and up to 30 to 60 at the best TiD point, against less than 0.2 expected background counts per run. The flat level is therefore real at the level of "counts were seen where none were expected", and loose at the level of its exact value.

## What it is constant in: not settled by the published data

The added signal was expressed per coulomb, per unit time, per incident deuteron and per unit of deposited beam energy. Which of these is flattest depends on how the published numbers are read: per unit time under one reading of the main figure, per coulomb under the other reading and under the paper's supplementary count table. The paper gives one current-against-voltage curve for the gun, not the current, charge or duration of each point, so inside this data set current, run time and beam energy all move together.

Short symbolic formulas fitted to the added signal pick out the host metal in every case and little else, and do not predict a held-out foil better than a plain constant does. With 24 low-energy points and two host metals there is no functional law to read off.

## Other laboratories

The same split was applied to each published data set in its own units. None of them is the same kind of measurement as UC Davis, so they are compared one at a time and not pooled.

- Titanium between 1.3 and 2.4 keV in one earlier experiment (Yuki 1997) follows the ordinary barrier shape, with no flat added signal.
- Palladium and oxide-faced palladium foils (Kasagi) show enhanced yields whose added part varies strongly with energy.
- Zirconium and titanium at the lowest energies in a recent detected-yield series (Czerski 2024) fall far more gently than the model, in the direction of an added signal.
- A loading experiment at 15 keV (Chen 2025) shows the rate rising 11 to 19 percent with loading.

## What this narrows

Stop searching in anything that changes the energy or the barrier of the collision between a beam ion and a target deuteron. That covers screening in all its forms, the motion of a single target deuteron in its well, the arrival transient and the threshold resonance at its measured strength: each would have to grow by sixteen orders as the energy falls and still leave a flat result. The signal behaves like a second source of reactions that the beam energy does not control.

Not covered by any of this: many deuterons moving together away from equilibrium. No node has built it. A source set by the state of the lattice and not by the beam energy is the kind of thing the flat signal asks for, and it is taken up in `g11_cmns` (collective part).

Still open, and not separable with what is published: whether that source scales with beam delivered or with time elapsed, and how much of the loaded-foil increase is deuterium content and how much is surface state or history.

## The measurement that would separate them

On one prepared foil at one fixed energy: run at several beam currents with the same total charge, and at one current with different total charge, recording collected charge and live time independently. A signal that goes with charge gives a rate proportional to current. A signal that goes with time does not. A dependence on history shows up between the equal-charge and equal-current arms. The per-run records due on 2026-12-01 may already contain enough to answer part of this.

## Lead's note

Regression archive SHA-256 f73abf37021ef0b01b95b2f7b0d8a177665436aee4e4460e17867394a5fdb34c; all 95 manifest hashes check. The lead reran the portable replay (`analyze.py --joined-input regression/analysis_rows.csv`): the residual, added-signal, loaded-ratio and model-shape tables are byte-identical and the two summary tables agree to 4 parts in 1e15. The symbolic searches were not rerun (PySR is not installed on the lead's machine); their formulas are the delivered ones. The count estimates above are the lead's arithmetic from figures in the paper's main text.
