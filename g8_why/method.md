# Plain-foil mechanism comparison, fixed transported profiles

The governing outcome is three sequential mechanism pages: reaction, driven target,
then beam. Nothing here changes established consumers or the saved deuterium
profiles. Collision work is separate.

`common.py` reads the exact-energy records through
`g4_beam.merit.api.reference`, including each separately stored native adsorbed
sheet. The same profiles, stopping law, 0.1-keV laboratory cutoff and 12-point
CSDA reaction quadrature serve every hypothesis. The baseline is smooth
Bosch–Hale S(E), 40 eV, 85/15 D+/D2+ **electrical-current** fractions and the
retained 5% FWHM source-spread assumption. No profile, sheet amplitude, layer,
site population or stopping parameter is fitted. Changes of incident species or
energy alone require their corresponding incident CSDA path in that same saved
nominal-run profile.

## Meaning of the numbers

`baseline.csv` separates raw modeled proton/neutron reactions per **true charged
coulomb** from the published curve normalized to each foil's own 6.5-keV-CM
reference. The latter is not an absolute measured yield or pointwise S(E).
The baseline readout retains high-voltage cup-focusing correction
`c(E)=Icup/Itrue`, the fixed earlier nuisance stitch `s=1.0413535612834366`, and
freezes p/n calibration at the smooth-baseline 4.5-keV value. This calibration is
an assumption, not an independent detector measurement. The primary comparison
uses the figure as plotted (`k=0`); `k=1` is the previously discussed additional
inverse-current reading. Both are saved, not combined as physical mechanisms.
A hypothesis-dependent 4.5-keV re-splice is explicitly labeled as a sensitivity,
never allowed to hide a primary predicted p/n branching change.

For a candidate channel factor `F_ch(E)=Ycandidate_ch(E)/Ybaseline_ch(E)`, the
fixed-calibration shape enhancement at proton energies is
`F_p(E)/F_n(6.5)`. Thus an observed normalized excess of 1 at the reference does
**not** force the physical factor there to 1. No new raw counts or run likelihoods
are reconstructed. The earlier Fig. 2/SI1 current-normalization discrepancy is
not evidence by itself of a duplicated current division.

### Absolute floor semantics

`readout_floats.csv` supplies all five energies and both current readings.
The following vectors are **conditional model requirements** in reactions/true C,
anchored to the smooth40 neutron yield at 6.5 CM, not absolute experimental rates.
Entries are ordered [.25,.5,1,2.5,6.5] CM keV; the last entry is neutron-channel.

| Foil | k0 reactions/true C | k1 reactions/true C |
|---|---|---|
| Pd | [4898,1065,1219,2082,511049] | [782,304,563,1207,511049] |
| PdD | [21604,7133,3840,9780,738682] | [3449,2038,1775,5671,738682] |
| Ti | [52355,48328,29265,57457,46351117] | [8359,13808,13526,33315,46351117] |
| TiD | [123659,40770,30699,180358,46611787] | [19744,11648,14188,104577,46611787] |

Actual SI1 **detected counts/cup C** at [.25,.5,1,2.5] CM are
Pd [1.075,.429,.778,1.705], PdD [1.535,.915,.760,2.494],
Ti [.768,1.344,1.220,3.170], TiD [2.865,1.773,1.980,15.246].
The supplied source is `S:B*0.007 counts/hour /(I_cup[A]*3600)`, not an
invented count history. Physical emitted proton reactions/true C are these
values times `f_col/eta_p`, with `eta_p` the unknown detector acceptance and
response, not total detector efficiency set to1. SI1 has no 6.5 datum.
The central values neither establish an energy-independent floor per C
nor a common floor across all four foils; uncertainties/run joins are unresolved.
The conditional model anchor and empirical detected rate must not be equated.

### Beam ratio parameters

`beam_neutral_requirements.csv` uses `neutral_rate_convention=0` for a neutral
D/true charged-unit ratio at every setting. With convention1, its solved ratio
columns and `beam_tail_fits.json.fitted_neutral_D_per_charge` instead name the
ratio **at I_INT**; local neutral D/charge is that parameter times
`I_INT/Itrue(E)`. They are not fractions of the total beam and can exceed1
without a mathematical contradiction; source/current/power constraints are
separate. Case identifiers containing `f3e4` in `beam_floats.csv` mean the
explicit test parameter **3e-4**, as fixed by `ftest` in `beam.py`, not3e4.
The13-keV test tail is a5%-FWHM Gaussian, not a measured Davis neutral spectrum.



## Reproduction

From the checkout root with NumPy/SciPy, use one thread:

```sh
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python g8_why/reaction.py
python g8_why/readout_floats.py
python g8_why/driven_target.py
python g8_why/beam.py
```

Each branch script saves its own compact floats and a targeted proof of the
number it produces. This is finite model computation, not a formal numerical
enclosure, a transport rerun or a new first-principles calculation. The scripts
have no Native/DFT or external-data dependencies. Source hashes accompany the
branch evidence. This directory stops at the requested mechanisms; it does not
start a new campaign.

## Delivery order and evidence

1. `reaction.md`: `reaction.py`, `reaction_source_notes.md`, `baseline.csv`,
   `reaction_floats.csv`, `reaction_retained.csv`, `reaction_shape_bound.csv`,
   `reaction_results.json`, `reaction_proof.json`, and unchanged `common.py`.
2. `driven_target.md`: `driven_target.py`, `driven_target_floats.csv`,
   `driven_target_budget.csv`, `driven_target_proof.json`; the shared absolute
   reporting pair is `readout_floats.py` / `readout_floats.csv`.
3. `beam.md`: `beam.py`, `beam_floats.csv`, `beam_neutral_requirements.csv`,
   `beam_wake.csv`, `beam_timing.csv`, `beam_optical_depth.csv`,
   `beam_tail_fits.json`, `beam_proof.json`.

Each page was completed and sent for recording before starting the next branch.
The original common helper and baseline remain unchanged after reaction.
All three targeted numerical consumers exited0. No broad suite was run.
The retained threshold law and ordinary driven heating fail the requested shape;
a genuinely energetic beam tail remains conditional on an unmeasured energy path.
These outcomes end this assignment, not initiate another mechanism campaign.
