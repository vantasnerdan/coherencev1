# Pass5: one physical yield merit and paper-axis comparison

The primary denominator is the energy- and run-bound **transported PdD profile
plus its separate adsorbate**, with the three causes refitted on all four foils.
Uniform pass4 targets are retained only as named closure references. The original
source is `853c37738caca76ba70cefd8f4df7350bc204937`; the archive also contains the
necessary shared `forward/consumer.py` refinement and its README update.

## API and units

```python
from g4_beam.merit.api import Target, merit, compare
candidate = Target('design', profiles, material='PdD', stopping_fn=stopping)
factors = merit(candidate, [.25, 1., 2.5, 6.5])  # dict: three cause arrays
rows = compare(candidate)                     # canonical long CSV-compatible rows
```

`profiles(E_cm_keV)` or an exact-energy mapping returns `(layers,nads_m2,state)`.
Each layer is `(zt_m,zb_m,n_host_m3,n_D_m3,Z,M_amu,n_metal_m3)`; the host inventory
excludes the separately counted D, but includes oxygen in oxides. The stopping
callback receives **lab keV** and the same SI layer. `columns=((area,profiles),...)`
models parallel island-covered/open columns; areas sum to one. `measurement`
contains explicit beam/current settings. `Ue_shift_keV` is a scalar or
`(E_cm,z_m)` callable added after the fitted cause. Treat a Target and its profile
snapshots as immutable after evaluation; the object caches numerical results.

`evaluate_fit(target,energies,cause,fit)` supplies raw physical/measurement yields
at one explicit transported-fit row. It never chooses a denominator; later
profile-likelihood consumers must evaluate numerator and reference with the
same fitted parameters. Its default-fit identity is exercised in `explicit_fit.log`.

`merit` is reactions/true charged coulomb of the design divided by the modeled
PdD value at the **same CM energy**, with the same cause and beam. It is not a
ratio of two self-normalized curves. `cause=` selects `screeningrise`, `surface`
or `neutrals`. Fixed extracted neutral flux and energy still react in the changed
target and therefore do not imply neutral merit one.

`compare` evaluates the14 paper energies0.25–6.5keV. Its columns distinguish:

- physical proton reactions/true charged C, ion and neutral contributions;
- charged and all-incident-D counts/C and the corresponding reaction probabilities;
- same-energy merit versus PdD and Pd, and same-target bare-BH enhancement;
- the **primary fixed-PdD6.5 curve**, physical self-R and stitched reported-R;
- literal Eq7 effective Ue at0/1/2 additional inverse-current powers, its log RMS,
  fit scope and reference-specific matched Ue ratios;
- actual SI1 counts/cup C and same-energy material ratios, model-bridged anchor
  markers, observed self-R, separately quoted paper Ue, and matched paper refits;
- cause fit chi2/dof/current reading, so poor reference fits remain visible.

The effective fit uses2.5<=E<6.5keV, with the6.5 reference fixed. Design Ue uses
PdD's observed energies/error weights; ratios to Pd use Pd's matched weights.
`self_R_current_power_j` names the self-normalized **reported** curve after the
same extra current correction applied to the digitized paper curve; `self_R`
alone is the physical proton-yield ratio. Quoted paper Ue is0.8/1.7/0.3/0.9keV
for Pd/PdD/Ti/TiD, each+-0.2keV, not our new refit values.

SI1 supplies absolute cross-material **same-energy counts ratios** only at0.25–3
keV. To place its markers on the primary modeled scale, compare multiplies each
ratio by modeled PdD(E)/PdD6.5 and writes the model bridge explicitly. Experimental
anchor/count cells at3.5–6.5keV are null. No measured0.45 PdD/Pd6.5 anchor exists.

## Kept physics and decisive omissions

All depth integration uses the shared consumer: molecular85/15 electrical-current
fractions, Gaussian spread, target-dependent stopping, explicit adsorbate, fixed
common Eq17 proton reaction background and BH neutrons. The fit accelerator
compiles the same consumer stopping nodes; final fit chi2 is evaluated by the
actual consumer, not interpolation alone. No free per-material normalization is
introduced. Transported references use actual2h/12h source histories; only missing
energies are newly generated, using the supplied BCA cache.

The three largest omissions are lateral film morphology/ballistic mixing,
qualified thin-film retention and surface reaction kinetics, and independently
resolved neutral/charge distributions. Within-snapshot evolution, straggling and
backscatter are dropped. Separate supplied snapshots carry time and erosion.
These are transport-model profiles, not measured D-depth reconstructions.

## Outputs and reproduction

`refit_summary.csv` gives all three causes at all three readings beside pass4.
`transported_material_ratios.csv` gives predicted PdD/Pd, TiD/Ti and TiD/PdD
beside actual SI1 ratios. `compare.csv`, `sensitivities.csv`,
`density_normalization.csv`, `sealed_reference_spread.csv` and
`pass4_ratio_closure.csv` are the numerical consumers; `merit.png` replots them.
`verification.json` and `path_parity.json` record actual finite numerical evidence.
All input meanings/status/source/consumer links are in `inputs_ledger.csv`.

From the repository root, with NumPy/SciPy/Matplotlib:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python g4_beam/merit/prepare.py
python g4_beam/merit/refit.py
python g4_beam/merit/refit.py --surface-alternative
python g4_beam/merit/report.py
# Replot existing outputs without transport or refitting:
python g4_beam/merit/report.py --plot-only
```

Existing supplied profiles are reused; no BCA rebuild, DFT, old pass3 replay or
embargoed Zenodo access occurs. The direct old/new consumer pilot agrees within
6.29e-6 relative. Uniform global D scaling with deliberately frozen stopping must
multiply absolute yield while canceling physical self-R; physical D stopping and
localized top5nm changes are reported separately, not conflated with that identity.

The stable target ID `ideal_sealed_0.7D_per_Pd` in this archive is a **density
control retaining the Pd adsorbed sheet**, plotted as “0.7D/Pd + retained
adsorbate”. It is not the actual ALD seal. `causes/` separately evaluates losing
that sheet and recomputes the headline at the density-profile fit rather than
dividing an old merit by a density multiplier.

Completed refit minima are51.4793 (screening rise),69.8902 (neutrals) and262.2449
(surface, poor), beside pass4's36.9248/46.3120/40.2913. The distinct5nm surface
start gives262.3363 and is not retained. The controlled global-density self-R
identity agrees to6.66e-16. An additional0.515/0.526/0.783/2.819keV gives10x at
0.25/1/2.5/6.5keV in the best screening-rise case. Actual SI1 material-ratio
disagreement remains visible; no per-material amplitude hides it.
