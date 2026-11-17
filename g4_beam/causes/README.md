# Pass5: discriminate the plateau causes and infer the missing surface density

This node consumes the accepted `merit/` archive (integrated as98efac2), without
replaying or changing its canonical m=1 reference. Its three cases are the actual
transported-profile refits: screening rise, surface layer and fixed-energy fast
neutrals. The canonical nearly full-sheet surface fit is poor (chi2=262.2449).
Removing the assumed Pd/PdD sheets and refitting instead gives surface chi2=44.2060:
that explicitly named **refitted_zero_sheet_surface** scenario must accompany the
canonical reference, not silently replace it. Cause ranking therefore depends on
unmeasured sheet occupation. No thermal-spike cause, free constant neutral floor,
new DFT or pre-embargo Zenodo access is introduced.

## The six responses

`predictions.csv` is the canonical compare() output for all14 paper energies and
all three causes, augmented with changed/base yields per incident charged C and
rates. It contains current, D2+ electrical-current fraction, incidence angle,
time into a run, foil temperature and near-surface D factors0.1/3/50.

- Current has both an immediate frozen-target response and a transported,
  beam-heating-adjusted response. Extracted neutral flux is held fixed while
  delivered charged current changes. Rates and yields/C are separate columns.
- Molecular fraction is an immediate beam switch: D2+ current supplies two D at
  half the atomic-D arrival energy. No long-time implantation distribution is
  invented for a composition not in the supplied BCA cache.
- Angle is an immediate straight-ray rotation. Depth is converted to distance
  along the ray, adsorbate optical thickness scales by sec(angle), and the
  fitted surface screening boundary retains its normal-depth thickness.
  Backscatter/reflection and detector-acceptance changes are not modeled.
- Temperature and time use actual retained-parameter Model histories, including
  the same3600s preloading and energy-specific2h/12h nominal runs.
- The density perturbation scales volumetric D only in the top10nm; the already
  nearly full outer sheet stays fixed. Host inventory is unchanged, D stopping
  is recomputed and no concentration clipping hides large diagnostic multipliers.

`laboratory_ratios.csv` selects three distinct knobs with the largest predicted
cause separation and states the measurable numerator/denominator. These are
proposed counts-per-independently-measured-charge ratios at fixed energy, not
new observations. Required beam, target, time, charge, detector and covariance
columns are in `zenodo_columns.csv`. No claim is made that those columns already
exist in the embargoed files.

## Density likelihood, not a new reference amplitude

The shared1D grid has13 multipliers0.1–100, applied to volumetric D over the top
10nm in BOTH Pd and PdD; Ti/TiD and all sheets are unchanged. At each point the
original cause/nuisance parameters are reoptimized at that cause's best current
reading. `density1d.csv` gives the profiled shape chi2 and the secondary frozen-m1
parameter response. The m=1 fit is reused. Each cause also has one m=1/sheet-zero
point, removing the sheet in both Pd foils; its fixed-parameter yield loss and
profiled fit are distinct. `density1d_summary.json` gives the minimum and sampled
ranges within delta chi2<=1 and4. These are finite-grid ranges, not interpolated
continuous confidence enclosures. Existing top5nm response values are reused in
`top5nm_reused_response.csv`; no second5nm profile-fit campaign is run.

The rising-screening2D grid has7×7 independent Pd/PdD multipliers0.3–30. Every
cell refits against the48 shape points plus the7 actual SI1 PdD/Pd ratios. For
that joint score the ratio log-error is the approved proxy
`sqrt(sigma_lnR_Pd(E)^2 + sigma_lnR_PdD(E)^2)` from the same-energy Fig2 bars,
assuming independent normalization errors; actual count/reference covariance is
not available and must replace this proxy when the per-run data open. Shape and
ratio components remain separate, so the result does not depend on concealing
that weighting choice.

`density2d_maps.png` shows, on the same grid, shape chi2 at the joint optimum and
the unweighted mean of the7 log(model/SI1) residuals. White outlined cells mark
shape chi2 within4 of the grid's shape minimum AND **all seven** modeled
same-energy ratios inside the exact observed min/max band0.97727–2.13333.
`density2d_region.csv` gives that conservative sampled region, including the
implied top10nm D/Pd window ratios at four energies and the maximum local ratio.
The geometric-mean-only band criterion is retained as an explicit secondary
column, not silently substituted for the primary region. The joint minimum is
reported alongside the two maps.

`inferred_baseline_sealed.csv` recomputes the physical0.7D/Pd numerator and the
inferred PdD reference at the SAME fitted cause parameters. It reports retained-
sheet and sheet-removed controls separately and does not divide old merit by a
bare density multiplier. These are instantaneous density controls, not a claimed
charge-integrated ALD device. The actual run-integrated design headline belongs
to the resolved erosion trajectories in `ald/`.

`material_ratios.csv` puts the original refit, each shared-density minimum and the
independent joint minimum beside actual SI1 PdD/Pd, TiD/Ti and TiD/PdD ratios.
There is no0.45 PdD/Pd6.5keV measurement or invented high-energy absolute anchor.

## Occupancy diagnostic and the second named reference

`sheet.py` profiles all three causes on a6×6 **linear** grid of independent
absolute Pd/PdD coverages0,.2,.4,.6,.8,1ML, with1ML=1.53e19D/m². Original
transported volumetric m=1 profiles stay fixed; Ti/TiD remain unchanged. Each cell
uses the same55-component joint likelihood and its cause's best charge reading.
This fixed-bulk/new-sheet inconsistency is an explicitly authorized diagnostic,
not a newly solved transport state. No further bulk grid is inferred from it.
`sheet_grid.csv`, `sheet_region.csv` and `sheet_summary.json` retain both score
components, the exact seven ratios, joint minima, accepted finite cells and
whether Pd is fuller than PdD. Maps use the same all-seven-ratios region rule
as the density maps, separately within each cause.

If a sampled occupancy region fits shape and ratio, one subsequent transport
boundary calculation must join bulk and sheet using actual boundary flux/kinetics;
patching an exported nAds value is not that calculation. Inferred parameters
remain inferred, and no barrier or material-specific sticking law is invented.

`references.py` implements the separate44.2060 reference without mutating the
closed merit API or a global cache. BOTH design and reference use that same
surface fit. It removes adsorbed D only on actually exposed Pd, independently in
each parallel area column; bulk D, film D, stopping, geometry, current and
material identity do not change. Covered oxide/cap columns do not acquire a Pd
sheet. Uncovered island/promoter fractions lose their exposed Pd sheet, and the
actual outer nanolaminate material determines its identity. An explicitly assumed
oxide adsorbate is not called a Pd adsorbate or silently erased. The transformation
is diagnostic, not a self-consistent boundary rerun.

## Decisive completed results

- The shared bulk scans prefer their sampled boundaries: rising and neutrals
  m=.1 (shape chi2=45.5231 and62.6930); surface m=100 (99.7561). Their sampled
  delta<=1 ranges are .1–.17782794, .1–.17782794 and100 only; delta<=4 extends
  the first two to .56234133. These are conditional finite ranges, not measured
  density confidence intervals.
- The independent rising grid has its joint minimum at(.3,.3):
  shape58.99394 + ratio15.24990 =74.24383. Its seven ratios2.320–3.410 remain too
  high. **No sampled cell** satisfies the requested shape-and-ratio region,
  even with the explicitly secondary geometric-mean band rule. More bulk D is
  not supported by this controlled scan.
- Removing sheets flips the surface ranking (262.2449→44.2060), but that44.2060
  fit uses shape alone; it does not establish cross-material counts agreement.
  For example, the separately JOINT-refitted zero/zero occupancy pilot gives
  shape88.04123 + ratio37.52351 =125.56474 and ratios4.339–6.367.
- All six knobs completed in393.34s,2workers/threads1:70 Model cases and714
  compare rows. The three distinct-knob laboratory ratios are density×50
  at.5keV (rising/surface/neutrals14.5407/1.01390/13.8010), D2 electrical fraction
  .85 at.25keV (3.07592/3.30731/1.00047), and half current with actual history
  at.5keV (1.00218/1.00016/2.00034). These are predicted changed/base physical
  proton counts per independently measured incident charged C, not measurements.
- The completed108-cell occupancy grid also has **no accepted cell for any
  cause**, including the geometric-mean alternative. Joint minima are rising
  (thetaPd,thetaPdD)=(1,.6), shape63.87299+ratio12.96737=76.84036; surface(0,.2),
  86.83761+38.23750=125.07511; neutrals(1,.4),86.84841+13.24812=100.09653.
  Rising/neutrals favor Pd fuller than PdD, but this does not solve the counts
  constraint. No additional grid or joined-boundary rerun is triggered.
- **Remaining gap:** no sampled target configuration joins shape and SI1 counts.
  The canonical comparison has model/SI1 excess2.55–6.83; the closer joint
  cells' geometric-mean ratios are1.91,2.69 and2.02 (roughly1.5–3×).
  Some other cells improve the ratio only by spoiling
  the shape. Thus this is not a claim that every individual discarded grid cell
  overpredicts every energy. The remaining candidate to check is the charge
  measurement: true charge per run and cup collection in each beam regime,
  assigned separately to g4_zenodo, not established here as the cause.
- Recomputing the physical sealed0.7 numerator AND inferred(.3,.3) denominator
  at the same rising joint fit gives no-sheet instantaneous merits
  136.5458/40.9755/22.2625/31.9153 at.25/1/2.5/6.5keV. The accepted-region
  spread is **null**, not an invented interval, because the region is empty.
  Under the separate refitted zero-sheet surface reference the corresponding
  density-control merits are52.5832/47.1932/45.1244/46.9447. Neither set is the
  charge-integrated ALD-device headline.
- The settled report exited0 in60.88s. `verification.json` records exact
  zero-sheet selfreference1 and area-composition error0 using supplied PdD/TiD
  columns as an algebraic witness, not an extra candidate. It confirms the
  original Pd sheet is unmodified and the non-Pd adsorbate is preserved.
  Maximum compiled/direct occupancy-grid chi2 difference is0.00110069; the
  final tables use the direct consumer values. All108 solver flags are success.

## Kept physics, omissions and reproducibility

Both the compiled fitting response and final raw yields use the accepted
consumer's stopping nodes and fixed common Eq17 proton/BH neutron reactions.
Final profile/joint scores are reevaluated by the actual consumer. No fitted
per-material amplitude is introduced. Input provenance is in `inputs_ledger.csv`.
The three principal omissions are resolved neutral/charge/angular distributions,
qualified evolving surface chemistry and thin-film retention, and the real
per-run count/normalization covariance. Transport retains its source IV-as-
incident-flux clock; true-charge and cup normalization are explicit in yields.

From the repository root with NumPy/SciPy/Matplotlib:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python g4_beam/causes/density.py --pilot-joint
G5_DENSITY_WORKERS=4 python g4_beam/causes/density.py
G5_RESPONSE_WORKERS=2 python g4_beam/causes/responses.py
python g4_beam/causes/sheet.py --pilot-sheet
G5_SHEET_WORKERS=4 python g4_beam/causes/sheet.py
python g4_beam/causes/report.py
python g4_beam/causes/report.py --plot-only  # no fitting or transport
```

Completed points are reused on restart; no valid fit is replayed for packaging.
`peak_RSS_KiB` is Linux lifetime `getrusage` highwater and may include an inherited
pre-exec image; it is not silently interpreted as current-job private allocation.
Explicit numerical threads are1; aggregate allocated heavy workers never exceed6.

Observed runtime: Python3.12.3, NumPy2.5.3, SciPy1.18.1, Matplotlib3.11.2.
The host reports assistant model `openai-codex/gpt-6.1-sol`; that identity is not
scientific evidence. Actual process exits/costs and output hashes are in
`commands.csv`; `SHA256SUMS` binds the frozen payload. Commands and gzip/floating
outputs may vary at the byte level across runtimes; hashes identify these actual
source-run bytes, not a promise of bit-identical regeneration.

The archive requires the original853c377 repository plus the accepted merit
archive (integrated98efac2); it does not overwrite the closed merit defaults.
The three `*_fits.json` tables are in the repository. `knob_profiles.json` (2.8 MB) is
kept only in the delivery asset `g5_causes.tar.gz` (draft release `g5_causes-delivery`);
`responses.py` rebuilds it. Its SHA-256 is listed in `SHA256SUMS`.
Figures replot from included finite tables without fitting or transport.
