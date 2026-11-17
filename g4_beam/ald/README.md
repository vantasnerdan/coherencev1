# Pass5 ALD designs: physical yield, lifetime and charge-integrated merit

This node builds on source853c377 and the accepted merit/causes deliveries
(integrated98efac2/80f3e87). It does not refit those closed causes, access embargoed
Zenodo data, run DFT/BCA, or claim experimental ALD performance. Source/parameter
status is in `inputs_ledger.csv` and `inputs/*.json`.

## Decision and delivered scope

The best resolved conditional candidate is the20nm TiD2 cap with retained D and
slow sputtering (nominal sputter times .1). Under `screeningrise`, its raw
charge-integrated ratios at .25/1/2.5/6.5keV are
210.292/17.925/36.855/37.987. The measured-TiD-counts proxy gives
40.925/8.504/9.439/null instead; the neutral scenario gives18.542 at2.5keV.
Above3keV only the raw model is available. Faster sputtering or lost cap D
removes the strong gain: retained/fast raw ratios are7.326/1.227/3.782/4.155;
lost/slow ratios are .111/.704/.179/.346. This identifies a retention/lifetime
material requirement, **not demonstrated tenfold ALD-device performance**.
Tenfold reach is energy/cause/scenario-specific, not an all-four-energy gate.

Delivery covers245 designs,10290 canonical target/cause/energy rows, both named
reference scenarios, and76 resolved four-energy cap/seal designs including the
nominal TiD2 case. All requested remaining transport cells are resolved.
`verification.json` confirms complete0..run-end curves, exactly reproduced
serialized ratio-of-integrals, and laminate inventory error5.45e-15.

Alternate-snapshot coarsening changes the resolved integral by at most5.31%;
this is a resolution diagnostic, not a certified error enclosure. Two shape-only
zero-sheet seal cases at2.5keV straddle10: Al2O3 .25nm/75% coverage gives
10.0679 versus9.8647 when coarsened; full Er2O3 .25nm gives10.1956 versus9.9095.
Treat those marginal tenfold classifications as resolution-sensitive. This does
not change the selected TiD2 material requirement or promote the zero-sheet
scenario to a joined physical boundary.

## Resolved cap decision

For retained20nm TiD2 at nominal sputter factor1, the charge-integrated results
are below. Each low-energy cell is **raw / SI1-TiD-counts-calibrated proxy**;
6.5keV has no measured anchor and shows raw model only.

| Cause | .25keV | 1keV | 2.5keV | 6.5keV raw |
|---|---:|---:|---:|---:|
| Rising screening |64.125 /12.480|3.270 /1.552|27.348 /7.004|28.495|
| Canonical surface |.364 /.056|.833 /.309|2.799 /1.617|15.772|
| Fast neutrals |24.781 /2.958|5.171 /1.447|27.590 /13.556|27.118|

Thus a nominal retained cap reaches the calibrated10x goal at.25keV in the
rising case and2.5keV in the neutral case, NOT under every cause/energy.
Retention and lifetime decide the reach: fast sputter×10 reduces canonical
rising raw merit to7.326/1.227/3.782/4.155; lost-D caps can attenuate PdD below1.
These are modeled conditions, not demonstrated ALD durability.

The calibrated largest verified sputter factors (next sampled failing factor)
are rising1.2287(1.3161) at.25keV and.03487(.03760) at1keV; the noneroding2.5keV
limit is only9.621. Neutral1keV requires.02598(.02792), and2.5keV allows
1.4096(1.5098); the noneroding.25keV neutral limit is7.450. At6.5keV, RAW-only
limits are3.4872(3.7598),1.8552(1.9870),3.2345(3.4872) for rising/surface/neutral.
The shape-only zero-sheet variant has its separate factors/results in the CSVs;
its stronger gains do not establish an eligible physical zero-sheet boundary.

Lifetime columns are distinct: `lifetime_nominal_h` and its low/high companions
state the source central inventory/sputter proxy and bracket. Cap
`scenario_sputter_lifetime_h` applies the selected sputter multiplier. Oxide,
promoter and BZY-front `model_soft_clearance_h` instead uses the actual
erosion-only law through the1pm cutoff; it is not the shorter linear inventory
time. These remain model/proxy clocks, not measured durability confidence bounds.
For the retained law `d_dot=-v*d/(d+lambda)`, that column is
`t_nominal*(d0-d_cut+lambda*log(d0/d_cut))/d0`, with `lambda=.02nm` and
`d_cut=.001nm`. It is an analytic model clock, distinct from a sampled numerical
clearance time and from a source uncertainty bracket.

## Completed designs3–6

All238 required transport cases are complete (42 surface-supply,84 promoter,
112 BZY), with0 unresolved cells. They feed112 design reports,4704 canonical
comparison rows and3840 erosion-comparison rows. The final checkpointed
report/erosion command exited0 in2019.10s without replaying those transport cases.

The following are **end-run** maxima in each family under `screeningrise`,
not run-integrated gains and not necessarily one common design across energies:

| Family | .25keV |1keV |2.5keV |6.5keV |
|---|---:|---:|---:|---:|
| D-bearing oxide |1.0000|1.0000|1.00598|1.00623|
| BaO/SrO/La2O3 promoter |1.0000|1.0000|1.00040|1.00073|
| BZY back/front endpoints |121.718|1.15659|19.2093|29.0700|
| Ideal loaded Pd/oxide laminate |158.815|46.8034|24.4693|35.5015|

Oxides/promoters returning close to1 at run-end have usually cleared; their
earlier effects remain in `remaining_erosion.csv`, not an invented zero response.
BZY's strong front-side cases require retained ideal hydration; its back-side
supply endpoints do not establish a room-temperature conductivity law. The
laminate gains require ideal .7D/Pd and the explicitly counted interface inventory,
with unproven ultrathin continuity. These endpoint results do **not** replace the
resolved cap/seal charge-integral headline. Other causes and the separately named
zero-sheet reference remain in `all_designs.csv`.

## Keeps, drops and three dominant omissions

The pre-code interface/choices are recorded in node entry1358, with the primary
reference amended in1374. This is a compact summary of the FINAL model, not a
claim that subsequent corrections were all present before coding.

**Keeps:** seven-field layered SI inventories and separate adsorbed sheets;
molecular beam/spread, host AND D stopping, fixed common nuclear reactions,
target-responsive neutrals, all three matched cause fits, actual Model histories,
area-weighted islands, host-inventory erosion, and separate physical charge
integrals. Every design is reported through the shared `compare()`.

**Drops:** straggling/backscatter, lateral coupling and morphology, a new
composition-specific implantation/BCA calculation, and invented microscopic
retention or promoter rate laws. It does not drop cap D, stopping or elapsed time.

The three dominant omissions are:
1. Three-dimensional closure/pinholes, percolation and ballistic mixing.
2. Measured thin-film D transport/retention, cap/Pd exchange and material-specific
   surface kinetics; frozen inventory endpoints are not a predicted outgassing law.
3. Independently resolved neutral/charge distributions and incident-charge/cup
   collection; fitted phenomenological laws and counts proxies do not settle them.

## Observable and references

The primary run result is

    integral I_true(t)*Y_design(t) dt / integral I_true(t)*Y_same-age_PdD(t) dt

with physical reactions per TRUE charged C on both sides, matching energy, age
and beam. It is NOT end-state merit, mean instantaneous merit, or a renormalized
self-R curve. Reported end merit, resolved peak and time merit>=10 accompany the
integral. The requested four energies are .25/1/2.5/6.5keV per-deuteron CM.
Canonical `compare()` tables additionally retain all14 paper energies, all three
causes, charge readings, self-R, literal Eq7 effective-Ue and actual SI1 anchors
(null above3keV). The primary overlay stays fixed modeled UC Davis PdD6.5keV.
No arbitrary material scale or measured6.5keV anchor is introduced.

References are explicit, not interchangeable:

- `canonical_full_sheet`: original transported m1 PdD with its actual source
  adsorbate, each of the three accepted transported cause fits. Best shape
  chi2: rising51.4793, surface262.2449, neutrals69.8902.
- `refitted_zero_sheet_surface`: original m1 bulk with exposed Pd adsorbate
  removed and the separately refitted surface parameters, shape chi2=44.2060.
  BOTH design and reference use that same fit. This is a shape-only diagnostic,
  not a joined transport solution or an explanation of SI1 count ratios.

The causes grids found no joint shape-and-counts region; no joined-boundary
rerun was triggered. The unresolved charge/run/cup-collection hypothesis belongs
to g4_zenodo. It is not silently absorbed into target density or a material scale.

`causes.references.without_exposed_pd_sheet` projects each parallel column:
remove only an actually exposed Pd adsorbate; leave bulk, film D, stopping,
geometry, material and beam unchanged. Full oxide seals/deuteride caps do not
retain a beam-side Pd sheet while present. Islands/promoters retain it only on
the actually uncovered/exposed Pd area. After clearance the actual Pd surface
returns. A separately assumed oxide adsorbate is not mislabeled as Pd.
`all_designs.csv` records exposed Pd area, retained Pd-sheet area/inventory and
outer-column identity. Outer nanolaminate material is stated separately.

## Measured TiD anchor and sputter requirement

UC Davis SI1 already reports TiD/PdD counts per collected cup charge ratios
of 6.11 at2.5keV and11.67 at3keV. These are measured counts/Cup-C ratios,
NOT independently calibrated incident-deuteron yields. The seven same-energy
anchors span.25–3keV; no anchor is extrapolated to3.5–6.5keV.

Cap rows retain their raw model merit and add `TiD_counts_calibrated_merit_PdD`:
raw merit divided by the cause-matched model/observed TiD/PdD ratio from
`merit/transported_material_ratios.csv`. This is a TiD-counts-calibrated PROXY,
not an observed cap performance. For non-Ti caps it additionally assumes transfer
of the TiD correction. The common PdD cause law used for design predictions is
retained; this empirical correction is not a new layer-specific screening fit.
`cap_TiD_calibration.csv` exposes every divisor and actual observed anchor.
The zero-sheet surface variant computes its OWN TiD/PdD ratio with that explicit
fit and sheet projection; it never borrows the canonical surface divisor.
Calibrating the yield does not rescale the assumed transport/erosion clock.

`cap_endpoint_spread.csv` places retained/lost-D and slow/fast-sputter integrated
results side by side. `nominal_cap_integrals.csv` adds retained20nm TiD2 at
sputter factor1. `cap_sputter_response.csv` and `cap_sputter_thresholds.csv` give
both raw and calibrated10x requirements from the SAME physical evaluations.
The calibrated threshold is the headline where SI1 covers the energy; above3keV
only the raw threshold is available. Factors multiply the input central sputter
proxy, not a measured TiD2 sputter yield. Rightmost sampled crossing brackets
have at most a10% factor ratio; the lower endpoint actually attains10x and the
upper endpoint fails. This is a finite numerical scan, not a global monotonicity
proof or a confidence interval. No successful sampled factor and no measured
anchor are separate statuses, not invented zero-valued thresholds.

The endpoint-model/SI1 divisor is held constant in age when applied to an
integrated run. It is not a measured time-dependent calibration, and does not
resolve the independently unknown charge collection or cap retention kinetics.

The noneroding-cap control checks the lower sputter limit. If it exceeds10 while
the positive grid does not, the positive bracket is extended below its initial
minimum instead of declaring the gain unattainable.

## Six requested families and their actual assumptions

1. **Al2O3/Er2O3 seals and islands.** Thickness .1/.25/.5/.75/1/1.5/2nm;
   permeability proxy3e11 times .01/100. At <=.5nm, projected covered fractions
   .25/.75/1 are separate parallel columns, not an assumed continuous monolayer.
   Actual Model transport/erosion is used. Er2O3 uses sourced dense-bulk formula
   density rather than inheriting Al2O3's density. Material permeability and
   ALD closure on Pd are unknown. The retained Model's oxide adsorption energy
   .2eV/Pd-site-density convention is a proxy, not a measured oxide monolayer.
2. **TiD2, ErD2, ScD2 and YD2 caps.** Thickness1/2/5/10/20nm; full/lost retained-D
   endpoints0/1 and host sputter scale .1/10. Each cap has actual host/D stopping,
   host inventory and erosion. Covered Pd adsorbate is removed, not moved onto
   the cap. Substrate follows the original matched PdD feed-state history;
   cap/metal exchange kinetics are not invented. Ti diffusion D0=3e-7m2/s,
   Ea=.60eV is a proxy, with .01/100 sensitivity for unsupported materials.
   These are conditional inventory endpoints, not measured retention bounds or
   a fitted outgassing law. Growth/loading/vacuum survival remain material-specific
   feasibility questions.
3. **D-bearing oxides.** Requested2–15at%D means D/(D+metal+O), not D/metal.
   For M2O3, D/M=2.5*f/(1-f). Full/lost retention endpoints reuse the actual seal
   transport and erosion, with changed film inventory AND stopping. D2O tracer
   incorporation is source-supported; beam-end retention is not measured.
4. **BaO/SrO/La2O3 promoters.** Coverage .1/.5 at assumed .3nm island height,
   one Pd-site/metal areal inventory. No work-function-to-barrier/screening law is
   guessed. Conditional1e-7/.5 loss endpoints apply exactly to the covered
   patch while promoter inventory remains; clean-Pd sticking returns at
   clearance. Uncovered Pd area is mixed once, not counted again through a
   thickness blend. Adsorption binding stays at the explicit Pd baseline.45eV,
   NOT a measured promoter energy. Neither quantity is a qualified
   material-specific barrier/rate. The explicit no-contact-resistance branch is
   used; actual film geometry still affects stopping and erosion.
5. **BZY20 back/front layers.** Thickness2/5/10/20nm. Backside supply endpoints
   blocked/baseline/double are not a predicted300K conductivity law. Backside
   chemical/electrochemical lifetime is unknown and not beam sputter lifetime.
   Front layers use dry/full ideal hydration0/.2D per formula; full hydration
   adds oxygen2.9->3.0. Front permeance and sputter are explicitly oxide proxies.
6. **Pd/oxide nanolaminates.** Pd1nm/oxide.2 or.5nm through150nm, unloaded or
   .7D/Pd plus1e19D/m2 at EACH interface, spread over .1nm Pd. The interface-only
   frozen inventory includes no additional outer adsorbate, even when the outer
   layer is Pd; this is not an equilibrium prediction for fresh Pd at1Pa.
   The zero-sheet scenario does not add a sheet to this design. Actual sequential
   layer stopping/erosion is evaluated. Continuity of1nm Pd, percolation, trap
   capacity, transport and retention are unproven; local D/metal is never clipped
   to hide an optimistic endpoint.

Pass4 Eads=.45eV, vacancy binding=.23eV, mtrap6 and clean sticking.5 remain marked
estimates until an eligible g5_recomb rate mapping exists. No converged new
barrier/rate has been inserted merely because a calculation or agent completed.
Fast neutrals remain fixed extracted flux/energy but physically react in and
lose energy through the changed target; there is no constant yield floor.

## Resolved trajectories and selection

`dense.py` pilots the full baseline/forward/I-O path. Selection is the UNION of
the canonical full-sheet original top10 peaks (PRIMARY), the previous joint-
reference top10, every coarse/reset10x disagreement, and all four cap endpoints
for each selected cap recipe; zero-sheet startup cannot crowd out lasting caps.
This gives75 resolved selected designs, plus the separately requested nominal
TiD2 row. Ties in the prior joint ranking use its coarse integrated score then
stable target name. Reset estimates are selection diagnostics, NOT bounds;
nonrerun integrated estimates are not published.
An original canonical peak is null when that candidate had no original trajectory;
it is not replaced by zero. Its separately cached joint-reference curve remains.

Coarse selection uses the32 actual saved Model ages at each energy; cap query
ages choose an available actual state, never an interpolated transport profile.
Only the selected ages extend the shared reference kernels. Completed kernels
checkpoint at minute scale; aligned coarse tails avoid redundant staggered grids.
An earlier universal-reference pilot timed out, retaining3279 valid four-fit age
kernels; none were discarded or replayed merely to package results.


Snapshots are approximately5min through actual clearance and >=1h afterward,
coarser to the run end, with additional near-clearance points. Oxide clearance
uses the actual soft erosion law d'= -v*d/(d+.02nm), down to the Model's1pm
boundary; cap clearance follows host inventory/sputter. Physical design/reference
Y are linearly interpolated separately. Time>=10 uses exact linear crossings of
Ydesign-10*Yref, not a sample count. Full0..run-end coverage is required; no absent
tail or0/0 reference is fabricated. A single shared age-reference history is not
132 target transport runs.

Model time retains its provided-IV-as-incident-flux convention. Conversion to
physical reactions uses an ASSUMED incident current from the fitted cup model;
both members share the same age/beam convention. We do not silently rescale histories when cup collection
changes. Lifetimes are source inventory/sputter proxies, not measured durability.
End-run depth CSVs retain all seven physical layer fields, area fractions and
separately counted adsorbed D; sheets are not smeared into volume D.

## Reproduce

From the repository root, with the accepted merit/causes data and the scientific
Python environment (Python3.12.3; NumPy2.5.3; SciPy1.18.1; Matplotlib3.11.2):

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
G5_WORKERS=4 python g4_beam/ald/designs.py --early --prepare-early
G5_WORKERS=4 python g4_beam/ald/designs.py --early
G5_DENSE_WORKERS=1 python g4_beam/ald/dense.py --pilot
G5_DENSE_WORKERS=1 python g4_beam/ald/dense.py --pilot-seal
G5_DENSE_WORKERS=4 python g4_beam/ald/dense.py --select
G5_DENSE_WORKERS=4 python g4_beam/ald/dense.py
G5_SPUTTER_WORKERS=1 python g4_beam/ald/cap_requirements.py --pilot
G5_SPUTTER_WORKERS=4 python g4_beam/ald/cap_requirements.py
G5_WORKERS=4 python g4_beam/ald/designs.py --remaining
G5_FINAL_WORKERS=4 python g4_beam/ald/finalize.py
python g4_beam/ald/plot.py
```

The accepted chronology is merit -> early1/2 records -> causes -> resolved ALD
integrals -> remaining3–6 -> final ALD archive. Cache checkpoints preserve actual
completed physical cases; start/end/RHS trial-progress logs diagnose a slow case
without changing equations. Trial RHS states are not accepted solver snapshots.
`--remaining-probe STAGE JSON_JOB` runs and saves one bounded selected stage case
before any full retry. A prior uninstrumented remaining run timed out at3600s;
it is not reported as completed science or evidence of a particular root cause.

Single-age raw forward evaluation explicitly opts into
`reference_E_cm_keV=None`: no unrelated6.5keV profile or normalized R is requested.
The default normalized API is unchanged. `raw_consumer_proof.json` records exact
physical-field parity for all3 causes/four energies, a valid single-energy profile,
raw physical zero, and retained normalized-zero rejection. This is a consumer
correction, not a new nuclear model, invented normalizer or fit replay.

Actual exits, costs and payload hashes are recorded with the delivery. Hashes
identify these actual source-run bytes; gzip times and floating-point execution
can vary across environments. Scientific gzip checkpoints are collaborator-only
draft-asset data; they do not authorise publication of private operational memory
or credentials. Archive publication is performed separately; no branch push or
merge is performed here.
The archive hash manifest covers the delivered payload. Diagnostic-log copies
redact operator-specific filesystem prefixes only; scientific source, numerical
data and accepted profile arrays are unchanged. Raw local diagnostic logs remain
available separately. Private plans, credentials, session history and disposable
Python bytecode are not delivery members. The nonrerun integral-estimate cache is
excluded; coarse/reset integral checks are retained only for the75 resolved
selection cases, plus the separately reported nominal cap calculation.

### Scoped producer corrections

At <=1pm the Model already omits the oxide. Its profile now also removes that
numerical residual from the Pd depth origin AND zc; densities, sheet inventory,
age and evolved Model solution do not change. The strict consumer validation
was not relaxed. `coordinate_migration.json` records the lossless cached-data
translation; maximum binary64 thickness-rounding change is2.71e-20m.
`transport_boundary_proof.json` exercises below/at/above cutoff using a retained
actual state and confirms untouched ordinary seal/reference paths. Those are
coordinate unit controls, not additional ALD predictions.

Current/charge fields are named `assumed_incident_current_A` and
`assumed_run_charge_C`. They are fitted/model charge normalization, not an
independently measured incident-current calibration. The transport evolution
still uses its stated provided-IV-as-incident clock; conditional absolute
reaction totals must not be presented as measured counts.

Cap endpoint companions include all retained/lost-D x slow/fast-sputter cases
for each selected material/thickness recipe, not just the optimistic top-10 rows.
Nanolaminate interface inventory is placed only at actual Pd/oxide contacts:
no interface sheet is hidden at the free beam-facing Pd surface, and the final
oxide/Pd-substrate contact is counted. Erosion counts film host layers, not that
substrate-interface inventory.

For a genuinely zero-yield target, `merit()` and `compare()` retain physical
merit0 and the fixed-PdD curve0. Self-normalized R, stitched reported R and its
effective-Ue fit are null where their target denominator is undefined; explicit
`self_R_defined`/`reported_R_defined` flags distinguish that from a measured zero.
There is no invented floor, missing row or exception-suppression fallback.
`zero_merit_proof.json` compares all2676 common numeric fields of42 positive
rows to the received merit-archive implementation exactly, exercises42 zero
rows, partitions raw/normalized caches and observes neutral target response20.632×.
The normalized low-level zero-reference rejection is unchanged.
Literal effective-Ue fits require a positive finite reported curve throughout
the requested weighted 2.5–6.5keV support. A zero/undefined point makes that fit
null; it is not silently removed to fit a narrower interval.

The transport implantation profile remains the existing shifted Pd BCA template,
including its overlayer first-Pd-cell routing approximation. Nuclear slowing in
the actual design composition is explicit in `Stopping`; this is not a new
composition-specific BCA calculation or a coupled cap/Pd exchange prediction.

The g5_recomb two-face result1681 supersedes the one-face point1675 and qualifies
mesh convergence of LOCAL NO-BEAM equivalents: gas Pd(theta0.781449ML;s0.0238816)
and bulk-fed PdD(theta0.999968ML;s5.06313e-10). Capture0.5 remains ASSUMED and bulk
v(z) is held at its actual7200s run-end profile. These are not universal sticking
inputs or an all-beam/material boundary conversion; marked pass4 values remain.
The earlier0.5ML local-harmonic mean binding1647 is not substituted as a
differential adsorption energy or barrier. No frozen reference is replayed.

## Files in the repository and in the delivery asset

The repository holds every file under 2 MB, plus gzip copies of `all_designs.csv` and `resolved_yields.csv` (`*.csv.gz`, decompress with `gunzip -k`). The other large files (the transport and profile caches, depth profiles and per-row compare tables) are in the collaborator-only delivery asset `g5_ald.tar.gz` on the draft release `g5_ald-delivery`. They are gitignored here, and `SHA256SUMS` lists all of them.
