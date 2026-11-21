# Site distributions, depth populations and binding-controlled inventories

## Decision

The continuous stationary tail fits the digitized UC Davis curves almost as well
as the transported sigmoid: **χ² 52.595 versus 51.479** on the same 48 points.
The discrete two-site law gives **59.455**. These approximate-error shape scores
do not identify a microscopic site, a posterior population, or calibrated counts.
A depth-dependent population can reproduce the central low-energy ratios with
small full-curve cost, but its depth, occupation and high-U value trade against
one another. Carry the **inventory-versus-U curves**, not a best-fit fraction
alone, to the site calculations. `site_classes.md` and its derived crossing table
contain the primary-source inventory comparisons and conditional part-(b) targets.

Authority is g7_micro(a), entries 2570/2573, extended by 2586 and 2592. No DFT,
pre-release Zenodo access, publication, repository commit or push was performed.
The frozen source main is `2ca008af1ad630056e65d2e7c19452c668bd3c87`;
`results.json` records SHA-256 identities of every consumed shared source/input.

## Exact model and meanings

Primary observations are `g4_beam/forward/errbars.csv`: Figure 2 normalized
responses with digitized approximate log errors, including a 5% floor. They are
**not raw calibrated counts**. All models use the same observed E<6.5 keV CM
mask, with the 6.5-keV reference fixed, and all four foils Pd/PdD/Ti/TiD. There
are 48 observations; PdD and TiD contribute 12 and 10. Low means E<2.5 keV CM;
high means 2.5≤E<6.5. No hidden interpolation fills missing observed energies.

`model.py` calls the actual `forward.consumer.csda_nodes`, source
`merit.reference` profiles and shared `g6_plateau.forecast.cross_section`.
It keeps the D–D proton threshold background (40 meV width, 115° interference,
7 fm nuclear parameter in `reaction_S`) and Bosch–Hale neutron baseline; stopping,
profiles and sheet inventories never change when screening changes. The beam is
85% D+ / 15% D2+ **electrical-current** fraction, two deuterons per molecular ion,
5% FWHM Gaussian spread and the existing five-node rule. Depth integration keeps
the shared 0.1-keV lab cutoff and stopping ODE. Layers are SI
`(zt_m,zb_m,n_host_m3,n_fuel_m3,Z,M_amu,n_metal_m3)`; the sheet is separate m⁻².
Cross sections are mb and physical yields are reactions per **true charged C**.

The plotted comparison retains high-voltage focus b_hi=1, b_lo=0; common current
readings 0/1/2 are compared, not selected separately for each family. Proton
signals below 4.5 keV are stitched to neutrons at 4.5 with one fitted common s;
neutrons at ≥4.5 are normalized to their 6.5 value. No per-material amplitude is
introduced. Detector/cup corrections enter only this reported observable, not
raw physical Y/C or D–T enhancements.

The laws are

- `two_site`: p(U)=(1−f)δ(U−0.040)+fδ(U−U_high);
- `exponential`: the same bulk atom plus weight f in
  exp[−(U−0.040)/scale]/scale above 0.040 keV;
- `sigmoid`: two U parameters per material, the original incident-energy law
  with **fixed external E0=1.472768567101136 keV nominal D–D CM**;
- `anchored_rise`: 0.040+positive_amplitude/(1+exp[(incident_lab/2−E0_m)/0.5]),
  two fitted parameters per material (amplitude and E0_m).

All have nine internally fitted parameters: two per foil plus the same stitch.
The external sigmoid knee was previously fitted on these very curves; freezing
it gives a **conditional matched-count comparison**, not independent evidence
for E0. Counting that earlier choice costs one further parameter (original k=10).
The anchored rise is a genuinely nine-parameter positive-slope control.

`shared` assigns the same population law to the already supplied sheet and foil.
`bulk40` leaves the same sheet at 40 eV while varying only occupied foil fuel.
Neither invents sheet sites or a sheet thickness. A fitted f is occupied-fuel
number fraction at each depth, not defect volume, available capacity or reaction
share. The kernel averages **cross sections/yields**, never σ at mean U. The
40-eV bulk is a declared comparison anchor motivated by prior cell calculations;
37–41 eV BO lowering at sampled R=0.1 Å does not independently establish an
energy-independent barrier shift over the nuclear turning region.

## Settled stationary fits

Power 1, shared sheet, full four-foil mask:

| Law | Total χ² | Low contribution | High contribution |
|---|---:|---:|---:|
| Fixed-knee sigmoid | 51.4793 | 15.5218 | 35.9575 |
| 40-eV anchored positive rise | 103.5418 | 25.9696 | 77.5722 |
| Two site | 59.4547 | 17.9995 | 41.4552 |
| Exponential | 52.5947 | 15.0380 | 37.5566 |

| Material | Two-site f | U_high (keV) | Tail weight f | Tail scale (keV) |
|---|---:|---:|---:|---:|
| PdD | 0.0161578 | 2.62576 | 0.312796 | 0.513766 |
| TiD | 0.0268854 | 1.72555 | 0.822539 | 0.294229 |

The tail weight is **not** a rare-site fraction: most of that continuous tail is
near the bulk anchor. Its far tail can dominate yield. Full four-foil parameters,
residuals and reported curves are in `results.json`/`fits.csv`, including controls.
All 28 optimizer runs returned success; this is finite optimization, not proof of
a unique global minimum.

Changing the common sheet convention is decisive: with bulk40, power-1 totals
are sigmoid119.242, anchored172.871, two-site174.974 and exponential167.137. At
power0/shared the corresponding totals are65.053/219.583/136.062/107.125; at
power2 they are146.992/144.273/133.066/141.478. Thus neither a site attribution nor
a stationary-versus-sigmoid conclusion is independent of measurement and sheet
assumptions. The fixed-E0 power0 sigmoid is not its separately optimized old fit.

`profile.csv` gives both U profiles with *all* other fit nuisances reoptimized,
and an explicitly labeled 2D f,U grid with other shapes fixed and stitch profiled.
For descriptive Δχ²≤5.99, the latter spans PdD f=.00750–.02689,
U=2.4977–2.8675 keV; TiD f=.01334–.04217, U=1.6504–1.8949 keV. These are grid
extents, **not confidence intervals or posterior support**, nor rectangular
independent ranges. Full U-profile extents on that threshold have the same
listed bounds. Neither admits the 10⁶-keV upper scan boundary. The exact grid
and thresholds 2.30/5.99/11.83 remain available; no statistical coverage is claimed
under digitization, ignored reference covariance and model misspecification.

At fixed local E, the retained Eq.7 factor turns over at E+U=B_G²≈986 keV for
D–D: increasing U without bound eventually *reduces* absolute yield. A high-U
boundary is not automatically evidence for a rare stronger site. Numerical
single-site ratio envelopes also show a stationary positive mixture cannot
reproduce every central sigmoid-like rise. For example a separately evaluated finite
.04–2048-keV grid gives power-1 PdD R(.25)/R(1)≤1.457 (shared) or .154 (bulk40)
versus central5.626; TiD≤2.605 versus4.028. This is not a certified continuous
bound. Large digitized errors make a good total approximate score compatible
with missing those central ratios; the exponential result must not be called
excluded on that basis.

A diagnostic high-energy-only fit leaves low-energy scores873.50/4887.33/497.34/
82.25 for sigmoid/anchored/two-site/exponential, respectively (all foils).
These are not external holdout likelihoods: E0 and several model choices already
used the same experiment. The low points materially constrain tail shape.

## Depth-resolved population and inventory curve

`depth.py` scans d=1,2,3,5,10,20,30,50 nm and U from .08 to10⁶ keV. It uses
`csda_nodes(...,break_m=d)` to split the integration exactly, f_s above d and
0≤f_b≤f_s below. The shared sheet uses f_s; bulk40 does not. It never infers
unmeasured phase populations from density or names in the transport model.
For fixed d,U, the central constraint R(.25)/R(1)=observed is **linear** in f_s,f_b;
`depth_region.csv` retains the complete feasible segment’s projected bounds and
the lowest whole-curve cost along it. Central-ratio feasibility is not a fit band.

| Material/sheet | d (nm) | f_s | f_b | U_high keV | Full one-foil χ² | s |
|---|---:|---:|---:|---:|---:|---:|
| PdD shared | 10 | .027989 | .001708 | 2.86564 | 2.9154 | .9584 |
| PdD bulk40 | 2 | .34434 | .002183 | 2.86564 | 3.1255 | .9660 |
| TiD shared | 20 | .018337 | .005790 | 1.83213 | 6.4805 | 1.4388 |
| TiD bulk40 | 20 | .018376 | .005815 | 1.83213 | 6.4718 | 1.4433 |

These are power-1 grid examples, with **one independently profiled stitch per
foil**, not an equal-parameter four-foil contest or a demonstrated common-stitch
solution. Extra depth freedom improves central PdD curvature, while the TiD
stitch differs appreciably. Both the entire curve cost and that relaxation must
remain visible before selecting a physical class.

`depth_inventory_curves.csv` profiles d,f_b,sheet and stitch on the power-1 grid.
It separates ratio-only feasible occupied inventory from the full-data
**descriptive Δχ²=4** band. Its ordinate is occupied special D per encountered
metal atom at exact nominal DD-CM profile key1, not available special sites:
**N_s ≥ N_special,occupied** for single occupancy. Blank full-band fields mean
no accepted full-curve solution at that U. Bands bracket objective crossings on
an81-node segment grid plus the optimum; they are finite-resolution descriptions,
not a calibrated likelihood or a proof of all disconnected branches.

The separate `full_band_volume_occupied_min_per_metal`/`...max...` columns do
**not** credit high-U fuel in the sheet to a bulk defect. They change class
triage: coarse annealed-Pd GB capacity ≤1.50×10⁻⁵/Pd is below the smallest full
band volume demand 2.48×10⁻⁵/Pd. The narrow Ti one-site-per-line-repeat scenario
also has no sampled intersection. A 0.5-nm Ti line cylinder survives only
2.866–3.205 keV, while the quoted Ti GB geometry survives 2.291–3.205 keV.
These are necessary-capacity tests conditional on the depth/ratio/error model,
not universal rejection of the classes; full source/conversion qualifications
and the retained vacancy/subsurface/second-phase cases are in `site_classes.md`.

`inventories.csv` uses clipped actual stopped paths, no nuclear weighting.
At key1 the encountered volume D/metal is .0103203 (PdD) and1.84971 (TiD);
the separate sheet contributes .00357409 and .00478510 per metal. The maximum
occupied D in the first1/3/10/20/50nm is respectively PdD
.00021462/.00064873/.00221741/.00458756/.01030515 and TiD
.00755833/.07032551/.28404891/.58389097/1.39764024 per encountered metal.
A constant fitted f cannot be equated to a fixed geometrical subsurface width.

## Binding-controlled populations: independent summary constraints

`occupancy.py`/`occupancy_results.json` use the same two-class Fermi law form,
not the same numerical parameters across reactions. A special capacity N_s
replaces ordinary sites in an explicit total capacity C; positive Δbinding raises
its equilibrium occupation. The chemical potential is fixed by loading x:
θ_b=u/(1+u), θ_s=Ku/(1+Ku), x=(C−N_s)θ_b+N_sθ_s,
f=N_sθ_s/x, K=exp(Δbinding/kBT). These assumptions do not establish beam-time
equilibrium or independently measure site capacity.

[Raiola2005 / EXFOR O1379](https://raw.githubusercontent.com/IAEA-NDS/exfor_master/master/exforall/O13/O1379.x4)
includes x=2.1/1.1/.26/.23/.20 at−10/50/100/150/200°C and
U<30/<50/250±40/295±40/290±65eV. The two limits are hard inequalities here;
the source supplies no confidence level. Capacity2.2/Ti, ordinary U=0 and a
5-keV-CM penetration proxy give descriptive χ²≈6.01 on the three detections,
predicting26.5/50/195.3/218.1/246.8eV. N_s trades from~.155/Ti at U_high=.308keV
to~1.23×10⁻⁵/Ti at30.47keV with almost unchanged objective. In the original
U≤100keV slice, the optimized-grid Δχ²≤1 ridge projects to
N_s1.36×10⁻⁶–.223/Ti, Δbinding.21–.5eV and U_high.224–100keV. The delivered curve
also exposes a mathematical extension to10⁴keV, reaching N_s≈8.67×10⁻⁸/Ti with
essentially the same minimum. That huge-U penetration-proxy extrapolation is
**not a physical asymptote or a part-(b) target**; finite bounds are scan boundaries,
not estimates. Ordinary40eV already violates the<30eV limit. The ordinary anchor
must not be imported blindly from
UC Davis. Pivot3.85–10keV changes the minimumχ²7.70–3.86; no original Ti curve
weights or full thick-target refit are claimed. Temperature, loading and phase
change together, so this is not an isolated loading experiment.

[Cvetinović2023](https://doi.org/10.1088/1742-6596/2586/1/012101) and its
[primary PLB paper](https://doi.org/10.1016/j.physletb.2023.137684) concern
**F+D**, about.296–.882MeV CM, with fitted U18.2±3.3keV hard rolled versus
3.2±1.9keV commercial soft. The ratio5.7 is of fitted shifts, not yields. Initial
total H(D)/Pd was .47/.70 using mixed85%D/15%H gas; both lost D under irradiation.
The soft foil is not a controlled annealed mate. Loading, cycling, chemistry and
preparation confound a causal cold-work interpretation. Distinct F+D penetration
(B_G≈379.984√keV, adiabatic reference2.19keV), not the D–D kernel, is used here.
At fixed capacities, binding and temperature, loading alone bounds the special
response-excess ratio by .70/.47=1.489, versus central16.5–20.2 across the stated
energy pivots. At500keV the soft central shift allows hard≤3.689keV, with relaxed
minimumχ²≥11.33. Sample-dependent capacity can interpolate but cannot identify
Δbinding/U_high: at Δbinding.2eV and U_high30keV, hard/soft N_s≈.2495/.02110 per
Pd (occupied fractions.5303/.0301). Their uncertainty-band curves and inventory
**change** are in the occupancy artifact, not a false universal site value.

[Huang2026 Fig.4](https://doi.org/10.1038/s41586-026-11238-4) supplies a physical
motivation for chemically selective populations: D-rich Al-rich L1₂ versus
D-poor adjacent Nb-rich D0₂₂ in Inconel718. A deliberate10s hold before quenching
accentuated loss from the lower-affinity phase; it is not an equilibrium partition
coefficient. The approximately−.28eV is an absolute solution energy in that alloy,
not fitted Pd/Ti relative binding. The no-new-defect TEM/positron controls were on
a separate (CoCrNi)94Al3Ti3 specimen. Minor vacancy/elastic/APT contributions were
not fully excluded. Neither local1–1.6at.%D nor this alloy’s site chemistry is a
measured second-phase capacity or a nuclear-screening value in Pd/Ti.

## Conditional D–T carry

`dt.csv` uses the existing plateau D–T Bosch–Hale kernel and fuel-recoil stopping;
**no D–D threshold resonance is transferred**. Freeze each PdD/TiD number profile
and its sheet, replace target D by T, and keep the fitted physical U law and
85/15 beam. This assumes unchanged chemistry, retention, electronic law and
relevant fuel site occupation under isotope replacement. It is not a fitted D–T
measurement. Lab1/2/5/13keV use exact source nominal DD-CM keys.5/1/2.5/6.5,
which must never be relabeled D–T CM. D–T CM instead uses m_T/(m_D+m_T) times lab.
Low-energy Bosch–Hale continuation and retained stopping remain model assumptions.

Shared-sheet examples below list raw mixed reactions/true C, screened/bare ratio,
and **yield-weighted** U=∑p U Y/∑p Y in keV. The exponential moment uses the same
tail quadrature, not the population mean or a yield-equivalent shift.

| Host law | Lab keV | Y/C | Screened/bare | U_yield keV |
|---|---:|---:|---:|---:|
| PdT two site | 1 | 1.17595e4 | 1.21207e9 | 2.62576 |
| | 2 | 4.37493e4 | 1.90710e4 | 2.62554 |
| | 5 | 6.70729e5 | 5.57811 | 2.10572 |
| | 13 | 1.60248e8 | 1.12475 | .285513 |
| PdT exponential | 1 | 1.39176e4 | 1.43449e9 | 4.06604 |
| | 2 | 3.91021e4 | 1.70452e4 | 3.66292 |
| | 5 | 6.63361e5 | 5.51684 | 1.80244 |
| | 13 | 1.77219e8 | 1.24387 | .383761 |
| TiT two site | 1 | 3.10188e4 | 5.27080e8 | 1.72555 |
| | 2 | 2.29538e5 | 6.77189e3 | 1.72513 |
| | 5 | 1.54569e7 | 3.78055 | 1.22359 |
| | 13 | 9.99264e9 | 1.12148 | .207227 |
| TiT exponential | 1 | 4.21003e4 | 7.15379e8 | 2.76968 |
| | 2 | 1.78161e5 | 5.25617e3 | 2.30363 |
| | 5 | 1.88290e7 | 4.60532 | .868998 |
| | 13 | 1.18982e10 | 1.33534 | .375636 |

The CSV also includes the separately refitted bulk40-sheet cases, exact bare
Y/C denominators and the two representative shared-sheet depth laws:

| Depth example | Lab keV | Y/C | Screened/bare | U_yield keV |
|---|---:|---:|---:|---:|
| PdT, d=10 nm | 1 | 3.29132e4 | 3.39239e9 | 2.86564 |
| | 2 | 8.45699e4 | 3.68653e4 | 2.86552 |
| | 5 | 8.91223e5 | 7.41185 | 2.44197 |
| | 13 | 1.59471e8 | 1.11930 | .301867 |
| TiT, d=20 nm | 1 | 3.03477e4 | 5.15676e8 | 1.83213 |
| | 2 | 1.87397e5 | 5.52864e3 | 1.83158 |
| | 5 | 1.14157e7 | 2.79213 | 1.10218 |
| | 13 | 9.57954e9 | 1.07512 | .123830 |

Their exact f_s/f_b/U parameters are the one-foil examples above, **not a
common-stitch selected model**. The actual depth moment passes the uniform-law
identity against the stationary mixed moment (tolerance2×10⁻⁶). Giant low-energy
ratios have tiny bare denominators; compare absolute Y/C too. Profile and law
variation is not a D–T uncertainty interval.

## Per-run consumer and reproduction

From the checkout root with NumPy2.5.4/SciPy1.18.1:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m g7_micro.analyze
# Narrow consumer for the newly authorized depth/inventory curve:
python -m g7_micro.depth
# Focused postprocessing of accepted depth cells, without refitting:
python -m g7_micro.depth --volume-bands-only
# Rebuild just D–T from the saved stationary/depth fits:
python -m g7_micro.analyze --dt-only
# Independent binding-law summary analysis:
python -m g7_micro.occupancy
python -m g7_micro.site_targets
# After release and an actual consistently masked/calibrated run join:
python -m g7_micro.score_runs runs.csv log_covariance.npz profiles.json laws.json
```

`score_runs.predict_run` returns raw reactions/true C on an actual supplied
profile. All original four-hypothesis CSV fields and covariance semantics in
`g6_plateau/README.md` remain mandatory. Add `material,lab_keV,isotope,channel,
D1_fraction,D2_fraction,spread_FWHM,profile_id,yield_convention`; convention must
be `reactions_per_true_C`. DD channel is p/n/total; DT is fusion. Profile JSON maps
profile_id to `nominal_DD_cm_key_keV,layers,nads_m2,state,source_id`; the exact key
must equal lab/2, even for D–T. Do not substitute nominal profiles for missing
actual source/history metadata. Law JSON maps model name to material-specific
`{family,a,b,sheet}` (a=f,b=U_high or scale in keV); a depth law instead carries
`{family:"depth",f_s,f_b,d_nm,U_high_keV,sheet}`. The laws conserve existing fuel.

The adapter calls the **unchanged real `g6_plateau.score_runs.score`** with the
identical joined rows, covariance and shared normalization design. Each new law
replaces the screeningrise prediction slot for its paired comparison, without
altering the frozen four-cause registry or declaring a fifth cause. NPZ IDs are
joined by identity, not row order. No independent masking, guessed charge,
invented exposure/counts or pre-release observations are used. Nonpositive or
sparse signals require the owning raw-count likelihood, not this positive-log
scorer. `g4_zenodo` retains ingest/mask/calibration ownership after2026-12-01.

The settled analyze CLI took156.94s, maxRSS222600KB with a1.5GiB virtual-memory
cap. Point40eV/3keV identities against the real shared DD forward and DT plateau
consumer agree within3.34e−16; f=0 identity is exact. Exact-versus-table χ² differs
by≤.002752; fitted exponential yields at96 versus192 tail nodes differ≤2.36e−14.
The actual adapter CLI was exercised on an explicitly labeled deterministic
arithmetic fixture (not experimental observations or raw-data simulation), with
reversed covariance IDs and a shared nuisance normalization. Its point3keV paired
χ²=911.695724535 exactly matches the unchanged scorer; f=0 χ²≈1.8e−30 and uniform
depth-bulk χ²=0. Depth uniform-population identities agree≤4.45e−16. The finite
exponential rule omits probability<1.000003e−10 outside its declared x window;
these numerical checks are not physical-error bounds.
The extended depth run took52.22s/159056KB; its accepted-cell volume-only
postprocess took45.17s/153836KB, and the final40-row D–T-only consumer took
6.48s/95912KB. No stationary or broad depth fit was replayed for those updates.

Actual Pyright stdio-LSP definition/references queries reached `predict_run` and
its cross-file consumers; affected diagnostics were clean after correcting the
workspace import scope. A positive Python AST query found the shared CSDA call.
The harness `xd://lsp` device was unavailable, so no nonexistent editor-device
success is claimed. No unrelated suite, build, formatter or old snapshot changed.

Artifacts: `model.py`, `analyze.py`, `depth.py`, `score_runs.py`; `results.json`,
`fits.csv`, `profile.csv`, `regions.csv`, `dt.csv`, `depth_results.json`,
`depth_region.csv`, `depth_inventory_curves.csv`, `inventories.csv`; the independent
`occupancy.py`, `occupancy_results.json`, `site_targets.py`, `site_targets.csv`
and `site_classes.md`.
The three dominant unresolved inputs are actual per-run measurement covariance
and charge, actual spatial/phase/retention populations including the sheet, and a
class-specific barrier law on the relevant approach geometry. The next useful
part-(b) comparison is a **site count plus its reaction-conditioned U target on
the curve**, not a required equality between one0.1-Å BO sample and a fitted U.

## Lead's note

This directory is part (a) of `g7_micro`, delivered as `g7-micro-a-2ca008a.tar.gz` (SHA-256 2ac20698cc42b4bb5c071bf16391589ccbe772115790603360c7e1445e75902b, 19 files) on the draft release for this node. The files are unchanged except for this note. The lead replayed every step in an isolated checkout at 2ca008a under NumPy 1.26.4 and SciPy 1.14.1, with ASE, GPAW and matplotlib blocked:
- all 24 stationary-fit χ² totals agreed within 0.0011;
- fitted fractions and screening values agreed within 0.05 percent, and so did every D-T row;
- `site_targets.csv`, `inventories.csv` and `occupancy_results.json` were identical.
The one large parameter difference is a diagnostic high-energy-only fit whose site fraction is about 2e-9, so its U value is not determined; its χ² is unchanged.

What it means, in plain terms:
- **A spread of site screening fits the measured curves.** A continuous tail of site values above the ordinary 40 eV fits the UC Davis curves about as well as the old rising-screening law (χ² 52.6 against 51.5 on 48 points). A single strong-site value fits somewhat worse (59.5): it puts 1.6 percent of the deuterium in PdD at about 2.6 keV and 2.7 percent in TiD at about 1.7 keV. Placing the strong sites in the top 10-20 nm reproduces the low-energy rise better than spreading them through the depth. None of these fits identifies a site type.
- **The data fix a curve, not a single site.** Fewer sites with stronger screening fit as well as more sites with weaker screening. In the titanium loading series the fit holds from 0.15 sites per Ti at 0.3 keV to about 1e-5 sites per Ti at 30 keV.
- **Some site types have enough room and some do not.**
  - Fine grain boundaries in polished foils, wide dislocation cores and grain-boundary bands in Ti have enough room for the share the UC Davis curves need.
  - Coarse grain boundaries in annealed Pd do not.
  - No usable count exists yet for vacancy clusters or second phases.
- **Cold-rolled Pd points to more strong sites.** If each site type binds and screens the same way in both foils, lower loading alone raises the strong-site response by at most 1.5 times. The hard and soft foils differ by 16-20 times in that response (F+D data), so the hard foil would hold more strong sites, as cold work makes. The two foils also differ in source, preparation and deuterium loss, so cold work is not isolated.

These keV values are effective shifts in the beam model, not the static lowering at 0.1 A. Comparing a computed site with them needs the barrier-penetration (action) calculation along the computed curve, not the 0.1 A value alone. Part (b), the DFT of the named site types, follows `g7_hydride/`.

First stage of part (b): `zr_calibration/`. An interstitial oxygen atom beside the target deuteron in ZrD2 changes the static lowering by less than 1 eV. Oxygen as a neighbour atom is not a strong site.

Second stage of part (b): `pd_oxide/`. A PdO sheet on a Pd slab lowers the static screening at the deuteron beneath it by 1.5 eV at 0.1 A and 5 eV at contact. The deuteron there sits one metal layer under the oxide (nearest O 4.2 A away); a deuteron placed in the interface itself (`pd_oxide/interface/`, nearest O 1.5 A) gives 38.7 eV at 0.1 A and 39.4 eV at contact, the same as the buried sites.
