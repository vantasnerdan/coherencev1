# Explaining the screening gap: candidate ledger and deciding measurements

Private scientific delivery. No publication, push, new DFT, embargoed per-run access, BO-converter start, or Fleet selection is part of this result.

## Decision

**Keep the material-site/profile explanation and a separately testable beam component at the front. Ordinary equilibrium target motion and the computed one-ion electronic response are much too small to supply the reported hundreds of eV. A genuinely pre-existing energetic population or a long-lived beam-prepared material state is a different, still parameter-dependent hypothesis—not energy borrowed from one cold recoil.** A D–D threshold term is a real channel-specific sensitivity, not an explanation that transfers to Li or F reactions.

The prominent full table is **`candidate_ledger.csv`**: 15 candidates, each with a numerical prediction or parameterized direction for all seven requested axes—charge, beam energy, temperature/loading/mobility, preparation, current/dose, metal/insulator, target D→T—and an assessment tied to source-row keys. **`experiment_rows_addendum.csv`** has 38 separately keyed published/derived records, explicit observable/normalization/error fields, and source qualifications. It corrects or extends the existing experiment ledger without editing its maintained generator. Blank error fields are unknown, not zero. `source_observables.md` is the consumed source specialist's working packet, with its CM conversion corrected here; it is not needed in the delivery payload.

| Family | What the numbers permit now | Decisive unresolved input |
|---|---|---|
| Strong sites / near-surface environment | Existing same-curve fits support distributions: tail χ²52.595 vs fixed-knee sigmoid51.479; two-site59.455 on48 digitized points. Conditional occupied PdD fraction .01616 at2.626keV, TiD .02689 at1.726keV. | Actual barrier law and occupied capacity; part(b) unknown. This is not a static keV-site measurement. |
| Equilibrium zero-point/thermal motion | At0.5keV lab: local bare-model rate log enhancement .0492 at300K; apparent shift .398eV. Tied harmonic300K gives .505eV. Actual frozen-PdD own-reference consumer gives only **+.276 / +.350eV** around its declared40eV comparison law. | Momentum-distribution/bound-state approximation, not a hidden keV kinetic scale. |
| Diffusion / stopping / evolving fuel profile / layers | Can change the energy-weighted yield without changing the elementary potential; temperature and impurity data make these compulsory shared parameters. | In-beam profile, species, stopping and charge at matched history. Not enough to declare a unique explanation. |
| Prompt cold recoil | Same binary pair gains no relative energy; general two-moving-D single-cascade bound is Erel≤Elab, not≤Elab/2. At nominal.25keV CM even maximal local two-D energy gives σ/σ(3keV-shift)=1.77×10⁻¹¹ under the stated model. | Encounter/path multiplicity and correlated distributions; local bound is not a complete thick-yield bound. |
| Pre-existing energetic D | At.5keV lab an independent all-hot Maxwell target with component θ=100eV gives apparent531eV; fhot=.00132 could mimic300eV locally, but even f=1 cannot mimic1 or3keV there. | Population, energy supply, spatial distribution and lifetime. None is inferred or invented. |
| Heated track / pulse | Conditional source narrow-resonance model gives ~10⁻¹⁷ reactions/incident D at1080K and its nominal exposure. Common source/UC-Davis track exposure needs a factor284 change under one-current-power normalization. | Actual T(r,t), density moment, hot screening, resonance/branch and track exposure; no absolute shared mechanism established. |
| Above-nominal neutrals / beam normalization | Real high-energy neutral route exists before deceleration; source-specific constraints matter. Figure/SI normalization comparison differs by4.68–4.89 at.25 vs3keV CM, but SI alone retains substantial low-energy excess. | Raw counts/live time/true charge, energy-resolved species and proton kinematics. No assertion of a duplicated current correction. |
| D–D threshold term | Frozen-background total-DD multiplier5.02 at.5keV lab; not carried to DT, Li or F. Published IPC/p is.42±.05 at10keV lab vs.060±.012 at20keV. | Entrance model, interference and channel response; measured branching is not fitted width ratio14. |
| One-ion / many-ion electrons | One2keV head-on encounter difference.22eV at.1Å, not time-step checked. At quoted flux,10nm radius/100ps overlap is3.46×10⁻⁷. | A longer-lived, populated state would need its own lifetime and dielectric/capture response. |

These are different observables, **not additive pieces of a 3keV budget**. A reported fitted U requires its forward model, channel, window, normalization and likelihood. The 40.55eV collision result is a potential lowering at R=.1Å, not a constant fitted U.

## 1. Governing observable and target-motion calculation

For projectile mass m1, target mass m2, α=m2/(m1+m2), β=m1/(m1+m2), laboratory projectile kinetic energy L, target kinetic energy K, and c=cos(angle between velocities), exact nonrelativistic relative kinetic energy is

`Erel = αL + βK − 2√(αβLK)c`, with stationary `E0=αL`.

The per-path reaction probability is `nD <σ(Erel) g>/vbeam`, where `g/vbeam=√(Erel/E0)`. Averaging σ alone describes a different observable. For an isotropic Maxwell momentum distribution K/θ has gamma shape3/2; the script integrates its radial distribution and c with48×48 Gaussian nodes, in the log domain. No Gaussian target-energy shift is substituted for this kinematics.

The transparent local comparison freezes astrophysical S:

`σ(E,U) = S/[√(E(E+U))] exp[−b/√(E+U)]`,
`b=31.397 (Z1Z2)√(μ/μDD) √keV`, `μDD=1.00705u`.

The apparent local U is the smallest nonnegative solution of

`ln F = ½ ln[E0/(E0+U)] + b(1/√E0−1/√(E0+U))`.

This is the retained Eq7 convention, not σbare(E+U)/σbare(E), a likelihood fit, or a nuclear-action calculation. `results.json` distinguishes `logF_sigma` and flux-weighted `logF_rate`. Masses use the repository's atomic-isotope convention. Relativistic and detailed charge-state corrections are not inferred.

**Gamow convexity and charge.** Put x=b/√E0. For small target velocity variance, the rate-weighted expansion is

`ln F_rate = βθ [b²/(4E0²) − b/(2E0^(3/2))] + O(θ²)`
`          = s²(x²/2−x)+O(s⁴)`, `s²=μθ/(2m2 E0)`.

Thus, in the low-energy convex regime, raw log enhancement grows approximately as `(Z1Z2)²`, but division by the screening sensitivity `(x−1)/(2E0)` yields

`Uapp ≈ (βθ/2) x(x−2)/(x−1)`.

Apparent U can therefore scale nearly linearly with charge too. **Charge scaling alone does not separate target motion from a potential shift.** At fixed relative velocity E∝μ, x∝Z1Z2/v; at fixed CM energy x∝Z1Z2√μ. Changing projectile/target orientation changes β and the target's momentum distribution. Neither fixed voltage nor a mixture of these experiments is a fixed-velocity test. Convexity is not universal at arbitrary E or across a narrow nuclear resonance.

### Local motion numbers

Entries are natural-log rate enhancement and apparent Eq7 U in eV, against U=0 with frozen S; not increments to the static40eV sample.

| State | .5keV lab / .25CM | 2keV lab /1CM |13keV lab /6.5CM |
|---|---:|---:|---:|
| Classical300K | .04919 / .3985 | .002981 / .1962 | .00006315 / .07256 |
| Classical473.15K (200°C) | .07744 / .6278 | .004701 / .3094 | .00009960 / .1144 |
| Harmonic300K, hbarωD=.047607eV | .06231 / .5050 | .003780 / .2487 | .00008006 / .09199 |
| All targets independently hot, θ=1eV |1.711 /14.45 |.1137 /7.525 |.002441 /2.805 |
| θ=10eV |9.947 /105.3 |1.014 /70.22 |.02424 /27.94 |
| θ=100eV |26.698 /531.0 |5.542 /498.0 |.22699 /269.1 |

The oscillator value follows `hbarω=(hbar c)√(Kii/mc²)` from saved PdD Kii=1.0917228eV/Å² and mass1.875612945GeV/c². Its source status is **provisional local well: absolute Pd spectrum fails**; it is an Einstein scenario, not measured dispersion. A harmonic Gaussian uses `θeff=(hbarω/2)coth[hbarω/(2kBT)]`, here .0327756eV, **not** kBT+hbarω/2. The30–100meV sensitivity scenarios remain in JSON. Convolving bound-state zero-point momentum with an on-shell cross section is an impulse approximation; it is not a derivation of a lattice bound-pair fusion rate. Strong binding correlations could require a different calculation.

For thermally diffusing D, the activation barrier is potential energy along a hop, not an extra projectile kinetic energy. A mobile fraction `exp(−EA/kBT)` reduces an explicitly mobile contribution; an unbiased selected thermal subset cannot exceed the all-target thermal result at the same T. Mobility can strongly change supply, profiles and sites even though its direct kinematic correction is small.

### Thin versus thick and the actual source normalization

The illustrative homogeneous thick calculation is

`Y(L,θ) ∝ ∫[.1keV,L] <σg>/vbeam /√l dl`.

This freezes fuel density, assumes pure D+, sqrt-energy stopping and a finite.1keV laboratory cutoff. It is neither a measured stopping law at zero nor the full experimental target. At.5keV lab its classical300K log enhancement is.05250 and yield-equivalent U=.4052eV; harmonic300K gives.06651/.5135eV. At13keV lab their U equivalents are.08079/.1024eV. Normalizing each curve to13keV changes log enhancement to.05241/.06639 at.5. These are not thin-kernel values relabeled as thick yields. Extending Eq7 with positive constant U and sqrt stopping to zero would produce a formal logarithmic low-energy singularity; no E→0 plateau theorem is claimed.

**Experiment-specific calculation:** `actual_source_profile_motion` instead consumes actual saved PdD energy-specific profiles and sheet, the shared CSDA stopping nodes,85/15 D+/D2+ **electrical-current** fractions (molecular component supplies two D),5%FWHM source spread assumption, .1keV lab cutoff, and the common proton/Bosch–Hale neutron kernel. Moving targets average that actual energy-dependent S(Erel), not frozen S. Below4.5keV nominal CM the proton increment is referenced to the6.5keV neutron increment; a fixed stitch cancels in the changed/base shape enhancement. Stopping and target profiles never change in this calculation. A pointwise constant-U curve inversion around the **declared40eV comparison law** supplies the following diagnostic, not a new fit to counts:

| Lab keV |300K log shape enhancement |300K ΔUapp eV |Harmonic300K log shape enhancement |Harmonic ΔUapp eV |
|---:|---:|---:|---:|---:|
|.5 |.026690 |.2762 |.033821 |.3500 |
|1 |.008369 |.2231 |.010609 |.2828 |
|2 |.002468 |.1833 |.003128 |.2324 |
|5 |.0003931 |.1372 |.0004984 |.1739 |
|13 |0 |not identifiable at reference |0 |not identifiable at reference |

At.5keV lab stationary/moving harmonic proton yields are5.5311×10⁻¹³ /5.7220×10⁻¹³ reactions/true C. At13keV stationary/moving neutron yields are7.38682×10⁵ /7.38757×10⁵. Tiny bare-like low-energy denominators and large ratios do not imply a measured absolute signal. This model is deliberately not the fitted1.4–1.7keV background, and these increments must not be added to that fit without recomputing it.

**Independent hot population:** for `Fmix=1+fhot(Fhot−1)`, the.5keV lab calculation with θhot=100eV requires fhot=.0013155 to mimic300eV locally. It requires1352 or3.60×10⁷ to mimic1 or3keV—impossible fractions within that particular distribution. θhot≤10eV cannot reach even300eV with f≤1. More energetic/non-Maxwell distributions are not bounded by those examples; they require an actual phase-space population, spatial overlap, persistence and energy budget. JSON supplies fixed target-K isotropic cases too, without claiming their occupation.

### Recoil energy conservation

For initially stationary equal-mass elastic D+D, the post-collision velocities are orthogonal and `|v'−u|²=|v|²`; the original pair's relative kinetic energy remains L/2. A recoil of energy K≤L hitting a **fresh stationary D** has Erel=K/2≤L/2. It is invalid to keep the projectile at full incident energy and simultaneously grant its recoil another full incident energy.

A general pair of two moving D produced in a cold single-ion cascade is less restrictive: relative energy is at most their total kinetic energy, hence **Erel≤L=2Enom**, with the lattice able to carry momentum. This is an energy bound, not a claim that an optimal counter-propagating pair is populated. At Enom=.25keV, the maximal local σ atErel=.5keV with a40eV comparison shift is3.016×10⁶ times σ(.25,.04), but only1.767×10⁻¹¹ times σ(.25,3keV). Matching that latter reference would need an additional5.66×10¹⁰ encounter prefactor before any path/flux accounting. At Enom=.5keV it needs9.32×10⁵. These frozen-S comparisons do **not** bound the full integrated cascade yield: independent bounds on multiplicity, fuel paths and reaction-weighted phase space are still absent. A preheated material or energy pooled from multiple ions is a separate model.

## 2. What the published charge, temperature and preparation data actually say

The addendum preserves original fit semantics. In particular:

- **Kasagi2002 DD:** Pd310±20(stat)±50(sys), PdO600±20±75, Fe200±15±45, Au70±10±40, Ti65±10±40eV. Detector was a50+200μm ΔE–E telescope at125°, not one90° Si. PdO's≈50 enhancement at2.5keV lab is a relative-yield/bare-relative-yield ratio normalized at10keV lab. [Primary DOI](https://doi.org/10.1143/JPSJ.71.2881).
- **Kasagi2004 Li+d:** PdLix1500±310 vs AuLix60±150eV, D lab30–75keV, thick alpha excitation functions and75keV normalization. Uniform Li profile, stopping and shifted bare cross section are assumptions. This is not the separate “anomalous alpha?” feature on a solid-Li slide, and no narrow peak energy/strength is extracted from that feature. Journal body unavailable; [author abstract](https://strathprints.strath.ac.uk/3020/) and author's slides provide the qualified anchor.
- **Cvetinović2023 hard/soft Pd:** Li+p2860±190eV is a relative hard/soft integrated-yield inference assuming the soft reference is weakly screened. F+p18700±1500eV comes from **1.53±.05 resonance-strength ratio**, with soft strength treated bare. F+d is18200±3300 hard and3200±1900 soft against an independent normal-kinematic CaF2 reference. These are different observables and references. Hard/soft totalH(D)/Pd≈.47/.70, not matched microstructure/occupation. The paper's empirical nearer-quadratic charge inference is not a universal law. [Full primary](https://scoap3-prod-backend.s3.cern.ch/media/files/75014/10.1016/j.physletb.2023.137684_a.pdf).
- **Kinematics correction:** actual atomic-isotope masses give F+p323.31keV CM→6417.99keV lab. The paper's “6417.9MeV” is a **unit typo**, not evidence for a6.466MeV resonance; the latter used integer mass numbers and is wrong on the resonance-width scale. The6410–6819keV lab window is322.91–343.51keV CM in this convention. F+d3089–9200lab is296.09–881.84CM. Nuclear-mass/charge-state refinement is distinct and much smaller; none is silently claimed.
- **Same reaction, different hosts:** original Cruz2005 Li7(p,α) at30–100keV proton lab reports Li2WO4 **185±150**, Li metal **1280±60**, PdLi **3790±330eV**. These original values are not the different expanded PoS recap values. A host-independent nuclear/charge-only model cannot describe this comparison; several material/transport mechanisms still can. [DOI](https://doi.org/10.1016/j.physletb.2005.08.036), author-institution abstract read, not full article.
- **Raiola2005 Ti:** at−10/50/100/150/200°C, D/Ti2.1/1.1/.26/.23/.20 accompanies U<30/<50/250±40/295±40/290±65eV. **Pt goes the other direction:**675±70eV at20°C→465±38 at300°C in the EXFOR author summary (separate Zeng report quotes±50 for the20°C value). C at200°C, loading.15 hasU<50. These distinguish a state/loading transition from universal monotonic thermal kinetic enhancement. [EXFOR O1379](https://raw.githubusercontent.com/IAEA-NDS/exfor_master/master/exforall/O13/O1379.x4) read directly. Limits have no supplied confidence level. Full IOP body was inaccessible to grounding; no claim to have read it.
- **Kowalska2025 O/C:** inferred U changes≈300→600eV under O implantation and≈600→400eV under C, but these are different starting series. The observable is the current-normalized **8/20keV lab thick-yield ratio**, not an8keV cross-section enhancement. Spot-to-spot D density systematics exceed statistics and were not numerically specified. [Primary](https://doi.org/10.3390/ma18061331).
- **Dubey2025:** publishedv2 TableI IPC/p=.42±.05 vs.060±.012 at10/20keV lab, a ratio7.0±1.6. ΓIPC/Γp=14 (at least10) is a fitted partial-width ratio, not those measured yield ratios. Photon/annihilation observables have their own errors and responses; arXivv1 photon errors differ. This is not a sub-keV Davis measurement or a bound-pair width. [Publishedv2](https://arxiv.org/pdf/2408.07567v2).

For transparent charge summaries, the reported U/(Z1Z2) values DD/Pd310, Li+d/PdLi500, Li+p/hardPd953, F+p/hardPd2078eV are not equal. U/(Z1Z2)² gives310,167,318,231eV. Neither sequence identifies a universal law because the target, orientation, energy window, resonance and reference differ. A fixed site's **linear-response** electrostatic interaction is bilinear in nuclear charge, while nonlinear electronic capture/rearrangement can scale differently.

The script's matched-Ecm5keV thermal300K control yields Uapp .0838eV DD, .1639 D+Li6, .4337 Li7+H,1.499 F19+H,1.973 F19+D. At fixed relative velocity (DD at5keV) these are .0838/.1333/.4645/1.538/1.465eV. At actual experimental kinematics the smooth constant-S proxies are smaller still: D+Li6 .046–.075eV; Li7+H .048–.122; F19+H .169–.175; F19+D .137–.246. This distinguishes correct projectile/target orientation from a false “p beam on F” motion calculation. F+p resonance strength requires its actual resonance convolution; the proxy is not a refit of that measurement.

The existing g7_micro binding-population analysis already explains an important sign: saturated rare-site inventory can make its occupied fraction rise as bulk loading falls, while absolute yield can still fall. Its Ti summary fit is degenerate in capacity/binding/Uhigh and uses a5keV penetration proxy, not the original thick-target likelihood. Ordinary40eV would violate the<30eV Ti upper limit if transplanted blindly. In F+d, loading alone at fixed capacities/binding bounds the special response-excess ratio by.70/.47=1.489, versus central16.5–20.2 under that analysis. It needs a different population/capacity or another effect, not an unqualified “cold work proved screening” statement.

Chen's +15±2% neutron-rate response to electrochemical loading at15keV DD CM supports real fuel/loading modulation; it does not independently measure ΔU. Thus it does not conflict with Ti's lower-loading/higher-fitted-U trend. [Chen2025](https://doi.org/10.1038/s41586-025-09042-7); direct intervention evidence supplied by Main, with source scope retained in the addendum.

## 3. Heated tracks, current overlap and electrons

The reused Czerski track model is a **hot-hot/mobile-trapped** reaction source, not the moving-target correction above:

`p_track = nD² exp(−EA/kBT) <σv>res(T,Uhot,ER) A`,
`A=V τ (Γp/Γ)`, not τ alone.

At T1080K, ER1eV, radius7fm, nominalA=1870nm³×100ps=1.87×10⁻³⁴m³s and branch1, the existing SI implementation gives p=2.377×10⁻¹⁷ Zr,2.481×10⁻¹⁷ Ti,8.649×10⁻¹⁷ Pd reactions/incident D using source Table1 densities and Uhot340/370/400eV. These are conditional computations, not newly measured yields. The source-detected plateaus are2.6×10⁻⁴,1.4×10⁻⁴,6.9×10⁻³ counts/(μA ks), respectively—.26/.14/6.9 detected counts/C, requiring acceptance before comparison to emitted reactions. The cold Ti fitted25±10eV is not the assumed370eV hot-track shift.

Sensitivity is large: for Zr at1080K, ER=.1/1/2eV gives1.191×10⁻¹³ /2.377×10⁻¹⁷ /7.247×10⁻²² per D at fixedA. At ER1eV, T800/1080/1500K gives1.836×10⁻¹⁹ /2.377×10⁻¹⁷ /1.019×10⁻¹⁵. Density enters quadratically, exposure linearly; no second amplitude rescaling is hidden. The source narrow-resonance convention ER1eV is not fixed by these plateau yields.

A useful scale check: uniformly heating1870nm³ of host density6×10²⁸m⁻³ from300 to1080K costs22.6keV with classical3kB heat per host atom alone. A2keV ion could heat at most165nm³ under those particular assumptions and100% deposition. This does not reject a nonuniform transient track or a hot-D subsystem; it means that peakT, volume and lifetime cannot all be read as independent uniformly hot quantities. Supply the reaction-weighted T(r,t), relevant heat capacity, density and deposition kernel before assigning an absolute rate. Pre-existing stored energy must be separately specified.

Retained actual shared fits (`g4_beam/czerski/README.md`) are already more discriminating than an arbitrary hot population: the same frozen-density, fixed-hot-screening model requires UC-Davis/Czerski A ratio **283.8** with one current power, **622.5** as plotted; even separate amplitudes leave material mismatches. In the separate48-point shape convention thermal χ²132.23 vs sigmoid36.92, layer40.29, fitted-neutrals46.31, common-nuclear/no-additional-floor181.45. These are conditional comparisons, not common across every other g7_micro table or statistical exclusions of all track models.

**Overlap:** for independent arrivals, `λ=Φ πr²τ`, `P(extra arrival)=1−exp(−λ)`. Using the author's estimatedΦ=1.1×10¹⁹D m⁻²s⁻¹, r10nm andτ100ps givesλ3.46×10⁻⁷. Radius1/100nm changes it by.01/100; lifetime10fs/1ps/1ns changes it by10⁻⁴/.01/10 relative to100ps. At10nm mean overlap1 requiresτ≈.289ms. An independent illustrative140μA pure-D+ beam uniformly filling a3mm-diameter **assumed spot** givesΦ1.24×10²⁰ andλ3.88×10⁻⁶ at10nm/100ps;85/15 molecular current would increase D-number flux by15%. The3mm extraction aperture is not a measured target spot, so these assumptions are not combined into a false measured intensity. An electronic10fs relaxation window gives another10⁻⁴ reduction. These numbers disfavor ordinary overlapping short-lived nanometre electronic excitations, not a measured long-lived defect/phase state. Sensor live time is not excitation lifetime.

The existing collision40.55eV static sample and.22eV single2keV encounter difference apply at.1Å to one head-on trajectory. The latter is **not time-step checked**; the vacant-path drag check does not supply that check. Frozen neighboring-D displacement changes static lowering by−.687eV, anchored encounter work by+.675eV; these are not extra screening terms. Separately, `g4_beam/screen/dynamic/README.md` gives contact PdD−Pd+13.58eV, matched semicore representation differences.112/.005eV, and low-q speed contrasts−5.55,+10.53,+11.42,+5.03eV at.5/2/5/13keV lab. They are different response constructions—not additive to the pair curve or encounter difference. No computed ordinary bulk closure here supplies the fitted keV scale.

## 4. Beam, analysis and isotope discriminators

The final Karahadian article fits Pd/PdD800±200/1700±200 and Ti/TiD300±200/900±200eV over2.5–6.5keV CM; its≈3keV PdD and2.2keV TiD values at.25CM are per-energy phenomenological inversions, not the same constant fits. Protons establish the lowest-energy floor; neutron evidence overlaps at higher energies. Beam-off electrochemistry gives background. [Final article](https://doi.org/10.1038/s41467-026-74421-1).

The85/15 source composition is an assumption, current comes from a separately measured cupI(V), and the≈5% energy spread was measured in a2.5–4.2kV retarding window only. Postbias at high voltage changes cup focusing; mere charge exchange at the low-voltage impact energy does not create an above-nominal tail. Czerski's≈3×10⁻⁴ pre-deceleration neutral control and D3 suppression are source-specific; Pd without deceleration lacks that particular high-energy-neutral pathway. [Czerskiv3](https://arxiv.org/pdf/2409.02112v3).

The saved direct Figure2/SI1 comparison at.25 versus3keV CM gives Fig/SI1 shape factors Pd4.84,PdD4.68,Ti4.89,TiD4.79, nearI(3)/I(.25)=145/38=3.82 times a common stitch. This is an **observable-normalization discrepancy**, not proof the paper duplicated a current division. SI1 per-current ratios Y(.25)/Y(3)=.181/.246/.038/.039 still retain excess relative to ordinary bare behavior. The retained common-proton uniform fits show the consequence of alternative measurement conventions, not four physical subtractions:

| U in eV |M0 as plotted |M1 one extra1/I reading |M3 1/I plus focus |
|---|---:|---:|---:|
|Pd |1247 |696 |254 |
|PdD |2039 |1400 |805 |
|Ti |977 |452 |31 |
|TiD |1315 |842 |450 |

The prior CSV and script define these alternative readings; no one is selected here without raw count/live-time/true-charge reconstruction. AtM1 changing effective Eq7 to the common nuclear proton model changes these fits by<10eV, not hundreds. A coherent energy offset has yet another signature: locally mapping σ(Eactual,.04keV) to σ(.25,3keV) requires an illustrative≈+7.16kV atomic-D lab error, not a small5% spread; a rare energetic tail must be constrained separately. This is a frozen-S pointwise sensitivity, not a thick-yield calibration bound.

**D→T:** same lab energy gives Ecm,DT≈.600L vs DD=.5L, and the nuclear Gamow/S factors change. For matched thermal distributions the incremental motion U_T/U_D≈.8; matched force constants in the zero-point limit give≈.65. At300K the tied oscillator calculation gives .377/.505≈.747 at.5keV lab, because it is not in the zero-temperature limit. These are incremental motion contrasts, not total fusion yield ratios.

Reuse the actual `g6_plateau/discriminators.csv` counterfactual: total DT/DD events/true C at **.5keV nominal lab** are **2.60 screening-rise,3.18 surface,55.6 fitted fast neutrals,2.54 fixed-lattice background**. They hold the Pd fuel number profile and sheet, fitted law and beam assumptions fixed; they are not measured, not transferred to other beams and do not include a computed isotope-dependent lattice state. DD threshold resonance is not transferred. Lowest-energy DT parametrization is extrapolated. T-specific profile/retention/site response would be required for an actual isotope experiment.

## 5. Ranked next discriminators — stop here

These are decision choices, not newly authorized experiments or jobs.

1. **Close the existing count→live-time→true-charge/reference join, and use the same runs' beam/proton phase-space information.** First reconcile SI1 and Figure2 under one raw Poisson likelihood with exposure/background, cup capture and shared reference/stitch. Do not log-transform sparse or nonpositive background-subtracted counts. Where available, join energy/species/neutral distribution and actual detector angle/absorber response. If the common4.7–4.9 shape factor is normalization, refit all mechanisms to that corrected observable; if it persists in independently delivered charge, an elementary material/beam source remains. Spectra can then separate high-energy beam and near-rest components: for DD at stationary target, proton first-order centroid boost is approximately `sqrt(3020keV * Ecm,trans) cos(theta)`. At illustrative135°, .25keV CM gives−19.4keV, while4–5keV energetic neutrals give−77.7 to−86.9keV, a58–68keV separation versus simulated35–39keVFWHM. **135° is not established for Davis;90° is a first-order blind spot.** SI1's.25keV PdD signal≈.21count/h implies only≈2.5–5 counts in12–24h, so no fabricated single-bin centroid precision: use a joint marked-energy count model if suitable released runs exist. This is the cheapest route because it uses the named normalization inconsistency and already recorded channels before commissioning new material physics.

2. **At one host/channel, separate instantaneous flux response from profile/site history and independently constrain the shallow fuel/chemistry.** Rapid current modulation or energy interleaving at matched cumulative dose, with measured species/energy/charge, distinguishes an immediate additive neutral term or nonlinear excited-state term from slower implantation/depletion/annealing. A fixed neutral flux while charged current halves predicts near2×yield/C; a co-scaling neutral source does not. Ordinary frozen material/site mechanisms predict≈constantY/C. Then matched depth/oxygen-carbon/temperature and total loading evidence decides whether the residual follows nD(z)/stopping or occupied special sites. Existing model near-surface×50, molecular-fraction and half-current ratios are parameterized forecasts, not observations. If the residual tracks measured profile alone, prioritize transport; if it tracks occupied defect capacity at fixed profile, prioritize the existing part(b) barrier/site calculation. Temperature alone is not discriminatory, and no new DFT is started here.

3. **Use a common-exposure channel/isotope contrast only after that observable is fixed.** Existing calibrated proton/neutron/IPC-sensitive ratios can constrain a DD-specific term: an energy-dependent IPC/p residual like the published Zr comparison promotes the channel mechanism, whereas a shared enhancement in unrelated Li/F reactions demands a material/beam component too. The frozen DT/DD forecast55.6 vs2.5–3.2 is a strong conditional separator of neutral-dominated and target-dominated floor models, but an actual T intervention needs the owner's separate authorization and T-specific retention/profile. F+p strength, F+d activation and DD branching must retain their separate forward models; they are not independent estimates of one universal U.

## Reproduction, evidence and delivery boundary

From the checkout root, with the existing scientific Python environment and NumPy/SciPy:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python g7_gap/motion_once.py
```

The script is one reproducible one-off, not a maintained general-purpose API. It writes `results.json` only in this directory. The final narrow consumer/model guard checks the cold identity against the actual shared kernel, equal-mass recoil conservation, the independently derived low-temperature expansion,48/72-node log results, and CSV key/shape coverage. Achieved cold/source relative identities are0 /4.44×10⁻¹⁶, elastic identity error4.44×10⁻¹⁶, representative quadrature log differences2.13×10⁻¹⁴, low-temperature expansion relative difference1.49×10⁻⁷. These are numerical/semantic checks of the declared models, not physical uncertainty enclosures or an experimental likelihood test. No broad suite, old fit, DFT job or campaign was replayed. Inputs carry SHA256 in JSON; code and sources remain graph-independent.

Delivery-ready whitelist: **README.md, candidate_ledger.csv, experiment_rows_addendum.csv, motion_once.py, results.json**. The corrected source working packet may stay private locally; no paper bytes, instruction snapshots, credentials, unrelated artifacts or per-run embargoed data are needed. The script records the packet's identity for local provenance but does not require it for calculations other than that optional hash. Marvin retains Fleet recording/selection and private draft-release delivery authority. There is no claimed publication, upload or independent peer audit of this new numerical model.

The remaining scientific uncertainty is concentrated: the actual normalized count/charge/profile/species consumer; occupied strong-site barrier laws; and the population/temperature/lifetime of any genuinely non-equilibrium state. The ledger and three deciding choices are the requested stopping point.

## Lead's note

Delivered as `g7-gap-9c0fb42.tar.gz` (SHA-256 f1dac110a1647178193f4f535feb1300c08809e7cf7ea8382dd9cf729c0dc3b8, 5 files), unchanged except for this note. The lead reran `motion_once.py` under NumPy 1.26.4; every number agreed with `results.json` to 2e-12.

What it means, in plain terms. All models give about 40 eV and experiments report 300 eV to 3 keV. The ledger tests each possible cause of that difference against the published data:
- **Ruled out as the main cause.**
  - Ordinary motion of the target deuteron (zero-point and thermal) adds under 1 eV.
  - A deuteron knocked on by the same beam ion cannot supply the energy.
  - Electrons excited by one beam ion add 0.2 eV, and the excitations of different ions overlap less than one time in a million.
- **Limited to D-D.** The threshold resonance multiplies the D-D yield about 5 times at the lowest energies. It cannot explain the large values in the lithium and fluorine reactions.
- **Still standing.**
  - The material: strong sites and the near-surface state. The same lithium reaction gives 185 eV in an insulator, 1,280 eV in lithium metal and 3,790 eV in PdLi, so the host matters.
  - The beam and its normalisation. The UC Davis figure and its supplementary table differ by 4.7-4.9 times in low-energy shape, and reading the charge differently moves the fitted screening from 2,039 to 805 eV for PdD.
  - A population of already-energetic deuterons, which would need an energy source nobody has specified.
- **Charge scaling does not decide.** Target motion and a real potential shift both grow roughly in proportion to nuclear charge.

Cheapest next discriminators, in order: rebuild the UC Davis yields from raw counts, live time and true charge when the per-run data open on 2026-12-01; separate instant beam-current effects from slow material history at fixed dose; then use a channel or isotope contrast.
