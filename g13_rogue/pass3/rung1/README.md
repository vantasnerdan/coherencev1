# Whole-foil in-band refinement — unpublished draft

Base: `b0a67c08661ec765c3916f1dad02c47f07cf3bdd`. This additive directory does not edit Pass 2. The single driver is `rungs.py`; rung-0 projection registration is included here as requested.

## Kept / dropped before the first refinement

Kept: the full equal-area 21.269-mm square surrogate of the 24-mm disk; 250-um thickness; free natural conditions on six faces; the identical seeded 8-by-8 rotated cubic columnar grains, pure-Pd tensor/mass baseline and fixed loaded D profile; physical 138-uA, 2.5-keV source and correlated D2+ parent recoil/sign; accepted portable addressed C samplers; causal pre-bin conditioning; beta 0 and .35; and exactly three beam shapes (uniform steady, accepted hot spots steady, accepted hot spots plus inherited stochastic bursts). Source-cell centers/weights remain the accepted 6-by-6 spatial quadrature. The two Gaussian hot-spot sigma values are 1.2 mm (FWHM 2.826 mm), on the inherited 3-mm Gaussian envelope, with gains/centers unchanged. They are assumed profiles, not measured beam structure.

Rung 0 is the accepted P5/P5/P2 space. Rung 1 is its genuinely nested P15/P15/P4 full vector space, integrated separately through all grains. A single global variational operator couples all spatial scales, both directions and every interface; no independent wave packets or local patch replace it. Its thickness polynomial is enriched, not a thin-plate-only propagation model. The physical eigenfrequency interval remains 5–250 kHz: 53 old versus 190 enriched retained modes. Increasing numerical spatial support **does not add a second 1%** inside this same physical frequency band.

The fixed support screen is minimum captured physical-volume L2 mass >=95% for real lateral plane waves, sampled at 37 directions and phases 0, pi/4, pi/2, separately applicable to each vector polarization. It gives lambda0=13.73147 mm; the factor-three target is 4.57716 mm and P15 reaches96.110%. This is an independent basis projection screen, not a proof of accurate free-foil transport for every mode, direction or polarization. The full-vector mass-weighted gradient minimum3.51919 mm and nominal250-kHz label remain separate. A representative near-screen eigenmode shifts -6.58% with squared mass overlap.729 after enrichment; that finite change is retained, not suppressed.

Dropped/not claimed: a measured mount, measured grain texture/PdD tensor, atomistic impact mechanics, unique microscopic launch phases, or continuum-converged transport. The 36-point source quadrature is not a resolved 20-nm impact footprint and may alias a truly fine beam map. Neither frequency enrichment nor spatial enrichment replaces beta by a lattice collision law.

## Source and comparison fixed before propagation

A signed vertical traction has a declared Gaussian lateral sigma20 nm and Gaussian time sigma0.5 ps, at each accepted species' recoil-depth first moment. These are dimensioned assumptions, **not measured lateral size/duration**; the CSDA stopping depth and transit time do not measure either. Its projection is normalized separately over the entire retained physical5–250-kHz band, to the historical quiet mean1% host-recoil self-energy. The same whole-band calibration, source waveform and physical parent marks apply at both spatial resolutions. At fixed parent mark/state, the old frequency-band self-energy is preserved, not each one of the old53 numerical amplitudes. Adding modes cannot multiply source energy by dimension.

The low/middle/high fractions remain `.01/.865 * (.1,1,10)` with quiet probabilities `.25,.74,.01`; beta can increase effective marked share. Coherent work contains the signed pre-existing-field cross term; the collective bin cap uses **all incoming host recoil**, including cross-sector work once. Quiet marked self-energy, realized marked fraction and actual coherent work are distinct quantities. No new disjoint frequency band is propagated in this first in-band refinement.

A one-component nanoscale Gaussian bulk pulse cannot in general reproduce the1% calibration within the whole recoil budget. `pulse_sensitivities.csv` states this explicitly rather than spreading1% across more frequencies. The permitted interpretation is a separate calibrated low-frequency effective component plus unknown high-frequency component; a bulk Gaussian illustration supplies dimensioned high-spectrum sensitivity, not a measured Lamb source. Every future added disjoint band would require its own nonnegative whole-source share with total <=host recoil, using the same fixed law. No higher-frequency fraction is silently inferred from the1% calibration.

The first comparison uses identical retained lifetime/range rules at both resolutions, then **substitutes**, rather than adds, the cited Cu mode-frequency loss proxy. The latter has zero separately declared mount loss. It is not calibrated Pd/PdD damping and is not added to a second grain-scattering law. The retained50-us/.2-m rule is a historical sensitivity, not rebranded as a material measurement.

Both loss comparisons start at zero, use the same235.2-us physical burn and100-us sample, 50-ns bins, and independent record IDs0,1. There is no claim this burn equilibrates the long-lived Cu-proxy modes. Matched total/time/site/species/sign addressed streams are preserved; feedback mark probabilities may legitimately diverge when pre-impact fields change. Equal totals alone are not proof of matching marks; draw digests are retained.

The frequency label is the retained **eigenfrequency interval**, not a hard Fourier filter on impulse accelerations or flux. Source ordering remains frozen within50-ns bins; sampling Nyquist is10MHz. There is no implicit resampling of the physical0.5-ps waveform as a resolved collision.

The predeclared height statistic is the maximum over time and **three fixed central2.5-mm-square x50-um-depth regions**, at the original depth slices0–50,25–75,175–225um. Each height is divided by that region's own positive time-mean elastic density, identically defined at both rungs. Exact polynomial cubature is split at grain boundaries. Local point energies at30,50,200um are separate. Spatial129-by129 maps reconstruct the same frame only to measure connected half-height widths; they do not turn an unobserved whole-domain space-time maximum into a measured one. Durations use adjacent above-half-height samples, with boundary censoring.

The width component is selected around a peak **inside the same central physical region**, not an unrelated remote maximum in that frame. Raw mapped widths are censored if they reach the foil edge or have fewer than two map intervals; a separate conservative band-support flag marks widths below half the rung's projection-screen wavelength. Flagged values are diagnostics, not estimates of a physically resolved hotspot width. Grain-dependent elastic density can also have sharp synthetic material contrast without a comparably narrow displacement wave.

A: impact on prior disturbance (>.25 canonical rms in strain or displacement), denominator physical charge-parent impacts. C: impact on prior signed +/-3 causal-rms strain/displacement, with overlapping unions/intersections retained. B: distinct-origin sectors with inward, separated-face flux and constructive energy at the accepted off-impact receiver; its denominator is receiver seconds, not parents. Standalone positive crests and negative holes are scored independently in the fixed physical regions. No extrapolated waiting-time law is inferred.

Exact modal damping is energy-budgeted. Conditional D conversion uses the **positive nonlocal quadratic** density formed by applying sqrt(modal energy-loss rate) to displacement and velocity before their elastic/kinetic forms. Its full-domain integral equals modal energy loss, but it is not measured local absorption or a measured D yield. Unit transfer and4.1-ps holding give conditional eV/D; no threshold claim is borrowed. Fixed n_D(z), physical point versus volume-averaged populations and elastic energy versus converted energy remain separate.

For conditional storage, the sampled production is treated as constant over the4.1-ps holding interval; production-rate times lifetime is the steady constant-production solution, not a time-resolved optical calculation. `leading_event_depth_energy.csv` evaluates **all three depths at the same selected event time**, with local and region-averaged populations/energies separate.

## Observed comparison and first support stop

The two matched records show a substantial **absolute** energy change after in-band enrichment, but not a loss-robust relative-tallness direction. Across the requested beam/beta/depth comparisons, the mean elastic-density rung1/rung0 ratio is4.47–15.98 and the ratio of observed maximum elastic eV/D is5.14–15.05. These are ranges of paired finite-record averages across arms, not population uncertainty bounds or a higher-band energy-yield measurement.

| 30-um region, beam / beta | Peak-to-own-mean ratio: rung1/rung0, retained loss | Same statistic, Cu proxy |
|---|---:|---:|
| Uniform steady /0 |0.978 |1.290 |
| Hot spots steady /.35 |1.378 |0.868 |
| Hot spots with bursts /.35 |1.286 |0.933 |

The direction reverses with the assumed damping in10 of18 beam/beta/depth comparisons. The representative same-band fields can become more concentrated in absolute energy without demonstrating a growing relative rogue tail. Two independent records do **not** resolve a population tail or a beam-shape winner; pointwise paired-record summaries preserve those record units and do not count ions/time bins as independent samples.

Width maps contract, but14 of24 enriched-event maps are below the conservative band-width screen;20 of24 old-event maps reach a foil boundary. Those censored widths do not establish a quantitative narrowing trend. Leading causes remain **unattributed**: signed cross-sector elastic work, prior A/C impact fractions and the separate remote B window are delivered as overlapping evidence, not forced double-impact/crest/hole assignments. Aggregation retains distinct source sectors but not unique within-bin parent trajectories.

At the same selected leading-event times, rung1 with retained loss gives the following **local** elastic/conditional-storage ranges across these finite beam/beta/record cases (eV/D, not uncertainty bounds):

| Depth | Local elastic eV/D | Conditional converted eV/D |
|---|---:|---:|
|30um |6.89e-10–6.30e-8 |3.59e-16–6.36e-15 |
|50um |1.90e-11–1.41e-9 |1.23e-17–1.47e-16 |
|200um |1.26e-11–6.69e-10 |6.85e-18–7.22e-17 |

These values consume the fixed local n_D, not a borrowed full-loading threshold. The Cu-proxy rows and separate volume averages are in `leading_event_depth_energy.csv`; neither model measures D conversion.

The common-band control matches all physical parent totals, all quiet addressed mark digests, all selected accepted C counts, and10 of12 accepted-control mark digests. In two hot-spot beta.35 record0 comparisons, the finite20-nm waveform/recomputed projection perturbs feedback marks: A parent counts differ by270,688,725 and464,607,821 (fractions0.00314 and0.00539 of the approximately86.13-billion parents); B occupancy differs by0.005 and0.0035. This is a disclosed changed-source/feedback control, not an exact replay claim. `accepted_consumer_comparison.csv` gives every row. Maximum source/propagation energy residual is1.56e-15 of incoming host energy. Realized quiet marked fractions are0.009982–0.010000; feedback fractions0.02250–0.03970. The Cu-proxy C counts depend strongly on the long-lived transient and its causal-rms initialization and are not stationary rare-event rates.

**Stop before rung2:** the next shared-screen target is1.52572mm, with kh1.03 and less than one old source cell per wavelength. A matrix-free global spectral or sparse hp representation could extend numerical support, and source quadrature could be hierarchically refined; neither supplies the missing physical damping. The already observed sign change between permissible loss assumptions is a **physical support boundary**, not a request for more solver precision or another low-band scan. No automatic rung2 solve is performed.

The first required measurement is **frequency- and branch-resolved guided Pd/PdD attenuation/lifetime at the relevant loading/temperature, with mount/edge loss separated**, across the represented5–250-kHz band and the branches entering the next shorter-wave band. Intrinsic absorption must be separated from grain redistribution before combining it with explicit grains. The acoustic launch's lateral spatial spectrum, duration and branch partition are a separate required source measurement; stopping depth/transit time cannot substitute. This does not require opening a micro/MD model.

## Reproduction and payload

From repository root, with NumPy/SciPy and C99 compiler for live draws:

```sh
export OPENBLAS_NUM_THREADS=5 OMP_NUM_THREADS=5 MKL_NUM_THREADS=5 PYTHONDONTWRITEBYTECODE=1
python g13_rogue/pass3/rung1/rungs.py compare
python g13_rogue/pass3/rung1/rungs.py readout
python g13_rogue/pass3/rung1/rungs.py package
```

The delivered `model.npz` freezes enriched eigenvectors/eigenfrequencies. `consumer_r0.npz` and `consumer_r1.npz` freeze source/projector/covariance/receiver/energy maps, bound to those exact mode bytes, so replay `compare` does not rebuild numerical inputs. The map files were reconstructed after the successful propagation from the same frozen modes and routines/runtime, without a second propagation; `numeric_consumers.json` states this provenance. `prepare` regenerates an operator and `freeze` generates numerical consumers as a separate input epoch; neither is required for the delivered replay. `readout` is compiler-free. Accepted Pass2 portability evidence remains the base acceptance, not a claim of a new cross-runtime replay. This draft has no commit, push, publication or claim promotion.

Payload guide:

- `projection_support.csv`, `mode_support.csv`, `support.json`: independent direction/phase screen, full-vector/polarization diagnostics, nominal and gradient labels, nested-operator and representative dispersion evidence.
- `records.csv`, `class_summary.csv`, `signed.csv`: A/C physical-parent counts/fractions/rates and overlaps; B receiver occupancy/connected rates; standalone signed crests/holes; independent-record uncertainty.
- `tails.csv`, `histories.npz`, `trend.csv`, `peaks.csv`: full sampled energy distributions, fixed-region statistics, adjacent durations, actual connected same-frame maps and all censor flags.
- `leading_event_depth_energy.csv`: same-event local/volume-average elastic and separately conditional converted eV/D at30/50/200um.
- `energy_ledger.csv`, `consumer_verification.json`, `accepted_consumer_comparison.csv`: actual recoil/work/loss budgets and exact versus nonexact common-band control observations.
- `pulse_sensitivities.csv`, `loss_law.csv`, `inputs_ledger.md`: model source fractions, incompatible single-pulse diagnostic, explicitly transferred/extrapolated loss and publication/measurement citations.
- `manifest.json`: exact root-relative payload SHA-256/byte manifest; `capsule.tar.gz` contains only regular files below this new directory.
