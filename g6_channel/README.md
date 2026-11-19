# g6_channel: resolved axial paths, encounters and drag

## Result to carry forward

**The best-supported gain candidate in the low-energy stopping domain is Pd(D/T), `<110>`,13keV lab,2deg polar RMS incidence at300K: about1.3–1.4x the random energy return in this coarse model.** The default256-particle D-T result is **1.420±0.060** (one particle standard error), with absolute nuclear energy return **1.886e-6**. D-D gives **1.341±0.055**, return **3.661e-9**. These are conditional model results, **not a validated material optimum**: the independent angular-scale check fails, and the electron-density sensitivity is larger than the particle error. No sub-degree collimation benefit is established; among the requested divergences,2deg gives the largest refined mean by admitting more fuel encounters.

At that D-T point the path-weighted fuel-site flux is1.198x random, total drag per path0.8625x random, and electronic drag0.9025x random. The1.420x yield is the actual joint integral, not a product of these means. Initial capture is91.4%; every captured sample eventually exits before the10eV terminal cutoff, with mean first-exit path181.5nm. **Initial capture is not the dominant restriction; channel survival/dechanneling and useful fuel overlap are.** The reported exit probability includes slowing-induced closure of the continuum allowed region, not only a thermal or fuel kick. It must not be read as a measured thermal dechanneling probability. The2keV best refined D-T point is1.411±0.083 at2deg, return5.110e-14 and mean first exit37.8nm; its gain is statistically tied with13keV, but its absolute return is much smaller.

The original32-particle full-matrix nominal absolute leader is Ti(D/T)2 `<100>`,200keV,2deg: D-T return0.002231 and channel/random factor1.114±0.078. It is **not distinguished** from the uncalibrated LiT random result0.002200, and both depend on unqualified high-energy electronic stopping. No global crystal/energy optimum follows. At that nominal200keV leader, another448x in nuclear return would be the unity-return specification; at the refined13keV PdT point it is5.30e5x. Geometry supplies an order-unity lever here, not that missing factor.

### Refined Pd `<110>` table

Particle standard errors, not confidence bounds or model errors:

| Fuel target | Lab keV | RMS divergence | Yield / random | Energy return | First-exit mean nm |
|---|---:|---:|---:|---:|---:|
| D | 2 | 0deg | 1.163±0.074 | 1.027e-16 | 46.2 |
| D | 2 | 0.5deg | 1.373±0.100 | 1.208e-16 | 44.6 |
| D | 2 | 2deg | 1.388±0.078 | 1.232e-16 | 37.2 |
| T | 2 | 0deg | 0.990±0.078 | 3.617e-14 | 46.1 |
| T | 2 | 0.5deg | 0.992±0.070 | 3.625e-14 | 45.4 |
| T | 2 | 2deg | 1.411±0.083 | 5.110e-14 | 37.8 |
| D | 13 | 0deg | 1.059±0.052 | 2.886e-9 | 203.5 |
| D | 13 | 0.5deg | 1.292±0.063 | 3.529e-9 | 199.4 |
| D | 13 | 2deg | 1.341±0.055 | 3.661e-9 | 178.1 |
| T | 13 | 0deg | 1.035±0.058 | 1.373e-6 | 205.5 |
| T | 13 | 0.5deg | 1.152±0.056 | 1.531e-6 | 200.8 |
| T | 13 | 2deg | 1.420±0.060 | 1.886e-6 | 181.5 |

`results.csv` is the complete **288-row current table**, including every material, axis/random, energy, divergence and isotope. It replaces only the Pd random/1102/13keV rows with these256-particle refinements; all other rows retain32 particles and their reported errors. Ratios are recomputed against the matching current random row, with `source_stage` and `random_reference_stage` explicit. `results.json` includes that table **and every curated raw stage**, so no normalized result hides its denominator.

### Dominant uncertainty and stopping point

Changing all superposed-atom electronic radial lengths to0.7x and1.5x, without changing the nuclear geometry or fitting a normalization, moves the13keV2deg D-T gain from **1.768±0.070 to1.322±0.058**. D-D moves from1.628±0.060 to1.188±0.043. This is a sensitivity range, not a justified physical enclosure. It identifies the bonded/local stopping density as the deciding missing input; another indiscriminate trajectory sweep would not settle it.

For the **actual best-case D-T consumer**,512-particle continuations at paraxial exit thresholds0.1 and0.3rad give gains1.3177±0.0380 and1.3260±0.0405, respectively: an observed **0.63% change**, below their roughly3% particle resolution. First-exit means are182.4/183.6nm. This supports the hybrid for this yield consumer at the achieved precision, not a1% global error certificate or a validation of every high-energy row. The alternative ensembles also show why the256-particle default1.420 should be summarized as about1.3–1.4, not an exactly established optimum.

The independent angular-scale comparison is **unmet**: across the original36 zero-divergence D crystal/axis/energy points the resolved critical-angle/Lindhard-scale ratio is0.120–0.706, with0/36 within20%. At the refined Pd11013keV point the angles are5.096deg versus10.763deg; at2keV8.449deg versus27.440deg. The observed capture fractions are separate trajectory observables. This disagreement leaves the continuum acceptance/dechannel model quantitatively uncalibrated; it is not repaired by defining one angle to equal the other. The current evidence supports the named conditional lever and the required density/acceptance input, not promotion to a validated experimental gain.

The 0/36 comparison above uses the characteristic scale \(\psi_1\), not a finite-energy acceptance angle. Andersen's Eqs. (28) and (33) give a separate regime check: solve \(\psi^2=\tfrac12\psi_1^2\log[1+(Ca/d\psi)^2]\) and take the smaller of that root and \(\psi_1\sqrt{\tfrac12\log[1+(Ca/\rho)^2]}\), with \(C=\sqrt3\), \(a=0.8853a_0 Z^{-1/3}\) and \(\rho=\sqrt2\,u_{\rm host}\). Using the delivered host spacings and displacement widths, the Pd `<110>` estimates are 9.08° and 5.47° at 2 and 13 keV lab, versus 8.45° and 5.10° from the projected potential (7% lower). The two estimates are within 20% for 4/9, 3/9, 1/9 and 0/9 of the delivered D-axis cases at 2, 13, 50 and 200 keV. The comparison is to a single-string continuum estimate, not a measured acceptance or a pass of the original all-case check. It supports the angular scale at the named low-energy Pd candidate but not across materials or at higher energies.

All runs were serial on one CPU with NumPy2.5.4/SciPy1.18.1. The full matrix took1209.33s; matched references13.21s; the four retained256-particle refinement/density runs232.96s; the actual-best-case hybrid comparison58.44s. Curated completed stages total1513.94s, peak RSS637.12MiB. The earlier deciding pilot took213.37s; its legacy output is not a consumer. No long sweep, DFT, build, lint or formatter was run.


This directory replaces independent mean stopping and mean fuel density with their **joint path integral**. It is a finite, coarse superposed-atom model, not a new DFT calculation. The governing node is `g6_channel`, definition1, accepted handoff `entry_2038`; the shared plan and pinned base are in `inputs.json`.

## Model and input meanings

The conventional cells contain four fcc host sites. PdD and LiD add four octahedral fuel sites (rocksalt); TiD2 adds all eight tetrahedral sites (fluorite). The orthogonal transverse/longitudinal supercells are `(010,001,100)`, `(1-10,001,110)` and `(1-10,11-2,111)`. Explicit enumeration counts projected rows and their line densities. These are **PdD1 and TiD2**, not the historical transport fills PdD0.70 and TiD1.70. The target is a uniform semi-infinite crystal with a flat ideal entry surface. No adsorbed sheet or transported loading gradient is silently imported.

Each atom supplies a normalized exponential effective electronic density `exp(-r/b)/(8*pi*b^3)`. The Fourier projection is `(1+b^2*k_perp^2)^-2`; convolution with the positional Gaussian supplies its Debye-Waller factor. Atomic lengths, lattice constants and their source/assumption status are in the inputs ledger. The host lengths are coarse Slater-scale estimates, not measured or PAW valence densities. Every species' local stopping multiplier has unit cell mean1; multiplying it by that species' actual bulk electronic stopping coefficient preserves the random-medium normalization without setting a yield ratio by definition.

The projected ZBL string potential is calculated from the same rows. In Fourier space each screened-Coulomb term contributes `4*pi*Z*e2*n(k)/(k^2+beta^2/a_ZBL^2)`. Its gradient steers the projectile. This **conservative confinement force is not drag**. The electronic loss is `sum_j n_j*Se_j(E)*rho_e,j(x,y) ds`. Independent nuclear kicks use the thermally broadened local nuclear density, a random impact disk and the existing ZBL scattering kernel. Every angular kick has its own kinematically associated recoil loss, separated into host and fuel columns; there is no additional mean nuclear stopping added on top. The continuum-plus-incoherent-kick split is an approximation, not an exact lattice MD decomposition. It omits correlations between successive displaced atoms and treats 300K Gaussian displacements as a stationary ensemble, not an anharmonic time-resolved phonon calculation.

The fusion probability is accumulated on **that same** `(E,x,y,z,ds)` path as `n_fuel*rho_n,fuel(x,y)*sigma(Ecm) ds`. Thus suppressing stopping while avoiding fuel sites does not automatically improve the yield. Reported fuel flux is the path-length-weighted local fuel density divided by its random density; the yield is the integral of the product, **not** the product of separately averaged flux and inverse stopping.

### Capture, channel exit and random continuation

The independent Lindhard reference is `psi1=sqrt(2*Z*e2/(E*d_hostrow))`. A different, geometry-resolving continuum criterion solves the isolated-string curvature condition `U_string''(rc)=8 E/d^2`, thermally enlarges `rc` to `sqrt(rc^2+(2u_host)^2)`, and finds the lowest total projected potential on that excluded-row boundary. Initial transverse phase-space points below that barrier are counted as captured. This gives a potential/curvature angle **and an observed finite-sample capture fraction**, not an angle defined to equal the Lindhard estimate. The angle comparison is only a scale comparison; low-energy closure and neighboring/fuel rows can make it disagree substantially.

First dechanneling is the first crossing of that barrier by an initially captured trajectory. It is not the first stochastic collision. The code continues the particle after this event: no active trajectory is paused to save runtime. A trajectory changes to homogeneous stochastic transport only when the allowed continuum region closes as it slows, or when it is above the barrier and its transverse direction exceeds the explicitly stated paraxial switching angle. That hybrid matching is a model assumption whose effect must be reported, not a proved instant loss of crystalline memory. The terminal boundary is an actual surface escape or10eV lab cutoff. A particle still captured at its terminal boundary contributes its **actual observed path** to `survival_restricted_mean_nm`; `captured_stopped_before_exit` states the count. The restricted mean is not an invented uncensored dechanneling length. `dechannel_length_observed_mean_nm` averages actual exits only.

The random reference is a dynamic, scattering trajectory ensemble with uniform incidence, not a CSDA integral relabeled random and not uniform angular incidence over the whole sphere. "Uniform incidence" means uniform transverse impact area and the same beam angular distribution. Random orientation removes the axial periodic field and local multipliers while retaining discrete angular scattering, associated nuclear recoil, electronic loss, surface reflection, chemistry and nuclear channels.

### Reaction and energy conventions

The beam remains deuterons in every row. Replacing target D by T changes recoil mass, harmonic displacement width and `Ecm=Mtarget/(MD+Mtarget)*Elab`; it does **not** change the lab velocity to a triton velocity. D-D keeps the existing forward model's proton-channel threshold-resonance term and Bosch-Hale neutron background; D-T uses the separate Bosch-Hale `D(t,n)alpha` rational fit, not a copy of that D-D resonance. `yield_primary` is DDp or DTnalpha, respectively; `yield_secondary` is DDn or zero. Energy return sums the isolated reaction-channel probabilities times their own Q values and divides by incident D lab energy. It is nuclear energy released, not a statement of harvestable electric efficiency. Screening is held at `Ue=0` for every numerator and denominator to isolate geometry; no fitted keV screening is imported.

The stochastic electronic stopping is the actual `g4_beam/transport/bca.py` Lindhard-Scharff implementation. Pd/Ti use its stated1.54/1.52 factors and H its0.75 factor. Li uses that implementation's unscaled LS default and ZBL nuclear force: **there was no existing calibrated Li target to reuse**. The Pd/Ti factors are low-energy estimates around5-10keV/u. Their50 and200keV lab predictions are uncalibrated velocity-law extrapolations; Li stopping is uncalibrated at every energy. The D-T fit is tabulated for0.5–550keV CM; the common continuation to10eV lab is an explicit sub-fit-domain extrapolation.

The scattering targets are stationary during each binary encounter: thermal **positions** and zero-point positional width generate the incoherent kicks; target thermal velocities are not separately sampled. Their leading speed scale is about1km/s versus438km/s for a2keV D beam, so this is a small relative-velocity approximation, not an omitted keV kick source. Collective target motion and correlated electronic response belong to the separate coupling/collision nodes. The paraxial conservative equation holds longitudinal energy fixed to leading order in `U/E` and angle; it is not an exact nonparaxial energy-conserving atomistic trajectory.

## Reproducibility and outputs

Use NumPy/SciPy on one CPU, with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`. The script imports only the pinned repository's `g4_beam/transport/bca.py` and `g4_beam/forward/model.py`; it does not build BCA caches, launch DFT, or alter transport profiles. Output JSON records source hash, RNG construction, particle count, spatial grid, maximum step, measured runtime and peak RSS. CSV contains the same scientific rows for downstream consumption. Each particle has its own adaptive step constrained by fractional electronic energy loss, transverse spatial resolution and an exponential optical-depth clock for each nuclear species. Homogeneous continuation is event-driven; crystalline integration additionally resolves the varying local densities. There is no fabricated end length on a budget expiration: failure to finish raises and emits no successful output.

The scope is three crystals, three axes plus random, lab energies2/13/50/200keV, polar RMS divergences0/.5/2deg,300K, D and target-T mass variants. Strained/layered crystals are not added. The three biggest missing pieces are electronic bonding/charge-transfer density, beam damage/evolving occupancy, and entry-surface disorder/chemistry. Thermally correlated row scattering and accuracy of the paraxial-to-random transition are additional explicit finite-model limitations.

## Deciding stages and independent reference

The first one-CPU TiD2 pilot used12 particles per selected row,96x96 transverse grids and maximum step0.08nm. It completed all8 selected rows in213.37s, with637.12MiB peak RSS. At13/200keV, all three aligned axes suppressed the mean yield relative to the dynamic random reference in this small pilot: the deciding issue was **fuel-site avoidance**, not an assumed stopping gain. Its broad particle error bars were retained; it was not used to claim an optimum. This justified completing the specified matrix cheaply rather than launching a long high-statistics trajectory sweep. Production uses an event-driven homogeneous continuation instead of the pilot's small Bernoulli steps; that optimization does not freeze surviving channels.

`references.json` compares512-particle independent dynamic random ensembles against **actual calls** to the pinned `g4_beam/transport/bca.py::run`, at2/13keV and historical homogeneous fills PdD0.70/TiD1.70. Both routes use the same species densities, uniform profiles,10eV cutoff, half-space escape convention and stopping extrapolation. The BCA range comparison is the implanted-only projected mean on both sides, not an unconditional endpoint mixed with a surviving-ion mean. Ratios of those ranges are1.0428/1.0535 for PdD0.70 and1.0185/1.0913 for TiD1.70. This is agreement within the declared coarse10% range tolerance, not a new stopping calibration. Poisson free paths versus the BCA driver's fixed mean free path remain a real transport-model difference.

The separately integrated existing forward `Target`, with exactly those **uniform** fills, the same two isolated D-D channels, `reaction='proton'`, `Ue=0` and10eV cutoff, gives dynamic/forward yield ratios1.0540/1.0347 (Pd) and1.1072/1.0737 (Ti). Thus the stochastic reference reproduces that thick-target comparison within the declared15% coarse tolerance; the largest observed discrepancy is10.72%, not zero. This does not reproduce the transported multilayer/adsorbate/fitted-screening measurement model, which is a different target. No Li agreement is claimed in the absence of calibrated Li range data.

The Lindhard angle reference is independently sourced from [Lindhard (1965), *Influence of Crystal Lattice on Motion of Energetic Charged Particles*](https://gymarkiv.sdu.dk/MFM/kdvs/mfm%2030-39/mfm-34-14.pdf) and [Andersen, *Notes on Channeling* (2018), Eq29](https://phys.au.dk/fileadmin/site_files/publikationer/Lecture_notes/Channeling_notes_2018.pdf). That reference is a characteristic string scale, not an exact acceptance angle in a multicomponent low-energy hydride. The curvature and thermal excluded-radius construction is a separate continuum model. A20% critical-angle agreement target is an honest comparison to report as **unmet** where it fails, not a reason to redefine either angle. The measured capture and dechanneling columns provide the observable consequence of that failure.


## Commands

From the repository root, with a Python environment containing NumPy and SciPy:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python g6_channel/channel.py --stage references --nions 512
python g6_channel/channel.py --stage production --nions 32 --grid 128 --step .05
```

The source's `--material`, `--axis`, `--isotope`, `--energy` and `--divergence` selectors keep a selected aligned case **and its independently simulated matched random reference**. `--density-scale`, `--vibration-scale` and `--switch-angle` expose the stated coarse-model uncertainty rather than hide it in a fitted normalization. Energies passed to `--energy` are lab keV; `--switch-angle` is radians. The unselected default production invocation covers all288 requested rows.

## Curated delivery and regeneration

The shareable set is exactly `channel.py`, `inputs.json`, `results.json`, `results.csv`, and this `README.md`. Intermediate per-stage files are not required for delivery. `results.json.raw_stages` preserves the selected runs with their full CLI settings and the common producing-source SHA256.

The focused stages are reproduced by:

```sh
python g6_channel/channel.py --stage refine --material PdD --axis 110 --energy 2 --nions 256 --grid 192 --step .025 --output g6_channel/refined-2.json
python g6_channel/channel.py --stage refine --material PdD --axis 110 --energy 13 --nions 256 --grid 192 --step .025 --output g6_channel/refined-13.json
python g6_channel/channel.py --stage refine --material PdD --axis 110 --energy 13 --nions 256 --grid 192 --step .025 --density-scale .7 --output g6_channel/density-narrow.json
python g6_channel/channel.py --stage refine --material PdD --axis 110 --energy 13 --nions 256 --grid 192 --step .025 --density-scale 1.5 --output g6_channel/density-broad.json
python g6_channel/channel.py --stage refine --material PdD --axis 110 --energy 13 --isotope T --divergence 2 --nions 512 --grid 192 --step .025 --switch-angle .1 --output g6_channel/best-exit-early.json
python g6_channel/channel.py --stage refine --material PdD --axis 110 --energy 13 --isotope T --divergence 2 --nions 512 --grid 192 --step .025 --switch-angle .3 --output g6_channel/best-exit-late.json
```

For the current table, key rows by `(material,isotope,axis,E_lab_keV,divergence_rms_deg)`, start from production, replace matching rows with `refined-2` then `refined-13`, and recompute each ratio and delta-method error against the matching random row. Alternative density/exit stages are sensitivities, not replacements. This is the complete aggregation rule saved in `results.json.assembly`; there is no fitted or hidden postprocessing. Scientific rows are deterministic for the saved settings/seeds; elapsed time and peak RSS are not byte-reproducible quantities.

