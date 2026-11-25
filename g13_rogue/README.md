# Continuously driven PdD waves: Pass 1

## Fixed-profile addendum: where the model sits

**The retained Pass 1 lattice, seed, stochastic outputs and contact replay describe a loaded-surface scenario absent from the fixed profile.** They place the full-loading 32-D cell (58.7223 D/nm³, `x=n_D/n_Pd=1`) at the beam track. They are preserved, not silently reinterpreted as a calculation in the transported foil. The following depth application reads the unchanged `g4_beam/merit/profiles.json` piecewise layers and accepted G11 stopping archive. The foil is **250 µm** thick; even its deepest PdD layers reach only `x≈0.61`, never the `x=1` calibration.

| Material / CM energy (keV) | x at 30 µm (12% of foil) | x at 50 µm (20%) | x at 200 µm (80%) |
|---|---:|---:|---:|
| PdD / 0.25 | 0.0115374 | 0.349625 | 0.610243 |
| PdD / 2.5 | 0.0107212 | 0.318560 | 0.608777 |
| PdD / 6.5 | 0.00971612 | 0.199051 | 0.603822 |
| Pd / 0.25 | 3.21731e-5 | 3.21731e-5 | 3.21731e-5 |
| Pd / 2.5 | 3.08955e-5 | 3.08955e-5 | 3.08955e-5 |
| Pd / 6.5 | 2.68888e-5 | 2.68888e-5 | 2.68888e-5 |

**30 µm is a dilute transition point, not a fully loaded region.** The samples are layer values, not interpolated loading. At 0.25 keV the PdD central D+ track ends at **23.5032 nm**: its traversed column `∫n_D dz/∫n_Pd dz=8.30593e-4`, with layer minimum/maximum `4.15309e-5/0.0136249`. At 2.5 keV it ends at **128.944 nm**, with column ratio **0.0116457** and minimum/maximum **7.78584e-5/0.0192581**: the beam-enriched mobile skin extends farther. Thus “about 1e-4 from 2–150 nm” is not an energy-independent description. These depths are only 0.00940%/0.0516% of the foil. `g13_rogue/output/depth_profile.csv` includes Pd and PdD point samples, containing-layer boundaries, actual SI densities and D+ column summaries; its 6.5-keV profile is exact, not a recomputed stopping track.

The stopping source is `g11_cmns/collective/part1/motion.py:stopping_split`: direct D recoil is `n_D*sn_zbl(e,1,M_D)` with SI densities (and the source cross-section conversion); it imposes **no recoil-transfer-energy cut**. Tracks stop at **100 eV incident energy**. No stopping split was rerun. The small direct-D channel therefore reflects this fixed track composition, not an absent D–D stopping term.

### Host waves reaching depth

The depth model keeps a thin near-surface source (24–129 nm versus the nearest 30-µm receiver), a uniform 1-cm² disk, its centerline receiver, isotropic downward-hemisphere emission, speed `v=4 km/s`, and energy survival `exp(-s/L_E)`. `L_E=5,30,50,100,250 µm` is a named **acoustic** attenuation sensitivity, not an optical lifetime measurement. Optical amplitude/energy lifetimes remain **8.2/4.1 ps**. No coherent phase preparation, acoustic anisotropy or measured host-to-optical conversion is supplied.

The source uses the accepted PdD host-recoil **90.7856/631.325 J per true C** at 0.25/2.5 keV separately from the chart's **38/138 µA true current**. Numerically J/C is eV per elementary charge: with `F_q=I_true/(e A)`, the launch power per area is `P=F_q E_host eta_host`. The 85/15 ion mixture is already in `E_host`; **do not multiply that power by 1.15 again**. Nucleus counts use `F_n=1.15 F_q` separately. The examples use one shared source fraction `eta_host=0.1` and one local conversion fraction `f_host→D=0.1`; they are sensitivities, not fitted efficiencies. All levels and active counts scale linearly with current; energy also scales with the source fraction, but impact count does not.

For `s_max=sqrt(z²+R_spot²)`, the scalar wave energy density is

`u(z) = P/v [E1(z/L_E)-E1(s_max/L_E)] = P/v ∫_z^s_max exp(-s/L_E) ds/s`.

It follows from each launch's `1/(2πs²)` hemisphere spreading and annular source area `2πr dr`. An extra angular cosine belongs to **plane-normal flux**, not this local scalar density. `u/n_D(z)` is acoustic energy divided by the local D inventory; it is neither coherent D-mode energy nor a bound on phase alignment. The CSV also gives division by the counterfactual `x=1` density to expose the denominator dependence.

The causal count alongside it is the declared survival-weighted **source-age/residence** model:

`N_n(z) = ∫_disk F_n ∫_(s/v)^∞ exp(-vt/L_E) dt dA`

`= (2π F_n L_E²/v)[(z+L_E)exp(-z/L_E)-(s_max+L_E)exp(-s_max/L_E)]`.

It counts weighted old launches whose wave has had time to arrive, not a finite acoustic packet duration, a number of modes or overlapping reaction events. Minimum direct arrival times at 30/50/200 µm are **7.5/12.5/50 ns**. The acoustic energy residence `L_E/v` is separate again (62.5 ns for 250 µm).

| CM keV | Depth µm | Active nucleus impacts, L_E=250 µm | Acoustic eV/local D | Local-conversion optical eV/D |
|---|---:|---:|---:|---:|
| 0.25 | 30 | 6.650e4 | 1.141e-11 | 7.484e-17 |
| 0.25 | 50 | 6.577e4 | 2.923e-13 | 1.918e-18 |
| 0.25 | 200 | 5.414e4 | 4.423e-14 | 2.901e-19 |
| 2.5 | 30 | 2.415e5 | 3.100e-10 | 2.034e-15 |
| 2.5 | 50 | 2.389e5 | 8.064e-12 | 5.290e-17 |
| 2.5 | 200 | 1.966e5 | 1.119e-12 | 7.343e-18 |

The last column uses **only local conversion of attenuated power**:

`E_opt,local/D = (v/L_E) u(z) f_host→D tau_E/n_D(z)`, with `tau_E=4.1 ps`.

At 2.5 keV and 30 µm, it is **3.166e-12** times the full-chart coherent-carrier reference **0.6423 meV/D**, and **1.997e-13** times the full-chart 300-K reference **10.1847 meV/D per polarization**. At 50 µm these ratios are **8.236e-14/5.194e-15**, and at 200 µm **1.143e-14/7.210e-16**. These are reference comparisons, **not thresholds calibrated at `x≈0.01–0.61`**, nor a probability that incoherent energy forms a carrier. The acoustic inventory must not be substituted for the optical column.

Attenuation changes the count and the stored energy differently. At **2.5 keV, 50 µm**, the same source gives:

| L_E µm | Active nucleus impacts | Acoustic eV/local D | Local optical eV/D |
|---:|---:|---:|---:|
| 5 | 9.713e-4 | 2.742e-17 | 8.993e-21 |
| 30 | 2.116e2 | 5.167e-13 | 2.824e-17 |
| 50 | 1.431e3 | 1.447e-12 | 4.746e-17 |
| 100 | 1.416e4 | 3.692e-12 | 6.055e-17 |
| 250 | 2.389e5 | 8.064e-12 | 5.290e-17 |

Longer acoustic survival increases the acoustic inventory; local loss power also contains `v/L_E`, so the converted optical level need not increase monotonically. `g13_rogue/output/depth_waves.csv` retains every depth/length combination through **250 µm**, with the two true-current/source values separate.

**First reflected return, separate high-reflectivity scenario:** a specular rear-face image replaces vertical distance `z` by `2H-z`, with **energy reflectivity 0.9** multiplying that path once. For 2.5 keV and `L_E=250 µm`, direct plus first echo gives acoustic levels **3.197e-10, 8.448e-12, 1.633e-12 eV/local D** at 30/50/200 µm, versus the direct values above. Echo/direct energy ratios are **0.03135/0.04764/0.4590**. At the rear face the interior limiting sum is 1.9 times direct. These are returns of the **same launches**, not independent extra impacts; separate weighted echo counts are retained, not added as an independent Poisson population. Only this one rear reflection is included; no front-face or repeated-reflection series or interference phase is assumed.

### What replaces the D band in the track

In the dilute-track approximation, D atoms have **local vibrations in Pd wells**, not the connected nearest-neighbour D hopping band used by the full-loading envelope. Rare occupied neighbour pairs are not zero: an iid occupied-site illustration gives mean `12x` occupied nearest neighbours per D (not a measured spatial correlation or a connectivity proof). Such finite pairs/clusters do not supply the extended carrier and modulation assumed in that instability formula. The full-band 0.6423-meV/D threshold consequently has no calibrated meaning for isolated deuterons; it cannot simply be divided by loading. Host-mediated longer-range effective coupling would be a different, uncalibrated mechanism. Even the deeper `x≈0.6` region needs its own disordered/nonlinear calibration; this addendum supplies the beam-to-depth energy link, not that missing Hamiltonian.

Reproduce only this addendum with `python g13_rogue/depth.py` (NumPy/SciPy). It writes the two depth CSVs and `g13_rogue/output/depth.json`, including exact consumed-source hashes. Its actual invocation checked piecewise dilution and column anchors, independent angular/disk quadrature, direct/reflected ordering, hemisphere energy balance, plane-flux normalization and charge/nucleus units once. No old stochastic, stopping, seed or contact producer was rerun for this depth result.

## Fixed-basis stochastic addendum (item 5)

This is an **additive refinement of the same loaded-surface counterfactual**, not a simulation of the dilute track or deep receivers above. Original `g13_rogue/output/` files remain byte-identical. The new ensembles and full threshold readouts are separate in `g13_rogue/output/statistics_addendum/`.

**Reproducible realization definition.** `g13_rogue/statistics_addendum.py` draws complex Gaussian variates in fixed original **site order**, using per-phase Philox raw streams and an explicitly ordered Box–Muller transform. It applies the principal Bose covariance root constructed in a fixed 32-wave FCC Fourier gauge (exact fourth roots of unity). Neither the initial draw nor thermal innovations are assigned to a solver's arbitrarily signed/rotated eigenvectors. Linear/nonlinear, ambient/packet, and radial strata share each phase's thermal stream; packet phases have independent stable case keys. `g13_rogue/output/statistics_addendum/basis.json` records covariance/root defects below 5e-16. A deliberately rotated degenerate eigenspace changes the old mode-assigned draw by **1.818**, versus **1.49e-15** for the site-root draw. This removes the identified realization ambiguity; library transcendental functions and floating-point evolution may still differ. **Cross-host bitwise trajectory equality has not been observed or claimed.**

**Paired closing excess.** All original source scenarios use **64 independent phases**; only the prespecified **2.5-keV, eta=0.1** case extends to **256**, with explicit phase and ensemble identifiers. At the retained **0.0932881-eV** inward-closing threshold and chart current:

| Quantity | Added 256-phase result |
|---|---:|
| Nonlinear minus linear Palm maximum excess | +75.60 per D/s |
| Paired phase-cluster standard error | 139.56 per D/s |
| Simultaneous category interval, coverage at least 95% | [-901.35, 1047.94] per D/s |
| Nonzero paired phases in the three radial strata | 41 / 32 / 31 of 256 |

The sign is unresolved. These are **arrival-centered packet-minus-paired-ambient maximum** differences, not a stationary arbitrary-window CDF or reaction rate. The finite-category interval bounds the signed {-2,-1,0,1,2} phase responses with simultaneous Clopper–Pearson intervals across strata/categories. It is **pointwise in energy threshold**, not a simultaneous band over the whole tail. Half the error budget is reserved at each of two focused looks; the 128-phase interval was [-1579.67,1667.96], and the preset absolute half-width of 1000 per D/s was reached at 256 (achieved 974.64), without choosing a favorable maximum or sign. For scale only, ambient maximum-survival divided by `32 D × 32.8 ps` is **3.76e8/8.71e8 per D/s** (linear/nonlinear); that window normalization is not an independent-event rate. Zero empirical errors elsewhere retain nonzero category bounds and one-phase resolution. No near-zero relative-precision claim is made.

**Does the roughly 1.6 tallest-closing ratio survive?** In the matched near-source 256-phase sample, the linear/nonlinear maxima are **0.14567/0.24323 eV**, ratio **1.670**. Thus the extreme ratio persists in this sample, but it is not a universal amplification factor: paired geometric-mean ratio is **1.310**, approximate phase interval **[1.262,1.359]**, and median maxima are **0.08999/0.11736 eV**. `g13_rogue/output/statistics_addendum/phase_maxima.csv`, `g13_rogue/output/statistics_addendum/comparison.json`, `g13_rogue/output/statistics_addendum/survivals.csv`, `g13_rogue/output/statistics_addendum/rates.csv` and `g13_rogue/output/statistics_addendum/current_maxima.csv` retain all sources, strata, thresholds, individual phase maxima and their distinct normalizations.

**Motion versus extra quartic score.** Closing is always **projected inward relative kinetic energy** `m max(c,0)²/4`; no local `chi*n²/2` energy is added to it. The solved nonlinear field has a frequency-dependent coordinate derivative. Decompose its signed closing speed as `c=c0+c_chi`, where `c0` uses the common linear `L` and damping clock on the **same** field, and `c_chi` is the explicit hardening-frequency term. At each nonlinear window's actual peak, the mean same-time kinetic decomposition is

`<K> = <m c0²/4> + <m c0 c_chi/2> + <m c_chi²/4>`

`= 0.088581 + 0.029387 + 0.002551 = 0.120519 eV`.

These terms are evaluated at the same peak, not separately maximized. With the hardening clock removed only from the diagnostic derivative, the paired maximum difference is **−2.347 ±1.340 meV** (phase SE), interval **[−5.368,0.673] meV**. Mean common harmonic local maxima are **0.10179 eV linear / 0.09084 eV nonlinear**, while the nonlinear branch's own local score is **0.10316 eV**. The larger derivative-based closing peak therefore comes predominantly from the **frequency-to-velocity term**, not extra quartic energy added to the closing score or demonstrated spatial concentration. This is a diagnostic decomposition, not a new evolution with `chi=0`. Some sampled action still exceeds the two-transition law's `n≤1` range (about **1.66% linear / 1.36% nonlinear** in the focused ensemble); the full continued tails remain visible, and coordinate support alone does not extend that calibration.

Reproduce the new statistics only (NumPy/SciPy/**Numba**) with:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 python g13_rogue/statistics_addendum.py --workers 4
```

Use `--readout-only` instead of `--workers 4` to consume the retained arrays without sampling. The producer and finalized array consumer completed in the statistics slice; the integrated replacement's actual `--readout-only` consumer also completed. Source/producer/consumer hashes are in `g13_rogue/output/statistics_addendum/inputs.json`; the readout revision is distinguished from the production revision. No original deterministic producer, stopping split, seed or contact replay was repeated.


## What this model keeps and drops (recorded before code)

This adds a driven-wave axis, not a new material or a changed deuterium profile. It keeps a single polarized deuterium vibration, the repository's local well and discrete excitation spectrum, signed neighbour forces on the occupied graph, amplitude loss with the focusing convention, transported impact energy and the beam-current readout. It compares nonlinear and linear evolution with identical dispersion, packet energy, random phases and thermal preparation. Travelling waves arriving from other impact positions are included explicitly through their causal footprint and attenuation, rather than counting only impacts on the observed pair. A separately prepared modulated carrier illustrates instability and its resulting finite-time hot spot; it is not identified with spontaneous beam preparation.

It drops new electronic-force calculations, changes to the foil loading/depth profile, later atom-level passes, pure electron-plasma ringing, and a fitted conversion from deposited energy to optical-mode action. That conversion will be reported as an explicit source-derived proxy and sensitivity, not a measured material constant. Thermal excitation excludes zero-point energy as a usable source. Screening in a collision does not supply the missing frequency derivative with electron density. Contact will consume an actual signed separation/velocity episode, with quantum preparation and spectral support retained, not a scalar peak assigned a reaction rate.

The deliverable contains only `g13_rogue/`. Source files elsewhere remain unchanged. Full equations, consumed-input ledger, output paths, finite sampling support, replay results and the three largest omissions are completed below with the capsule.

## Original loaded-surface scenario: computed sizes

The [110] restricted well gives first/second transition energies **42.6473/47.9886 meV**, hence hardening **5.34132 meV per quantum**. At the engine-equivalent 8.2-ps amplitude lifetime, the lowest sampled full-loading threshold is **0.01505 quanta/D**, **0.6423 meV/D**, peak displacement **0.03827 Å**, at modulation wavelength **11.29 nm**. This is a coherent-carrier threshold, not a threshold on incoherent thermal energy.

| Loading/source preparation | Modulation wavelength (nm) | Threshold (meV/D) |
|---|---:|---:|
| Full loading | 11.29 | 0.6423 |
| 0.875 spread, clean / near probe | 12.76 / 12.76 | 0.6857 / 0.6852 |
| 0.75 bunched, clean / near probe | 8.319 / 6.925 | 0.8428 / 0.8257 |
| 0.875 bunched, near probe | 8.319 | 0.7284 |

The partial closures already shift some sites **0.159–0.262 Å**, beyond their **0.03-Å** measured small-probe support. Those threshold numbers are transferred-law contrasts, not calibrated partial-loading nonlinear predictions. They do not establish the observed loading dependence of a foil reaction signal.

One explicit pumped carrier starts at **0.12326 quanta/D (5.2972 meV/D)** with a **12.321-nm** modulation. Its finite evolution reaches **69.7466 meV/D**, equivalent longitudinal width **2.649 nm**, and **970 fs** contiguous half-peak duration. This is not a freely produced three-dimensional beam hot spot. Its peak action **1.4954** crosses the two-transition law's declared low-action range; the first sampled crossing is **7190 fs**, with 10-fs reporting resolution. That law's `n=1` local-energy boundary is **45.3180 meV/D**. The replayed 7601-fs closing episode therefore uses that explicitly marked Kerr continuation even though its actual displacements remain inside the 0.8-Å radial extent. The simultaneous inward kinetic peak is **0.0859065 eV**, at **R=2.88941 Å**; the separate minimum is **2.34631 Å**. They are not one state.

The normalized 300-K bath carries **10.1847 meV/D per retained polarization**, or **0.24383 thermal quanta/D**, excluding vacuum energy. For the `eta=0.1` source scenario at the chart-derived current, the stationary source-balance sizes are:

| CM energy (keV) | True current (µA) | Direct optical level (eV/D) | Acoustically converted D level (eV/D) |
|---|---:|---:|---:|
| 0.25 | 38 | 4.873e-17 | 6.397e-15 |
| 2.5 | 138 | 1.157e-15 | 2.944e-14 |
| 6.5, scaled-partition proxy | 140 | 3.096e-15 | 7.879e-14 |

Those are conditional source-power/loss balances, not a measured coherent occupation. They do not give a probability that the thermal population aligns into the pumped carrier. The input-driven prediction is not energy-flat: changing recoil partitions and track depth changes the supplied level. Above 2.5 keV its extension is the named partition/geometry proxy. The observed flat added signal and rising loading are therefore comparison targets, not parameters fitted here.

Across the stored finite phase/bath ensemble, the largest local energies are **0.14034/0.15143 eV** (linear/nonlinear), and the largest inward pair energies **0.11467/0.18358 eV**. These aggregate extrema combine preparations and do not by themselves measure an enhanced beam tail. All sampled amplitudes pass the 0.8-Å coordinate diagnostic, while **1.54%/1.57%** of site-times cross the two-transition low-action range. The whole continued distributions remain in the tables. `g13_rogue/output/current_maxima.csv` adds current-weighted arrival-centered maximum excesses and **12-phase cluster standard errors**, including the paired nonlinear-minus-linear comparison. For example, at 2.5 keV, `eta=0.1` and the 0.09329-eV closing threshold, the nonlinear single-packet maximum excess is **−43.6 ±29.4 per D/s** (empirical phase SE). That signed difference is not a negative physical event rate. This retained phase uncertainty, rather than the largest single sampled value, controls any comparison of those maximum tails; unseen bins remain unresolved.

The actual radial replay gives **log10 proton chance per normalized pair/pulse = −114.248…−110.330**, and either ordinary DD branch adds 0.29364 dex. The separately reported **peak proton hazard** is **10^(−101.032)…10^(−96.822) per pair/s**; it is not an integrated probability or a beam formation rate. Changing only the retained energy headroom changes exposure by **1.565–1.680 dex**. Quantum/classical RMS radius mismatch is **0.114–0.172 Å**; only **47.4–50.8%** of the exposure lies on the closing half of the episode. Thus these are finite-band, explicitly prepared conditional chance sizes, not a converged material-pulse probability. The unchanged 0.0625-fs integrator was exercised; no new time-step convergence claim is made. Full occupancies, signed excesses and branch probabilities are in `g13_rogue/output/contact_sizes.csv` and `g13_rogue/output/contact.json`.

## Equation and meaning of its coefficients

One complex amplitude `a_i` carries action `n_i=|a_i|^2` in quanta for the [110] polarized D vibration. Its deterministic equation, in eV, Å and fs, is

`i hbar da_i/dt = sum_j L_ij a_j + chi |a_i|^2 a_i - i hbar a_i/tau_A + F_i(t)`,

with `L_ii=Delta01`, `L_ij=-(hbar^2/2m)/Delta01 * J_ij` off diagonal. The corresponding coherent-state energy is `a* L a + chi sum_i n_i^2/2`. The printed local energy is the positive on-site diagnostic `Delta01 n_i + chi n_i^2/2`, not a unique partition of the interacting Hamiltonian. The linear comparator sets `chi=0` but retains the same `L`, forcing spectrum, packet harmonic-energy normalization, phases and thermal-bath innovations. Its local energy uses its own linear Hamiltonian. A different energy tail can therefore reflect both changed motion and anharmonic energy at the same action; the stored action support and full maximum samples separate those readings.

This is a rotating-wave envelope with the declared harmonic transition-matrix-element closure `|u01|^2=(hbar^2/2m)/Delta01`; it is not an exact projection of every anharmonic coordinate matrix element. The same closure fixes both hopping and action-to-displacement normalization, rather than using one normalization for instability and another for the pair history.

`Delta01=E1-E0` and `chi=(E2-E1)-(E1-E0)` come from a one-dimensional nuclear eigenproblem on the accepted [110] ray. They are discrete transitions, **not** the G10 SSCHA free-energy Hessian. This scalar restriction is not G10's full three-dimensional excitation spectrum. The coherent-state quartic is the normal-ordered Kerr model: its level energies are `m Delta01 + chi m(m-1)/2`, so its first two transitions match by construction and its classical action-frequency shift is `chi n`. Higher transitions, coordinate leakage outside 0.8 Å, and the independent classical action-frequency integral are all retained in `g13_rogue/output/transitions.csv` and `g13_rogue/output/action_frequency.csv`. Action above one quantum is a visible continuation of that two-transition law, not newly calibrated large-amplitude spectroscopy. No vacuum half quantum is counted as deposited or thermally usable energy.

The signed bulk dispersion is `epsilon(k)=Delta01-(hbar^2/2m)/Delta01 sum_delta J_delta cos(k.delta)`. Each original folded bond is split among its equally short periodic images; its weight is not multiplied by infinite-lattice shell counts. For a carrier `k0` and modulation `p`, put `d=[epsilon(k0+p)+epsilon(k0-p)]/2-epsilon(k0)`. The conservative sideband growth is `sqrt[-d(d+2 chi n)]/hbar` when the square-root argument is positive, and zero otherwise. Subtract `1/tau_A` for net amplitude growth. The odd sideband difference gives propagation, not growth. For hardening `chi>0`, only a negative-curvature sideband with `d<0` can focus; a band minimum does not acquire instability merely by being a band edge. At fixed negative `d`, growth exceeds loss when `n>[d^2+(hbar/tau_A)^2]/[-2 chi d]`.

`g13_rogue/output/growth.csv` gives the full wavelength/action/loading/loss grid, including stable Gamma rows, threshold energy in meV/D, amplitude in Å, and source/low-action support flags. Partial-loading rows preserve the actual occupied graph and all measured local forces in an unrelaxed cubic equilibrium closure, then take its explicitly labelled translational average for this dispersion. That average is not a disorder eigenmode. Only the measured source-axis curvature is changed; unprobed/transverse curvatures and quartics retain G10's stated transfers. A threshold outside the partial source's 0.03-Å displacement range is marked there rather than treated as measured partial-loading nonlinear physics.

### Loss convention repaired locally

The focusing README calls 4100 fs an amplitude-loss time, but its actual engine multiplies velocity by `exp(-dt/4100 fs)`. That is viscous drag `-v/4100 fs`: in the underdamped harmonic limit the **energy** lifetime is 4.1 ps and the **amplitude** lifetime is 8.2 ps. The main envelope uses `tau_A=8200 fs`; the growth table also retains 4.1 ps and 2.05 ps as stronger-loss sensitivities. Stationary injection–loss balance uses `tau_E=tau_A/2`, never `tau_A` twice. Acoustic attenuation is separately assumed, not measured by either the optical lifetime or that label.

## Input ledger and source interpretation

Exact byte hashes and sizes of consumed files are in `g13_rogue/output/inputs.json`. No new force data or foil profile is produced.

| Input | Repository owner and use |
|---|---|
| Local Kii, signed shell forces, lattice positions | `g3_loading/x1/dense_Pd/k3/final-x1-dense-readout-20261114T0213Z.json`; `g10_quantum_site/model.py` supplies its existing occupied graph. Local Kii is used once; `Kii-sum(J)` is diagnostic, never subtracted again. |
| Partial occupied indices, original probe axes, K and all F0 | The five actual partial source paths enumerated by `g10_quantum_site/model.py`, byte-bound in `g13_rogue/output/inputs.json`. Source-shell means are multiplied by that occupied graph once. Missing source-shell means retain G10's explicit full-loading closure, not invented zeros. |
| Wide well energies and outward slopes | `g10_quantum_site/bulk-refit-080/results/potential.json`, with source bytes in `g10_quantum_site/bulk-wide-wells/bulk-wells.json` and near-origin K/K4 in the local readout. The scalar callback retains the accepted Hermite interpolant, the original even/signed-axis convention, and its declared quadratic exterior. Direct E/slope anchor assertions run when it is built. This avoids the upstream module's missing hidden overlay without changing that module. |
| Beam current, charge and composition | `g12_residual/inputs_ledger.md` and `g4_beam/forward/model.py:I_cup,collected_fraction`. The 85/15 D+/D2+ current mixture supplies 1.15 nuclei per elementary charge, not 1.15 electrical currents. |
| Deposition and track depth | Fixed PdD records in `g11_cmns/collective/part1/stopping.json`, with the SI ODE definition in `g11_cmns/collective/part1/motion.py`. Electronic, host-recoil, direct-D recoil and unresolved cutoff energy remain distinct. |
| Optical loss | The actual velocity update in `g11_cmns/collective/focusing/engine.py`, interpreted as above. |
| Electron term | `g6_response/README.md` and `g6_response/inputs-ledger.json`; the absent phonon frequency-density derivative remains null. |
| Conditional reaction consumer | `g11_cmns/collective/contact/bulk_pair.json`, its fresh matched spectrum in `g13_rogue/input/contact_spectrum.json`, and the unchanged `g11_cmns/collective/focusing/pulse/replay.py`. |

**Track coordinate correction:** Part 1's field named `terminal_depth_A` stores the unconverted SI depth returned by its ODE; here it is read as metres and multiplied by `1e9` for nm. Its `timeline.depth_A` is already correctly in Å and is divided by 10 only for the agreement assertion. The central D+ depths are about 23.50 and 128.94 nm at 0.25 and 2.5 keV CM. These are used as an explicit uniform-depth deposition proxy, not an inferred range measurement. The local cell has side 0.8168 nm, footprint `(0.8168 nm)^2`, and exactly 32 occupied D: `density * area * depth =32`. A packet's local share includes the cell-depth/track-depth factor once. Assertions check the converted terminal/timeline agreement and cell normalization. Earlier Part 1 results are not assessed or modified.

## Continuous beam and waves arriving from elsewhere

`g13_rogue/output/drive.csv` keeps cup current, assumed cup/true-current ratio, true current, nuclei flux and each energy partition separately. The non-fitted focusing chart uses `b_hi=1,b_lo=0`: above 2.75 keV CM the inferred true current is 140 µA even when the cup current is larger. The nominal 150-µA scale is shown separately in `g13_rogue/output/causal_overlap.csv`; it is not substituted for the current-versus-energy input. At fixed true charge, a single-packet contribution is independent of current; an ambient time rate per charge scales as `1/I`, and a two-impact term scales as `I` per charge. Peak flux is not an integrated reaction probability.

The saved energy-resolved tracks exist at 0.25 and 2.5 keV. The 6.5-keV display is explicitly a proxy: retain the 2.5-keV partition ratios and depth, scale resolved energy to `2 E_CM -115 J/C`, and keep the cutoff unresolved. It is not a recomputed 6.5-keV stopping curve. Mode conversion is uncalibrated. The two linked scenarios `eta=0.001,0.1` assign `eta` of direct-D recoil to optical packets, `eta` of host recoil to travelling acoustic waves, and `eta` of that acoustic energy to the local D receiving channel. They are sensitivities, not best material values. Electronic stopping and the unresolved cutoff are never silently reassigned to coherent vibrations.

Remote impacts are an incident two-dimensional Poisson process. Acoustic transport uses a 4-km/s expanding annular sheet, uniform through the declared track depth, with energy attenuation `exp(-r/L_E)`. At first arrival `t=r/v`, its energy is divided by `max(cell area,2 pi r cell side)`, so a cell receives a fraction of the travelling shell energy rather than an entire impact. Conversion into the D envelope and subsequent graph propagation use the same equation as the local packet. Three radial **area** strata retain near-source, intermediate and remote arrivals, with their radii, annular areas, delays and weights in `g13_rogue/output/packet_geometry.csv`. The receiver clock starts at that first arrival; different delays preserve the homogeneous Poisson arrival rate. Midpoint area quadrature is finite geometry resolution, not an exact Green function or a calibrated foil boundary condition. The sampled radius ends at five energy attenuation lengths; the omitted annular energy fraction is `exp(-5)=0.006738`, while the nonlinear tail response to it is unresolved.

For a hard causal interior that counts an impact only after its wave can reach the point, and until source age `L/v`, the mean count is `pi F L^3/(3v)`. Counting a disk for a full `L/v` instead is a different residence model, three times as large. Survival-weighting the causal interior gives `2 pi F v^2 tau^3`. These counts measure active arrivals, not energy at a point and not a coherent closing projection. `g13_rogue/output/transport.csv` therefore also gives the separate energy-balance level `F E_packet tau_E/(D density * track depth)`, with named acoustic lifetime sensitivities, plus conversion into D. A shell/GF energy density and an interior residence count are not interchangeable.

## Full tails, thermal reference and rate semantics

`g13_rogue/output/survivals.csv` retains every energy threshold, not just a percentile: local and pair-closing occupancy, per-window largest-local and largest-closing survivals, and upward threshold crossings. `g13_rogue/output/maxima.npz` retains every phase realization's raw maxima. Zero sampled exceedances mean **unresolved**, with one-count maximum-survival resolution `1/12`; no zero is promoted to a physical tail cutoff. Site-times and neighbouring bonds are correlated and are not independent sample counts.

The 300-K ambient preparation is a harmonic Bose complex-amplitude ensemble with `n_B=1/[exp(epsilon/kT)-1]`. Its discrete linear OU bath and the nonlinear comparator use identical noise draws and damping. The nonlinear bath is a specified stochastic classical-envelope model, not an exact nonlinear quantum canonical distribution. The linear baseline is normalized by its modal Bose covariance and each incoming packet by its quadratic energy in the same `L`. The source has random global phase; neither that phase nor total packet energy is substituted for the signed nearest-neighbour closing projection. The nonlinear Hamiltonian includes additional quartic source work beyond the shared harmonic packet normalization; this is not a measurement of conversion efficiency.

The source normalization is harmonic energy, not a hidden equality of nonlinear Hamiltonian work: `g13_rogue/output/packet_geometry.csv` also retains the mean/min/max **actual injection work** in each comparator. Cross work with the thermal field can have either sign for an individual prescribed phase. The undepleted recoil reservoirs and uncalibrated conversion are not identified by that phase-dependent work.

`g13_rogue/output/rates.csv` separates the ambient background, the **single-arrival Palm excess** (the actual conditional packet-minus-paired-ambient history integrated over time and weighted by incident area/flux), and the finite-memory probability of two or more arrivals. Excesses may be negative; they are not probabilities. No beam is artificially accelerated. The `0.1,1,10` current multipliers change the physical Palm weights, not simulated arrival spacing. The unknown overlap-tail coefficient stays null: the Poisson probability is not an upper bound on nonlinear amplification or a substitute for an overlap initiation rate.

`g13_rogue/output/current_maxima.csv` is deliberately a different rate observable: one maximum per **arrival-centered** receiver window, with ambient maximum survival and a signed Palm excess per D/s. It is not the CDF of an arbitrary stationary laboratory window, and it is not the multiple-crossing count in `g13_rogue/output/rates.csv`. Zero cluster SE with no sampled response is a finite-resolution row, never a physical zero. `g13_rogue/readout.py` derives this readout and the contact sizes from the preserved arrays without repeating a simulation.

The window lasts four amplitude lifetimes. An isolated linear packet then has residual energy fraction `exp(-8)=0.00033546`; a nonlinear post-window tail is not assigned zero. Initiation rates are **2-fs-grid upward crossings per D per second**, with each unordered nearest-neighbour bond counted once and divided by the D count, not a fitted independent-event frequency or a nuclear rate. Occupancy is a time/site or time/bond fraction. The local energy Hamiltonian and signed pair relative kinetic energy are different observables. `g13_rogue/output/statistics.json` reports coordinate and low-action support fractions; distributions include the visibly labelled Kerr continuation rather than suppressing their large samples.

## Prepared hot spot and actual contact replay

The bounded seed calculation uses the bulk [110] dispersion on a 256-layer slice with transverse uniformity. A weakly modulated carrier is held by an explicit coherent pump; this supplies the uniform background needed by the linear stability calculation. It is a preparation-specific longitudinal concentration, not a prediction that the sparse beam prepares that phase or a three-dimensional localized volume. `g13_rogue/output/seed.csv`, `g13_rogue/output/seed_field.npz` and `g13_rogue/output/seed.json` retain the entire evolution, energy/participation/longitudinal width and sampled half-peak lifetime, including censored holds and source support.

The physical displacement mapping is `u_i=2 sqrt[(hbar^2/2m)/Delta01] Re(a_i)` in Å. Its velocity is the derivative of the **same solved field**, including carrier phase, damping and prescribed forcing. For a nearest-neighbour [110] pair two layer intervals apart, `R=R0+u_j-u_i`, `c=-dR/dt`, and `K_rad=m c^2/4`; receding samples have zero **inward** energy, not an invented inward sign. The actual deterministic history around its strongest complete closing/rebound episode is stored at 0.5-fs resolution in `g13_rogue/output/histories.npz`, with the pair-hash-bound manifest `g13_rogue/output/histories.json`.

The unchanged replay consumes those actual R/c samples, finds signed outer/inner/outer turns, retains both radial joins and both flagged-row treatments, regenerates and checks the exact spectrum's energies/contact gauges, and propagates its declared normalized spatial `l=0` Gaussian preparation in both 1.5/4.5-quantum headroom bands. It preserves positive exposure, once-subtracted signed excess, contact occupancy and proton/either-branch probability separately. Output is `g13_rogue/output/contact.json`. Band dependence, the quantum/classical trajectory mismatch and any crossed envelope support are part of that conditional result; neither a scalar peak, a peak-times-duration estimate nor a physical active-pair population is substituted.

The consumer also prints its historical one-pulse-per-incident-D conversion. That is **not** an inferred formation rate for this pumped seed, and is not combined with the beam peak flux. A physical foil reaction rate and the collision-screening route remain unknown. An electron-density term would require `g_n=d(hbar omega_D)/d ln(n_e)` in eV at specified polarization/loading/strain, together with the density ripple per acoustic strain. The dielectric/drag tables do not provide it. Its coefficient is null, never zero or fitted back from collision screening. Pure electron plasma ringing is noted but not built.

## Reproduction and reading the figure

From the accepted base plus this directory, with NumPy, SciPy and Matplotlib:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python g13_rogue/run.py --spectrum g13_rogue/input/contact_spectrum.json
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python g13_rogue/readout.py
python g13_rogue/package.py
```

The first command produces every scientific table/figure above and invokes the existing pulse CLI. The packaged spectrum is a fresh solve bound to the **current** pair descriptor, not the older historical contact result with another input hash. To regenerate just that input:

```sh
python g11_cmns/collective/contact/contact.py --pair g11_cmns/collective/contact/bulk_pair.json --output g13_rogue/input/contact_spectrum.json
```

To consume only the saved new history:

```sh
python g11_cmns/collective/focusing/pulse/replay.py --pair g11_cmns/collective/contact/bulk_pair.json --spectrum g13_rogue/input/contact_spectrum.json --histories g13_rogue/output/histories.json --trace-dir g13_rogue --output g13_rogue/output/contact_replay.json
```

`g13_rogue/output/waves.png` shows signed loading-dependent dispersion, the prepared seed's full maximum-energy history, phase-resolved maximum survivals against the matched linear baseline, and the signed Palm pair-initiation excess versus threshold. The figure's selected current/energy is not the whole tail: the CSV files retain all simulated thresholds, currents, energies and source-conversion scenarios. `g13_rogue/output/result.json` is the compact numerical readout.

### Observed invocation and package

The successful settled invocation was `python g13_rogue/run.py --spectrum g13_rogue/input/contact_spectrum.json --reuse-solved`, with the thread environment shown above. It preserved the already computed transitions and deterministic seed after the local coordinate interpretation repair, produced the drive/transport/full-tail/current tables and figure, and naturally completed the **actual unchanged pulse CLI with exit 0**. The final `python g13_rogue/readout.py` consumed those saved results and produced the current-maximum uncertainty and reaction-size tables. No broad test suite, upstream force run or later pass was run.

`g13_rogue/package.py` builds `g13_rogue/capsule.tar.gz` as an exactly scoped replacement-directory overlay on accepted commit `09b1b0779ed83f014d3e9ab21328fcd106f2a1bb`. It supersedes original draft asset `793110859` (archive SHA-256 `c82650415c8aa2188251cc918c89aec0960cb3b805e0eb9eebc78b95404a2ea8`), not its scientific evidence. Apply only listed `g13_rogue/` members; if that directory has subsequently changed, compare those members before replacing it and preserve unrelated additions. `g13_rogue/manifest.json` and `g13_rogue/SHA256SUMS` bind its source/input/output members; archive hash and authenticated unpublished-release/asset identities remain external delivery facts, avoiding a self-hash cycle. Nothing outside this directory is packaged. The retained result does not depend on publication, merge or a later pass.

## Biggest three omissions

1. A measured, geometry-resolved conversion from each recoil/electronic channel into optical and acoustic D-driving modes, including actual acoustic attenuation and foil boundaries. This dominates the physical beam-preparation uncertainty.
2. A calibrated nonlinear/disordered vector wave Hamiltonian over the sampled action and angular range. This pass uses a scalar two-transition law, transferred partial-loading closures and a finite stochastic bath; its continued tails are not new force data.
3. The material's quantum pair preparation/population and phonon frequency-density derivative. Conditional contact replay does not identify those missing quantities or a physical reaction rate.

Later atom-level checks, matched electronic response and other materials remain listed only; none is started by this capsule.
