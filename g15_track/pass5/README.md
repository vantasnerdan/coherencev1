# G15 Pass 5 — signed whole-foil conditional receiver

This additive draft consumes exact selected Pass4 members on base `f6010129f11ba8010815b78a67b3b6bed1db7194`. It does not rerun or modify accepted producers, start DFT, restart stopped campaigns, commit, push, merge or publish. `selected_inputs.json` binds the actual inherited files; the two accepted manifests are preserved with their exact hashes. The draft contains conditional model results, not a calibrated physical tail exclusion or a nuclear contact probability.

## Governing observables

The mechanical projection is **volumetric strain** `e = tr(epsilon)` and tensile-positive mean normal stress `sigma_h = K_bulk e`. The returned strain sampler is the declared isotropic projection `epsilon = e I/3`, not a reconstructed shear tensor. A negative trough, a complex trace-strain intensity hole, and a physical total-energy hole are different observables; the complex projection diagnostic is not relabelled as total acoustic energy.

For the saved eight complex trace-band amplitudes at each actual incident-probe coordinate/time, `R = |sum b_band|² / sum |b_band|²` measures **cross-band projection cancellation**, not total-energy depletion. `probe_projection_holes.csv` gives the 0.1, 0.5 and 1 cutoffs under the normalized target-column-weighted probe exposure, with conditional cluster errors. This is neither a time-only fraction nor an event/contact probability; the raw coordinates, times, weights and intensities are retained in `probe_projection_samples.csv` and the same band arrays feed the signed covariance consumer without recomputation.


For one electrical parent, all represented bands and retained direct/face/edge copies are summed into one signed receiver kernel `h` before forming Campbell cumulants:

`kappa_n(x,t) = integral rate(s) E[h(x,t-s,mark)^n] ds dmark`.

The source uses current accepted complex mechanical mean, late LTE and G16 pressure amplitudes. Square-root recoil variance is not a coherent phasor. A catalog of actual centered-Poisson recoil draws reuses marks across frequency, direction and copies; two molecular constituents share the electrical parent and echo delay while their collision innovations remain separate. The finite 64-draw recoil approximation is priced against the exact diagonal source second moment. Signed source work is retained; receiver projection does not add a second carrier or a heating source.

The AC packet is the causal derivative of a finite Hann window with an integer number N≥2 of carrier cycles. For `z=wavelet exp(-i omega t)`, its integral and `integral(z²)` vanish, `integral(|z|²)=1`, and the carrier Fourier coefficient is positive real. Consequently `sqrt(2) Re(A z)` has isolated-port real-flux energy `|A|²` for either launch quadrature. This does not remove cross-port or reflected-path interference. Static thermoelastic strain/permanent deformation are separate. Both the initial unipolar packet and the later complex-L2-only gamma packet are superseded; the gamma age law is only a positive importance proposal with its explicit density. The exact finite support is used for Cox convolution and the union of retarded local birth windows before drawing physical parents.

All eight current source ports are represented, through the 0.25–2.5 THz inherited band. Old G13 positive integrated frequency energies provide a **shape sensitivity**, not native pressure calibration or unknown per-node native phases. Port energies and signed source phases remain the current Pass4 ones. The named bandwidth alternatives are 0.5, 1 and 2 times that retained shape width. Bulk phase/group dispersion, the retarded D susceptibility, and reflected piecewise depth-loading integrals remain distinct.

The misleading historical track field `terminal_depth_A` is already in metres; only timeline `depth_A` rows are in Å. The exact 2.5-keV terminal ranges are 128.944460425 and 81.871005481 nm. Both the source centroid and unresolved source/Fresnel support retain this distinction.

The force transform is `Ftilde = integral F exp(+i omega t - i k.r)` and reconstruction is `exp(-i omega t)`. Source coefficients stay **unconjugated**. The outgoing 3-D radiative strain projection supplies **+i**, and the positive-slope circular guide Green supplies `exp(+i pi/4)`; the guide receiver profile already contains its spatial derivative. The guided eigenvector gauge is fixed explicitly and its known normal-impulse phase uses the conjugated source mode. A unit-modulus source-basis map preserves accepted per-node launch energy and signed work. Residual pressure/recoil/LTE guide-force phase remains a declared scenario, not a calibrated modal projection. Only the four supported in-plane source directions are interpolated for guides. The omitted nonradiative `1/(kR)` correction is separate from finite-source phase error.

For a mass-normalized guide, `integral rho |u|² dz = 1`, so `h = integral |div(u)/omega|² dz` has units Pa⁻¹: there is **no extra division by foil thickness**. Fifteen-node Gauss–Legendre integration of the actual degree-14 modal polynomial on each constant-D layer gives `h` and `n_eff = integral n_D |div(u)/omega|² dz / h`. The same `h n_eff` is used once in the real phase, group delay and passive loss, including the live exterior variance; the receiver evaluates the same modal polynomial. The uniform-density reduction and zero-real/retained-loss control are recorded in `passive_controls.json`.




## Two receiver consumers

`sampling.py` integrates a positive whole-beam radial/age mark measure. Large near jumps stay finite Poisson; dense weak jumps use a Gaussian closure with an absolute-third-moment diagnostic. The complete parent kernel, not individual copy powers, supplies cumulants. The tabulated Bennett expressions apply only to the finite retained recoil/cubature model. They are not bounds on omitted source-interior physics. Random-time samples spanning two physical seconds give time fractions; they are not a turning-point census or a resolved acoustic upcrossing rate.

`live.Strain` supplies arbitrary actual coordinates and physical times to the optical packet consumer. It uses a whole-beam radial exterior variance and finite Gaussian spectral carrier realization, plus globally named uncompressed local Poisson parents. Actual coordinates enter carrier phase even below a nanometer; there is no propagation patch and no interpolation of an rms CSV into THz cycles. This **conditional, locally stationary exterior closure** does not reconstruct higher-return covariance or a measured phase chart. `matched_controls.csv` separately compares complete-parent, independent-band-parent and direct-only receiver laws over the same duration. `nested_receiver_maps.csv` contains the declared signed projection at nested coordinates from the whole-foil grid to 0.125 nm.

`low.py` retains the accepted 303 mass-normalized modes and exact source-surface projection rows, without the large predecessor basis or a producer replay. Prompt stopping and fast native-pressure/optical-stress impulses drive an actual Cox-convolved mean, including turn-on and bursts. Its same-parent Campbell covariance gives a correlation-independent innovation bound; no fictitious microscopic low-mode phase is inserted. FV/surface interpolation is explicitly declared, with the separate first 125-nm layer. `low_accepted_window.csv` and the selected modal states preserve the inherited full-force finite windows. Those windows do not bound the unclosed thermoelastic/slow-force continuation over seconds.

Pass4's `2.5054719912489504e-8` and warm `0.00556844834632475` are **fractional optical frequencies**, not strain. Every receiver comparator uses the evolved fixed-T Gamma pole at the same T/U, including `dQ/dK`; the full signed strain/stress is retained separately. Frozen-Q sensitivity is a named comparison, not the bare Gamma denominator. Actual packet transport and yield use the finite-q consumer rather than assuming the two native allocations have equal finite-q sensitivity.
The inherited `response.covariance` key `theta_per_linear_strain_eV` already contains `g/3`; numerically it multiplies the volumetric trace in this declared isotropic projection. It must not be divided by three again. The quiet operator input `linear_strain` remains trace/3.


## Depth and optical contracts

`medium.build(out,path4)` produces `medium.npz` and physical ledgers. `Medium(bundle).sample(ec,beam,xyz_m,mount='nominal',time_s=2)` returns `(T_K,Ue_eV_per_formula,gradT_K_per_m,gradUe_eV_per_formula_per_m)` for scalar or batched coordinates. `time_s=inf` is the independent mean steady state. Finite fields retain the named original Cox history at the stated heat time resolution. Mount alternatives are 1e3, 1e4 and 1e5 W/m²/K edge transfer. Tagged single-parent electron perturbations are separate, not added to the mean background a second time.

`medium.optical` evaluates the same native quartic/G16 landscape, including background strain `tr(epsilon)/3`, and returns angular frequency, group velocity and ray wavevector derivative. `medium_codec.json` describes a lossless uint32-delta/byte-shuffle/LZMA representation of the already calculated binary32 depth innovations. Every decoded history was compared bit-for-bit; compression does not recompute heat or discard finite-time arrays.

`meetings.build(out,path4,medium.sample,strain.sample)` uses distinct electrical parents, actual packet covariances, coordinates, times, flight depths, shared molecular echoes and matched controls. The geometry-event, live event-response and incident-probe Palm proposals have separate counts and denominators. A chronological pair is counted once; reflected residence is not a new parent. The event-conditioned anchor mean prices strain in the larger geometry quadrature; the live optical response/probe rows use pointwise signed strain. The smaller live design is a measured runtime choice, not a change in the geometric rate weights.

## Yield last, with replacement

`yield_consumer.py` preserves the original 8×8 Maxwell/BH quadrature and uses the full trajectory-relative speed `|v_probe-v_target|`. The flux denominator remains the original incident speed, not the shifted target-relative speed. Connected differences use log ratios and `expm1`. The accepted `yielding.Kernel.difference` is an independent covariance-control consumer.

The final per-incident-D addition is replacement heat plus signed ordinary acoustic covariance response plus the first Palm term plus the connected second Palm term:

`thermal + acoustic + integral [F(parent)-F(quiet)] dPalm1 + integral [F(joined)-F(parent1)-F(parent2)+F(quiet)] dPalm2`.

Thermal replacement, ordinary acoustic response and both optical Palm terms use the same independent incident-probe time/depth/Cox weights over **0–2 s**. The quiet sheet and own 6.5-keV reference stay unchanged. The total Monte Carlo error is estimated from their sum within each probe, with retained-parent temporal clustering; it is conditional on the fixed field chart and does not include source/phase model uncertainty. The former 2.83e-7 pilot used a two-second endpoint heat budget and is not the exposure mean.


`joined - additive` separately isolates the retained pair dynamics. Existing Pass4 native additions are replaced, not added again. The unchanged sheet is included quietly once; original stopping/loading and the raw own 6.5-keV reference remain fixed. Signed first-order acoustic covariance reads the same passive susceptibility/loss port; no second-order mean heating credit is inserted. Supported struck/flight channels remain the separate exact inherited `prior_inflight_yields.csv`. A numerical encounter or a Gaussian covariance is not inferred nuclear contact.

The two rate comparators are U=0 and U=0.04 **keV** (the named 40-eV shift), not a fitted new screening energy and not the separate historical 40-meV threshold convention. Neither is equated with `Ue_eV_per_formula`. Phase-only frozen/shuffled/strain-off connected controls subtract their own matched one-parent velocities at the same transported population and work; they do not remove every physical source.

## Computed native Gamma/X allocation

Axiom's additional selected input is exact commit `f3114eb77af3343d6203a00efbd59e588cdbeb7f`; its three SHA-bound follow-up JSON files are under `native_Gamma_X_inputs/`. For `F=-Kii u+sum J u_neighbor`, the measured increments obey `Delta diagonal=(Delta K_Gamma+Delta K_X)/2` and `Delta opposing J=(Delta K_X-Delta K_Gamma)/2`. The latter opposes the diagonal change at Gamma. It is **not** a 75:25 convex mixture of the old finite-q endpoints.

`native_Gamma_X_inputs/producer_native_gammax.py.txt` preserves the exact completed producer bytes (SHA `4c443411a6d913aceb1538395312f69c3554c53eba64b1a6110c303defab0d89`). The saved native meeting reports hash this immutable evidence directly. The maintained module subsequently corrected only scalar-onsite indexing in the tagged snapshot observer; no launch, propagation or meeting result was replayed or re-labelled as a new producer version.


The declared representative is `Delta K(q,U)=Delta K_Gamma(U)+A(U)[1-F1(q)]`, with `A=3(Delta K_X-Delta K_Gamma)/4`, `F1(Gamma)=1` and `F1(X)=-1/3`. All constrained opposing change is assigned to shell 1 while aliased shells are held fixed. This reproduces the longitudinal Gamma/X data but is not a unique shell fit or measured material group velocity. Gamma uses the exact low-U crossover contrasts; below the first measured A knot its linear interpolation is a named closure. Mesh 12 and quadratic low-A alternatives are reported. The sampled 1000–1500 K sign change is not a mesh-resolved crossover bracket.

`native_Gamma_X` is a **frozen-incoming-source propagation/meeting scenario**: accepted onsite packet birth states and archived inner histories, original pressure and the same room G3 base/heat are retained. The common medium's finite-q stiffness, gradients, phase and transport are changed and signed joint work is priced. `coupling_frozen_launch` supplies the matched second endpoint on that same incoming source. `frozen_launch_controls.csv` uses within-probe paired errors; the original endpoint runs are not replayed or silently re-labelled. Zero-U sound susceptibility is identical, so its exact saved source charts and matching probe rows are reused, not injected as another carrier.

Only the 2.5-keV archives contain retained spatial electronic U snapshots. `native_Gamma_X_tagged_sensitivity.csv` applies the new law to those actual local values and the selected frozen covariances; its signed instantaneous potential contrasts are not a regenerated launch or a donor-work time integral. Full local `U(t,r)` is absent elsewhere. It is neither inferred from projected DeltaK nor replaced by Umax. Source-feedback and microscopic-contact conclusions therefore remain unresolved.

## Replay and delivery

Use the selected Python environment with BLAS/OpenMP/MKL threads fixed to one. After extracting **both selected archives** into an empty directory:

```sh
python g15_track/pass5/run.py --consume
```

This reads saved medium, packet and signed acoustic-probe arrays; it does not rebuild scientific producers. `cold_bundle_consume.json` records fresh extraction, exact selected-input checks, consumer code identity and byte-identical final outputs. `remote_delivery.json` records authenticated private-draft inventory, target base, download and hash verification. No public release is created.

The same entrypoint runs the bounded joined case (`--case`), the new frozen-launch native case (`--native-case`), and all case stages (`--all --stage receivers|meetings|probe|native|native-observers|assemble|controls|maps|report`). Energy stages accept disjoint selectors; `native` requires one. `native-observers` prices the retained local-U and refraction contrasts once before the final report. Assembly joins saved arrays and exact invariant source aliases without re-running an endpoint. No physical job is authorized beyond two hours.

## Readout files

- `receiver_statistics.csv`, `band_cumulants.csv`, six `E*_receiver.npz` histories: signed strain/stress and the separate exact Gamma comparator, over two physical seconds. Sparse PDF samples are not a realized crossing census.
- `source_projection.csv`, `passive_controls.json`: signed launch work, exact source variance versus finite recoil catalog, and isolated phase/loss controls.
- `finite_track_contrast.csv`: one resolved primary-stopping correction, added as cached + finite − far **before squaring**. It does not correct native pressure, recoil or LTE by analogy; its G12 values are kernel changes in m², not fractions of the residual.
- `medium_profiles.csv`, `medium_refraction.csv`, `nested_receiver_maps.csv`: actual depth/whole-foil fields, mount alternatives and nested x/y/depth receivers.
- `meetings.npz`, `meetings.json`, `meeting_summary.csv`: separate event geometry, live response and incident-probe Palm weights, coordinates, times, covariances, coherent velocities and carrier ledgers.
- `yield_outcomes.csv`, `decision_summary.csv`: all six chart energies (0.25, 0.5, 1, 1.5, 2 and 2.5 keV), three beams, native allocations, retained closures and rate comparators. Component and summed within-probe/cluster errors remain separate from physical model uncertainty. The exact inherited struck/flight table remains under `inherited/frozen/`.
- `native_Gamma_X_constraints.csv`, `native_Gamma_X_transport.csv`, `native_Gamma_X_refraction.csv`, `native_Gamma_X_maps.csv`: exact measured contrasts, named finite-q continuation and changed depth/whole-foil observer. `native_Gamma_X_receiver.npz` holds only the recomputed Gamma columns and aliases the unchanged signed source histories.
- `nativeGX_meetings.npz`, `nativeGX_meetings.json`, `nativeGX_acoustic_probe.npz`: third-allocation and matched frozen-launch-control states, with exact source-row aliases and no duplicate acoustic injection.
- `frozen_launch_controls.csv`, `frozen_launch_outcomes.csv`, `native_Gamma_X_tagged_sensitivity.csv`: paired consumer contrasts, control component uncertainties and the explicitly limited retained local-U sensitivity.
- `probe_projection_holes.csv`, `probe_projection_samples.csv`: separately defined complex cross-band trace-projection cancellations; not physical total-energy holes.
- `parameter_sensitivity.csv`, `omitted_scales.csv`, `parameters.json`: conditional bandwidth/temperature/quadrature sensitivity and unresolved-support scenarios. A sufficient far-field radius is not an observed interaction radius or a nuclear opportunity probability.

## Receiver evidence and interpretation

The completed endpoint grid has 432 two-second marginal laws, with 256 random-time samples per law. Its largest sampled absolute trace is `4.66933e-6`; the maximum absolute matched optical-frequency modulation is approximately `1.01e-5` (2.5 keV, fast Cu, coupling, hot steady, receiver 2). This is `402.235` times Pass4's wave comparator, and about **0.181% of the separate `0.00556845` warm comparator**, not a computed same-case heat ratio. These are finite sampled model extrema, not physical maxima or warm-tail exclusions.

The nominal hotspot cubature refinement from 12,800 to 102,400 marks changes the second cumulant by −34.6% and the fourth by a factor 8.33. Its rare tail is therefore **not converged**. The dense weak-parent standardized third-absolute-moment diagnostic is about 0.001–0.0015 in that contrast, but it does not repair the unresolved discrete near tail or native source interior. The tabulated time fractions and finite-catalog Bennett expressions remain conditional diagnostics; they are not acoustic crossing rates or universal probabilities. The live common-field packet consumer is a separately declared exterior/local closure, not silently equated with that marginal PDF.


## Three deciding omissions

1. Source-interior/Fresnel/nonradiative elasticity, residual guide-force and native pressure angular/within-band phases, and the rare recoil tail beyond finite sampled marks. `omitted_scales.csv` includes `kr=10` radiative support and source-size scenarios, geometric opportunity and a coarse energy cap. These are not measured interaction volumes or pointwise threshold exclusions.
2. The isotropized per-branch bulk reference, radiative Hankel guide approximation, first-order reference-temperature susceptibility, unclosed low thermoelastic continuation, higher-return/mode-conversion covariance and nonstationary exterior phase correlations. The finite first-return parent quadrature and locally stationary live Gaussian exterior are separate, not silently equated.
3. Close-wavevector optical exchange, microscopic phase-dependent covariance/contact, and backreaction beyond the leading Gaussian RWA/eikonal pair calculation. Opposite ray-work demand is priced as a limitation, not credited as unexplained heat.

These remain model uncertainties. An empty finite overlap sample is not a rate limit; a small conditional yield is not a proof that omitted physics is absent. The accepted Pass4 effective/rate-weighted opportunity volume and a realized tail are not the same quantity.

## Conditional incident-probe result

The final consumer contains 648 energy/beam/allocation/branch/rate rows: 162 main meeting cases plus 54 matched frozen-launch controls. The largest modeled addition is `7.63547e-7` of the G12 residual, with conditional Monte Carlo standard error `6.48902e-7` (2.5 keV, long retained, onsite, hot bursts, p branch, 40-eV comparator). It is dominated by replacement heat; errors include within-probe covariance and physical clustering, not the unclosed physical source/phase uncertainty.

Across the same grid, the largest absolute signed ordinary-acoustic contribution is `1.54997e-12` of the residual and the connected pair contribution is `4.09242e-18`. On the explicitly common frozen incoming source, the native-Gamma/X versus onsite propagation change is at most `3.71326e-20` of the residual. These conditional differences do **not** price regenerated source feedback: the retained local-U snapshot sensitivity changes a sampled optical frequency by as much as `0.00791868` relative to the onsite case. Full `U(t,r)` and corresponding source/donor-work feedback are missing, so the small fixed-source yield is not a physical exclusion or a contact-probability result.
