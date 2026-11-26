# Signed whole-foil sources — unpublished draft

Base: `fe514a174214b01558c228fc9f8982ef3791fb06`. Only this directory and its dependent `../rung2/` are new. The accepted `../rung1/` and `../../pass2/replication/` files remain authoritative and unchanged.

## Kept / dropped before propagation

Kept: accepted whole equal-area square foil, 250-µm thickness, full rotated cubic tensor and 8×8 seeded grains; frozen rung0 and rung1 operators and fixed physical observation regions; electrical-parent counts, accepted source profiles, random signed quiet 1% 5–250-kHz launch and causal beta 0/.35 feedback; both loss laws as separate alternatives. Added normal **inward** momentum from each electrical parent, and one-sided deposited heat diffusing through the thickness with Pd heat capacity/conductivity, coupled through the full eigenstrain-stress virtual work. Modal load signs are projection signs, not absolute values. All sources start together; static deformation is excluded from dynamic A/C and peak fields. Six beam cases use even/steady, hot-spots/steady, hot-spots/bursts, and the earlier **even** spatial ripple profile at 100 kHz, a beam-coupled nearby mode, and a literal between-modes frequency. The same frequencies and physical counts/exposure are used at both later spatial resolutions.

Dropped / conditional: no measured mount, damping, PdD thermal constants, absorbed-electron fraction or impact footprint is available. The primary dynamic thermal model neglects lateral diffusion over the short record (diffusion distance is reported against the source-cell and beam scales), but solves thickness diffusion spectrally, retaining insulated-face DC growth. A separate flange-isothermal lateral conduction solution supplies bounded static bias; it is not called a measured flange. Normal static deformation is the elastic shape about the center of mass with the rigid force/torque balanced by a distributed inertial reaction. Free-foil rigid translation is reported separately, never as bending. No atomistic impact resolution, new material fit, new random frequency band, or measured D-energy conversion is claimed.

The primary statistic uses 16 independent record IDs, paired by electrical parents and addressed marks. Record-level mean/SD and paired intervals, not bins or parents as replicates, decide whether a peak-to-own-mean change is resolved. Finite-observation response is not a settled resonance measurement. The producer and input ledger specify the exact thermal boundary, energy partition, damping integration and subtraction of forced bias.

## Equations and fixed readout

For mass-normalized mechanical modes, the common producer advances
`q_i'' + gamma_i q_i' + omega_i² q_i = Q_i[T]`
with pressure and random parent impulses added to velocity. It uses the exact viscous affine step for a bin-constant generalized thermal load. The two inherited `gamma_E` prescriptions remain alternative amplitude-envelope energy rates; this is not a claim that viscous instantaneous energy loss equals the old rotational model's constant fractional energy decay at every phase. Both randomized-only and combined controls use this same declared propagator, at all spatial rungs.

`Q_i = alpha (C11+2 C12) ∫ div(u_i) ΔT dV` is exactly the cubic isotropic-eigenstrain contraction `∫ epsilon_i : C : (alpha ΔT I) dV`; it is not a second external pressure. Cosine thickness modes solve diffusion with insulated front/back boundaries, including the undamped uniform temperature mode. The cumulative stopping tracks separately project electronic/host/direct-D heat, with cutoff energy conditionally deposited at its terminal depth. The 64-mode temperature field resolves the mechanical thickness moments and ~9-µm 100-kHz thermal skin, not the original nanometre collision. The actual finite-record lateral diffusion distance is recorded; lateral conduction of the **transient** is omitted relative to the 1.2-mm beam structure, not declared zero in actual Pd. The separate steady flange calculation uses lateral conduction and its explicit ideal isothermal edge.

The heat reservoir receives incident retained energy minus actual random and pressure work; thermal mechanical work is removed once after the affine step. Negative work returns energy. Modal damping leaves the tracked reservoir and is not credited again as fresh incident heat. This is a leading linear thermoelastic/work-coupled heat model, not a calibrated temperature-dependent coupled thermodynamic material law.

The stored dynamic field is displacement minus the instantaneous forced elastic bias, with the slow bias velocity separated. The thermal bias includes uniform-temperature growth; a steady pressure bias does not appear as vibrational strain. First-bin A/C are explicitly empty, and every later mask and feedback probability is evaluated before that bin's new marks. All signed component and origin-sector cross terms are retained.

Here “dynamic” means the field about the instantaneous elastic equilibrium, including finite thermoelastic phase lag and common-onset transients; it is not a claim that every sample is a freely oscillating mode. A time-independent mean spatial load alone gives a relaxed static deformation, not sustained growth of material vibrations.

The mechanical calculation is the force/torque-balanced elastic component in the CM frame. Distributed rigid inertial reaction cancels the unrepresented translations/rotations, doing no work on the mass-orthogonal elastic modes. The reported free-body acceleration is the alternative laboratory rigid motion, not added to the bend or to an elastic modal peak. The closed numerical heat/work ledger is for this stated balanced convention; a genuinely freely drifting body's tiny rigid kinetic energy is outside it, rather than silently being called elastic energy.

The primary height statistic remains the leading event over the same three 2.5-mm-square ×50-µm physical regions, divided by its selected region's own positive time mean. Fixed-region tails are also retained separately. Six beams × two beta × two losses × randomized-only/combined give 48 arms per spatial model, each with 16 independent records. Spatial maps/widths are the preselected expensive adjunct for records 0 and 1 in **every** arm; all 16 records have full sampled tails, durations/censoring, signed crests/holes and same-event 30/50/200-µm energies. Pointwise paired-ratio intervals are reported; a 95% interval wholly within [0.9,1.1] means only compatibility with less than a 10% change at this finite-record resolution. An interval crossing one is not called a nil effect.

Near-null fixed signed projections (operator norm below 1e-10 of the same observable's projector-family norm) have undefined normalized crest/hole statistics; `verification_r*.json` lists these explicitly. This does not remove any source arm, fixed-region energy history or same-event depth readout.

The conditional eV/D columns preserve the accepted positive sqrt(gamma)-weighted elastic-plus-kinetic quadratic, multiplied by 4.1 ps and divided by the frozen loaded-D population. Under the declared viscous dynamics this is a **cycle/envelope storage-transfer proxy**, not instantaneous viscous heating; the actual work/loss ledger uses the propagated mechanical energy change. Elastic eV/D is reported separately and neither quantity is measured conversion.

Pressure parents are aggregated into 50-ns impulses. `pressure_sampling_r*.csv` checks the **actual pre-kick** harmonic transfer and constant-mean fixed point, separately from the bin-constant thermal forcing. For rung1 at the three declared drives and both losses, the maximum complex transfer difference is 0.0428%; the stationary pressure-displacement sampling error is at most 0.0508% per retained mode, and the even-mean-pressure stroboscopic kinetic artifact is at most 1.76e-24 J. This bounded binning artifact is not interpreted as growing material vibration. The separate known-mode affine thermal-forcing check differs by about 1.55% in complex transfer, predominantly the half-bin sampling phase; this numerical phase is not a physical thermal diffusion delay.

`summary_r*.csv`, `contrasts_r*.csv` and `verification_r*.json` are consumer readouts. Per-arm `records`, `tails`, `peaks`, `signed`, `energy`, `same_event` tables and `histories.npz` retain actual independent-record evidence. `transfer_r*.csv` distinguishes settled harmonic gain/phase from the fixed finite exposure; a literal midpoint between modes is not labelled an antiresonance. `static_r*.csv` states the CM/inertial-reaction mechanical and ideal flange-thermal assumptions.

## Reproduction

From repository root with NumPy/SciPy and a C99 compiler, use `python3 g13_rogue/pass3/rung1b/signed.py run --rung R --beam B --beta X --loss L`, adding `--random-only` for the original-source control. R is 0 or 1 here (2 uses the dependent capsule), B is 0–5 in the stated beam order, X is 0 or .35, and L is `retained` or `Cu_proxy`. Defaults select all 16 records. Frozen `signed_r*.npz` consumers avoid numerical reassembly. `summarize --rung R`, `proof --rung R` and `package --rung R` are separate readout/check/packaging actions, not an implicit wider campaign. `prepare` deliberately creates a new numerical-input epoch and is not needed for replay.

## Representative signed clips

[Retained loss](clip_r1_hot_bursts_b035_retained.png) and [Cu proxy](clip_r1_hot_bursts_b035_Cu_proxy.png) are actual 200-frame APNG trajectories for preselected hot-spots/bursts, beta .35, solo record 0. Rows are u_z and epsilon_xx at 30 µm; columns are random, momentum, thermal and combined. **Each panel has an independent fixed sampled-clip physical scale**: equal colors do not mean equal amplitudes. Sampling is 0.5 µs over the 100-µs observation at 0.332-mm map spacing; playback is slowed to 20 fps. `clips_r1.json` authenticates all frames and verifies identical solo/ensemble source marks for both losses. Annotation-only reflow preserved every original physical-map pixel in all 200 frames without replaying propagation. The retained-loss clip's sampled displacement maxima are 349 nm random, 3.18 pm pressure and 0.511 nm thermal; Cu gives 1.071 µm, 6.25 pm and 0.917 nm. These are clip-sample maxima, not ensemble peak statistics. Viewing APNGs requires an animation-capable image viewer/browser, not the scientific runtime.

The complete P15/P15/P4 biased static solve includes the soft modes below 5 kHz and balances all six rigid modes. For the ideal flange-thermal/full-retention illustration, the source gives about 7.26 K (even) or 8.63 K (hot) maximum steady rise; center-minus-near-corner elastic deflection is about 23.16 nm or 16.30 nm, overwhelmingly thermal rather than pressure. These are **conditional static displacements**, not part of the dynamic A/C pre-field. Free-body incident acceleration would be 1.559e-3 m/s², not a static bend.

At the 1.525719-mm screen, a homogeneous full-cubic reference gives lowest free-face branches about 491 kHz along a crystal axis and 588 kHz along the diagonal. P4 versus P8 thickness shifts are below 3.3e-7 for the first three branches in the two inspected directions. This supports the thickness representation in that reference, **not** a continuum-convergence claim for the synthetic grained foil. The 1.526-mm **projection support** is not relabelled as an actually populated 1.526-mm propagating band: the comparison deliberately keeps 5–250 kHz and adds no disjoint random-band energy.

## Finite-record result: preselected hot spots / bursts

The following are **paired per-record ratios**, 16 independent records per row; brackets are pointwise 95% Student-t intervals. “Height” is the leading event divided by its selected region's own mean. The full six-beam results, SDs and fixed-region contrasts are retained in the CSV readouts.

| beta / loss | Combined / randomized-only height: mean (SD), interval | Same-sign rung1 / rung0 height: mean (SD), interval |
|---|---|---|
| 0 / retained | 0.999748 (0.000780), [0.999332, 1.000163] | 1.33061 (0.28952), [1.17633, 1.48488] |
| 0 / Cu proxy | 0.999913 (0.001528), [0.999098, 1.000727] | 1.15931 (0.27408), [1.01327, 1.30536] |
| .35 / retained | 0.999407 (0.020460), [0.988505, 1.010309] | 1.26429 (0.34191), [1.08210, 1.44648] |
| .35 / Cu proxy | 0.990467 (0.021029), [0.979261, 1.001672] | 1.08748 (0.30864), [0.92302, 1.25194] |

For this case the additional signed-source intervals are inside the predeclared ±10% window; this is not a universal nil-effect claim or a model-uncertainty interval. The same-sign spatial-rung increase is resolved for retained loss at both beta values and for Cu at beta 0, but not for Cu at beta .35. Pressure-plus-thermal **diagonal** energy shares at these leading events are at most 5.88e-6; they are not exclusive attribution probabilities, and the actual cross terms and algebraic signed-only bounds remain in the tables.

In these records, mean actual random work is 0.845–3.004 mW versus 690 mW incident. Mean thermal work is about −0.49 to −0.81 µW: this includes signed work against the total trajectory and can return energy to the heat reservoir; it is not negative thermal component self-energy. The per-record source ledger, not a sum of separately driven arm energies, closes the budget. A/C rates in the summary are **electrical-parent strikes per second** (counts divided by actual exposure); B initiation rates count distinct-origin receiver windows, not individual impacts.

Same leading-event time within each record, beta .35, fixed center point (not a newly selected depthwise maximum), mean across 16 records:

| Loss | Depth, µm | Elastic eV/D quotient | Conditional 4.1-ps proxy, eV/D |
|---|---:|---:|---:|
| retained | 30 | 2.984e-8 | 3.401e-15 |
| retained | 50 | 6.555e-10 | 8.260e-17 |
| retained | 200 | 3.672e-10 | 4.595e-17 |
| Cu proxy | 30 | 1.350e-7 | 7.809e-16 |
| Cu proxy | 50 | 2.986e-9 | 1.782e-17 |
| Cu proxy | 200 | 1.600e-9 | 9.542e-18 |

The frozen depth-dependent D population changes the quotient denominator; these depth contrasts are not pure wave attenuation and neither column is a measured energy yield. `same_event` tables also retain the region-average values and all individual records. The separate flange-static illustration uses full incident retention as a heat ceiling before the small mechanical-work subtraction; it is not added as another heat source to the dynamic ledger.

## All six matched beam cases

Each of rung0 and rung1 has the complete 48-arm ×16-record matrix (768 record readouts per rung), with identical total electrical-parent counts for each paired record. Maximum numerical heat/work balance residuals are 2.10e-14 and 1.95e-14 of incident energy.

Across all 24 rung1 combined/randomized-only leading-height contrasts, paired means range **0.99047–1.01447** and every pointwise 95% interval lies within [0.9,1.1]. Four intervals lie above one: “small at this resolution” is not “zero,” and these are not simultaneous intervals. Same-sign rung1/rung0 means range **1.08748–1.33061**, with 22/24 pointwise intervals above one. At these leading events, the algebraic upper bound on signed-only elastic energy share is **5.97e-6**; random-law feedback and normalization are not an exclusive direct-source energy attribution.

Combined-arm width maps use 24 representative records per loss (records 0/1 in every beam/beta case). Retained-loss x/y spans are 1.83–4.49 / 1.83–5.32 mm; Cu gives 1.66–7.15 / 1.66–4.82 mm. None reaches a foil edge, but **17/24 and 18/24** respectively carry the conservative retained-band width-censor flag. Sampled half-height temporal spans are 2.10–16.60 µs (retained) and 1.40–20.75 µs (Cu); **44/192 and 39/192** combined-arm leading events are observation-boundary censored, so those are not complete physical-duration distributions.

The even-profile frequency comparison is pinned at 100 kHz, the appreciably pressure-coupled **98.626752-kHz** rung1 mode (signed normal participation about +0.479 kg^-1/2 in the frozen gauge), and **100.475986 kHz**, the literal midpoint to the stable 102.325221-kHz mode. It is not an antiresonance selection. `frequencies.json` preserves the selection, and the same frequencies—not a new rung2 tuning—are used throughout.

Raw modal thermal-force phases in `transfer_r*.csv` inherit the frozen eigenvector sign. A gauge-invariant relative force phase is `arg(Q_thermal / F_pressure)`; the signed pressure participation must not be discarded. The 100-kHz temperature phase itself is about −0.832 rad in the inspected 128-mode surface diagnostic. A projected force's extra sign/geometry phase is not acausal heating. Settled unit-peak transfer and actual finite, exposure-normalized beam trajectories are separate quantities.

