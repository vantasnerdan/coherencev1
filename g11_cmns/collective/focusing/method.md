# Bulk PdD focusing: equations, selection and reproducibility

## Scientific contract and inputs

This is a classical, finite-energy, frozen-Pd Hamiltonian comparison in fully occupied bulk PdD, not a new electronic calculation, pump-field model or finite-temperature ensemble. The given 0.6 meV/D means an already prepared pumped energy budget; no guessed pump efficiency or continuing field is added. The other budgets are 4.8 and 25 meV/D. Each launch has total excess energy N e above the classical x=v=0 reference. The quantum ground energy is **not** included in this classical budget.

The physical D mass is inherited from `g10_quantum_site/model.py`; the kinetic coefficient is in eV fs²/Å². The actual original register has D25 and D63 as rows 12 and 31, zero periodic image, n=[110]/√2, reference R=2.88782409436586 Å. All cells retain those physical target rows. The 2³ conventional FCC cell has 32 D; tiling its geometry gives 4³/256 and 8³/2048. Lattice constant is 4.084 Å. `input_manifest.json` binds every imported prior input by path, SHA-256 and bytes; no previous large state arrays are copied into this capsule. Its revision labels are supplied provenance, distinguished from byte-equality joins actually checked by the packager.

## One energy and its analytic force

For site displacements x_i,

`V = Σ_i V_080(x_i) − Σ_(i<j) J_ij x_i·x_j + Σ_nn [W(q_ij) − J_nn q_ij²]`,

`q_ij² = max(n_ij·x_i,0) max(−n_ij·x_j,0)`.

`V_080` is the complete accepted .80 Å local well, including its analytic cubic angular gradient, reproduced by evaluating the same Hermite coefficients. Below .03 Å it is the old regular polynomial, with K=0.590658682890108 eV/Å² and F0=0. This campaign applies that accepted one-site closure at every D. The prior equilibrium one-site reweighting did not itself validate such an all-site dynamical model.

The old scalar folded graph is preserved, not refitted. For each folded neighbour, enumerate **all equal shortest periodic images** and divide the original J by that image count. The resulting translation-invariant finite-range stencil is applied on the actual FCC register. Refolding at N=32 exactly restores every original J element, including signed distant shells. At larger N this is one explicit lift of the folded model, not independently measured infinite-crystal force constants. Time stepping uses edge lists (63,488 undirected harmonic edges at N=2048), not a dense matrix or all-pairs distance calculation.

The native five-row connected W is positive and has positive signed dW/du: it adds approach cost, not attraction. On the symmetric isolated cut x_i=q n and x_j=−q n, the original harmonic bond contributes **+J_nn q²**; replacing it with W therefore gives exactly `2 V_080(q n)+W(q)`, without double counting. Other displacement rays are an assumed off-path closure, not native angular information.

At q<.12 Å use

`W(q)=J_nn q²+A(q/.12)^4+B(q/.12)^6`,

where `A+B=W(.12)−J_nn(.12)²` and `4A+6B=.12 W′(.12)−2J_nn(.12)²`. These two coefficients are determined, not fitted. Above .12 the native energy/derivative Hermite interpolant is used. Thus the correction has no quadratic Hessian and is C1 at a one-site limit, retaining the independently measured harmonic limit. A naïve first Hermite interval would instead create a negative Γ curvature; that inconsistent promotion is not used. The angular correction vanishes when either site has no positive inward projection.

The scalar small-signal curvature is `H(k)=K−Σ_d J_d cos(k·d)`. `consumer_proof.json` records the complete finite-cell FFT spectra and direct edge-force equivalence. The angular weights are nonnegative (Schur's inequality gives w100≥0, q≥9p gives w110≥0; w111≥0) and sum to one. The proof script evaluates all derivative roots and endpoints of each native radial cubic minus K r²/2, and of W−Jq². Their nonnegativity plus the positive finite-cell spectrum establishes V≥0 on the declared domain and hence **K_rel≤N e** for the conservative supported trajectories. This is a finite-model analytic/numerical argument, not a formal outward-rounded error enclosure or a bound on the real material.

## Launch families and bounded selection

See `patterns.md` for the exact preparations. All candidate metadata, selected launches and time series are retained, not only the winning peak.

* **Converging waves:** up to 16 periodic longitudinal reciprocal modes with opposite phases about the target. A reciprocal-space null constraint sets both target velocities to zero. The harmonic target adjoint proposes phases; the full nonlinear energy/force evolves the launch.
* **Collisions:** two compact, approaching localized packets, target initially silent. Their kinetic participation is eight D in this register. This deliberate preconcentration is reported, not called an extensive distributed state.
* **Reverse spreading:** start with all N e as outward target relative kinetic energy, evolve for 1200 fs or to the last supported step, reverse all velocities and propagate back. This is a genuine energy-conserving time reversal with nonzero displacement phases, not a velocity-only approximation. Exact recovery of N e is a reachability control. It is not collective evidence unless the actual reversed launch is spread and not target-preloaded. Both kinetic IPR and positive-site-energy IPR are reported; the latter uses K_i+V_080(x_i) and explicitly excludes signed bond/correction energies. A turning-point kinetic IPR alone cannot establish spreading.
* **Direct optimization:** target velocities and initial displacements are zero. Harmonic adjoints propose three focus times from a 20-fs grid on 120–2280 fs, separated by ≥160 fs. The exact capped linear amplitude optimizer has `max |v_i|≤2 sqrt[N/(N−2)] sqrt(2e/m)`, excluding target rows and enforcing the fixed total budget. Every proposed launch is evaluated nonlinearly. Four ±12° polarization refinements and four amplitude/phase mixtures of two proposals are evaluated nonlinearly; choose the best supported target peak. This is a bounded local search in a declared distributed family, not a global optimum.

The comparison takes the highest supported pre-exit target peak among the declared candidates, preferring genuinely spread candidates when any exist. Headline eligibility requires launch positive-site-energy participation ≥N/4 and target positive-site-energy fraction <5%; all failures remain in the data. For velocity-only wave/direct states this is the ordinary kinetic IPR. The whole-cell maximum is measured **on that same chosen trajectory**, not separately optimized. Its radial/separation-at-peak columns are null rather than mislabelled target values.

## Dynamics, support, hold and controls

Velocity Verlet has dt=.5 fs and a fixed 2400-fs scoring horizon. Reversal uses its actual preparation time. Every integration step checks all site radii and all nearest-neighbour separations; the trajectory stops on the first |x_i|>.80 Å or R_ij<1.28782409436586 Å. The first crossing sample is never scored. A smooth local continuation is evaluated solely to make that crossing step defined, never as physical post-exit evolution. The support box is an assumed modelling extent, **not angular validation**. The native .48/.80 connected rows carry the original PAW-overlap flags; no flag is converted into an electronic error bound.

`candidate_supported` means the scored peak precedes any crossing, while `trajectory_supported` means no crossing in the full requested horizon. Both are explicit. Peak total relative kinetic energy is `m|v_i−v_j|²/4`; the same-time radial value is `m[(v_i−v_j)·Rhat]²/4`, with signed inward speed separately recorded. Neither is an incident beam energy. The hold is the contiguous sampled interval around the selected peak at ≥half that peak. Endpoint-censored intervals are lower observed holds, not completed durations; dt=.5 fs gives a finite-grid endpoint resolution. A later independent pulse is not added to the hold.

Two matched harmonic counterfactuals use the *same selected nonlinear direct launch* at N=2048 and e=.6/25 meV, with the same support checks and horizon. They distinguish dispersion from finite-amplitude defocusing/clipping without optimizing a new control. Only one damping sensitivity is run: N=2048, e=25 meV, same launch, `dv/dt=−v/4136.03694 fs` implemented by a symmetric split. Its energy decay is intentional, not conservative drift. No loss fit or sweep is performed.

One settled proof runs the integrated 72-row consumer, force/sign/refolding/spectrum checks, exact native-cut role, positivity argument and the **single** dt=.25-fs refinement of the selected N=2048/25-meV direct launch. That refinement also samples the updated local wall and K r²/2 at the actual most-displaced site's direction at every step. It answers whether time resolution changes the largest-energy conclusion and whether the sampled wall is harder than harmonic. Other sound trajectories are not replayed.

The selected .5/.25-fs peaks are .98026129798/.82574185077 eV at 2367.5/1563 fs, with 18.5/18.75-fs holds; neither trajectory exits support. Their signed inward speeds are +.13707240/−.12580597 Å/fs: the latter is receding. This is a changing late-pulse selection, not clipping, and is not a converged high-energy anchor or an incoming-state identification. A provisional 3% assertion initially aborted before serializing the comparison; the same authorized configuration was repeated only to recover its numbers, with no new optimization or refinement level. Energy drift decreases .00026554→.00006527 but does not certify the long-time peak. Both exact rows, ±30-fs local-window comparisons and the saved fine trace are retained.

Production diagnostic arrays have the original ten columns in `results.json`; the saved single refinement additionally has local and harmonic V at the most-displaced site. The final reusable engine also emits those last two diagnostic columns on reproduction; they do not alter forces, scoring or support. Initial positive-site energy diagnostics are computed directly from the saved launches, not by replaying trajectories.

## Quantum/G12 and extrapolation boundary

`rate_mapping.md` defines the read-only adapter. For illustration only, assign the classical K_rel as **excitation ΔE above each retained radial model's own quantum ground**, so the level comparison uses E_ground+K_rel, not K_rel as absolute eigenenergy. The two adjacent retained stationary l=0 levels are shown for every original contact variant and both joins. Their normalized contacts in Å⁻³ map to emitted p/s, p/true C and dimensionless PdD .25-keV k0/k1 **for one continuously active normalized pair** using the existing coefficient and current convention. These endpoint diagnostics are neither a population nor rigorous bounds on a launched quantum state. No hold×stationary-rate event yield is inferred, and no detector/foil rate or high-energy contact continuation is made.

Energy-budget lower bounds for 100 eV and 1 keV assume perfect concentration. Finite-size power-law slopes and optional constant-efficiency scenarios are separately labelled extrapolations; three finite sizes do not establish a limit or absence of saturation. The retained contact excitation ends near 7.9 eV, so 100-eV/1-keV G12 remains unresolved even if a classical scaling law is assumed.

## Reproduction

Use an existing scientific Python environment with NumPy, SciPy, Numba and Matplotlib; no package install or environment mutation is part of this capsule. From this directory:

```sh
python run.py --run
python analyze.py
python verify.py
python analyze.py
python package.py
```

The second analysis only incorporates the settled proof into the page; it does not simulate. `launches.npz` stores only the original campaign's selected launches; `trajectories.npz` stores 2-fs diagnostic samples, while scoring and support use every .5-fs step. Capsule packing is whitelist-only with a fresh safe extraction and byte equality; the parent owns actual cross-host delivery.

## Ordered continuation on this same page

The added `continuation.py` does not import or call the contact/rate adapter. The old comparison, arrays and contact tables are retained unchanged. It carries out these decisions in order:

1. At .6 meV/D extend **only waves and direct** to 16,384 and 131,072 D, at .5 fs, 2400 fs and the original support censor. `LargeLattice` preserves the original tile order, shortest-image stencil, target and exact finite-cell dispersion. The wave preparation still has at most **16 reciprocal modes**; this launch-subspace cap is not a truncation of the lattice dispersion. Three wave trials and exactly the original 3+4+4 direct trials are evaluated. Candidate declaration order and first-wins ties are preserved; refinements follow the three direct proposals. Old small-N results are reused, not simulated again.
2. At N=2048/25 meV/D reselect the direct launch at **.25 fs**, using the same three harmonic proposals, four phase rotations and four mixtures. The base winner and hence the dependent trials are chosen anew at that timestep. This is a fine-step selected result, **not a convergence certificate** and not the already archived same-launch comparison.
3. Only after that result, at N=2048 and 4.8/25 meV/D evaluate the full direct search at .25 fs with the original .8-Å/1.28782409436586-Å support stop lifted, separately under two declared artificial exteriors. Every step still checks the original support; all enter/reentry transitions are retained. Selection uses the full model-valid extended trace, not its supported prefix. `peak_after_original_exit` is distinct from the actual peak/neighbourhood being inside the extent after a reentry.

For `delta=r-.8>0`, the soft local rays use `E(.8)+S(.8)*delta+K*delta²/2`. The stiff rays use the accepted `wide_quartic` exterior: add `K4_ray*(1.5*.8²*delta²+.8*delta³+delta⁴/4)`. Both retain the same angular interpolation and its analytic gradient. Below .8 Å the radial functions, angular closure and forces are unchanged, not refitted.

For `delta=q-.8>0`, both pair exteriors use `W(.8)+W′(.8)*delta+C*delta²/2`, retaining the exact same projected `q` and force chain rule. Soft `C=max(2 Jnn,.1 eV/Å²)` and stiff `C=max(W″(.8−),4 Jnn,4 Csoft)`. Thus both curvatures are positive, soft is no smaller than `2 Jnn`, and neither silently extrapolates the final cubic. Actual numerical coefficients are saved in `continuation_results.json`. These are C1 **model choices outside the computed electronic range**, not new PdD force data.

An extended branch stops separately on actual NN separation <.3 Å, reversal of the reference-axis NN lineage, nonfinite energy or >1% energy drift. These are declared model/numerical failures, not material thresholds or physical peaks. A failing branch is not preferred over a complete model-valid branch. Original support crossing, model failure and finite-horizon censor are separate fields.

The retained geometry includes actual Cartesian displacement vectors for all D simultaneously at each selected peak, ±10-fs neighbourhood, largest radius, smallest NN separation and first original exit, plus the all-candidate largest-radius and smallest-NN states. JSON names the lattice row, original D index, tile translation, unit crystal direction, actual target and minimum-NN vectors, actual R and effective symmetric q; `continuation_states.npz` retains full simultaneous x/v arrays. An off-axis state is never replaced by invented collinear coordinates.

`continuation_verify.py` is the single settled changed-consumer check: exact on-domain force equality, shared on-domain trajectory equality, synthetic C1 boundary/outer-force checks, and an actual selected-state force gradient plus the full extended trace's measured energy drift. It does not reselect again, run another timestep campaign or claim convergence. The selected traces retain every integration step. Raw nonselected candidates retain all-step summary/extrema/support-transition records and compact 2-fs traces with peak/extrema/crossing samples, labelled in the JSON; `peak_index` indexes the saved trace and `peak_integration_step` names the original all-step position.

Each case declares its trace columns: original supported trajectories retain the engine's 12 diagnostics; extended trajectories retain the first 10 physical/support/energy diagnostics only. No missing local-wall diagnostic is represented by a numerical zero. The all-candidate crossing states additionally retain full simultaneous vectors at each such candidate's largest-radius step, including the noncollinear phase-2/3 states.

After unpacking the continuation capsule at its root, enter `g11_cmns/collective/focusing/` and use the existing scientific environment:

```sh
python continuation.py probe
python continuation.py size16
python continuation.py size32
python continuation.py fine
python continuation.py outside
python continuation_verify.py
python continuation_package.py
```

The capsule includes exact prior input bytes at the imported root-relative paths, including the original hidden staging paths; no baseline source identity is rebound. `continuation_manifest.json` hashes the source/data whitelist, excluding itself; the archive hash is reported separately to avoid a circular identity. No push, publication, damping sweep or electronic/contact calculation is part of this continuation.


## Actual nonlinear launch feedback

`nonlinear.py` adds a new search; it does not rewrite the accepted eleven-candidate banks, engine, input curves, old comparisons or contact/rate adapter. The cases are N=2048 at .0048/.025 eV/D and dt=.25 fs, and N=16384 at .0006 eV/D and the original dt=.5 fs. Every scoring horizon is 2400 fs. The saved accepted launches seed each objective independently. The 16k mixture is reconstructed from its stored proposal times and **stored** base rankings, not by repeating or changing the old bank selection.

The admissible launch has x(0)=0, both target velocities zero, `sum_i |v_i|²=2 N e/M`, and

`max_i |v_i| ≤ CAP = 2 sqrt[N/(N−2)] sqrt(2 e/M)`.

The existing capped linear optimizer gives the nearest Euclidean point on this fixed-norm/site-cap intersection: minimizing `|u−z|²` at fixed `|u|²` is exactly maximizing `z·u`. Its null-field tie convention remains unchanged. Kinetic IPR, realized total energy, target preload, sphere residual and cap residual are saved for the actual selected launch. No target displacement or velocity is preloaded.

There is also a direct projected-step witness. For feasible v, `z=v+alpha*g` and its nearest feasible u, comparison with the feasible competitor v gives `2 alpha g·(u−v) ≥ |u−v|² ≥ 0`. This is an algebraic statement about the proposed direction, not proof that the finite nonlinear objective improves. The actual propagated acceptance test remains authoritative, especially with an approximate force-Jacobian pullback and a moving peak time.

Three **different** full-trajectory readouts are maximized on each launch: total `Krel=M|v_i−v_j|²/4`, strictly inward `Kin=M max(0,c)²/4`, and approach `Rref−R` (equivalently minimize actual R), where `c=(v_i−v_j)·Rhat`. Receding radial kinetic energy is never counted as inward. Each selected state reports simultaneous Krel, Kin, signed c, R, time, all-site extent and its contiguous half-objective hold. The closest approach and energy maximum need not occur together; a half-approach hold is not an energy-pulse duration.

### Discrete gradient and actual acceptance

The adjoint differentiates the executed nonlinear split, not the harmonic proposal kernel. Write `h=dt/2`, `k=h/M`, `d=exp(−dt/tau)` (d=1 without loss):

`w=v+k F(x); x'=x+h w+h d w; v'=d w+k F(x')`.

For endpoint duals a_x', a_v', first add `k DF(x')ᵀ a_v'` to a_x'. Then `a_w=d a_v'+h(1+d) a_x'`, `a_x=a_x'+k DF(x)ᵀ a_w`, `a_v=a_w`. The force derives from the same energy, so its locally differentiable Jacobian is symmetric; centered directional differences of that **actual force** supply the Hessian-vector product. Adjacent kicks at an identical x share the pullback. The adjoint uses checkpoint replay. At a force-Hessian knot a unique derivative is not asserted.

`adjoint.endpoint_gradient` returns the unprojected initial-velocity gradient and simultaneous endpoint x/v. Its smooth soft/stiff force branch does **not** confer support. The search removes the sphere-normal component, projects trial launches with the original cap and checks their actual nonlinear trajectories. An optional fixed-time prefix prefilter can reject a nonimproving proposal; every improving proposal is completed to the same horizon or its first original exit, rescored on all steps, and retargeted to its new actual maximum. Joining prefix and remaining propagation preserves the original split and elapsed physical time.

The first unsupported sample is never scored. A later support exit does not invalidate an earlier supported peak, but its truncated hold and horizon remain explicit. Soft/stiff exteriors are available only as separately labelled artificial scenarios; none of their post-exit states may enter the supported ranking. The optimizer does not silently change forces to obtain a derivative.

Accepted histories and feasible rejected steps distinguish a measured local plateau/backtracking stall from a bounded stopping decision. The declared meaningful relative-gain threshold and the larger-step plateau challenge are saved with the trials. A search without a demonstrated local stop is reported as **best observed, no demonstrated plateau**, never as convergence or global optimality. A supported edge can be an active constraint; it is not a material ceiling.

Final delivery reranks **every observed feasible conservative trial**, including rejected steps, across all three objectives. A separate `nonlinear_select.py` recovers a genuinely winning unsaved launch from the exact stored seed, accepted-prefix gradients and trial scale, then completes its actual trajectory. It does not repeat the original candidate bank. Independent search histories remain intact. A cross-objective winner is labelled best observed without an objective-matched plateau unless the active inward radial K pulse and endpoint gradient are identical to the total-K case; only then can the same rejected-step evidence be reused.

### One loss value, no condition drift

Exactly one new amplitude decay, tau=4100 fs, is inserted once per step between the two half-drifts. The original seed under loss, conservative optimum propagated under the same loss, and newly optimized launch under loss share the same initial total N e, cap, dt, horizon and support. Post-hoc `exp(−2t/tau)` scaling of the conservative kinetic trace is shown only as a different, non-dynamical comparison; it does not damp the nonlinear trajectory or supply its R. No thermal, electronic or rate input changes. The older 4136.03694-fs sensitivity remains part of the immutable prior result, not the decay used in this new comparison.

The actual-loss comparator is propagated from the **delivered cross-selected conservative K winner**, not accidentally retained from a smaller independent K-search finalist. Loss-search seed lineage remains explicit when those launches differ. A staged comparator is reused only for an identical launch, dt, tau, horizon and frozen-model fingerprint; this does not substitute post-hoc scaling for propagation.

The settled `nonlinear_verify.py` exercises the actual changed nonlinear consumer and saved objective readouts; `nonlinear_package.py` verifies an exact fresh archive roundtrip and cold extracted saved consumer. These are finite-model implementation witnesses, not timestep convergence, new electronic support or a physical emission/rate calculation.

```sh
python nonlinear.py N2048_4.8meV --no-prefilter
python nonlinear.py N2048_25meV
python nonlinear.py N16384_0.6meV
python nonlinear_select.py
python nonlinear_verify.py
python nonlinear_package.py
```

`nonlinear_manifest.json` separates this delivery from the exact necessary external replay bindings. `nonlinear.tar.gz` is one unpublished local capsule; unrelated old results, rate data and prior archives are not copied into it. Physical times and timesteps above define the calculation.

Selected finals reference their exact arrays in the corresponding case NPZ rather than copying a large trajectory for each objective label. `nonlinear_selection.npz` carries only the recovered unsaved trial states and the matched-loss comparator traces; `selected_result` is the shared resolver used by verification and the extracted consumer.

In a fresh capsule without the intentionally omitted original outputs, use `python nonlinear_verify.py --replay-only`: this replays the new selected states and changed consumer, while explicitly not rechecking the absent old-output bytes. The recorded source-workspace proof retains that separate byte-preservation result.

Exact selected-state reproduction uses the archived launch arrays and the verifier. Restarting the optimizer is not required and can produce a different bounded-search endpoint; the recorded settings, histories and stopping reasons identify this result.

The new commands have all necessary replay bindings in the delta capsule. Earlier sections, their original commands and their figure/data links refer to the unchanged accepted baseline tree; overlay this capsule there to render the complete same page.


## Inward-first continuation at the pumped budget

`inward.py` retains .0006 eV/D, x(0)=0, both target velocities zero, the existing site-speed cap, dt=.5 fs and the 2400-fs horizon. At N=131072 it reconstructs **only** the stored `direct_mixture_-0.25`: the stored best base at 2260 fs minus .25 times the first other base at 1940 fs, followed by the original capped projection. The 2100-fs base and its stored ranking are retained in the provenance; the old candidate matrix is not propagated again. One actual matched initial trajectory checks the reconstruction against its saved all-step trace. The initial objective is its strongest positive-closing pulse at 1781.5 fs, not a new harmonic focus time.

At N=16384 the initial state is the exact previously selected inward launch, resolved from `nonlinear_selection.json` and its referenced arrays. Its .2946611352565502-eV pulse at 1551.5 fs is the starting point, not a replacement bank or the older .12483784-eV receding total-energy maximum.

The scalar optimized on every trajectory is `Kin(t)=M max(c(t),0)^2/4`, where `c=(vi−vj)·Rhat=−dR/dt`. The discrete adjoint targets the strongest observed closing sample, then each feasible projected trial is propagated with the original nonlinear force. A prefix ending at that sample can reject a trial early; an improving prefix must be completed to 2400 fs or the first original support exit and rescored on every retained integration step before acceptance. A rejected prefix is not evidence about unobserved later pulses or another basin. This is a local search, not a complete search over feasible launches.

The first unsupported sample is excluded. The score, total relative K, signed c, actual R, half-inward-energy hold and simultaneous extent are kept together; the minimum R is a separately selected event. Full-horizon support and an otherwise valid peak preceding a later exit are distinct recorded facts. No selected supported score is allowed to use a post-exit exterior. Consequently a support-limited branch is not evidence for a material focusing ceiling.

One actual N=131072 feasible small-step directional witness compares the native fixed-time inward gain with the adjoint prediction at the restored strongest closing pulse. The settled extracted consumer recomputes the saved all-step readouts and launch constraints; it does not repeat the optimization.

After each conservative search, the exact winner is propagated with the single amplitude decay tau=4100 fs inside the original split. Its strongest inward pulse is the matched passive-loss comparator. At 16k the exact previous loss-designed launch is also retained and compared by its inward pulse, not its possibly larger receding total K; the better feasible start seeds the loss-in-dynamics inward search. No post-hoc exponential is substituted for actual lossy geometry.

Search histories distinguish a measured .2%-gain local plateau, challenged by a larger feasible step, from a bounded best-observed result without a demonstrated plateau. Neither is a global optimum or timestep/horizon convergence result. The **selected initial and final launches and their full scored traces**, including support crossings, are retained in `inward_N*.json/.npz`; accepted histories and rejected-trial readouts are separate records, not full state arrays for every intermediate iterate.

The cap is a deliberate distributed-preparation/local-energy constraint, **not a measured maximum material velocity**: each off-target site has at most four times the mean free-site launch kinetic energy, implying participation `IPR≥(N−2)/4`. Active-site counts are computed from the saved launch arrays without additional dynamics.

For the delivered 131k conservative result, only **one verified actual update** is selected; it is not a local stopping certificate. The first feasible step is scaled by the measured projected slope, rather than restarting a long large-step ladder, and the raw gradient plus simultaneous endpoint/support witness are retained before line search. The corresponding 4100-fs loss feedback uses only already available conservative launch increments, projected at the same fixed energy and cap. Each tested direction is scored in actual lossy dynamics; the searched signs and best-observed status are explicit in JSON. A lossy result is called reoptimized only when those measured scores select a different launch from the matched passive-loss comparator. No global search or new force surface is implied.

## Conditional contact on the six prior inward episodes

The [pulse method](pulse/method.md) and [exact results](pulse/results.json) retain the conservative-inward and previous loss-designed launches for all three earlier cases, selecting the largest **positive** signed closing pulse on each launch. The actual sampled minimum R and first receding sample are separate records; sign-bracket interpolation supplies the enclosing outer–inner–outer episode, not an invented quantum state.

An explicitly prepared normalized l=0 packet is propagated in each archived joined radial operator with a spatial exterior drive `−F(t)f(r)`, where `F=μ R̈+U′(R)`, `f′=1` along the supplied classical path and `f=0` near the Coulomb cusp. This makes the path a classical solution without claiming its quantum mean follows it. Both joins and both original flagged-row treatments remain separate. The two declared finite energy bands use the same initially projected packet; they are a sensitivity comparison, not confidence intervals or a measured material energy distribution.

The positive total exposure is `Ξ=∫ρ(t)dt` in Å⁻³ s. Its once-subtracted signed excess over the same undriven ground is recorded separately. Snapshot occupancy integrates the declared Coulomb-model density inside 1 fm and is not a first-passage probability. Proton pulse yield is `1−exp(−Ap Ξ)`, with the unchanged proton coefficient; the either-ordinary-DD-branch result is separate. No peak×hold rate or classical shell-visit count enters any of these quantities.

For a normalized state confined to the declared finite span, Cauchy–Schwarz gives `ρ≤Σnρn`, hence `Ξ≤Δt Σnρn`. The numerical [class ceilings](pulse/class_bounds.json) use the archived contact slopes; this exact finite-space inequality does not prove the material occupies that space or enclose uncertainty in the adopted potential. The observed radius mismatch and contact sensitivity to the wider band leave the physical pulse probability unresolved.

The fixed G12 conversion retains emitted p/s, emitted p/true C and k0/k1 under **one such pulse per incident D**, with the same 85/15 weighting, proton branch and own .25-keV PdD n65. It is distinct from the unmodified old continuously active normalized-pair stationary table. All conditional examples, both bands and once-subtracted excesses remain available rather than selecting a favorable join or cutoff.

## Combined unpublished pulse/size capsule

The new root-relative `pulse_size.tar.gz` is an exact overlay on the accepted baseline. `pulse_size_manifest.json` binds the owned pages, new results, supplied pulse source and the unchanged necessary model/input bindings; its own hash and the archive identity are separate to avoid circularity. The original stationary table and force sources are not rewritten.

From this focusing directory in the existing scientific environment:

```sh
python pulse_size_consumer.py
python pulse_size_package.py
```

These inspect the saved selections and package their exact bytes; they do not repeat optimization, the six pulse histories, the old candidate bank or DFT. The packager performs safe fresh extraction, exact archived-byte comparison and one cold extracted saved-consumer check of the launch constraints, same-time native readouts, finite-class ceilings and proton/G12 units. The 131k native endpoint and feasible directional witness are retained in its case JSON/NPZ. This is finite-model implementation evidence, not timestep convergence or identification of the material's quantum state.
