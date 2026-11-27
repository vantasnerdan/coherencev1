# Pass 3: direct electron-driven curvature and breathing widths

This pass changes the even site curvature, not the sign of a force or a cold displacement seed. The electronic carrier histories, deposited energy, fixed foil profile, three partition/loss closures and three beams are inherited from Pass 1. The six comparison coefficients are +/−0.001, +/−0.01 and +/−0.1 fractional curvature per eV of electronic excitation per formula-volume, with a zero control.

## Direct-curvature result

The specified scan changes the matched variance by at most +0.2241% (rms width +0.1120%). This is the native-axial cage estimate at 2.5 keV, D₂⁺, fast closure, neighbors on and −0.1/eV. The corresponding bulk maximum is +0.2048% variance. The largest local energy excess is 4.921 meV, in the 0.25-keV D₂⁺ fast bulk case. Every scanned curvature remains positive; the smallest is 0.979894 eV/Å². The exact D-column normalization differs from the accepted input by at most 6.7e−16. The event energy defect is at most 2.44e−7 eV, with a 4.03e−7-eV bound after positive retarded remapping; actual carrier stocks remain nonnegative and no uncertainty violation is resolved.

The 2026-11-27 frozen-epoch bulk example at 2.5 keV, D₂⁺, fast closure, neighbors on has net signed work −0.192616 eV for +0.1/eV and +0.193208 eV for −0.1/eV. Their even-in-coefficient average, 2.9562e−4 eV, is a parity diagnostic: it neither proves that all odd work is recoverable nor identifies the even part with surviving squeezed energy. At 1106 fs the actual covariance-energy changes are +0.00142114 and −0.00142113 eV, while the ordinary optical stock is 0.108464 eV. The changed optical energy dissipated over this window is −0.155505 and +0.155953 eV, routed once through the existing 20% thermal-electron / 80% high-Pd ports; the remaining carrier changes and signed foil port are in the tables.

The local source-clock correction matters: at the same 0.0157139-eV/formula deposit and source area, the selected native interval is 1.0606 fs rather than a 21.9418-fs spatial-cell transit. Its weak local beta at 0.1/eV is 0.0007302 rather than 0.0003702. This is finite-pulse covariance breathing, not a maintained parametric gain.

## What is propagated

The local Hamiltonian is p²/(2m) + q K(t) q/2. A positive electron internal-excitation density, nonthermal plus thermal electrons before the stiffness-work debit, gives

    delta K(t) = K_local * coefficient * U_e(t) per formula-volume.

The volume is a³/4 from the same frozen lattice epoch. It is not 1/n_D in a partially occupied layer. Actual occupied D counts and source energy counts are different weights. This energy coordinate is internal excitation, not electron temperature or thermal electronic free energy. The original Pass 1 generator carries electrons, D optical population, all Pd bands, recoils, stopped structural energy, local lattice flow and the three foil exits.

Eight depth cells and six radial annuli integrate the deposited source vectors. Each sample's pulse starts at the original track's arrival at its depth plus the 10 Å/fs lateral front, and lasts its original local source interval capped by the accepted 2-nm local deposition kernel divided by local projectile speed. That duration is not the full spatial-cell transit. Formula-volume energy still uses the full spatial cell, not a smaller native-row volume. Occupied-D weights use the exact piecewise profile column at cell faces, not midpoint density. D₂⁺ has two half-energy constituents at the actual 0.74 Å/v axial delay, one electrical parent; the sibling observes only the first constituent's prior field.

The frozen Pass 1 operator has a 73.3014-fs local period. Its actual FCC neighbor law is sampled on a 4×4×4 reciprocal cell, with identical frequencies grouped without changing their weights. The independent control sets that neighbor law to zero. Each local reciprocal-mode covariance is propagated by an exact constant-curvature symplectic step. Its position, momentum and symmetrized q-p entries are retained. Fourier characteristics reconstruct relative covariance with 2[1−cos(k·r)], removing common translation. `field_pairs_*.npz` contains spatial projections of the local anomalous-variance field, including cross-projections. These are ordered-width observables, NOT global normal-mode covariances. The local dynamics is a Gaussian Wigner/leading-gradient closure, not an eight-site motif tiled through the track.

For bulk, the G10 300-K quantum moments are q²=0.0269276534 Å² and (p/ħ)²=16.6396110 Å⁻² per axis, with symmetrized q-p=0. They already contain zero-point motion. The exact `joint` files are single-site marginals and an excitation/kinetic quadrature, not an intersite density matrix. For the original six-D cage, the finite native cut gives K_self=(E_A+E_B)/u² and J_close=W/u². Its COM and relative normal modes have K_self−J_close and K_self+J_close; their 300-K quantum priors provide the site and relative widths. The relative observable excludes the COM mode exactly. These are finite-secant axial estimates, not a measured full cage Hessian; `cage_axial_transfer.csv` retains their cut-dependent spread and separate bulk-width comparison.

Those G10 moments are not stationary under the frozen G3 quadratic operator. The local effective-law initial-value problem starts with the native 300-K state at the causal local pulse onset; every driven state and its quiet comparison use that same local clock, state, ordinary optical inflow and loss. The pre-arrival increment is zero. The free native reference is never damped to zero. Ordinary Pass 1 optical energy enters as positive phase-random covariance noise and is already inside the propagated state, not added again at the yield consumer. Changed optical inflow replaces the duplicate delta optical reservoir. The explicit Gaussian bath map is C -> exp(−dt/T1) C + [1−exp(−dt/T1)] C_native_quiet(t), plus the positive optical source. It preserves Gaussian normalization, positivity and uncertainty. T1 is the energy lifetime; separate extra pure dephasing does not become lost energy.

Whole-footprint readout integrates the retarded age inside each unchanged cell, with source/profile breakpoints and positive quadrature weights. It does not replace a whole cell by a single phase. The numerator and matched constant-phase norm use the identical weights and causal support; no interpolation value before the local front is admitted. The reported common captured window needs no late-time extrapolation. The same remapped second moment drives the reciprocal acoustic source.

The phase fraction is undefined when its matched norm is zero; the table leaves it blank. Only analytically empty pre-arrival support is set to zero. Positive-support values, however small, are retained and distinguished from a printed floating-point scale. The resolved maximum is 0.999752; the measured triangle excess is zero. This is phase alignment of a width field, not a phase fraction of particles or a large amplitude.

## Retained global normal modes

`four_mode.py` separately propagates four genuine global Gaussian modes on the frozen 2.5-keV D₂⁺ fast footprint, for both environments, neighbor on/off, every comparison coefficient and the same pulse-off reference. Bulk retains axial cos/sin partners at ±π/L and ±2π/L; the cage retains cos/sin at ±π/L in each native COM/relative sector. The functions span the full occupied footprint. The physical mass Gram is the integral of m_D n_D phi_a phi_b over volume; the geometric volume Gram is also saved. Orthonormalization precedes projection of the accepted G3 bond-gradient Hamiltonian and actual retarded K(t). Frozen 4-fs driver samples are causally interpolated onto 1-fs exact finite constant-Hamiltonian steps; the optical-population interpolant preserves positive reference inflow. The coupled 8×8 position/velocity covariance uses that same reference optical noise and the matched native-reference bath.

The native moment values are the same as above. This conditional joint Gaussian has one global preparation/quiet clock; it is not an exact reconstruction of the local-age IVPs or a measured G10 intersite state. `whole_mode_pairs.csv` gives the actual real and imaginary parts of delta< a_mu a_nu >; `band_*.npz` retains the full covariances, orthonormal functions, Grams, frequencies, driver and signed work. The minimum quantum-uncertainty eigenvalue is +0.02634, zero-drive covariance difference is exactly zero, and the retained-system energy defect is at most 2.18e−14 eV. No new electron debit or energy reservoir is introduced: this is a conditional retained-band consumer of the same frozen drive, not an extra term to sum into the source ledger.

At −0.1/eV with neighbors on, bulk peak anomalous Frobenius norm is 1.3450e−4 with mode mixing versus 6.4380e−5 when off-diagonal projected-drive entries are omitted. The first 1/e crossing after the post-deposition peak is 11 fs in both cases; the local width-field readout gives 332 fs. For the cage the corresponding norms are 3.7987e−5 versus 3.4928e−5, crossing times 303 fs in both cases, versus 336 fs for the local field. Mixing changes the retained amplitude, not appreciably this crossing time. The global band and local field are different observables/joint-clock closures; their difference is not attributable solely to mixing. These are first-crossing diagnostics, not fitted decay constants: the electronic drive still has 1.65%/2.84% of its peak matrix norm at window end. Four modes are 0.0218% of the footprint's axial mode count and do not bound the omitted spectrum or establish a forcing-free lifetime.

## Finite pulses and the conserved work

For a weak pulse returning to zero, the finite pair-creation amplitude is

    |beta| = |(Omega/2) integral epsilon(t) exp(−2i Omega t) dt|
           = |(1/4) integral epsilon_dot(t) exp(−2i Omega t) dt|.

A 15-fs rise with 100-fs decay and fractional amplitude 1e−4 gives exact |beta|=8.85950e−6, versus 2.5e−5 for an ideal sudden step. A 500-fs rise gives 2.90451e−7, not zero. The preserved exact finite-propagation check agrees with its weak Fourier prediction. `controls.json` also compares a native local source clock with the rejected spatial-cell transit clock at unchanged source area and energy/site.

The code projects the signed K history at twice its actual reference optical frequency. The maintained resonant Mathieu criterion, h_threshold=2 kappa/Omega, is a reference for a sustained narrowband modulation; neither the transient peak nor a finite-window coefficient is called sustained gain. The actual normalized covariance amplification, anomalous covariance, phase-aligned second moment, post-source persistence and energy above the matched reference are also measured. First moments stay zero. Ordered breathing is not a single coherent displacement, and a Gaussian covariance alone is not evidence of entanglement.

At each curvature change, work is exactly 1/2 Tr(delta K C_qq), in eV with physical mass and Å/fs units. Signed work debits or credits the actual electron channels. Negative return work is not discarded. The finite oscillator excitation is separate from the reversible initial-state curvature offset. Optical damping returns the signed changed energy to the existing 20% thermal-electron / 80% high-Pd destinations. Every correction then follows the same positive carrier generator and its foil port. `phase_summary.csv` carries event work, returned work, stored excitation, onward dissipation, actual carrier positivity, uncertainty and the independent balance defect. Its diagnostic minima/maxima for curvature, work fraction and conservation are conservative bounds over the seven branches of that physical source case.

## Pairs, foil and yield

`pair_instantaneous.csv` gives the actual time-conditioned Gaussian separation distribution, absolute pumped and matched-quiet probabilities and their signed difference at every supplied native separation, with relative kinetic energy and finite inward-shell flux. Near-track and median-strength sites remain separate with their positions. Negative signed-separation tails are explicitly reported.

The largest signed instantaneous shell changes are:

| Conditional environment, 2.5 keV D₂⁺ fast, −0.1/eV | Separation | Pumped probability | Same-time no-curvature-change probability | Signed excess |
|---|---:|---:|---:|---:|
| Bulk, neighbors on | 2.407824 Å | 0.146594864 | 0.146549494 | +4.5370e−5 |
| Original six-D cage, neighbors off | 1.726075 Å | 0.203289427 | 0.203164351 | +1.2508e−4 |

Both are at the near-track sample r=0.5657 Å, z=51.1694 Å and t=44 fs; the conditional relative kinetic energies are 0.07519 and 0.08070 eV. These are finite-shell, not nuclear-contact, tails. The matched reference keeps ordinary Pass 1 beam heating and removes only the curvature change; it is not a beam-off distribution. At low energy the ordinary population already makes broad tails: in the 0.25-keV bulk path example 511/512 draws leave the native potential chart. The trajectory-minimum CDF therefore retains censoring bounds rather than inventing a far-tail continuation.

`pair_trajectory_minima.csv` is different: 512 paired phase-space draws are propagated through the native directional pair potential, sharing the same explicit reference-bath noise and ordinary optical inflow. Each draw keeps its phase through the whole pump/lifetime window. A path is censored when it exits the native cut; its unknown future mass is retained as an upper probability, not assigned a contact distance. Native cage forces and bulk finite-cut energies determine this semiclassical nonlinear sensitivity. Single-site native non-Gaussian tails and the basis-8 comparison are separately retained in `native_quiet_tails.csv`.

The trajectory initial draws use paired scrambled-Sobol quadrature, not independent measured events. The zero-count IID reference scales are not QMC or Cox confidence intervals. Initial/noise-map support censoring does not determine a unique continuous crossing energy: the path kinetic diagnostic is the recorded deterministic-exit contribution per censored draw. The complete time-conditioned kinetic value at a supported separation is in `pair_instantaneous.csv`.

`whole_foil.csv` counts one maximum per physical electrical parent, across all six energies, three closures, three beams and both environments. Its compatible overlap uses the actual spatial support and selected time window. `parent_windows.csv` contains the actual center/hotspot Cox event marks and shared-parent ordering. `strength_and_support.csv` keeps occupied-site strength quantiles distinct from event counts and exposes the part of a Gaussian distribution outside the finite well chart.

The beam-on-target consumer uses the actual local CSDA projectile speed and the breathing momentum covariance at the encounter age. Struck D uses the accepted four-transfer recoil rays, their local slowing and actual arrival clocks. The same-parent sibling is separate from other-parent pumping. The native sheet is retained once in the unchanged reference and is not itself pumped. The ordinary target-motion baseline is subtracted once, not added to another zero-point estimate. `yield_mixture.csv` joins both species and distinguishes per incident D, per electrical parent and per true charge, beside the unchanged G12 residual and its 12-h/2-h exposure boundary. The exact relative-flux kernel at the largest computed covariance change is compared with the small-change consumer in `consumer_checks.json`.

The p/static40 column is the physical G12 comparison here. Neutron and bare columns are separate model alternatives; an envelope over those controls is not an observed proton/neutron mixture or a fitted join.

The largest actual p/static40 joined change is 1.42652e−29 per incident D, or 1.64049e−29 per electrical parent, at 2.5 keV, fast cage, neighbors off, +0.1/eV. It is 1.11666e−14 of the unchanged G12 residual. All three beams were consumed: their maximum fractions are 1.1166607e−14 (even), 1.1166587e−14 (hot steady) and 1.1166498e−14 (hot bursts). The six-setting envelope is 5.91e−30, 9.15e−24, 1.65e−18, 1.68e−16, 2.48e−15 and 1.12e−14 of that residual. Thus this specified curvature scan does not account for it in the supplied non-neighbor consumer; this conclusion uses the actual covariance/time/velocity runs, not a mean-energy assertion.

Compatible local-overlap upper estimates reach 1.37e−5, 6.58e−5 and 3.70e−4 for the respective beams. The corresponding longest-window marks are 536, 499 and 427 parents with zero observed prior overlap; those finite Cox samples are not a universal zero-overlap bound. Some conditional parents exceed 0.1/1 meV of site excess, none exceed 10/100 meV in the scanned branches. Counts are maxima per electrical parent, not sampled sites; each table carries its actual 12-h/2-h exposure.

Neighbor-pair inward crossings and their quiet-subtracted changes are in `neighbor_consumers.csv`. Unfitted contact-transmission sensitivities 0, 1e−12 and 1 are carried separately; none is added to the joined non-neighbor yield as a measured nuclear probability or a fitted explanation of the residual.

## The same event's acoustic side

The accepted strain slope is dK/d(linear hydrostatic strain), not a volume derivative. The reciprocal hydrostatic stress change is

    delta sigma_hydro = n_D * (dK/d epsilon_linear) * delta<r²> / 6.

The equilibrium structural strain and a propagating strain are different outputs. A causal elastic mode with transit R/c and an explicit outgoing-radiation port is driven by that signed stress. Its work equals stored elastic energy plus outgoing elastic energy. That work is reallocated from the covariance energy, not added as a second electronic input; the calculated stiffness-backreaction ratio is reported. The saved G13 same-band finite-foil Green coefficients then carry the new outgoing source through the two accepted loss alternatives. No G13 producer is run and no already counted thermoelastic/recoil launch is added again. The clamped-ion electronic pressure needs a same-state strain derivative of electronic free energy; K(U_e) by itself does not supply it.

The largest conditional propagating volume strain is 2.8481e−10 (1-keV D₂⁺ fast bulk, neighbors on, −0.1/eV), versus the same event's quasi-static 1.1153e−8 structural estimate. Peak reciprocal pressure is 2153 Pa; acoustic work is 4.7186e−12 eV, stored plus outgoing, with no energy cap activated. The maximum feedback/direct stiffness ratio is 2.11e−6. This is a conditional Born readout with explicit covariance-energy reallocation, not a self-consistent electronic-pressure/elastic solution.

## Inputs, files and replay

The base is `2d85230f2580efbee4d7de5e8aece2c9540b3c1c`. `source_index.json` binds all consumed frozen inputs by SHA-256; `inputs_ledger.md` gives their meanings. With Python, NumPy and SciPy, from that checkout:

    python g15_track/pass3/run.py

This regenerates the source covariances and every actual consumer. To consume the supplied evidence arrays without replaying valid source dynamics:

    python g15_track/pass3/run.py --consume

Three energy-group evidence archives contain all `cov_*.npz`, local-field `field_pairs_*.npz`, retained global `band_*.npz` and physical-parent windows, each file once. Each archive extracts independently into the same root and has a distinct `manifest_<group>.json`; every manifest shares a source-epoch SHA-256 derived from the input-binding digest and the sorted hashes of all frozen source covariance arrays. Reference covariance is binary64; small signed local differences are stored separately in binary32 rather than rounded into the large quiet baseline. Variance ratios are <q²>/<q²>_matched; take the square root for rms width. The main archive contains source, pages, decisive tables, input bindings and governing plan. External `delivery.json` and `remote_delivery.json` record archive hashes and actual unpublished-draft byte readback. No branch is pushed or release published.

The bounded retained-band consumer alone can be reproduced from supplied evidence with `python g15_track/pass3/four_mode.py`. It does not replay the source, pair paths, beams or yields.

The remaining physical uncertainties are the matched hot-electron q/occupation/temperature/internal-energy curvature join, the quantum six-D cage state and transverse Hessian, and continuation from finite native shells to contact together with actual cage abundance. Those are inputs, not fitted reaction terms. The supplied native non-Gaussian prior and finite support are not replaced by claims about a Gaussian far tail.

The measurement that would most directly distinguish the coefficient branches is a same-state, q-resolved electronic-temperature/occupation Hessian together with internal excitation per formula-volume; time-resolved diffuse-scattering width and momentum-correlation breathing would distinguish a finite filtered pulse from an assumed sudden or maintained modulation. No applicable G16 slice was supplied at source freeze.
