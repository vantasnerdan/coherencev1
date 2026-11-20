# Resolved lattice-to-reaction coupling

## Levers and present part-(d) result

The new axis is the **joint target displacement/velocity distribution**, not a multiplier by the number of coherent sites. The material-facing levers are relative (not common-mode) displacement, target velocity along the beam, effective charge/response spectral weight, and loss per occupied site. The electronic displacement response remains open.

For the filled PdD geometry, the conditional screened point-charge spectator calculation gives `|delta U| <= 0.0639229 eV` **only if every site's instantaneous displacement is bounded by 0.0432 Angstrom**. An RMS is not such a bound. Pure uniform Gamma displacement of all D, including the struck target, gives exactly zero relative spectator change. The 0.0639 eV envelope is not an electronic-polarization bound and is not a measured locked-minus-unlocked screening shift.

At fixed charge density, curvature and loss, locked and unlocked labels have identical linear susceptibility: their computed contrast is zero. A deliberately generous **full ionic response versus no ionic response** benchmark gives 1.65992 eV at zero speed and:

| D lab energy (keV) | Full ionic addition (eV) | Native target yield factor, ionic benchmark | Positive spectator-envelope factor |
|---:|---:|---:|---:|
| 0.5 | 0.00517950 | 1.00050243 | 1.00626549 |
| 1 | 0.00258811 | 1.00009976 | 1.00246692 |
| 2 | 0.00129364 | 1.00001894 | 1.00092614 |
| 3 | 0.00086234 | 1.00000713 | 1.00051736 |
| 5 | 0.00051736 | 1.00000211 | 1.00024855 |
| 8 | 0.00032333 | 1.00000071 | 1.00012783 |
| 13 | 0.00019897 | 1.00000024 | 1.00006525 |

These are absolute physical yield/C ratios on the **same transported PdD target and outer adsorbed sheet**, with a declared 40 eV electronic baseline, the native D+/D2+ electrical-current fractions 0.85/0.15, 5% spread, stopping, and fixed proton threshold term. They are projections, not a refit of the closed bulk screening table. `response.json.projections` also contains independently normalized `R(E)/R(13 keV lab)` factors; these must not be confused with absolute yield gains. The benchmark is not added a second time to an already ionic-screened baseline.

### Targets, stated once

At the collision reference **2 keV lab**, hold the response poles/kernel fixed and use `U_ion=0.00129364396 eV` as a unit ionic-weight reference:

- **Effective charge:** an added 0.5/2 keV needs approximately **622/1243 times** the reference effective charge if the response weight scales as charge squared; exactly `Z_locked/Z_unlocked=sqrt(1+delta_U/U_ion)` for a unit unlocked reference.
- **Oscillator strength:** it needs additional weight **3.87e5/1.55e6 times** that unit reference. A larger coherent domain does not supply this weight. Charge or pole changes require a new response calculation and sum-rule accounting; this is not an extrapolated physical prediction.
- **Collective polarization:** the required contact-potential change is **500/2000 V** for a +e projectile. For the *same frozen spectator response shape only*, it corresponds to **7.82e3/3.13e4 times** the 0.0639229 eV reference projection. There is no bound here on a different host electronic response.

The linear ionic benchmark saturates with domain size; **no finite domain reaches 0.5–2 keV in that fixed-strength model**. This is conditional nonexistence in this model, not a restriction on a changed electronic/nonlinear/Floquet response. The point-charge triangle envelope reaches 500/2000 eV only at imposed displacement caps 1.329/1.406 Angstrom, approximately 46/49% of the nearest-neighbour spacing, outside the small-amplitude scenario.

### Fourth cause for the plateau table

| Cause | Size for beam fit | Near-surface D x50 | D2+ current fraction 0.85 | Half current | Distinctive signature |
|---|---|---|---|---|---|
| Lattice state | Ionic spectator envelope plus **open electronic term**; equal-pole linear contrast zero, full ionic benchmark above | Change volumetric D only in top 10 nm, retain outer sheet; state contribution also follows loading/coherence | Change constituent energies and count two D per molecular charge; no independent state multiplier for fixed pump/T | Frozen physical yield/C unchanged; a beam-maintained state may relax and needs pump/transport history | Onset/loss, hysteresis or relaxation correlated with loading, temperature, pump threshold and phase-correlation length, rather than instantaneous current alone |

Machine-readable rows are `response.json.projections`, keyed by `intervention` and `E_lab_keV`; `fourth_cause` records the state signature and counterfactual missingness. No DT beam or LiD lattice response has been computed here. The approximate same-D ionic benchmark `delta_U_eV ~ 0.002588/E_lab_keV` can be used only as a labelled fixed-response counterfactual, not as an isotope-recalculated prediction. The DD threshold interference is not transferred to DT.

The collision handoff requests a **paired frozen electronic calculation**: undisplaced PdD versus an imposed Gamma, x-polarized hydrogen-sublattice displacement, **0.0768 Angstrom per D**, with the metal host fixed, at **2 keV lab**. This is the updated primary intervention, rounded from the cold pumped PdD record maximum. The original **0.0432 Angstrom** intervention remains a budget-conditional linearity check; the archived preliminary part-(d) result refers to that original control. This shape is taken from the corrected NiH Gamma mode; it is not a verified uniform locked PdD domain. Include the struck target in the common displacement and reference its collision coordinate consistently, so moving the target does not fake screening. The displacement-induced electronic difference must be reported separately from the ionic spectator term. The parent owns that handoff; this directory does not run or fabricate it.

Forward: the paired calculation is in `g6_collision/`. At 0.0768 A the displaced cell changes the lowering by −0.69 eV at R 0.1 A, and by −0.01 to −1.07 eV over R 0.1-1.6 A.

## Settled matched-case results: (a)–(c)

These are pump-on/off contrasts, **not three verified locked/unlocked phases**. Corrected cold NiH is Gamma-led and locked across its finite box; cold PdD is patchy, and the 300 K partially loaded PdD pair remains unlocked. The phasor amplitudes are 0.0636121, 0.0379761 and 0.0091872 Angstrom, respectively; the last is a thermal finite-window phasor, not the total thermal RMS. Native-step observation covers 20.7266 ps for NiH and 5.77436 ps for PdD. No corrected warm NiH record is substituted from the older epoch.

### Beam motion

Native **proton-channel, fixed-threshold** physical yield/C factors, on identical transported PdD profiles, 40 eV baseline and beam definitions:

| Pump contrast | Beam direction | 0.5 keV lab | 2 keV lab | 13 keV lab |
|---|---|---:|---:|---:|
| NiH H velocities transplanted into DD | parallel x | 1.00771594 | 1.00078684 | 1.00003844 |
| same | transverse z | 1.00016733 | 1.00002601 | 1.00000191 |
| cold PdD, x=1 | parallel x | 1.00078791 | 1.00007827 | 1.00000341 |
| same | transverse z | 1.00001732 | 1.00000269 | 1.00000020 |
| 300 K PdD, x=.875, neither locked | parallel x | 1.00041226 | 1.00002298 | 0.99999782 |
| same | transverse z | 0.99997965 | 1.00000098 | 1.00000068 |

`beam.csv` also retains separate neutron and total-DD rows, the no-threshold background, and normalized shape factors. These are kinematic changes, not additional screening. The named 0.5 keV NiH-transplant contrast ranges **1.00627–1.00838** over eight matched contiguous temporal blocks through the actual target consumer. That is observed finite-time variability, not an independent-sample confidence interval. Central-difference velocity-squared attenuation at the mode is 0.96784 (NiH) or 0.98945 (PdD), an additional declared discretization limitation; no spurious precision is inferred from the printed digits.

### Closest approach and the conditional pair consumer

| Matched pair | Smallest sampled separation off -> on (Angstrom) | Median per-pair minimum off -> on (Angstrom) |
|---|---:|---:|
| NiH cold | 2.63892249 -> 2.60204406 | 2.63892250 -> 2.63866449 |
| PdD cold x=1 | 2.88518596 -> 2.85941272 | 2.88673633 -> 2.88404597 |
| PdD 300 K x=.875 | 2.02792551 -> 2.00177485 | 2.39752509 -> 2.40227260 |

There are 3072 sampled joint pairs for full occupancy and 2677 for partial occupancy, with the same identities in each matched contrast. Near-uniform Gamma motion gives large lab velocity but little relative compression. `pairs.csv` reports a width only when it exceeds the accumulated-second-moment resolution; the nearly motionless unpumped NiH width is unresolved, with an observed-support upper bound **2.36e-8 Angstrom**, not a subtraction-generated width. Raw joint minima/approach arrays are retained. Every observed finite-start support lies outside the 0.287993 Angstrom bare turning-radius maximum, so its reduced bulk barrier fold is one.

The **actual lattice pair rate is not identified by these records**. The properly specified quantum-potential comparison gives:

| Host geometry and independent k=Kii/2 | Isotope | Achieved quantum mean (Angstrom) | log10(rate / s^-1), background | Conditional on/off rate ratio |
|---|---|---:|---:|---:|
| NiH geometry, k=.856123 eV/Angstrom^2 | DD | 2.733683 | -126.02244 | 1 |
| same | DT | 2.733112 | -137.14366 | 1 |
| PdD geometry, k=.295329 eV/Angstrom^2 | DD | 3.034213 | -122.41411 | 1 |
| same | DT | 3.032968 | -133.16079 | 1 |

All four use the explicitly conditional .541 Angstrom screened potential described below, not a measured material pair Hamiltonian. The ratio is one because the same normalized ground state/potential is assigned to each pump label, **not** because an actual pumped pair contact has been measured unchanged. The DT reduced mass changes its penetration even though its nuclear coefficient is larger. The positive-energy DD resonance does not identify a bound-state rate. The consumed `harmonic()` solver SHA-256 is **99007768cb574ae3e612338dabbdfd17c10f598f0d8ee321985ca1dd4110bffc**.

### Pump target

| Cold pumped case | Gross drag (eV/site/s) | Harmonic Omega drag (eV/site/s) | Explicit finite-window pump work (eV/site/s) | Present Omega/gamma | Required Omega/gamma for .0176 eV/site/s, fixed velocity |
|---|---:|---:|---:|---:|---:|
| NiH | 4.6750e9 | 4.7385e9 | 4.7917e9 | 105.97 | 2.815e13 gross / 2.853e13 harmonic |
| PdD patch | 3.6950e8 | 2.9095e8 | 5.1114e8 | 130.77 | 2.745e12 gross / 2.162e12 harmonic |

These are per **occupied site**, not per coherence domain. The finite-window PdD pump-work estimate exceeds the independently reduced drag; transients and estimator/discretization differences are not resolved into a steady balance. They specify the loss target at fixed observed motion; existence/stability along that continuation is not asserted. The DD heat-budget requirements are larger by 17.6/3.65. No nonzero periodic dynamic lock persists pump-off at positive loss in this model.

At 300 K, gross drag is **3.6302e10 -> 3.6338e10 eV/site/s**, whereas thermal injection is **3.7503e10 eV/site/s** in both records. The matched gross-drag increment is 3.63e7 eV/site/s; finite-window explicit pump work is **-2.86e7 eV/site/s**. These do **not** establish negative pump demand, useful extraction, or a steady thermal balance. A thermal/noisy un-locked record cannot be assigned the cold locked-state pump budget or interpreted from its tiny Omega phasor alone.

### Record-specific spectator scenarios

The original .0432 Angstrom cap does not bound these pumped records. `record_limits.json` recomputes the nonlinear triangle estimate from actual observed maxima and case-specific lattice spacing. With the declared .541 Angstrom screening scenario, the cold PdD off/on envelopes relative to the undisplaced geometry are **.00316874/.123080 eV**; their summed absolute contrast envelope is **.126249 eV**, not a signed predicted contrast. The warm pair's much looser envelopes are **8.93483/8.21576 eV**, giving **17.1506 eV** after discarding every cancellation. These bound only the saved discrete records under that charge model, not an unobserved thermal tail or the material.

For NiH, .541 Angstrom is only a shared screening sensitivity, since no NiH screening length is identified here: its cold on-record envelope is .445309 eV under that scenario, **not a material bound**. Actual sampled nearest-neighbour relative displacements remain much smaller than twice global site maxima. A pure common Gamma component cancels exactly.

For the updated **imposed** .0768 Angstrom PdD snapshot, the same spectator scenario gives **.122996 eV**; the .0432 control gives .0639229 eV. Both are recomputed, not linearly scaled. No electronic-response extrapolation is allowed before a returned paired calculation and, for amplitude scaling, a returned linearity comparison.


## Reduced estimators and exact conventions

### (d) Frozen charge and collision-frequency response

The spectator envelope sums `V(R-2A)-V(R)` over FCC sites, `V(r)=C exp(-r/lambda)/r`, `C=14.399645478 eV Angstrom`, `lambda=0.541 Angstrom`, excluding the struck site. This is a triangle estimate for a **difference** in contact potential; constant spectator potentials are not a screening shift. Its positive value is projected only as a sign-optimistic envelope; a real signed electronic or ionic contrast is not inferred.

The conditional ionic model uses `n=4/a^3`, `a=4.084 Angstrom`, `kappa=1/lambda`, `qmax=(6 pi^2 n)^(1/3)`, one unit charge per D, `omega_p^2=n e^2/(epsilon0 m_D)`, and `Omega=2 pi*10.06374152 THz`. Loss is the actual epoch value `gamma=2/tau_amp`, `tau_amp=4.13603694 ps`, not `1/tau_amp`. Put `w=q^2/(q^2+kappa^2)`, `Omega_L^2=Omega^2+omega_p^2 w`. The contact-energy addition is

`delta_U(v)=(2C/pi) integral dq w omega_p^2 w Re<1/(Omega_L^2-(q.v)^2-i gamma q.v)>_angles`.

For `p=sqrt(Omega_L^2-gamma^2/4)`, the angular average is `Re[atanh(qv/(p+i gamma/2))/(p qv)]`; its zero-speed limit is `1/Omega_L^2`. SI q and velocity enter the frequency denominator; integration q is in Angstrom^-1. The finite-domain comparison replaces the lower endpoint by `2 pi/L`, **not the oscillator weight by N**. It is an explicit IR-cutoff comparison, not an exact finite-crystal mode calculation. Angular averaging matters; using a single head-on resonance overestimates the effect.

The native `forward_yields`/`csda_nodes` consumer and its stopping/profile integration are unchanged. A temporary, restored reaction-kernel adapter adds the response at **each local slowing-down energy**, including each molecular constituent, instead of applying a nominal-energy shift to the whole path.

### (a) Moving target and starting separation

For each measured vector target velocity `u` and beam speed `v_b`, use `g=|v_b e_axis-u|`, `E_rel=m_D g^2/4`, and

`effective_sigma(E_lab)=mean[reaction_sigma(E_rel,U_e)*g/v_b]`.

`g/v_b` is the relative-flux weight per incident beam fluence. The full joint longitudinal/transverse velocity samples are retained; independent marginal tails are never multiplied. Fold this kernel through the actual target integrator, 0.5–13 keV lab, with pump-parallel x and transverse z directions, separate DD proton/neutron channels, smooth background (`eq7`) and the existing proton-only threshold term (`proton`). Reaction probabilities and normalized shape ratios are separate columns. The NiH reference contains **H**, not D: its saved H velocities are explicitly transplanted into a DD reaction kernel as a kinematic sensitivity. This is not a NiD or DT trajectory prediction. PdD uses D trajectories directly.

In homogeneous bulk, translating an encounter by the target displacement does not change an impact-parameter-integrated cross section: the impact-parameter integration is translation invariant. Starting a WKB segment at a lattice-neighbour reference adds no tunnelling exponent if the entire observed displacement support remains outside the outer turning point. The conservative observed support is `r_start >= a/sqrt(2)-max|u_position|`; compare it with the largest bare Coulomb turning radius `C/50 eV=0.287993 Angstrom` at the consumer cutoff. This reduction, rather than an invented `lambda(mean distance)` multiplier, handles the bulk starting-separation distribution. Surface displacement/density correlations and a beam channel's spatially nonuniform flux are not identified by these periodic bulk records; the transported target is held fixed.

### (b) Joint nearest-neighbour approaches and quantum contact

Select occupied starting sites once with a fixed seed; retain each occupied pair in all six nonduplicate FCC nearest-neighbour orientations. Stream

`r_ij(t)=|a CELL d + u_j(t)-u_i(t)|`, `rmin_ij=min_t r_ij(t)`.

Arrays retain pair identities, reference orientation, closest distance/time and backward-difference relative speed at that minimum, plus per-pair first/second distance moments and the pooled distance histogram. For minima at the first observation, the speed has no preceding observation and its stored zero is **unavailable**, not a zero collision speed. The inward radial-flux histogram uses `max(0,-delta r/delta t)` and is a classical approach diagnostic, not a fusion probability. Reference orientation is retained; an instantaneous nuclear spin/angular-momentum state is not measured. Uniform Gamma motion cancels exactly from `u_j-u_i`, even when lab velocities and host-relative displacements are large.

Short records do not determine high-barrier nuclear tails. The empirical distribution of per-pair minima is **not** the stationary separation density; neither can fix a quantum contact without a state/potential. Actual lattice pair rates are left explicitly unidentified, not assigned zero and not evaluated at `rmin`.

A proper conditional comparison instead consumes `g6_pairlaw.script.pairlaw.harmonic(lambda, r_center, k, isotope)`, the normalized radial s-wave ground state of

`V(r)=C exp(-r/lambda)/r + k(r-r_center)^2/2`.

Choose an explicit illustrative potential: `lambda=0.541 Angstrom` for both geometries, `r_center=a/sqrt(2)`, `k=Kii/2` (the relative stiffness of two independent fixed host springs). This screening is a shared scenario, not an inferred NiH electronic value. `r_center` is the **external well centre**, not the computed mean and not an instantaneous separation. The same potential is assigned to both pump labels, so its conditional rate ratio is exactly one. This is not evidence that the actual pair rate is one; changed conditional preparation or electronic confinement remains to be supplied. The rate is `A*|psi(0)|^2` per normalized pair per second; DD background includes both branches and the peer's allowed-spin convention, and DT has its own reduced mass/coefficient. A bound-state threshold-resonance matrix element is not determined by positive-energy `reaction_S`, so that rate is explicitly unavailable. The obsolete two-anchor interpolation is not consumed now that the conditional radial model exists.

### (c) Per-site loss and pump budget

Use actual mass and loss: `gamma=2/tau_amp`, `P_drag/site=m gamma <|v|^2>`. **Never divide this by domain site count.** Langevin injection is `3 gamma k_B T` per occupied site for the three velocity components. At 300 K this bath contribution must not be called pump power. `power.csv` reports gross drag, thermal injection, and independently sampled explicit Hamiltonian pump work.

The modulated term is `H_pump=h K_Gamma cos(2 Omega t) x^2/2`, so

`P_pump/site=-h K_Gamma Omega mean[x^2 sin(2 Omega t)]`.

Native-dt positions give central-difference velocity estimates. For a harmonic component their squared attenuation is `sinc(Omega dt)^2`, recorded explicitly; stochastic unresolved kinetic power is not an exact OU-work measurement. The native fitted vector phasor `A` supplies a separate harmonic estimate `P_Omega=m gamma Omega^2 |A|^2/2`. Neither finite-window pump work nor a short phasor is claimed as an exact steady-state total-power balance.

Define ratio 4 as `R4=Omega/gamma`, and also report `f/gamma=R4/(2 pi)` to remove the frequency convention ambiguity. Holding the measured velocity distribution fixed,

`P_drag(R4)=m Omega <v^2>/R4`, `R4_required=R4_current P_drag/B`.

For `1e-9 events/site/s`, `B_DT=17.6e6*1e-9=0.0176 eV/site/s`; the illustrative DD heat budget is `B_DD=3.65e6*1e-9=0.00365 eV/site/s`, **not** the DT Q value. At finite T the gross-drag requirement is a reference in which the pump supplies all kinetic loss, not an inferred net pump requirement. Harmonic-only requirements and matched pump excess are distinct. This continuation does not promise a locked solution at the new damping.

Pump off, at zero temperature, `dH/dt=-m gamma sum_i |v_i|^2`: positive drag cannot support a nonzero periodic locked oscillation. Static vacancy displacements are not pump-free dynamic lock; finite-temperature noise supplies thermal motion, not the imposed coherent pump.

The most sensitive nuclear lever is **ratio 7 plus the independently specified confinement**, not domain N. For pumping, ratio 4 controls the per-site cost and ratio 5 the threshold; relative phase matters for pair compression. The ratio table reports coupling/well ratio, anharmonicity in both supplied directions at the zero-point scale, vibration/thermal energy, both loss conventions, linear Mathieu pump threshold, vacancy fraction, and the declared screening/physical nearest-neighbour spacing. The vacancy-induced shift is a spatial input law, not silently replaced by one scalar. No physical material screening length is inferred from a measured pair separation.

## Reproduction and provenance

From the repository root, using an interpreter with NumPy/SciPy:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python g6_couple/response.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python g6_couple/trajectories.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python g6_couple/couple.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python g6_couple/record_limits.py
```

Run serially. The native corrected `g3_lattice/rerun/part2/lattice.py` observer is used without upstream edits. Only needed matched trajectories are reproduced: corrected NiH cold x=1 pump off/on, PdD cold x=1 off/on, PdD 300 K x=.875 off/on. Within a contrast, seed, occupancy, input epoch, physical duration, noise law and dynamics are matched; observation stride is changed to one native time step. The NiH run has 110592 sites; PdD has 32768 sites. Saved summaries do not contain the joint raw records needed by these estimators, which is why these bounded reruns are necessary. No new material or long campaign is introduced.

`response.json` is the frozen preliminary fourth-cause handoff. `trajectories.json` lists per-case JSON and compact NPZ arrays; JSON binds parameter and engine hashes and retains actual case configuration. `beam.csv`, `pairs.csv`, `power.csv`, `conditional_pair.csv`, `record_limits.json` and `summary.json` are reduced consumers. The native rest-kernel agreement is **1.55e-15 relative**, and the actual named pumped/unpumped contrast changes that same consumer by **0.7716%** at 0.5 keV lab. After a near-zero width subtraction was found unresolved, only the affected pair summary was repaired from saved data; its real consumer check distinguishes unresolved unpumped width from the .0006056 Angstrom pumped radial width. Unaffected beam/quantum evidence and all successful trajectories are reused. No build, lint, formatter or broad test suite was run.

| Input | Authority / epoch | Status and consumer |
|---|---|---|
| NiH springs, quartics, quantum variance | corrected `updated/ni_epoch/parameters.json` | Adopted mixed source epoch; H mass preserved; not promoted to a new convergence claim |
| PdD stiffness/loading/bias/variance | corrected `updated/parameters.json`, native loading and vacancy-bias laws | Existing adopted finite model; spatial laws retained; no new first-principles certification |
| Decay times | same parameter files | Physical damping is an uncertain model input; gamma uses twice inverse amplitude lifetime |
| Transported profiles and adsorbed sheet | `g4_beam.merit.api.reference` | Unchanged within every motion/screening contrast; no refit |
| Beam/reaction/stopping semantics | actual `forward.consumer`, `forward.model` | Lab/CM, charge fractions, two molecular constituents and branch conventions retained |
| Conditional quantum contact | `g6_pairlaw/script/pairlaw.py` | Explicit potential, normalized state and source hash; not inferred from classical minima |
| Displacement-induced electronic response | paired collision snapshot | Open; no fabricated value |

## Three largest omissions

1. **Displacement-induced host electronic/nonlinear collision-frequency response**, including changed effective charge or spectral weight: requires the paired frozen electronic calculation.
2. **Prepared many-body quantum pair density/contact**, isotope-specific state and bound-state nuclear threshold-resonance matrix element: a classical closest-approach record does not supply these.
3. **Pump delivery/efficiency and backreaction, long-time thermal/damage/transport evolution, and high-barrier rare tails** beyond short finite periodic classical records.
