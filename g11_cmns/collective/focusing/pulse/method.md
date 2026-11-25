# Conditional radial contact during an inward pulse

This is a causal example for a **single normalized tagged pair in the spatial ℓ=0 sector**, not a reconstruction of the material's quantum state or a foil rate. No electronic input, original trajectory velocity, join, or stationary readout is refitted.

## What the classical trajectory establishes

The source is the selected full nonlinear conservative **inward** launch and the previously **loss-designed** launch (`loss_K`, 4100 fs amplitude loss in actual motion) for each of the three cases. On each launch find the largest **positive** stored signed closing speed `c = −dR/dt` before the first original lattice support exit. The episode begins at the preceding outer turn, includes the ensuing inward turn, and ends at the first outer rebound. A cubic spline of the recorded radii, clamped to the two bracket-end velocities, places the three zeros of its derivative inside the corresponding archived sign brackets. Actual smallest stored `R`, its index, and the index of the first receding sample are recorded separately. Loss histories are not obtained by damping the conservative paths. Both 25-meV/D episodes end before their first original lattice support exits.

Along a selected classical path the instantaneous quantity `Ecl(t)=μ c²/2+U(R(t))` is not the energy of a freely moving radial particle. `dEcl/dt = Fenv dR/dt` if an extra radial force `Fenv` makes the trajectory respond to the joined `U`. In the trajectories this energy varies across each episode. Nor can even a specified mean energy fix quantum contact: it determines just a diagonal energy moment, whereas `u′(0)=Σ an un′(0)` depends on occupations **and** coherent phases. A classical pair radius and signed velocity are not tomographic data; assigning them to quantum mean radius/current would at most add two moment conditions, not supply the radial state or the expectation of the force over its width. Even the spatial driving field is not identified by the trajectory: a smooth perturbation whose spatial derivative vanishes on the path preserves that classical solution but generally changes the finite-width quantum packet. An inward classical visit is not, by itself, a quantum absorption/hitting event.

For example, the real states `(u0±u1)/√2` have the same mean energy and zero instantaneous current, but contacts `(√ρ0±√ρ1)²/2`. This is not a model of the launch; it isolates information missing from an energy assignment.

## An explicit preparation, force and spectrum

Use the accepted radial operator `H0=−t d²/dR²+U(R)`, `t=ℏ²/(2μ)`, `u(0)=0`, and both existing joins crossed with both connected-row treatments ([contact method](../../contact/method.md)). Units are eV, Å and fs, so the relative mass in the classical equation is `M=ℏ[eV fs]²/(2t[eV Å²])≈104.35 eV fs²/Å²`. The following **externally prescribed** radial force makes the interpolated path a classical solution:

```
M R̈ = −U′(R)+F(s) f′(R),    F(s) = M R̈(s)+U′(R(s)).
H(s) = H0 − F(s) f(r).
```

Here the radial coordinate of the quantum operator is `r` and `s` is time. The exterior coupling has `f=0` below 1.0 Å and `f′=1` above 1.4 Å, including every selected path. In the transition `f′=6z⁵−15z⁴+10z³` with `z=(r−1.0)/0.4`; integrate it from 1.0 Å for `f`. Thus the time-dependent potential leaves the Coulomb cusp unchanged. “Exterior” refers to this **spatial drive profile**, not the lattice search's exterior continuation. The 25-meV/D pulses cross the 1.600 Å radial splice. That accepted splice is only C⁰; exact inverse forcing has a corresponding finite time jump there, which is treated by separate integration segments. No temporal smoothing across these jumps, pre-turn switching protocol, or physical claim about the environment is assumed.

At the preceding outer turn prepare a real, reflected reduced radial Gaussian centered there, with density standard deviation `σ=√(t/ℏω)` before projection. The local nominal width uses `ℏω=√(2t U″(R0+))=√(t K_one-body)` from the archived outward small-displacement curvature; neither the packet width nor its preparation is an inference from the trace. For each history the maximum of stored `K_rad+U(R)` across its sign-bracket samples and all four variants sets a common energy **reference** `Eref`. Project and normalize the Gaussian in the levels `En≤Eref+1.5ℏω`; use exactly the same prepared state for the wider propagation `En≤Eref+4.5ℏω`. The offsets are declared headroom (half a local quantum plus one or four more), not a derived population law or material bound. The cut is in absolute eigenenergy on the existing lattice energy zero, not in excitation above the ground. Within each support solve `iℏ ȧ=(diag(En−E0)−F f_projected)a`. Only retained eigenvectors are rebuilt on the archived grid; the archived Coulomb-matched contact magnitudes are used with the origin-positive basis gauge checked against a fresh match for each retained vector. The projected Gaussian's retained norm, mean radius and current are reported, not forced to reproduce the classical path.

Advance the retained state with two unitary factors per step. With `H1,H2` at the two Gauss times and `α,β=(3±2√3)/12`, apply `exp[−ih(αH1+βH2)/ℏ]` first and `exp[−ih(βH1+αH2)/ℏ]` second. This is the fourth-order commutator-free rule on smooth time segments; no shared order claim is made at higher-derivative breaks of the tabulated force. Steps are at most 0.0625 fs and never straddle the join-force jumps. Simpson integration uses the positive contact density on each segment. The upper-band conservative 25-meV/D, 1.600 Å three-row case tests the exposure and sampled peak contact against halving this step.

## Three different observables and a finite-class ceiling

- **Snapshot occupancy:** `Pr(r<1 fm at s)=∫₀^(10⁻⁵ Å) |u(r,s)|² dr`. The coherent cusp amplitudes use `un(r)=un′(0) r[1+Cr/(2t)+…]` and the energy-specific series for the adopted near-origin potential `C/r+b0+b2 r²`. Integrate the finite-radius density (not density without a volume). The reported maximum is over sampled times of this model pulse, not the probability of ever entering that region.
- **Exposure:** `ρ(s)=|Σ an(s)√ρn|²` in Å⁻³, and `Ξ=∫ρ(s) ds` in Å⁻³ s, with times converted from fs. Integrate positive density across the actual time-dependent episode, separately accumulating before and after the inner geometric turn; neither a peak times a classical hold nor a log-average is used. The signed excess is `Ξ−(tf−t0)ρ0` with the same undriven joined ground subtracted **once**. An excess is not a total pulse probability.
- **DD reaction:** for the archived proton-only hazard strength `Ap=5.81220367098898×10⁷ Å³/s`, the proton-only hazard gives `Pp=1−exp(−Ap Ξ)≈Ap Ξ`; for either ordinary branch add the same-convention `An/Ap=53.701/55.576`, not a second angular factor. If the channels are treated as competing causes of the first reaction, the proton probability is `Ap/(Ap+An)` times the either-branch probability, indistinguishable on the printed log scale. These are conditional probabilities per pair per pulse, with no multiplication by the number of simulated D sites.

For the finite-radius integral write `un(r)=un′(0) r Σj bj (r/a)^j`, `a=1 fm`, `b0=1`, negative indices zero. The regular solution gives `j(j+1) bj = (Ca/t) b[j−1] + ((Ureg(0)−En)a²/t) b[j−2] + (b2,reg a⁴/t) b[j−4]`, where `Ureg=U−C/r` and `b2,reg` is its quadratic coefficient. The implementation retains through order nine and integrates the squared coherent sum by 12-point Gauss quadrature on `[0,a]`.

For every normalized state in a declared retained space, `ρ≤Σn ρn` by Cauchy–Schwarz. The instantaneous contact maximum is attained by the coefficient vector proportional to the positive origin-slope vector. Therefore `Ξ≤(tf−t0) Σn ρn` and the same monotone hazard law gives a finite-class probability ceiling. At finite radius, the analogous snapshot ceiling is the largest eigenvalue of the positive matrix with entries `∫₀^(1 fm) um(r) un(r) dr`, evaluated with the same cusp quadrature. These are statements for a specified support, not an enclosure of either the electronic curve or physical-pulse probability. The [standalone class artifact](class_bounds.json) needs only the archived spectrum and sign windows for the contact ceiling.

The calculation isolates the additional levels in the upper-band solution before squaring and integrating. Their isolated positive contribution divided by the full coherent exposure diagnoses contact sensitivity despite small norm; it need not be ≤1 because omitted interference can be destructive. The sampled integrated current is `⟨i[H0,r]/ℏ⟩`; the predicted `d⟨r⟩/ds` in the finite space includes the extra projected `−iF[f_projected,r_projected]/ℏ` correction. Neither is substituted for an unverified classical match.

## Observed conditional readouts

The times below are the interpolated outer turn, stored maximum closing-speed sample, interpolated inner turn, and interpolated first outer rebound. The minimum radius is a **stored** value, not the radius of the first receding sample. For example, the conservative 16k minimum is 1.924150293 Å at 1569.0 fs; the first receding sample at 1569.5 fs is already farther out. Across all six windows the spline derivative differs from the archived velocities by at most 3.0×10⁻⁵ Å/fs at the samples.

| D sites, meV/D, launch | Outer / closing peak / inner / rebound (fs) | Peak inward radial K (eV) | Stored min R (Å) |
|---|---:|---:|---:|
| 2048, 4.8, conservative | 660.310 / 678.00 / 693.691 / 727.966 | 0.569113 | 1.705783 |
| 2048, 4.8, loss-designed | 665.921 / 683.75 / 700.089 / 735.255 | 0.482500 | 1.772555 |
| 2048, 25, conservative | 1334.635 / 1351.50 / 1365.441 / 1396.125 | 0.969086 | 1.458947 |
| 2048, 25, loss-designed | 1004.487 / 1021.75 / 1036.942 / 1069.030 | 0.752524 | 1.583653 |
| 16384, 0.6, conservative | 1530.918 / 1551.50 / 1569.077 / 1607.268 | 0.294661 | 1.924150 |
| 16384, 0.6, loss-designed | 1715.417 / 1734.50 / 1754.525 / 1794.101 | 0.238535 | 2.058895 |

Each cell below is `log₁₀ Ξ` [Å⁻³ s], lower → upper energy band. The last two join columns use only the three unflagged connected rows with tangent continuation; the first two joins include the two flagged rows. These are separate declared model choices, not a confidence interval.

| History as above | 1.287824 Å, five | 1.600 Å, five | 1.287824 Å, three | 1.600 Å, three |
|---|---:|---:|---:|---:|
| 2048, 4.8, conservative | −111.429 → −110.635 | −111.569 → −110.810 | −109.168 → −107.292 | −110.970 → −109.257 |
| 2048, 4.8, loss-designed | −113.521 → −111.979 | −113.607 → −112.117 | −111.329 → −109.295 | −112.962 → −111.076 |
| 2048, 25, conservative | −102.648 → −101.036 | −102.620 → −100.060 | −99.182 → −98.293 | −101.497 → −98.867 |
| 2048, 25, loss-designed | −107.565 → −107.071 | −108.934 → −107.332 | −105.065 → −104.304 | −106.895 → −106.084 |
| 16384, 0.6, conservative | −116.131 → −115.515 | −116.182 → −115.577 | −113.527 → −113.284 | −115.034 → −114.818 |
| 16384, 0.6, loss-designed | −117.594 → −116.937 | −117.632 → −116.989 | −115.127 → −113.672 | −116.586 → −115.238 |

At fixed join, flag treatment and energy band, the prior loss-designed episodes reduce exposure relative to the respective conservative episodes by 1.31–2.16, 4.92–7.27 and 0.39–1.60 dex in the three cases. On the other hand, the wider band at fixed launch increases exposure by 0.22–2.63 dex. Halving the time step in the selected splice-crossing control changes `log₁₀ Ξ` by −2.7×10⁻⁶ dex and sampled `log₁₀ ρ_max` by −1.2×10⁻⁶ dex. Thus neither apparent agreement of some time histories nor numerical time-step stability makes the projected-state exposure cutoff-independent.

Across all variants the RMS radius mismatch is 0.132–0.391 Å, the RMS mean-velocity mismatch is 0.0133–0.0338 Å/fs, and the closing share of exposure is 3.0–92.5%. In the 16k loss-designed, 1.287824 Å three-row upper-band example, the current at the classical closing peak is −0.0489 Å/fs while the imposed speed is −0.0676 Å/fs; the mean radius reaches its minimum 2.313 Å at 1756.523 fs, rather than the imposed 2.059 Å at 1754.525 fs. Peak contact is later still, at 1786.548 fs, when the packet's mean current is outward (+0.0280 Å/fs). This variant's added levels carry at most 3.23% of the norm, yet their isolated term is 103.5% of the integrated coherent total; destructive interference accounts for the excess over 100%. Only 3.0% of its total exposure occurs on the classical closing leg. This rules out assigning all contact to the instant of highest inward kinetic energy or to the closest radius.

Across the 48 declared cases, `log₁₀` maximum snapshot occupancy inside 1 fm spans −118.63 to −98.86 and proton probability per pair per pulse spans −109.87 to −90.53; including both ordinary DD branches adds 0.294 dex. All signed excesses over the matched undriven ground are positive in these examples, but neither the near-origin occupancy nor the total/excess probability identifies the physical material-pulse chance.

## Proton comparison and replay

At **one preparation per incident D** and equal per-D yield in the two modeled beam components, multiply the small positive total proton-pulse yield (or the separately positive excess) by `(0.85+2·0.15)/e` for emitted p/true C, and by the same 38 µA model true current for emitted p/s. The dimensionless `.25-keV PdD` model additions use the unchanged own `N65=434741.4188244709/C`, `f_cup=1`, and for `k1` the `238/38` factor. See the [fixed conversion contract](../../contact/readout_contract.md). No repetition or formation probability is measured. An assumed fraction multiplying this event incidence would multiply the one-pulse-per-D readout by that fraction. These are emitted protons and model ordinate additions, not detected counts or a comparison to a measured rate.

Base-10 ranges below span the named joins, flag treatments and energy bands, with **total** pulse probability for the one-pulse-per-incident-D conversion. Add 4.4202 dex to emitted p/s for emitted p/true C, and 0.79679 dex to `k0` for `k1`.

| Pulse history | Max snapshot P(r < 1 fm) | Proton per pair/pulse | Emitted p/s | PdD `k0` addition |
|---|---:|---:|---:|---:|
| 2048, 4.8, conservative | −112.46…−107.97 | −103.80…−99.53 | −89.37…−85.09 | −90.59…−86.31 |
| 2048, 4.8, loss-designed | −114.51…−110.09 | −105.84…−101.53 | −91.41…−87.10 | −92.63…−88.31 |
| 2048, 25, conservative | −103.53…−98.86 | −94.88…−90.53 | −80.45…−76.09 | −81.67…−77.31 |
| 2048, 25, loss-designed | −109.59…−105.10 | −101.17…−96.54 | −86.73…−82.10 | −87.95…−83.32 |
| 16384, 0.6, conservative | −117.00…−114.17 | −108.42…−105.52 | −93.98…−91.08 | −95.20…−92.30 |
| 16384, 0.6, loss-designed | −118.63…−114.35 | −109.87…−105.91 | −95.43…−91.47 | −96.65…−92.69 |

The adjacent legacy concentrated preparation's one-event-per-D readout is log₁₀ emitted p/s `−100.99` to `−98.79` over the four branches (`k0`: `−102.20` to `−100.01`). The [original stationary comparison](../README.md#conditional-quantum-contact-and-g12) has same-launch fine-step `k0` spans `−88.31…−85.26` and `−87.32…−84.23` for the lower and upper levels (emitted p/s `−87.10…−84.04` and `−86.10…−83.01`); the original coarse-step `k0` spans are `−86.31…−82.19` and `−85.26…−81.16`. Those are the different **continuously active normalized pair** hypothesis. The rows are unmodified, not stationary approximations to the pulse exposure, and not physical rates.

Replay the committed bulk case from the repository root without a lattice run:

```sh
python g11_cmns/collective/contact/contact.py --pair g11_cmns/collective/contact/bulk_pair.json --refine --output <levels file>
python g11_cmns/collective/focusing/pulse/replay.py --pair g11_cmns/collective/contact/bulk_pair.json --spectrum <levels file> --output <pulses file>
```

Both the spectrum and the `--histories` manifest bind the exact [site-pair description](../../contact/method.md#site-pair-input) by SHA-256. The default `bulk_histories.json` names the archived six traces. A new manifest has schema `g11-radial-histories-v1`, the descriptor's `pair_sha256`, and `histories`: each row names `case`, `launch`, an `archived_trace` NPZ relative to `--trace-dir`, its `trace_key`, `scored_samples`, and `first_original_exit_fs` (null if absent). The existing trace array layout supplies time in fs, radial energy in eV, separation in Å, and signed closing speed in Å/fs at columns 0, 2, 3, and 9, respectively. Each trace must include a complete supported closing/rebound episode. The initial packet width uses the combined outward well curvature; the entire chosen path must stay in the drive's unit-slope region. The retained energies must cover both selected energy bands.

Pulse replay checks regenerated eigenenergies against the matched spectrum to \(5\times10^{-9}\) eV, allowing the observed eigensolver portability difference. The independent \(2\times10^{-5}\)-dex contact-gauge check is unchanged. The committed [`results.json`](results.json) is a historical result, not overwritten by the replay commands.
