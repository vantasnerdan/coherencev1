# Conditional radial readout and size extrapolation

## Classical peak → two archived eigenlevel readouts

`K_rel = μ|v_i − v_j|²/2`, with `μ=m_D/2`, includes transverse motion. The focusing table reports the radial component separately. **Assign** the supported peak `K_rel` as `ΔE` above the radial ground; in variant `j` the lookup compares it against `E_j,n − E_j,0`. This is a convention, not an inferred radial population, angular projection, scattering or beam energy. No potential energy is added and no ground-state contact is subtracted. Keep target-bond and anywhere-in-cell peaks separate.

Use the nearest retained energy below and above the assignment; a match returns the same level twice. Quote the actual endpoints and their logs, not a log-interpolated state. Every archived contact sequence examined is increasing; neither that fact nor the two endpoint rates bounds the contact of a coherent preparation or general mixture. A selected peak must lie on the supported lattice branch, strictly before any first wall exit. If a peak is outside lattice support or there is no upper retained level in a join, report **unresolved**, not zero and not the last level carried forward. Supported variants can still be shown individually when another variant is unresolved.

[Accepted contact spectra](../contact/results.json) (`ed4ec3bcd24bccd654ff42e9e6e2aabcbbfc3d65823fcb6761bcddf4147a5bef`) are used directly, with no contact recomputation. The cutoffs below are the last actual eigenenergies **relative to each own ground**, not the nominal 8-eV absolute solver cutoff.

| Join Rj (Å) | Connected treatment | Levels | First excitation (eV) | Last excitation (eV) |
|---:|---|---:|---:|---:|
| 1.287824 | Five incl. flagged | 116 | 0.041969610 | 7.908783 |
| 1.600000 | Five incl. flagged | 116 | 0.041969610 | 7.906308 |
| 1.287824 | Three + tangent | 117 | 0.041967368 | 7.951881 |
| 1.600000 | Three + tangent | 117 | 0.041967368 | 7.969219 |

The `across_join_endpoint_span` is the min/max **at each side separately** across the two joins and two connected treatments, not a confidence interval or a combined energy sensitivity. The three-row cases bridge unmeasured distances; the span is a comparison of named model choices, not an error bound. There is **no supported G12 value at 100 eV or 1 keV** even if a size law returns a finite cell size.

## What the G12 numbers mean

Use the unchanged [proton-only transform](../contact/g12_readout.py) and [fixed .25-keV CM PdD row](../contact/readout_contract.md). `Γp = A_p ρ` with `A_p = 5.81220367098898×10⁷ Å³/s`, `ρ = |u_R′(0)|²/(4π)` in Å⁻³. The archived beam-ordinate conventions give `Yp = Γp/I_true` in emitted p/true C, `ΔR_k0 = Yp/(f_cup N65)` and `ΔR_k1 = (238/38) ΔR_k0`, where model `I_true=38 µA`, `f_cup=1` and `N65=434741.4188244709`/C. These are positive model additions to the selected ordinate (negative logs do not imply a negative addition); they are not residuals, emitted neutron-inclusive rates or detected counts. Each log describes **one continuously active, normalized stationary spatial ℓ=0 radial pair**, not an excess over an already active radial ground pair. The chart conversion does not move from .25 keV to the assigned `ΔE` or to the extrapolation targets.

**No focusing event rate.** The contiguous classical half-peak duration is not an interval of known quantum contact. It cannot be multiplied by a level contact and used as an event. A separate event hypothesis requires `Ξ=∫(ρ_prepared(t)−ρ_matched_ground(t))dt`, the same radial Hamiltonian and a prepared state/contact history, plus a declared event incidence for conversion to per-second/per-charge estimates. The focusing dynamics supply none of those occupancies, rates or detector efficiencies.

Lookup discrimination, **not a measured focusing peak**: assigning `ΔE=0.05 eV` selects levels 1/2 (about 0.04197/0.08759 eV) in all four variants. Searching unshifted energies would incorrectly select 0/1. Across the four variants, lower/upper endpoint `log₁₀(p/s)` spans are [−102.2285, −100.0361] and [−100.7147, −98.4739]. The corresponding `log₁₀(ΔR_k0)` spans are [−103.4465, −101.2541] and [−101.9327, −99.6919]; `k1` adds 0.79679336 dex. All four output logs match the coefficients from the two archived [Part 1 rate rows](../contact/sources/part1-rates.csv) to ≤1.7×10⁻¹⁴ dex in this comparison. At `ΔE=7.91 eV` only the two three-row spectra have both endpoints; the across-join span is unresolved. These checks use the existing table only.

## 100 eV and 1 keV: energy floor versus size scenario

At fixed `e` per D, `Ne ≥ T` is necessary to put `T` into a single pair from the supplied excess energy, with no other energy source. Equality would concentrate all of it. The lower bounds need no fitted exponent and are not designs for a launch.

| Desired pair energy T | Minimum N at 0.6 meV/D | Minimum N at 4.8 meV/D | Minimum N at 25 meV/D |
|---:|---:|---:|---:|
| 100 eV | 166,667 | 20,834 | 4,000 |
| 1 keV | 1,666,667 | 208,334 | 40,000 |

| Desired pair energy T | Minimum e at N=32 (eV/D) | N=256 (eV/D) | N=2048 (eV/D) |
|---:|---:|---:|---:|
| 100 eV | 3.125 | 0.390625 | 0.048828125 |
| 1 keV | 31.25 | 3.90625 | 0.48828125 |

Use `conditional_size_scaling(rows)` only for a single distributed launch pattern, bond scope, `e`, loss and pre-peak/censor and selection policy across cells. Do not fit a preloaded pair, a best-pattern envelope with changing winners, or an unsupported peak. It returns two different arithmetic scenarios:

- **Constant observed efficiency** at the largest supported size `N_*`: `η_* = K_*/(N_*e)`, `N(T) = T/(eη_*)`. This is an explicitly assumed linear continuation; even a single observed `η_*` is not evidence that it persists. Integer cells must be at least the ceiling of the scenario result.
- **Finite-size fit**, only when all three comparable supported peaks grow: fit `log10 G = y0 + α(log10 N − x0)`, `G=K_rel/e`; return `log10 N(T) = x0 + [log10(T/e) − y0]/α`. Examine both adjacent slopes and the fit residual; a fitted `α` alone is not a prediction. A count below the energy floor violates the model budget, and a superlinear law cannot persist indefinitely. With a nonincreasing step or too few supported sizes, no growing power-law inversion is produced.

Three cell sizes cannot identify a saturation limit: flattening or falling efficiency may suggest loss of focusing, but an asymptote remains unresolved. Both scenarios leave the calibrated lattice branch and/or contact range at the target energies; they do not give a supported N for achieving those energies or a G12 event/foil estimate there.
