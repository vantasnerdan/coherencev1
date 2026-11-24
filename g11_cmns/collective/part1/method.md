# Part 1 — concentration is not uniform translation

**At exactly the same initial energy, changing the spatial pattern changes relative D–D motion.** On the actual 32-D PdD G10 graph, the 0.037968-Å preparation gives:

| Cold preparation, 2-ps history | Closest pair separation (Å) | Peak relative kinetic energy anywhere in that history (eV) | Smallest positive-energy participation ratio | Longest interval with ≥25% on one site (fs) |
|---|---:|---:|---:|---:|
| Uniform Γ control | 2.887824 | <1.1×10⁻³¹ | 32 | 0 |
| Finite-q, (2π/8.168,0,0) Å⁻¹ | 2.860309 | 0.0006055 | 16 | 0 |
| Localized seed | 2.804337 | 0.0027067 | 1 | 412.5 |
| Localized seed, energy matched to **twice** the Γ amplitude | 2.735134 | 0.013712 | 1 | 1065.5 |

The first three share **0.00484409 eV total** excess potential energy, not a per-site energy or a COM amplitude assigned to each pair. The last has its own higher, exactly matched Γ reference. These are deliberately prepared patterns, not the archived pumped state. Crucially, all three nonuniform cold closest approaches occur at the initial turning point: **their contemporaneous relative momentum/kinetic energy is zero**. The later kinetic peaks must not be combined with those minimum distances to manufacture a rate. Every realization's joint closest-distance/time/pair/momentum/energy record is retained in `realizations.csv`; `history.csv` retains the time-dependent joint observable.

The full Hamiltonian is the saved `g10_quantum_site/model.py` K tensors, F₀, positions and reciprocal folded J, with its cubic anisotropic quartic. The Γ identity is exact for identical sites: equal displacements and velocities cancel in pair differences. A 10%-transverse finite-q seed, parametric h=.08 drive, actual saved **41.62-meV clock** and model loss remain explicit; its peak relative energy is 5.45×10⁻⁶ eV in the archived-amplitude-sized small-graph preparation, with no quarter-energy concentration over 2 ps. This is not a general absence of parametric growth. The bare G10 Γ clock is **20.460 meV**, distinct from the saved zero-point-Hartree clock; no zero-point correction was silently stacked onto G10 dynamics.

**Saved pump geometry matters.** The archived complex x-channel summary has q=(−0.048078,−0.048078,+0.048078) Å⁻¹, wavelength **7.5453 nm**, carrier weight .14234 and Γ coherent fraction **.02886**. The vector phasor spatial RMS is .037968 Å, not a common Γ peak or time RMS. Its full per-site complex field/global phase is absent from the delivered packet. For an ideal *single carrier* at that actual q, the NN relative peak amplitude would be ≤.007443 Å and kinetic peak 1.156×10⁻⁵ eV; that illustration is not a reconstructed trajectory. The 8.168-Å G10 cell cannot represent this long Bloch carrier. Neither scalar amplitude nor a changed q is substituted for the missing field.

### Thermal, one-ion and one captured D

Eight finite classical canonical chains at 77/300 K use the actual potential and **physical Maxwell momenta**, not PIMD fictitious velocities. Classical mean kinetic energy is 3kT/2 (.00995/.03878 eV); G10's quantum equilibrium values are .03669/.05161 eV, so these references are not interchangeable. A thermal 300-K realization reaches 1.975873 Å at 1139.25 fs with contemporaneous ΔK=.043345 eV. Only about 49% of full-model site-times are within the .24-Å radial diagnostic; the .875 partial model, retaining its measured K/F₀/J and transferred closures, has about .13% within .03 Å. Their smaller-distance extrema are quartic extrapolations, not measured-well predictions.

The fixed-profile 0.25-keV PdD ledger deposits **294.1453 electronic +90.7856 host recoil +.06916 direct-D recoil J/true C**, and leaves **115 J/C unresolved at the 100-eV-per-D transport cutoff**. At 2.5 keV the components are 4251.1947+631.3248+2.48054 J/C. `stopping.json` retains separate D⁺/D₂⁺ central-track depth/time histories and all four foil splits. Direct D recoil is not all stopping. A specified .1% share of direct recoil reaches the modeled cell, with electronic/host-to-D efficiency 0 or .001 and impulsive versus four-pulse timing. The actual impulse work, including the thermal cross term, is recorded. Nearest-site assignment and pulse timing are **model sensitivities**, not measured recoil sequences; leftover deposition is not injected twice. Excess-energy participation subtracts a simultaneous same-initial-state/no-impulse reference. At efficiency .001 the impulsive/split cases reach 1.967852/1.995464 Å; their excess concentration and same-site retention are in the realization records.

The capture scenario uses **one actual deleted site, original index 5**, in the .875 bunched cell: 28 existing D plus one tagged stopped D. Its added well/bonds are an explicit x1 closure, not a relaxed capture calculation or double occupation. Existing atoms begin in the partial thermal state; the new tracer starts at the vacant register with a 300-K Maxwell momentum, **not 100 eV**. Its closest tagged pair is 1.931404 Å at 467.5 fs, ΔK=.044253 eV, p_rel=3.03894 eV·fs/Å. The conditional equilibrated tagged comparison has minimum 2.027537 Å. Post-capture/equilibrium tagged 2-Å counts are **2/0** over eight 2-ps realizations, all beyond radius support. This is a tagged-event scenario, not new bulk D density.

### Reaction consumer and adjacent data

The **2-Å observation shell is not a nuclear/scattering matching surface**. Cold cases have 0 sampled crossings over 992 pair·ps each (one-crossing resolution 1/992 pair⁻¹ps⁻¹, not a universal zero). The 300-K full/partial cases have 2/2; η=.001 impulsive/split tracks have 1/1, all outside measured radius support. No compact-contact or valid incident-flux map is supplied. The existing target-conditioned [0.1,1.6]-Å single-ray BO samples do **not** overlap the supported lattice chart: even the generous .24-Å two-site radius leaves R≥2.407824 Å, a **.807824-Å gap**, and the preparations differ. See `radial_matching.md` for actual coordinate/force matching and the no-double-counting rule for incoming connected cuts.

Thus **physical emitted rates of these trajectories are unresolved**, not zero. `rates.csv` puts this unknown beside every relevant G12 A_k and independent SI1 count proxy, with conditional numbers for Γ_p=A_pρ_contact: A_p=5.81220×10⁷ Å³/s per normalized allowed-spin s-wave pair. For one event's excess integrated summed contact Ξ [Å⁻³s], Y_p/C=4.17185×10²⁶ Ξ (before a conditional capture fraction); for one stationary active domain's summed contact C_d [Å⁻³], Γ_p/s=5.81220×10⁷ C_d and Y_p/C=Γ_p/I_true. Active domains are a separate specified multiplier. These are **conversion coefficients, not rates inferred from the sampled shell crossings**.

At .25 keV, the adjacent k0 comparison is:

| Foil | Conditional stationary ΔR per C_d, one active domain | G12 A₀ at this setting | Independent SI1 detected/s |
|---|---:|---:|---:|
| Pd | 5.0844×10⁶ C_d | .01656 | 4.0833×10⁻⁵ |
| PdD | 3.5182×10⁶ C_d | .05052 | 5.8333×10⁻⁵ |
| Ti | 5.6091×10⁴ C_d | .00195 | 2.9167×10⁻⁵ |
| TiD | 5.5779×10⁴ C_d | .00458 | 1.0889×10⁻⁴ |

Each ΔR uses **its own model neutron 6.5-keV denominator**, primary p/n join1 and retained focusing c(E), not detector acceptance. k1 multiplies model readout by I_cup(6.5)/I_cup(E); the observed ordinate is unchanged. SI1 acceptance and absolute observed reference are unknown, so counts are not subtracted from emitted rates. Ti(D) columns are unit conversion only: no matched Ti(D) force graph has been supplied. No profile/loading or fit is changed.

**Three largest omissions:** (1) matched joint inner/contact dynamics, including angular state and reaction normalization, not supplied by wide one-site cuts alone; (2) wide/connected force support, especially partial-cell transverse/other-site and scalar-bond closures; (3) measured stopping-to-D/capture preparation and active driven-domain population/loss, plus a matched Ti(D) graph.

**Evidence and reproduction:** `results.json`, `consumer_proof.json`, the complete realization/history tables, split stopping and conditional readout table. Run `python g11_cmns/collective/part1/motion.py` in the installed scientific environment. `--refine` reuses the unaffected cold cases and regenerates the paired excess, tracer and crossing observables.

The beam discriminators remain the plateau p/n ratio, proton line shape, loading/temperature changes and an applied vibration with measured spatial phase; they distinguish preparation, relative motion and nuclear branching rather than treating larger COM amplitude as a reaction count.
