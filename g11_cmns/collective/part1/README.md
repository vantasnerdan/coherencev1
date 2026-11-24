# Part 1 — relative motion and concentration

**Equal mode energy distinguishes concentration from uniform translation.** Actual G10 full-PdD K/F₀/J/positions, 32 D, 2-ps histories:

| Cold preparation | Closest R (Å) | Later peak relative kinetic energy (eV) | Minimum energy participation ratio | Longest ≥25%-on-one-site interval (fs) |
|---|---:|---:|---:|---:|
| Uniform Γ control | 2.887824 | <1.1×10⁻³¹ | 32 | 0 |
| Finite-q, (2π/8.168,0,0) Å⁻¹ | 2.860309 | .0006055 | 16 | 0 |
| Localized seed | 2.804337 | .0027067 | 1 | 412.5 |
| Localized, twice-amplitude energy sensitivity | 2.735134 | .013712 | 1 | 1065.5 |

The first three share **.00484409 eV total**, matched to Γ peak displacement .037968 Å; the stronger row has its own matched Γ reference. **Every nonuniform cold minimum occurs at the initial turning point, with relative momentum zero.** The later kinetic peak and minimum R cannot be combined into a reaction event. All cold rows stay inside the .24-Å radial diagnostic; the cubic-direction/bond closures remain model assumptions. Joint per-realization R/time/pair/μΔv/μ|Δv|²/2 and histories are retained, not averaged independent extrema.

Γ is a control, **not the saved pump**: its actual carrier q=(−.048078,−.048078,+.048078) Å⁻¹ has wavelength7.5453 nm, carrier share .14234 and Γ fraction .02886. Full complex phase-field bytes are absent; the 8.168-Å force graph cannot represent that Bloch wavelength. The saved .037968-Å spatial phasor RMS and **41.62-meV model clock** are preserved separately from G10's bare Γ clock20.460 meV. A finite-q-seeded h=.08 driven case gives peak ΔK=5.45×10⁻⁶ eV over2 ps; this finite case does not settle parametric growth in other preparations.

**Thermal/one-ion/captured cases.** Eight finite classical canonical chains use physical Maxwell momenta, not PIMD velocities; classical 77/300-K kinetic energies .00995/.03878 eV are distinct from G10 quantum .03669/.05161 eV. At300 K the full/partial models reach1.975873/1.986494 Å, outside their .24/.03-Å support. Only49%/.13% of site-times lie inside those diagnostics. The PdD .25-keV stopping split is **294.1453 electronic +90.7856 host recoil +.06916 direct-D recoil J/true C**, with115 J/C left unresolved at the transport cutoff. Explicit .001 transfer/share and impulse/timing sensitivities preserve that ledger; paired excess energy subtracts a simultaneous no-impulse reference. No100-eV cutoff energy is assigned to vibrations.

One conditional capture puts **one tagged D at actual vacant index5**, not an occupied site or new bulk density. Its x1-transferred added well/bonds are not a relaxed capture calculation. Postcapture minimum tagged R=**1.931404 Å at467.5 fs**, contemporaneous ΔK=.044253 eV, p_rel=3.03894 eV·fs/Å; equilibrated tagged comparison2.027537 Å. Tagged2-Å counts are2/0 over eight2-ps realizations, all beyond support.

**Rate boundary.** The2-Å observation shell is neither nuclear contact nor a scattering match. Cold rows have0 sampled crossings over992 pair·ps each: resolution1/992 pair⁻¹ps⁻¹, not universal zero. All thermal/track/capture crossings are beyond support. Existing [0.1,1.6]-Å single-ray BO samples have a **.807824-Å gap** from the generous supported lattice domain and a different preparation. `radial_matching.md` states what connected pair data can repair, including subtracting the already-present +J_ABu²; one-body wide wells do not close the inner/contact bridge.

**Physical trajectory fusion rates remain unknown.** Conditional normalized-contact conversion is Γ_p=5.81220×10⁷ρ_contact/s; event-integrated excess contact Ξ [Å⁻³s] gives Y_p/C=4.17185×10²⁶Ξ, while stationary summed contact C_d [Å⁻³] gives Γ_p/s=5.81220×10⁷C_d **per explicitly active domain**, divided by I_true only for conditional/C. These are not rates inferred from shell crossings. Example .25-keV k0 columns:

| Foil | Conditional stationary ΔR | Adjacent G12 A₀ | SI1 detected/s |
|---|---:|---:|---:|
| Pd | 5.0844×10⁶C_d | .01656 | 4.0833×10⁻⁵ |
| PdD | 3.5182×10⁶C_d | .05052 | 5.8333×10⁻⁵ |
| Ti | 5.6091×10⁴C_d | .00195 | 2.9167×10⁻⁵ |
| TiD | 5.5779×10⁴C_d | .00458 | 1.0889×10⁻⁴ |

`rates.csv` supplies all adjacent k0/k1 emitted/C/time coefficients, each own modeled n6.5 denominator, join1 and focusing once. SI1 acceptance is unknown; counts are not calibrated/subtracted as emitted rates. Ti(D) has no matched motion potential. Profiles/loading and the G12 baseline are unchanged.

**Three omissions:** matched inner/contact state and flux; wide/connected force support and partial closures; measured transfer/capture/domain preparation and loss, including Ti(D).

**Evidence:** `results.json`, `realizations.csv`, `history.csv`, `stopping.json`, `rates.csv`, `consumer_proof.json`; derivations and replay in `method.md` and `radial_matching.md`. Beam discriminators: plateau p/n, proton line shape, loading, temperature and applied vibration with measured spatial phase.
