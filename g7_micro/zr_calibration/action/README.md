# Clean versus O-adjacent ZrD₂: D–D barrier action

**Decision:** the supplied static base model does not produce a strong O-adjacent site. The O-adjacent site has *less* penetration than the clean site at every evaluated energy for all three specified contact continuations. This is a result for one ordered ε-ZrD₂ proxy with one interstitial O per 8 Zr, not experimental Zr phase identification or an action-level convergence claim. A paired large-cell qualification and full-curve electronic-grid/k-point bounds are absent. The Pd/PdO comparison remains separate.

| D–D CM (keV) | Clean Eq.7 equivalent (eV) | O-adjacent Eq.7 equivalent (eV) | ΔG (O minus clean) | exp(ΔG), O/clean penetration |
|---:|---:|---:|---:|---:|
| 0.05 (CSDA cutoff) | 38.83 | 37.85 | −0.58093 | 0.55938 |
| 0.25 | 41.51 | 39.97 | −0.15190 | 0.85907 |
| 1.25 | 43.02 | 40.69 | −0.02410 | 0.97619 |
| 4 | 44.43 | 41.89 | −0.004581 | 0.99543 |
| 6.5 | 45.27 | 42.67 | −0.002242 | 0.99776 |
| 10 (reference only) | 46.22 | 43.55 | −0.001183 | 0.99882 |

Nominal chord-Hermite rows shown. Across 0.05–6.5 keV CM and the three sampled interpolants, O−clean U_eq is **−2.63 to −0.85 eV**. Both individual shifts are energy-dependent; the stated comparison does not use a constant `U_high` or map to an occupied-site inventory. At 4 and 10 keV CM (the 8 and 20 keV lab points of [Kowalska et al.](https://doi.org/10.3390/ma18061331)), the nominal site-only ratios above imply a local penetrability double ratio of 0.9966, not their measured thick-target, current-normalized proton yield ratio. Their inferred fitted shift rose by about 300 eV under O implantation, whereas the modeled static contrast is small and negative. This disfavours the tested local environment as a stand-alone explanation; oxide/interface states, profile changes and measurement effects remain unseparated. The larger fractional change at 50 eV is at a *local* collision energy, not a comparable change in counts or in the measured plateau.

## Input and action convention

Producer private draft `453206951`, exact map asset `773449526` SHA-256 `2d4f74195bfb19f82b1f07be3de2047b19eeeddaeaef6eb8ade309f6ee87e5be`, source checkpoint `ae619f4bb2b6161c181d7d08278b9a84438275f2` and Fleet entries `2641`, `2667`, `2800`. Uses the hash-bound compact base pair and endpoint-only comparisons, not the multi-gigabyte GPW replay. Both base geometries passed the supplier's retained-site/force criteria. Identical incoming ray and cell; each quartet and contact fixes its own relaxed spectators. Both are q+1 probes with no added electron. Total BO uses the supplier's mapped Hae internal-energy convention; electronic/multipole images remain. The clean/O signed contact corrections are −41.655/−39.233 eV. The h/k checks are only at contact and 0.1 Å; the O large geometry was stopped and the clean large contact SCF did not converge. There is no numerical interval for the physical action contrast.

The reader uses `g6_pairlaw.script.material_pair.read_interaction` and the pinned forward model at `9c0fb426dfe9498e1049eb05c95daeed9b590205`. In Å and eV, write \(V=C/R+W\). Use PCHIP for the regular part on the six positive radii. Join W₀ and the first positive point at 0.1 Å with Hermite (chord slope nominal, zero slope for sensitivity) or linearly. The contact derivative is unknown; the alternative shapes are not an uncertainty enclosure. Verify that any interval supporting a turning point is monotone in the numerical potential, with the omitted tail below the query energy within the supplied domain. For both potentials integrate the positive square-root barrier from the same 4 fm inner radius, ending at the respective turning radii. Subtract under the integral in log radius to obtain `G = A_bare − A_site`. No out-of-support estimates are used. Invert the actual Eq.7 log cross-section ratio on `E+U ≤ B_G²` with the same S(E), including the flux prefactor; distinguish unsupported, suppression and missing root. The channel's S(E) cancels in this ratio; no resonance parameter or thick-target profile is refit.

## Reproduce

Supply the exact map as `inputs/delivery-map.json`, base curves under `inputs/curves/`, the two contact sidecars and `paired-comparisons.json` at `inputs/`, as specified by the map. From the pinned repository root with this directory present:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m g7_zr_action.convert inputs g7_zr_action/results.json
```

The actual compact assets reproduced 402 actions across energy/shape cases. Bare gain was zero; analytic finite-radius constant-shift gain agreed within 5.3×10⁻¹³; Eq.7 reconversion for both nuclear channels matched log gain within 1.5×10⁻¹⁴. The swapped real contact, an unsupported energy, and an Eq.7 no-root request exercised the respective refusal/status branches. These are numerical checks of this interpolation model only. No new Native calculation, producer requalification, public release or source push.
