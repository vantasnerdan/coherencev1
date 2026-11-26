# TiD₂ pair closing work and conditional steady contact

This is the **constrained, fully occupied cubic-fluorite** fingerprint supplied in `g10_quantum_site/tid2-fluorite/`, not an assignment of a titanium foil's phase. Native energies and forces are frozen-ion electronic calculations. The radial contacts are a conditional model using a **borrowed PdD inner electronic shape**, not measured Ti contact or a physical Ti/Pd reaction-rate ratio.

## Geometry, common register and closing work

Original zero-based **A=D8**, **B=D19**, zero image `[0,0,0]`, A→B **+[100]**, and **R₀=2.211551631258434 Å**. The actual saved positions obey xA=xA₀+u n, xB=xB₀−u n and **R=R₀−2u**; u is displacement **per site**. The zero and every A/B/AB control share reference, cell, original species/index order and electronic settings. `run.py` binds and reads the actual point records, rather than obtaining derivatives from the exported all-atom mixed force.

With each material's own undisturbed electronic energy E₀ as zero,

- V_A=E_A−E₀, V_B=E_B−E₀;
- W=E_AB−E_A−E_B+E₀;
- total closing work U=E_AB−E₀=V_A+V_B+W;
- dU/du=−F_AB,A·n+F_AB,B·n;
- dW/du=−(F_AB,A−F_A,A)·n+(F_AB,B−F_B,B)·n;
- d/dR=−½ d/du. The source `delta_compression_path_force_eV_A` is **−dW/du**, not dW/du and not a projection of `delta_forces_eV_A`.

All native pair rows remain primary:

| u per site (Å) | R (Å) | V_A (eV) | V_B (eV) | Total work U (eV) | Connected W (eV) | dU/du (eV/Å) | dW/du (eV/Å) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| .12 | 1.971551631 | .040003246 | .040003146 | .085782877 | .005776485 | 1.431632477 | .106817106 |
| .24 | 1.731551631 | .156971355 | .156971567 | .342237679 | .028294756 | 2.829357207 | .280087060 |
| .48 | 1.251551631 | .581935703 | .581936190 | 1.321095701 | .157223808 | 5.213873082 | .849528475 |
| .60 | 1.011551631 | .859882976 | .859883320 | 2.016737846 | .296971550 | 6.538485065 | 1.686125533 |

These are total electronic-cell increments for one registered moving pair, not per-ion energies. Positive W says closing costs more than the sum of the two frozen one-body moves; it alone does not define a free-space repulsive pair potential.

**B is not silently mirrored.** Independent B(−[100]) companion points exist at .12/.24/.48/.60 Å and are used unchanged. A's + and −[100] records are also kept separate, with signed energy derivatives. Full species-preserving periodic geometry verifies midpoint inversion exchanging D8↔D19; the independent companion energies differ from A by at most 4.87×10⁻⁷ eV, with maximum inverted full-force mismatch 9.01×10⁻⁵ eV/Å. Only B's unsampled .80-Å inward knot and B's outward ray are supplied by this explicit inversion of the corresponding A records. They are modeled symmetry completions, not independent B calculations. The separate site twofold rotation checks A's ±[100] convention without averaging them.

## Outer model and common-distance PdD comparison

`pair.json` is the existing **`g11-site-pair-v1`** input. Each site's signed ray uses native energy/slope Hermite knots, including the .60-Å A/B companions. The origin slope is set to zero by ideal crystal symmetry, not fitted to the small raw residual force. The connected interpolant has W(0)=W′(0)=0. Its reported origin curvature is only the **origin-to-.12-Å tangent/Hermite curvature**, not a measured Hessian or a full neighbour graph. There is no extra J u² or coupling matrix.

For separation beyond R₀ the explicit assumption is W=0. Each signed one-body ray continues past its final knot with that endpoint's value/slope and a positive quadratic curvature inferred from its first Hermite segment. The exterior confinement is modeled, not a measured infinite tail.

`result.json` contains both the five **native PdD** rows from `g10_quantum_site/pdd-nearest-pair/pair-cut.csv` and **common-distance interpolated** comparisons at R=1.971551631258434, 1.92782409436586, 1.73, 1.60 and 1.28782409436586 Å. The latter evaluate the actual all-row G11 outer models, within both native connected-distance ranges, and are not new electronic calculations. At each material's R₀ the work is zero; absolute electronic energies between different materials are never compared. Native Ti's 1.25155-Å row and native PdD's 1.28782-Å row are not relabeled as the same separation.

At declared common separations the all-row outer interpolation gives:

| R (Å) | Ti work (eV) | Pd work (eV) | Ti W (eV) | Pd W (eV) |
|---:|---:|---:|---:|---:|
| 1.971551631 | .085783 | .277776 | .005776 | .064586 |
| 1.927824094 | .119937 | .316031 | .008396 | .073669 |
| 1.730000000 | .344436 | .538128 | .028513 | .126925 |
| 1.600000000 | .552351 | .733226 | .050630 | .177311 |
| 1.287824094 | 1.227947 | 1.390325 | .142315 | .376384 |

Ti's total closing work is lower at each of these common distances, but most remains one-body work. The declared near-origin connected curvature is d²W/du²=**.626583529 eV/Å²**; its limited tangent/interpolation meaning is unchanged.

## What the contact observable means

The supported consumer is **`SitePair` + `solve` in `g11_cmns/collective/contact/contact.py`**, unchanged. It solves a single half-line, pure spatial l=0 D–D Hamiltonian with reduced mass m_D/2. The fixed normalization is ∫|u_R|²dR=1 and **ρ(0)=|u_R′(0)|²/(4π)**, in Å⁻³. A contact density is not a point probability. Ground, 77 K and 300 K are stationary radial quantum preparations; the thermal values use linear contact densities and Boltzmann weights via log-sum-exp, not averages of their logarithms. These are equilibrium ensembles, not a driven state or the full orientational/lattice ensemble.

No Ti inner electronic cut is supplied by the Ti source or the existing G6/G11 inputs. At each of the two accepted joins, **Rj=1.28782409436586 and 1.600 Å**, the conditional construction is

U(R)=U_lat(R) for R≥Rj;

U(R)=V_G6(R)+U_lat(Rj)−V_G6(Rj) for R<Rj.

Both materials use the same G6 PdD occupied/vacant total interaction samples and **C=14.399645468602673 eV Å** Coulomb coefficient. The existing G11 consumer interpolates the electronic remainder V_G6−C/R and uses its existing smooth-core continuation below .10 Å. The additive energy shift is recorded separately at every join; no force fit, blend or Ti screening refit is made. The recorded slope mismatch is U′_lat−V′_G6; the **radial force mismatch is its negative**. The common borrowed electronic shape, including its shallow interpolation shelf, is not new Ti evidence.

The primary uses all four Ti connected knots. The explicitly labeled sensitivity omits **only** the shortest two connected knots and uses the .24-Å endpoint tangent beyond them, leaving the same one-body rays. Its last native connected distance is 1.73155 Å, so neither join has native connected-distance overlap in that sensitivity. In this descriptor `connected.flagged` is a **qualification marker for missing short-pair PAW checks, not a claim of demonstrated overlap**. Neither treatment is an error bound; the all-row values are not discarded.

`result.json` retains the contact logs, additive shifts, derivative/force mismatch sign convention, refined contact-weighted support, final-log numerical changes and actual native-knot consumer agreement. Refinement uses the existing G11 choices h=.001→.0005 Å, Rmax=9→10 Å, Ecut=8→10 eV, barrier match .15→.25 eV and ODE rtol=2×10⁻⁹→3×10⁻¹⁰. Contact-numerator energy support and the top-1-eV contact fraction are checked separately from the partition-function tail. The cusp propagation is matched to a resolved normalized barrier amplitude; tiny eigenvector contact tails are not used directly.

The **matching PdD spectrum is reused**, not regenerated for custody: descriptor SHA-256 `7bcd74342a1ff22599eb0f98e0f8ed55aa5a722d450c40331de4a3c7c0163966`, original spectrum SHA-256 `5f40e1dc03e2127a1f29930de197c263a06e672b68338bceddb0ad0b008573f6`. `bulk-reference.json` retains only relevant stationary model/readout rows and those bindings. Its separate historical refinement retains the historical input hash; that old hash is not relabeled as the new descriptor. No Ti G12 normalization, capture weight, active-pair population, foil rate or emitted-proton conversion is inferred.

## Retained steady-contact result

All-row primary values below are **log₁₀[ρ(0)/(Å⁻³)]**, not reaction rates:

| Join (Å) | Preparation | Ti | Pd | log₁₀(Ti/Pd) |
|---:|:---|---:|---:|---:|
| 1.287824094 | ground | −104.728452 | −111.714162 | 6.985710 |
| 1.287824094 | 77 K | −104.728445 | −111.674167 | 6.945723 |
| 1.287824094 | 300 K | −103.187596 | −106.720420 | 3.532825 |
| 1.600000000 | ground | −104.588201 | −111.714082 | 7.125881 |
| 1.600000000 | 77 K | −104.588194 | −111.674465 | 7.086271 |
| 1.600000000 | 300 K | −103.011951 | −106.752558 | 3.740607 |

The Ti outer-wall model raises this conditional contact by about **7 decades cold** and **3.5–3.7 decades at 300 K**, under the common borrowed inner shape. It does not establish that physical Ti has this contact advantage.

The join is energy-continuous but deliberately not force-matched:

| Ti row treatment | Join (Å) | Inner additive shift (eV) | Radial force mismatch F_lat−F_inner (eV/Å) |
|:---|---:|---:|---:|
| all four | 1.287824094 | .810058098 | +2.502041830 |
| all four | 1.600000000 | .733546769 | −3.283483059 |
| omit short connected knots | 1.287824094 | .758179474 | +2.244621676 |
| omit short connected knots | 1.600000000 | .729634194 | −3.345039901 |

Omitting the two shortest **connected** knots raises Ti ground/77-K/300-K log contact by **.421953/.421953/.531218** decades at the short join and **.041366/.041366/.055344** at 1.600 Å. This changes the contact magnitude, not the direction of the conditional Ti/Pd comparison. The primary two-join spread is .140250 decades cold and .175644 at 300 K; these sensitivities do not bound the absent Ti inner-electronic contribution.

The single settled replay exercised the actual G11 consumer and reproduced all four native Ti total energies, connected values and signed potential derivatives to maximum absolute discrepancy **5.56×10⁻¹⁶** in their respective eV or eV/Å units. The largest final-log change under the declared numerical refinement was **4.06×10⁻⁵ decades** across both joins and both row treatments. For the all-row primary at 300 K, **99.9% of the contact numerator** lies below 1.381/1.464 eV at the short/long join; only 6.63×10⁻⁹/1.38×10⁻⁹ lies above 2 eV. The refined top-1-eV contact fraction (9–10 eV) is below 5.30×10⁻¹¹¹ for either primary join. At 77 K the ground level supplies over 99.9% of contact; the reported zero high-energy fractions are binary64 underflow, not exact proofs of zero. The largest retained-state outermost-.10-Å grid probability is below 9×10⁻⁷⁷ at the coarse outer boundary. These are finite numerical/support observations, not an electronic-model uncertainty enclosure.

## Three largest omissions

1. **Ti inner electronics are absent.** This is Ti outer-wall sensitivity under borrowed bulk PdD inner electronics, not a comparative physical material reaction rate.
2. **The actual material state and dynamics are not this one-path equilibrium model.** Frozen fully occupied cubic fluorite, pure l=0, imposed exterior confinement, no bath/relaxation/angular population, no driven preparation and no established foil phase.
3. **Pair electronic numerical uncertainty is not quantified.** The only 4→6 mesh comparison is the −.48-Å [111] single-site point and its own zero: energy-rise change .0023921746124671017 eV and maximum force-response change .010781986651845514 eV/Å. It is not a pair mesh comparison or a whole-curve bound. The two shortest pair rows have no separate PAW-overlap check; their overlap status and associated error remain unknown, not zero and not definite overlap.

## Replay and retained files

From the repository root, with NumPy and SciPy:

```sh
python g10_quantum_site/tid2-site/pair/run.py
```

This runs the actual existing G11 radial consumer, reconstructs `pair.json` from bound native records and writes `result.json`. It does not rerun DFT or PdD, edit G11 or start later work. The retained files are `run.py`, `pair.json`, `result.json`, `bulk-reference.json` and this page. Source and consumer hashes are retained in the result; no nested capsule is needed.
