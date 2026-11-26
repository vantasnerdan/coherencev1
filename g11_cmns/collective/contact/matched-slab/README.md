# Covered PdD(111): matched D26/D28 steady contact

Unpublished replacement for the provisional borrowed-bulk-W calculation. Base is exactly `6d795cc90fbf62cfca46d286447d1e46262831e4`. Historical source, driver, output and prose are unchanged except the single forward pointer in `g11_cmns/collective/contact/covered-slab.md`.

## Stationary contact result

All density columns are **log₁₀[ρ/(Å⁻³)]**; the last two columns are log₁₀ density ratios. All five primary native knots, including flags, are retained.

| Rj (Å) | Preparation | Matched slab | Accepted bulk PdD | Accepted TiD₂ | Slab/Pd | Slab/Ti |
|---:|:---|---:|---:|---:|---:|---:|
| 1.287824094 | Cold | −102.565325 | −111.714162 | −104.728448 | +9.148837 | +2.163123 |
| 1.287824094 | 77 K | −102.563788 | −111.674167 | −104.728441 | +9.110379 | +2.164652 |
| 1.287824094 | 300 K | −98.298891 | −106.720420 | −103.187592 | +8.421529 | +4.888701 |
| 1.600000000 | Cold | −104.993319 | −111.714082 | −104.588161 | +6.720763 | −0.405158 |
| 1.600000000 | 77 K | −104.992062 | −111.674465 | −104.588154 | +6.682403 | −0.403908 |
| 1.600000000 | 300 K | −101.864843 | −106.752558 | −103.011949 | +4.887714 | +1.147106 |

**The supported long-join result:** the matched slab outer model is +6.72 decades above bulk PdD cold and +4.89 at 300 K, under the common borrowed inner shape. It is **.405 decades below TiD₂ cold but 1.147 above at 300 K**; there is no temperature-independent slab/Ti ordering. The short-join numbers are the explicitly requested **conditional quadratic-A/B + tangent-W continuation figures**. Moving to the long join lowers slab contact by 2.428 decades cold and 3.566 at 300 K: the short answer must not be promoted to a measured material contact.

### Distinct sensitivities, not error bars

| Variant | Join (Å) | Cold log density | 77 K log density | 300 K log density |
|:---|---:|---:|---:|---:|
| Coarse 2×2×1, matched A/B/W | 1.287824094 | −103.339569 | −103.338381 | −99.033306 |
| Coarse 2×2×1, matched A/B/W | 1.600000000 | −105.608592 | −105.607605 | −102.389857 |
| Mesh3 terminal-inward-curvature continuation | 1.287824094 | −102.437988 | −102.436438 | −98.080102 |

The coarse comparison changes **all A/B/W duties and their own origin forces**, not W alone; it retains the same outward normal proxy and projection rule. Relative to primary, coarse shifts cold/77/300 K by −.774245/−.774593/−.734415 decades at the short join and −.615273/−.615543/−.525013 at the long join.

The one targeted continuation contrast changes only inward one-body curvature **beyond u=.48**, from default K_A/K_B=1.353295/3.124932 to terminal native Hermite K_A/K_B=.790096/.229434 eV/Å². Native points and derivatives, origin, outward proxy and W tangent are unchanged. It shifts short-join cold/77/300-K contact by +.127337/+.127351/+.218789 decades. The long-join Hamiltonian is exactly unaffected because that extension lies below its join and is replaced by the common inner shape. This contrast does not reverse the short-join ordering, but is neither an uncertainty bound nor evidence for the unmeasured connected tangent.

**Old → new, long join:** the unchanged provisional near-exterior calculation gave cold/300-K logs −104.303054/−101.028372; matched replacement gives −104.993319/−101.864843, changes −.690265/−.836471 decades. This replaces inward A/B **and** the bulk W stand-in, not an isolated-W counterfactual.

Primary additive inner shifts are .416148063/.632674106 eV at short/long join. Lattice-minus-inner **radial derivative** mismatches U′_lat−V′_G6 are −1.372470733/+4.006043155 eV/Å; corresponding radial-force mismatches are their negatives. The finite non-force-matched joins are deliberate, not hidden in the contact normalization.

### Numerical support of the new observable

The actual consumer reproduces all native mesh3/mesh2 total energies and dU/du to maximum discrepancy **5.56×10⁻¹⁶** in eV or eV/Å. New slab calculations compare h=.001→.0005 Å, Rmax=9→10 Å, Ecut=3→4 eV, forbidden match height .15→.25 eV and ODE relative tolerance 2×10⁻⁹→3×10⁻¹⁰. Largest final log-contact change across all five model/join cases and all temperatures is **2.98×10⁻⁶ decades**.

At primary 300 K, 99.9% of the **contact numerator**, not merely the partition function, is below 1.180/1.386 eV at short/long join. The peak contributing states are .864/1.126 eV, far above the .03424-eV ground. Fractions above 2 eV are 1.02×10⁻¹²/5.78×10⁻¹⁰; contact fractions in 3–4 eV have logs −26.07/−22.71. Across all variants the largest top-1-eV contact fraction is below 4.0×10⁻²³. At 77 K the contact support is near the ground (.104 eV contains 99.9% for primary). These are finite support/convergence observations, not rigorous infinite-spectrum or electronic-error bounds; zero reported cutoff changes mean floating-point resolution, not exact zero tails.

The independent explicit-basis/native proof gives maximum 5.56×10⁻¹⁶ eV or eV/Å discrepancy, with zero outward-knot discrepancy. Its separate 28-level long-join solve (h≈.000750 Å, Rmax=9.5 Å, Ecut=2 eV, u=R exp(z) origin propagation) reproduces the named cold/77/300-K results within **3.06×10⁻⁶ decades**; the retained primary spectrum quantifies the omitted >2-eV contact contribution. No independent new electronic calculation is claimed.

## Registered model

A is **D26, subsurface, moving toward +axis**; B is **D28, adsorbed, moving toward −axis**. They are not equivalent wells or exchanged labels. The shared frozen reference has R₀=2.4370398315432267 Å and axis A→B=(0.00545223290348, 0.00315290979902, 0.99998016596139). Each ion moves u=(R₀−R)/2. Native A, B and AB controls have the same geometry, atom order, electronic settings and energy zero within each mesh. The recorded DFT energies are force-consistent electronic free energies.

U_lat(R)=V_A(u)+V_B(u)+W(u), with V_A=E_A−E₀, V_B=E_B−E₀, W=E_AB−E_A−E_B+E₀. The independent path derivative is dU/du=−F_AB,A·axis+F_AB,B·axis; dW/du subtracts the corresponding single-site derivatives. The source `connected_dE_du_eV_A` is **dW/du**, neither dW/dR=−½dW/du nor an all-atom mixed-force projection. No extra J u² or bulk W is added.

`g11_cmns/collective/contact/matched-slab/pair-m3.json` retains all five primary 3×3×1 A/B/W knots at u=.06/.12/.24/.36/.48 Å, with native force slopes. In particular A(.06)=−.000301425439 eV is preserved. The actual-axis mesh3 origin slopes are A=−.041599493679, B=+.021804708141 eV/Å: the reference is **not a radial stationary minimum**. Connected W(0)=W′(0)=0. Four shortest knots retain their original D–Pd PAW-overlap flags; no primary point is removed and no electronic-error bound is inferred.

### Signed exterior and mixed-mesh origin

The inward rays use actual-axis matched controls, not normal-cut substitutes. The outward rays are the **separate original outward normal wells**: D26 uses normal −.12/−.36 Å and D28 normal +.12/+.36 Å, each energy relative to the original accepted **mesh2** zero, not a mesh3 absolute-energy subtraction. Map q=|z|/|axis_z|, with dE/dq=−F_z dz/dq. The omitted transverse response is an assumption (last outward point has only .00227 Å transverse displacement).

Original normal-zero forces and mesh3 actual-axis zero forces differ. To satisfy the existing G11 C1 signed-ray contract without changing measured outward knots, each outward origin derivative is set to the negative of that descriptor's matched inward derivative. The first outward Hermite interval bridges that derivative to the original mesh2 .12-Å value/slope. This is a **modeled near-reference mesh/direction bridge**, not a new mesh3 outward force measurement. Coarse uses its own actual-axis mesh2 origin forces and the same outward source knots, energy reference, projection and curvature rule; hence the outward first interval is not pointwise identical to the primary. It would be inconsistent to hold both that interval and different native origin forces fixed while claiming C1.

Beyond the final knot each signed one-body ray has the existing endpoint quadratic form E_end+s_end t+K t²/2. Default K is the original signed ±.12 normal-force secant multiplied by axis_z², positive for each site. W is **zero outward** and tangent-continued after u=.48 inward. No outward W data or additional forces are invented. These tails confine the selected radial model, not a measured infinite-distance slab.

### Shared inner interaction and joins

Use the unchanged `SitePair`/`solve` in `g11_cmns/collective/contact/contact.py`. No G11 library edits or `attach_readout` conversion are made. At R≥Rj use U_lat; below use V_G6(R)+U_lat(Rj)−V_G6(Rj). The common **bulk PdD G6 occupied/vacant total** interaction includes C/R, C=14.399645468602673 eV Å. Its electronic remainder is PCHIP on the accepted .10–1.60 Å samples, with the existing smooth electronic continuation below .10 Å. This is not a slab short-range electronic measurement. Energy matches additively, with a reported derivative mismatch; there is no force fit or blend. The old inner interpolation shelf is inherited, not a new slab well.

The long join 1.600 Å has u=.418519915772 Å and is interpolated inside all five native A/B/W curves. Last native R=1.477039831543 Å. The short join 1.28782409436586 Å has u=.574607868589 Å: **.094607868589 Å per D beyond the last sample**. Its requested number is conditional on quadratic one-body plus tangent-connected continuation, not an unqualified material number.

## Normalization and evidence

Contacts are log₁₀[ρ(0)/(Å⁻³)] for one normalized pure spatial l=0 radial pair, ∫|u_R|²dR=1, ρ(0)=|u_R′(0)|²/(4π), reduced mass m_D/2. Density is not point probability. Cold is the ground eigenstate; 77/300 K are canonical **radial quantum** ensembles, summed in linear density with Boltzmann weights via log-sum-exp. Neither a full lattice/orientational ensemble nor a physical pair population is supplied. No proton/foil rate, G12 factor or Ti G12 label is inferred.

Accepted bulk PdD and TiD₂ stationary results are reused, **not rerun**, from `g10_quantum_site/tid2-site/pair/bulk-reference.json` and `g10_quantum_site/tid2-site/pair/result.json`. Their exact descriptor and common inner samples are checked; Ti rows use its retained refined observables. Numerical discretization changes on the new slab models are recorded separately from coarse electronic-mesh and short-continuation sensitivities.

`g11_cmns/collective/contact/matched-slab/result.json` retains spectra, final-log refinement changes, both joins, all three preparations, explicit contact-numerator energy support and the unchanged old/new comparator. `g11_cmns/collective/contact/matched-slab/proof.json` records an independent finite proof: direct native Cartesian-force energy/derivative reconstruction versus the actual G11 outer model, explicit Hermite basis, signed zero/C1 and outward-knot checks, and a separately assembled long-join Hamiltonian with a different regular-Coulomb ODE variable. This is finite numerical evidence, not an interval or electronic/continuum error bound.

## Source custody and replay

Accepted pair SHA-256 `e3760bfca4a7d8d78627ad248589d2252813250dd7083d2f394873c6f0ac9036`; signed wells SHA-256 `5b1d19e247bcaf546d3e38fe387b65ef7da78e4db26379688d69331075ea94f7`. Authenticated original asset **794388700**, release **458472584**, is 89555 bytes, archive SHA-256 `2734ee7808f3a3ce6f452292fc017dfdd804d4e63adeb1e4d0e3ceb0fb321049`; all **49** content sums checked, and accepted pair bytes matched. No upstream assets are repackaged here. The single 4×4×1 quartet at R=1.717039831543 Å is preserved: ΔW=−2.266024050 meV and ΔAB=+5.115599092 meV versus mesh3. It is one point, not a replacement curve, mesh sweep or error bound.

From the repository root at the exact base plus this capsule, with Python, NumPy and SciPy:

```sh
python3 g11_cmns/collective/contact/matched-slab/run.py
python3 g11_cmns/collective/contact/matched-slab/verify.py
sha256sum -c g11_cmns/collective/contact/matched-slab/SHA256SUMS
```

The first command exercises the real existing G11 consumer and writes only the three new descriptors and new result. The second writes the new finite proof without calling `solve` or the producer. Original PdD/Ti and provisional results are read only. Deterministic scientific values should be compared with the reported finite numerical tolerance across NumPy/SciPy versions; hashes bind the delivered bytes, not a promise of cross-build floating-point identity.

The three largest physical omissions remain: **slab inner electronics**, **actual angular/lattice population and dynamics**, and **unquantified electronic/PAW plus modeled outward/short-join continuation uncertainty**. The bounded result answers the matched-path stationary-model comparison only. No pulses, envelope, lattice propagation, new electronic calculation, fit to nuclear data, publication or claim promotion is part of this draft.

## Lead's note

Archive SHA-256 6755fc7b7d21cb2a74ff8bb674ac9dc8af0f6b57112ea47d3bf56052f8d3b0ef, all 9 file sums checked as delivered (`SHA256SUMS` lists this page before this note was added). The lead reran `run.py` and `verify.py` from scratch on a clean copy of the repository at commit 6d795cc (NumPy, SciPy): the three descriptors are identical and every contact figure reproduces to all printed digits; the only differences are in quantities at the 1e-75 level.
