# Oxide-cap beam-fit mechanisms: Pd first

## Decision

Within the commissioned **40 eV true-shift, fixed-profile CSDA model**, an exposed D-rich layer cannot reproduce the 310 or 600 eV uniform-fit reference curves. At 2.5 keV lab its optimistic infinitesimal-sheet limit is **4.136 times** the normalized uniform-40 yield; the 310/600 references require **9.899/61.580 times**. No finite thickness or density ratio reaches either full curve in that exposed frozen-stopping family. The sheet's best whole-curve projection is **239.85 eV**, with **29.06% maximum fractional shape residual**: a fitted number is not a successful shape match.

A hypothetical, constant actual beam-energy gain can imitate the shift very well, but needs about **−650 V at the collision**, not a demonstrated oxide voltage. Beam-driven positive charging has the opposite sign, and a conservative cap field cannot leave a lossless energy gain at a collision in grounded buried Pd. The cited film evidence does not establish the required negative supply, irradiated electrical isolation or in-cap D inventory. **Neither mechanism, alone or together, is established as the explanation of the measured oxide contrast.** This is not an all-profile/all-transport no-go. The Pd physical-reach condition for extending this calculation to Zr with O is unmet, so **no Zr yield calculation was performed**. In particular, no exact ZrO2 cap or 8/20-keV yield ratio is invented.

Authority is Axiom's Fleet `Coherencev1/g7_micro/part-a`, exact `entry_3185`; the existing part-a selection is untouched. This private calculation does not change the retained static ~38–39 eV results, start a new native calculation, or refit the original experimental data.

## Observable, equations and units

Kasagi measured D(d,p)T relative yields from **2.5–10 keV D lab = 1.25–5 keV CM**, bracketing each setting with 10-keV measurements on the **same prepared foil**. The compared Pd and PdO/Pd/Au foils were different prepared targets. Reported fits are 310 ±20(stat) ±50(sys) and 600 ±20 ±75 eV, respectively—not measured local potentials or depth profiles. Source links and exact caveats are in [material_source_notes.md](material_source_notes.md).

We use the existing `g4_beam.forward.model.reaction_sigma(..., ch='p', reaction='eq7')`: Bosch–Hale D(d,p)T S(E), with the commissioned Eq7 screened kernel

\[
\sigma(E_c,U)=\frac{S_{BH}(E_c)}{\sqrt{E_c(E_c+U)}}\exp[-31.397/\sqrt{E_c+U}].
\]

The default consumer's **later common-resonance** option is deliberately not inserted into a comparison whose source baseline is Bosch–Hale. No optional resonance sensitivity is needed for the physical-reach decision. The beam here is D+, monoenergetic, with no unsupported molecular fraction or energy-spread parameter from the unrelated default consumer beam.

For a frozen composition and D profile,

\[
Y(E_0)/C=e^{-1}\int n_D(z)\sigma(E_{lab}(z)/2,U)10^{-31}\,dz,
\quad dE_{lab}/dz=-S_{stop}/1000,
\quad R(E_0)=Y(E_0)/Y(10\,\mathrm{keV}).
\]

Here sigma is in mb, stopping in eV/m, lab energy in keV and densities in m⁻³. `curves.csv` stores physical reactions per **true incident coulomb**, its own-foil 10-keV ratio, the fitted uniform ratio, and residual. Every monitor sees exactly the same profile and, for offset cases, the same offset law. There is **no free amplitude, stitch, detector correction, or unshifted substitute monitor**.

Fit U alone by equal-weight least squares in **log R** at 2.5, 3.0, …, 9.5 keV lab; 10 keV is the exact reference and has identically zero residual. `log_rms` and `max_fractional_shape_residual=max|R/Rfit−1|` report shape mismatch. These declared synthetic-grid weights are not the authors' statistical weights. A lower-bound U=0 fit, notably the buried-cap curves, is a failed shape representation, not evidence of vanishing physical screening.

The energy-integral specialization uses the **same consumer Pd/D stopping** and `dx=1000 dE/S`; cumulative path inversion is the constant-composition specialization of `csda_nodes`. Its 0.1-keV-lab cutoff matches that consumer. The final script checks homogeneous Pd, a D-rich exposed Pd step, an exposed in-cap step with the exact density ratio, a buried step behind PdO, and actual `forward_yields` reference normalization. No project-wide suite is run.

**Source convention is not erased:** Kasagi used its earlier adjusted reaction curve and Andersen–Ziegler stopping and describes approximately 50× enhancement over its own dotted bare curve at the lowest energy. Our Eq7/Bosch–Hale+ICRU/ZBL reference is a commissioned comparison, not their raw-data refit. `results.json` separately gives the exposed sheet, uniform-40 and uniform-600 ratios against this model's **own U=0 curve**, so the 61.580× number above is not mislabeled as the measured 50×/bare value.

## Depth-profile map and reach

`map.csv` spans t = 0.05, 0.1, 0.3, 1, 3, 10, 30, 100 nm and r = 1, 2, 5, 10, 30, 60, 100, 1000, 100000 (in-cap t≤30 nm). **r always means nD_layer/nD_bulk**, including in PdO; it is not a ratio of D/host stoichiometries. Bulk nD is a declared **0.01 D/Pd × 6.80×10²⁸ Pd/m³**, not a measured near-surface value. This is within the scale of the source's yield-equivalent densities, but those densities assume the very uniform profile being questioned.

Four geometry/model columns distinguish:

- `Pd_frozen`: exposed 0–t region with enhanced reaction density, bulk stopping held fixed. This permits the positive-mixture analytical envelope.
- `Pd_D_feedback`: same profile, with D stopping included at its local density.
- `in_30nm_PdO`: exposed 0–t **D-bearing oxide**, then D-free oxide to 30 nm, then uniform Pd bulk. It requires accessible D *inside the cap*, not merely underneath it.
- `under_30nm_PdO`: D-free oxide from 0–30 nm, then the rich t-nm Pd layer, then uniform bulk. A buried interface is not relabeled as an exposed sheet.

PdO has explicit **Pd+O**, not an effective Z. Its formula-unit density comes from the sourced crystalline cell (a=0.3043 nm, c=0.5336 nm, two PdO units), about 4.05×10²⁸/m³. Pd and D retain the consumer ICRU/ZBL convention; O uses the repository's estimated LS electronic plus ZBL nuclear stopping, summed atom by atom. Oxygen's total stopping is also multiplied by 0.5/2 as a **matched numerical sensitivity**, not an experimental confidence interval or stoichiometry change. The actual irradiated film's phase, porosity and stopping are unmeasured here.

For frozen stopping, write the yield as B+(r−1)L_t. At each incident energy,

\[
R=\frac{B_{10}R_B+(r-1)L_{t,10}R_{L_t}}{B_{10}+(r-1)L_{t,10}},
\]

an exact convex mixture with **reference-yield weights**, not arbitrary density weights. The exposed-sheet limit is the favorable boundary. Equivalently, for the shared-depth kernel ratio h(z)=sigma(E_low(z)/2)/sigma(E_ref(z)/2), decreasing per-species functions s_j(E) d(log sigma)/dE imply decreasing h for any nonnegative Bragg mixture. The script records negative finite-grid increments for Pd, D and O over 0.1–10 keV lab. This supplies a model-scoped sufficient-condition diagnostic, **not an interval-certified universal theorem** about real transport. Even the lower-energy point alone exceeds the favorable boundary by factors **2.394 (310)** and **14.890 (600)**; `t_max_nm` and `r_required` are therefore **null/unreachable**, not zero or an unsearched large value. Letting t→0 and r→∞ does not remove that shape boundary.

The finite map includes local D-stopping feedback. On its hypothetical ≤1 D/Pd benchmark subset the best exposed-Pd grid point is **t=1 nm, r=100**, Uapp=**205.64 eV**, with **25.36%** maximum shape residual (frozen stopping: 208.15 eV, 25.42%). These are maxima over the declared grid, not optimized physical capacities. Ratios above the single-site benchmark remain explicitly labeled sensitivities. One octahedral site/Pd and separate electrochemical D/Pd≈0.7 measurements are useful scales, **not measurements of this oxide interface or universal bounds**. The oxide has no independently established D-per-formula capacity; the ≤1 label there is only a hypothetical comparison. An uptake barrier does not prove D accumulates beneath it.

For the explicit in-cap geometry, its hypothetical ≤1 D/formula-unit grid subset reaches **209.48 eV** at **t=1 nm, r=30** (D/formula≈0.504), with **24.35%** maximum shape residual. Halving/doubling O stopping moves that point to **204.71/215.23 eV**, not 600 eV. The deliberately excessive-inventory in-cap grid maximum is **230.24 eV**. None of these is measured D retention. The same-model exposed sheet is **6.217× its U=0 reference**, while uniform-600 is **92.562× that reference**; these are distinct from the source's approximate 50× comparison.

**Direct use of the saved action converter.** At the two measured lab endpoints, use the nominal chord-Hermite gains in `g7_micro/pd_oxide/action/results.json`. For an exposed sheet, the D(d,p)T ratio with either computed barrier is `sigma_BH(1.25,0) / sigma_BH(5,0) * exp(G_site(1.25) − G_site(5))`, with energies in keV CM. Relative to the same uniform-40 reference at 2.5 keV lab, this gives **4.355×** for the bare saved site's barrier and **4.170×** for the coated saved site's barrier, instead of 4.136× for exactly 40 eV. This check uses the computed gain at both endpoints, not the 0.1 Å lowering; it does not move the measured-contrast decision or change the assumed 40 eV in the map.

### Paths versus reaction depths

`depths.csv` gives every 0.5-keV setting. In the present straight-ahead model, Pd path to 0.1 keV is **82.2→198.5 nm** from 2.5→10 keV, while its mean reaction depth is **3.86→15.05 nm**. These quantities are not interchangeable, nor are CSDA path and a scattered/projected stopping range.

Kasagi's stated ~5/15-nm effective reaction depths and statement of >70% of 10-keV yield within 5 nm are retained as **author reports**, not stopping ranges or extra data used to fit this calculation. Our homogeneous model instead gives **25.8%** of its 10-keV yield in the first 5 nm (71.8% at 2.5 keV). Thus those source statements cannot jointly be silently treated as a verified depth law for this consumer.

`oxide_paths.csv` reports the actual 30-nm cap traversal, exit energy and full CSDA path for all energies and all three O-stopping scales. In this straight-ahead sensitivity every setting reaches the buried interface, but with substantial loss; that is **not measured transmission** or a claim that all real ions reach it. If an actual path length exceeds the listed stopping path, access is blocked and an under-cap enhancement cannot contribute; scattering/roughness cannot be replaced by assigning the 5/15-nm reaction depths to the cap. The nominal buried curves suppress low-energy yield and do not fit the enhanced reference. No experimental PdO range is inferred from a Pd-metal surrogate.

## Constant energy offset and electrical feasibility

For D+, a collision at potential phi relative to the beam reference has delta_E_lab=−e phi; the first scale for extra 300 eV CM is therefore about **−600 V**. Solving the actual normalized kernel, including its shifted monitor and prefactor/stopping, gives:

| True → apparent U (eV) | Actual lab gain (eV) | Collision phi (V) | \|field\| across 30 nm (GV/m) | Max shape residual |
|---|---:|---:|---:|---:|
| 40 → 340 | 650.57 | −650.57 | 21.69 | 0.801% |
| 310 → 610 | 667.38 | −667.38 | 22.25 | 0.172% |
| 310 → 600 | 644.76 | −644.76 | 21.49 | 0.172% |
| 40 → 600 | 1231.74 | −1231.74 | 41.06 | 0.749% |

21.69 GV/m is **216.9 MV/cm**, not 21.69 MV/cm. This is a required DC stress, **not a measured PdO breakdown threshold**. A conservative electric field over a cap that returns to grounded metallic Pd gives zero net lossless collision-energy offset in that metal: entry acceleration is undone before the buried collision. Charging can instead alter stopping, transmission, implantation, reflection and where reactions occur; those effects require a transport/potential profile, not a fictitious retained constant energy gain. In-film reactions can sample a nonzero local potential but require the unmeasured in-cap D and voltage distributions. A bias on the entire holder/foil, or a misreferenced accelerator/Faraday assembly, is a different, testable global energy offset. A neutral is not accelerated by the static field; where it neutralizes and its previous energy history would matter.

**Sign and current are decisive.** D+ arrival and escaping secondary electrons drive positive charge unless an independent negative electron supply/other current overwhelms them; metal neutralization/contact leakage can instead prevent accumulation. Positive potential retards D+. The source's same nominal lab-power condition implies I∝1/E, so a fixed Ohmic leakage path would normally produce a larger *positive* potential at low energy, not a constant negative offset. Even a wrong-sign negative voltage artificially proportional to this current needs about −170 V at the monitor and −678 V at 2.5 keV to project to 340 eV; its shape mismatch is **53.7%**, unlike the constant-offset curve. The corresponding positive-voltage control suppresses the normalized lowest yield. A deliberately monitor-only positive voltage could suppress the denominator and fake enhancement, but requires a different, multi-kV monitor bias and the opposite energy ordering; its necessary lowest-point values are stored in `results.json`.

Kasagi's 4-mm spot and 60 µA **conditioning** current at 10 keV imply **only conditionally**, if the measurement monitor had that same current, I=240→60 µA and mean J=19.10→4.775 A/m². The monitor current was not independently reported; peak J and local temperature need not equal footprint means. With all of this current as Ohmic leakage, holding a 650.57-V magnitude over 30 nm requires **rho=1.14×10⁹–4.54×10⁹ ohm m**, vertical resistance **2.71–10.84 Mohm**, and sheet charge **0.1920 epsilon_r C/m²**. A leakage fraction f requires rho larger by 1/f. These are conditional magnitude requirements; they do not supply the negative current needed for the desired sign.

As a **different, engineered dielectric comparator only**, [Falidas et al., Materials 15, 8325 (2022), §3.2/Fig.9](https://pmc.ncbi.nlm.nih.gov/articles/PMC9738140/) measure room-temperature voltage-ramp breakdown **2.05–2.65 MV/cm = 0.205–0.265 GV/m** for ALD **Al2O3-doped ZrO2 3D TiN/MIM capacitors**, without a D beam; leakage is <1 µA/cm² at approximately +3 V in those stacks. The present required 216.9 MV/cm is about **82–106 times larger**. This is not a PdO breakdown threshold, not the oxide in O-implanted Zr, not an in-beam resistivity, and not a Zr fusion-mechanism extension.


The cited PdO proxy is *in-plane* conductivity roughly **4–11 S/cm** at 40°C for thermally oxidized Pd/quartz films: rho≈0.0009–0.0025 ohm m. It is not irradiated through-plane resistivity of the 170–230-K, 30-nm PdO/Pd contact; no extrapolated thin-film rho, insulating premise, or breakdown threshold is claimed. No applicable PdO dielectric strength or D+-induced secondary yield was found. ZrO2 is a different dielectric, not a substitute material match. The capacitor formula leaves epsilon_r explicit. See cited methods, geometry, temperature and phase-change caveats in the source packet.

The reported **13 eV energy error** is specific beam-calibration information, not a measurement of oxide phi. The **0.97±0.03 suppressed-cup/target current ratio** is a 3–10-keV current diagnostic across targets, not a 2.5-keV voltage probe, resolved PdO secondary yield, or trapped fraction. Its central value would correspond to an effective escaping-electron ratio ~0.031 only under a restrictive current-balance model. It must not be added again to the measured target current.

### Together, and the next build decision

An unlimited exposed D sheet is already the optimistic density-shape limit; combining it with an attractive constant potential still needs at least **−774.63 V** (magnitude **25.82 GV/m** across 30 nm) to reach the lowest-point 600-eV reference. The corresponding optimistic boundaries for 310/340 eV are −213.74/−270.94 V. These precede the rest-of-curve fit, finite D inventory, in-cap location and return to grounded metal. `results.json` records these **necessary lowest-point boundaries**, not fitted mechanisms. Conversely, a monitor-only positive potential needs **+2.369 kV for 310→600** (+2.819 kV for 40→340), with the opposite current ordering from usual constant-power Ohmic charging. A profile can bias extracted U substantially, while an actual global negative energy offset can imitate the curve closely; neither supplies a defensible PdO material-design branch with the available evidence. This keeps the Zr extension closed without excluding time-dependent profiles, non-CSDA transport, other reaction baselines or a separately measured applied bias.

**One ranked decisive measurement:** on the **same prepared foil**, keep independently measured incident D energy fixed and step beam current up/down while simultaneously measuring shallow D areal density/depth (calibrated NRA/ERDA or an equivalently depth-resolved method), local surface potential and the incident/landing-energy reference. Use a grounded bare-Pd control and reversible current ordering, control beam footprint and spot temperature, retain suppressed-cup current accounting, and bracket yield with the same shifted-or-unshifted physical monitor rather than assuming it. A yield change that tracks potential/landing energy at fixed shallow D separates charging from density redistribution; a change tracking shallow D at fixed energy/potential identifies profile evolution. No unexplained amplitude renormalization should hide either response. This measurement addresses the missing physical branch before designing a thicker oxide or transferring unrelated loading capacity.

## Files and reproduction

- `oxide_beam.py`: sole numerical implementation and narrow integrated consumer checks.
- `map.csv`: four-geometry apparent-U map, r semantics, inventory labels and shape residuals.
- `curves.csv`: selected absolute per-charge and own-monitor curves, uniform fits and residuals.
- `depths.csv`, `oxide_paths.csv`: every energy's distinct reaction/path/access quantities.
- `results.json`: decisive reach, voltage, current, combined-boundary, parity and explicit Zr-gate values.
- `material_source_notes.md`: directly consumed cited grounding packet (separate source owner).
- `README.md`: this decision record.

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python g7_oxide_beam/oxide_beam.py` from the checkout using its existing NumPy/SciPy environment. The final run completed in **16.04 s**, one thread; no Native/DFT, network, new package or broad test suite. The maximum specialization/`csda_nodes` relative discrepancy was **5.025×10⁻⁶**; actual `forward_yields` normalized-R discrepancy was **5.389×10⁻⁷**. These are numerical parity errors, not experimental/model uncertainty. The script reads the established forward model from the supplied detached base; it is not a standalone duplicated forward pipeline. Exact check errors are in `results.json`. External Fleet recording/private upload belongs to the parent; this delivery closes the scientific implementation and shared-plan scope only.
