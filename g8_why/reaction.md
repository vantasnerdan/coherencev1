# Reaction — the retained threshold pole does not explain the plain foils

**Mechanism (two sentences).** A near-DD-threshold ⁴He 0⁺ amplitude changes the energy-dependent nuclear strength and interferes with the broad proton amplitude, without moving any target D. A host could affect its entrance penetration or threshold detuning, or competing electronic decay, but a common intrinsic reaction law must first survive each foil's different fixed response integral.

**What the data require.** At nominal CM energies **[0.25, 0.5, 1, 2.5, 6.5] keV**, the required *own-reference normalized excess* `Q=Robs/Rsmooth40` is:

| Foil | Q, as plotted (fixed cup-focus/p:n readout) |
|---|---|
| Pd | [5.23e16, 3.77e9, 5.35e4, 2.99, 1] |
| PdD | [2.29e17, 2.49e10, 1.63e5, 12.68, 1] |
| Ti | [2.08e17, 2.47e10, 7.98e4, 2.14, 1] |
| TiD | [4.90e17, 2.09e10, 8.40e4, 6.77, 1] |

These are **not pointwise S measurements or absolute measured yield factors**. For local-energy multiplier `m_ch(ε)`, `F_m,ch(E)=∫K_m,E,ch(ε)m_ch(ε)dε / ∫K_m,E,ch(ε)dε`, with the immutable CSDA/sheet/source response `K`; the primary fixed p/n calibration gives `Q=F_p(E)/F_n(6.5)` below the join. Absolute scale is unidentifiable from own-reference data. Different kernels prevent general cancellation of a common *energy-dependent* S; a constant common multiplier does cancel. `baseline.csv` preserves raw modeled reactions/true C separately: e.g. PdD proton baseline is [9.43e−14, 2.87e−7, .02362, 771.2, 7.54e5]. The additional-inverse-current interpretation gives Pd/PdD/Ti/TiD Q(.25)=[8.34e15, 3.66e16, 3.31e16, 7.83e16]; all five values are in the CSV. The previously reported Fig.2/SI1 normalization discrepancy remains unresolved, not a proven double correction.

**Predicted versus required.** The existing literal Eq.17 model (Γp=40 meV, phase115°, radius7 fm, ER→0; BH neutrons) predicts these **physical proton Y/C factors**, not normalized Q:

| Foil | Fp at the five CM energies |
|---|---|
| Pd | [5.861, 1.928, 1.085, .9410, .9612] |
| PdD | [5.862, 1.930, 1.086, .9411, .9601] |
| Ti | [6.410, 2.086, 1.120, .9422, .9583] |
| TiD | [6.410, 2.086, 1.120, .9422, .9584] |

Neutron factor=1; thus the physical reference proton factor is not forced to1. Total-DD factors at .25 are 3.47/3.47/3.75/3.75, **not** the earlier g7 factor5.02, which used a fitted-high background/reference. Formal one-strength fits with everything else held give equivalent Γp coefficients common27.83 eV, or per foil31.02/69.81/7.345/50.68 eV; χ²=30233 for48 shape points, versus the common positive-basis candidate399. These coefficients are failed fits, **not inferred physical widths**: the E² denominator reduction cannot independently determine Γtotal, and such large Γp would invalidate its narrow-pole premise. Alternative-current common coefficient16.80 eV still gives χ²=27634; hypothesis-dependent re-splicing also fails and is only a reporting sensitivity.

A decisive splice-free result is stronger: for **any nonnegative pole-strength coefficient and any phase** in this E⁻² term, `|I|≤(P+wB)/√3` bounds proton Y(.25)/Y(2.5) by [4.70e−14,4.19e−14,3.07e−15,3.09e−15]. The alternative-current data require [.648,.608,.251,.189]; direct prior SI1 comparisons are consistent with these order-one shapes. Retuning that pole separately by foil cannot close this 1.4e13–8.2e13 shape gap. A 12-bin positive common S/BH basis over .05–7.5 keV was also tried: χ²399/48 (alternative279), versus separately fitted sum40.9 (alternative26.8), with numerical log-sensitivity ranks only4–7/11 and low-energy coefficients as large as1e17–1e21. This is a nonunique, basis-limited inverse illustration, not an exclusion theorem for every possible S(E), nor a measured host-dependent nuclear law.

**What must be true.** A different reaction explanation needs a specified entrance/threshold response capable of the weighted factors above and the observed loading/material dependence. Smooth exit-channel shifts of100 eV change MeV p/n phase space only by order3e−5; changing a genuine ≤1 eV total width only in the off-resonant denominator changes sampled ε≥50 eV strengths by ≤1e−4, and extra IPC competition alone suppresses that term. Neither supplies the missing factors. An eV-scale pole shift can matter near that pole, but is not a demonstrated keV-beam mechanism here.

**Non-yield discriminator.** Measure matched p/n and IPC spectral branching versus energy and foil, without re-calibrating away a branch change. Dubey's thick-target IPC/p=.42±.05 at10 keV lab and .060±.012 at20 keV lab are measured ratios; fitted ΓIPC/Γp=14 is not either branch measurement and is not a UC Davis sub-keV datum. Sources and convention qualifications: [reaction_source_notes.md](reaction_source_notes.md), [Czerski2022](https://doi.org/10.1103/PhysRevC.106.L011601), [Czerskiv3](https://arxiv.org/pdf/2409.02112v3), [Dubey publishedv2](https://arxiv.org/pdf/2408.07567v2).

**Exactly three largest omissions:** (1) unresolved same-run charge/reference/detector likelihood and p/n calibration; (2) a microscopic host-dependent entrance/pole/partial-width calculation; (3) unique sub-keV S reconstruction beyond this finite basis and assumed pole convention.

**Finite proof:** shared smooth40 raw consumer parity1.11e−15 relative; positive-basis partition and retained resonance replay pass (`reaction_proof.json`). No transport, profile or established consumer was changed. **This retained reaction branch does not explain it.**
