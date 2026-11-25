# Covered PdD(111): provisional D28–D26 steady contact

The **1.600 Å join** gives a useful comparison **within the borrowed bulk-connected model**. The **1.288 Å join** needs too much unmeasured one-body continuation to quote a rate there. Neither result is a surface measurement or a prediction of foil emission.

## The reported join

Logs are base 10. Contact is the spatial \(l=0\) density at zero separation for one normalized pair (Å⁻³), not a dimensionless probability. Emitted p/s assumes one continuously active pair. G12 \(k_0\) is the unchanged **PdD 0.25-keV** comparison, not detected counts. The slab values use the near-curvature exterior described below and all five bulk-connected points; the three-retained-point branch changes the 1.600 Å slab logs by only +0.026/+0.036 at ground/300 K. This small difference does not bound the error from borrowing the bulk curve.

| Preparation | Bulk log contact | Slab log contact | Slab log emitted p/s | Bulk → slab log G12 \(k_0\) |
|---|---:|---:|---:|---:|
| Ground | −111.71 | −104.30 | −96.54 | −105.17 → −97.76 |
| 300 K | −106.75 | −101.03 | −93.26 | −100.21 → −94.48 |

At this join, changing only the rest separation **2.88782 → 2.43704 Å** raises the ground/300-K contact by **+10.70/+8.46 decades**; then replacing the bulk wells with the two slab wells at that same separation subtracts **3.29/2.73 decades**. The effects are ordered counterfactuals at fixed connected term and join, not unique additive material contributions. The inward slab wells are **stiffer** than the bulk [110] well on these cuts: at 0.36 Å per site the costs are 0.1932 eV for D28, 0.0755 eV for D26, and 0.0561 eV per site in bulk. D26 is softer than D28, not softer than this bulk reference.

## What limits the short join

| Join \(R_j\) | Displacement per site \(u\) | Beyond signed normal cut | Spread across the three one-body continuations, ground / 300 K |
|---|---:|---:|---:|
| 1.600 Å | 0.41852 Å | 0.05851 Å | 0.31 / 0.43 decades |
| 1.288 Å | 0.57461 Å | 0.21460 Å | 3.41 / 5.05 decades |

The short-join thermal spread exceeds the 3.66-decade nominal well effect. Its ground spread is also comparable with that effect. No short-join slab rate is selected. With the two flagged connected points excluded, the contrast is similar (short-join 3.46 / 5.22 decades, long-join 0.31 / 0.43). The 1.600 Å result is still extrapolated and conditional—not an error bar or a bound on the unknown surface curve. Matched one-site controls along the actual pair cut would replace the unmeasured wall as well as support the connected subtraction.

## Registration, model and replay

- Common relaxed covered slab at [13d3e76](https://github.com/vantasnerdan/coherencev1/commit/13d3e76), authenticated source asset **792605823**, delivered `site-wells.json` SHA-256 `3ad2ed3b2cdf6d54bcc14408bbb3ed7e06bd945a13e508c41eac270e6162c89d`. The reported D28→D26 vector is **0.361°** from the negative surface normal. Normal-cut energies were evaluated at the axial normal component and slopes multiplied by that projection; the unmeasured transverse displacement at the last point is 0.00227 Å. The origin uses the measured residual force, not an imposed stationary zero. For the exact accepted signed source bytes use the delivered archive; the version on main has only path strings redacted.
- The registered reduction uses \(u=(R_0-R)/2\) and \(U_{\mathrm{lat}}=V_{28}(u)+V_{26}(u)+W_{\mathrm{bulk}}(u)\). The signed one-body Hermite curves pass through 0, ±0.12 and ±0.36 Å with their measured slopes; positive quadratic exteriors match the last value and slope. The three continuations use terminal curvature (floored at one fourth of the near secant), the symmetric ±0.12 Å slope secant, or four times that secant. They are sensitivity examples, not measured exteriors. The bulk connected curve is a stand-in, not reinterpreted as sampled slab energy. At each original join, shift the same old bulk inner interaction in energy to match values; neither its geometry nor its force match is established for the slab.
- Ground and canonical 300-K contact use the accepted normalized radial solver and both original join/flag choices. Dropping the highest 1 eV of the retained 3-eV spectrum changes the 300-K log contact by at most \(1.5\times10^{-10}\) dex; the contact-weighted peak is below 1.06 eV. This tests the thermal sum, not the exterior model or a new spatial-grid convergence claim.

From this directory: `python surface_steady.py PATH_TO_SITE_WELLS_JSON --output NEW_RESULTS_JSON`. [Saved results](surface_steady_results.json) include the exact three descriptors, the shorter-distance control and raw sensitivity rows, not extra slab rate recommendations. There is no vacancy pair or bound vacancy-centered state in the supplied cell. No pulses, lattice runs, or new force calculations were performed.
