# Pd and Pd/PdO: D–D Eq. 7 equivalents on the delivered base frames

**Decision:** the saved coated frame does not show the large positive static shift sought in this comparison. It gives lower penetration than the saved bare frame at every sampled energy and contact shape. This is a property of the two supplied frames. It is not a resolved sign for equilibrated Pd/PdO: the coated frame is partly relaxed, with no geometry-error bound; no paired h/k check or cell/slab/vacuum control was delivered.

| D–D CM (keV) | Bare Pd Ueq (eV) | Pd/PdO Ueq (eV) | Coated − bare (eV) |
|---:|---:|---:|---:|
| 0.05 (CSDA cutoff) | 39.00 | 37.37 | −1.63 |
| 0.25 | 43.49 | 39.94 | −3.55 |
| 0.50 | 44.77 | 40.44 | −4.33 |
| 1.00 | 45.74 | 40.96 | −4.78 |
| 1.25 | 46.03 | 41.14 | −4.88 |
| 2.50 | 46.96 | 41.84 | −5.12 |
| 4.00 | 47.70 | 42.45 | −5.25 |
| 5.00 | 48.10 | 42.79 | −5.32 |
| 6.50 | 48.64 | 43.25 | −5.39 |
| 10.00 (separate extension) | 49.67 | 44.15 | −5.52 |

Nominal chord-Hermite results. Across the sampled 0.05–6.5 keV CM window, including zero-slope Hermite and linear W, coated-minus-bare Ueq is **−5.45 to −1.63 eV**. Maximum spread due to these three shape choices at a common energy is 0.49 eV. Neither range bounds numerical or geometry error. Different observables: the local nominal coated/bare penetration ratio is 0.379 at 0.05 keV and 0.706 at 0.25 keV, despite a few-eV shift on the very steep cross section; neither is a thick-target yield ratio.

The [Kasagi et al. comparison](https://doi.org/10.1143/JPSJ.71.2881) fits about 310 eV for Pd and 600 eV for PdO/Pd/Au from 2.5–10 keV D lab thick-target proton yields normalized at 10 keV lab. On those 1.25–5 keV CM endpoints our nominal site-only coated/bare penetration is 0.9510 and 0.9932. This comparison is opposite in direction and much smaller in equivalent shift; it does not identify the sampled sites, reconstruct the prepared target, or compare calibrated counts. There is no constant `U_high` or occupied-site inference.

## Source, conventions and replay

Scientific inputs in `g7_micro/pd_oxide/capsule/` at `a31c95ebb9e722fce5c0937d43a2cafab7172133` were compared byte-for-byte to the producer's draft `453374618`, asset `774098833` and its bound direct inputs. The same six actual radial samples occur in each original and common-support projection. Each signed contact is bound to the original curve. The adapter uses only bound metadata to populate the earlier reader's text provenance fields; arrays and energy gauge are unchanged. The two compared frames have the same cell, ray and ray precommit. The bare 15Pd+D frame passes the producer's force criterion (free fmax 0.0360 eV/Å); the coated 19Pd+4O+D frame adds a PdO sheet and is partly relaxed (free fmax 0.4217 eV/Å). Geometry error is unbounded for both. No nuclear-charge or site-retention inference follows from the sheet addition.

The unmodified [Zr action kernel](../../zr_calibration/action/convert.py) (commit `e996395214680cd038a7cf9b208d0a5b088e51e2`) computes the finite-4-fm WKB gain `G = A_bare_Coulomb − A_site`, using PCHIP in regular W at positive R and the pinned short-range contact continuations. Only the increasing branch of the actual Eq. 7 cross-section ratio is inverted, at fixed S(E); q+1 supplies no added electron. Negative coated-minus-bare ΔG denotes weaker penetration than the bare *slab*, not suppression relative to unscreened Coulomb. At all evaluated points both sites have positive G relative to Coulomb.

From the checkout root, with NumPy/SciPy:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m g7_micro.pd_oxide.action.convert \
  g7_micro/pd_oxide/capsule g7_micro/pd_oxide/action/results.json
```

The actual inputs ran 402 paired energy/shape cases. Bare gain was zero, analytic finite-radius constant-shift log gain agreed within 5.3×10⁻¹³, and Eq. 7 log-ratio roundtrip with both nuclear channels agreed within 1.5×10⁻¹⁴. The wrong real sidecar was refused. No Native calculation, full producer replay, equilibrium claim, publication or source push.
