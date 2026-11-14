# Final NiH one-point quartic: 100 only

Supersedes the coarse **100** quartic used with the harmonic-control bundle. The new-profile **111 quartic was not measured**. Axiom's entry896 allowed it only if the first actual Native unit cost stayed around100 core-min; the declared strict continuation gate was100. The first cost106.6883261666666666666666667 core-min, so production stopped after100. Do not silently advertise a fully matched-profile two-direction NiH well or cubic quartic tensor.

## Model and actual data

Original Ni32H32 geometry, a=3.732 A, source H25. Exactly one additional self-consistent calculation moves this H by+0.24 A along100, other nuclei fixed, restarting the same actual zero checkpoint as the signed harmonic pair. Electronic profile, PAWs, communicator and force-change stopping are unchanged: PW500/PBE/Ni10/H, k3 gamma-centered, nonmagnetic neutral, Fermi0.1 eV, density1e-4/force-change1e-4/energy1e-7, six genuine k-point ranks, no point-group operations. The original harmonic source epoch remains byte-identical; four extension files have a separate bound identity.

Fixed K2=1.712245063783307 eV/A² from the independently measured±0.03 A pair. Convention:

    V(u)-V(0) = K2*u²/2 + K4*u⁴/4
    F(u) = -K2*u - K4*u³

The one-point joint fit minimizes squared energy residual plus squared u-scaled force residual, the existing `analyse.fit_direction` convention. It does not fit K2 again.

| Actual endpoint quantity | Value |
|---|---:|
| u (A) |0.24|
| E(u)-E(0) (eV) |0.062244150897356576|
| F_parallel(u) (eV/A) |-0.6349102666443648|
| K4_joint (eV/A⁴) |16.16569693679319|
| K4_energy-only (eV/A⁴) |15.590631100980586|
| K4_force-only (eV/A⁴) |16.201638551531477|

Energy-only versus force-only K4 differs3.7712695%. Joint endpoint residuals, predicted minus actual: +0.00047698260685641314 eV and+0.0004968568821420716 eV/A. The u⁴ coefficient is K4/4=4.041424234198297 eV/A⁴, not K4. Full raw64-ion forces and actual positions are in `dense_Ni/k3/quartic/q100_u024.json`; readout in `dense_Ni/k3/quartic/readout-100-final-20261114T1001Z.json`.

One positive amplitude is not a multi-amplitude Taylor-convergence check or a quantum spectrum. Negative amplitudes use ideal cubic inversion symmetry; no negative large-amplitude point was calculated. K2 is the adopted finite-step force coefficient, not an independently extrapolated u→0 value. No new111 coefficient or all-direction potential is implied. The old111 coefficient1.12836572 remains coarse PW400/Gamma data if a consumer elects to retain it; any such combined two-direction model stays explicitly mixed-profile.

## Consume and reproduce

Choose a new output path:

```sh
python nih_quartic_readout.py --direction 100 --output dense_Ni/k3/quartic/readout-recomputed.json
```

The actual full reader reproduced identical parsed output in a GPW/PAW-free cold folder. A100 record placed at the111 path is refused, not used as111. The actual joint coefficient agrees with both the qualified legacy fitter and exact-rational normal equations; this algebra check is not a physical-accuracy certificate. Source, point, checkpoint and readout hashes are retained. GPW/PAW bytes and raw Native logs are excluded; exact active-criterion excerpts are provided instead.

The base harmonic/coupling scalars are unchanged from entry893; no new harmonic or coupling pair was run. `nih_quartic_profile.py` binds the extension to that source; `nih_quartic.py` reuses its actual calculator/common zero; `nih_quartic_readout.py` owns the one-point fit; `run-nih-quartic.sh` is the MPI entry. Their separate identity preserves the already accepted harmonic bundle.

## Final stop

The100 unit succeeded and stopped09:59:07Z. Exact CPUUsageNSec=6401299570000, retained as a string with InvocationID/MemoryPeak. Combined charged prior+Native2872.399736049966666666666668 core-min plus20 conservative preparation/review reserve stays inside3000. Remaining107.600263950033333333333332 was not reallocated: first cost exceeded the continuation gate,111NOTRUN, and no further Native calculations follow. Decimal costs use28-digit precision; raw integer counters are lossless. See `receipts/nih-last-quartic-cost-stop-20261114T1001Z.json`.
