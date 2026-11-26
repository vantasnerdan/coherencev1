# Primary-point independent-record addendum

Unpublished additive replication of accepted base `6512dbc91202bb413f6bb4d188afeafcb3c69a65`. The accepted producer, source normalization, profiles, operators and delivered outputs are unchanged. Only this directory is added. This is not another material/physics calibration or an extension of the scan.

## Model and sampling unit

Exactly 18 arms: nine spatial/timing beams, each at beta 0 and .35; PdD-profile 2.5 keV, 138 µA, intrinsic **energy** lifetime 50 µs, **energy** attenuation length .2 m, free-edge/default reflectivity 1, eta .01. The retained modal band is 5–250 kHz, source/pre-bin step 50 ns, exposure **100 µs per independent record**. Each arm has the same fixed 32 independent record IDs 0–31. Physical ions, sites and time bins are not independent replication units.

The exact accepted `live.simulate` is called through an additive RNG/input adapter, not a second wave implementation. Every record starts from zero, burns for the accepted five maximum modal energy lifetimes (4,704 steps, 235.2 µs), then measures 2,000 steps. Field masks and feedback use only the preceding field. The causal-prior rms, stop/count conventions, collective-bin energy cap, species and two-recoil parent semantics remain the accepted ones. Common record totals are Poisson at I/e, split separately over burn and measurement; all 18 matched arms share each record's charge. Each parent still has an actual sampled low/middle/high mark and sign, not a mean mark or Gaussian surrogate.

Space 0/1/2 means uniform / 3-mm Gaussian / Gaussian plus the accepted two hot spots. Time 0/1/2 means steady / .7-amplitude 100-kHz ripple / ripple plus Markov ON/OFF burst (gain 9, mean ON 2 µs, OFF 18 µs), normalized over each charge-conditioned interval. These are hypotheses, not measured beam profiles.

`inputs.npz` is the canonical, generated **numeric input**, produced once by accepted `coarse`/`live` routines in the source runtime. It freezes the actual 8×8 synthetic columnar-grain assignment's retained eigenvectors, eigenfrequencies, attenuation map, signed normalized source columns, reference covariance, spatial/receiver/region projectors and all 32 beam-intensity paths. The scalar Legendre basis degree is (5,5,2), with accepted mass normalization. `source_ledger.json` binds exact source/input hashes and runtime. Replaying uses this file; regenerating it is a distinct numerical-input epoch and is not required for replay.

## Portable addressed draws, version 1

`portable.c` vendors the NumPy **2.1.2** BTPE/inversion binomial and PTRS/product Poisson implementations, with their log-gamma helper and state structure; see `NUMPY_LICENSE.txt` and upstream SHA in the source ledger. They are compiled as standalone C99 with `-O2 -ffp-contract=off -fno-fast-math`; the installed NumPy sampling implementation is never called. Multinomial allocation is sequential conditional binomial, with the remaining count assigned to the last category. No parent-count thinning or stochastic-mark mean replacement is used.

Each scalar transform starts a fresh rejection substream. The ASCII prefix is `g13-pass2-primary-portable-v1|` followed by stage, `|`, and decimal record ID; its SHA-256 first eight bytes, little-endian, are XORed with `step*262144+cell` and passed through the 64-bit SplitMix finalizer. Each uniform word adds `0x9e3779b97f4a7c15`, then uses that same finalizer; the upper 53 bits divided by 2^53 give U in [0,1). The exact finalizer/constants and arithmetic are in the two small sources.

Stages are `total`, `time`, `site`, `species`, `sign`, `low`, `high`, `beam`. Total/time step indexes are the burn/measurement interval; site/species/sign/marks use the full physical bin index. Flattening is C order: site (x-major 6×6 grid), species (D+, correlated D2+ parent), sign (+,−), then mark (low,middle,high). A multinomial component i gets subkey `mix(key XOR mix(i))`. The accepted low-then-conditional-high split is preserved. Rejection in one component cannot consume another component's stream. Spatial/timing/feedback labels are deliberately absent from matched source/mark keys: totals, meaningful beam paths and draws are common random numbers, while distinct record IDs remain independent pseudorandom records. This is comparison coupling, not additional experimental replication.

All binomial probabilities, including sequential multinomial conditionals, are rounded to the nearest 2^-40 (ties upward). The absolute probability change is at most 2^-41 per conditional Bernoulli parameter. This is a declared fine numerical approximation to the accepted continuous probability law, not exact real-arithmetic identity, nor a bound on feedback-sensitive observables. It avoids sensitivity to ordinary last-digit differences without asserting that every possible rounding boundary is protected.

`portability.json` records the **observed** stage: the actual delivered space2/time2/beta.35 batch, records 0–7 including burn, reproduces every record field, actual low/high integer draw digest, quantized mark-probability digest and energy ledger exactly in NumPy 2.1.2/SciPy 1.14.1 versus NumPy 2.5.5/SciPy 1.18.1 on one platform. A separate one-record replay also matches signed-field bytes; a consistent alternating modal-sign change leaves those class counts and signed field identical. Record 0's integer draws/counts also match between solo and eight-record execution, but its quantized probability digest differs with batch shape: draw identity is not identity of every floating field/probability digit. This does not assert an unobserved 18-arm second-runtime replay or universal bitwise C-libm/BLAS/platform independence.

The original reported event-table divergence remains authoritative but its initial cause is **not established**. A separate fixed-n,p probe matched raw words and transformed draws in both runtimes. A local conditional-p perturbation demonstrably changed accepted rejection consumption and later draws. Thus equal parent totals or a stable raw generator alone never proved equal live marks. The addressed design removes that hidden cross-cell stream-consumption effect; frozen numeric inputs address the independent operator/basis boundary.

## Tables, units and uncertainty

- `records.csv`: all 576 arm/record observations; physical totals, actual mark totals/digests, A and C signed/union parent counts, fractions and impact rates; B occupancy, complete-event count, censoring and rate.
- `summary.csv`: all 18-arm means, **sample SD across 32 records**, n and pointwise 95% Student-t intervals (31 df), for A/C fractions/rates and B occupancy/rate.
- `contrasts.csv`: within-record differences from **uniform/steady at the same beta**, matched beta .35−0 for every beam, and feedback-by-beam difference-of-differences. Intervals use the 32 paired differences, not independent-arm or parent-level errors.
- `signed.csv`, `signed_summary.csv`: existing central 2.5-mm-square, 0–50-µm signed epsilon_xx/u_z time-occupancy and complete-event summaries at ±3 and ±5 causal-prior rms. These are standalone signed wave excursions, not impact-mark probabilities.
- `energy.json`, `ensemble_verification.json`: actual accepted source/propagation ledger, exact arm/record coverage and parent/mark conservation.

A is impact into a pre-existing disturbance (.25 quiet canonical rms in either component); C_any is impact into the union of ±3 causal-prior-rms strain/displacement excursions. Signed C components overlap; do not add them to C_any. A and C overlap. Their rates have physical-parent seconds as denominator. B is off-impact distinct-sector inward-flux plus constructive interference; its occupancy denominator is receiver time, not parents. Its accepted `initiation_rate_s` counts **complete connected excursions**, excluding boundary-censored events. Never add B's receiver rate to an impact rate.

Intervals are approximate, pointwise finite-record intervals, not simultaneous multiple-comparison bounds or model/discretization uncertainty. Full paths are the statistical units: within-record temporal correlation and all matched-arm correlation are preserved by the record-level contrasts. Exposure is short relative to some field/beam variation; the fixed burn/rms initialization is retained, not proved to reach a stationary law. A zero record count means no event observed in that record, not a zero inferred population rate. Constant finite samples retain mean and SD but leave inferential intervals empty; no Poisson zero-count upper bound is imposed on correlated wave excursions.

The practical precision criterion fixed before seeing the ensemble is **95% paired C_any-fraction halfwidth ≤0.5 percentage point** for each beam-minus-uniform/steady comparison. N=32 is fixed, with no optional-stopping extension. Unmet precision or intervals crossing zero remain unresolved; selection of a largest observed arm is descriptive, not an unbiased optimized effect.

## Results and next decision

Entries below are **mean ± sample SD across 32 independent records**, in percent. The columns use their different denominators as defined above. Full signed-C, event-rate and uncertainty rows are in the CSVs.

| Space | Time | beta | A impacts % | B receiving time % | C_any impacts % |
|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 92.144 ± 4.201 | 34.333 ± 14.095 | 0.055 ± 0.142 |
| 0 | 1 | 0 | 92.006 ± 4.758 | 34.962 ± 13.713 | 0.105 ± 0.309 |
| 0 | 2 | 0 | 90.257 ± 5.954 | 33.398 ± 13.678 | 0.158 ± 0.395 |
| 1 | 0 | 0 | 92.204 ± 3.640 | 33.566 ± 11.951 | 0.097 ± 0.348 |
| 1 | 1 | 0 | 92.049 ± 3.570 | 34.288 ± 12.779 | 0.160 ± 0.479 |
| 1 | 2 | 0 | 90.194 ± 6.438 | 33.933 ± 12.043 | 0.363 ± 0.979 |
| 2 | 0 | 0 | 91.771 ± 4.616 | 33.159 ± 10.784 | 0.061 ± 0.152 |
| 2 | 1 | 0 | 91.339 ± 4.124 | 34.089 ± 12.097 | 0.217 ± 0.832 |
| 2 | 2 | 0 | 90.230 ± 4.851 | 34.005 ± 12.390 | 0.210 ± 0.656 |
| 0 | 0 | .35 | 96.408 ± 2.803 | 35.530 ± 10.621 | 1.512 ± 2.047 |
| 0 | 1 | .35 | 95.893 ± 3.844 | 33.705 ± 11.225 | 1.354 ± 1.837 |
| 0 | 2 | .35 | 95.078 ± 4.515 | 31.998 ± 13.097 | 1.487 ± 2.852 |
| 1 | 0 | .35 | 96.270 ± 2.528 | 33.562 ± 8.906 | 1.382 ± 1.605 |
| 1 | 1 | .35 | 96.082 ± 3.062 | 33.698 ± 11.328 | 1.222 ± 1.593 |
| 1 | 2 | .35 | 95.112 ± 4.388 | 32.163 ± 11.625 | 1.700 ± 3.059 |
| 2 | 0 | .35 | 95.775 ± 2.899 | 33.606 ± 9.043 | 1.248 ± 1.954 |
| 2 | 1 | .35 | 94.936 ± 3.536 | 32.752 ± 9.219 | 0.924 ± 1.535 |
| 2 | 2 | .35 | 94.949 ± 3.062 | 33.534 ± 13.489 | 1.054 ± 1.995 |

Selected **paired** contrasts, in percentage points, with pointwise 95% intervals:

| Contrast | A | B occupancy | C_any |
|---|---:|---:|---:|
| Space2/time2 − uniform/steady, beta .35 | −1.459 [−2.405,−.513] | −1.995 [−6.633,2.642] | −.458 [−1.268,.352] |
| Uniform/steady, beta .35 − 0 | 4.265 [3.260,5.269] | 1.197 [−2.311,4.704] | 1.457 [.735,2.178] |
| Space2/time2, beta .35 − 0 | 4.719 [3.336,6.102] | −.470 [−3.839,2.898] | .844 [.304,1.384] |
| Space2/time2 feedback-by-beam interaction | .455 [−.781,1.691] | −1.667 [−5.017,1.683] | −.613 [−1.349,.123] |

All nine paired feedback effects on C_any are positive at the **pointwise** interval level. In contrast, **none of the 16 beam-minus-uniform/steady intervals excludes zero for either B occupancy or C_any**. A clumping advantage for B/C is therefore unresolved, not established absent. The positive/negative C components are separately tabulated under symmetric thresholds; their finite-record imbalance is not a damage or sign-bias law.

The predeclared C_any halfwidth target is met in **10/16** beam contrasts; the six unresolved precision targets are beta .35 at space/time 0/2, 1/1, 1/2, 2/0, 2/1 and 2/2. The widest halfwidth is **1.119 percentage points**. No sample-count extension or optional stopping was used. All 576 record/arm observations conserve marks/parents and matched charge; maximum relative source/propagation energy residual is **3.77e−16**.

**Decision:** the conditional feedback intervention is distinguishable in these finite records; a B/C beam-shape winner is not. Do not select a favorable single realization or rank clumping arms as resolved. Any further precision allocation is an owner decision outside this completed fixed ensemble, not an automatically started wider scan.

## Reproduction

Extract at repository root on the exact accepted base. Python needs NumPy and SciPy; a C99 compiler is needed for production/replay. Use the supplied `inputs.npz` rather than regenerating the operator:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python g13_rogue/pass2/replication/replicate.py run --batch 0
python g13_rogue/pass2/replication/replicate.py run --batch 1
python g13_rogue/pass2/replication/replicate.py run --batch 2
python g13_rogue/pass2/replication/replicate.py run --batch 3
python g13_rogue/pass2/replication/replicate.py summarize
python g13_rogue/pass2/replication/replicate.py portability --output g13_rogue/pass2/replication/representative.json
python g13_rogue/pass2/replication/replicate.py package
```

`prepare` transparently regenerates canonical numeric inputs by calling the accepted routines, but it is not the default replay path. `summarize` reads saved arrays/tables without importing the producer or Numba; it consolidates and removes its own batch shards. `package` creates a deterministic root-relative archive with only the additive filenames and a per-file byte/SHA-256 manifest. No accepted source/output is included as an archive member. The delivery receipt records authenticated download/manifest verification; the release remains an unpublished draft. No commit, push, merge, wider scan or claim promotion is part of this addendum.
