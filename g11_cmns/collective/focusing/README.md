# Can a distributed bulk PdD launch focus onto one D–D pair?

**Finite-model result, unpublished.** The target is the actual D25→D63 [110] bond; every launch uses the same total N e. Wave/direct launches start with zero target velocity and zero displacement. Reverse recovery and localized collisions are shown, but are not silently promoted to collective preparation.

For the same distributed direct launch at N=2048, 25 meV/D, the two resolutions give **0.83–0.98 eV**, both above the wave's 0.157 eV. This two-result spread is **not an error bar** or converged peak. Only 1.61–1.91% of the 51.2-eV budget reaches the target; no common incoming state or saturation is established.

## Same-budget comparison

At **dt=.5 fs**, each cell is **gain K_rel/e; contiguous half-peak hold (fs)**. † is initially preconcentrated/not headline-eligible; ‡ marks an edge-censored hold. The final direct anchor has the visible .25-fs sensitivity below. Absolute energies, target radial component, separation, crossings and same-trajectory whole-cell maximum are in [comparison.csv](comparison.csv).

| D count | meV/D | Waves | Reverse | Collision | Direct |
|---:|---:|---:|---:|---:|---:|
| 32 | 0.6 | 7.11; 29.5 | 32; 14.0†‡ | 7.77; 30.5 | 14.9; 27.5 |
| 32 | 4.8 | 5.2; 27.5 | 32; 12.0†‡ | 7.53; 28.5 | 17.3; 24.0 |
| 32 | 25 | 9.01; 22.5 | 32; 10.0†‡ | 9.32; 22.5 | 18; 20.0 |
| 256 | 0.6 | 15.9; 27.5 | 256; 12.0†‡ | 43.2; 26.5† | 35.6; 25.5 |
| 256 | 4.8 | 13.6; 25.5 | 256; 9.0†‡ | 70.1; 21.0† | 31.9; 22.5 |
| 256 | 25 | 5.68; 22.5 | 256; 4.5†‡ | 0.0282; 1.5†‡ | 29.6; 18.5 |
| 2048 | 0.6 | 12.8; 28.0 | 2.05e+03; 9.0†‡ | 569; 20.5† | 221; 24.0 |
| 2048 | 4.8 | 6.98; 24.5 | 2.05e+03; 3.5†‡ | 0.0759; 1.0†‡ | 99.5; 19.0 |
| 2048 | 25 | 6.29; 22.5 | 2.05e+03; 1.5†‡ | 0.00206; 0.0†‡ | 39.2; 18.5 |

Direct launch positive-site-energy participation is 827/2048, target preload 0%. Collision participation is 8: its eight-site packet gets more energetic as N grows, not more extensive. Reverse states have positive-site participation 2–9.44 and target site-energy fraction 8.7–1e+02%; exact recovery is a reachability control wherever that launch is poorly spread.

## What limits the result?

3/36 selected trajectories hit the .80-Å site or 1.287824-Å nearest-bond boundary; every quoted peak is strictly before any exit, and censored holds remain labelled. The same supported-prefix rule is used for every family. The model's extent is a box of assumed support, not measured angular coverage.

At N=2048/.6 meV the same-launch harmonic peak is 0.472 eV before its site-boundary exit at 703 fs, versus 0.132 eV for the nonlinear trajectory. At 25 meV the harmonic trajectory exits at 150 fs, long before the nonlinear peak: its 0.0235-eV supported prefix is **clipping, not a same-horizon nonlinear advantage**. The one fixed 4136-fs amplitude-loss sensitivity changes the coarse high-energy direct peak to 0.316 eV (-67.7%). No damping sweep was used.

The same launch at .25 fs peaks at **0.826 eV** at 1563 fs, versus 0.98 eV at 2367.5 fs at .5 fs: -15.8%, so **the largest-energy peak is not numerically converged**. Holds are 18.5/18.75 fs; energy drift is 0.00027/6.5e-05. Both exact same-launch rows and their pre-peak exits are retained. The updated local wall costs 1–3.46 times the old harmonic wall at the actual most-displaced site's sampled directions above .03 Å (maximum radius 0.7858 Å). The .25-fs trajectory reaches the finite horizon without a support exit; neither a near-end peak nor a first boundary crossing establishes saturation. The matched high-energy harmonic trajectory exits early, so its smaller supported peak is clipping, not a same-horizon nonlinear advantage.

| Same launch, no reoptimization | .5 fs | .25 fs |
|---|---:|---:|
| Peak Krel = radial K (eV) | 0.9802612980 | 0.8257418508 |
| Peak time / half-peak hold (fs) | 2367.5 / 18.5 | 1563 / 18.75 |
| Signed inward speed (Å/fs) | +0.13707240 (inward) | -0.12580597 (receding) |
| First support exit | none by 2400 fs | none by 2400 fs |

This is late-pulse/peak-selection sensitivity, **not changed clipping**: the fine trace near the coarse peak (±30 fs) reaches only 0.188 eV. The two opposite signs must not be treated as the same incoming state.

## Energy floors first; extrapolations second

Perfect-concentration floors for 100 eV / 1 keV are N ≥ **166,667 / 1,666,667** at .6 meV/D, **20,834 / 208,334** at 4.8 meV/D and **4,000 / 40,000** at 25 meV/D. At N=2048 the per-D floors are **48.83 / 488.28 meV**, before inefficiency. Neither target is simulated.

**Scaling is finite and preparation-dependent, not an asymptotic limit.** Direct gain slopes are between 32→256 and 256→2048 D at .5 fs. The 25-meV anchor is not converged; the table exposes both timestep-based arithmetic assumptions rather than claiming a fitted high-energy size law.

| anchor meV/D | two coarse finite-size exponents | 100 eV / 1 keV N if efficiency stays fixed | eV/D at N=2048 if gain stays fixed |
|---:|---:|---:|---:|
| 0.6 | 0.42, 0.877 | 1.55e+06 / 1.55e+07 | 0.453 / 4.53 |
| 4.8 | 0.293, 0.547 | 4.29e+05 / 4.29e+06 | 1.01 / 10.1 |
| 25 | 0.24, 0.135 | 2.09e+05 / 2.09e+06 (.5 fs); 2.48e+05 / 2.48e+06 (.25 fs) | 2.55 / 25.5 (.5 fs); 3.03 / 30.3 (.25 fs) |

Those N and e values are **assumption-only extrapolations**, not simulated reachable states; measured efficiency is already energy-dependent, and the largest-energy step sensitivity limits its numerical anchor. `results.json` also retains the finite-size power-law scenarios and their incompatible slopes/residuals. No three-point flattening or growth establishes saturation or its absence. The fixed 2400-fs horizon, nonlinear wall/defocusing and model exits remain distinct.

## Conditional quantum contact and G12

All 36 target peaks have a [retained per-level lookup](results.json), with both joins and both original flagged-row treatments. The illustrative assignment is **ΔE=classical K_rel above each quantum model’s own ground**, not equality to an eigenenergy or an inferred quantum state. Adjacent l=0 eigenlevels are separately reported; their endpoint span is not a rigorous wavepacket bound.

| Base-10 logs, stationary-level assignment only | .5-fs peak: lower / upper endpoints | .25-fs peak: lower / upper endpoints |
|---|---:|---:|
| contact (Å⁻³) | [-92.86,-88.73] / [-91.81,-87.70] | [-94.86,-91.81] / [-93.87,-90.78] |
| emitted p/s | [-85.09,-80.97] / [-84.04,-79.94] | [-87.10,-84.04] / [-86.10,-83.01] |
| emitted p/true C | [-80.67,-76.55] / [-79.62,-75.52] | [-82.68,-79.62] / [-81.68,-78.59] |
| PdD .25-keV k0 | [-86.31,-82.19] / [-85.26,-81.16] | [-88.31,-85.26] / [-87.32,-84.23] |
| PdD .25-keV k1 | [-85.51,-81.39] / [-84.46,-80.36] | [-87.52,-84.46] / [-86.53,-83.43] |

These are **one continuously active normalized-pair** units, not detector counts, a foil rate or a peak×hold event yield. Preparation probabilities, radial populations and coherence are unknown. The retained excitation support ends at 7.906–7.969 eV; no 100-eV/1-keV contact or rate is extrapolated.

**Three largest omissions:** (1) the off-ray local/connected angular force surface is not measured; (2) dynamic Pd/electrons, finite-temperature disorder and physical damping are absent except the stated one sensitivity; (3) generating the prescribed phases and converting a classical focus into occupied quantum radial states is not supplied.

The closest connected alternate-H difference is −0.0047422856 eV, with signed slope difference +0.0413149205 eV/Å; this finite-cutoff comparison is not a PAW error bound.

![Finite-size gains and selected target histories](focusing.png)

Equations, selection, support, resource controls and portable commands: [method.md](method.md). Exact staged capsule bytes and whitelist: `capsule.json`, `ARCHIVE_FILES.txt`, `SHA256SUMS`; prior-input assumptions: `input_manifest.json`. No prior contact or .80-well calculation was rerun.
