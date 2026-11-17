# Clean Pd(111) sheet thermodynamics and recombination boundary

Ordered neutral PBE/Pd10/H1/PW300 adsorption phases on three-layer Pd(111), with the bottom two Pd layers clamped and top Pd/D relaxed. Deuterium masses enter the vibrations. Subsurface D, chamber C/O, larger reconstructions and nonadiabatic capture dynamics are omitted.

## Qualified sheet inputs

`input-ledger.json` binds the actual 1, 0.5 and 0.25 ML inputs. Full-mobile harmonic thermochemistry is retained at 1 ML. Lower coverages use three local-D modes, calibrated by the measured full-minus-local 1 ML correction. Configurational entropy is added once, outside those vibrational free energies.

The fitted excess free energy per site is

`g_ex(theta,T) = -(B1(T)+w(T))*theta + w(T)*theta**2`.

The chemical potential uses its derivative, not the mean binding:

`mu = kT*ln(theta/(1-theta)) - B_diff`, where `B_diff = B1+w-2*w*theta`.

At 305.033/301.405 K, `B1=0.1167729183/0.1211904141 eV` and `w=0.0707613153/0.0707124778 eV`. The maximum mean-energy fit residual is 1.995 meV/D. The quarter-ML ledger uses its tight-spectrum-center electronic energy; the original loose-endpoint convention and measured 10.85 micro-eV/D correction are retained separately.

## No-beam operating states

Intrinsic capture `s0=0.5` is an explicit assumption, not a calculated capture probability. Chamber D2 pressure is 0.03666355 Pa. The loaded case has reservoir `x_back=0.8` and holds vacancies at the actual 7200 s reference-run end; vacancies are not allowed to thermally fill forever.

| No-beam state | Sampled-route factor | θ barrierless (ML) | θ sensitivity (ML) | D2 release barrierless (m−2 s−1) | D2 release sensitivity (m−2 s−1) |
|---|---:|---:|---:|---:|---:|
| Pd gas, 305.033 K | 0.0239513157 | 0.7814494103 | 0.7814494103 | 6.5814738e19 | 1.5763496e18 |
| PdD bulk-fed, 301.405 K | 0.0228991847 | 0.9999681778 | 0.9999815638 | 7.7504425e20 | 1.7750944e19 |

Neither loaded result reaches `x+xt=0.95 D/Pd` in the top 100 nm, so the overflow-on comparison is not triggered. Beam implantation, beam vacancy production, sputter erosion and adsorbate stripping are off; overflow is omitted in this primary thermal readout. The original bulk isotherm, Darken transport, finite electrochemical supply and traps at the saved vacancy profile remain.

The original Model's public `grid_both=True` option resolves both foil faces. The common kinetic factor cancels from gas-equilibrium coverage, but not from the supplied steady state. Back-anchored shooting resolves the stiff finite-current back boundary without changing the original operator or invariant gates. `route-sensitivity-comparison.json` contains the mesh/residual evidence and the exact baseline-observable replay.

`physical-result.json` and `coverage-state-ledger.json` retain the barrierless baseline. `ridge-sensitivity.json` is a separate scientific metadata derivative, bound to its original execution. `elbow-result.json` preserves the final route's original execution identity. Metadata repacking is explicit; actual GPW/geometry bytes and scientific arrays are unchanged.

## Detailed-balance conversion

Let Gamma_exact be isotope-correct molecular impingement per Pa and Gamma_model the original transport normalization. With standard pressure P0=1 Pa:

- `R_des = Gamma_exact*s0*P0*theta**2*exp(-2*B_diff/kT)`.
- `R_ads = Gamma_exact*s0*p_D2*(1-theta)**2`.
- `J_D = 2*(R_des-R_ads)` counts D atoms, not molecules.
- `s_eff = (Gamma_exact/Gamma_model)*s0*(1-theta)**2` is the **state-local** coefficient in the original boundary.

Do not copy either table value into a universal sticking default. `coverage_boundary.py` uses the original Model's kT and molecular-flux normalization explicitly. The chemical-potential, mean/differential and molecule/atom distinctions have independent symbolic/numerical checks.

The barrierless columns use the assumed capture without an electronic-ridge multiplier. The optional `--elbow-readout` and `--elbow-readout-sha256` arguments apply `exp(-max(E_ridge,0)/kT)` to **both** adsorption and desorption, then solve the loaded steady closure again. This is the sampled route's factor upper bound, not a saddle free energy or a calculated capture probability. An incomplete route cannot supply it.

That optional run is a **sensitivity**, not a replacement for the material/state ledger. The sampled ridge is restricted-route PBE/PW300/10x5x1 data without peak k/slab qualification. `transport.py` and its current surface inputs remain unchanged.

## Reproduce the physical readout

The CLI takes explicit source, operating-point, material-ledger and saved-vacancy inputs; it does not infer a workstation layout. In the integrated `g4_beam/transport/recomb` directory, choose an available unreserved CPU (9 in this example):

```bash
ledger_sha=$(sha256sum input-ledger.json)
ledger_sha=${ledger_sha%% *}
taskset -c 9 python run_sheet_coverage.py \
  --ledger input-ledger.json --ledger-sha256 "$ledger_sha" \
  --transport-source ../transport.py \
  --operating-points coverage-operating-points.json \
  --profiles inputs/run-end-profiles.json \
  --profiles-sha256 e2e91077595de602b815ae86f45c1a675d6a553153272d0e1295b72e2a5a9943 \
  --profile-revision 98efac29bcea7effa13eb7495537eff43548858e \
  --capture-probability 0.5 --cpu 9 --output coverage-reproduction
```

The pinned transport source and its existing BCA table are required. New BCA generation is forbidden by the reader. The output directory must not exist. A wrong source digest, missing coverage, inconsistent electronic reference, unstable local mode or unconverged physical result is refused, not replaced by an estimate.

For the separate sensitivity, add `--elbow-readout elbow-result.json --elbow-readout-sha256 1c0a9c27b53f6414e64499c035d90f46ec51cd4311c3a51e5cb46b39e6b8800b` and use a new output directory such as `ridge-reproduction`.

`native-force-banks.json` indexes 57 actual compact GPW/geometry records: the 40 material force records and all 17 corrected-route frames. Payload locations are explicit in its v2 schema. All records were consumed after portability conversion. The verifier resolves paths relative to its manifest and also works from a foreign working directory.

```bash
bank_sha=$(sha256sum native-force-banks.json)
bank_sha=${bank_sha%% *}
python verify_native_banks.py --manifest native-force-banks.json \
  --sha256 "$bank_sha" --output native-bank-check.json
```


## Independent mathematical fixture

`sheet_proof.py --help` describes the bounded symbolic/numerical fixture for the derivative sheet law and frozen-vacancy closure. It takes explicit pinned transport, operating-point and saved-profile inputs, and supplies no lower-coverage material energies. Its mathematical checks do not substitute for the physical ledger or Native banks.

## Connected clean-route calculation

`elbow_route.py` constructs a continuous Cartesian gas-to-qualified-fcc path, with a nested nine-knot/17-frame mesh. `run_elbow_path.py` records actual compact GPW banks; `elbow_path_results.py` checks them against the exact saved frames and diagnoses unsampled peaks using force-consistent interpolation. A low connected path is an upper-bound argument for non-activation within the specified resolution; a positive sampled maximum is not a global activation barrier. No capture probability follows from existence of a path.

The final corrected route is complete. No further path refinement or peak k/slab control was performed.

| Route quantity | Excess above clean + gas (meV) |
|---|---:|
| Relaxed central knot 8 | 24.3957358 |
| Highest actual frame 14 | 98.0912974 |
| Force-Hermite diagnostic maximum, not DFT | 102.5460259 |

The table's kinetic factors use the highest actual frame, not an interpolated electronic energy. This sampled fixed-coordinate route does not establish the global minimum activation barrier.

```bash
elbow_sha=$(sha256sum elbow-input-ledger.json)
elbow_sha=${elbow_sha%% *}
python analyze_elbow_bundle.py --manifest elbow-input-ledger.json \
  --sha256 "$elbow_sha" --output elbow-reanalysis.json
```

`portability-proof.json` records exact reproduction of all 17 scientific rows, the two-state sensitivity, and refusal of authentic banks paired with the wrong target curve. `native-final-cost-ledger.json` retains exact Native counters and conservative tail reserves. `SHA256SUMS` binds the published payload.

## Open items

- Matched k-mesh and slab-thickness qualification at the computed route peak. Existing adsorption-phase controls do not qualify that peak.

Global activation/barrierlessness and the deferred oxide/donor/permeation cases remain explicit in `pass5-open-items.json`.

## Files in the repository and in the delivery asset

The 57 GPAW `.gpw` banks under `data/` (about 184 MB) are kept only in the delivery asset `g5_recomb-delivery-20261117t225836z.tar.gz`, on the draft release of `vantasnerdan/coherencev1` (SHA-256 `d65ed48e…`). They are gitignored here, and `SHA256SUMS` lists them. Every other file is in the repository.

The reproduction commands pin the pass 4 `transport.py` (SHA-256 `383aad30…`, commit 853c377). Pass 5 later changed the solver in `transport.py` (commit 849ae84) without touching its equations, so the reader rejects the current file by digest. To rerun, extract the pinned version with `git show 853c377:g4_beam/transport/transport.py`. The lead reran both coverage readouts this way, and they reproduce the reported coverages and the ridge factors exactly.
