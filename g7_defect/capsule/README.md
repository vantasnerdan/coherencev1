---
id: coherencev1-g7-defect-screening
description: Static target-conditioned Pd/D BO point curves separating loading and a Pd-vacancy cell
schema: pdd-defect-screening-capsule-v2
status: point-qualified-cell-unqualified
entrypoint: reader.py
---

# Defect and loading: static BO point curves

At R=0.1 A, increasing loading raises the perfect-cell lowering by **1.08864x** but changes the vacancy-cell value by **0.98951x**, while introducing the vacancy changes low/high-loading values by **1.02760x/0.93403x**—small, mixed-direction point changes, not a threefold enhancement, with cell-size convergence unestablished.

These are conditional point estimates, not a certified error bar resolving the small contrasts or a physical no-go for other defects. The sampled values do not reach the approximately 1 keV target. E was not triggered and was not computed.

## Lowering table

Lowering is isolated unit-charge C/R minus the mapped total BO interaction, in eV. Full-precision values, total pair interactions and source bindings are in `results/readout.json` and `curves/`.

| R (A) | A: low-loading Pd | B: reused PdD | C: low-loading vacancy | D: high-loading vacancy cell |
|---:|---:|---:|---:|---:|
| 0.1 | 37.248098 | 40.549800 | 38.276102 | 37.874734 |
| 0.2 | 34.408719 | 35.908817 | 35.004010 | 35.942697 |
| 0.3 | 30.806499 | 32.417775 | 30.707338 | 30.831452 |
| 0.5 | 23.964135 | 23.981516 | 23.847728 | 23.912123 |
| 1.0 | 13.884146 | 13.738790 | 13.743712 | 13.667282 |
| 1.6 | 9.016885 | 9.180974 | 8.931257 | 8.822521 |

## Actual configurations and force provenance

All new geometries were relaxed at fixed inherited lattice a=4.134141684668574 A, with all ions free. D means the hydrogen nucleus with deuterium mass; static BO energies do not introduce a mass-dependent electronic potential.

| Case | Inventory | Geometry raw fmax (eV/A) | Vacancy-centred target distance (A) | Geometric trapped n |
|---|---|---:|---:|---|
| A | 8 Pd, 1 D, no Pd vacancy | 0.04550734 | not applicable | not applicable |
| B | Original pristine PdD baseline, reused without Native runs | inherited qualification | not applicable | not applicable |
| C | 7 Pd, 1 D; one missing Pd | 0.07855682 | 0.34933119 | 1 at radii 1.3/1.5/1.7 A |
| D | 7 Pd, 8 D; fixed inventory after removing one Pd | 0.06395311 | 2.37680058 | 0 at radii 1.3/1.5/1.7 A |

D did **not** retain an occupied VH_n complex in this relaxed local basin. Its tagged target moved away from the vacancy; it was neither pinned nor silently replaced by a different D. Its curve is the actual high-loading vacancy-cell target-conditioned interaction, not evidence for a stable saturated trap, a saturation maximum, or universal absence of trapping.

Axiom2518 admitted a 0.1 eV/A force gate. C uses the actual completed standard-H PAW step36. Its original 80-step run failed the earlier 0.05 gate; the final-minus-selected free-energy drift was **-0.2053987822 eV**, and the original failed record remains unchanged. This drift and the residual forces are provenance, not a certified force-to-screening error bound. D required one fresh-SCF, small-step BFGS refinement from its saved step75 positions. No FIRE continuation was run.

## Hamiltonian, geometry fallback and mapping

- Non-spin-polarized PBE, Pd10, Fermi-Dirac smearing width 0.1 eV. New A/C/D BO grids are actual 24x24x24 with Gamma sampling; C h0.16 uses actual 28x28x28. B retains its separately accepted conventional-cell k2/h0.20 baseline.
- The approved geometry fallback was **standard-H PAW relaxation**, with the Native PAW wavefunction predictor. Every occupied/vacant host and probe energy was then independently recomputed with **Hae** at the fixed chosen positions. No standard-H PAW energy enters a Hae quartet.
- Hae BO uses documented `experimental.reuse_wfs_method='keep'` (`NullWfsMover`) only as an initial-guess mechanism; each electronic state is fully SCF-converged. The earlier attempted all-electron ionic trace used keep but failed relaxation; it is not the current geometry or a successful all-electron force qualification.
- V_pair = E_probe,occupied - E_probe,vacant - E_host,occupied + E_host,vacant. Hosts are actual neutral states. Vacant means removal of only the tagged target D and its electron, with the same background frozen. The charge+1 probe adds a nucleus and **no additional electron** relative to each host.
- Use internal electronic energies E_free - minus_TS. Nuclear regularization is replaced by the accepted **same-zero-G periodic point-minus-Gaussian** block, Hae alpha1=10 Bohr^-2, pair width sqrt(5)/Bohr. No arbitrary offset, far-point fit or bare-only image subtraction is applied.

Probe position is target - R*direction. A/B use Cartesian [110]. C uses (-0.4319765216, 0.6576289795, -0.6171875); D uses (0.2620186206, -0.9380915073, -0.2265625). C/D rays maximize sampled background-nucleus clearance over R0.1–1.6 in the vacancy-side hemisphere: 26 crystal directions plus128 spherical samples before filtering, not a global angular optimum. Full coordinates and clearances remain in the Native records.

## Checks and the explicit cell waiver

The representative C h check holds positions and ray fixed: h0.20 lowering38.27610245 eV versus h0.16 lowering38.70349476 eV, absolute difference0.42739230 eV, below max(5 eV,10 percent)=5 eV.

The 27-site all-ion relaxation was explicitly stopped on Axiom's revised direction. **No converged cell-size check was obtained.** `point_qualified=true` records the actual force-qualified A/C/D electronic results and the C grid check; `cell_qualified=false` remains false.

The included `sensitivity/` result gives a R0.1 lowering of **38.38087628 eV**. It uses an exact seven-Pd/one-D minimum-image patch from C embedded into a 27-site cell, with nineteen far-field Pd left at ideal, unrelaxed positions. Its four actual electronic energies are independently SCF-converged; no forces were computed or ions relaxed for this embedded point. It is not a passed cell-convergence gate, even though its numerical difference falls below the original tolerance.

The one-cell PAW comparison is against the **unconverged** earlier AE trace, after vacancy registration with the same seeded atom labels:

| Position-provenance metric | Difference (A) |
|---|---:|
| Tagged D displacement | 0.62990908 |
| Pd RMS displacement | 0.27666399 |
| Pd maximum displacement | 0.52097740 |

This is not agreement between two converged setup geometries. The source-bound provenance artifact records the original failed AE force and the numerical global-origin sensitivity.

## Read and verify

`python reader.py CAPSULE_DIRECTORY` verifies inventory, ORIGINAL/EXPORTED bindings and actual saved neutral-energy/electron/geometry/Green calculations, reproduces the point table and C h result, and checks any optional sensitivity without promoting it to convergence.

`python portability_proof.py ARCHIVE OUTPUT_DIRECTORY` uses the real extracted archive with GPAW, ASE, SciPy and plotting imports blocked. Its negative variants test a byte alteration and a rebound extra-probe-electron error. No prepared state, U(0), rate, drag, core continuation or far-tail completion is supplied. Material-pair files retain the `coherencev1-target-conditioned-pair-v1` schema and exact reused B bytes.

## Three omissions

1. **Converged cell size and sampling:** no converged two-cell qualification; any larger-cell result is unrelaxed sensitivity only, with finite k sampling and vacancy density still unresolved.
2. **Alpha-Pd volume:** the common expanded PdD lattice was retained; alpha-Pd near a=3.89 A was not computed.
3. **Model and ensemble:** static non-spin-polarized PBE/PAW-geometry/Hae-BO only, without a thermal ionic ensemble, ZPE, dynamics or fusion rates; one fixed-inventory local basin is not a global basin search or saturation capacity.
