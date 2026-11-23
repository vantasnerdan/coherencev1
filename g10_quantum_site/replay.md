# Bulk replay and handoff

Use the already installed i-PI 3.3.1 / ASE 3.29 / SSCHA + CellConstructor 1.6.2 environment, with NumPy 2.5.5 and SciPy 1.18.1. The saved input sources belong to pinned repository revision `57609d0dcfd85a51bd5f4a671f311a2d289dc050`; `output/input_inventory.json` binds the actual consumed files and coefficients. No GPAW, private checkpoint or new electronic-structure calculation is used.

From the repository root, with `PYTHON` naming that environment's interpreter:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$PYTHON" g10_quantum_site/deliver.py exact
"$PYTHON" g10_quantum_site/deliver.py coupled
"$PYTHON" g10_quantum_site/deliver.py derive
"$PYTHON" g10_quantum_site/deliver.py verify
"$PYTHON" g10_quantum_site/deliver.py capsule
```

`exact` and `coupled` are numerical producers, not prerequisites for consuming saved tables. `derive` consumes their persisted outputs and the existing actual-reader case; `verify` is the final changed-consumer check without rerunning trajectories or eigensolves. Required full coupled references are D in PdD and H in NiH, each at 300 and 77 K. H/T marked `Dhost` are selected-site mass tracers in the D-isotope host, not homogeneous PdH/PdT. The separately named genuinely sampled H-host reference is a different physical case. Metal geometry and Born–Oppenheimer inputs stay unchanged.

The exact finite-bead tracer change of variables is `x' = xc + sqrt(m_D/m_I)*(x-xc)`. It leaves the selected ring's spring action invariant; its constant Jacobian cancels in normalized expectations. Remaining log weights are `-(beta/P)*sum[V(x',host)-V(x,host)]`, including the actual neighbour field once. Momentum uses transformed paths, normalized importance weights and the new isotope mass in its free-end factor. Host coordinates stay coupled, not frozen; their isotope remains D. Frame/block participation and displacement-dependent overlap are exposed, not assumed adequate in a low-overlap tail.

Saved full-site trajectories are reduced to the selected site's unwrapped bead paths and actual neighbour field, sufficient to replay one-body momentum and conditioned-coordinate estimators. Complete physical graphs, occupied original indices and signed force axes remain in the inventory. Partial/strain output is a separate analytic Gaussian model-contrast table, not a fabricated anharmonic distribution or a measured strained local well.

## Table entry points

- `output/tables/*exact_position_band.csv` and `*exact_momentum_band.csv`: both radial meanings, integrated tails, three directional densities, fit-continuation min/max and finite-basis comparison.
- `output/tables/*coupled_position_band.csv` and `*coupled_momentum_band.csv`: non-Gaussian i-PI coordinate/momentum estimators, continuation bands and explicit uncertainty/detection columns. Signed Fourier estimates are not clipped or renormalized.
- `output/exact/PdD_300K_joint_g8_joint.json` (also 77 K): list-valued `states` containing scalar `excitation_eV`, eV `kinetic_eV` arrays and **joint** positive `probabilities`, summed over all states/nodes to one. Finite radial-bin quadrature width is recorded; it is not a factorized excitation–kinetic ansatz. Supply the list to `BoundSpectatorKernel.fold('PdD', energies_cm_keV, states, removal_eV=...)`. Coupled momentum has no known joint excitation measure.
- `output/tables/*g9_normal111_band.csv`: plane-integrated signed coordinate marginal in Å⁻¹, plus separately labelled on-axis line density in Å⁻³. Positive s points along declared inward +[111]/√3; reverse it for an opposite surface. [100]/[110] context is included. This is a bulk frame, not a supplied slab orientation or target offset.
- `output/tables/g9_mapping_example.csv` imports the actual marginal and applies `z_target=Z0+s`, `arrival_shift=s/v` with explicitly declared example Z0 and speed. No g9 source file changes.
- `output/consumer_cases.json`, `output/verification.json`: actual reader/import result and the scoped final verification.
- `output/coupled/sscha_reference.json`: installed SSCHA local Hessian and gradient, force-rounding change. `gaussian_reference.json` is separately the coupled analytic Gaussian variational reference, not a tail replacement.

Every position distribution carries its computed-cut endpoint. For full loading it is 0.24 Å on the available positive cuts; [110] remains a cubic prediction and negative large displacements an even-extension assumption. Partial loading has only ±0.03 Å measured on its named signed axis, with transferred quartics and assumed transverse closure. Strained Γ data are collective: their local-well conversion requires frozen J, and no strained local large-displacement support is implied. The Fourier tables compare potential continuations over the whole k range; a coordinate endpoint is not a momentum endpoint.

The fit-continuation band spans energy-only, force-only and joint fits; it is a sensitivity band, not a complete uncertainty enclosure. Coordinate quadrature, finite basis, bead count, trajectory correlation, Fourier window and continuation effects remain separate. A zero-count histogram or a signed noisy Fourier value does not resolve a physical tail.

## Archive whitelist

`ARCHIVE_FILES.txt` enumerates **every** included relative path, one per line, and `SHA256SUMS` binds those bytes (excluding the checksum file itself). The allowed set is exactly these sources/documents and numerical output tree:

```text
g10_quantum_site/model.py
g10_quantum_site/exact.py
g10_quantum_site/coupled.py
g10_quantum_site/conditioned.py
g10_quantum_site/coupled_momentum.py
g10_quantum_site/coupled_momentum_notes.md
g10_quantum_site/deliver.py
g10_quantum_site/README.md
g10_quantum_site/replay.md
g10_quantum_site/tails.png
g10_quantum_site/neighbours.png
g10_quantum_site/output/ [only enumerated .json, .npz, .csv files]
g10_quantum_site/ARCHIVE_FILES.txt
g10_quantum_site/SHA256SUMS
```

No run logs, bytecode, temporary files, private paths, other directories, upstream input duplication or surface/vacancy production belong in the capsule. `capsule` produces local `g10_quantum_site/g10_bulk.tar.gz` and prints its SHA-256. Upload/draft-release publication is the receiving owner's step; this directory does not push or change main.
