# Replay: frozen PdD refit through 0.48 Å

This is an unpublished, scoped scientific capsule, not a replacement of the original G10 archive. Extract its repository-relative members into an empty directory. `ARCHIVE_FILES.txt` is the exact whitelist; `SHA256SUMS` binds every payload member. The archive itself is listed separately in `capsule.json` because it cannot hash itself.

The input archive supplied by Giuliano is asset 789081289, SHA-256 `d488fd0e11c5b2c82cd8a438e6174e49a56353b627bd7054641c7c51f133a74f`. The exact selected signed JSON is `inputs/bulk-wells.json`, SHA-256 `000612c601ee3cce30fc6c904ba527607ef3680ca2cec4d54ecf2d674d364a82`. `inputs/old-readout.json` preserves the prior joint coefficients, geometry, occupied-index registration, local K and neighbour shell means. Original near-point/probe bytes, provider semantics, and archived positive quantum baseline are included. The raw wider archive is not re-packed here; the exact consumed JSON is. Native and inversion-derived rows remain distinguishable. No NiH data are used even though the signed provider export also contains NiH rows.

## Environment and commands

Run from the extracted repository root using Python 3.12 with NumPy, SciPy and Matplotlib. The author's interpreter is `python`; no environment installation or network access occurs. Set one thread for each numerical runtime:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
PY=python
D=g10_quantum_site/bulk-refit-048
```

The numerical resolutions and settled replay commands are recorded in `results/summary_main.json`, `results/summary_check.json`, `results/summary_isotopes.json` and the final page. Recomputing overwrites only this refit directory's results.

```sh
$PY -B "$D/solve.py" --npar 10 --order 48 --tag main
$PY -B "$D/solve.py" --npar 12 --order 64 --tag check --compact
$PY -B "$D/solve.py" --isotopes H T --npar 10 --order 48 --tag isotopes --compact
$PY -B "$D/coupled_sensitivity.py" quadratic quartic
$PY -B "$D/deliver.py"
$PY -B "$D/verify.py"
```

The coupled replay deliberately includes the original `model.py`, `coupled_momentum.py`, exact old readout path, trajectory summary and both retained closed-path NPZ files. It does not require i-PI, SSCHA, a new trajectory, electronic-structure tools, or the original checkpoint. Saved neighbour fields are actual co-beaded joint samples, not a mean-field replacement.

The settled run used Python 3.12.3, NumPy 2.5.5 and SciPy 1.18.1. `results/verification.json` records the consumer check: full 3D position/momentum normalization, exact table readback, signed-plane marginal normalization/symmetry, the recreated archived positive baseline, force/splice evidence and consumed coupled witnesses. It does not rerun coupled sampling or claim independent electronic/physical validation.

## Data and operator contract

`well.energy_force(x, variant)` accepts arbitrary vectorized Cartesian `x (...,3)` in Å and returns `(local_energy (...), local_force (...,3))` in eV and eV/Å. Variants are explicit: `old`, `quadratic` (primary), `quartic` (alternative). The old source has zero F0 and the new fit leaves that, K, geometry and neighbour graph unchanged. This API is local only: do not add or subtract neighbour K twice.

The canonical isolated operator is `H=−ℏ²∇²/(2m)+V(x)`, with distinguishable localized nuclear mass and no exchange/hopping, fixed PdD electronic surface and cell. H/T change only the isotope mass in this same surface. Temperatures 77/300 K change the ensemble only. Four Cartesian parity representatives with multiplicities 1/3/3/1 restore all eight parity sectors, not one one-dimensional ray. Gauss-Hermite quadrature assembles the full nonseparable local potential. Every saved density is a positive canonical sum of squared eigenfunctions or their unitary Fourier transforms. The momentum coordinate is `k=p/ℏ` in Å⁻¹, not classical velocity and not additional incoming energy.

`*_distributions.npz` contains angle-averaged density, per-radius measure, integrated exterior probability, and [100]/[110]/[111] on-ray densities for position and momentum. Main D files also save full positive-octant 3D coordinate/momentum grids. Reflection reconstructs the full distribution exactly under this even model; cubic permutations provide redundant equality, not independent physical measurements. Spectrum NPZ files retain the actual eigenvectors/energies and basis length, hence replay the distribution at arbitrary vectors. Signed [111] `P(s)` is a plane-integrated marginal; it is not the on-ray density. Positive s uses +[111]/√3; no physical slab orientation/offset is inferred from a bulk cell.

`tables/*position.csv`: radius Å; spatial density Å⁻³; radial measure `4πr²ρ̄(r)` Å⁻¹; tail dimensionless. Momentum analogues: k Å⁻¹, density Å³ and radial measure Å. Their maxima and one-millionth crossings differ. Off-ray values are the declared cubic closure, not measured coverage of an r≤.48 Å ball. Every exterior prediction remains a continuation even after numerical convergence.

`coupled_position.csv` is the tagged site's genuine joint bead marginal with shell/bin measures and paired contiguous-block SE. `coupled_characteristic_111.csv` is the midpoint displaced-path characteristic of the signed momentum component. No Fourier inversion or coupled far-tail contour is claimed. Reweighting changes one chosen site's local well in the old joint 32-site host; reduced saved paths cannot reconstruct an all-cell change. Frame/block ESS measures retained participation, not a guarantee of unseen overlap or independence. Existing bead discretization and residual time correlation remain separate limitations.

The two new continuations agree exactly through .48 Å and match E and radial slope there; the softer quadratic uses bare K curvature, while the quartic retains the old positive ray-quartic curvature and higher Taylor terms. Their spread is a structural model sensitivity, not a statistical confidence interval or a DFT-error estimate. Numerical differences under the stated basis/quadrature comparison are recorded separately. Finite-comparison contours are labelled resolved only when the candidate crossing is stable and the local difference is below 20% of 10⁻⁶ of the peak; this is not a rigorous enclosure.

No density, source profile, S factor, arrival dynamics or yield is fitted. No 0.80 Å well data, connected pair cut, surface/vacancy response, damping or full all-cell refit is supplied. The optional direct g8 joint-excitation consumer is not asserted by this capsule.
