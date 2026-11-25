# NiH: one signed [111] cut through 0.80 Å

**Decision.** For physical H at 77 and 300 K, replacing only the old [111] ray softens the two tested cubic completions: isolated MSD rises by 0.43–0.71% / 0.53–0.88%, kinetic energy falls by 0.28–0.48 / 0.24–0.42 meV and the 0.48 Å exterior probability rises by 3.6–5.8% / 4.2–6.7%. The tagged finite-bead lattice estimates agree in sign. Varying the unmeasured angular spread changes the size of these shifts; changing only the continuation beyond 0.80 Å is immaterial to the moments. This is a contrast between specified models, not a physical angular bound.

## Input and model

The positive [111] rows are static frozen-neighbour PBE free-energy rises at the signed-export settings (500 eV, Γ-centred 3×3×3, unpolarized, 0.1 eV smearing), relative to the material's own zero. The negative rows are inversion-derived, with the same energies and opposite forces **on the same fixed +[111] axis**, not independent calculations.

| Positive s (Å) | New E (eV) | New fixed-axis F (eV/Å) | Old E (eV) | Old fixed-axis F (eV/Å) |
|---:|---:|---:|---:|---:|
| 0.24 | 0.04933334 | −0.41217157 | 0.05024857 | −0.42653734 |
| 0.48 | 0.19679311 | −0.80523798 | 0.21222522 | −0.94666585 |
| 0.80 | 0.50052896 | −0.98812778 | 0.66346307 | −1.94751930 |

Reference: `g10_quantum_site/model.py:48–73`, `full('Ni','joint')`. The G3 NiH dense-cell on-site `K=1.7122451 eV/Å²` is not collective Γ; the 32-H graph is unchanged (`ΣJ=0.6759926 eV/Å²`). The independent late G3 [100] joint energy/force coefficient is `q100=16.1656969 eV/Å⁴` (fixed harmonic term; old −/+0.03 Å probes and +0.24 Å displacement context, no new [100]); [111] used `q111=1.1283657 eV/Å⁴` from `g2_refit/computed_inputs.json` (entry_401; coordinate decision entry_458; older settings per `g3_loading/README.md:72–82`). With `B=(3q111−q100)/2`, [110] has the old cubic **prediction** `q110=(q100+B)/2=4.8876985 eV/Å⁴`, not a new or independent [110] well. Input hashes and all six signed anchors are in `support.json`.

Use cubic Hermite in `u=r²` through the three positive anchors, with slopes `dE/du=−F/(2r)` and `(E,dE/du)=(0,K/2)` at zero. It is nonquartic, exact in energy and radial slope on all six signed anchors, and has old Cartesian near-origin Hessian `K I`. Outside 0.80 Å continue by `E(.80)+E′(.80)(r−.80)+C(r−.80)²/2`; `C=K` normally and `C=8K` for the `w_stiff_exterior` sensitivity. Both splices are C¹; no smooth high-k asymptote is inferred. For `l_i=x_i²/r²`, `w=27 l₁ l₂ l₃ ∈ [0,1]`; the alternatives are `V_old+w (E_new111−E_old111)` and `V_old+w²(E_new111−E_old111)`. The weight and its gradient vanish on coordinate planes, so the old [100]/[110] energies and forces remain exactly. The increment is O(r⁴) at the origin. The old [111] quartic is the minimum over directions, so one can write each potential as a nonnegative old angular excess plus a convex blend of two positive confining rays. The Hermite ray itself is increasing on all intervals (`min dE/du=0.61758 eV/Å²`). No angular observations away from this ray fix the blend.

## Numerical readout

`observable_table.json` is an array of shared rows. `support.json` separates source/probe data, numerical comparisons and paired uncertainty. Isolated results use the same Cartesian parity-reduced Hermite Galerkin operator as the accepted 0.48 Å calculation, but actual proton mass, NiH old quadratic/quartics, and the NiH oscillator scale (0.2029207 Å). The archived old `.npz` stores readout grids, not eigenvectors; the legacy Hamiltonian was reproduced on its own grid before comparing sharper operator moments. Its old integral-grid MSD is about 2.10×10⁻⁴ Å² above the operator moment, comparable to the smallest refit shifts; the 10-to-12 parity basis change in the old operator MSD is below 10⁻¹¹ Å². Compare matched results, not the published coarse integral against a new operator moment.

The independent energy differentiation agrees with the signed forces within 3.6×10⁻⁶ eV/Å, with left and right slope estimates differing by at most 6.6×10⁻⁶ eV/Å at the probe step. Reversing the sign creates a discrepancy of ~1.98 eV/Å; off-axis force also matches Cartesian differentiation. The combined 10→12 basis, 42→52 GH and readout refinement for each angular closure changes MSD/kinetic by at most 1.1×10⁻⁷ Å²/eV, `P(r>0.48)` by 2.5×10⁻⁷ absolute and `P(r>0.80)` by 3.6×10⁻¹⁰ absolute. This is empirical discrimination, not a numerical enclosure.

Isolated `P(r>0.80)` rises from 2.66×10⁻⁷ to 4.21–4.92×10⁻⁷ at 77 K, and from 4.69×10⁻⁷ to 7.99–9.51×10⁻⁷ at 300 K for `w²` and `w`. At fixed `w` the stiffer continuation reduces the new tail by 0.072×10⁻⁷ / 0.183×10⁻⁷ (about 1.5% / 1.9%); the completion spread is about 10× / 8× larger. This tail probes an assumed continuation, not probabilities inferred directly from DFT. The stiff exterior moves MSD by less than 5×10⁻⁸ Å².

`distributions.npz` keys begin with `<variant>_<temperature>K_` (`old`, `w`, `w2`, `w_stiff_exterior`). They contain both canonical densities built from nonnegative sums of squared eigenfunctions (Fourier phases `(-i)^n`, `k=p/ℏ`, in Å⁻¹). The spherical density is the angular average with the original 3D density units, the radial measure is `4π r²` or `4π k²` times that average, and `*_tail_to_window` integrates only up to the saved window (1.2 Å / 48 Å⁻¹). The line columns sample the original volumetric density along [100], [110], [111]; the signed [111] **plane-integrated** position marginal separately integrates two transverse coordinates, in Å⁻¹. Neither is an arrival-direction distribution. The radial measure peaks are near 0.208/0.212 Å and 5 Å⁻¹ at 77/300 K for all these models. At 0.48 Å, the volumetric line-density ratio [111]/[100] changes from 2.47 to 2.64–2.67 at 77 K and from 2.63 to 2.82–2.87 at 300 K. The isolated momentum curves are not a tagged lattice momentum readout and high-k tails are not certified.

## Tagged finite-P result

Reweight the old primitive trajectory by `exp[−β mean_beads(V_new(x)−V_old(x))]`, changing exactly one tagged on-site potential with the same mass, springs, graph and original co-beaded neighbour field. The kinetic column is the tagged centroid-virial scalar including that field, not the old all-site mean or a momentum distribution. Pair the estimates on nine contiguous time blocks (block-influence standard errors, not confidence bounds). Effective participation remains ≥798 of 801 frames and nearly nine of nine blocks. For `w` the paired MSD shifts are +0.000467 ± 0.000056 Å² (77 K, P=128) and +0.001049 ± 0.000200 Å² (300 K, P=32); `w²` gives +0.000295 ± 0.000038 and +0.000663 ± 0.000124 Å². The corresponding `w` shifts in `P(r>0.48)` are +0.000659 ± 0.000108 and +0.001929 ± 0.000550. No saved bead reaches 0.80 Å, so tagged 0.80 Å probabilities are **null**, not zero; the saved ensemble cannot distinguish the exteriors. The archived equal-time path reweight does not supply a tagged momentum distribution, all-cell refit, continuum-P extrapolation, or dynamics.

## Replay

Use a local environment with NumPy and SciPy (used: 2.5.5 and 1.18.1). From a checkout or extracted capsule, set paths for the prior G10 assets and the verified signed provider export:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
python path/to/nih/replay.py --reference-root path/to/reference-root \
  --signed-export path/to/bulk-wells.json
```

The reference root follows original repo-relative paths: `g10_quantum_site/{model.py,exact.py}`, `g10_quantum_site/bulk-refit-048/{solve.py,well.py,inputs/{old-readout.json,bulk-wells.json}}` (shared solver imports); the three coefficient-source JSON files in `support.json:model.source_sha256`; the two old `g10_quantum_site/output/exact/NiH_{77,300}K_joint.npz` reference grids; `g10_quantum_site/output/coupled/{trajectory_summary.json,NiH_77K_paths.npz,NiH_300K_paths.npz}`. Replay never modifies these inputs or reruns trajectories. It writes only `observable_table.json`, `support.json` and `distributions.npz` next to the script.
