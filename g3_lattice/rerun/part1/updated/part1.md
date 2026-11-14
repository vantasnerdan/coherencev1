## Result

Adopted fixed-k3 Pd local-well refit (762): local 0→1 energies are **74.034667 meV (PdH)** and **48.883453 meV (PdD)**; the coupled Γ response poles are **62.264359/39.826501 meV**. The H/D ratios are 1.514514 locally and 1.563390 collectively. The 55.8/37 meV neutron band descriptors are not identified Γ modes, so this comparison is not a Γ validation.

The frozen-local-variance Hartree + actual PartA contrast construction gives tensile static onsets, at 0/300 K, **7.4565/8.7699% (PdH)** and **5.2704/7.0506% (PdD)**. These are extrapolations beyond the computed ±1% family, not SSCHA or established material phase boundaries; zero-strain stiffness remains positive.

The joint 80-proton, three-cell central vacancy-force fit is **C1=0.07073357956, C3=−0.0001342605521 eV/Å**. Physical half-cell image cancellation makes C2 unidentified; **C2=0 is an explicit assumption**, not a measured zero. Vector-RMS calibration error is **0.00412955 eV/Å** (component RMS 0.00238420); leave-one-cell-out vector RMS for b/c/d is **0.00298770/0.00863463/0.00425596 eV/Å**. These are calibration/held-out errors on these motifs, not bounds for arbitrary arrangements.

Completed measured-F0 versus fitted-F0 template comparisons give static displacement-field relative differences **2.44/3.55/3.05% (b/c/d)**. The absolute separated-Ω field differences are **0.00063145/0.00036450/0.00027716 Å**. Phase-aligned relative differences are 2.04/72.50/9.53%; the large c/d relative errors concern very weak Ω signals and must not be read without their absolute sizes. “Exact template” means the measured reference force field, not a measured full Hessian or a DFT-relaxed minimum. This static signal is **vacancy-driven displacement**, not additional spontaneous ordering.

The archive contains all quantum and bias numeric summaries, the original scientific input readouts and runnable sources, and the completed three-cell comparison JSON. No GPAW checkpoint/archive, raw trajectory bulk or Native log is included.

## Assumptions
Local cubic anharmonic wells plus same-profile scalar J; static closure freezes local variance and additional elastic ΔΣJ, while the vacancy law uses periodic image-aware central forces with C2=0 and explicit scalar/orientation transfer.

## Script
With Python and NumPy, from the extracted packet root:

```sh
(cd pd_x1_quantum && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python pd_x1_refit.py)
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python fit_vacancy_bias.py \
  --case new_inputs/loading_b_793 --case new_inputs/loading_c_827 \
  --case new_inputs/loading_d_850 --out vacancy_bias_reproduced.json
```

`updated/part1_summary.json` binds the completed comparison to the actual 30-case batch; the comparison is a reused result, not another simulation.
