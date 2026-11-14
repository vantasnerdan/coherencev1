---
description: Registered-site PdD loading columns with vacancy-facing probes and explicit bond normalization
config_schema: loading-pd-config-v2
result_schema: loading-pd-local-v2
geometry: loading_geometry.py
producer: loading_pd.py
consumer: loading_readout.py
---

# Loading-column protocol

The remaining Pd and D nuclei occupy the registered, uniformly scaled rocksalt sites. These are not relaxed defect configurations. The common zero sample retains all forces, including nonzero forces on the selected D sites. Electronic H species represent D; force constants are mass-independent.

The reference follows the approved Schirber–Morosin77K beta-PdD relation a(x)=4.084 exp[0.044(x-1)] A. The x=1 endpoint is author-extrapolated; .875 is model interpolation; .75 is an authorized extrapolation below the abstract's measured range. Pointwise lattice uncertainty is unavailable, not zero. The electronic PBE free-energy calculation uses Fermi smearing0.1eV, not a77K electronic temperature, and includes no nuclear zero-point motion.

## One registered axis per probe

| Probe | Axis in original cubic coordinates | Positive displacement |
| --- | --- | --- |
| clean originalH25 | (1,0,0) | +100 |
| near originalH1 | (-1,0,+1)/sqrt(2) | Toward vacant originalH5 |

Nearest D-site vacancies lie on110-family rays. The near axis is therefore not100 or111. Each config records the actual unit vector, label, target vacancy and periodic image. The deterministic choice is the lowest original-index vacant nearest neighbour; all three commissioned configs select H5. This is a requested vacancy-facing sample, not a measured principal soft direction. The .875 spread/bunched near probes have2 vacant NN; the .75 near probe has4. The clean probes have0.

One common zero and actual -0.03/+0.03A pairs for clean and near give five SCFs for the .875 spread and .75 bunched cases. Axiom's revised sampling plan (Fleet798) commissions only the zero plus near pair for .875 bunched: three SCFs, explicitly no clean pair. No reduced-loading quartic or second-axis sample is implied. Never reconstruct missing signs in vacancy cells. The fitter requires the declared probe set; pass `--near-only-config x0875-bunched4` for that intentional omission, not for missing data in a full case.

## Native and frozen interfaces

Run only inside an admitted bounded Native/systemd unit with the existing all-exit/accounting path; these commands do not admit resources themselves. Use the installed GPAW-capable Python through GPAW_PYTHON, single-thread BLAS, genuine k-point-only MPI, one pinned rank/core, and actual world/rank checks.

```sh
bash run-loading.sh zero --config loading-configs/x0875-spread4.json --mesh 3 --ranks 6
bash run-loading.sh probe --config loading-configs/x0875-spread4.json --mesh 3 --ranks 6 --probe clean
bash run-loading.sh probe --config loading-configs/x0875-spread4.json --mesh 3 --ranks 6 --probe near
python loading_readout.py --config loading-configs/x0875-spread4.json --mesh 3 --probe clean --output NEW.json
```

PW500/PBE/Pd10, nonmagnetic, Fermi0.1eV, nbands-20, Davidson2, Mixer(.05,5,50), max80; energy1e-7/density1e-4/force-change1e-4. Spatial point-group reduction is disabled. Every displaced state reads the same actual zero checkpoint with a full self-consistent density. No frozen-density or private-density transplant. The mesh is campaign-fixed k3, not a new mesh-convergence claim.

Frozen reading needs matching source/config/point JSON files and NumPy/ASE, but no private GPW/PAW bytes. Native restart requires exact checkpoint bytes. Complete points publish atomically; orphan checkpoints require explicit disposition. Native exit success is not numerical-accuracy acceptance: retain energy/force curvature consistency, force-evenness and full-ion acoustic diagnostics before adopting a model.

## Force-column and bond convention

For the recorded unit axis d, the actual central column is Phi_(i alpha,jd)=-(F_i alpha(+h*d)-F_i alpha(-h*d))/(2h). Kii_dd is its source vector dotted with d. Offsite pulling J_dd is minus the corresponding neighbour vector dotted with d. The near column must not be labeled Jxx. Raw Cartesian vectors and the axis remain available; scalar/isotropic extension is a model closure, not a measured full tensor.

Each finite-cell shell retains occupied sum, occupied multiplicity, full-reference distinct-site multiplicity, mean per existing bond, original indices, vacancy membership and nearest-image offsets. The fitting/kernel convention is mean_pulling_Jdd_eV_A2: actual occupied sum divided by occupied multiplicity. Apply lattice occupancy n_i*n_j once. Feeding the occupied sum as a full undiluted kernel would dilute twice. Empty-shell means are null, not zero-strength invented bonds. Mref times the mean is only a full-reference-equivalent diagnostic; farther folded shells are not unfolded infinite-FCC coordinations.

Kii minus the sum of all actual occupied J is an H-only heterogeneous column sum, not a global Gamma eigenvalue or an all-ion acoustic zero. Source-image terms remain folded in the source diagonal. No acoustic correction is imposed.

## Loading versus additional elastic strain

The fitted spacing argument is bulk a(x), coupled to composition along this registered path, not an independent strained EMT local-spacing knob. Loading coefficients already include that isotropic reference change. For PartA's absolute4.084A-reference tensor law A, let s=a(x)/4.084 and E_absolute=s*(I+E_loading)-I. Additional elastic contrast is A(E_absolute)-A((s-1)I), not minusA(0). This makes the extra correction zero at the loading bulk reference. The source-bound reference/sensitivity [receipt](receipts/loading-strain-reference-contract-20261114T0010Z.json) is Fleet693.

Transferring this collective x1 strain contrast into local Kii at vacancy loading requires the explicit frozen-neighbour-sum approximation; it is not independently measured here. Retaining the new x1 quartics at reduced loading is likewise a model assumption, not loading-quartic data.

## Reproducing a completed-case figure

The figure consumer reads the actual point files bound by the two frozen readouts. It explicitly removes the measured origin-force tilt for display; a positive tangent Kii does not make the registered origin a stationary minimum. Lines are central harmonic guides only inside the measured displacement window, not loading-quartic fits.

```sh
python plot_loading_case.py loading_Pd/x0875-spread4/k3/readout-clean-20261114T0421Z.json loading_Pd/x0875-spread4/k3/readout-near-20261114T0500Z.json --output loading_Pd/x0875-spread4/k3/registered-probes.png
```

The first complete case has actual GPW/PAW-free replay and a wrong-spacing refusal in `receipts/loading-clean-admission-20261114T0422Z.json`. The independent diagonal-axis projection, wrong-axis sensitivity and frozen near replay are in `receipts/loading-near-admission-20261114T0501Z.json`. Nonzero origin forces are retained in the data, not silently removed from the physical model.

For a commissioned near-only case, pass its one actual near readout to the same figure command. No duplicate clean curve or synthetic point is supplied.

## Anchored descriptive law

`loading_fit.py` fits `y=y_x1+c_epsilon*(a_reference/4.084-1)+c_v*n_vac_NN` with the measured x1 anchor fixed, for Kii_dd and each matched-shell mean J_dd. The five reduced-loading probes are b clean/near, c clean/near, and d near only. Equal-weight fits retain signed residuals, scaled-design conditioning and leave-one-configuration-out predictions when identifiable; these are descriptive errors, not confidence intervals or new SCF/mesh bounds.

Axiom852 authorized one post-contrast alternative, selected after observing the spread/bunched difference: add only `c_pair*vacant_NN_pair_count`. This count is the number of unordered pairs among the source's vacant nearest-neighbour sites that are themselves nearest neighbours, using minimum images in the registered2a cell. The actual near counts are0 for .875 spread,4 for .75 bunched,1 for .875 bunched; clean and x1 counts are0. `--motif` selects that three-slope alternative; no interaction or further term is fitted.

The two actual saved models are `loading-law-base.json` and `loading-law-motif.json`. `loading-law-comparison.json` contains both residual/holdout arrays and independent rational normal-equation checks. Kii RMS changes from.04501 to.04276eV/A2, but maximum configuration-holdout error worsens from.10830 to.40442eV/A2. Reduced training residual is not evidence of improved prediction. Some sparse outer-shell holdouts are ill-conditioned or unidentifiable. No universal two-descriptor law, predictive accuracy guarantee or preferred variant is asserted. The initial v1 source/output epoch remains in `evidence/loading-affine-v1-before-motif-20261114T0755Z/`; it is not the current API.

The current v2 API is `loading_law.evaluate(n_vac_NN, a_reference_A, law, vacant_NN_pair_count=None, extrapolate=False)`. Motif evaluation requires the actual pair count explicitly; the baseline has no pair term. Coefficients are a named `coefficients_eV_A2` mapping in each observable model. Unidentified coefficients and out-of-hull queries are refused; explicit extrapolation returns the affected names. Domains are separate two- or three-dimensional sampled convex hulls, not rectangular validation boxes. No stiffness clipping or invented empty-shell bonds. The spacing coordinate is the composition/reference-spacing path, not extra EMT strain. Native origin forces remain separate data; tangent fitting neither centers the potential nor measures loading quartics.
