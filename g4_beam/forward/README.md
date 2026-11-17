# Forward yields and the knob-table consumer

This delivery uses exactly the wide-bounds layer, `sigmoid_tt`, fixed-energy fast neutrals (fitted energy, with 5 keV nested), Czerski thermal spike, and a no-additional-plateau control. `sigmoid_xs` exists only as the independent Eq7 curve check. Standalone fitted-width resonance and arbitrary free floors are removed.

## Reaction, independent route, and corrections

`sigma` is the unchanged Eq7/Bosch–Hale cross section. `reaction_sigma` adds the COMMON proton nuclear background from Czerski v3 literal Eqs14/17: Gamma_p=40 meV, phase115 degrees, channel radius7 fm, E_R->0. It is not another plateau candidate. Neutrons remain Bosch–Hale because their resonance width/phase are unspecified; `reaction='mirror'` is the single explicitly separated sensitivity, never a fitted width. Substitution of the printed Eqs15–16 into Eq2 appears to differ from printed Eq17 by a factor pi; the literal Eq17 convention tied to the source's fitted constants is retained, without a compensating width.

`second_route.py` NEVER imports `model.py` or `fits.py`: it implements ICRU/ZBL stopping, adaptive lab-energy quadrature, Gaussian/molecular averaging, independently integrated angular-first-moment layer geometry, current normalization and curve roots. `second_route_reference.json` freezes only the required original comparisons. The original Eq7 evidence is in `second_route.csv`; the changed common-reaction evidence is in `second_route_common.csv`, with an SI evaluation of Eq17 independent of the maintained keV/fm expression.

Comparison tolerances are 0.005 keV for roots, 0.01 keV for bulk Ue, and 0.05 absolute chi2 for one-material frozen-parameter comparisons. Adaptive integrals request relative tolerance2e-8; these comparison tolerances are not experimental uncertainties. Exact R=1e-4 cross-section roots (bare/Pd/Ti/TiD) are **1.90728/1.17925/1.63879/1.08496 keV**, not the quoted raster-row1.92/1.20/1.66/1.10. The original curve check targeted a row above the bottom axis, not exactly1e-4. Exact sqrt-stopping thick-target roots are2.13027/1.39512/1.86303/1.29637; actual-material ICRU/ZBL single-D roots are2.13999/1.41243/1.87662/1.32006. No Eq7 parameter was retuned to a pixel position.

Eq7-effective independent Pd/TiD bulk fits give1.26842/1.33751 keV as plotted,0.69549/0.84143 with one cup-current power. `bulk_effective.json` preserves that unaffected branch. `bulk_screening.csv` gives every material's effective and fixed-reaction bulk Ue, Jacobian uncertainties, the one-sigma shift flag, and the mirrored-neutron sensitivity. These use the same2.5–6.5keV fit and shared proton/neutron stitch. In the fixed-reaction result PdD under M0 moves by-0.161keV (1.09 old sigma); other listed shifts do not exceed one old sigma. Mirroring neutrons moves bulk Ue by-0.0061 to+0.0240keV and the6.5keV physical neutron reference by-3.32% to-4.34%, at matched baseline Ue.

The comparison now exposes **all six hypotheses M0–M5 (24 material rows)** in both conventions, reusing the completed fits. The two-current-power M3/M5 rows give effective/common Ue(Pd/PdD/Ti/TiD)=0.253475/0.254375,0.804722/0.805181,0.030363/0.031073,0.448796/0.450141keV. The complete comparison still has only the PdD/M0 one-sigma warning above. Mirrored-neutron sensitivity is computed only for M0/M1; other rows are explicitly `not computed`, with JSON null/blank CSV values, not zero or inferred alias results. `bulk_comparison(results, mirror_sensitivity=completed_rows)` reuses matched existing sensitivity rows; its maintained 24-row JSON/CSV write and Boolean parse were exercised in0.00228s without repeating a fit.

The maintained log grid now covers32keV lab: the allowed fitted neutral energies no longer silently clip at the old16keV endpoint. Below the integration cutoff `Y_rows` returns zero rather than negative extrapolation. Legacy extra-check drivers/outputs and obsolete candidate inversions are removed; an archive overlay must apply `removed_paths.json`.

## Physical callable

```
from g4_beam.forward.consumer import forward_yields
out = forward_yields(E_cm_keV, profiles, Ue_keV=Ue,
                     candidate="none", measurement={"efficiency": efficiency})
```

`profiles` is an exact-energy mapping or callable returning `(layers, nads_m2, state)`. A layer has **seven** SI fields `(z_top_m,z_bottom_m,n_host_m3,n_D_m3,Z,M_amu,n_metal_m3)`; n_host includes nonmetal atoms in oxide effective hosts. Bind each profile's nominal energy, condition, temperature and time; a transport time index is not an energy index. Profiles are frozen, not evolved or rerun here. `Ue_keV` may be scalar or callable `(E_cm,z_m)`. Alternatively supply `potential=(r_m,total_V_eV)`.

The WKB nuclear inner boundary is the explicit FIXED model radius **4fm =4e-15m =4e-5Angstrom**, distinct from the7fm resonance channel radius. The total-potential table MUST cover that radius and the outer turning point (at least50eV CM support). A0.01–3Angstrom reporting-only table is rejected for missing the inner domain; a table missing the outer turning point is also rejected. No guessed short-range extrapolation is performed: interpolation extrapolation is disabled, and the action always starts at4fm, even when the first supplied sample lies at a smaller radius. Providers must supply the full physical inner domain, not merely a plotting/reporting window.

This path uses the s-wave WKB action ratio at that same fixed4fm radius with the fixed nuclear S-factor, not an inferred constant Ue. The bare pre-exponential is retained, a distinct documented potential approximation.

The result separates `physical_yield_per_C` (reactions/TRUE charged coulomb, including explicitly modeled neutrals), `counts_per_C` (efficiency-weighted counts/CUP coulomb), and `reported_signal`/`R` (possibly another current power). Sigma mb converts by1e-31m²; incident deuterons per true coulomb carry1/e. D2+ current contributes TWO deuterons. The old transport helper's arbitrary `mb*D/nm²` quantity is NOT calibrated yield/C: its conversion would be1e-13/e. A zero reference is rejected, not replaced by a fake fallback.

Parameters for `thermal` must explicitly provide `A`, `E_R_eV`, `T_track_K`, `E_A_eV`, and `Ue_hot`. The fitting convention is ER1eV,T1080K and source hot Ue(Pd/Ti/Zr)=0.400/0.370/0.340keV. A is an EFFECTIVE FITTED track/resonance product in m³s at that convention, not a predicted lifetime or absolute rate. When changing explicit ER/T, keep A fixed as the convention-bound coefficient; do not also absorb those changes into A a second time. The rate is N0² exp(-EA/kT)<sigma v> A, not a fitted per-material additive floor. It averages N_D² over the traversed lattice track, excluding the separately counted adsorbate.

`profile_thermal.py` supplies the same track-density definition to fits, using the actual exported baseline profiles and log-linear interpolation of their positive energy-dependent density moments. Exactly one A is fitted across ALL materials within a dataset; see `../czerski/README.md` for the cross-dataset result.

Default stopping supports ICRU Pd/Ti/H. For oxide/alloy effective hosts the caller must supply a stopping function; unsupported Z is NOT silently guessed. `consumer_example.py` supplies the original transport/BCA LS/ZBL kernel explicitly, including the real Ti oxide Z10.7. It exercises all16 supplied profiles, every confirmed candidate, and both supplied Pd TF-potential endpoints. Inputs are the exact0aab89af99e2bad8ace81c4b413c033cf87812e2 exported baseline bytes (`inputs/profiles_baseline.json`, source SHAf06a1df5fcc98e204e621080e038569a171a19da7fe0cf900de8b3f5291aa51d), not a transport replay. Screening inputs are labelled fitted/TF MODEL estimates supplied for use until computed screening arrives, not new DFT.

The final actual-profile command exited0. Its uniform-target integral limit differs by9.19e-5 relative, the surface-sheet unit identity agrees to the reported binary64 resolution, detector/current changes leave physical yield/C unchanged, and a zero reference is rejected. A separate coupled supplied-potential/surface-layer call exercises the extended WKB energy domain without replaying unaffected cases; its output is in `consumer_log.txt`. For supplied Pd TF lambda0.57Angstrom plus a3keV surface shift, E_cm0.5/6.5keV gives R0.0699026/1 and2.06924e5/2.96017e6 physical reactions/C.

The targeted `--potential-domain` command exited0 in26.83s on the actual saved-profile API. Both the0.01–3Angstrom reporting-only table and an incomplete outer table were rejected with explicit domain errors. Full-domain recovery succeeded; extending the supplied samples from4fm down to1fm while retaining the physical4fm boundary changed yield/C by only2.69e-5 relative (interpolation-grid difference). The domain evidence is retained alongside the unaffected checks in `consumer_verification.json`.

### Pass5 transported-profile fit acceleration

The pass5 merit consumer reuses this same kernel. Stopping is now solved once
per immutable seven-field profile and incident energy, independently of the
nuclear/screening parameters. `csda_nodes` returns 12-point Gauss reaction nodes
on segments no longer than10nm; a fitted surface boundary is split exactly.
Both production yields and the compiled fit responses consume those nodes.
The previous coupled stopping/yield ODE is not retained as a second production
integrator. `parameters.additional_Ue_keV` is an explicit scalar or `(E_cm,z_m)`
addition after the selected cause, used for the tenfold-screening sensitivity.

The actual PdD profile pilot agrees with the previous ODE within6.29e-6
relative over the no-plateau, surface, screening-rise and neutral cases
(`../merit/path_parity.json`). This is numerical parity, not material-input
uncertainty. Transported cause refits and fixed-PdD-reference comparisons live
in `../merit/`, separately from the retained uniform pass4 fit artifacts.

## Reproduce and consume

From the repository root with NumPy/SciPy/Matplotlib installed:

```
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python g4_beam/forward/second_route.py
python g4_beam/forward/fits.py
python g4_beam/forward/second_route.py --common
python g4_beam/forward/analysis.py
python g4_beam/czerski/fit.py
python g4_beam/forward/consumer_example.py
python g4_beam/forward/consumer_example.py --potential-layer
python g4_beam/forward/consumer_example.py --potential-domain
```

The original Eq7 second route is already computed and need not be replayed for a changed common reaction. `fits.csv/json`, `bulk_screening.csv`, `abs_yield.csv`, `consumer_outputs.csv` and `consumer_verification.json` are the numerical consumers. `fig_fits.png` and `replot.csv` show only the confirmed mechanisms. `fits_log.txt` retains the observed final-serialization failure; its NumPy boolean was fixed and the complete numeric bulk CSV was reserialized without repeating any fit. Analysis and consumer commands run separately afterward. Recorded scientific runtime: Python3.12.3, NumPy2.5.3, SciPy1.18.1, Matplotlib3.11.2, explicit single-thread BLAS/OpenMP.

The original range table is an angular-first-moment approximation, not a replacement for the supplier's BCA ranges. Subatomic fitted layer thicknesses remain a screening/amplitude degeneracy, not a literal resolved ALD layer. The three principal omitted physical effects are target evolution, angular/straggling corrections to the supplied-profile yield path, and an independently measured neutral-energy/charge-collection distribution. No target-time reconstruction, new DFT, old campaign replay or embargoed data download is performed.
