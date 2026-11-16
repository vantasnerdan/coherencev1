# Primary Czerski / EXFOR fits and shared-material result

## Data and observable

`data/README.md` binds the exact [arXiv:2409.02112v3](https://arxiv.org/pdf/2409.02112v3) source, digitization and [IAEA EXFOR O1084003](https://github.com/IAEA-NDS/exfor_master/blob/2193cdb7066421fb7931b696d6bab219f9cc8b03/exforall/O10/O1084.x4). All 162 extracted vector marker centres remain in `figure_points.csv`; the complete 19-point Ti EXFOR table is preserved. The primary-source extraction owner transferred these files before integration.

Czerski yields are detected counts/(microampere kilosecond), i.e. 1000 detected counts/C. They are not total emitted reactions/C without detector acceptance. Energies are already per deuteron; E_cm=E_lab/2, with NO second division by D2/D3 multiplicity. EXFOR is normalized at **6.45 keV LAB = 3.225 keV CM**, not 6.5 keV CM. No Karahadian cup/current hypothesis is applied to these data.

Panels2a/3a/4a/7a are yields. Panel3b is a dependent enhancement ratio; panels2b/4b/7b are proton-peak CHANNELS, not extra independent yields. `objective_membership.csv` records exactly which source rows enter the yield objective. Fig7's overlapping high-energy D1 points are not counted again; its distinct below3keV CM extension is retained. All original points remain available. `enhancement_comparison.csv` and `spectral_comparison.csv` model the dependent observables without duplicating their likelihood.

The magnetically selected beam still has source-described neutral components before deceleration: Zr D2 uses approximately10keV lab neutrals, Fig7 D1 approximately20keV, and D3 approximately4keV. These remain separate beam states. Fitted Zr neutral fractions are measurement nuisance parameters common to every candidate, starting from the stated ~3e-4 estimate. They are not new material candidates. Pd had no deceleration and does not acquire an invented above-nominal neutral floor. Ti D3 residual neutral contamination is unspecified and neglected in the clean-series comparison, not measured zero.

Statistical error bars in the Czerski figures were NOT digitized. EXFOR has four stated statistical percentages and a common8% compiler digitization resolution, with unspecified normalization correlations. Therefore the primary/shared objective is an explicitly labelled balanced log-residual score, **not an experimental chi2, AIC, significance or confidence interval**. For each dataset its sum of squared log residuals is divided by its point count, so datasets have equal objective weight. `rms_factor = exp(sqrt(mean(log(prediction/observation)^2)))` is a descriptive mismatch factor. Only the Karahadian columns also report that dataset's existing digitized-error chi2 convention.

## Shared implementation and normalization

The physical candidate set is exactly layer (wide bounds), sigmoid_tt, fixed-energy fast neutrals, thermal spike, plus no-additional-plateau control. The COMMON FIXED proton nuclear reaction background is Czerski v3 Eqs14/17 with40meV,115degrees,7fm; it is not a fitted-width plateau candidate. It is applied to every source and Karahadian material. The neutron branch remains Bosch–Hale, with its unspecified resonance handled only in the separately requested sensitivity in `../forward/bulk_screening.csv`.

Direct yields use the same numerical thick-target machinery as `forward/model.py`. For these sources `SourceTarget` uses the authors' electronic sqrt(E) stopping approximation. The effective absolute coefficient is calibrated with the Table2 M values in their high-energy asymptotic convention, using the stated detector geometry. `normalization.csv` exposes that coefficient, M, acceptance and the independent Table1-composition density used in the thermal rate. M is a fitted approximate-curve normalizer, NOT an independently measured density or efficiency. Retained Bosch–Hale energy dependence and the full Eq7 integral need not give the same fitted cold Ue as the authors' approximate Eq23; fitted values here do not replace their measured-fit report of25±10eV for cold Ti. The source's370eV hot Ti value is an inferred thermal-model input, not that cold measurement.

Layer depth uses the same angular-first-moment definition as the original forward fits. Very small layer widths are unresolved screening/amplitude degeneracies, not literally identified physical layer thicknesses. EXFOR supplies no D-depth measurement; its normalized direct comparison uses the source Ti stopping convention. Thermal EXFOR comparisons condition on that same homogeneous Ti density, rather than inventing a measured EXFOR density. The callable can replace this assumption when a profile exists.

## Thermal candidate: one fitted A, not a free floor per material

The implementation is the actual source Eqs7–9/21:

```
N_mobile = N0 exp[-EA/(kB T)]
rate = N0 * N_mobile * <sigma v>
yield_per_incident_D = rate * A
A = effective V*tau*(Gamma_p/Gamma), in m^3 s, fitted at ER=1eV convention
```

The narrow-resonance reactivity, including sqrt(ER), exp(-ER/kBT), screening penetration and the7fm radius, is implemented in SI in `../forward/thermal.py`. Source hot Ue is FIXED at Zr0.340, Ti0.370, Pd0.400keV; T_track1080K and ER1eV define the fit convention. EA is Zr0.413, Ti0.680 and Pd0.150eV. Density enters as **N0²**, not N0 alone and not an all-mobile N²/2. Czerski density comes from Table1's stated hydride compositions and mass densities; Table2 M is never used as thermal density. Karahadian density moments come from the ACTUAL supplied end-of-run transport profiles, including distinct Pd/PdD/Ti/TiD conditions. The common straight-ahead track-average N_D² and energy interpolation are declared in `../forward/profile_thermal.py`.

Exactly ONE A is fitted across ALL materials within each dataset. No per-material amplitude or hot-screening fit remains. Separate-dataset and single-universal-A fits use the same fixed hot screening, barriers and density convention. A absorbs unresolved track/branching and the ER convention; it is fitted, NOT an absolute-rate prediction, lifetime or microscopic identification. In the callable A, ER, T_track and EA are explicit. Vary ER/T through the reactivity once, keeping the convention-bound fitted A fixed rather than rescaling it a second time.

## Deciding results

Individual primary-dataset descriptive RMS mismatch factors (thermal is constrained across all three materials, not independently fitted per row):

| Candidate | Czerski Zr | Czerski Ti | Czerski Pd | EXFOR Ti alone |
|---|---:|---:|---:|---:|
| No additional plateau |2.416|2.757|2.271|1.120|
| Wide layer |1.301|1.187|1.186|1.120|
| sigmoid_tt |1.336|1.175|1.186|1.120|
| Fast neutrals, source beam constraints |2.416|2.757|2.271|1.120|
| Thermal, shared A |1.318|1.315|1.926|1.120|

The EXFOR no-plateau case drives additional plateau strength toward zero (thermal A~1.0e-47m³s, compared with Czerski A=1.3266e-34m³s). Layer and sigmoid have effectively bare no-plateau limits there. With shared Ti physics between Czerski and EXFOR, mismatch factors are layer1.561/1.973, sigmoid1.676/2.448, and no-additional-plateau/neutral3.120/2.932. A single thermal A across all Czerski materials plus EXFOR gives aggregate Czerski/EXFOR factors1.932/1.635, with cold Ti Ue0.02699keV. That cold value must not hide the other datasets' remaining mismatch.

Pd can share the layer or sigmoid parameters between Czerski and current-corrected Karahadian with factors1.293/1.369 or1.226/1.314, respectively. The corresponding Ti shared fits remain poor once the zero-plateau EXFOR data are included: e.g. sigmoid gives Czerski/EXFOR/Karahadian1.347/3.273/3.296. These are conditional model residuals, not statistical exclusion of mechanisms or of target-condition effects.

The cross-material thermal result is decisive under the supplied frozen-density and track assumptions:

| Karahadian normalization | A_Czerski (m³s) | A_Karahadian (m³s) | A_K/A_C | Separate-A RMS factors C/K | One-A RMS factors C/K |
|---|---:|---:|---:|---:|---:|
| As plotted M0 |1.8019e-34|1.1216e-31|622.5|1.544 /2.261|3.188 /4.465|
| One cup-current power M1 |1.7708e-34|5.0263e-32|283.8|1.535 /1.970|2.660 /3.448|

These rows share cold material parameters in the joint fit and use one amplitude per DATASET, not one per material. A universal A does not reproduce both curve sets under this model. Even separate amplitudes do not make the fixed-hot-screening, common-track description a precise simultaneous fit of all materials. Allowing unrelated per-material floors would conceal rather than resolve this result. This is not a universal physical no-go: the track approximation, actual D distributions and source absolute normalization remain physical assumptions, and EXFOR's density is not measured here.

`thermal_shared.json` contains the full fitted products, objective and grouped residuals. `cross_dataset.csv` is the complete per-material/per-dataset table with each `global_scope` explicit. The best Karahadian-only thermal fit under the six existing measurement hypotheses has chi2=132.23 (M4), compared with sigmoid36.92, layer40.29 and fitted-energy neutrals46.31; fixed5keV gives49.65 and the common-reaction no-additional-plateau control181.45. These are the same48 points and same existing Karahadian error convention. The fixed nuclear reaction term itself is now in every row.

The dependent proton-channel check is useful but not an independent extra likelihood. With calibration fitted only above3keV CM, low-energy Pd channel RMS is2.75 for the constrained thermal model versus5.61 for a direct layer; Zr D3 improves from1.89 to1.04channel, whereas Fig7's neutral-dominated data need their neutral component (thermal1.35 versus direct/neutral0.89channel). This distinguishes component kinematics without claiming an established causal identification from yield curves alone.

## Reproduce

From the repository root:

```
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python g4_beam/czerski/fit.py
# Re-render existing fitted outputs without redoing fits:
python g4_beam/czerski/fit.py --plot-only
```

Final corrected source fit exited0 in55.61s, single numerical thread. `primary_fits.png` and `spectral_check.png` were rendered and inspected; `replot.csv`, the original marker coordinates, enhancement and spectral tables retain the numerical plotting inputs. `forward/fig_fits.png` is the third delivered figure. No PDF, author-request raw data, invented error bars, new plateau candidate, per-run reconstruction, DFT job or old campaign is included.
