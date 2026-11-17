# transport: near-surface deuterium with beam damage

`transport.py` is a 1D finite-volume model over the full 250 um foil. Diffusion follows the gradient of the D chemical potential on Pd-D (alpha/plateau/beta) and Ti-D (alpha/plateau/delta) isotherms, so hydride fronts move. It carries vacancy traps (Pd: 0.23 eV, 6 D per vacancy, 1% saturation), vacancy diffusion and supply, surface desorption against 0.037 Pa chamber D2, sputter erosion, a chemisorbed D monolayer, an optional overlayer (sputtered, regrown, with its own permeance and D capacity), a back-side fugacity limited by electrochemical supply, and beam heating. Knobs are arguments to `Model(case)`: `overlayer`, `d_ol`, `j_rel`, `dT`, `strain`, `alloy`, `s_front`, `vsat`, `no_gas`, `eta_abs`, `v_init`.

`bca.py` is a binary-collision Monte Carlo (ZBL potential, scaled Lindhard-Scharff stopping, NRT vacancy counting). The supplied cache is required by the updated runner: no automatic multi-core table rebuild. Historical selected RustBCA comparisons are retained in `lit/`; they are not a new all-target4% validation. The new sourced LSS/ZBL calculation is an infinite-medium projected moment, not the cache's half-space implanted-only mean; Ti metal and TiD2 must not be conflated. `run_all.py` and `post.py` define the original metrics; `export_profiles.py` writes `baseline/profiles_full.json`, the full end-of-run seven-field profiles that the forward model reads; `checks/run_updated_knobs.py` retains the original estimate rows and recomputes only changed Pd rows.

Results:
- Independent integration of all16 original SI profiles reproduces all48 five/twenty/seventy-nm D/metal windows exactly. Bare adaptive-CSDA TiD/Ti yield ratios are0.9919–1.0057 over the four supplied energies: beam-loaded Ti and TiD both approach D/Ti~1.84 in the70nm window.
- Closed-form dilute Pd Sieverts gives D/Pd3.22021e-5 at300K and0.03666355Pa chamber D2, under the supplied empirical isotherm.
- Pd250um constant-D **dilute permeation time lag** is176.981s; an independent diffusion-eigenmode moment agrees. This is not nonlinear hydride-front equilibration or a promise of full back-side loading within minutes.
- For2nm **nonregrowing** TiO2/Al2O3, the tested particle/atom balances give95% removal0.472–1.633h. Assumed water regrowth pins native TiO2 at~0.92–1.23nm and prevents95% clearance in these cases. These are model sputter-inventory times, not measured permeability lifetimes.

Drops: 3D microstructure, blisters, ballistic mixing, vacancy clustering, run history, electrochemical kinetics, surface chemistry beyond one sticking coefficient, angular scattering in the yield integral.

Three biggest omissions: beam-side surface chemistry (carbon, sulfur, real adsorbed D), microstructure of the 200 nm sputtered film, run history and lateral beam profile.

The model retains Eads0.45eV, Eb0.23eV, mtrap6 and sticking0.5 by the lead's current-results decision; uncontrolled electronic endpoints and unfinished paths are not parameter replacements. Other estimates include vacancy saturation1%, oxide permeances/capacities, most alloy transport kinetics, Ti vacancy migration0.5eV and electrochemical absorption5%. The controlled conditional rerun, recombination barrier/TST rate and C/O cases were not computed.

## Vacancy outcome

Only n=0/1 were tested. The8-Pd pilot gives first-D electronic binding+0.130945eV relative to a matched bulk interstitial. The last27-Pd trial gives0.242437eV arithmetic but remains unrelaxed (force0.169324eV/A versus0.08 target); neither is image/k/cutoff qualified. They are display-only comparisons beside0.23eV, not finite-temperature parameter replacements or evidence of six-D capacity. The private immutable original retains all17 findings/dispositions and actual failure/source/checkpoint records.

## Clean Pd(111): current-results stop

| Quantity and sign convention | eV/D | Status |
| --- | ---: | --- |
| Retained positive Langmuir Eads | +0.45 | Finite-temperature model estimate |
| Computed1ML mean adsorption energy relative to half D2 | -0.465699 | Uncontrolled0K electronic endpoint |

The electronic row uses the shared H/D Born–Oppenheimer surface, not isotope-specific ZPE or a gas free energy. Its opposite sign gives a positive electronic binding magnitude; it is not a coverage-differential finite-temperature Langmuir parameter. Inputs are PBE/Pd10,300eV,k6x3x1, a fixed3.890A three-layer1x2 slab with clamped bottom atoms. The k/slab/vacuum/cutoff and lower-coverage controls were not run.

The28 saved CI-NEB bands end at step27 with fmax0.514006eV/A versus0.05 target. User-directed termination is not numerical completion. No converged barrier, harmonic prefactor, TST rate, controlled conditional transport rerun or C/O result is claimed. Raw unfinished path energies are not barriers. See `surface/surface-public-row.json`, `surface/surface-public-inputs.json` and the complete18-finding `surface/surface-findings.json`.

## Updated data/parameter contracts

Profiles contain seven-field SI layers `(z_top,z_bottom,n_host,n_D,Z,M,n_metal)`, adsorbed D/m² and state. Updated case JSON saves the exact end-state `profile_SI` for the forward consumer, not the driver's sampled four-column plotting profile. `yield_per_ion` returns raw mb*D/nm²; multiply by1e-13 for probability. `yield_per_charge` is raw yield per charged particle; the original runner already multiplies by1e-13/e for physical DDp reactions per true incident C. Detector efficiency and Faraday-cup bias are separate.

Electronic differential adsorption may enter an explicitly conditional0K scenario only with matched half-H2 reference and missing ZPE/vibrational/configurational/gas corrections labelled. A barrier/prefactor never replaces sticking. The updated case validator refuses wrong binding/adsorption references and automatic sticking/capacity overrides. Original estimates remain a distinct scenario.

## Pd-Rh loading-side estimate

[Losiewicz et al.(2020), DOI10.3390/ma13010162](https://doi.org/10.3390/ma13010162) reports Pd80Rh20 H-electrosorption at298.15K: H/metal0.002/0.4/0.8 at RHE potentials+0.050/−0.033/−0.150V. A conditional reversible-Nernst fit supplies a300 = 0.0472524, b300 = 0.602170 and dH = −21.3069kJ/mol while retaining the specified Pd entropy/kinetic/surface assumptions. Three fitted parameters/three source points provide no holdout validation. H→D and temperature extrapolation are assumptions; the source is not direct gas-pressure or DFT data.

This20%Rh foil is not the ordered25%Rh screening cell; no causal d-band derivative is inferred from combining them. All four estimate-only Pd/PdD by0.25/2.5keV alloy rows and exact SI profiles are saved. At Ecm0.25keV, the Pd-Rh row gives top5nm D/metal0.00094873 (0.3053×Pd baseline) but fixed-fitted-Ue yield0.9587×baseline because the monolayer remains nearly saturated. These rows do not consume the unqualified vacancy trials.

The second-route range source is [Lindhard/Scharff/Schiott1963, pp24–25](https://gymarkiv.sdu.dk/MFM/kdvs/mfm%2030-39/mfm-33-14.pdf). Its infinite-homogeneous-medium boundary condition is explicit in the published checks; it does not replace the conditional implantation cache.

## Reproduce the baseline alloy rows

Run serially on one allocated core with NumPy/SciPy installed. The saved BCA cache is required and included; no table rebuild or Rust toolchain is used. Outputs must name a fresh directory. The shared `g4_beam/screen/screen_io.py` comes from the preceding screening commit.

```sh
python g4_beam/transport/checks/run_updated_knobs.py \
  --reference-runner g4_beam/transport/checks/reference/run_all.py \
  --original-table g4_beam/transport/checks/reference/knobs.csv \
  --parameter-sets g4_beam/transport/checks/inputs/baseline-parameter-scenario.json \
  --pdrh g4_beam/transport/checks/results/Pd80Rh20-conditional-isotherm.json \
  --output pd-rh-baseline-reproduction --cpu-budget-s 120
```
