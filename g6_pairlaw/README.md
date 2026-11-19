# Pair-at-rest response, pass 6

This replaces a two-anchor interpolation with normalized radial contact calculations for (a) a scaled molecular BO **model** and (b) a screened repulsion in a harmonic confinement. The material levers are screening length, actual mean pair separation, and confinement stiffness. These three inputs are independent: a small screening length alone does not hold two nuclei close together. Part (c) has not run; no material total pair potential has been supplied.

## Model and rate definition

The calculation keeps relative nuclear kinetic energy, a real central potential, the normalized radial ground state, the Coulomb cusp, and the nuclear short-range coefficient. It drops anisotropic/many-body lattice dynamics, electronic nonadiabaticity, and resolved nuclear spin/channel dynamics. It is a conditional static pair response, not a claim about a prepared lattice or power gain.

For isotope reduced mass mu, solve `-t u'' + V u = E u`, `t=hbar^2/(2 mu)`, `u(0)=u(infinity)=0`, `integral |u|^2 dr=1`. Length is Angstrom and energy eV. The three-dimensional s-wave is `psi=u/(sqrt(4 pi) r)`, so `rho_contact=|u'(0)|^2/(4 pi)` in Angstrom^-3. Each rate is `A rho_contact` per normalized pair per second; no site count, encounter count, or detector efficiency enters it.

For distinguishable incoming plane waves, `sigma=S(E)/E exp(-2 pi eta)`, `C_eta^2=2 pi eta/(exp(2 pi eta)-1)`, and

`lim_(E->0) sigma v / C_eta^2 = S(0)/(pi alpha mu c)`.

This matches the Coulomb prefactor, rather than multiplying a Gamow exponent by an arbitrary vibration frequency. The code converts keV barn to J m^2 and m^3/s to Angstrom^3/s explicitly. For DD the measured unpolarized cross section includes the identical-boson factor 2 and spin weights 1/9 and 5/9 for s-wave S=0,2. A **normalized allowed-spin statistical mixture** has weights 1/6 and 5/6; its A is therefore 3/4 of the beam expression. For a nine-state initially random spin ensemble with its odd-spin part excluded from this s-wave contribution, the coefficient instead is 1/2 of the beam expression (2/3 of reported DD rates). Individual polarized S=0 and S=2 rates cannot be recovered from the summed Bosch-Hale factors. Do not apply an extra identical-pair factor of 1/2 to a rate already per pair; that factor belongs in counting unordered pairs in a population.

DD includes both `D(d,p)T` and `D(d,n)3He`: S(0)=55.576 and 53.701 keV barn, respectively, from the existing Bosch-Hale branch polynomials in `g4_beam/forward/model.py`. Those branch fractions are retained only as low-energy background effective coefficients, not a full spin-resolved nuclear model. DT uses the measured-data fit quoted in [Kamimura, Kino and Yamashita, arXiv:2112.08399v3, Fig.4](https://arxiv.org/html/2112.08399v3):

`S_DT(E) = (26 - 0.361 E + 248 E^2)/(1+((E-0.0479)/0.0392)^2) MeV barn`, E in MeV.

Only its smooth E=0 limit is consumed. This avoids extrapolating the Bosch-Hale DT rational fit outside its low-energy domain. The entrance is predominantly nuclear spin 3/2, with unpolarized weight 2/3; reported DT maps use the **unpolarized** coefficient. The pure quartet coefficient is 3/2 times larger, and is used for the state-matched muonic comparison. DT's established broad positive-energy 5He structure is already in the fitted S; a new DD threshold resonance is not assigned to DT.

## Potentials and the meaning of distance

(a) `V_m(r)=m E_h [E_BO(m r/a0)+1]`: full two-electron ground singlet BO curve, both length and energy rescaled; the nuclear reduced mass is not rescaled. All 82 entries of [Pachucki, arXiv:1007.0322v2, Table 4](https://arxiv.org/html/1007.0322v2) are in `source-data-if-needed/pachucki-table4.csv`. Energies include nuclear repulsion; subtract `1/R` before interpolating the smooth electronic part. At R<0.1 bohr use `E_BO=1/R+E_He+a R^2`, with the infinite-mass nonrelativistic helium limit -2.9037243770341196 hartree and a fixed by the first tabulated point. A second-point core and source-table thinning are explicit controls. Above 20 bohr an unused van der Waals continuation exists; every consumed molecular calculation ends at 15 bohr and has a negligible outer tail. The input table's 10^-15 energy precision is NOT claimed for interpolation or rate accuracy.

The electronic length `a0/m` is a scale label, not a fitted Yukawa screening length. Although m=206.768262 uses a muon/electron mass ratio as an illustrative endpoint, this is still a scaled **two-electron** Hamiltonian. A dd-mu or dt-mu molecule contains **one** muon and has nuclear-muon mass-polarization/rearrangement dynamics. It is not this endpoint. No coefficient is fitted to force agreement.

(b) `V(r)=C exp(-r/lambda_s)/r + k (r-r_center)^2/2`, `C=e^2/(4 pi epsilon0)`; k in eV/Angstrom^2. The radial s-wave probability determines `mean_A=<r>` and `rms_A=sqrt(<r^2>)`. `center_A` is the centre of the *external harmonic term*, not the potential minimum or achieved separation. Repulsion and zero-point motion push the achieved mean outward; even a zero centre generally has a finite minimum mean. A threshold not reached at centre>=0 is recorded as `no_root_nonnegative_center_rate_too_low`, never a negative/unphysical centre or a fabricated separation.

**Ratio 7 is lambda_s / a_lattice, not lambda_s / <r>.** Here the independent reference lattice-site spacing is explicitly 2.78 Angstrom (the existing closest relaxed D-D reference); the caller must replace this denominator with its actual physical lattice spacing. The scaled-family column `ratio7_electronic_scale` uses a0/m as a model proxy only, not an inferred material screening length. This distinction prevents a compressed pair from redefining the material's lattice spacing.

## Muonic comparisons, without forbidden anchor calibration

For like-state controls the script uses the full nonadiabatic Coulomb three-body nuclear contact, not a fitted muonic rate: DD ground L=0 contact 2.43871205e-6 in muon atomic units from [Frolov, arXiv:1111.7056v3](https://arxiv.org/html/1111.7056v3), and DT ground L=0 contact 8.870847(1)e-7 from [Frolov and Wardlaw, arXiv:1008.3010v4, Table 8](https://arxiv.org/html/1008.3010v4). Convert by `(a0/206.768262)^-3`. These are source-calculated contacts, not contacts computed by this radial script. Their point-Coulomb, nuclear-force-free meaning is preserved. The DT quartet contact approximation is compared with the 1.15e12/s ground J=v=0 coupled-channel rate of Kamimura et al.; disagreement measures short-range nuclear and factorization model error as well as coefficient approximation.

The commonly quoted DD-mu value is not a valid s-wave anchor. [Wu, Cui and Kamimura, arXiv:2508.12783v2](https://arxiv.org/html/2508.12783v2), Eq.22 and Tables 3-4, identifies the measured `(J,v)=(1,1)` effective rate `(3.81+-0.15)e8/s` as `lambda_11+Gamma_dex`, with Gamma_dex about 0.2e8/s. Fusion is **p-wave**, nuclear S=1. Its roughly 3.61e8/s fusion part needs a derivative-contact/partial-wave nuclear operator; `|psi(0)|^2=0` for an isolated l=1 relative state. Source p-wave calculations give lambda_11 from 1.75e8 to 5.05e8/s across the reported nuclear inputs (not an s-wave uncertainty interval). `source_comparison.json` separately reports the scaled-BO/L=0 contact mismatch and any arithmetical comparison to the nonlike experimental rate, explicitly labelled. Muons here are diagnostic literature controls, not proposed material ingredients.

## Resonance: what can and cannot transfer

The existing `reaction_S` carries the literal positive-energy scattering Eqs.14/17 convention of [Czerski et al., arXiv:2409.02112v3](https://arxiv.org/pdf/2409.02112v3), proton width 40 meV, phase 115 degrees, channel radius 7 fm, ER->0; neutrons remain Bosch-Hale. `on_shell_resonance_control.csv` evaluates that exact expression at **positive** CM energies. Its E^-2 resonant term does not define a finite constant S(0). Inserting the bound eigenenergy (whose zero can be shifted) into it, or substituting a vibration energy without a matching nuclear model, is unjustified. A scattering fit/phase alone does not identify a bound-state matrix element, the necessary energy dependence, and its nuclear boundary conditions. The physical inner point-contact convention here is also distinct from the forward consumer's fixed 4 fm WKB boundary and the resonance's 7 fm channel radius.

Accordingly there is no claimed physical DD/DT resonance-transfer curve or measured relaxation factor. Threshold tables include background F=1 and **sensitivity controls** F=1e3,1e6 multiplying the constant contact coefficient. These are chosen exploratory multipliers, not bounds on the unknown resonance and not confidence intervals. They answer exactly how much extra short-range strength would relax a target if a future matched nuclear model supplied it. Their finite ranges cannot enclose an unidentified amplitude. No mirrored neutron width, no DD width transferred to DT, and no forbidden anchor tuning is used.

## Numerical method, consumer and cost

Set x=ln(r), u=sqrt(r) phi. The generalized eigenproblem is `-t phi_xx+t phi/4+r^2 V phi=E r^2 phi`; whitening by r gives a symmetric tridiagonal matrix on a uniform log grid. The well eigenvector is normalized with `integral u^2 r dx=1`. For a narrow remote well the eigensolve starts 14 local oscillator widths to its left, preserving well resolution instead of undersampling it to include a remote origin. The full origin-to-well interval is always propagated separately. Tiny eigensolver tails are numerical noise: the origin amplitude is obtained by the regular outward Riccati system `y'=Q-y^2`, `z'=y`, with y=phi'/phi and z=ln(phi), using LSODA's automatic stiff integration on x. It is matched at a resolved inner-forbidden-region amplitude. The origin starts with `u=r[1+C r/(2t)+O(r^2)]`; no `u(r_min)/r_min` is mistaken for the normalized contact. Rates remain in log10 space even below underflow. The first fixed-step Numerov attempt failed at a very stiff high barrier (denominator pole); the stable Riccati representation repairs that failure rather than suppressing it.

WKB is a second barrier-propagation route matched to the **same direct well solution**, not an independent eigenvalue calculation. In log coordinates Q=1/4+r^2(V-E)/t, it propagates `phi ~ Q^-1/4 exp(integral sqrt(Q) dx)` from the Coulomb-origin limit Q->1/4. `wkb_minus_direct_dex` reports its approximation error rather than tuning a prefactor. The direct rate is the consumer. Two grid sizes (6500,11000), core/source controls, changed matching height, and exterior extension on representative extremes price the deciding numerical error; they are not statistical intervals or rigorous enclosures.

`harmonic_map.csv` supplies both isotope responses over lambda_s=0.01,0.03,0.1,0.4,0.54 Angstrom; k=10,100,1e4,1e6 eV/Angstrom^2; and centres=0,0.02,0.05,0.1,0.25,0.75,2.78 Angstrom. Very stiff compact cases are lever controls, not measured material wells. `harmonic_thresholds.csv` holds all four requested contours, with no-root status and minimum achieved mean. `scaled_thresholds.csv` holds the analogous BO-family targets. The reusable call is `harmonic(lambda_s, r_center, k, isotope)`; import `script/pairlaw.py` by path. Its arguments cannot silently mean a fixed mean separation.

`scaled_fit.json` gives a closed form in x=2 ln(m)/ln(206.768262)-1, degree-5 Chebyshev, with 16 interleaved held-out direct solutions. `harmonic_fit.json` gives degree-12 closed forms in x=2*r_center/0.75-1 at each *fixed* isotope/lambda/k. Coefficients multiply T_j(x). Each has 12 independent interleaved holdouts, an observed error, and a `usable_max_0p15dex` flag. Consume only flagged slices, only within the stated domain, and never extrapolate/interpolate across lambda/k without an added error check. These finite holdout errors are not uniform proved bounds. Thresholds are solved from the direct radial callable, not from an unchecked fit.

Reproduce the settled scientific/consumer run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python g6_pairlaw/script/pairlaw.py > g6_pairlaw/output/run.log
```

`summary.json` records actual runtime/environment, source-table hash, and the reference material-facing DD/DT calls. Assertions check the log-grid tolerance, achieved direct contour rates, unchanged ratio-7 denominator, and the distinction between centre and achieved mean. No external potential is downloaded or invented at runtime. No build/lint/formatter or unrelated campaign is needed.

## Inputs ledger

| Input | Source/settings | Status and reader |
|---|---|---|
| BO potential | Pachucki Table 4, all 82 samples, hartree/bohr | Source precision high; interpolation/core controls in numerical output; radial model |
| DD branch S(0) | Existing Bosch-Hale p/n coefficients | Effective unpolarized beam intercepts; spin/exchange conversion above; DD coefficient |
| DT S(0) | Kamimura Fig.4 measured-data fit, MeV/barn | Smooth low-energy extrapolation, not a new nuclear fit; DT coefficient |
| muonic nuclear contacts | Frolov/Frolov-Wardlaw full Coulomb L=0 states | Literature values only, correctly converted muon units; mismatch controls |
| DD threshold term | Existing forward reaction_S, exact positive-energy convention | Not identified for bound states; on-shell control only |
| Harmonic k, centre, screening | Declared independent parameter grid | Conditional model inputs, not converged first-principles material values |
| Physical lattice-site spacing | Independent explicit reference 2.78 Angstrom | Never pair mean; ratio-7 label only |
| Material total pair potential | None delivered | Part (c) not run |

## Settled result

**Material-facing target:** in this model, 1e-9/s needs either lambda_s=0.01 Angstrom with k=10 eV/Angstrom^2 and achieved mean DD/DT separation 0.522/0.523 Angstrom (ratio 7=0.00360), or an extremely stiff k=1e6 eV/Angstrom^2 well holding 0.0451/0.0446 Angstrom at lambda_s=0.54 Angstrom (ratio 7=0.194). These are two different confinement routes, not a universal separation law. A physical resonance relaxation is not identified; an assumed millionfold constant-contact gain relaxes the scaled-family 1e-9/s mean from 0.0781/0.0768 to 0.1178/0.1125 Angstrom (DD/DT), a sensitivity control only.

The final command exited 0 in **98.09 s wall time** (97.12 s inside the script), single-thread NumPy 2.5.4/SciPy 1.18.1/Python 3.12.3. One earlier 10.06 s run exposed the fixed-Numerov instability described above; it did not produce a completed material result. The retained `output/run.log` is the settled successful run.

### Source comparisons and numerical scale

| Quantity | Resolved result | Meaning |
|---|---:|---|
| D2, m=1, DD ground s-wave background | 1.510e-64/s | 0.503 of 3e-64/s; reproduced within factor 2.0 without fitting |
| Scaled BO m=206.768, DD background | 2.918e11/s | Not the observed J=1 muonic state |
| Scaled BO / true DD-mu ground L=0 contact | 17.55 | Hamiltonian/model mismatch, not a calibration parameter |
| True DD-mu ground L=0 contact times DD coefficient | 1.663e10/s | Conditional s-wave contact estimate, not a prediction for measured p-wave DD-mu |
| DT-mu ground quartet, literature contact times A | 9.625e11/s | 0.837 of the state-matched 1.15e12/s coupled-channel result; agreement within factor 1.20 |
| Scaled BO / true DT-mu ground L=0 contact | 20.91 | Two-electron endpoint again overcompresses the contact |
| Scaled BO DT quartet / coupled-channel rate | 17.50 | Endpoint does not reproduce the real one-muon molecule |

The scaled DD endpoint is arithmetically 808 times the J=1 measured fusion part, but this is **not a like-state validation** and is not used as one. The decisive conclusion is to retire the two-anchor physical interpolation, not tune its coefficient until two different Hamiltonians/states line up.

Across the nine representative molecular/material extremes, changing 6500 to 11000 log-grid points shifts the direct rate by at most **4.81e-5 dex** (0.0111%). BO source thinning changes it by at most 0.00238 dex, and the united-atom core control by 0.00155 dex. Matching-height shifts are below 3.91e-6 dex; the largest outer-domain shift is 1.27e-5 dex, including a very stiff remote well. WKB/direct differences are **0.000269 to 0.10775 dex** in those cases, at worst a factor 1.282; WKB therefore agrees within 30% on the explicitly checked extremes. These are observed numerical differences, not a rigorous total uncertainty. Model/state identification and nuclear factorization dominate the usable physical precision.

### Scaled-family targets

Baseline, no new threshold-resonance transfer. Ratio 7 below is explicitly the **electronic-scale proxy**, with the independent physical denominator 2.78 Angstrom.

| Isotope | Target /s | m | Achieved mean /Angstrom | a0/m /Angstrom | Ratio-7 proxy |
|---|---:|---:|---:|---:|---:|
| DD | 3e-13 | 7.891 | 0.10036 | 0.06706 | 0.02412 |
| DD | 1e-9 | 10.236 | 0.07808 | 0.05170 | 0.01860 |
| DD | 4e-9 | 10.739 | 0.07456 | 0.04928 | 0.01773 |
| DD | 1 | 24.287 | 0.03427 | 0.02179 | 0.00784 |
| DT | 3e-13 | 8.134 | 0.09690 | 0.06506 | 0.02340 |
| DT | 1e-9 | 10.348 | 0.07676 | 0.05114 | 0.01839 |
| DT | 4e-9 | 10.815 | 0.07356 | 0.04893 | 0.01760 |
| DT | 1 | 22.621 | 0.03629 | 0.02339 | 0.00841 |

Thus the earlier 0.04–0.05 Angstrom pair target is not recovered as a unique 1e-9/s mean spacing: the actual BO-family mean is about 0.077–0.078 Angstrom, while its electronic length is about 0.051 Angstrom. Confusing these two distances produces an apparently similar ratio-7 target with a different physical meaning.

An immediately usable closed form is `log10(rate)=sum_(j=0)^5 c_j T_j(x)`, `x=2 ln(m)/ln(206.768262)-1`, **1<=m<=206.768262**:

| Isotope | c0 | c1 | c2 | c3 | c4 | c5 | Held-out max error /dex |
|---|---:|---:|---:|---:|---:|---:|---:|
| DD | -15.22584583 | 35.33703150 | -10.57785354 | 2.25777754 | -0.37015055 | 0.04805728 | 0.00604 |
| DT | -15.98106619 | 38.57087912 | -11.57947939 | 2.47382240 | -0.40495309 | 0.05274512 | 0.00658 |

This is a fit to calculated rates, not to either anchor. Its largest observed held-out error is 1.53%. For the harmonic family 38/40 fixed-parameter closed forms meet the declared 0.15-dex holdout criterion; their largest error is 0.1013 dex. The DD/DT lambda=0.01, k=1e6 slices are explicitly unflagged (0.1507/0.1671 dex); use the direct solver there. No failed fit is silently presented as a validated response.

### Harmonic material contours and consumer consequence

All separations below are **achieved means**, not harmonic centres. Both tables use F=1. Full centre/mean/stiffness/screening and F controls are retained in `harmonic_thresholds.csv` (480 rows: 204 roots, 276 no-root/domain cases).

| lambda_s /Angstrom | k /(eV/Angstrom^2) | Isotope | 3e-13/s mean | 1e-9/s mean | 4e-9/s mean | 1/s mean |
|---:|---:|---|---:|---:|---:|---:|
| 0.01 | 10 | DD | 0.67741 | 0.52186 | 0.48878 | no root |
| 0.01 | 10 | DT | 0.66584 | 0.52305 | 0.49330 | no root |
| 0.54 | 1e6 | DD | 0.05087 | 0.04509 | 0.04404 | 0.02678 |
| 0.54 | 1e6 | DT | 0.04999 | 0.04464 | 0.04367 | 0.02798 |

At lambda_s=0.03 and k=1e4, zero centre achieves a minimum mean 0.07792 Angstrom (DD) or 0.07764 (DT), but rates only 3.60e-10 and 2.22e-10/s: **there is no nonnegative-centre 1e-9/s root**. F=1e3 controls move that target to means 0.09630/0.09301 Angstrom; F=1e6 to 0.11867/0.11359. These quantify the strength a future nuclear model would need, not what a physical DD resonance presently delivers.

For the explicit reference lambda_s=0.54, external centre=2.78, k=10, the solver achieves mean 2.78673 Angstrom and **log10 rate=-212.031 DD / -231.418 DT** per second. Ratio 7 remains 0.19424: dividing by the achieved pair mean or its instantaneous closest approach is not allowed. A material that can supply only this screening scale and site-centred soft confinement needs a changed pair-confinement mechanism; screening length alone does not create the compact state. This conditional comparison is not a computed PdD/vacancy potential or a result on a driven lattice.

A trajectory histogram of nearest approach, relative speed and orientation does not specify a normalized ground-state contact. Part (c) reads a material interaction without assigning a rate to a collision scan.

### Part (c): mapped target-conditioned pair interaction

From the repository root, run `python -m g6_pairlaw.script.material_pair curve.json`. The input JSON has these fields:

| Field | Meaning |
|---|---|
| `schema` | `coherencev1-target-conditioned-pair-v1` |
| `observable` | `target_conditioned_BO_total_interaction` |
| `R_A`, `V_pair_total_eV` | Paired ascending radii and total interaction energies, including bare repulsion. The producer first subtracts the vacant-cell BO energy from the occupied-cell BO energy and supplies the mapped isolated interaction with its stated energy zero. |
| `reference` | `source`, `configuration`, `energy_zero`, `periodic_to_isolated` (provenance strings); `bare` must be `isolated_unit_charge_coulomb`. |

The reader reports only at supplied radii. It computes `V_bare = C/R`, the signed electronic correction `V_pair_total − V_bare`, and the positive-lowering convention `V_bare − V_pair_total`. Comparisons require the same physical reference and target-relative path. The reader does not subtract an arbitrary far-point constant, turn raw `U_phi` into an interaction energy, or mix the potential with electronic energy loss. Keep the undisplaced and displaced profiles separate; a 2 keV head-on snapshot does not supply arbitrary-angle or energy-dependent screening.

`--isotope DD` or `--isotope DT` additionally solves a **specified prepared radial ground state**. Supply a `prepared_state` object with `preparation`, `inner_continuation_source`, and `outer_boundary` (provenance strings), `V_confinement_eV` at every listed radius (explicit zeros if the well is already in the interaction), and `regular_potential_at_zero_eV` for the *combined* potential after subtracting `C/R`. The reader joins these inputs by interpolating the regular part with a smooth, zero-slope origin, adds `C/R` exactly once, and uses the last radius as the outer wall. It rejects profiles without a resolved single interior well or with a non-negligible outer tail. Those checks do not prove the proposed origin continuation, preparation, or physical far boundary; they must come from the input model. Report masses and states separately under isotope substitution. Without such a state there is no material rate. No mapped material curve has been supplied or run here.
