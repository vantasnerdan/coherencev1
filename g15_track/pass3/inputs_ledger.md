# Pass 3 inputs and model

## Kept and dropped before code

Keep the frozen Pass 1 positive carrier generator, fixed depth profile and adsorbed sheet, two half-energy molecular siblings with their real axial delay, all six energies, three loss/partition closures and three Cox beams. The new work is an even quadratic Hamiltonian perturbation. Keep quantum position and momentum widths, their symmetrized correlation, signed reversible work, native finite pair cuts, local retarded clocks and dispersion-resolved whole-footprint second moments. Subtract an identically prepared quiet calculation, not a zero-width state.

Use a Gaussian Wigner/local-density closure: an entire reciprocal-cell quadrature propagates each local covariance, and its Fourier characteristic gives the neighbor-relative law. Full-footprint projections of its local anomalous-variance field measure ordered width, not global normal-mode covariance. A separate bounded four-global-mode Hamiltonian consumer supplies actual retained mode-pair covariance on the same frozen driver. The source varies on nanometres, not a tiled eight-site motif. A directional relative normal form uses each native pair cut to test tails. Drop profile evolution and million-site nonlinear trajectories.

## Frozen authority

Accepted checkout: `2d85230f2580efbee4d7de5e8aece2c9540b3c1c`. Governing plan: `coherencev1-g15-pass3-direct-stiffness-squeeze-20261127`; pinned template `research-evaluation` revision `14c5716f18cd5b1b8ee77624079ed31a0879c64388a625f4db88162dda02faf4`. Only this directory is additive. No G13 rerun, DFT, commit, push, merge or publication.

| Input | Source | Meaning and consumer |
|---|---|---|
| Carrier energy and source clocks | `../pass1/inputs.json`, `../pass1/run.py` | Exact accepted positive transfer generator, stopping/source intervals and closure lifetimes. Electron energy is excess internal nonthermal plus thermal excitation, not electron free energy or temperature. |
| Site volume | frozen `lattice.a_A^3/4` | eV per PdD formula-volume, independent of partial D occupancy; physical D counts come separately from fixed D(z). |
| Six signed coefficients and zero | Axiom Pass 3 scope | Fractional curvature per eV per formula-volume: ±0.001, ±0.01, ±0.1; zero. |
| Frozen lattice epoch | Pass 1 lattice snapshot plus `g3_lattice/coherence/lattice.py` | Use one provisional local-clock epoch throughout, including its zero-point-renormalized operator; never splice the corrected 39.8-meV pole into it. |
| Reciprocal strain slope | frozen `collective_strain_slope_eV_A2` | Derivative of the same quadratic Hamiltonian, not clamped-ion electronic pressure. |
| Bulk quantum state | `g10_quantum_site/output/exact/PdD_300K_joint.npz`, `exact_summary.json`, `PdD_300K_joint_g8_joint.json`, basis8 archive | True 300 K position and momentum moments already contain ZPM. Symmetrized q-p is zero by time reversal. Exact radial tails and basis sensitivity remain separate from Gaussian tails. |
| Bulk coupled state | `g10_quantum_site/output/tables/all_quantities.csv` | Coupled 300 K variance and kinetic uncertainty are a sensitivity, not an independent added half quantum. |
| Bulk pair | `g10_quantum_site/pdd-nearest-pair/pair-cut.csv` | Conditional symmetric closing cut; finite native potential support. |
| Six-D cage pair | `g10_quantum_site/pd-vacancy-r4-six-D-pair/pair.json`, `native/reference.json` | Exact original relaxed cage, not surface-vacancy basin. Native axial harmonic transfer supplies an explicit quantum-width estimate, not a measured cage state. |
| Cross section, profile, sheet, exposure and denominator | `g12_residual/inputs_ledger.md`, Pass 1 frozen inputs/tables | Bosch–Hale and common 40-eV comparator; one adsorbed sheet, no same-data fit. 12 h through 2 keV and 2 h at 2.5 keV. |
| Macroscopic exits and temperatures | accepted Pass 1 tables/finite-foil arrays | Retained unchanged reference, plus explicit new stiffness-work delta and its actual onward carrier destinations. |

The full equations, numerical closure, source checks and measured result are in `README.md` and generated tables. Native short-distance energies are conditional shell energies, not nuclear reaction probabilities.

## Local clock and normalization

The spatial quadrature has eight depth cells and six radial annuli. It integrates every original source-channel deposit over the spatial cell. Occupied-D weight uses the exact accepted piecewise D(z) column between its two cell faces, times annulus area; formula weight is volume/(a³/4). The local pulse starts at the recorded sample-depth arrival plus the lateral front. Its duration is the original containing interval capped by a 2-nm local deposition kernel divided by actual projectile speed, not the full cell's transit. Thus no coarse cell energy is divided by a native row's smaller volume or launched at one longitudinal time.

The rejected construction used the full spatial-cell transit time as this local pulse width; it low-pass filtered the fast 2-Omega signal. The corrected replay replaces that clock only. `controls.json` records both clock/filter amplitudes on identical spatial weights, source area and energy per formula-volume. The finite-pulse analytic check and governing carrier normalization are unaffected and reused.

The reciprocal-cell covariance is a leading-gradient Gaussian Wigner closure: opposite-momentum anomalous correlations retain dispersion phases and neighbor characteristics. Its spatial field projections are not a native G10 intersite state or a global normal-mode covariance. The separately orthonormalized four-mode consumer retains ±q through real cos/sin partners, uses the physical D-mass and geometric-volume Grams and the exact accepted G3 bond-gradient law, and reports its conditional Gaussian prior, matched global quiet clock, finite support and omitted-spectrum limitation. Its diagnostic work is not added or debited a second time. The nonlinear pair sensitivity uses the finite cut, shared Gaussian reference-bath noise, ordinary optical feed and paired paths censored at support.

The native initial state and matched quiet operator share a causal local effective-law clock. Before the local pulse their increment is zero. Full-footprint phase and pressure integrate the retarded local age inside each unchanged cell, splitting at original source/profile boundaries. The complex numerator and constant-phase control use exactly the same positive physical D weights and causal support. Positive linear remapping preserves the local budget identity; the measured per-formula error supplies the reported whole-footprint bound.

The cage transfer uses tagged D0/D3 in the original relaxed six-D cage. K_self=(E_A+E_B)/u² and J_close=W/u² define axial COM/relative normal modes K_self−J_close and K_self+J_close, with separate 300-K quantum priors. The relative coordinate cancels COM rather than evaluating a bulk lattice character at a nonlattice cage separation. The three finite displacements give secant-width uncertainty; residual forces, mirrored opening branch and all-points nominal PAW overlap remain part of this estimate. Transverse isotropic transfer and bulk shell dispersion outside the pair are conditional controls, not a measured six-site cage Hessian. The separate coupled G10 local moment and block uncertainty are in `controls.json`; no neighbor paths are inferred from one site's imaginary-time beads.

At fixed electronic temperature the native curvature is a derivative of Mermin free energy, while deposited excitation is a change of internal electronic energy. A finite-temperature Gamma relative mode is not automatically the local partial Hessian. A native G16 result is applicable only with its actual q/branch, occupation, temperature, energy-per-formula reference and supported interval; the six comparison coefficients remain distinct.

Experimental variance oscillations establish the mechanism in other systems, not its PdD coefficient: [Johnson et al., Bi, PRL102,175503](https://pubmed.ncbi.nlm.nih.gov/19518793/), [Trigo et al., Ge, Nature Physics9,790](https://doi.org/10.1038/nphys2788). Model-assisted hot-electron hardening of Au is a scale comparison only ([Descamps et al.2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10857355/)).
