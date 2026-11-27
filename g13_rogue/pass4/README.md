# Pass 4 — conditional whole-foil acoustic transport

Unpublished, additive to base `650a2ebd57c35f371ecd23c4db975ebd6e3cc55a`. Authority: exact entries 4171/4173. This directory owns the producer and its readout; accepted source files are read-only dependencies. No global 475–625-kHz extension, atomistic or independent-patch calculation.

## Kept / dropped before execution

Keep the equal-area square whole foil, 250-µm thickness, the accepted 8×8 grain orientations (seed 173), 1-cm² beam, 0.85/0.15 electrical D+/D2+ parents, 138 µA and 0.690 W. D2+ is one electrical parent, with coherent constituent deposition. Keep both accepted resolved 5–250-kHz loss alternatives and their sunk beam records. Above that band, seven disjoint decades end at 2.5 THz. A first common-band impulse/receiver calculation must decide whether a bulk Green kernel can replace guided resolved propagation. It is not assumed to do so. Retain the resolved kernel if it cannot.

Use a deterministic receiver/image Green construction, not independent local patches. Low thickness branches are obtained from a free-face one-dimensional thickness spectral problem, not from isotropic bulk transport. At shorter wavelengths use crystal-direction Christoffel phase velocities, a stated dispersive acoustic cap, and group velocity for transit/lifetime. Grain transmission/reflection and face/edge budgets conserve packet energy. Missing measured high-frequency Pd losses are conditional lifetime scenarios, not an extension presented as measurement. The Cu 32 dB/m/MHz law is only an explicitly extrapolated analog scenario. No added grain-loss coefficient duplicates explicit interface scattering.

The three largest omissions are measured loaded-foil high-frequency attenuation/texture, microscopic lateral/time launch shape, and microscopic energy transfer from host acoustics to D. They are retained as conditional model uncertainties, not blanket exclusions of any band. Linear one-way transport has no backreaction onto the accepted slow field. Coarse frequency quadrature precedes any refinement.

## Equations fixed before execution

The incident normal impulse is `p_s=sqrt(2 m_s E_parent)`, with `m_s=m_D,2m_D` and the same 5-keV energy per electrical parent. Track transit (0.327–0.408 ps), terminal depth (82–129 nm), assumed lateral Gaussian radius (20/60 nm) and assumed acoustic duration (0.5/2 ps) are different inputs.

For mass-normalized acoustic modes, `qdd+omega² q=F`. A parent contributes the complex free-oscillation impulse

`A_j = p_s e_j,z F_p(k,omega) - i beta (k·e_j) Q_s F_Q(k,omega) / [rho c_p (D k²-i omega)]`,

where `beta=alpha(C11+2C12)`, `D=kappa/(rho c_p)`, the track Fourier transform retains deposition position/time, and the transverse/temporal source transform is Gaussian. This is thermoelastic work from deposited heat, not an electronic acoustic fraction. In bulk, the classical oscillator density of final states gives `dE_j=|A_j|² k² dk dOmega/[2 rho (2pi)³]`; the angular/polarization projection and momentum/thermal cross term are retained. For a free slab the equivalent expression is `dE_j=|A_j|² d²k_parallel/[2(2pi)²]` with thickness eigenfunctions normalized by `integral rho |u_j(z)|² dz=1`; source and receiver depth projections use those same eigenfunctions. Thus no arbitrary energy normalization is attached to built launch. Fourier heat diffusion is a pure-Pd continuum transfer at mesoscopic lengths only; its atomic-scale continuation is conditional and must be labelled.

An assumed comparison reserves the unchanged historical 1% low-band launch and assigns one additional 10%-of-remaining-host-recoil broad Gaussian acoustic envelope above 250 kHz, shared across all seven bands, with the rest retained as random motion/heat. There is no 1% per high band. The built and assumed scenarios are alternatives, not additive incident power. Every species/branch/band work is debited once, with heat and unrepresented spectrum as explicit remainders; an over-budget result is rejected, never clipped.

Energy attenuation length is `L_E=10/[ln(10) alpha_dB/m]` and `tau_E=L_E/v_group`. At each interface `R+T=1`; free faces have unit energy reflection, and the declared edge retention/escape completes the budget. Transit uses group speed; carrier phase uses phase wave vector.

At a receiver, paths retain their electrical-parent identity. Distinct parents in a declared finite packet envelope and participation `(sum w)²/sum(w²)` are different outputs and neither determines coherence. Sparse realizations sum signed complex vector amplitudes from their outset, paired with a parent-phase-shuffled control. Shared images and parents at neighbouring receivers determine spatial support; independent time traces cannot define width. Dense results may use the compound-Poisson packet distribution and correlated moments, not a universal Rayleigh law. Finite-record zeroes, exposure/importance weights and censoring remain visible.

Local elastic energy density divided by the frozen `n_D(z)` is a quotient in eV/D, not a measured microscopic D excitation. The accepted 4.1-ps dissipation/storage proxy is reported separately. Thermal reference is 300 K; population denominators are never fitted to the acoustic output.

## Realized scope and replay

The two launch sensitivity cases change footprint **and** duration together; they do not identify separate effects. The implemented signed field is a mode-diagonal three-branch complex envelope, not a reconstructed microscopic strain tensor. Energy is its cycle-averaged elastic envelope. The nearest-grid pre-existing label is earlier-parent-only energy, not a sub-grid before-arrival measurement.

The 12-energy-length importance core is not the full arriving-parent census. The latter uses first-five-image/eight-sigma front support within a declared finite history; no claim of an exact all-history count is made. Uncensored modeled-path energy and omitted-contribution expectation bounds are separate. Neither participation nor the count alone determines coherence.

Cu uses `exp(-alpha_E*path_length)`; Al uses `exp(-transit/tau_E)`. Analog coefficients are assigned as conditional absorptive rates in the heat/residence ledger, not asserted to be deconvolved Pd intrinsic attenuation. No additional empirical grain-scattering rate is added. Explicit grain reflection is conservative redistribution and remains a receiver remainder.

Accepted dependencies are fixed to commit `2afd47ffb85df8e373abb30dde577620e450c251`; the additive code base remains `650a2ebd57c35f371ecd23c4db975ebd6e3cc55a`. At a checkout already containing that accepted rung2 directory, defaults are root-relative. Otherwise extract the exact commit's `g13_rogue/pass3/rung2/` into a separate dependency directory with read-only `git archive`, and pass `--reference-dir DEPENDENCIES/g13_rogue/pass3/rung2`. `dependencies.json` binds numeric inputs and consumed tables by SHA-256; the capsule does not duplicate those accepted bytes.

`READOUT.md` gives the ordered portable CLI. `refine --reuse-spectrum` consumes the included, source-epoch-checked selected spectrum; omit that flag only when rebuilding the same two deciding spectra. `check --reuse-structural` reuses the unchanged source/guided/interface evidence while exercising the corrected length/time join and actual signed records. No old low-band matrix is rerun. `LEAD.md` is the decision-ready positional readout; archive manifests and the separate delivery receipt bind actual delivery bytes.
