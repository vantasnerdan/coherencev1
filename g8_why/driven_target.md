# Driven target — ordinary heating cannot sustain the required reacting tail

**Mechanism (two sentences).** If the ordinary40eV barrier remains, the low-setting floor needs a relative-energy tail of order keV, not a small correction at nominal energy. Beam deposition could prepare targets for later ions, but a cold0.5keV-lab cascade has at most0.5keV total relative energy available: a keV tail must pool several ions' energy or draw on a separately supplied reservoir.

**Distribution and power.** Keep every D/profile/sheet fixed and change only momentum: `fD=(1−f)Mcold+fMhot`, Maxwell scale θ and mean kinetic energy3θ/2. The actual rate weight is `<σ(Erel)g>/vbeam`, with `Erel=L/2+K/2−√(LK)cosα`, integrated through the frozen CSDA response. A harmonic sensitivity uses `θcold=(ℏω/2)coth(ℏω/2kBT)`, provisional PdD ℏω=.047607eV; transplanting it to Ti is **not** a measured Ti spectrum. Saved300–381K temperatures are conditional transport results, not measured beam-on distributions.

For a uniform momentum mixture over the actual traversed D column including its sheet, `f≤ηPdep τ/(ND,column·3θhot/2)`, η≤1. The assumed1cm² footprint and cup-focus convention give D flux **[2.73,4.88,7.90,9.91,10.05]e18m⁻²s⁻¹** and deposition **[146,602,2073,6741,18039]W/m²** at **[.25,.5,1,2.5,6.5]CM keV**. This is an optimistic family-specific bound, not a theorem for arbitrary phase space; recoils never receive extra energy while the original projectile retains its full incident energy.

**Predicted versus required.** All-target thermal+ZPM physical proton factors span [1.040–1.044,1.012–1.014,1.0034–1.0039,1.00059–1.00072,1.000097–1.000135]. At θhot=100eV, even **all deposition to D with τ=100ps** predicts:

| Foil | Physical proton Y/C factor at the five energies |
|---|---|
| Pd | [1.254,1+1.25e−5,1+3.62e−8,1+1.13e−9,1+2.56e−10] |
| PdD | [1.251,1+1.14e−5,1+2.92e−8,1+7.47e−10,1+1.05e−10] |
| Ti | [1.0149,1+3.13e−7,1+4.96e−10,1+9.15e−12,1+1.27e−12] |
| TiD | [1.0149,1+3.13e−7,1+4.96e−10,1+9.17e−12,1+1.27e−12] |

Own-reference shape additionally divides by the changed neutron6.5 factor: its reference is1, whereas physical factors above are not forced to1. Required normalized excesses at.25 are5.23e16–4.90e17, or8.34e15–7.83e16 under alternative current; all five paired comparisons are in `driven_target_floats.csv`.

Even f=1 at100eV fails below1CM: the minimum required fractions, leaving the6.5 reference unexcited, are4.7e6–3.3e7 at.25,2.6e4–1.7e5 at.5 and221–804 at1. **Fractions>1 are outside this model.** A harder1keV Maxwell scale (mean1.5keV, not one cold cascade) needs f(.25)=[1.52,6.62,.231,.547] for Pd/PdD/Ti/TiD, alternative [.243,1.06,.0369,.0873]; all-power minimum lifetimes are [41,181,679,1605]s, alternative [6.55,28.9,108,256]s. Infeasible fractions stay infeasible. At100ps this harder scenario gives only Fp(.25)≈[1.27e5,1.27e5,3.05e4,3.05e4], still orders short.

**What could sustain it?** At10nm, independent-ion overlap is8.6e−8–3.2e−7 for100ps; order-one overlap needs.32–1.17ms (prior quoted-flux result.289ms), still nowhere near the uniform-population energy budget. Even assigning all deposited energy to electrons for1ps bounds the fraction above1eV by5.4e−13–9.1e−12, conditional on one active electron/host; the near-equilibrium part can be Fermi-distributed, but its nonthermal shape and dielectric reaction response are unknown. A long-lived populated state with demonstrable transfer is required, not supplied here.

Times are material/condition-specific: [Paulatto et al.](https://arxiv.org/html/1411.5628v1) predict **stoichiometric PdH at295K**, optical HWHM13.9–46.7cm⁻¹, corresponding to nominal Lorentzian energy-decay scales.057–.192ps but strongly non-Lorentzian spectra—not a measured driven-PdD/Ti lifetime. [Laser-pumped Pd nanocrystals](https://doi.org/10.1038/s43246-025-01063-z) show picosecond electron-to-lattice transfer/tens-of-ps strain redistribution;120ps is an acoustic **period**, not a lifetime or this foil's preparation. Adopted PdD homogeneous loss remains unknown. Storing100eV in the tested K≈1.09eV/Å² well would require13.5Å amplitude versus the tested.0768Å (.0032eV/site), outside the intact-mode model.

**Non-yield observable/scaling.** Beam-on/off target-momentum or vibrational/electronic spectra with resolved decay identify stored energy. One-cascade reactions scale∝I (constant Y/C); an unsaturated state injected∝I and hit by later ions gives rate∝I² (Y/C∝I); thermal/collective onset can add temperature dependence and delayed relaxation. This is distinct from one molecule's coincident constituents. [README](README.md#absolute-floor-semantics) separates all five conditional model reactions/true C from actual SI1 detected counts/cup C; neither establishes a common flat physical floor.

**Exactly three largest omissions:** (1) measured spot, deposition and joint target/electron phase space; (2) material-specific lifetimes, transfer efficiency and nonlinear electronic response; (3) within-cascade correlations and reaction multiplicity beyond the uniform impulse mixture.

**Finite proof:** cold identity exact;48×48/72×72 changed-consumer differences≤1.61e−14, exit0 (`driven_target_proof.json`). **Ordinary heating does not explain it; existing evidence supplies no sustained pooling mechanism.**

