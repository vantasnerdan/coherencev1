# Pass-four electronic screening

Scope: lead decisions 1162 and 1228 on Fleet node g4_screen, not the superseded seventeen-state/full-Fourier preregistration. The JSON computation-status fields are authoritative.

## Result and status

Delivered checkpoint: **14/15 DOS cases, all9/9 small-q RPA cases, and both neutral adiabatic cycles**. PdD0.5 DOS is **not computed**: original and bounded recovery terminated before a completed DOS checkpoint. No guessed/interpolated number is supplied. This accepted partial scope follows lead decision 1228; `not_computed_materials` is explicit in the JSON.

The primary static proxy slope, using occupied **Pd-projected** d-centres across Pd−2/0/+2% and the two alloys, is **−8.3547eV/eV**, conditional interval **[−22.3040,+5.4259]**. The neutral two-host adiabatic slope is **−0.06631eV/eV**, observed-control span **[−1.35314,+1.22052]**. Neither sign is resolved under its stated allowance. These are associations, not universal/causal derivatives.

Changing the coordinate to the occupied all-metal d-centre changes the pooled central static slope to+0.6194, interval[−16.3150,+17.5664]. That is definition sensitivity, not positive evidence. Strain-only Pd gives−2.3764, interval[−17.6111,+12.8103]. Geometry, chemistry and projection change together; do not interchange these paths.

## DOS

N(EF) below is **total states/eV per metal atom**, not per cell or per D. Centres/widths here are all-metal PAW d projections, occupied[-20,0]eV, Gaussian0.1eV; width is RMS. The JSON also supplies per-cell normalization, element-resolved moments and D-projector N(EF) per D.

| Case | N(EF)/metal | d centre/eV | d RMS width/eV |
|---|---:|---:|---:|
| Pd−2% |1.937082|−2.684492|1.633261|
| Pd−1% |1.996982|−2.550022|1.558481|
| Pd0 |2.055141|−2.423867|1.487566|
| Pd+1% |2.112425|−2.305210|1.420256|
| Pd+2% |2.168696|−2.193430|1.356328|
| Pd3Rh L12 |1.848451|−2.466938|1.450973|
| Pd3Ag L12 |1.357087|−2.514947|1.473386|
| PdD0.25 |1.334321|−2.244792|1.472389|
| PdD0.5 |not computed|not computed|not computed|
| PdD0.75 |0.515728|−2.458479|1.598287|
| PdD1 |0.496742|−2.684474|1.630624|
| Ni |1.275890|−2.225807|1.347567|
| NiH |0.727075|−2.213837|1.674474|
| hcp Ti |0.856514|−1.263564|0.768266|
| fluorite TiD2 |2.452073|−2.588121|2.293976|

The primary slope instead uses the Pd component consistently: Pd3Rh−2.547382eV and Pd3Ag−2.074253eV. Do not mistake the centroid of an added Ag/Rh d band for a shift of the Pd band.

## Direct low-q RPA and ratio7

Fit q²[epsilon_M(q)−1]=kappa²+regular*q², lambda=1/sqrt(kappa²), using three commensurate finite-q points. epsilon_M=1/epsilon_inverse(0,0), including local fields. **C/lambda is an explicit Yukawa-model proxy, not a full-RPA near-origin potential.**

| Case | lambda/A | Ue Yukawa/eV | lambda/DD | lambda/metal | radial force factor |
|---|---:|---:|---:|---:|---:|
| Pd−2% |0.390771|36.849303|0.144965|0.144965|0.007974|
| Pd0 |0.397479|36.227469|0.144504|0.144504|0.007822|
| Pd+2% |0.403550|35.682445|0.143834|0.143834|0.007605|
| Pd3Rh |0.404211|35.624044|0.145611|0.145611|0.008190|
| Pd3Ag |0.478445|30.096754|0.169307|0.169307|0.018798|
| PdD1 |0.541196|26.607102|0.185133|0.185133|0.028869|
| Ni |0.375779|38.319444|0.154524|0.154524|0.011558|
| NiH |0.451455|31.896082|0.175083|0.175083|0.022198|
| TiD2 |0.435312|33.078879|0.196843|0.139189|0.037811|

DD means the actual ordered D sublattice when present; otherwise the registered ideal octahedral fcc probe sublattice, not real D in a pure metal. Force factor=(1+d/lambda)exp(−d/lambda), **not** a longitudinal pair Hessian or a second screening of already-DFT forces.

## Neutral nonlinear electronic cycle

Fixed spectator and probe ions; neutral host+DDfar versus host+He at siteA, same total valence electrons. SiteA is the body octahedron (six Pd). SiteB is an edge octahedron (six Pd in pure Pd, four Pd/two Rh in L12). This site choice is not an alloy-site average.

| Host | neutral Ue_ad/eV | excess over neutral vacuum/eV | observed control change/eV |
|---|---:|---:|---:|
| Pd |42.253359|−9.262907|0.112697|
| Pd3Rh |42.261550|−9.254716|0.046245|

GPAW reports energies relative to reference atoms. Raw host cross-species differences are not physical totals. The actual H/He GPWs and original PAW XML agree on2Href−Heref=53.74937188eV; restoring that offset gives the totals above. Neutral vacuum2H(doublet)−He(singlet)=51.51626597eV. The negative metal excess is **not** negative screening relative to bare D+. Do not add either cycle value to the independent RPA/Yukawa proxy.

Four/eight-metal periodic-cell+sampling changes and8/10A atomic-box changes give0.267%/0.109% observed changes relative to total Ue_ad. This is a measured control, not a rigorous infinite-dilution/all-electron error bound. Pd/Rh core response remains frozen. The neutral fixed-N result is not an ionic-beam limit.

## Settings, qualification and failures

- PBE/new-PW GPAW26.7.0, ASE3.29.0, NumPy2.5.3, SciPy1.18.1, Python3.12.3; serial rank/domain/band, one thread. Bulk PW500eV, gamma-centred12³, Fermi width0.1eV. RPA96 saved states, response cutoff50eV, eta0.005eV, static direct point integration. Actual geometry/q/settings accompany each number.
- Pd hydrostatic reference is empirical a=3.890A, not relaxed PBE Pd. Alloys/hydrides are separately relaxed ordered phases; Ni/NiH use spin-polarized calculations. Stationarity is not global phase/phonon stability.
- Pd-only tetrahedron N(EF) k12/24/32=2.055141/2.045181/2.103378; requested5% check passes. Matched-cutoff k12/24 lambda0.397479/0.394438A differs0.771%. Band/cutoff/eta and per-case q-window controls are in `controls/`; worst delivered q-window change is7.51%.
- The±10% per-host lambda allowance is transferred from Pd plus individual q-window checks. Its propagated slope intervals are deterministic conditional ranges, **not** statistical confidence or independent all-material convergence. PAW projected-centre uncertainty and responsive-core/model error are not included.
- B terminated at120min: Ag andPdD0.25 checkpoints survived; incomplete relaxation/SCF did not. DOS-only cases had unnecessarily inherited96-state full-Hamiltonian diagonalization. Source now saves all converged SCF states for DOS, requires >1eV top-band margin, and leaves RPA96 unchanged. Actual PdD0.25 tetrahedron N(EF) is identical for96 versus29 retained states; retained top-band minimum is+9.5111eV.
- Original reused Pd0 DOS lacked element moments. Full reanalysis of the exact saved GPW supplies them, bound as an explicit supplement; no Gaussian/element alias. Wrong Pd/Rh projection, H2 vacuum reference and omitted energy-zero offset are refused. Reference-gauge invariance and new producer sidecar were exercised on actual saved outputs.
- Biaxial strain, large-q/Fourier tails, Ti-metal/partial-loading RPA and extra nonlinear hosts were explicitly dropped by1162, not represented as zero/completed results.

## Files and reuse

`scalars.json` contains observables, actual settings/geometry, slopes, neutral endpoints/calibration, qualification and missing-case status. `qualified-cases.json` preserves the complete admitted case records; `input-ledger.json` identifies actual scientific checkpoints even when their owning batch timed out. Native CPU counters remain exact strings; the budget ledger includes failed/terminated work and the conservative20core-min reserve.

Modules separate geometry/relaxation, ground/DOS, response, length/cutoff, saved-output admission and scalar synthesis. `collect_screen.py --case` admits completed individual checkpoints; default complete-tier collection refuses missing states. `--allow-incomplete` is explicit and lists missing materials. `screen_report.report` can resynthesise the scalar associations from the published qualified records without GPAW.

Executed source/input manifests preserve original versions; current producer improvements do not re-pin historical Native records. Admitted DOS spectra and direct-RPA matrices are included for numeric readback. Large GPWs, trajectories and private operational logs are not public attachments; their original byte identities/retained-bank scope are explicit. Recomputing Native requires the recorded engine/PAW versions and new invocation/source bindings. This archive is not a claim that unarchived wavefunctions were independently reproduced elsewhere.

Source: [Yan et al., PAW dielectric response, arXiv1104.1273v1](https://arxiv.org/abs/1104.1273). Raw PAW and installed-engine source hashes, analytic TF/Yukawa/cutoff checks and actual saved-output receipts bind the implementation.

Pass 5: the static contact potential at finite G with the Pd 4s/4p shells open, and the electron-speed closure, are in `dynamic/README.md`.
