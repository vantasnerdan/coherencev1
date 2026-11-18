# Pass 5: enumerated 8-Pd pair capsule

This is the complete **single-cell-size** 8-Pd scientific deliverable: twelve final, all-mobile neutral Pd10/H1 PBE cells (four references, four axis and four tetrahedral-start vacancy states), seven saved-density readouts, matched electronic binding, a conditional ionic calculation, finite-temperature source-Eb diagnostics, and three selected conditional yield results. No unfinished 27-Pd endpoint is included. This is not a global configuration/spin search, a qualified vacancy capacity, or first-principles keV screening.

## Results and limits

`scientific-result.json` is an explicitly portable locator extract of the frozen report. Original identity: 53,593 bytes, SHA256 `39d3c516b29f5e036e90a2aabf256114170a44c237815f6201e04944220de98f`. Only locators and consequent extract identities change; scientific numbers do not. `source-provenance.json` preserves original byte counts and SHA256 alongside staged identities for every admitted source. These extracts are **not original executions**. Original producer identities remain bound to exact historical execution-source bytes; public adapters are separate, import/locator-only migrations. No execution is silently re-pinned to the migrated adapters.

Lowest energies among the enumerated local basins select axis D1, tetrahedral-start D2, axis D4 and axis D6 (plus empty vacancy). Binding per added D relative to matched bulk-interstitial insertion:

| Interval | eV/D | Meaning |
|---|---:|---|
| 0→1 | 0.1793562752951381 | single increment |
| 1→2 | 0.2302005440942949 | single increment |
| 2→4 | -0.014379834208964937 | block mean |
| 4→6 | -0.2753018133291967 | block mean |

The axis-only 2→4 block is +0.02732602580698007 eV/D; its 4→6 block is -0.2753018133291967 eV/D. The axis electronic plateau around four D is not a proved global or finite-T capacity: D3/D5, configurational degeneracies, hydride-background changes and differential ionic free energies are absent. Source **Eb=0.23 eV and mtrap=6 remain unchanged**; mtrap6 remains unqualified and must not silently become four.

The finite-T ledger uses the actual 7200 s PdD run-end profile, T=301.40496912917996 K and source mu0=0.22604589791756174 eV. Each electronic block separately fills the original source Eb slot, assuming omitted differential ionic free energy is zero. The alternative equilibrium curves are not a simultaneous occupancy partition, a newly evolved vacancy/trap profile, or a replacement trajectory. The frozen maximum front-odds identity residual is 3.1086244689504383e-15.

All seven midpoint density grids meet the declared five-percent refinement criterion. Source PAW valence density is **not conduction susceptibility**. One D has no distinct intracluster D–D pair; periodic self-images do not create one. The selected D2/D4/D6 spectator shifts relative to the isovolume defect-free Pd8D8 reference are respectively +0.053828604140771716, -0.1739506965798038 and -0.1941659756966383 eV. Their alternative local-valence TF-proxy shifts are -10.365087132163225, -3.4467416540340636 and -3.517991740105746 eV. **Spectator and TF closures are never summed.** The spectator term is a frozen leading coalescence limit, not a solved three-body trajectory.

The yield adapter assigns every transported D to one pair: **100-percent-pairing sensitivity, not a measured pair fraction**. Material profiles, stopping, D inventory, adsorbate and fitted cause parameters remain unchanged. Maximum merit across the declared cases is 1.000286983367914. Literal charge readings and paper-weighted Eq7 fits carry the preexisting phenomenological cause baseline; their keV effective Ue values are not first-principles screening predictions. `results/unchanged-merit-baseline.json` is unchanged exact baseline evidence from integrated commit `7ec4b53c1a9608d7d5f9ae1fd266a89625d93f8c`.

The ionic calculation is an infinite homogeneous local linear drift-diffusion unit-charge model, U/U_static=asinh(s)/s, s=v/(2 D_collective kappa). At 0.25 keV CM its saved-profile cell maxima are 4.115861020629336e-5 eV in the top 100 nm and 0.10721302580300492 eV across the foil. The BZY10 ideal static estimate is 24.190357368148142 eV at the explicit 300 K assumption, 0.03 H/formula unit and cited cubic lattice 4.206 angstrom. No room-temperature mobility or finite-speed oxide value is invented. The below-100-eV trigger is **model-only**, not physical exclusion of ionic mechanisms. Atomic-q, strong-field, correlation, hopping and inertial applicability is unestablished. The complete admitted ionic review (including empty findings, all assessments, limits and dispositions) is in `reviews/ionic-complete-review.json`; its original identity is recorded in provenance. No PDFs or private execution logs are shipped.

## Execution-source dependencies

`_upstream/` contains four **immutable exact-byte snapshots**, owned by canonical Substrate issues236/237/238, not a second maintained implementation. They are execution-source dependencies, not a fork or compatibility shim; math is unchanged. Source commits:236=`81f9bc778ccf160d17957d52dbd8c2f3df209a2b`,237=`15c6881f`,238=`c673d15b0d050e206d78e432858e31dd3510670f`; proposals240/239/241 respectively remain institutional records, not execution prerequisites. Exact file hashes are recorded in provenance. Public imports use only `g4_beam.pairs` and its private `_upstream` package. Historical producer snapshots under `execution-source/` retain their original imports as evidence, not public entry points.

`setups/Pd.10.PBE.gz` and `setups/H.PBE.gz` are the exact original setups for offline GPW density readback. GPW and geometry trajectory bytes are unchanged. Primitive, input, execution, density and derived-result JSON files with rewritten locators are explicitly identified as extracts. Original cgroup administrative strings are omitted. The executed transport source, BCA dependency and tables are frozen exact bytes under `execution-source/transport/`; the actual run-end profile is exact bytes in `inputs/`.

## Reproduction (saved data and pure evaluation only)

Execution environment: Python 3.12.3, NumPy 2.5.3, SciPy 1.18.1, ASE 3.29.0, GPAW 26.7.0 and SymPy 1.14.0. The actual acceptance run records these again in `checks/dependency-versions.json`. Set `PYTHON` to that environment's interpreter. No SCF is needed or authorized by these commands. From an integrated repo containing `g4_beam/merit` and this directory:

```sh
cd g4_beam/pairs
export PYTHONPATH=../..
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export GPAW_SETUP_PATH=setups
mkdir -p checks/reproduced
"$PYTHON" -m g4_beam.pairs.verify_saved_banks --output checks/reproduced/banks.json
"$PYTHON" -m g4_beam.pairs.pair_binding_results --root data/axis --root data/tetra --reference data/reference --output checks/reproduced/binding.json
"$PYTHON" -m g4_beam.pairs.binding_transport_ledger --binding checks/reproduced/binding.json --profiles inputs/run-end-profiles.json --transport execution-source/transport/transport.py --selection-label 'Final lowest enumerated8Pd basins; source-Eb equivalents only, not physical capacity or a new trajectory' --output checks/reproduced/finite-T.json
"$PYTHON" -m g4_beam.pairs.run_ionic_estimate --profiles inputs/run-end-profiles.json --transport execution-source/transport/transport.py --output checks/reproduced/ionic.json
"$PYTHON" -m g4_beam.pairs.ionic_proof --material checks/reproduced/ionic.json --output checks/reproduced/ionic-proof.json
"$PYTHON" -m g4_beam.pairs.pair_density_results --point data/tetra/mono8T-D2/primitive.json --output checks/reproduced/D2-density.json
"$PYTHON" -m g4_beam.pairs.beam_pair_results --cluster-density checks/reproduced/D2-density.json --reference-density results/density/Pd8D8-reference.json --selection-label 'Lowest enumerated8Pd D2 basin' --output checks/reproduced/D2-yield.json
# Expected nonzero exit: authentic D2 bank under a declared D4 target.
"$PYTHON" -m g4_beam.pairs.pair_density_results --point checks/refusal/mono8-D4/primitive.json --output checks/reproduced/must-not-exist.json
```

For a split staging tree, add the read-only integrated baseline repo root after `../..` in PYTHONPATH. The actual acceptance runner records that root symbolically as `$INTEGRATED_REPO`; no machine-specific path is needed in this capsule. `checks/acceptance.json` records actual commands, stdout/stderr identities and exits. The scientific comparison is numerical, not byte equality of changed locator metadata or process timings. The refusal fixture is an unchanged portable authentic axis-D2 primitive placed under the D4 requested target name, not a forged force bank. `verify_saved_banks` checks all twelve bank bytes against report bindings, GPW forces, energies, coordinates, isotope and geometry metadata. Acceptance also consumes the actually selected tetrahedral-start D2 saved density through the changed public path and runs integrated `compare(Target)`.

## Administrative final diagnostic and sealing

The sole final bulk27 numerical discriminator is **closed**. It converged at unchanged PBE/PW300/k222/width0.1 and tight SCF tolerances after a combined zero-initial-moment,157-band and point-symmetry strategy. Actual rawfmax=6.880640629410049e-5 eV/angstrom; energy per Pd differs from the8-Pd crystal by5.7130326764909256e-5 eV, below the preregistered1e-4 target. `final27-diagnostic.json` is an authenticated readback locator extract; exact experimental source and the input extract are under `execution-source/final27/` and `inputs/final27-discriminator.json`. Its original GPW/primitive remain retained evidence, not an additional published material bank. This numerical crystal point **does not qualify D1/D2 binding size convergence**.

`native-final-cost-ledger.json` records all seven actual counters as strings and a conservative60-second tail per unit:227.8346696667 Native core-min. With the20-core-minute pure reserve, the conservative whole cost is247.8346696667<800. **No further27-Pd Native is authorized**, regardless of unused allocation. Protected P254 resources remain untouched.

The complete independent scientific review is in `reviews/pairs-complete-review.json`: empty findings/issues, ten assessments and eleven explicit limitations, supported within the declared closures. Its raw original is130425bytes, SHA256`518a875cb9c49e7d8c13fc27d422bae3c0e80fbaf8d6cebfe888395a205f4b15`; only source locators are made portable. Saved-data acceptance above remains valid at unchanged scientific/model/consumer bytes. The single private draft asset is for the designated integration owner.

## Files in the repository and in the delivery asset

The 12 GPAW `.gpw` banks under `data/` are kept only in the delivery asset `g5_pairs-delivery-20261118t033027z.tar.gz`, on the draft release `g5_pairs-delivery` of `vantasnerdan/coherencev1` (SHA-256 `b0308210…`). They are gitignored here, and `SHA256SUMS` lists them. Every other file is in the repository.

The lead reran `pair_binding_results`, `run_ionic_estimate` and `beam_pair_results` from the saved data, the last from the delivered `checks/D2-density.json`. The binding and yield outputs reproduce exactly, and the ionic estimate agrees to within 2e-16. `pair_density_results` and `verify_saved_banks` need GPAW to read the banks.
