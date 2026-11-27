# G16 item 1: hot electronic occupations in fixed-volume PdD

Fleet `g16_hot_electrons` definition 1, handoff 4234. This is the two-atom material input for the G15 track model, not a calculation of nonthermal electron motion, trapping, or foil heating.

## What the calculation gives

At 300 K the collective stiffness is 0.182970 and 0.176704 eV/Å² on the 12³ and 16³ meshes. The centered D remains force-free by inversion symmetry in this homogeneous cell (maximum raw D force across both meshes <5×10⁻²⁹ eV/Å; the largest combined harmonic force/inversion residual is 1.74×10⁻⁵ eV/Å). This is not a bound on the force at a defect or in a spatially nonuniform electronic excitation.

At the G15 thermal-only peak of 2.33819×10⁻⁵ eV per model D site, the matched-energy interpolations give ΔK of +1.510×10⁻⁴ and −2.875×10⁻⁴ eV/Å². At 1.76354×10⁻⁴ eV/site, hypothetically thermalizing all simultaneously sampled electronic energy, the values are +2.200×10⁻⁴ and −1.388×10⁻³ eV/Å². **Both comparisons reverse sign across the meshes.** Keep both as numerical sensitivities rather than choosing a signed material coefficient. The two differences are not error bounds on the continuum mesh limit.

From 3,000 K upward the two meshes agree on stiffening at the specified tolerance. On the 16³ mesh, K rises to 0.2531, 0.5793 and 1.8846 eV/Å² at 3,000, 10,000 and 30,000 K, with respective excitation energies 0.07368, 0.83029 and 7.42270 eV/site. This energy range is well above the modeled G15 peaks. The corresponding fixed-volume total pressures are 7.900, 15.929 and 71.830 GPa (7.144 GPa at 300 K, of which only the difference is the thermal pressure change). Pressure contrasts also change sign across the two meshes at 310 and 350 K.

The 12³ [100] quartic fit changes from 11.24 to 14.33 eV/Å⁴ between 300 and 30,000 K; the [111] fit changes from 3.41 to −0.19 eV/Å⁴. The latter is a finite-interval fit with an energy cross-check, not a negative harmonic stiffness or proof of a new minimum. Maximum holdout energy residual is 0.000377 eV across the fits. No isolated-site Hessian, defect force or partial-loading result is inferred. Thermal occupations are a coarse stand-in for the first nonthermal femtoseconds.

## Saved output and replay

`response.csv` is the G15 input with two rows per temperature; `response.json` retains full diagnostics and the SHA-256 and byte count of all 110 native points. `track-response.json` gives each mesh’s independent bracket at the two matched G15 energies. `inputs_ledger.md` contains proposed downstream ledger rows, not edits to the foil model. Signed displacements, occupation settings and band-tail checks, forces and stress are in `native/`; `sources/local_space_time.csv` is the G15 producer’s original table.

From this directory with Python, NumPy and ASE, replay without performing another electronic calculation:

```sh
mkdir replay
python readout.py --root native --output replay/response.json
python export_table.py replay/response.json replay/response.csv
python verify_table.py native replay/response.json replay/response.csv
python track_join.py replay/response.json sources/local_space_time.csv replay/track-response.json
python export_ledger.py replay/response.json replay/track-response.json replay/inputs_ledger.md
```

## Model and sources

Use the G2/G3 reference: primitive rocksalt PdD (conventional a = 4.084 Å), fixed positions and strain, Pd.10/H PAW, PW500/PBE, neutral nonmagnetic 11-valence cell, Γ-centered k meshes 12³ and 16³ with the symmetry of each actual displaced cell. Scalar Fermi–Dirac occupations at 300, 310, 350, 1,000, 3,000, 10,000 and 30,000 K. Ions are at neither a fitted temperature nor ongoing relaxation. Keep the energies, signed forces on both ions, and undisplaced full stress from native output, resolving the numerical residual separately. Check the largest occupation of the highest represented band; no acceptance if greater than 1e−6 per spin.

The G15 pass-1 `local_space_time.csv` at `2d85230` (299,176 bytes, SHA-256 `73fa1794047daf82b3dc25247e280b281dba8c6a6421c74ce87694926df18bed`) gives a maximum sampled thermal-electron excess of 2.33819e−5 eV per stoichiometric model D site, or 1.76354e−4 if the simultaneously sampled nonthermal excess is added to the thermal one before maximizing. These scalings multiply the volume density by the fixed 17.029316176 Å³ per PdD formula unit; they do not assign a physical electron temperature to a nonthermal packet or change the foil loading. G2 `primitive.py` SHA `4aabcab9887dddeffc551e3b264b5ef0de9575ba20dfebc9aac5aef283167b06`; calculator `117779df0a231bba1269850b3053386eba691b51b1f428b8207708c5842bfb2a`. G15 pass-2 `LEAD.md` pins the stiffness versus directed-force distinction (`95389b41aeb57ffe3c1bf6ed831712e0548f6ea204dc19051bb7349e18f14d8d`).

## Quantities and verification

- Use the force-consistent electronic Helmholtz free energy F and the native −TS to reconstruct U = F−(−TS). Retain the GPAW extrapolated smearing energy as an independent consistency check, not U. Report ΔU at undisplaced geometry relative to 300 K on the same mesh in eV per formula unit (= per D site). Report signed pressure as −trace(stress)/3, conversion to GPa and ΔP at fixed volume.
- The primitive D displacement is repeated in all cells with Pd fixed. Opposite [100] displacements at ±0.01 Å give the signed collective Γ stiffness, the q=0 sum of the same-sublattice Hessian. The two-atom cell cannot independently separate the single-site partial Hessian or neighbour coupling. At the centrosymmetric undisplaced site the directional thermal bias force must vanish, within the reported raw force residual; nonzero changes of restoring force at finite displacement are distinct. Recover odd forces, even energies and total forces before assuming these symmetries. Compare ±0.02 Å at 300, 350 and 30,000 K. Target force cancellation and inversion within 2e−5 eV/Å and low-T mesh/step resolution 2e−4 eV/Å² or 10% of the resolved ΔK, whichever is larger. Keep the actual differences when the signal is below that scale.
- At k12, obtain [100] and [111] even quartic curves using ±0.12 and ±0.24 Å. Fit K4 from raw antisymmetrized forces with the common cubic K2 and check against the symmetrized free-energy differences. A two-amplitude fit describes that stated interval, not a universally transferable polynomial. Show residuals and the two separate amplitude values. The two-mesh numerical comparison is for ΔK; quartic and pressure have their own actual validation status.
- Compare temperatures at the same cell and the same mesh, then assess mesh sensitivity of the *difference* across the 12³/16³ pair, not the absolute baseline alone. Report secants in Δ(K/K300)/ΔU without silently interpolating the small track-energy range when bracket or precision fails. Compare force curvature with the central second difference of F, accounting for finite-step quartic contributions. An old 0.1 eV smearing result is not a cold reference.
