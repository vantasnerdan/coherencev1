---
id: g7-hydride-finite-radius-result
description: 'Qualified same-code LiD and atomic-D BO comparison in neutral and bare probe branches'
schema: hydride-screening-readout-v1
status: mandatory-data-qualified
entrypoint: reader.py
---

# LiD anion and same-code atomic-D reference

At R=0.1 Å, LiD gives 40.8523 eV (neutral probe) and43.6614 eV (bare probe), or1.00870/1.07806 times the stated PdD40.5 eV reference and0.80644/1.18051 times the corresponding same-charge atomic-D gas calculation; the PdD reference itself used a bare probe, so its neutral comparison is a charge bracket, not a matched-charge host contrast.

| R (Å) | LiD q0 lowering (eV) | LiD q1 lowering (eV) | Atomic D q0 lowering (eV) | Atomic D q1 lowering (eV) |
| --- | --- | --- | --- | --- |
|0.1|40.85227168|43.66143687|50.65759995|36.98521627|
|0.2|36.36904125|39.13738401|47.36082944|34.91723424|
|0.3|32.39111940|35.11990397|41.66384436|31.04193702|
|0.5|23.28740293|25.65010222|33.65945465|25.75488820|
|1.0|11.75922849|12.56878502|20.49210008|16.89775339|
|1.6|6.36940642|6.07128709|12.46658810|11.55685425|

These are finite-R static target-conditioned BO lowerings, not experimental fitted U0, reaction rates, prepared-state penetrabilities, or a universal material bound. LiD's nominal hydride-anion picture does not assign an exact atomic electron population in this calculation. The neutral/bare bracket changes the electronic ensemble; it must not be collapsed into one projectile charge.

## Definition and matched charge references

All host/probe energies are new PBE/Hae references on the same frozen standard-H PAW-relaxed background. Internal energy is free energy minus the stored negative T*S term at Fermi-Dirac width0.1 eV. No standard-H PAW energy enters the BO quartet.

V_pair = E_probe,occupied - E_probe,vacant - E_host,occupied + E_host,vacant, followed by the accepted same-zero-G periodic point-minus-Gaussian nuclear mapping. Lowering=C/R-V_pair, C=14.399645468667815 eV Å. Hae alpha1=10 Bohr^-2; the pair Gaussian width is sqrt(5)/Bohr. No fitted offset, far-point subtraction or bare-only image correction is applied.

Hosts are neutral. The vacant reference removes only the selected targetD nucleus and ONE electron, not a separately relaxed vacancy or a forced two-electron anion removal. Probeq0 brings one nucleus and one electron; probeq1 brings the nucleus without an electron. F host electron counts are16/15; q0 probe counts17/16 and q1 counts16/15. G hosts1/0; q0 probes2/1 and q1 probes1/0. The gas vacant host is exactly the zero-nucleus/zero-electron Hamiltonian reference, explicitly analytic, never a missing-SCF fallback. The gas one-H zero-electron q1 probe was an actual completed Native calculation.

The periodic charged states use the compensating background of the Native solver. Both probe references have the SAME q and cell, so a common q² background/gauge constant cancels in their difference; the neutral hosts match as well. The source-bound symbolic check includes a wrong-charge sensitivity. This does not cancel geometry-dependent electronic/multipole images: those are limited by the actual F cell check, not declared absent.

## Geometry and method

F starts from the actual g6_gap PBE PW500 eV/k4 hydrostatically relaxed LiD primitive structure: cell off-diagonal components2.01433373614844 Å (rocksalt a4.02866747229688 Å), original stress max0.0001524125 eV/Å³. Default Li.PBE has one valence electron and frozen1s²; standard H has one. This inherited PBE volume is not a new FD volume optimization.

All ions are independently relaxed with standard-H PAW in each finite FD calculation. Base/fine F are2³ primitive repeats(8Li8D), cell control3³(27Li27D); actual rawfmax values are1.3263e-12,9.99615e-6 and2.0642e-12 eV/Å, respectively, below0.1. G is an isolated atomicD reference in a10 Å periodic cubic box. Every approach is head-on along the stored[111] unit vector; F uses the open tetrahedral side and the actual target/path/clearance is preserved.

The common family is non-spin-polarized PBE, FD0.1 eV and Gamma sampling, h0.20 Å with h0.16 Å control. Odd-electron gas states therefore use the same fractional-spin ensemble as the metal calculation, not a separate spin-polarized atomic benchmark. Standard-H geometry uses the PAW predictor; Hae energy SCFs use keep/NullWfsMover and fully reconverge after each move. Actual grid/band/SCF/predictor records are in the saved states.

## Actual F numerical controls at0.1 Å

| Control | q0 absolute lowering difference (eV) | q1 absolute lowering difference (eV) | Tolerance (eV) |
| --- | --- | --- | --- |
|h0.20→0.16|0.82682593|0.91766458|5|
|primitive repeat2→3|3.15721905|0.12999965|5|

Both controls pass for BOTH charge branches using actual relaxed positions. There is no inherited Pd-defect cell waiver or unrelaxed far-field embedding. These representative checks are not an exact isolated limit or a separately exhaustive k-space convergence study.

## Evidence and scope

All four six-R curves and both-charge controls passed independent saved-energy/electron/geometry/source/standalone-state/Green replay. The qualified readout is77715 bytes, SHA2567e6e5c798d922493d837050b97072139ba8241305168f9484b4c0d18ec19c840. Scientific input/source/state bytes remain unchanged; consumer corrections admit the valid empty package marker, actual numeric Nv1 and actual one-point Gamma count without weakening their hash/physics checks.

The positive-R table does not specify a value at R=0. Any separately computed united-atom contact sidecars must be matched to their own charge/energy convention and bound to the exact curve SHA; neither V_total(0) nor a zero contact slope is implied by these data.

Largest remaining method limitations: common nonspin/static-PBE/FD ensemble rather than a thermal ionic/ZPE/dynamical calculation; finite periodic cells and representative controls rather than an exact isolated limit; one ideal ordered LiD environment and atomic-D reference rather than a microstructure/configuration/experimental-gas ensemble. No unconsumed experimental paper supplies a dependent claim.
