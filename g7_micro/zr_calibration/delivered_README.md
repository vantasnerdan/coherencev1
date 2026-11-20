---
id: g7-micro-zr-base-hk
description: Source-bound static BO and regular contact contrast for clean and O-adjacent ordered epsilon-ZrD2
status: paired-h-k-qualified-oxygen-cell-omitted
entrypoint: export.py
scope_authority: entry_2782
---

# Zr: interstitial-O contrast

In this ordered epsilon-ZrD2 proxy at base h0.20/Gamma, adding interstitial O next to the target D reduces the R0.1 lowering by 0.587715 eV and makes the signed contact W0 2.421958 eV less negative.

| Method | Clean R0.1 lowering (eV) | O-adjacent R0.1 lowering (eV) | Clean signed W0 (eV) | O-adjacent signed W0 (eV) |
|---|---:|---:|---:|---:|
| base h0.20 Gamma | 39.103402 | 38.515687 | -41.654941 | -39.232983 |
| h0.16 Gamma | 39.895059 | 39.361347 | -41.905798 | -41.261403 |
| h0.20 k222 | 37.388326 | 37.378212 | -39.920012 | -38.077843 |

Both paired h/k changes pass the commissioned5eV criterion at R0.1/contact. This is not statistical significance, full-curve refinement, an isolated-cell limit or a paired O-cell qualification. The exact paired changes and their actual source/control map are in paired-comparisons.json. The qualified=true field is scoped explicitly to the revised paired h/k commission, with paired_cell_qualified=false.

## Three omissions

1. **Paired O large-cell check omitted by Axiom2782.** Its all-position geometry stopped at its registered CPU cap after step2 reached fmax0.2856681470219027eV/A (about0.29), above the required0.05. No accepted O large geometry or O-cell BO exists; no continuation or force waiver was bought. No O finite-cell uncertainty bound is inferred.
2. **Qualified contact derivative omitted.** m0 is null, not zero: the exported internal energy is free minus negativeTS and no matched entropy/force derivative was commissioned. The nominal contact continuation uses the chord; zero-contact-slope Hermite and linear-W forms are shape sensitivities, not an uncertainty enclosure.
3. **Full-R refined/control curves omitted.** The completed h/k controls have R0.1/contact only, not full refined positive-R curves; action-level refinement sensitivity is not bounded by these endpoints. No completed cell-control case is included. The two base curves retain all six actual positive radii.

## Clean-cell availability

This snapshot contains the six complete base/h/k cases. The clean large-cell run completed four finite-R quartet legs, then failed with KohnShamConvergenceError in the united-He SCF at21:28:59Z; no converged clean-cell W0/full case exists. The current exporter requires W0 for an included cell case, so its finite-R-only result is omitted without parser rework, as Axiom instructed. Raw converged legs/checkpoints and the failed contact log are preserved. No new Native or completed supplement is promised; any later admitted clean-cell shift is a single-clean-site result, never the missing paired O-cell check.

## Physical and numerical convention

The target-conditioned q+1 incident nucleus brings no added electron. Spectators, cell and incoming ray are frozen per actually relaxed frame; hosts are neutral and differ by the target H electronic nucleus/one electron. Isotope mass is D. Zr12/O6 core-frozen PAW geometry is all-position force-qualified; BO Hae and public Hae-derived Z2 contact use the established internal-energy, absolute-energy and same-zero-G point-minus-Gaussian mapping. Raw source/Ne/state/energy/geometry/operator proofs are replayed, not replaced by analytic or synthetic states. Contact is finite signed W0=lim(V-C/R), not singular V(0), positive screening magnitude or Eq7 Ueq. Negative W0 means lowering.

The material is an ordered stoichiometric epsilon/FCT proxy with one O per8Zr, not a reconstruction of the experimental alpha/delta mixture or a unique microscopic site assignment. Reported retention is a geometric criterion, not global stability. No full thick-target yield, fitted experimental shift, universal multiplier or dynamical-screening conclusion is claimed.

## Replaying and provenance

With Python and NumPy installed, from this directory run:

```sh
python export.py replay .
```

The reader needs NumPy and stdlib, not GPAW/ASE/SciPy. Frozen producer files are hashed/parsed data, never executed. Scientific comparisons use1e-9 absolute and zero relative tolerance. The archive retains the actual bound GPW/source/setup bytes. GPAW User log redaction carries ORIGINAL/EXPORTED bindings; replay does not generate bytecode or mutate the strict inventory. Curves are unchanged target-conditioned-pair-v1 arrays; the two separate contact sidecars are SHA-bound to their corresponding curves. Isolated NumPy replay and actual byte/Ne/W0/migration-diagnostic refusal receipts accompany the private draft.

## Wider part(b)

The Zr delivery is the revised2782 scope only. Pd/PdO and all later site classes are not reported complete. After this delivery the next authorized work is a smaller Pd/oxide geometry/inventory and short measured-SCF-cost pilot, with a proposal before production. The original whole deadline remains23:00Z; no source/main push or public release is authorized.
