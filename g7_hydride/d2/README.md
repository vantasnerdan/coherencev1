# Optional D2 finite-R supplement

This separate combined-provenance capsule adds the two G2 molecular-D2 curves. It does not replace the original four F/G curves integrated at `156dd2b`, nor their contact-sidecar bindings. The original four curve files remain byte-identical at their original locations. The F/G copies in this new capsule have identical numerical arrays but augmented combined-readout provenance; integrate only the new G2 curves under a separate path.

## Actual molecular result

A standard-H PAW ionic relaxation in a 10 Å periodic box takes one move from the existing 0.74 Å initialization to bond length 0.7523205149682946 Å, raw maximum force 0.009280931754715655 eV/Å (<0.1). Both nuclei carry D masses. The ray is [111] approaching target0 from outside the bond; the other D stays frozen in every target-conditioned quartet.

| R (Å) | Neutral q0 lowering (eV) | Bare q1 lowering (eV) |
| --- | ---: | ---: |
| 0.1 | 33.7502923206 | 43.0736316626 |
| 0.2 | 30.6855480050 | 39.0736471479 |
| 0.3 | 26.6515426042 | 34.1608095610 |
| 0.5 | 20.3652314588 | 26.3598370695 |
| 1.0 | 11.5705048725 | 14.5415166073 |
| 1.6 | 7.6934163778 | 8.2684342796 |

These are static target-conditioned BO lowerings, not a fitted experimental molecular-gas U0. Relative to the same-code atomic-D reference, molecular bonding changes the neutral and bare branches differently; it is not a universal charge-independent gas correction.

## Method, provenance and scope

Original common PBE, nonspin FD0.1 electronic ensemble, h0.20 Å, Gamma1, 10 Å periodic box. All BO legs use Hae with alpha1=10 atomic units; only geometry uses standard H PAW. Internal energy is free-minus-negativeTS. Target-only neutral removal leaves one H and one electron. Neutral probes add one electron; bare probes add none. Both molecular branches use actual Native occupied/vacant-host and probe states, with full source/geometry/charge/electron/energy replay. The same zero-G periodic point-minus-Gaussian mapping gives V_pair_total; lowering=C/R−V_pair_total. No fitted far-point or atomic offset is used.

The frozen Native epoch is unchanged from the accepted mandatory route, SHA256 `3a530cd7bd8bb167037e59f0375424e6860a5b17aae539af10b8e358f66b50c6`. Molecular run-result is 3053 bytes, SHA256 `04c18513eb9e2bda3de776c67831b963c5312a00c22c13fb131834edb468b173`. Combined readout is 90953 bytes, SHA256 `eccb07724afea37fcca0e6ff236dec201c80936c0106dff8e175a4983d70b4d4`.

The original two F charge-state h and relaxed-cell controls remain passed; they are not a new D2 h/cell convergence test. G2 has one frozen molecular orientation, no rotational/thermal/ionic ensemble, no exact isolated-electronic-limit claim and no independently converged molecular-gas screening benchmark. No G2 contact anchor or derivative was requested. No universal material or rate conclusion follows.

Actual G2 cost 465399343000 CPU ns =7.75666 core-min, below its admitted15 core-min cap. All Native work including pilot/contact failures totals2169877271000 CPU ns =36.16462 core-min; all owned units and watches are settled. Original whole400/100reserve/18:00 allocation is unchanged. P254 resources were not borrowed.

## Reader and proof

```
python reader.py CAPSULE_DIRECTORY
python portability_proof.py ARCHIVE NEW_PROOF_DIRECTORY
```

The reader recomputes all six curves and both F controls from raw/source JSON. Floating comparisons use absolute1e-9 and zero relative tolerance; structures/discrete fields are exact. The separate archive proof must bind the actual archive SHA256 and exercise isolated NumPy/stdlib-only replay, byte tamper and hash-rebound extra-q1-electron refusal. No GPAW/ASE/SciPy or Native execution is needed for replay. Both new G2 curves are also checked with the pinned2645b19a material-pair CLI, without inventing isotope or prepared-state flags.

## Privacy

Scientific JSON/source is unchanged. Only approved GPAW User/CWD/library-location metadata and dataset/GPW machine path prefixes are transformed in text logs; basenames and ORIGINAL/EXPORTED hashes/counts remain. No GPW, PDF, trajectory, private launch/watch output or provider/session payload is included. Draft assets remain unpublished; integration is by Axiom, not a direct Root main push.
