# LiD / atomic-D signed contact anchors

Separate supplement to the unchanged finite-R delivery, archive SHA256 `7b553071782e30084cfc67a3e473aff06c7bf2972aa6dbbea7e4dc6c9f4a23f6`, integrated as `156dd2b`. The four copied positive-R curves are byte-identical to that delivery. No R=0 row is inserted into their arrays.

## Observable and measured grid sensitivity

`W0 = lim[R→0](V_pair_total(R) - C/R)` is the signed regular contact term. Negative W0 denotes lowering; it is not the singular total potential at zero and is not an experimentally fitted screening energy.

| Host / probe | W0, h0.20 (eV) | W0, h0.16 (eV) | Absolute h difference (eV) |
| --- | ---: | ---: | ---: |
| LiD / neutral q0 | -47.2127459993 | -43.9928323873 | 3.2199136120 |
| LiD / bare q1 | -50.0929582346 | -46.7542643528 | 3.3386938818 |
| Atomic D / neutral q0 | -57.5258376217 | -53.3046219386 | 4.2212156831 |
| Atomic D / bare q1 | -41.9202506296 | -38.9855516960 | 2.9346989337 |

All four meet the registered `max(5 eV, 0.1*abs(fine W0))` criterion. These are two-point grid sensitivities, not extrapolated h→0 values. Physical positions, species, cell, masses, target and ray remain fixed between each pair; independent relaxation/force/GPW provenance may differ. A 0.01 Å moved-geometry negative was refused for both hosts. Standard-H PAW geometry and common nonspin PBE, FD0.1 occupations, Gamma-only Hae BO references are inherited from the finite-R calculation. Fine gas neutral-host energy/log is the actual completed q1 fine host, reused byte-exactly by q0.

## Absolute operator convention

An ordinary He PAW dataset is not substituted for two Hae nuclei. GPAW's H-only all-electron source is preserved under each raw source tree. The public prebuilt BaseSetup path constructs a faithful Z=2 operator with the same alpha1=10 and alpha2=300 (atomic units), zero atomic reference, Delta0/vbar doubled and M quadrupled. The actual built radial arrays, source assembly and meaningful wrong-M/vbar/atomic-offset sensitivities are independently checked by `verify_contact_mapping.py` and `contact_replay.py`. The coincident-two-H Native input guard is not bypassed and no synthetic coincident SCF is claimed.

For each charge branch, `Q_He = E(He+spectators) - E(H+vacancy spectators) - E(occupied host) + E(vacant host)` uses actual source-bound internal energies, `free-minus-negativeTS`, with the same charge and electron counts as its finite-R quartet. q0 adds one probe electron and q1 adds none.

`W0 = Q_He + B_native - B_source + regular_images` follows the source-functional identity and the same zero-G point-minus-Gaussian mapping as finite R. `B_native=68.65793304732587 eV`, `B_source=68.65793360383724 eV`; their -5.5651136677e-7 eV unit-constant residual is retained, not silently set to zero. There is no fitted atomic, far-point or host offset. Analytic regular-limit series and numerical approach-to-limit checks are retained in each sidecar.

## Scope and derivative

The cells are the original finite periodic LiD repeat2 and atomic-D 10 Å box. Matched-charge gauge/background cancellation does not remove all electronic/multipole images. Existing finite-R LiD cell-size checks are not contact-cell convergence checks. Neither an exactly isolated electronic limit nor a universal material bound is claimed.

`derivative_eV_per_A` is null: no matched internal-energy derivative/entropy-force calculation was qualified. No zero-slope or symmetry assumption is made.

## Saved-output reader and proof

With Python and NumPy, from this directory:

```
python contact_export.py --replay .
```

The reader recomputes all four sidecars from archived Native states, operators, original finite-R sources and curves. Structures/discrete values compare exactly; floating values use absolute tolerance 1e-9 and zero relative tolerance. It performs no Native calculation and needs neither GPAW, ASE nor SciPy.

To prove this actual archive, including isolated `python -I -S` replay with forbidden imports actively blocked and source-aware negative variants:

```
python contact_proof.py CONTACT_ARCHIVE.tar.gz NEW_PROOF_DIRECTORY
```

The separately published proof binds the actual archive SHA256. Its negative fixtures are modified copies only, never positive or fabricated Native data: a byte change, an extra raw q1 electron after standalone/outer hash repair, and a 1 eV wrong W0 after sidecar hash repair must refuse.

## Privacy and provenance

`contact-sidecars.json` binds the scientific sources, sidecars, curves and log inventory. Scientific JSON/source is byte-exact. Only approved GPAW User/CWD/library-location metadata and dataset/checkpoint machine path prefixes are removed from text logs; dataset/GPW basenames remain. Every log retains ORIGINAL and EXPORTED byte counts/SHA256 and transformation counts. Missing original private header bytes cannot be reconstructed from redacted logs; the export receipt attests that operation and the archive hash binds the exported inventory. Unknown private path forms are not silently scrubbed. No GPW, PDF, trajectory, account/session/provider payload or launch/watch administrative output is included.
