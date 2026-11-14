# NiH x=1 k3 numerical-profile control

Three actual self-consistent calculations: zero and independently measured source-x displacements -0.03/+0.03 A. Original Ni32H32 rocksalt geometry, a=3.732 A, source H index25; all other nuclei fixed. No new quartic, relaxation, principal-axis tensor or mesh-convergence claim.

Profile: nonmagnetic neutral PBE, PW500 eV, gamma-centered 3x3x3 k mesh, Ni10/H PAWs, Fermi width0.1 eV, Davidson2, no point-group operations, density1e-4/force-change1e-4/energy1e-7 electronic stopping. Six genuine MPI ranks, k-point-only decomposition. Both displaced points restart the same fully self-consistent zero checkpoint. Original PW400/Gamma/density1e-6 and new PW500/k3/density1e-4/force-change1e-4 differ in several controls: this is not an isolated k-point derivative.

## Result

Units eV/A²; changes relative to the original source-bound coarse coupling.

| Quantity | Original | New | Change |
|---|---:|---:|---:|
| Local Kii | 1.467942406443983 | 1.712245063783307 | +16.6425% |
| Sum of all actual non-source H pulling Jxx | 0.2342034912155344 | 0.6759926184905987 | +188.6347% |
| Same-matrix H-only q=0 stiffness | 1.2337389152284486 | 1.0362524452927084 | -16.0072% |

All three cross Axiom's approximately10% update threshold. The new nearest-neighbour mean pulling Jxx is0.04046277575362452 across12 distinct periodic-image-folded sites; its sum0.4855533090434942 is NOT the all-shell sum. Preserve negative as well as positive J. Full actual H columns, all64 ion response vectors, folded shell multiplicities and original coarse column are in `dense_Ni/k3/readout-20261114T0917Z.json`.

Force-derived Kii1.712245063783307 versus energy-derived1.7144299353175383 differs0.1274401%. Raw unpatched all-ion acoustic residual/Kii is2.1435157%; no acoustic sum-rule repair. Direct64-ion force-column checks agree within1e-12 algebraically; this is not a1e-12 physical-accuracy claim. The actual full reader reproduced identical parsed output in a GPW/PAW-free cold folder and refused a minus record substituted for the plus record before writing output.

The H-only q=0 value is Kii-sumJ with metals fixed in this fully occupied cubic cell. It is not the coupled all-ion acoustic mode, an anharmonic quantum gap, or an interface-wide spectrum. Old directional NiH quartics remain coarse PW400/Gamma inputs; they were neither refit against the new Kii nor measured at the new profile. Combining those quartics with the new harmonic/coupling values is an explicitly mixed-profile model, not a newly converged full NiH potential.

## Consume

From this extracted bundle, choose a new output filename:

```sh
python nih_readout.py --output dense_Ni/k3/readout-recomputed.json
```

The reader needs source files and the saved profile/three point JSONs, not GPWs or PAW bytes. It validates source identity, original NiH geometry, actual signed positions, genuine communicator layout and common-zero restart provenance. `coupling.atoms` and `coupling.shells` use the original coupling schema; `original_coarse_coupling` and `relative_changes` retain the comparison. Do not substitute pure-metal Ni spacing or truncate the actual column silently.

`nih_profile.py` owns identity/admission; `nih_dense.py` produces Native points; `nih_readout.py` consumes saved numerical data; `run-nih.sh` is the MPI entry. The numerical source profile remains unchanged after production. PAWs and GPWs are intentionally not in the public bundle; hashes are retained.

## Cost and stop

Both Native units succeeded; all planned Native stopped at09:15:28Z. Exact integer counters are preserved as strings. Ni control cost253.4381018333333333333333334 core-min; whole current combined allocation, including charged prior Pd and all failed/rejected Native epochs, is2765.7114098833 core-min plus the conservative20 core-min preparation/review reserve, inside3000. Remaining214.2885901167 is not permission for another family or quartic. See `receipts/nih-cumulative-final-20261114T0920Z.json` and its previously accepted ledger chain. Decimal costs are derived at28-digit precision; raw nanosecond counters are lossless.
