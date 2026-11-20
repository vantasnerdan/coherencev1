# g6_response: actual low-Q screening and electronic-drag outputs

This capsule contains only the actual hosts enumerated in `scope-status.json`. Missing hosts are unmeasured, not zero. It is NOT completion of the original high-Q/four-host deliverable or of a missing adopted host. `table.csv` and host records come from saved actual frequency-dependent microscopic RPA matrices. False numerical eligibility flags and failed convergence results are retained, not hidden.

Keeps: actual PBE/Pd10 or Ti PAW Kohn-Sham bands; full local-field inversion; physicalQ=q+G; complete declared q rules; k/band/G controls; Pd-only eta/frequency controls; exact piecewise-linear loss integration and unmodified loss signs.

Drops: high-Q25/Angstrom short-range screening, now open; NiH/Pd3Ag additions; fitted homogeneous dielectric response; ICRU as a drag substitute. LiD stopping/settings is independently owned by g6_gap. Pd4s/4p shells are frozen: an explicit omitted response, not an asserted zero.

C/lambda is a labelled low-Q Yukawa proxy, NOT screening at0.5/0.1/0.04Angstrom or a contact extrapolation. Host records retain the stopping curve, low-Q head samples, all individual control changes and a coupled empirical error estimate. A qualified static fit does not qualify the drag or merit. ICRU is an untuned separate30percent comparison using actual density/stoichiometry.

Replay an included host with NumPy/SciPy in a fresh output path:

```
python readout.py --directory data/Pd --output replay-Pd.json
```

The three biggest omissions are high-Q short-range electronic response, nonlinear projectile charge transfer/trajectory correlations, and finite-temperature disorder/ionic response. The frozen Pd semicore omission is additional model scope. These are not restored by a successful numerical check. The candidate material lever is strong low-Q static response per unit accessible Doppler-frequency loss; no short-range host ranking follows.

For the first controlled Pd family, screening passes the declared low-Q tests but stopping fails G200-to400 convergence by up to33.97percent. G400 has12–15percent of its computed stopping above physicalQ7/Angstrom at the diagnostic energies. Neither a silent2kF cutoff nor a zero omitted high-Q tail is justified. The saved-output-only repair forecast is a cost estimate, not a computed correction or convergence claim. No other host is fabricated.

Repository copy (lead's note): the 33 saved response matrices `data/Pd/point-*.npz` (about 120 MB) are not committed. They are in the delivery archive `g6-response-lowq-delivery.tar.gz` (SHA-256 1be3de744451f7938d44bd3dced26854849643e80e8c84f1f35e77991a384c10) on the draft release `g6_response-delivery`; their individual hashes are in `manifest.json`. Extract them into `data/Pd/` to replay. The lead's replay from those matrices reproduced every physics value in `Pd-result.json` within 5e-13 relative.
