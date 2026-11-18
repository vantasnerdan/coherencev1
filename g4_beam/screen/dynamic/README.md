# Static screening with the Pd inner shells, and the electron-speed effect (g5_dynamic)

## Result

PdD minus Pd in the static contact potential is **+13.58 eV** in the main calculation (N3, 1200 eV plane waves, 18-electron Pd setup). Opening the Pd 4s and 4p shells adds almost nothing. Within this fixed-lattice PBE/RPA model, the inner shells do not supply the hundreds of eV in the fitted loading contrast.

| Contact potential (eV) | Pd | PdD |
|---|---:|---:|
| U10, PW500 | 17.195 | 30.858 |
| U18, PW500 | 17.307 | 30.863 |
| U18, PW1200 | 17.308 | 30.890 |
| U18 − U10 at matched PW500 | 0.112 | 0.005 |

All rows use k12, G200 and N3, with 160 and 156 empty bands. The last row is the difference between two PAW setups at matched settings. It includes the changes in representation and ground density, so it is not an isolated core contribution.

At N2, raising the response cutoff from G200 to G400 adds 2.96 eV to the loading contrast. That is a sampled sensitivity, not a bound on the omitted G, and the 1 eV / 10 % precision flag stays false. The largest potential over every sampled profile is 33.9 eV (PdD, G400). Every main N3 U18 profile is below 31 eV.

| Fitted loading benchmark | Rows | Range (eV) |
|---|---:|---:|
| Transported cause × charge refits, all readings | 9 | 213.3–839.6 |
| Cause parameters selected for the targets (subset of the above) | 3 | 379.4–609.5 |
| Pass 4 uniform thick-target bulk fits | 6 | 550.8–792.1 |
| Direct Eq 7 curve fits | 9 | 891.0–1977.5 |

These are fitted values with different definitions. They are not potentials to add to U. `results/fit-reference.csv` holds all 24 rows with their selection flags.

## Electron speed

The ten-valence low-q closure (Huke et al., arXiv:0805.4538) gives PdD minus Pd at four laboratory energies:

| Laboratory energy (keV) | Pd (eV) | PdD (eV) | PdD − Pd (eV) |
|---:|---:|---:|---:|
| 0.5 | 31.7 | 26.2 | −5.55 |
| 2 | 14.2 | 24.7 | 10.53 |
| 5 | 8.7 | 20.1 | 11.42 |
| 13 | 5.3 | 10.4 | 5.03 |

Through the merit consumer (`results/beam/`), the Pd profiles and stopping stay unchanged. The PdD fitted family is the common calibration, and only the Pd-minus-PdD speed shift is added at the nominal incident speed. The 42 rows cover all three causes and both charge readings. Merit relative to PdD is 0.646–1.017. That is neither a bare enhancement ratio nor an isolated speed-effect ratio. The speed closure is not rebased onto the new static U and is not applied to the bound 4s/4p electrons.

## What the model keeps and leaves out

It keeps neutral, nonmagnetic, ordered, fixed-lattice PBE with RPA response, and the static contact potential at 11 radii from 0 to 1 Å with full local fields. The speed step uses one isotropic intraband velocity and a frozen shift at the incident speed.

It leaves out the nonlinear D-D collision, charge capture, ionic response, disorder, surfaces, response beyond RPA, and full frequency and local-field dynamics. No profile exceeded 100 eV, so the real-time step was not started.

Compute: 502 of 800 core-minutes, including failed attempts.

## Files

- `results/static-profiles.csv`: 16 models and controls at 11 radii.
- `results/screening-report.json`: the matched activation; the q, G, band, k, plane-wave and reference controls; the benchmark comparison.
- `results/speed10.json`: four laboratory energies, the full paper grid and the Pd velocity-grid control.
- `results/beam/`: the 42-row Target comparison.
- `results/fit-reference.csv`: the 24 benchmark rows.
- `results/atomic-admission.json` and `proof/`: diagnostics of the 18-electron setup, and checks of the speed closure and the quadrature.
- `data/static/registry.json`: the saved response banks, read without GPAW.
- `ledger.json`, `capsule.json` and `static-source-provenance.json`: source identities and the file map.
- `COMMANDS.md`: the readout, GPAW and source-check commands.

## Files in the repository and in the delivery asset

These are kept only in the delivery asset `g5_dynamic-delivery-20261118t104657z.tar.gz`, on the draft release `g5_dynamic-delivery` of `vantasnerdan/coherencev1` (SHA-256 `10262d46…`). They are gitignored here:
- the 13 `.gpw` files;
- `data/static/banks/`, the saved response matrices (34 MB);
- `data/atomic/Pd.18.PBE`, the generated 18-electron setup (2.3 MB);
- the third-party copies used only by the optional source checks: `data/primary/yan-paw-response.txt` (paper text), `data/speed-proof-sources/huke-page10.png` (a paper page) and `data/speed-proof-sources/chi0_drude.py` (GPAW source).

`SHA256SUMS` lists every file, including these.

The lead reran the six readout steps of `COMMANDS.md` on the asset tree with NumPy 2.5.4 and SciPy 1.18.1. Four files match byte for byte: the atomic diagnostics, the static profiles, and the fit-reference and report CSVs. The speed, Target-comparison and JSON report outputs agree to within 2e-15 relative; the JSON reports also record the hashes of their regenerated inputs. From the repository alone, the speed, Target-comparison, fit-reference and report steps run, and the report step reads the committed `results/static-profiles.json`. The atomic and static steps need the asset.

Sources: Yan et al., PAW dielectric response, arXiv:1104.1273; Huke et al., arXiv:0805.4538.
