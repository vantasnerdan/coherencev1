# g2_input: GPAW inputs for the ordering model

Computed for NiH and PdD at full loading. GPAW 26.7, PBE, plane waves at 500 eV, Fermi smearing 0.1 eV, metal atoms held fixed, hydrogen displaced by plus and minus 0.01 Angstrom for stiffness. Not re-run on a second machine; `mesh_readout.py`, `collective_readout.py` and `report.py` rebuild the tables from the saved results.

## Collective stiffness (all H or D moving together) against k-mesh

| Material | k-mesh | Stiffness (eV/A^2) | Harmonic energy (meV) |
|---|---:|---:|---:|
| NiH | 6 | 1.174 | 69.8 |
| NiH | 8 (accepted) | 1.050 | 66.0 |
| NiH | 12 | 1.028 | 65.3 |
| NiH | 16 | 1.017 | 65.0 |
| PdD | 6 | 0.421 | 29.5 |
| PdD | 8 | 0.221 | 21.4 |
| PdD | 12 (accepted) | 0.181 | 19.4 |
| PdD | 16 | 0.175 | 19.1 |

Raising the cutoff to 600 eV changes the stiffness by -5.6 percent (NiH) and +1.8 percent (PdD). The harmonic energies sit 26 percent (NiH) and 48 percent (PdD) below the measured 89 and 37 meV: the wells are strongly anharmonic, and the harmonic value is not the vibration energy. PdD is very sensitive to sampling; the k = 6 value used in the first loop was more than twice the converged one.

## Response to strain

Change in stiffness per unit strain, as a fraction of the unstrained stiffness (central difference over plus and minus 1 percent):

| Material | Hydrostatic | Uniaxial, along the motion | Uniaxial, across the motion |
|---|---:|---:|---:|
| NiH | -22.0 | -23.9 | +1.2 |
| PdD | -113 | -80.9 | -15.9 |

PdD at +1 percent hydrostatic expansion has a negative collective stiffness: -0.013 eV/A^2 at k = 12 and -0.018 at k = 16.

## Anharmonic (quartic) coefficients, eV/A^4

`V(u) - V(0) = K2 u^2/2 + K4 u^4/4`.

| Case | Along [100] | Along [111] |
|---|---:|---:|
| NiH, one H displaced in a 64-atom cell | 18.3 | 1.13 |
| PdD, one D displaced in a 64-atom cell | 9.83 | 1.37 |
| NiH, collective, unstrained | 16.2 | 1.35 |
| PdD, collective, unstrained | 11.2 | 3.29 |
| PdD, collective, +1 percent | 9.88 | 3.26 |

The 64-atom values use a single k-point and 400 eV; the same cell gives a single-site stiffness of 1.47 (NiH) and 1.09 (PdD) eV/A^2 with summed neighbour coupling 0.23 and 0.06. A properly sampled PdD supercell run is queued.

## Cost

164 core-minutes of GPAW in total (NiH 80, PdD 85), including failed and stopped attempts.

## Files

`convergence.json`, `adopted-strains.json`, `collective-wells.json`, `wells.json`, `model-inputs.json`: results. `primitive.py`, `converge.py`, `strain.py`, `collective.py`, `well.py`: the GPAW runs. `primitive_*/`, `collective_*/`: per-run outputs.
