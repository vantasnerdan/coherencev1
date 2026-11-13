# g0_scales: measured inputs for Ni/H, Ti/H and Pd/D

    python3 scales.py

| Material | a (Angstrom) | H/D site | H or D per metal | Density (m^-3) | Bare ionic plasma energy (eV) | Optic vibration (meV) | Hop barrier (eV) |
|---|---:|---|---:|---:|---:|---:|---:|
| Ni/H | 3.732 | octahedral | 1 | 7.70e28 | 0.240 | 89 | 0.37 |
| Ti/H (delta TiH2) | 4.454 | tetrahedral | 2 | 9.05e28 | 0.261 | 149 | 0.49 |
| Pd/D | 4.084 | octahedral | 1 | 5.87e28 | 0.149 | 37 | none |

Sources are in `scales.py`. Not found by lookup: conduction-electron density, screening length, proton effective charge, and a PdD hop barrier. The first three are what `g2_input/` sets out to compute.
