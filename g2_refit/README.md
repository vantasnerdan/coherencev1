# g2_refit: the ordering model with computed inputs

Inputs are GPAW results from `g2_input` (in `computed_inputs.json` and `adopted-strains.json`). Run:

    PYTHONPATH=../g0_model python3 anharmonic.py computed_inputs.json > anharmonic.json
    PYTHONPATH=../g0_model python3 refit.py computed_inputs.json anharmonic.json adopted-strains.json > refit.json

## The oscillator

Each proton or deuteron sits in its own anharmonic well, `V = K u^2/2` plus a cubic-symmetric quartic, and is coupled to its neighbours by the computed force constants. The well comes from displacing one atom in a 64-atom cell. Its levels are found by direct diagonalisation in three dimensions.

| Material | First excitation, single site | With neighbour coupling | Measured (neutron) |
|---|---:|---:|---:|
| NiH | 96.0 meV | 90.8 meV | 89 meV |
| PdH | 82.3 meV | 80.6 meV | 55.8 meV |
| PdD | 55.7 meV | 54.5 meV | 37 meV |

NiH agrees with measurement to 2 percent. Pd is about 45 percent too high. The Pd H/D ratio, 1.48, is close to the measured 1.51 and away from the harmonic value 1.41. The likely cause of the Pd miss is that the 64-atom well was computed with a single k-point: the converged collective stiffness of PdD is 0.18 eV/A^2, against 1.03 from that cell.

The second excitation is not evenly spaced and splits by direction (NiH: 107 and 117 meV along the first motion, 91 meV across it). That uneven ladder is what saturates a pumped mode.

## Static branch: ordering under expansion

The converged collective stiffness falls steeply as the lattice expands:

| Linear strain | NiH (eV/A^2) | PdD (eV/A^2) |
|---|---:|---:|
| -1% | 1.295 | 0.397 |
| 0 | 1.050 | 0.181 |
| +1% | 0.834 | -0.013 |

Zero-point and thermal motion stiffen the well, so the collective mode softens to zero later than the bare numbers suggest. In a mean-field estimate (well shape and mean-square displacement held fixed, stiffness extrapolated linearly beyond 1 percent):

| Material | Ordering strain at 0 K | at 300 K |
|---|---:|---:|
| NiH | 8.3% | 8.5% |
| PdH | 3.8% | 4.0% |
| PdD | 3.0% | 3.5% |

Deuterium orders at a smaller strain than hydrogen, by 0.76 of a percentage point at 0 K and 0.50 at 300 K. Just beyond threshold the displacement is along the cube diagonal, about 0.12 Angstrom. At zero strain and full loading nothing orders.

## Pumped branch

Mode frequencies from this well: 21.9 THz (NiH), 19.5 THz (PdH), 13.2 THz (PdD). In the rocksalt structure a pump field cannot modulate the well stiffness at first order; the allowed route is second order in the field, so the pump carrier sits near the mode frequency, not at twice it. The pump threshold in stiffness modulation is 0.8 to 1.8 percent with the placeholder loss. Above threshold the mode locks in one of two phases. Loss, the optical coupling and the domain size are not computed.

## What the model needs next

1. The single-site well and neighbour coupling for PdD at converged sampling, then a self-consistent treatment of the whole lattice.
2. The phonon dispersion, for how far a pumped mode travels.
3. The mode's real linewidth. The one value found, calculated for PdH at 295 K (Paulatto et al., arXiv:1411.5628), would put the loss about ten times above the placeholder.
4. The second-order optical coupling and how far light penetrates.
5. Partial loading, vacancies and deuterium arrangement, on a lattice built from these wells.
