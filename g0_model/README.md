# g0_model: which ingredients give an ordered state

Coupling protons to the transverse field alone gives no ordered state (`g0_stability/`). This script keeps that coupling as a stable baseline and adds five ingredients as switches, testing all 32 combinations for NiH, TiH1.84 and PdD:

1. screened dipole-dipole attraction between proton displacements
2. a quartic (anharmonic) term in the well
3. coupling to the conduction-electron response
4. drive and loss
5. a gradient (stiffness) term

Run:

    python3 domain.py > domain.json

## Result

Only dipole attraction together with the quartic term gives a bounded ordered state by itself: a static collective displacement of the protons, appearing when the attraction exceeds the well stiffness. Whether a real hydride reaches that point is decided by how strongly the metal's electrons screen the attraction.

| Material | Screening length | Loading threshold | Ordered at reference loading | Displacement | Energy per H or D |
|---|---:|---:|---|---:|---:|
| NiH | 0.56 Angstrom (Thomas-Fermi) | 9.5 | no | | |
| TiH1.84 | 0.61 Angstrom (Thomas-Fermi) | 18.3 | no | | |
| PdD | 0.59 Angstrom (Thomas-Fermi) | 6.0 | no | | |
| NiH | 5 Angstrom | 0.55 | yes | 0.22 Angstrom | -22 meV |
| TiH1.84 | 5 Angstrom | 2.4 | no | | |
| PdD | 5 Angstrom | 0.29 | yes | 0.43 Angstrom | -104 meV |

With drive and loss switched on the protons follow the drive as a forced response: 21.5 THz (NiH), 36.1 THz (TiH1.84), 8.9 THz (PdD), extending about 2 to 2.6 nm from a localised drive and ringing down in 0.9 to 4 ps once it stops. This is a response, not a self-sustained state.

Quartic alone, electron coupling alone and the gradient term alone do not order. Dipole attraction without the quartic term runs away and is rejected.

## What this leaves open

- The ordered state found here is static. A self-sustained oscillating state needs gain, which this script does not have.
- The gradient term gives a correlation length but does not pick a domain size.
- The quartic coefficient, stiffness, damping and the 5 Angstrom screening length are placeholders. Screening and the well shape are the inputs to compute next.
