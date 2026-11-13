---
description: Centred cubic PdD zero-pressure PBE control and matched collective curvature
data: cubic_Pd/result.json
producer: cubic_spacing.py
verification: receipts/cubic-frozen-law-comparison-20261113T2345Z.json
---

# PdD reference spacing control

| Centred cubic reference | a (A) | collective K Gamma (eV/A²) |
| --- | ---: | ---: |
| Schirber–Morosin77K author-extrapolated x1 input | 4.084 | 0.1810596906 |
| Zero mean stress in this PBE calculation | 4.134262354 | -0.0547109881 |

Both curvatures are actual two-atom collective D displacements, ±0.01A, PW500/PBE/Pd10, k12, nonmagnetic, Fermi smearing0.1eV. The centred cubic zero-pressure point has mean stress6.5284e-6eV/A³, below the0.0005 criterion. This is electronic force-consistent free energy without nuclear zero-point motion, not a new experimental lattice measurement. The4.084 input is compressed relative to this PBE reference; none of the existing4.084 inputs has been replaced.

The fitted PartA law, unchanged, predicts -0.0578204517 at the independently computed PBE cubic spacing: absolute tensor-entry error0.0031094636. It supports the isotropic-stretch branch locally; it does not validate arbitrary mixed strains or the extrapolated EMT mismatch sites.

The relaxed-film value Kzz=-0.180665958 at fixed in-plane4.084 belongs to a different anisotropically strained centred geometry. Neither that number nor the negative cubic Gamma curvature is localKii, a globally relaxed D well, or a nuclear quantum-state result. Compare localKii, neighbour sums and Gamma only with matching geometry/operator/electronic profile.

The saved force-pair consumer recomputed both signed curvatures and verified every original result source binding. The control cost357276125000ns =5.9546020833core-min is new PartA work; the reused4.084 stress/cost is not charged again. Public numeric/source files only: no PAW/GPW files, primary-paper copies or private paths. Native rerun: in the PartA source tree with its normal dependencies/setups, use a fresh output directory and run: python cubic_spacing.py.
