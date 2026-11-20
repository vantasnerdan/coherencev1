# g7_hydride: screening at a deuteron in LiD, and the same calculation for an isolated D atom

The question: in a saline hydride the hydrogen sits as an anion carrying its own electrons. Does that give an incoming deuteron more static screening than the electron sea of a metal? An isolated D atom is computed with the same code and settings, so lattice and gas are compared within one method.

## Result

Static Born-Oppenheimer lowering at R = 0.1 A, by the `g6_collision/` method, with two charge states of the incoming deuteron:

| Target | Neutral probe (eV) | Bare probe (eV) |
|---|---:|---:|
| D on its anion site in LiD | 40.9 | 43.7 |
| Isolated D atom, same code | 50.7 | 37.0 |
| PdD (`g6_collision/`, bare probe) | | 40.5 |

At 0.1 A the hydride anion gives the same static lowering as the metal electron sea: 1.01 times PdD with a neutral probe and 1.08 times with a bare probe. Against the isolated atom with the same probe charge it gives 0.81 times (neutral) and 1.18 times (bare). Lowering against R from 0.1 to 1.6 A for all four cases is in `capsule/README.md`. The curves are exported in the `g6_pairlaw` material-pair schema under `capsule/curves/`.

What it means:
- Metal, insulating hydride and isolated atom all give 37-51 eV at 0.1 A. Static screening at this separation is set by the electrons right around the colliding pair, and depends little on the host. The anion's own electrons add no more than the metal's electron sea.
- This agrees with the published ab initio cluster value for LiH, 40.5-44.7 eV (Huke et al., Eur. Phys. J. A 35, 243 (2008)). It also agrees with the measured gas-like screening in insulators (Raiola et al. 2002).
- The several hundred eV to keV measured in metals is not a property of an ideal static site of any host type tried so far. It has to come from what ideal static cells leave out: the distribution of real sites near the surface (`g7_micro/`), local chemistry such as oxygen, or the collision dynamics.

## Checks

- Grid spacing, h 0.20 to 0.16 A, at 0.1 A: +0.83 eV (neutral probe) and +0.92 eV (bare probe).
- Cell size, 2x2x2 to 3x3x3 primitive repeats, both relaxed: 3.16 eV (neutral) and 0.13 eV (bare).
- All four changes are within max(5 eV, 10 percent).
- Forces after relaxation are below 1e-5 eV/A in all three LiD cells.
- The structure starts from the relaxed PBE LiD of `g6_gap/`.

## Three omissions

1. A static, non-spin-polarized PBE calculation: no thermal or zero-point motion of the ions and no collision dynamics.
2. Finite periodic cells with representative controls, not an exact isolated limit.
3. One ideal ordered LiD environment and one atomic reference, not a microstructure or an experimental gas.

## Contact values

`contact/` holds the value at contact (R = 0), which the `g7_micro/` converter needs to turn these curves into a barrier-penetration shift. Each value is the regular part of the pair interaction as the two nuclei meet, computed with a helium nucleus on the target site. Lowering at contact, fine grid (h 0.16 A):

| Target | Neutral probe (eV) | Bare probe (eV) |
|---|---:|---:|
| D on its anion site in LiD | 44.0 | 46.8 |
| Isolated D atom, same code | 53.3 | 39.0 |

Each lies about 3 eV above its 0.1 A value, and the picture is unchanged: the hydride and the isolated atom screen alike at contact. Changing the grid from h 0.20 to 0.16 A moves each value by 2.9-4.2 eV, within the registered tolerance. The slope at contact was not computed, so the converter treats the shape between 0 and 0.1 A as a stated sensitivity.

## Lead's note

`capsule/` is the delivered archive `hydride-screening-delivery.tar.gz` (SHA-256 7b553071782e30084cfc67a3e473aff06c7bf2972aa6dbbea7e4dc6c9f4a23f6, 203 files), unchanged. The GPAW logs in it had machine user and path lines removed before export; their ORIGINAL and EXPORTED hashes and the transform counts are in `capsule/manifest.json`. Replay with `python capsule/reader.py capsule` from this directory (NumPy only). The lead's isolated replay under NumPy 1.26.4, with SciPy, ASE and GPAW blocked, qualified all four curves and both controls, and refused a curve altered by 1 eV. Compute used: 23.5 core-min of a 400 core-min budget, including a failed pilot publication.

`contact/` is the separate supplement `contact` archive (SHA-256 6f9429941ab6586621320d9947b5e93454b34c49e6e9b140e34f0160d609c2e3, 593 files), unchanged. Its four copied curves are byte-identical to those in `capsule/`. Replay with `python contact/contact_export.py --replay contact` (NumPy only). The lead's isolated replay passed and refused a contact value altered by 1 eV. Compute: 4.9 core-min.
