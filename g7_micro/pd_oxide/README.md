# Palladium with an oxide sheet: does the oxide make a strong site beneath it?

Part (b) of `g7_micro/`, second stage. Kasagi et al. (2002) report the largest screening for Pd covered by PdO. The zirconium stage (`../zr_calibration/`) showed that a single oxygen atom beside the target does nothing. This run tests an oxide sheet: a three-layer Pd slab with a deuteron in its second plane, bare and with a PdO sheet (4 Pd + 4 O) on top, in the same cell. Same pair calculation as before, bare incoming deuteron, same approach direction in both slabs.

## Result

The oxide sheet lowers the static screening at the deuteron beneath it by 1.5 eV at 0.1 A and by 5.0 eV at contact. It does not raise it.

| R (A) | Bare slab lowering (eV) | Oxide-coated lowering (eV) | Coated minus bare (eV) |
|---:|---:|---:|---:|
| contact | 44.8 | 39.8 | -5.0 |
| 0.1 | 39.4 | 37.9 | -1.5 |
| 0.2 | 35.3 | 34.4 | -0.9 |
| 0.4 | 26.9 | 26.3 | -0.6 |
| 0.8 | 16.3 | 16.3 | 0.0 |
| 1.2 | 11.4 | 11.4 | 0.0 |
| 1.6 | 8.8 | 8.8 | 0.0 |

Contact values are the lowering at R = 0 (the delivered signed W0 with its sign reversed).

Where the site sits: the deuteron is buried in the metal, 1.9 A below the top metal layer and 4.4 A below the oxide sheet. Its nearest oxygen is 4.2 A away and its nearest Pd 1.9 A. This run therefore tests a site one metal layer under the oxide. It does not test a deuteron at the oxide/metal interface plane or inside the oxide.

What it means:
- A deuteron under a PdO sheet screens like every other static site computed so far: about 38-39 eV at 0.1 A. The oxide changes it by a few eV, downward.
- Two static routes for oxygen are now closed at this level: an oxygen neighbour atom (Zr) and an oxide sheet one metal layer above the site (Pd). The measured oxide and oxygen effects are hundreds of eV. A site in the interface plane itself, or inside the oxide, is the one static oxide case not yet computed.
- What remains for the oxide measurements is what a static cell cannot hold: a different deuterium population or depth profile under the oxide, charge trapped in or on the oxide under the beam, or the measurement normalisation (`g7_gap/`).

## Carried through the barrier

`action/` (Marvin) turns both curves into what a beam fit would read, with the converter from `../zr_calibration/action/`:

| D-D energy, centre of mass (keV) | Bare slab (eV) | Oxide-coated (eV) | Coated minus bare (eV) |
|---:|---:|---:|---:|
| 0.25 | 43.5 | 39.9 | -3.6 |
| 1.25 | 46.0 | 41.1 | -4.9 |
| 5 | 48.1 | 42.8 | -5.3 |
| 10 | 49.7 | 44.2 | -5.5 |

The oxide-coated site reads 1.6-5.5 eV below the bare site at every energy and for all three ways of joining the curve to its contact value. Kasagi et al. (J. Phys. Soc. Jpn. 71, 2881 (2002)) fit about 310 eV for Pd and 600 eV for PdO on Pd. The computed contrast is opposite in direction and about a hundred times smaller.

## Checks and their limits

- The lead ran the bundled verifier (`python3 capsule/export_pd.py verify capsule`, NumPy only, with GPAW, ASE and SciPy blocked). It replayed all values from the saved energies, and it refused a contact value altered by 1 eV.
- The lead recomputed the lowerings from the delivered pair energies.
- The bare slab is relaxed (forces below 0.04 eV/A on free atoms). The coated slab is only partly relaxed: 0.42 eV/A at the step used. Its relaxation hit the compute cap and a restart failed. No bound on the geometry error is claimed.
- A second relaxation frame was started to show how much the unfinished relaxation moves the result. It did not finish (one of its four states timed out), so no spread between frames exists.
- No k-point or grid check was run in this stage; both slabs are at one k-point and h 0.20 A. The zirconium stage found these checks move the contrast by under 2 eV.

## Three omissions

1. Coated geometry not fully relaxed, and single k-point on a 6.1 A metal cell with no paired check.
2. One thin model: three metal layers, one oxide sheet, one target site in the second plane. No slab-thickness or vacuum check, and no site at the oxide/metal interface itself.
3. Static calculation with a bare probe. No charge on the oxide, no beam-driven state.

## Lead's note

`capsule/` is the delivered archive `g7_micro-pd-base-20261121.tar.gz` (SHA-256 97d55c55be5bb9151b5310a8ff1b0b480933279b425deba8bd4191e731895bb2), unchanged; no GPAW state files or raw logs are in it. The delivered `capsule/README.md` carries the producer's full statement of scope. `action/` is Marvin's delivered converter archive (SHA-256 b9d32d3836202cff59f529371f63c16d908274f86bcfe9941345b44618881793), unchanged. The lead's replay on `capsule/`, with GPAW and ASE blocked, reproduced `action/results.json` byte for byte. `supplement/` is the delivered geometry and second-frame record (SHA-256 47ce009c9dc0d3f22397219c20186dc4e16cc352f6288323f7344d3db257b9ae), unchanged. Compute for this stage: 465 core-min.
