# g6_collision: one deuteron approaching a hydrogen site in PdD

Lump 1 of pass 6. The beam model shifts the collision energy by one fitted constant U_e (0.25 to 2.04 keV from the fits). This directory computes, from first principles, how much the host electrons lower the barrier between an incoming deuteron and a deuteron on its site in PdD. It also tests whether moving the deuterons in step (the lattice state, lump 6) changes that lowering.

## What was run

- Real-time TDDFT (GPAW) of one deuteron at 2 keV lab moving head-on along [110] toward an occupied octahedral D site in relaxed PdD. A matched run aims at the same site with its D removed. The difference between the two runs is the target-conditioned interaction, and the vacant run gives the drag.
- Born-Oppenheimer reference states (electrons relaxed at each frozen position) at seven separations R = 0.1, 0.2, 0.3, 0.4, 0.5, 1.0 and 1.6 A, for both the occupied and the vacant site.
- Two lattice states. The **baseline** is relaxed PdD. In the **displaced** state, every D (the struck target included) is shifted in step by 0.0768 A along [100], with Pd fixed. This is a frozen snapshot of the collective mode at the amplitude the lattice runs give for cold pumped PdD (`g3_lattice/`, `g6_couple/`).
- Lowering is reported as positive: the bare Coulomb repulsion C/R minus the computed pair interaction V_pair. V_pair = E(probe, occupied) − E(probe, vacant) − E(host, occupied) + E(host, vacant), mapped from the periodic cell to an isolated pair.

## Results

Screening lowering against bare C/R, eV:

| R (A) | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 1.0 | 1.6 |
|---|---|---|---|---|---|---|---|
| Baseline | 40.55 | 35.91 | 32.42 | 27.54 | 23.98 | 13.74 | 9.18 |
| All D displaced 0.0768 A | 39.86 | 35.19 | 31.35 | 27.40 | 23.80 | 13.71 | 9.17 |
| Change | −0.69 | −0.72 | −1.07 | −0.14 | −0.19 | −0.02 | −0.01 |

- **Electron excitation during the pass** at 2 keV is 0.22 eV at R 0.1 A, 0.005 of the lowering. This is the dynamic share at one speed; no speed scan was run. The node's second batch (LiD, Pd without D) was to run only if the dynamic share exceeded 20 percent, and at 2 keV it is 40 times below that, so the batch does not run.
- **Drag on the head-on path** (vacant site) is 2.25 eV/A in the baseline and 2.12 eV/A in the displaced cell.
- Both curves are exported in the `g6_pairlaw` material-pair schema (`curves/baseline_BO_pair.json`, `curves/gamma_BO_pair.json`) and pass that reader at commit 2645b19.

## Checks

- Grid spacing, h 0.20 to 0.16 A: the R 0.1 A lowering moves by −0.45 eV.
- Cell size, 2x2x2 cell with one target vacancy: the R 0.1 A lowering moves by +0.25 eV.
- Both are within the declared tolerance of max(5 eV, 10 percent), so the original cell stands for all four production runs.
- Force-work identity along each trajectory: residual at most 0.012 eV.
- Time step: a half-step rerun of the baseline vacant trajectory is running. It checks the drag only, with acceptance within 10 percent. Encounter excitation is reported from the production step and is not time-step checked.

## Three largest omissions

1. One speed (2 keV lab) and one head-on direction. The displaced cell is a frozen snapshot with an imposed amplitude, not a computed locked state.
2. Separations 0.1 to 1.6 A only. There is no extrapolation to R = 0 and no rate here; the pair rate comes from `g6_pairlaw/` reading the exported curves.
3. Encounter excitation is not time-step checked, and drag is for one finite head-on path.

## Lead's note

Repository copy: the 30 saved density and kernel arrays (`*.npz`, about 10.4 MB) and the 44 GPAW text logs (`*.txt` other than `README.txt`, about 2.4 MB) are not committed. They are in the delivery archive `collision-delivery.tar.gz` (SHA-256 8aaaf01cae63199cacd3caeddd9cb34a0b7860e1cf6165bcd47e75715da50cbe) on the draft release `g6-collision-pass6-20261120`, and their hashes are in `manifest.json`. The reader requires exactly the manifest files, so extract the archive into a fresh directory and run `python reader.py <that directory>` to replay. This README is a lead's addition and is not in the manifest. The lead's replay (NumPy 1.26.4, isolated Python, SciPy, ASE and GPAW imports blocked) reproduced the delivered output to 6e-13 eV and refused a curve altered by 1 eV.

What it means:
- The static first-principles lowering in PdD is about 40 eV at 0.1 A. It is settled to under 0.5 eV in grid and cell size. The electron response during a 2 keV pass adds 0.2 eV to it.
- Moving every D in step by 0.077 A changes the lowering by −0.01 to −1.07 eV across the curve, slightly reducing it. At this amplitude the lattice state does not supply the keV screening the beam fits ask for.
- The fitted 0.25 to 2 keV must therefore come from something a bulk PdD cell does not contain. The candidates still open are the near-surface causes carried in `g6_plateau/`: surface layer, fast neutrals and rising screening. Swapping D for T in the target separates them: at 0.5 keV lab, fast neutrals give 55.6 times the yield and a target cause gives 2.5 to 3.2 times. `g6_couple/` part (d) carries how the change scales with displacement amplitude.
