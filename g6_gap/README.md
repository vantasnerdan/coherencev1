# LiD low-energy stopping and conditional LiT return

HAND-ON: within the declared ordered-PBE/retained-ZBL model, yes: there is a nuclear-dominated low-energy window, up to 1.881–1.901keV lab on q4. At 1keV nuclear drag is 84.21%; break-even needs 52.33b (51.54–53.13b observed-control band), versus the supplied 200b TiT1.7 reference. This is not measured LiT or a physical gain verdict.

## Replay

Requires Python and NumPy only; no GPAW, ASE or SciPy. From this directory run:

    python replay_gap_loss.py --output replay-results

Use a fresh output directory. It verifies input hashes, contracts every saved loss array, reconstructs the matched seven-control envelope, and produces:
- `tables/stopping-band.csv`: D-in-LiD electronic stopping, 0.1–13keV; q4, low/high observed envelopes, source ICRU and inherited sqrt references separately. The labelled nuclear column is D-on-Li6T ZBL, not an isotope-swapped LiD measurement.
- `tables/G-band.csv`: D-on-Li6T G over 0.5–200keV, bare and Eq7 Ue0.84/1.4keV; q4/low/high electronic loss with unchanged ZBL nuclear stopping and the ICRU-law baseline. Lower drag means upper G. The inherited Target cutoff remains 0.1keV.
- `tables/decision.json`, `control-analysis.json`, `verification.json`: crossing, break-even, bridge and exact replay receipts.

The q4 half-source-ICRU stopping boundary is 5.000–5.035keV lab; the observed seven-control envelope brackets it within 4.781–5.303keV. This is an operational factor-two boundary, not exact zero stopping. Finite-eta tails remain. Full stopping accuracy is NOT met everywhere: the observed envelope reaches 1.652 times max(15%|Se|,0.01eV/angstrom) near 0.717keV.

The literal Li+H compilation overlap misses 30% at -37.43%; the separate untuned Geant4 model is -12.85%. Neither is normalized. The primary G table uses the explicit ICRU model above the computed 25.1keV range; all physical returns are conditional. Literal-compilation alternatives and inherited-law comparisons remain separately labelled in the full return JSON files.

## Scope and provenance

Relaxed ordered PBE LiD, frozen Li1s core, one occupied band, full microscopic inverse before diagonal loss contraction. LiD electronic response is reused for LiT only at the same BO volume, without zero-point/isotope corrections. The error band is an executed-settings L1 envelope, not statistical coverage or an infinite-basis bound; its electronic lower endpoint is constrained nonnegative.

Saved arrays retain every sampled Q and frequency needed by the delivered observables. `unbundled-raw-matrices.json` lists original full-inverse NPZ hashes and sizes, retained by the producer, NOT included or represented as local files here. This transport scope was approved in Fleet2225; no multipart assets. Replaying these projected observables does not rerun GPAW, invert epsilon again, or prove continuum convergence.

Sources: Campillo/Pitarke/Eguiluz cond-mat/9808113v1; Yan et al. arXiv1104.1273v1. D-T coefficients are Bosch/Hale, Nuclear Fusion32(1992), pp620–621 Eqs8–9/TableIV; below 0.5keV centre-of-mass they are explicitly model extrapolation. The software ICRU/Ziegler comparator is frozen Geant4v11.3 G4BraggModel, not an unseen original ICRU table. ZBL is the unchanged Coherencev1@47f0902 forward model. Literal IAEA compilation rows preserve unavailable Ep1995/Ph1953 primary-edition gaps.

`manifest.json` hashes all bundled inputs and code. Private original paths were removed from derived metadata; original hashes remain provenance, not hashes of a rewritten original.

## Lead's note

Repository copy: the saved loss arrays `data/*/loss-*.npz` (about 26 MB), the three full return-curve files `tables/{q4,low_drag,high_drag}/return-curves.json` and `data/expected-q4-returns.json` (about 7 MB each) are not committed. They are in the delivery archive `g6-gap-delivery.tar.gz` (SHA-256 3f94d3494a9e22c9b21aa52a82f0b2cc8d137dc723ccfc44411e5be7028dc19e) on the draft release `g6_gap-delivery`, and their hashes are in `manifest.json`. Extract them into this directory to replay. The lead's replay (NumPy 1.26.4, SciPy blocked) reproduced every stopping curve, both control envelopes and every q4 energy-return value with zero difference, and `tables/` byte for byte.

What it means for the energy return:
- The gap leaves the best return where it was: 0.362 percent bare and 0.379 percent at U_e 1.4 keV, both near 180 keV, the same as the ICRU-law baseline. Drag there is unchanged. (The 0.47 percent quoted in the node brief was a rough lead estimate with a simpler high-energy stopping law; this table supersedes it.)
- Below 13 keV the return rises 1.4 to 2.9 times over the ICRU-law baseline, most near 2 keV.
- At 1 keV lab, break-even needs 52 b in LiT against 206 b in TiT1.7 on the same footing (ICRU law and ZBL from `g4_beam/forward/model.py`). In the Eq 7 form D-T reaches 52 b at 1 keV lab with U_e of about 107 keV. TiT1.7's 206 b lies above the largest value the Eq 7 form gives at that energy (159 b), so in a metal hydride no screening value reaches break-even at 1 keV.
- The remaining drag in LiT at 1 keV is 62 percent nuclear stopping on lithium, 22 percent nuclear stopping on tritium and 16 percent electronic. A gapped host in which tritium is the only scatterer would have a break-even floor of 11.6 b at 1 keV, reached in the Eq 7 form at U_e of about 42 keV. Solid T2 and other insulating tritides go on the list.
