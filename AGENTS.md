# AGENTS.md

Working rules for agents on this repository. They come from the first campaign (passes 0 to 3) and a review of what slowed it down. Read this before starting work.

## Goal

Describe the ideal lattice for a coherence domain, a region where many hydrogen nuclei move in step, and then find a way to build it. No single metal is the target. Palladium and nickel are reference points because measurements exist for them.

This is a building project. Published theory is an input. Do not frame work as testing, validating or refuting anyone, and do not add gates, disclaimers or sections that narrow claims. Nothing here is evidence that excess heat or nuclear effects exist.

## Method: one model, two layers

Keep two layers apart.

- **Lattice layer.** Knows nothing about which material it is. It takes a handful of numbers and says whether the hydrogens lock in step. `g3_lattice/` holds it, with the single-site oscillator in `g2_refit/` and the reruns in `g3_lattice/rerun/`.
- **Material layer.** Turns a real compound into those numbers. `g2_input/`, `g3_loading/` and `g3_skew/` hold the first-principles workflow for NiH and PdD.

A new idea is either a new point in the lattice layer's input space or a new axis. Decide which before doing anything.

### Speak in ratios

Describe every case by these ratios, not by the material's name:

1. neighbour coupling over well stiffness
2. departure of the well from a simple spring, at the size of the zero-point motion
3. vibration energy over thermal energy
4. vibration frequency over loss rate
5. pump depth over pump threshold
6. fraction of empty sites, and the shift one empty site causes compared with the coupling
7. screening length over lattice spacing

The ideal is a region in this space. The map of that region, with a fitted boundary law, is the shared baseline. NiH and PdD are its two known points.

### New point or new axis

- **New point** (cheap): strain, isotope, a different metal, alloying, loading, temperature, double-bottomed wells (negative stiffness). Get the fingerprint and place it on the map.
- **New axis** (extends the lattice layer; needs a stated reason): slow long-range modulation of stiffness (for example near a Curie point), long-range coupling in insulating hydrides, four-cornered sites with a different neighbour geometry.
- **Outside the map**: how the pump is delivered (beat notes, antennas). These do not change the lattice.

### Fingerprint first

Place a candidate on the map with the cheapest numbers available, in this order: published measurements (neutron vibration energies give the stiffness), then a two-atom-cell calculation (minutes; it agreed with the 64-atom cells to between 1 and 11 percent here). Only a candidate near the target region gets large cells or long runs.

The loss rate and the pump coupling are unknown for every material. Work that pins either one improves the whole map and outranks any single new material.

### Coarse to fine

Model the whole chain in its cheapest faithful form, then add resolution where the output is most sensitive. Size work by compute, not by calendar.

- "Cheapest" means the cheapest faithful coarsening of the thing being modelled. Before any code, write down what the model keeps and what it drops.
- Each pass names the observable that would settle its question and measures it in the earliest model that can, even crudely.
- Each pass lists the three biggest things the model leaves out.

## Rules for computed inputs

Coarse-to-fine applies to physics resolution. It does not excuse an unsettled number.

- **Converge before anything consumes it.** Any first-principles number gets at least a two-point check of the dominant numerical setting (for a metal, the k-point mesh) before another model reads it.
- **Cross-check by a second route.** On-site stiffness minus summed neighbour coupling must agree with the two-atom-cell collective stiffness within a declared tolerance, at matched settings. Use any such identity as a standing check.
- **Agreement with experiment is not a convergence check.** The first NiH well matched measurement to 2 percent with a coupling three times too small.
- **Keep an inputs ledger.** One table: each consumed number, the settings it was computed at, its convergence status, and what reads it. Re-read the row when a number becomes load-bearing.
- **Relax positions.** Several soft or negative stiffness values here are saddles because atoms were held on ideal sites. Compute wells where the hydrogen actually settles.
- **Who decides what.** The lead sets budget, priority and the accuracy needed. The specialist chooses the numerical settings that deliver it.

## Rules for records

- A finished result is a record on its node and a commit, at the time it is computed. Messages carry decisions and pointers only. Batch reruns, not reporting.
- When a result supersedes an earlier one, edit the earlier README in the same commit with one line pointing forward. Keep the top-level README current.
- Keep checks that test the number being produced. Drop ritual: refusal tests, receipts, cost ledgers to many digits.
- No compute accounting. No core-minute budgets, splits, reserves or admission pilots, and no CPU, memory or invocation figures in messages or pages. Run what the result needs; stop a job that is stuck or out of proportion to the question and say why in one sentence.
- One message shape: the result, one line of assumptions, the decision needed if any.

## Rules for scope and reporting

- Define the stopping point as a deliverable and fix it when it is requested. New ideas go on a list, not into scope.
- When work is added to the critical path, say how much waiting it adds and offer the option of shipping now with the open item marked.
- Spend compute in a stated order of value.
- A question from the person directing the work is a question. Answer it. Offer the work in one line; do not queue it.
- Report to that person only on results, decisions or problems, in plain prose, from recorded results. Do not relay interim numbers. Keep one table of current best numbers and update it.

## Repository

- Public. No PDFs, credentials, machine paths, personal names other than cited authors, or email addresses.
- One writer. Specialists deliver hash-bound archives; the lead verifies hashes, reruns what it can, and commits.
- Each directory holds a small script, its output and a README. `p*` directories are the first loop and are kept for the record. `g*` directories are the current model.
- NumPy and SciPy cover everything except the first-principles runs (GPAW), the EMT geometry (ASE), and two first-loop tools (PySR, PySINDy).

## Open items carried forward

1. Map the ratio space with the existing lattice simulation and fit the locking boundary. This also explains why NiH locked and PdD stayed patchy.
2. Relax the hydrogens in representative strained and partially loaded cells.
3. Measure or estimate the loss rate of the hydrogen mode, and the stiffness modulation per unit of drive.
4. Matched NiH and PdD lattice runs, at larger size and at 300 K.
5. Converge the 64-atom k-point mesh (PdD collective value is 11 percent above the two-atom result) and recompute the NiH 111 quartic.
6. Second-neighbour vacancy force, in a cell that does not cancel it.
7. A first-principles check of one soft interface region with hydrogen present.
