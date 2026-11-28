# AGENTS.md

Working rules for agents on this repository. Read this before starting work.

## Goal

Describe the ideal lattice for a coherence domain, a region where many hydrogen nuclei move in step, and then find a way to build it. No single metal is the target. Palladium and nickel are reference points because measurements exist for them.

This is a building project. Published theory is an input. Do not frame work as testing, validating or refuting anyone, and do not add gates, disclaimers or sections that narrow claims. Nothing here is evidence that excess heat or nuclear effects exist.

## Correction, 2026-11-28

Until this date the method did three things:

- It reduced PdD to a few numbers taken from DFT cells of 2 to 64 atoms: site stiffness, neighbour coupling, and pair energies at chosen wave vectors and distances.
- It fed those numbers to a lattice model that does not know which material it is, and built every later model on top of that.
- It used no interatomic potential, and it kept models to NumPy so the lead could rerun them on one machine.

The person directing the work: "we have spent more time trying to narrow and toyify the models and simulations to the point that the results are not useful for larger dynamic predictions, so by trying not to waste time we wasted a week."

The rules below replace that method. The earlier directories stay as the record. When the full simulation supersedes a result that rests on the reduced chain, that result's README gets one line pointing forward.

## Method: simulate the real material whole

- **The base model is the field's standard full simulation of the real material.** For PdD that means:
  - a large 3D block of Pd and D atoms at the loading in question, every atom free to move;
  - forces from an interatomic potential checked against DFT;
  - a standard code (LAMMPS or equivalent);
  - the deuterium's quantum motion from path-integral MD;
  - for an impact, three further pieces: electronic stopping, the short-range repulsion, and the electrons as a field coupled to the atoms (two-temperature MD).

  Vibrations, their dispersion and decay, collective motion and close approaches come out of the run. They are not chosen in advance.
- **Continuum fields carry the scales beyond the block** (a whole foil, seconds of beam). They take their properties from the block and exchange waves with it, as in `g15_track/`. They do not replace the block.
- **First principles builds and checks the potential.** It also covers what a potential cannot, such as electrons at high temperature (`g16_hot_electrons/`) and the fast pass of a deuteron (`g6_collision/`). It does not supply a handful of numbers to a reduced model.
- **Do not narrow.** Do not reduce the material to a few parameters and build a new model on them. Do not narrow a question to a patch, a pair, a site or a single process when the full simulation holds all of them.

### Speak in ratios

Place every case on the map by these ratios, read from its full simulation:

1. neighbour coupling over well stiffness
2. departure of the well from a simple spring, at the size of the zero-point motion
3. vibration energy over thermal energy
4. vibration frequency over loss rate
5. pump depth over pump threshold
6. fraction of empty sites, and the shift one empty site causes compared with the coupling
7. screening length over lattice spacing

The ideal is a region in this space. NiH and PdD are its two known points. The lattice layer (`g3_lattice/`, `g2_refit/`) is kept as a tool for reading these ratios. It does not stand in for the simulation.

### Changing the case

A new idea is one of two things:

- **A change to the simulated system:** loading, strain, isotope, temperature, defects, surface, alloying, another metal, or a drive.
- **A process the simulation does not yet hold.** Add it to the simulation.

Decide which before doing anything. How the pump is delivered (beat notes, antennas) is a drive on the same simulation. Run Pd first. Run other metals only if the Pd result proves useful.

### Measurements check the simulation

Set published measurements beside the simulation's own output: lattice constant against loading, neutron phonon dispersion and vibration energies, elastic constants, diffusion and isotherms. Measurements check the simulation. They do not replace it.

The loss rate of the hydrogen vibration and the way a drive couples to it are outputs of the full simulation: anharmonic decay in the block, and electron coupling in two-temperature runs. Pinning either one improves the whole map.

### Coarse to fine

- **Coarse means a coarse run of the full simulation:** a smaller block, fewer shots, a shorter run or a lighter potential. It never means a reduced model of a different kind.
- **Refine by enlarging, lengthening and repeating the same simulation** where the output is most sensitive. Size work by compute, not by calendar.
- **Write down what the simulation keeps and what it drops** before writing code. Each pass names the observable that would settle its question and lists the three biggest things it leaves out.

## Rules for computed inputs

- **Converge before anything consumes it.** Any first-principles number gets at least a two-point check of the dominant numerical setting (for a metal, the k-point mesh) before another model reads it.
- **Check the potential on the configurations the simulation actually visits.** Compare DFT energies and forces on thermal snapshots, defects and close-approach frames, and the measured properties above. Say which numbers the potential misses and by how much. Where a published potential misses a number the result rests on, fit a machine-learned potential to DFT data.
- **Cross-check by a second route** wherever an identity exists, within a declared tolerance and at matched settings.
- **Agreement with experiment is not a convergence check.** The first NiH well matched measurement to 2 percent with a coupling three times too small.
- **Keep an inputs ledger.** One table: each consumed number or potential, the settings it was computed at, its convergence status, and what reads it. Re-read the row when a number becomes load-bearing.
- **Relax positions.** Several soft or negative stiffness values here were saddles because atoms were held on ideal sites.
- **Who decides what.** The lead sets priority and the accuracy needed. The specialist chooses the numerical settings that deliver it.

## Rules for records

- A finished result is a record on its node and a commit, at the time it is computed. Messages carry decisions and pointers only. Batch reruns, not reporting.
- When a result supersedes an earlier one, edit the earlier README in the same commit with one line pointing forward. Keep the top-level README current.
- Deliver the input decks, potential files and scripts that produced the result, so the lead reruns it with the same codes. Do not rewrite a simulation as a reduced NumPy reader for the sake of a rerun.
- Keep checks that test the number being produced. Drop ritual: refusal tests, receipts, cost ledgers to many digits.
- No compute accounting. No core-minute budgets, splits, reserves or admission pilots, and no CPU, memory or invocation figures in messages or pages. Run what the result needs. Stop a job that is stuck or out of proportion to the question, and say why in one sentence.
- One message shape: the result, one line of assumptions, the decision needed if any.

## Rules for scope and reporting

- Define the stopping point as a deliverable and fix it when it is requested. New ideas go on a list, not into scope.
- When work is added to the critical path, say how much waiting it adds and offer the option of shipping now with the open item marked.
- Spend compute in a stated order of value.
- A question from the person directing the work is a question. Answer it. Offer the work in one line; do not queue it.
- Report to that person only on results, decisions or problems, in plain prose, from recorded results. Do not relay interim numbers. Keep one table of current best numbers and update it.

## Repository

- Public. No PDFs, credentials, machine paths, personal names other than cited authors, or email addresses.
- One writer. Specialists deliver hash-bound archives. The lead verifies hashes, reruns with the same codes, and commits.
- Each directory holds its scripts or input decks, its output and a README. `p*` directories are the first loop and are kept for the record. `g*` directories are the current work.
- Tools: use the standard codes the question needs.
  - LAMMPS, or an equivalent MD code, with a checked interatomic potential, for atomistic dynamics.
  - Machine-learned potentials (for example MACE) where a published potential falls short.
  - GPAW for first principles; ASE and phonopy.
  - i-PI for quantum nuclei, and the rest of the toolset in `tools_marvin/` (SSCHA, NuHamil, QuTiP, quimb).
  - NumPy and SciPy for analysis and continuum fields; PySR and PySINDy.

## Open items carried forward

1. `g17_pdd_block/`: the full 3D PdD block with loading as a parameter, then deuteron shots into it.
2. The loss rate of the deuterium vibration and its coupling to hot electrons, read from the block.
3. Join the block to the foil-scale fields of `g15_track/` so that waves and heat pass both ways.
4. Rebuild on the full simulation the results that rest on the reduced chain, in order of value:
   - the deuterons' close approaches (`g10_quantum_site/`, `g15_track/`);
   - the pumped state and collective motion (`g3_lattice/`, `g13_rogue/`, `g15_track/`).

The earlier list is superseded wherever the block answers it: the ratio map from the lattice layer, relaxing strained cells, the 64-atom k-mesh, the second-neighbour vacancy force, matched NiH and PdD lattice runs, and the soft interface region.
