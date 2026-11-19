# coherencev1

A coarse end-to-end model of the coherence-domain hypothesis for hydrogen-loaded metals (Ni/H, Pd/D), refined in passes.

The aim is a model that maps apparatus-level inputs to measured outputs, so there is something to engineer against. It does not set out to show the effect exists, and nothing here is evidence that it does.

## Method

Model the whole chain in its cheapest form first, then add resolution where the end-to-end output is most sensitive. Passes are sized by compute.

- **Pass 0, negligible compute:** a minimal domain model (matter coupled to one field mode, a few ODEs), looked-up material scales, and one-line placeholders for the rest of the chain.
- **Pass 1, light compute:** sweep the pass 0 model, compress it into closed-form laws (PySR, SINDy), rank the unknowns by sensitivity.
- **Pass 2, moderate compute:** replace the top-ranked placeholders that a tool or a lookup can supply with real values; re-run and see how far the output moves.
- **Pass 3 and on, heavy compute:** converged runs where the ranking points.

Passes 0 to 2 are complete. Each directory holds a small script, its output and a README.

| Directory | What it is |
|---|---|
| `p0_scales/` | Looked-up material scales for Ni/H and Pd/D |
| `p0_domain/` | Minimal domain model and the diamagnetic-term comparison |
| `p0_chain/` | End-to-end chain: material in, heat power density out |
| `p0_bridge/` | Per-site rate the bridging step would have to deliver |
| `p1_sweep/` | 10,000-run sweep and 200 full trajectories |
| `p1_pysr/` | Symbolic regression on the sweep |
| `p1_sindy/` | Reduced equation for the order parameter |
| `p1_rank/` | Sensitivity ranking of every unknown |
| `p2_input/` | Heat capacity by lookup, bare gap by one GPAW calculation |
| `p2_refit/` | The chain and ranking with those two values in place |

## What one loop found

- **The minimal model orders only if the diamagnetic term is evaded.** With that term off, a coherent ordered state appears above a density threshold, and Ni/H and Pd/D at loadings 0.6 to 1.0 sit above it. With the term at the value the sum rule fixes, there is no threshold at any density and the end-to-end output is zero. The model carries the term as a knob `eta` from 0 to 1; a threshold exists only for `eta < f tanh(1/(2 theta))`.
- **The model gives a threshold, an order parameter and an ordering energy per site. It does not give a locked frequency, a domain size or a lifetime.** The ordered state is static, the model has no geometry, and the ordered equilibrium does not decay.
- **The bridging requirement.** Working back from claimed excess-power densities for Pd/D, the bridging step would have to deliver about 3e-13 to 4e-9 events per coherent site per second. One nuclear event equals the ordering energy of 1e8 to 2.5e10 sites.
- **What the output depends on.** The per-site rate (about 10 decades), the thermal time, and the diamagnetic knob lead the ranking. None of these is computed by any tool used here. Steady heat depends on thermal time and heat capacity only through their ratio.
- **No bursts.** Across 10,000 sampled runs every ordered case settles to a steady output.
- **Method checks.** Symbolic regression recovered the known threshold law exactly once the target was re-expressed without its singularity, and failed on the raw singular target. It found a compact approximate law for the coherent fraction and no reliable law for heat. The reduced order-parameter equation reproduces end states of the full model to a fraction of a percent near threshold.
- **Folding in real inputs.** Heat capacity by lookup collapsed a hundredfold placeholder range to about 20 percent and dropped it to the bottom of the ranking. A computed bare gap 27 percent below the measured value moved the output by 4 percent with the diamagnetic term off, and switched it off entirely in the case that had been just inside the threshold.

## Loop 2: a model that orders

The first loop's model was too simple: two-level matter and one field mode. It differs from the published coherent-plasma work as shown below (Gamberale and Modanese, arXiv:2205.14015; Gamberale, arXiv:2209.06139; arXiv:2403.15487; Modanese, arXiv:2608.21452). That work is an input here, not a target. `g0_stability/` shows that coupling to the transverse field alone, treated consistently, gives no ordered state in either model, so the second loop looks for the ingredient that does: dipole-dipole interaction between protons, anharmonic wells, coupling to the electron plasma, drive and loss, and a gradient term that gives the domain a size.

| | This repository | Published theory |
|---|---|---|
| Matter | two-level systems | harmonic oscillators with about 11 bound levels in a 4.8 eV well |
| Field | one mode | all propagation directions and both polarizations |
| Diamagnetic term | a free knob `eta` | kept in full; ordering is attributed to the all-directions mode sum |
| Energy scale | 0.037 eV optic phonon | 0.41 eV dressed oscillator, proton plasma frequency 0.2 to 0.3 eV |
| Gap per particle | 0.1 eV or less | 0.37 to 1 eV, far above thermal energy |
| Electrons | absent | a second coherent plasma; capture uses the electron density at the proton |
| Nuclear step | free per-site rate, steady, no drive | electron capture with a 782 keV threshold, rate proportional to N squared per domain, needs pumped excitation |
| Domain size | none | 1 to 10 micrometres, 1e10 to 1e13 protons |
| Material | Pd/D central case | Ni/H and Ti/H; Pd/D is stated to be a different mechanism |

One point of agreement: the published paper's single-direction comparison (its Eq. 32) shows no instability, which is what this model finds at `eta = 1`. The two-atom GPAW cell is a periodic crystal and is adequate for a phonon frequency, but that frequency is not an input the published theory uses.

## Pass 3: computed wells, loading, strain and a lattice of up to a million sites

Work after the table above, in directories `g0_*` to `g3_*`, each with its own README:

| Directory | What it is |
|---|---|
| `g0_stability/`, `g0_model/`, `g0_scales/` | Field coupling alone gives no order; a small model that tests which ingredients do |
| `g1_rank/` | Sweep, coarse laws and ranking for that model |
| `g2_input/` | Hydrogen wells from GPAW: two-atom cells, strain derivatives |
| `g2_refit/` | Anharmonic oscillator on the first computed wells (inputs since superseded) |
| `g3_loading/` | Corrected single-atom wells, partial loading, a loading law, NiH check |
| `g3_skew/` | Ni on Pd mismatch pattern; PdD well under shear and thin-film strain |
| `g3_lattice/` | Lattice simulation: loading, arrangement, strain and pump; phase coherence; loss |

Status: the single-atom wells in `g2_input` and everything computed from them in `g2_refit` and `g3_lattice` used settings that turned out to be too coarse. `g3_loading` holds the corrected wells. The oscillator and lattice results are being rerun on them; until that lands, treat the numbers in `g2_refit` and `g3_lattice` as provisional.

## Pass 4: a second observable on the material layer

| Directory | What it is |
|---|---|
| `g4_beam/` | Forward model of low-energy D-D beam-fusion yields in Pd, PdD, Ti and TiD (Karahadian et al. 2026), near-surface deuterium transport with beam damage, measured screening potentials. Outside the lattice map; shares ratio 7 (screening length over lattice spacing) and the near-surface loading with it. |

## Pass 6: resolving the lumped inputs of the beam and pair models

Pass 6 replaces seven lumped inputs with resolved quantities: the fitted screening shift, the low-energy excess, random-medium stopping, unrelated sources for screening and drag, gapless stopping in insulators, the missing lattice-to-beam coupling, and the pair-at-rest rate law. Directories appear as the nodes deliver.

| Directory | What it is |
|---|---|
| `g6_pairlaw/` | Fusion rate of a pair at rest from a normalized radial solution in a real potential, for D-D and D-T. A scaled D2 Born-Oppenheimer family, and a screened repulsion inside a harmonic well (screening length, well stiffness and well centre as separate levers). It reproduces the D2 molecule rate within a factor of 2 without fitting, and retires the earlier two-point interpolation. |
| `g6_couple/` | Part (d): the locked lattice state as a cause of the keV screening. The deuterons' own (ionic) response adds 0.005 eV at 0.5 keV lab, falling with speed. A frozen-spectator estimate gives 0.064 eV for an imposed 0.0432 A displacement cap, while the lattice runs reach 0.077 A (cold pumped PdD) to 0.72 A (warm PdD). The locked-minus-unlocked contrast is zero at fixed charge and poles. Reaching 0.5-2 keV would need about 600-1200 times the effective charge, so the open route is the host electrons' response to the locked state. Parts (a) to (c) to follow. |
| `g6_plateau/` | What each of four causes of the low-energy excess (rising screening, surface layer, fast neutrals, lattice state) predicts for D-T, LiD and a pair at rest, plus a keyed per-run scorer for the Zenodo data. Swapping D for T in the target at 0.5 keV lab gives 55.6 times the yield if fast neutrals cause the excess, and 2.5-3.2 times if the target does. |

## Reproducibility

NumPy and SciPy are enough for everything except `p1_pysr/` (PySR), `p1_sindy/` (PySINDy) and `p2_input/frozen_phonon.py` (GPAW). Those three were run on one machine only; all other outputs were reproduced on a second.
