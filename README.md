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

## Reproducibility

NumPy and SciPy are enough for everything except `p1_pysr/` (PySR), `p1_sindy/` (PySINDy) and `p2_input/frozen_phonon.py` (GPAW). Those three were run on one machine only; all other outputs were reproduced on a second.
