# p2_refit: the chain with two real inputs in place

The pass 2 inputs from `p2_input/` replace two placeholders in the chain: heat capacity 4.1e6 J/(m^3 K) (range 3.5e6 to 4.2e6) and the computed harmonic gap 0.02707 eV in place of the measured-peak value 0.037 eV. Everything else is held where it was, then the long-time heat and the sensitivity ranking are recomputed. Run from the repository root:

    python3 p2_refit/refit.py

It reads `inputs.json`, `p0_chain/pass0.py`, `p1_rank/sensitivity.py` and `p1_rank/ranking.json`, and writes `refit.json`.

## What was held fixed

Mode energy 0.037 eV. The near-threshold case keeps its original `eta = 0.5312`; it is not re-centred on the new threshold. Per-site rate 9.1e-11 /s, `Q` = 25 MeV, heat fraction 0.8 and the time constants are unchanged. The computed gap comes from the fully loaded two-atom cell and is inserted into the loading 0.8 model.

## How far the output moved

Long-time heat power density, W/m^3:

| Case | Before | Heat capacity only | Heat capacity and gap | Change |
|---|---:|---:|---:|---:|
| `eta = 0` | 5.984e6 | 7.549e6 | 5.725e6 | -4.3% |
| `eta = 0.5312` (was just inside threshold) | 1.056e4 | 1.427e4 | 0 | switched off |
| `eta = 1` | 0 | 0 | 0 | none |

With the lower gap the oscillators are effectively hotter (`theta` rises, `tanh(1/(2 theta))` falls to 0.480), and the case that sat just inside the threshold is now outside it: `1/rho_c^2 = f w - eta` becomes negative, so there is no threshold at any density for that `eta`. At fixed density the critical `eta` moves from 0.537 to 0.403. This is the same threshold law evaluated at a new input, not a new form.

The `eta = 0` temperature moves from 2295 K to 1696 K, still far outside the range where the material inputs apply.

## Domain quantities

| Quantity | `eta = 0`, before to after | `eta = 0.5312`, before to after |
|---|---|---|
| `rho_c` at 300 K | 1.276 to 1.443 | 3.472 to none |
| ordering energy at 300 K, meV per site | -108.6 to -81.3 | -0.263 to 0 |
| ordering energy at the thermal fixed point, meV per site | -51.8 to -36.3 | -0.0086 to 0 |

## Ranking after the replacement

Decades of long-time heat over each unknown's range, `eta = 0` case:

| Unknown | Before | After |
|---|---:|---:|
| `k_site` | 9.69 | 9.67 |
| `tau_th` | 3.27 | 3.29 |
| `eta` | 2.75 | 1.82 |
| `f` | 1.12 | 1.61 |
| `z_eff` | 1.12 | 1.61 |
| mode energy `E_c` | 2.14 | 1.10 |
| `Q` | 1.04 | 1.02 |
| heat fraction | 0.47 | 0.46 |
| loading | 0.25 | 0.26 |
| species (mass only) | 0.23 | 0.24 |
| bath temperature | 0.05 | 0.08 |
| `Cv` | 1.69 | 0.06 |
| bare gap | 0.19 | pinned |

Heat capacity falls from sixth to near the bottom: its earlier weight came from the width of the placeholder range, not from the physics. The gap shows no span because it is pinned to the one computed value; that does not mean it is settled, since the calculation is 27 percent below the measured peak.

In the case that was near threshold, every unknown except `eta` now leaves the output at zero, so no span can be stated for them.

## Coarse laws

No new fit was made. The exact threshold law recovered in `p1_pysr/` is re-evaluated at the new inputs. The feature functions of the reduced equation from `p1_sindy/` are re-evaluated in `coarse.py`, but these material points lie outside the range that equation was fitted on, so those values are diagnostics and not predictions. One new number supplies no residual to fit a correction to.
