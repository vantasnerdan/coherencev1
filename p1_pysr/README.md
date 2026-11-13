# p1_pysr: symbolic regression on the sweep

PySR 2.7.2 run on `p1_sweep/sweep.csv`. Copy `sweep.csv` beside the script and run:

    python3 pysr_laws.py

Each target was searched once: 20 s budget, serial, seed 20261113, operators `+ - * /` and `sqrt exp log tanh`, maximum size 18. Errors below are training errors. The fronts in this file were produced on the author's machine and have not been re-run elsewhere; a time-limited search can stop at different expressions on another machine.

Supplied, not discovered: `w = tanh(1/(2 theta))` and `G = f rho^2 / (1 + eta rho^2)`, the two coordinates the static equilibrium depends on.

## Threshold

The known law is `rho_c = 1/sqrt(f w - eta)` when `f w > eta`.

**Fitting `rho_c` directly: not recovered.** Best expression at complexity 16 has MSE 4.9 and a maximum error of 79. `rho_c` diverges as `f w` approaches `eta`, so a few near-singular rows dominate the squared error.

**Fitting `1/rho_c^2` with identical data and settings: recovered exactly.**

| Complexity | Expression | MSE |
|---:|---|---:|
| 1 | 0.254 | 4.5e-2 |
| 3 | `f - eta` | 3.3e-2 |
| 4 | `tanh(f) - eta` | 2.4e-2 |
| 5 | `f*w - eta` | 2.4e-33 |

Lesson: a target with a singularity has to be re-expressed as a regular quantity before a symbolic-search failure means anything.

## Coherent fraction C (ordered rows)

| Complexity | Expression | MSE |
|---:|---|---:|
| 5 | `tanh(log(G w))` | 9.2e-3 |
| 7 | `tanh(1.673 log(G w))` | 9.8e-4 |
| 9 | `tanh((theta + 1.400) log(G w))` | 4.0e-4 |
| 10 | `tanh((theta + 1.549) tanh(log(G w)))` | 2.6e-4 |
| 14 | nested tanh and log | 1.8e-4 |

A compact approximate law. The maximum error stays near 0.08 to 0.14, so it does not replace the exact equilibrium.

## Ordering energy, dU / (E_c r) (ordered rows)

| Complexity | Expression | MSE |
|---:|---|---:|
| 5 | `0.227 (1.306 - G)` | 1.5e-2 |
| 7 | `0.239 (theta + 0.945 - G)` | 5.8e-3 |
| 14 | `0.220 (exp(theta + 0.0429 (-0.213 - G) G) - G)` | 4.0e-3 |

Rough. Some of these give the wrong sign near onset and none should be used in place of the exact value.

## Long-time heat, as a fraction of `n Q h k_site` (ordered rows)

| Complexity | Expression | MSE |
|---:|---|---:|
| 4 | `exp(-1.438 - thermal_load)` | 5.4e-2 |
| 8 | `exp((-diamag - thermal_load - 4.371)/coupling)` | 3.6e-2 |
| 13 | `sqrt(exp((-thermal_load - 4.965)(diamag + theta + 0.316)/coupling))` | 2.5e-2 |

`thermal_load = (tau_th/Cv) n Q h k_site / T0`, `coupling = f rho^2`, `diamag = eta rho^2`. The search picks up the thermal load, which contains the ratio `tau_th/Cv`, but no compact law fits: errors in power reach 9e9 W/m^3. No reliable global heat law was found in this budget.
