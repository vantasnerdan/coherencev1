# p1_sweep: sweep of the pass 0 chain

10,000 seeded runs of the chain over its dimensionless groups and placeholder parameters, plus 200 full five-variable trajectories for equation discovery. Run inside this directory:

    python3 sweep.py > sweep-run.json

Writes `sweep.csv` and `timeseries.npz`. Seed 20261113. About one CPU-minute.

`pass0.py` here is the chain from `p0_chain/` with one numerical change: the slow fuel and temperature equations are integrated in coordinates that resolve a very small coherent fraction without subtracting two nearly equal temperatures. The equations are unchanged and the `p0_chain` example reproduces.

## Inputs sampled

- `rho_bare` 0.1 to 10 (log), effective charge `z_eff` in (0, 1], with `rho = z_eff * rho_bare` entering both the coupling and the diamagnetic term.
- `r` 0.3 to 3 (log), `theta` 0.05 to 2 (log), `eta` 0 to 1, `f` 0.1 to 1, `kappa` and `gamma` 1e-3 to 1 (log).
- `k_site` 1e-20 to 1 /s (log), `Q` 1 to 25 MeV, heat fraction 0 to 1, `Cv` 1e5 to 1e7 J/(m^3 K), `tau_th` and `tau_refill` 1 to 1e6 s (log). Escape fraction fixed at 0.1.
- Reference density 1e28 /m^3 and mode energy 0.1 eV convert the dimensionless model to W/m^3 and kelvin.

Rows: 200 trajectory design points near threshold, 7,840 plain draws, 1,960 extra points paired on either side of the threshold surface `eta = f*tanh(1/(2 theta)) - 1/rho^2`.

## Outputs per row

`ordered`, `rho_c`, `threshold_margin`, equilibrium `q`, coherent fraction `C`, `dU_Ec`, `dF_Ec`, `growth_per_tau`, and from the chain `heat_long_W_m3`, `T_long_K`, `n_long_over_n0`, `C_long`, `behaviour`. The long-time value is taken at 20 times the longer of the two slow time scales. There are no columns for frequency, size or lifetime because the model does not define them.

## Result

1,615 rows are ordered and all of those reach a steady output. 8,385 are below threshold and give zero. No oscillating or decaying outputs occurred. The energy ledger residual is at most 2e-16.

`timeseries.npz` holds, for 100 runs above and 100 below threshold, the time, the five state variables, the exact right-hand side, and the group values. These trajectories were chosen where one slow real mode is well separated from the others, so they suit a reduced equation near threshold and do not cover the dynamics in general.

Outputs can differ in the last digits between machines.
