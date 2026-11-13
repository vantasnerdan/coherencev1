# p1_sindy: reduced equation for the order parameter

PySINDy 2.1.0 fit on the 200 full five-variable trajectories in `p1_sweep/timeseries.npz`. Copy `timeseries.npz` beside the script and run:

    python3 sindy_reduction.py

The result below was produced on the author's machine and has not been re-run elsewhere.

## Result

    dq/dtau = 0.970 alpha q + (1.343 alpha - 0.964 beta) q^3

`q` is the field amplitude and `tau = omega_c t`. The coefficients are functions of the dimensionless groups:

    w     = tanh(1/(2 theta))
    A     = (f w - eta) rho^2 - 1                      distance from threshold
    D     = 2 kappa + f rho^2 w gamma / (r^2 + gamma^2)
    B     = (2 f^2 rho^4 / r) (w - (1 - w^2)/(2 theta))
    alpha = A / D,   beta = B / D

The order parameter grows when `alpha > 0`, which is the same threshold the full model has, and saturates at `q^2 = 0.970 alpha / (0.964 beta - 1.343 alpha)`.

## What was supplied and what was fitted

`alpha` and `beta` were supplied as candidate features, derived from the slow-response expansion of the pass 0 equations. That expansion suggests `dq/dtau = alpha q - beta q^3`. Six candidate terms competed (`q`, `q^3`, `alpha q`, `alpha q^3`, `beta q`, `beta q^3`); one sparse fit, threshold 0.05, run once. It kept three terms with coefficients close to but not equal to the expected ones, and added the `alpha q^3` correction. So the fit confirms the form and measures the coefficients; it did not find the form unaided.

## Scope

Fitted only near threshold (margin 0.002 to 0.015 on either side), where one slow real mode is well separated from the fast ones, and only after the fast variables have relaxed: 72,036 samples. It is a local reduction, not a model of the oscillatory field and spin dynamics in general.

## Reduced equation against the full model

Integrated from the first post-relaxation point to the end of each trajectory. These four runs were also in the fit.

| Run | Side of threshold | Full q at end | Reduced q at end | Largest error along the way |
|---:|---|---:|---:|---:|
| 0 | below | 4.9e-12 | 2.3e-12 | 1.2e-7 |
| 1 | above | 0.03253 | 0.03259 | 7.0e-4 |
| 100 | below | 1.1e-12 | 2.2e-12 | 1.3e-7 |
| 101 | above | 0.09218 | 0.09222 | 7.0e-3 |

Growth, saturation and decay are reproduced. For run 101 the largest transient error is about 7.6 percent of the saturated amplitude, so the reduction is coarse for transients and close for end states.
