# p1_rank: sensitivity ranking of the unknowns

One-at-a-time variation of every unknown in the pass 0 chain, ranked by how far the long-time heat power density moves. Run from the repository root:

    python3 p1_rank/sensitivity.py > p1_rank/ranking.json

## Setup

Central case: Pd/D, loading 0.8, 300 K, mode energy equal to the gap (0.037 eV), `f = 1`, `k_site = 9.1e-11 /s` (the bridging number), `Q = 25 MeV`, heat fraction 0.8, `tau_th = tau_refill = 1000 s`, `Cv = 3e6 J/(m^3 K)`.

Two cases, because with `eta = 1` the output is zero and nothing can be ranked:

- `eta = 0`: central heat 5.98e6 W/m^3 at 2295 K. That temperature is far outside the range where the material inputs apply.
- `eta = 0.531`, just inside the threshold (`eta_crit = 0.537`): central heat 1.06e4 W/m^3 at 303.5 K.

The long-time value comes from the chain's stationary balance, solved with the chain's own equilibrium function. It gives no time behaviour.

`D` is the number of decades between the largest and smallest positive heat over the sampled range. Zeros are recorded separately as a switch.

Kind: (a) computable with an existing tool, (b) measured, published or set, (c) open physics or unspecified geometry.

## Ranking

| Rank | Unknown | Kind | Range | D, eta = 0 | D, near threshold | Can switch output off (eta = 0 / near) |
|---:|---|---|---|---:|---:|---|
| 1 | `k_site` | c | 1e-20 to 1 /s | 9.69 | 8.48 | no / no |
| 2 | `tau_th` | c | 1 to 1e6 s | 3.27 | 4.46 | no / no |
| 3 | `eta` | c | 0 to 1 | 2.75 | 2.75 | yes / yes |
| 4 | mode energy `E_c` | c | 0.01 to 1 eV | 2.14 | 0.50 | yes / yes |
| 5 | loading | b | 0.6 to 1 | 0.25 | 2.05 | no / yes |
| 6 | `Cv` | a/b | 1e5 to 1e7 J/(m^3 K) | 1.69 | 1.96 | no / no |
| 7 | bare gap | a/b | 0.037 to 0.065 eV | 0.19 | 1.80 | no / no |
| 8 | `f` | a | 1e-3 to 1 | 1.12 | 0.05 | yes / yes |
| 9 | effective charge `z_eff` | c | 1e-3 to 1 | 1.12 | 0.15 | yes / yes |
| 10 | `Q` | c | 1 to 25 MeV | 1.04 | 0.25 | no / no |
| 11 | species (mass only) | b | D or H | 0.23 | 0.94 | no / no |
| 12 | heat fraction | c | 0 to 1 | 0.48 | 0.07 | yes / yes |
| 13 | bath temperature | b | 300 to 600 K | 0.06 | 0 | no / yes |
| 14 | `tau_refill` | c | 1 to 1e6 s | 0 | 0 | no / no |
| 15 | escape fraction | c | 0 to 1 | 0 | 0 | no / no |
| 16 | `gamma` | c | 1e-3 to 1 | 0 | 0 | no / no |
| 17 | `kappa` | c | 1e-3 to 1 | 0 | 0 | no / no |

## Reading it

- The per-site rate and the diamagnetic knob lead, and no tool here computes either.
- Steady heat depends on `tau_th / Cv` only, not on the two separately.
- The losses `kappa` and `gamma` drop out of the long-time output because the domain is treated as instantly equilibrated. That is a property of this closure, not evidence that losses are physically irrelevant.
- A span of zero with a switch (bath temperature near threshold) means only one sampled point gave positive output.

Restricted to what a tool or a measurement can supply, in order: loading, `Cv`, bare gap, `f`, species, bath temperature.
