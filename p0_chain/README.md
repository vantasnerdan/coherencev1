# p0_chain: end-to-end chain, pass 0

One call from prepared material to heat output. The domain model from `p0_domain/` followed by one-line placeholders for the remaining rungs. Run:

    python3 pass0.py chain > chain.json

`pass0.py` here supersedes the copy in `p0_domain/`: the diamagnetic switch is now a continuous knob `eta` in [0, 1], with `d = eta * rho^2 / 4`. `eta = 1` is the sum-rule value; `eta < 1` is a diagnostic, not a claimed physical mechanism. A threshold exists only when `eta < f * tanh(1/(2 theta))`, at `rho_c^2 = 1 / (f * tanh(1/(2 theta)) - eta)`.

## Rungs

- **Preparation:** species, density from loading and lattice constant, temperature, bare gap, free mode energy.
- **Bridging:** event rate per volume `R = n * C * k_site`. `C` is the squared transverse polarization of the equilibrium domain state, zero when there is no order. `k_site` is a free per-site rate.
- **Nuclear event:** energy `Q` per event. No reaction is chosen.
- **Energy disposal:** heat `h * Q * R`, radiation `(1 - h) * Q * R`.
- **Product fate:** one escape fraction.
- **Feedback:** `dn/dt = (n0 - n)/tau_refill - R` and `dT/dt = heat/Cv - (T - T0)/tau_th`. The domain is treated as instantly equilibrated.

## Free parameters

| Parameter | Example | Range | Basis |
|---|---|---|---|
| `eta` | 0 and 1 | 0 to 1 | diagnostic knob |
| `f` oscillator strength | 1 | 0 to 1 | sum-rule bound |
| `kappa`, `gamma` losses | 0.15, 0.1 | 1e-3 to 1 | guess |
| mode energy `E_c` | equal to the gap | 0.01 to 1 eV | guess |
| `k_site` | 1e-9 /s | 1e-20 to 1 /s | placeholder until the bridging number |
| `Q` | 1 MeV | 1 to 25 MeV | no reaction chosen |
| heat fraction `h` | 0.8 | 0 to 1 | definition |
| escape fraction | 0.1 | 0 to 1 | definition |
| heat capacity `Cv` | 3e6 J/(m^3 K) | 1e5 to 1e7 | unknown for the loaded hydride |
| `tau_refill`, `tau_th` | 1000 s each | 1 to 1e6 s | unknown |

## Example

Loading 0.8, 300 K start, mode energy equal to the gap, values from the table above.

| Material | eta | t (s) | T (K) | C | heat (W/m^3) |
|---|---:|---:|---:|---:|---:|
| Ni/H | 0 | 0 | 300 | 0.971 | 7.66e6 |
| Ni/H | 0 | 1000 | 1791 | 0.800 | 6.31e6 |
| Ni/H | 0 | 10000 | 2076 | 0.675 | 5.33e6 |
| Pd/D | 0 | 0 | 300 | 0.994 | 5.99e6 |
| Pd/D | 0 | 1000 | 1500 | 0.877 | 5.28e6 |
| Pd/D | 0 | 10000 | 1807 | 0.751 | 4.52e6 |
| both | 1 | any | 300 | 0 | 0 |

With `eta = 1` the output is zero for both materials. With `eta = 0` self-heating lowers the order and the output settles; no bursts appear at these values. The temperatures reached are far outside the range where the reference cells and gaps apply, so they show the placeholder's behaviour, not a material prediction.

Energy ledger: heat plus radiation equals events times `Q` to machine precision. Charge, nucleon and lepton ledgers cannot be closed until a reaction is chosen.
