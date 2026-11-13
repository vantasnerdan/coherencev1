# p0_domain: minimal domain model

N two-level oscillators (protons or deuterons) coupled to one field mode, mean field, counter-rotating terms kept. Five ODEs. Run:

    python3 pass0.py domain > domain.json
    python3 pass0.py a2 > a2.json

This is a hypothetical model for locating a threshold. It is not evidence that the effect exists.

## Equations

Time in units of `1/omega_c`; `q + i p` is the field amplitude per root oscillator; `s` is the spin vector.

    q' = p
    p' = -(1 + 4d) q - g s_x - 2 kappa p
    s' = (4 g q, 0, r) x s - gamma (s - s_eq(q))
    s_eq = -tanh(B / (2 r theta)) (4 g q, 0, r) / B,   B = sqrt(r^2 + 16 g^2 q^2)

Groups: `rho = omega_p/omega_c`, `r = E_gap/E_c`, `theta = kB T/E_gap`, losses `kappa` and `gamma`, oscillator strength `f`; `g = rho sqrt(f r)/2`. The mode energy `E_c` is a free input. The diamagnetic term is a switch: off `d = 0`, on `d = rho^2/4`, which is the value the sum rule fixes, not a fitted one.

## Result

- **Diamagnetic term off:** ordering appears when `f rho^2 tanh(1/(2 theta)) > 1`. At `theta = 0.25`, `f = 1` the threshold is `rho_c = 1.0185`.
- **Diamagnetic term on:** the condition becomes `f rho^2 tanh(1/(2 theta)) > 1 + rho^2`, which cannot hold for `f <= 1`. There is no threshold at any density.

Example with H, `E_c = E_gap = 0.1 eV`, 290 K, `kappa = 0.15`, `gamma = 0.1`, `f = 1` (rates in units of `omega_c`):

| rho | n (m^-3) | seed growth, off | order parameter q, off | ΔU (meV/site) | ΔF (meV/site) | seed growth, on |
|---:|---:|---:|---:|---:|---:|---:|
| 0.80 | 8.52e27 | -0.1057 | 0 | 0 | 0 | -0.0976 |
| 0.98 | 1.28e28 | -0.1012 | 0 | 0 | 0 | -0.0924 |
| 1.05 | 1.47e28 | +0.1030 | 0.1849 | -0.568 | -0.104 | -0.0906 |
| 1.30 | 2.25e28 | +0.4418 | 0.5221 | -8.643 | -6.618 | -0.0846 |

## What the model does not give

- **Frequency of the ordered state:** zero. The order is static, not a locked oscillation.
- **Domain size:** undefined. The model has no geometry. The vacuum wavelength at `E_c` (12.4 µm at 0.1 eV) is a reference length only.
- **Lifetime:** none. A bath-stabilized equilibrium does not decay in this model.

## Limits

Single mode, dipole coupling, mean field, two-level truncation, phenomenological thermal bath, bare nuclear charge and mass, no host screening, no geometry, no drive.

## Sources

- Emary and Brandes, Phys. Rev. E 67, 066203 (2003), arXiv:cond-mat/0301273: non-RWA Dicke Hamiltonian and critical coupling.
- Nataf and Ciuti, arXiv:1006.1801, Eqs. 8, 11, 13: diamagnetic term and the sum rule.
