# p0_scales: material reference scales for Ni/H and Pd/D

Looked-up values and SI arithmetic. No simulation. Run `python3 material_scales.py > scales.json`.

| Material | x | a (Å) | n_metal (10^28 m^-3) | n_H/D (10^28 m^-3) | E_pi (eV) | λ_pi, vacuum (µm) | E_gap (eV) | ρ* | g*, f=1 | θ at 300 K* | θ at 600 K* |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Ni/H | 0.6 | 3.732 | 7.695 | 4.617 | 0.1862 | 6.658 | 0.089 | 2.092 | 1.046 | 0.290 | 0.581 |
| Ni/H | 0.8 | 3.732 | 7.695 | 6.156 | 0.2150 | 5.766 | 0.089 | 2.416 | 1.208 | 0.290 | 0.581 |
| Ni/H | 1.0 | 3.732 | 7.695 | 7.695 | 0.2404 | 5.158 | 0.089 | 2.701 | 1.351 | 0.290 | 0.581 |
| Pd/D | 0.6 | 4.084 | 5.872 | 3.523 | 0.1150 | 10.777 | 0.037 | 3.109 | 1.555 | 0.699 | 1.397 |
| Pd/D | 0.8 | 4.084 | 5.872 | 4.698 | 0.1328 | 9.333 | 0.037 | 3.590 | 1.795 | 0.699 | 1.397 |
| Pd/D | 1.0 | 4.084 | 5.872 | 5.872 | 0.1485 | 8.348 | 0.037 | 4.014 | 2.007 | 0.699 | 1.397 |

`E_pi` is the plasma energy of the hydrogen species treated as a bare charge +e. Starred columns use one illustrative choice of mode, `E_c = E_gap` = the local-vibration peak. That choice is not the model's convention; the mode energy `E_c` is a free input.

## Mapping to the domain model's groups

With `E_c = hbar*omega_c` and `E_gap` in eV:

- `rho = E_pi / E_c`
- `r = E_gap / E_c`, detuning `= r - 1`
- `theta = kB*T / E_gap`
- `g = rho * sqrt(f*r) / 2` for oscillator strength `0 < f <= 1`
- diamagnetic coefficient under the full sum rule: `d = rho^2 / 4`
- vacuum wavelength at the mode: `1.239841984 / E_c[eV]` µm
- losses `kappa` (field) and `gamma` (matter), each over `omega_c`: guessed range 1e-3 to 1, no measured basis.

## What these numbers are and are not

- Fixed-cell loading scenarios. The lattice constant and vibration energy are held at the near-stoichiometric reference value for every loading and temperature.
- Bulk NiH at x = 0.6 to 1 is a high-pressure or metastable state, not ordinary nickel at 300 to 600 K.
- One octahedral site per metal atom, so site occupancy equals x under that model only.
- Vacuum wavelengths are reference lengths. They are not domain sizes and do not show that a metal supports such a mode.
- Host electron plasma energies (optical fits, unloaded metal): Ni 15.92 eV, intraband 4.59 to 4.93 eV; Pd 9.72 eV, intraband 5.58 eV.
- `kB*T` = 0.025852 eV at 300 K and 0.051704 eV at 600 K.

## Sources

1. Antonov et al., Europhys. Lett. 51, 140 (2000), p. 142: NiH a = 3.732 Å at 100 K; optical peak 89 meV, shoulder to 130 meV at 5 K.
2. Antonov et al., Phys. Rev. B 80, 134302 (2009), pp. 1 and 3: stoichiometric PdD a = 4.084 Å; peak 37 meV, shoulder to 65 meV.
3. Rakić et al., Appl. Opt. 37, 5271 (1998), Tables 1 to 3: host optical-fit parameters.
4. NIST CODATA 2022 constants.
5. Antonov et al., J. Alloys Compd. 580, S109 (2013): NiH phase equilibrium.
