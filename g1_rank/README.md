# g1_rank: laws, ranking and measurable predictions for the ordering model

Uses `g0_model/domain.py` unchanged. Run:

    PYTHONPATH=../g0_model python3 rank.py > rank.json

## Static branch (dipole attraction plus anharmonic well)

Screening length the metal would need for the protons to order at stoichiometric loading:

| Material | Critical screening length | Thomas-Fermi estimate |
|---|---:|---:|
| NiH | 1.55 Angstrom | 0.56 Angstrom |
| PdD | 1.09 Angstrom | 0.59 Angstrom |
| TiH2 | none: does not order even unscreened | 0.61 Angstrom |

Closed forms: `r = k0 - x J1(lambda)`, with `J1` proportional to `(1 + d/lambda) exp(-d/lambda)`. Threshold loading `x_crit = k0/J1`; displacement `q0^2 = -r/beta`; energy `-E0 r^2/(4 beta)`; correlation length `sqrt(c/curvature)`.

Ranking, decades moved across each input's range: screening 5.0, quartic coefficient 2.0, loss 2.0, choice of vibration peak 1.6, stiffness 1.0, loading 0.7. Everything else is below 0.3.

## Pumped branch (gain)

A proton optic mode with gain `G`, loss `kappa` and saturation settles into a self-sustained oscillation once `G > kappa`. It does not need weak screening. At `G = 2 kappa` and 300 K:

| Material | Frequency | Linewidth | Gain length | Decay after pump stops |
|---|---:|---:|---:|---:|
| NiH | 23.5 THz | 7.6 GHz | 5.2 nm | 1.6 ps |
| TiH2 | 39.5 THz | 20.5 GHz | 4.2 nm | 0.9 ps |
| PdD | 9.7 THz | 2.8 GHz | 6.2 nm | 4.1 ps |

What is not known is how any physical pump turns into gain. With a placeholder gain cross-section of 1e-20 m^2 the thresholds are far out of reach: current density above 1e14 A/m^2 (the electron drift must exceed the mode's phase velocity before any gain is possible), hydrogen flux near 1e32 per m^2 per second, temperature gradient above 1e10 K/m. A steady bulk current cannot pump this mode. The gain coefficient is the missing microscopic quantity.

## What an experiment could see

- Static branch: a proton vibration that softens as loading approaches the threshold, then a static distortion, with a mode that hardens as the square root of the excess loading.
- Pumped branch: a pump threshold, a narrow line at the frequencies above whose strength grows with pump excess, and picosecond decay when the pump stops.
- Isotope: at the same host, the D frequency is close to the H frequency divided by root 2, while the static threshold and ordering energy do not change.

## Inputs to compute next

The screened proton-proton force constant at finite wavevector (decides the static branch), and the sign and size of current-induced gain on the proton mode (decides the pumped branch).
