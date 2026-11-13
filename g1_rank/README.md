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
| NiH | 23.5 THz | 8.4 GHz | 4.8 nm | 1.6 ps |
| TiH2 | 39.5 THz | 21.4 GHz | 4.1 nm | 0.9 ps |
| PdD | 9.7 THz | 3.4 GHz | 5.3 nm | 4.1 ps |

## Which pump (`pumps.py`)

    PYTHONPATH=../g0_model python3 pumps.py > pumps.json

| Pump | Threshold | Reachable |
|---|---|---|
| Steady or pulsed current | electron drift must exceed the mode's phase velocity: above 1e14 A/m^2; a 1e13 A/m^2 pulse reaches 3 to 10 percent of that | no |
| Hot protons arriving during loading (PdH) | arrival flux above 6e30 per m^2 per second for one atomic layer, equivalent to about 270 bar of H2 striking with every atom sticking hot; about 5000 bar for a 5 nm layer | no, and heating a harmonic mode does not invert it |
| Resonant light at the mode frequency | none: a forced response at any strength. At the zone-centre frequencies (21.5, 36.0, 8.9 THz for NiH, TiH2, PdD) an internal field of 2 to 11 MV/m gives a 0.02 Angstrom amplitude | yes, as a probe |
| Parametric drive at twice the mode frequency | well stiffness modulated by 1.6 to 1.8 percent (a strain near 0.3 percent) | yes, with existing mid-infrared pulsed lasers |

The parametric drive is the one carried forward. The mode locks to half the pump frequency (pump at 46.9 THz for NiH, 78.7 THz for TiH2, 19.4 THz for PdD), with a threshold, saturation and two possible phases. What is still needed is how strongly a pump field or strain changes the well stiffness in the loaded hydride.

## What an experiment could see

- Static branch: a proton vibration that softens as loading approaches the threshold, then a static distortion, with a mode that hardens as the square root of the excess loading.
- Pumped branch: under a pump at twice the mode frequency, a threshold in pump strength, then oscillation at exactly half the pump frequency that saturates, with picosecond decay when the pump stops.
- Isotope: at the same host, the D frequency is close to the H frequency divided by root 2, while the static threshold and ordering energy do not change.

## Inputs to compute next

The screened proton-proton force constant at finite wavevector (decides the static branch), and how much a pump field or strain changes the proton well stiffness (decides the pumped branch).
