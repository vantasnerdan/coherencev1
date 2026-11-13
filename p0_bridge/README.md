# p0_bridge: the bridging number

The chain run backwards: from a reported excess-power density to the rate per site that bridging would have to deliver. It is a conditional requirement, not a measured nuclear rate and not a verdict on the claims. Run `python3 bridging_number.py > bridging_number.json`.

## Formula

    rate per coherent site = P / (C * n * Q)

`P` reported power density, `n` number density of H or D, `Q` energy per event (assumed fully converted to heat), `C` the fraction of occupied sites that are coherent. Every rate below is for `C = 1` and scales as `1/C`.

## Result

Reference value: **9.1e-11 per site per second**, for 21.4 W/cm^3, `Q` = 25 MeV, Pd/D at loading 1.0. Envelope over both claims, `Q` = 1 to 25 MeV and all tabulated densities: **3.1e-13 to 3.8e-9 per site per second**. The envelope is a parameter range, not a measurement uncertainty.

| Claimed P (W/cm^3) | Q (MeV) | Events / (m^3 s) | Rate per site per s, all densities | Rate per site per s, Pd/D only |
|---:|---:|---:|---:|---:|
| 0.095 | 1 | 5.93e17 | 7.7e-12 to 1.7e-11 | 1.0e-11 to 1.7e-11 |
| 0.095 | 25 | 2.37e16 | 3.1e-13 to 6.7e-13 | 4.0e-13 to 6.7e-13 |
| 21.4 | 1 | 1.34e20 | 1.7e-9 to 3.8e-9 | 2.3e-9 to 3.8e-9 |
| 21.4 | 25 | 5.34e18 | 6.9e-11 to 1.5e-10 | 9.1e-11 to 1.5e-10 |

## Energy per site against energy per event

The domain model supplies an ordering energy of 1 to 10 meV per site with the diamagnetic term off, and none with it on.

| Q | Q / (ordering energy per site) |
|---:|---:|
| 1 MeV | 1e8 to 1e9 |
| 25 MeV | 2.5e9 to 2.5e10 |

This is arithmetic: the number of sites whose ordering energy adds up to one event. It is not a domain size, a participant count or a mechanism.

## The claims

Both power densities are claims from Fleischmann, Pons and Hawkins, J. Electroanal. Chem. 261, 301 (1989), p. 304, Table 1, for Pd/D. They are normalized to the full electrode volume. If only a fraction of that volume is active, the local power density and the required rate rise without an upper bound, because the source gives no lower bound on the active fraction. They are not independently verified calorimetry.
