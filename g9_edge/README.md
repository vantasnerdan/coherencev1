# g9_edge — the finite arrival edge

**Model:** a charged deuteron is absent for t<0, appears at z=0 at t=0, and enters z>0 at fixed incident speed. Identical gates in the bulk and entry branches isolate the surface term; this is conditional on that atom-source convention, not a simulated approach from vacuum. Pure-Pd bulk RPA is a **conditional proxy** for the four unchanged Pd/PdD/Ti/TiD transported profiles, never measured foil or hydride response. No store coupling is fitted.

**Electron result:** `output/curves.csv` gives induced potential energy q_D phi (not half the self-polarization energy), at the ion and 0.1, 0.5, 1 and 2 Å ahead. `peaks.csv` separates the initial screening deficit from the first signed opposite excursion, with peaks, times and finite lobe widths for all five CM speeds. The reference is the frozen/adiabatic response of the **same causal spectral realization** at the same geometry; the separately labelled fully formed moving **bulk** wake has infinite source prehistory and is not silently substituted for it. `edge_windows.csv` reports a source-moment electronic clock, 1 Å transit/period ratio, and the initial-plus-first-opposite-lobe window; it is not an excitation lifetime.

At 0.25 keV CM and R=0.1 Å, the entry potential starts **+25.02 eV** above frozen, then reaches **−16.82 eV** at 0.1135 fs (opposite lobe width **0.07264 fs**); the identically gated bulk overshoot is **−2.886 eV**. The realized electronic moment is 32.95 eV, T=0.1255 fs. Its sparse-grid static mismatch is retained, not normalized away.

| CM keV | transit/T | first-window end fs | pair barrier deficit / next peak eV | raw p yield / smooth40: Pd, PdD, Ti, TiD |
|---:|---:|---:|---:|---|
| .25 | 3.640 | .1468 | −13.62 / +6.444 | .2773, .2820, .7588, .7587 |
| .5 | 2.574 | .1469 | −13.72 / +6.441 | .6183, .6247, .9506, .9505 |
| 1 | 1.820 | .04118 | −13.74 / +.0953 | .8458, .8512, .9914, .9913 |
| 2.5 | 1.151 | .05136 | −13.75 / +.6923 | .9639, .9674, .9992, .9992 |
| 6.5 | .7139 | .1526 | −13.75 / +4.697 | .9924, .9948, .99993, .99993 |

These are **conditional gated-source ratios**, not material rankings. Merely replacing the exact t=0 sheet sample by the separately computed formed bulk wake gives Pd ratios **1.0172, 1.0108, 1.00535, 1.00216, 1.00062**; the sheet/prehistory choice dominates the low-energy result. The prepared-store release remains an unknown separate term, not zero inferred from this column.

**Receivers/stores:** the field is differentiated at a stationary target (z=2.042 Å) plus its twelve fcc neighbors, never along a co-moving probe. At .25 keV CM the target's peak impulse is 0.1139 eV·fs/Å, kick energy 31.1 µeV; the simultaneous 13-site total peaks at 0.365 meV. `fixed_sites.csv` contains their signed time histories and pair relative energies. Target impulses are causal-minus-frozen; the projectile subtracts formed bulk drag already in CSDA. `stored_energy.csv` retains all requested stores. Cold pumped g3 PdD holds 0.602 meV/D, 19.7 eV across its **whole 32,768-D box**; 1.377 nm phase length is not a population. Impact-accessible store release, coupling and receivers remain unknown; equilibrium cohesion/clouds/ZPE are not excited energy.

**Rate accounting:** `yields.csv` keeps raw p/n yields, settled outside-edge and transient inside-edge subtotals, and correctly normalized shares. Species, spread, transported profile, exact t=0 sheet, stopping and smooth40 Bosch–Hale remain fixed. At each local collision energy/depth, use the sharp first window, signed pair barrier difference **inside 1 Å**, and separately computed relative kinetic correction **before 1 Å**: no duplicate field work or averaged U. Lowest nominal R_turn=0.0576 Å; halving it changes the pair correction by at most .279 eV. Details and the ballistic first-spacings clock z/v_incident are explicit in [replay.md](replay.md).

**One non-yield observable:** the conditional 3.02-MeV proton-line recoil shift along the impact axis. For the central target at .25 keV CM it ranges **−14.54 to +16.89 eV** during the chosen edge; `fixed_sites.csv` gives its speed dependence, not an emission efficiency.

**Exactly three largest omissions:** (1) microscopic nonlocal half-space and nonlinear projectile charge-state/physical-prehistory dynamics; (2) impact-accessible excited-store population and its release/receiver coupling; (3) self-consistent moving nuclei and time-dependent tunnelling.

**Numerical/model sensitivity:** unchanged-profile denominator agrees to 8.9e−16. Step .001→.0005 fs changes pair correction .0213 eV and selected local rate .00256. Local-image kcap1→2 Å⁻¹ changes a lowest-speed near-ion curve up to **14.41 eV** and total yield ratio .00641; this is not exact nonlocal-surface ab initio response. Final replay used 94.2 CPU s, 534 MiB RSS; total including the corrected first replay about 210.3 CPU s, single-threaded. Equations/replay: [replay.md](replay.md); checks: `output/check.json`. Stop here: no new campaign or profile changes.
