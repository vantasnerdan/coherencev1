# What the existing radial bytes do and do not join

The retained `g6_collision/curves/baseline_BO_pair.json` supplies seven target-conditioned **total interaction** samples at R=0.1–1.6 Å. Its source is a frozen ordered conventional PdD cell, a head-on [110] projectile checkpoint from a 2-keV track, and occupied/vacant neutral-host subtraction. It is not an isotropic pair potential, a prepared bound pair, or the joint path of two initially lattice-bound deuterons. No supplied origin/contact or exterior continuation belongs to that file.

The present full G10 register has nearest D–D spacing 2.8878240944 Å. Even the generous radial diagnostic |u_A|,|u_B|≤0.24 Å leaves R≥2.4078240944 Å: **at least 0.8078240944 Å separates this diagnostic domain from the largest existing radial sample**. The smallest cold separation actually sampled here is 2.735134 Å, farther away. All 2-Å observation-shell crossings are outside that full-well diagnostic; partial cells have only ±0.03-Å measured probe-axis support. There is therefore neither distance overlap within the currently measured lattice domain nor a common full coordinate/preparation chart with the single-ray BO samples. A scalar 40-eV Eq.7 shift specifies neither missing force surface.

The forthcoming same-geometry connected pair cut can extend a distinct, properly registered approach path. For original PdD indices A=25, B=63, zero-image [110] unit vector n and displacements u_A=+u n, u_B=−u n, R=R0−2u:

- W(u)=E_AB−E_A−E_B+E_0.
- W′(u)=−n·(F_AB,A−F_A,A)+n·(F_AB,B−F_B,B). The all-index mixed-force difference is a different observable.
- The old Hamiltonian already contributes +J_AB u² on this path. Add **W(u)−J_AB u²**, with force derivative W′−2J_AB u, or replace that bond consistently. Do not add W and the old bilinear term twice. Derivatives with respect to approach distance d=2u are half the u derivative.

A joint cut reaching u≈0.644 Å would reach R=1.6 Å, and u=0.8 Å reaches R≈1.288 Å. Those are coordinate reach estimates, not received force data or permission to splice unmatched preparations. Wide *one-body* wells alone do not supply that joint electronic response. A supplied connected path can repair the pair force/energy observable only on its actual geometry/support; it cannot by itself supply the remaining three-dimensional angular/contact state, continuum incident flux, nuclear-scale continuation, mode population or capture fraction.

## Rate meaning

The recorded 2-Å flux is a joint r,p observable, not a scattering matching surface. For a genuine incoming ensemble one would need its inward phase-space flux **and** a transmission/contact map for the same inner Hamiltonian. Only at a valid scattering match may an isotropic impact-parameter ensemble be contracted with the ledger's separate smooth σ_BH40(E)v. Instantaneous lattice ΔK is not that asymptotic incident energy.

For a normalized bound allowed-spin s-wave pair, the existing `g6_pairlaw` coefficient instead gives Γ_p=A_pρ_contact, A_p=5.81220367099×10^7 Å³/s (random nine-spin ensemble: 3.87480244733×10^7 Å³/s). No arbitrary attempt frequency enters. `rates.csv` reports conditional linear conversion coefficients for an unknown *event-integrated excess* contact Ξ and an unknown stationary summed contact C_domain. These coefficients are **not rates of the sampled trajectories**. Per-event contributions subtract the simultaneous matched equilibrium/no-impulse trajectory once; a stationary domain rate uses an explicitly supplied active-domain count and divides by true current only as a conditional per-coulomb conversion.
