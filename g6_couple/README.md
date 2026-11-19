# Lattice-state fourth cause: part (d)

Run from the repository root with NumPy/SciPy:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python g6_couple/response.py
```

This is a reproducible reduced ionic calculation, not a host electronic-response calculation. It uses the actual `g4_beam.forward.consumer` integrator and unchanged transported PdD profiles/outer adsorbed sheet, a declared 40 eV electronic baseline, native stopping, 0.85/0.15 D+/D2+ electrical-current fractions, 5% spread and the fixed proton threshold term. Local slowing-down and molecular constituent energies are used in the ionic response. No fit parameters or upstream files are changed.

## Result and machine interface

`response.json.projections[]` is keyed by `intervention` and `E_lab_keV`; energies are 0.5, 1, 2, 3, 5, 8, 13 keV **lab**. It reports absolute physical yield/C factors and separately normalized R factors. Controls are `baseline`, `near_surface_D_x50` (volumetric D only in top 10 nm; outer sheet unchanged), `D2_current_0.85` (electrical fraction; two D per molecular charge), and `half_current_frozen`.

- The signed locked-minus-unlocked linear-response contrast is **zero when charge/density/poles/loss are identical**. Coherent site count does not multiply oscillator strength.
- Full ionic addition versus no ionic response is **1.65992 eV at zero speed**, and **0.00517950/0.00129364/0.000198969 eV** at 0.5/2/13 keV lab. This benchmark is not a measured locking contrast. Its native physical-yield factors are **1.00050243/1.00001894/1.00000024**.
- If each site's instantaneous displacement is capped at **0.0432 Angstrom**, a screened point-charge spectator triangle envelope gives **|delta U| <= 0.0639229 eV**. This is neither an RMS-derived cap nor a signed screening prediction nor an electronic-polarization bound. The optimistic positive-envelope physical-yield factors are **1.00626549/1.00092614/1.00006525** at 0.5/2/13 keV lab. Pure Gamma motion relative to a comoving target cancels exactly.

The spectator model sums `V(R-2A)-V(R)` on FCC shells, excluding the target, with `V=C exp(-r/lambda)/r`, `C=14.399645478 eV Angstrom`, `lambda=.541 Angstrom`, `a=4.084 Angstrom`. It differences the spectator potential; a constant background is not counted as a new Ue.

The ionic RPA model uses `n=4/a^3`, `kappa=1/lambda`, `qmax=(6 pi^2 n)^(1/3)`, one unit charge per D, `Omega=2 pi*10.06374152 THz`, and **gamma=2/(4.13603694 ps)**. Put `w=q^2/(q^2+kappa^2)`, `Omega_L^2=Omega^2+omega_p^2 w`. Integrate `(2C/pi) w omega_p^2 w` times the angular real average of `[Omega_L^2-(q.v)^2-i gamma q.v]^-1`. The analytic average is `Re[atanh(qv/(p+i gamma/2))/(p qv)]`, `p=sqrt(Omega_L^2-gamma^2/4)`; the static limit is `1/Omega_L^2`. Domain size only sets a lower q cutoff, not an N-fold oscillator weight. No finite domain reaches 0.5–2 keV at fixed strength/poles in this model; a changed electronic/nonlinear response is a separate open mechanism.

## Targets and fourth-cause row

At **2 keV lab**, 0.5/2 keV additional screening requires, with the response kernel fixed:

- Effective charge: **621.7/1243.4 times** a unit unlocked charge if weight is charge squared.
- Additional oscillator weight: **3.865e5/1.546e6 times** the full unit-charge ionic reference.
- Polarization: a **500/2000 V** contact-potential change; for the *same frozen spectator projection shape only*, **7.822e3/3.129e4 times** its .0639229 eV reference. This gives no bound on a different electronic response.

These are targets, not permissible physical extrapolations; new charge/poles require new response and sum-rule accounting. Exact target numbers are `response_lever_targets_at_2keV_lab`.

| Cause | Size | Density/current/isotope signatures |
|---|---|---|
| Lattice state | Ionic envelope **plus open electronic term** | Loading/temperature/pump-threshold and correlation-length onset, hysteresis or relaxation; top-10nm D x50 affects volumetric inventory, not the sheet; D2+ acts through constituent energies/counts at fixed state; half current leaves frozen physical yield/C unchanged but a beam-maintained state may relax |

`fourth_cause` also records missing DT/LiD branch predictions. The approximate **same-D-model** benchmark `delta_U_eV ~ .002588/E_lab_keV` may be used as a labelled fixed-response counterfactual only. It is not an isotope-recalculated DT response; the DD threshold term does not transfer to DT.

The collision handoff is a paired undisplaced/displaced electronic calculation: impose the corrected NiH Gamma x-polarized shape on PdD, **.0432 Angstrom per D**, metal host fixed. This is not a verified uniform locked PdD domain. Include the target in the common D displacement and reference the collision coordinate consistently. No electronic response is fabricated before that paired result exists.

## Three largest omissions

1. Displacement-induced host electronic/nonlinear response, including changed effective charge or spectral weight.
2. Quantum conditional pair preparation/contact and bound-state nuclear resonance transfer; a classical nearest-neighbour distance does not fix these.
3. Pump delivery/backreaction, long-time thermal/transport evolution and finite-record/high-barrier tails.

Parts (a)–(c) are separate matched-trajectory consumers and do not enter this standalone part-(d) archive. This preliminary archive supplies the missing plateau dependency now; it does not claim completion of the whole node.
