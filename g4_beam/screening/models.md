# Ue models: formulas, material inputs, reach, and knob sensitivities

Conventions. Ue in eV unless stated. e^2/(4 pi eps0) = 14.40 eV·Å. D+D Gamow energy
E_G = 986 keV. The screened cross section is sigma(E) = S(E)/sqrt(E(E+Ue)) exp(-sqrt(E_G/(E+Ue)))
(paper Eq. 7; same form in Czerski et al. 2004 EPL 68:363, Kowalska 2025, Lipoglavšek 2026; the older Bochum
analyses use the enhancement factor f = exp(pi eta Ue/E), which agrees for Ue << E).
The constant-shift form holds only when the screening length lambda is much larger than the classical
turning radius r_c = e^2/E, which is the same as Ue << E (lambda/r_c = E/Ue). In the paper's fit window
this ratio is 1.5-3.8 for PdD (Ue = 1.7 keV), and at E_cm = 0.25 keV with Ue = 3 keV it is 0.08. So the
per-energy Ue of 2-3 keV is not a screening length; it is a re-expression of the yield floor
(see section 8). All numbers marked "computed" come from `analysis.py` or `density_vs_ue.py` in this
directory; "estimate" means my own estimate with the reasoning given.

## Summary table

| Model | Ue formula | Material inputs | Pd value | Reaches 0.3-3 keV? | Main knobs (per unit) |
|---|---|---|---|---|---|
| Adiabatic / atomic | difference of electron binding energies, separated vs united atom | none (atomic) | 20-40 eV (gas definition) | No | none |
| Debye-Hückel plasma (Bonomo 2003; Raiola 2004, 2005) | Ue = 2.09e-11 sqrt(n_eff rho_a / T) eV (SI rho_a, T in K) | n_eff from Hall coefficient, atomic density, T | 318 eV (n_eff=1), 800 eV (n_eff=6.3), computed | Only with n_eff = 28 (1.7 keV) to 89 (3 keV) at 293 K, or T = 21-65 K: unphysical | dlnUe = 0.5 dln n_eff + 0.5 dln rho_a - 0.5 dln T |
| Thomas-Fermi with band DOS | Ue = e^2 k_TF, k_TF^2 = e^2 N(E_F)/(eps0 Omega) | N(E_F), atomic volume | 77 eV (N=2.3/eV/atom, estimate) | No | dlnUe = 0.5 dln N - 0.5 dln Omega |
| Free-electron TF | Ue = 42.5 r_s^-1/2 eV | valence electron count, volume | 25 eV (1 e/atom) to 37 eV (10 e/atom), computed | No | Ue ~ n^(1/6) |
| Self-consistent dielectric + cohesion (Czerski 2004, 2006; Huke 2008) | U_pol from eps_v(q) eps_c(q) with LFC and TF charge form factor; U_coh from ZBL lattice potential | valence density, core polarizability, lattice geometry | 133.8 eV = 90.6 + 43.2 (Czerski 2004 Table I) | No (factor 2-6 below HV metal data) | weak; r_s and host Z |
| Dynamic (Lifschitz-Arista) | U_dyn^2 = U_ad^2 [1/2 + (v_F^2-v^2)/(4 v_F v) ln|(v+v_F)/(v-v_F)|] | Fermi velocity | factor 0.996 (0.25 keV) to 0.88 (6.5 keV), computed | No (reduces Ue) | v_F |
| Threshold 0+ resonance + screening (Czerski 2011-2026) | sigma = flat part(E+Ue) coherently plus narrow d+d resonance in 4He | Ue (as above); resonance parameters are nuclear | fitted Ue drops to ~100 eV for clean Zr when resonance included | Apparent Ue can rise at low E; plateau needs a separate mechanism | Ue, D density |
| Thermal spike in ion tracks (Czerski arXiv:2409.02112v3) | plateau rate ~ N_0 N_0 exp(-E_A/kT_track) <sigma v>(T_track, Ue, resonance) | D density, D diffusion barrier E_A, heat capacity, conductivity, e-ph coupling, Ue | plateau 1000x Zr per D (their Table 2) | Produces a flat floor; read as per-energy Ue it gives 2-3 keV | dln(plateau)/dE_A = -1/kT_track ~ -11 per eV; dln(plateau)/dUe ~ +0.06 per eV at thermal energy |
| Kim-Zubarev / quantum tail (JJAP 45 L552 2006; Coraddu 2004) | rate averaged over a Galitskii-Yakimets momentum distribution with power-law tail | deuteron-host collision rate (host density x elastic cross section), T | not given as a Ue | Claimed large enhancement; Zubarev 2008 (arXiv:0805.4348) argues the averaging is not justified | collision broadening |
| Defects, lattice position, impurities (Kowalska 2025; Cvetinović 2023; Lipoglavšek 2010) | phenomenological: Ue ~ sqrt(m*/m) (TF), or Ue set by local DOS at the H site (Knight shift) | vacancy density, D trapping, O/C impurities, stress | Zr: ~100 -> 600 eV with O; Pd hard vs soft foil: d+d-equivalent ~0.3-2 keV vs ~0.36 keV | Yes in the data (JSI inverse kinematics, Z-scaled); no predictive formula | m*/m = 36 needed for 100 -> 600 eV (estimate) |
| Deuteron-density dynamics (Huke 2008) | not a screening model: depth profile of D changes the yield ratio | D depth profile vs time | factor 2 near-surface excess = +0.28 to +0.44 keV apparent Ue in a 2.5-6.5 keV fit, computed | Yes (as apparent Ue) | factor 4 = +0.6 to +1.0 keV |

## 0. Adiabatic limit and gas value

Ue_ad is the gain in electron binding energy between the separated reactants and the united atom
(Assenbaum, Langanke, Rolfs 1987 Z. Phys. A 327:461). For d+d it is quoted as 27 eV (Raiola, PoS NIC-IX 169,
2006), 20 eV for gas and 39 eV for the Bochum metal analysis (Vesić 2024 Particles 7:834, Table 1). Measured gas
values: 25 ± 5 eV (Greife et al. 1995 as quoted by Raiola 2006; Vesić prints 15 ± 5) and 13.4 ± 0.6 eV from the
Trojan Horse method (Vesić 2024, Table 1). No material input. This is the floor every solid-state model adds to.

## 1. Debye-Hückel plasma model (Bonomo et al. 2003 NPA 719:C37; Raiola et al. 2004 EPJA 19:283; Raiola 2006 PoS)

R_D = sqrt(eps0 k T / (e^2 n_eff rho_a)) = 69 sqrt(T/(n_eff rho_a)) m; Ue = Z_p e^2/(4 pi eps0 R_D)
= 2.09e-11 sqrt(n_eff rho_a / T) eV. Bochum extension for hydrogen content (Raiola 2006, Eq. 2.2):
Ue(T) = (1/4.78e10) sqrt(n_eff(T)(1-y(T)) rho_a / T) for y <= 1, and 0 for y > 1, with y the hydrogen solubility.

Inputs: n_eff = 1/(R_H e rho_a) from the Hall coefficient (Hurd 1972 tables; Pd remeasured by Bochum to give
n_eff = 6.3 ± 1.2), rho_a, T.

Numbers (computed, Pd, 293 K): n_eff = 1 gives 318 eV with R_D = 0.045 Å; n_eff = 6.3 gives 800 eV with
R_D = 0.018 Å. The Debye sphere then holds 1e-5 to 3e-5 electrons, so the classical continuum picture does not
apply (Huke 2008 PRC 78:015803 make the same point: Debye screening needs kT above the Fermi energy, ~1.8e5 K
for Ta). Reaching 1.7 keV needs n_eff = 28.5 at 293 K (or T = 65 K with n_eff = 6.3); 3 keV needs n_eff = 89 or
T = 21 K.

Knobs: dlnUe = 0.5 dln n_eff + 0.5 dln rho_a - 0.5 dln T. +1% linear strain: -1.5% from rho_a alone.
293 -> 600 K: x0.70. Hydrogen loading: the Bochum form predicts Ue falls with D content and vanishes at
y >= 1, which is the opposite sign to the paper's Pd -> PdD rise.

Status: Huke 2008 rank tests on the Bochum table give Ue vs n_eff rho_s = 0.49 (p = 0.013, n = 25), falling to
0.12 (p = 0.42) when the omitted In, Sn, Sb, Pb, Bi are added and 0.02 (p = 0.91) when the sign of R_H is kept;
the log-log slope is 0.41 ± 0.04 against the required 0.5 with an unacceptable chi^2 (Q ~ 2e-13). The T^-1/2 test
rests on one Pt point (Pt 675, 530, 530, 465, 480 eV at 20-340 C; Raiola 2005 J. Phys. G 31:1141). Radioactive
decay rate predictions from the same model were not confirmed (Vesić 2024 Sec. 2.5).

## 2. Static linear response: Thomas-Fermi, Lindhard/RPA, Czerski dielectric theory

Thomas-Fermi with the true DOS: k_TF^2 = e^2 N_V(E_F)/eps0, N_V = N(E_F)/Omega. Ue = e^2 k_TF/(4 pi eps0).
Free-electron TF: k_TF = sqrt(4 k_F/(pi a0)), giving Ue = 42.5 r_s^-1/2 eV (the HTML rendering of Huke 2008
shows r_s^-1/6; its prefactor 2[9/(4 pi^2)]^(1/6) e^2/a0 = 42.5 eV matches the r_s^-1/2 form).

Czerski self-consistent dielectric theory (Czerski et al. 2004 EPL 68:363; 2006 EPJA 27:83; Huke 2008 Eq. 8-9):
V(r) = (e^2/(2 pi)^3) ∫ 4 pi phi_1(q) phi_2(q) / (eps_v(q) eps_c(q) q^2) e^{iqr} d^3q -> e^2/r - U_pol,
eps_v = 1 - v(q)P(q)/(1 + v(q)G(q)P(q)) (Lindhard P, local-field correction G), eps_c for core electrons,
phi(q) = 1 - z + z q^2/(q^2 + k_TF^2) with z = 1. Cohesion term U_coh = potential-energy gain of two d vs one He in
the field of 12 neighbours (ZBL potential, 1 Å). Ue = U_pol + U_coh.

Czerski 2004 Table I (theory vs experiment, eV): Al 96.6 vs 191 ± 12; Zr 112.4 vs 295 ± 7; Pd 133.8 vs 296 ± 15;
Ta 135.7 vs 302 ± 13; C 14.3 vs -20 ± 5. Fitted free-electron polarization: U_pol,ff = Z1 Z2 (250 ± 20)/r_s^1/2
eV (Huke 2008), which is ~6x the free-electron TF value and gives 216 eV for Pd with 10 e/atom (computed).
Low-energy limit: the effective Ue at E -> 0 is 0.58 of its high-energy value (polarization part), 0.78 including
cohesion (Czerski 2004).

Static screening lengths (computed in `analysis.py`; N(E_F) inputs are estimates, see report):

| | lambda_TF (DOS) | Ue | lambda/d_MM | lambda_TF free s-p | Ue |
|---|---|---|---|---|---|
| Pd (N=2.3) | 0.19 Å | 77 eV | 0.07 | 0.57 Å | 25 eV |
| PdD0.7 (N≈0.4) | 0.48 Å | 30 eV | 0.17 | 0.53 Å | 27 eV |
| Ti (N≈0.9) | 0.33 Å | 44 eV | 0.11 | 0.59 Å | 24 eV |
| TiD2 (N≈1-2.5 per Ti) | 0.22-0.35 Å | 41-65 eV | 0.07-0.11 | 0.51 Å | 28 eV |

Ue = 0.3, 0.8, 1.7, 3 keV need lambda = 0.048, 0.018, 0.0085, 0.0048 Å. In free-electron TF, Ue = 100 eV already
needs 270 e/Å^3 (Pd valence density is 0.68 e/Å^3) because k_TF grows only as n^(1/6).

Knobs (estimates): TF-DOS dlnUe = 0.5 dlnN - 0.5 dlnOmega. With canonical d-band scaling (width ~ Omega^-5/3,
N ~ Omega^5/3): +1% linear strain gives about +1% in Ue. Pd -> PdD0.7: N falls ~2.3 -> ~0.4 and Omega rises 11%,
so Ue_TF falls ~60% (77 -> 30 eV). Pd90Ag10: N -15 to -20% (rigid band, d-holes vanish near 60% Ag), Ue -8%.
Pd90Rh10: N +5 to +15%, Ue +3 to +7%. Pd90Pt10: N about -5%, Ue -2%. Free-electron TF and the Czerski fit:
Ue ~ n^(1/6), so +1% linear strain gives -0.5%.

Status: consistent with gas and insulators, the Ti values (≤30 to 65 eV in most accelerator runs, 25 ± 10 eV
UHV) and the defect-free UHV Zr value of 105 ± 15 eV (fit with threshold resonance). Below the HV metal data by
a factor 2-6 and below keV values by 1-2 orders. It predicts PdD < Pd; the paper finds PdD > Pd.

## 3. Dynamic screening

Lifschitz-Arista form applied by Czerski (2006) and Huke (2008 Eq. 10). For a beam deuteron on a resting d in
Pd (v_F = 1.46e6 m/s for one s electron per atom, computed): U_dyn/U_ad = 0.996, 0.985, 0.960, 0.880 at
E_cm = 0.25, 1, 2.5, 6.5 keV. The correction slightly raises the apparent energy dependence of Ue but at the
10% level; it cannot produce keV values. Input: v_F (or the plasma frequency) from DFT.

## 4. Threshold resonance plus screening (Czerski group)

A narrow 0+ single-particle d+d resonance in 4He near threshold, added coherently to the flat S-factor
(S_0 = 57 keV b). Width Gamma_d = 2 P(E) hbar^2/(mu a^2) |theta|^2; strength ~ Gamma_p. Early fits used
E_R = 10 eV, Gamma = 10 eV (Targosz-Ślęczka et al. 2011 APPB 42:795); later: proton width 40 meV, phase 115°,
total width 200-300 meV, decay mainly by internal e+e- pairs (Czerski arXiv:2409.02112v3; Dubey et al. 2025
PRX 15:041004, partial pair width at least 10x the proton width). Kılıç et al. 2017 (EPJ Web Conf. 154:01024) tie
the resonance position to the local screening energy. Effect on Ue: in UHV Zr, including the resonance lowers the
fitted Ue from ~300-500 eV to 105 ± 15 eV, near the dielectric value (~110 eV). In TiD2 the screening fit gives
25 ± 10 eV and the resonance supplies a factor 10 at E_d = 1.5 keV. Lipoglavšek & Cvetinović 2026 (PRR 8:043070)
found no e+e- signal of this resonance in spontaneous fusion from D-loaded Pd. Material inputs: only Ue and the
D density. Consistent with Czerski's own UHV Zr and Ti data; untested elsewhere.

## 5. Thermal spike plateau (Czerski arXiv:2409.02112v3, 2026)

Yield floor below E_d ~2-6 keV in ZrD2, TiD2 and PdD0.05 with protons emitted from a resting centre of mass.
Model: each beam deuteron makes a track (~1 nm wide, central T ~2900-4100 K, average 1080 K, lifetime
100-1000 ps); mobile D density N = N_0 exp(-E_A/kT); rate R = N_0 N <sigma v>(T) with narrow-resonance rate
and Ue entering through exp(-sqrt(E_G/U_0)). Material inputs and their published values: D/M (Zr 2, Ti 2, Pd 0.05),
E_A (Zr 413, Ti 680, Pd 150 meV), Ue (340, 25 [370 needed in the hot track], 400 eV), heat capacity and
conductivity (similar for all three). Pd plateau per unit D density is ~1000x Zr: x16.9 from E_A and x61 from
Ue. Sensitivities (computed from their formula): dln(plateau)/dE_A = -1/kT_track ≈ -11 per eV;
dln(plateau)/dUe = 0.5 sqrt(E_G) Ue^-3/2 ≈ 0.06 per eV at Ue = 0.4 keV. Read as a per-energy Ue, any floor
R(0.25 keV) between 1e-3 and 1e-2 gives 2.3-3.2 keV (computed, section 8). Czerski v3 states that the Pd
plateau level reported by the American group agrees with theirs and that the Ti level is much larger, and
suggests ion-source energy spread could inflate fitted Ue.

## 6. Kim-Zubarev and the quantum tail

Coraddu et al. 2004 (arXiv:nucl-th/0401043) and Kim & Zubarev 2006 (JJAP 45:L552): collisions in a dense medium
broaden the energy-momentum relation, and the Galitskii-Yakimets momentum distribution acquires a power-law tail
(~p^-8). Averaging sigma over this distribution gives large low-energy enhancements in metals. Material inputs:
host atom density and deuteron-host elastic cross section (collision rate), temperature. It is not expressed as
a Ue and gives no material ranking that I found. Zubarev 2008 (arXiv:0805.4348) concludes that the results "do not
support the idea of averaging the fusion reaction cross-section over the momentum distribution". Kim's other
models (Bose-Einstein condensation of deuterons in traps) predict rates, not Ue. For our layers, the natural input
is the D momentum distribution from our own single-site oscillator (zero-point kinetic energy ~3/4 hbar omega),
which is narrow and Gaussian, so the tail would have to come from collisions, not the well.

## 7. Defect, non-equilibrium site and impurity models

- Kowalska et al. 2025 (Materials 18:1331; UHV, ZrD2): O implantation raises Ue ~300 -> ~600 eV, C lowers it
  ~600 -> ~400 eV; vacancies made by the beam (PAS, Kowalska 2023 Materials 16:6255). Mechanism stated: electrons
  localized at vacancy-D-impurity complexes, larger band curvature, larger effective mass; "according to the
  Thomas-Fermi model the screening energy is proportional to the square root of the effective electron mass".
  The authors say microscopic DFT is needed. Implied: 100 -> 600 eV needs m*/m ≈ 36 (estimate).
- Cvetinović et al. 2023 (PLB 838:137684; JSI): gas-loaded Pd, hard (cold-rolled, H/Pd 0.47) vs soft (annealed,
  0.70). Hard foil: 1H(7Li,α) 2.86 ± 0.19 keV, 1H(19F,αγ) 18.7 ± 1.5 keV, 2H(19F,p) 18.2 ± 3.3 keV; soft foil
  3.2 ± 1.9 keV for 2H(19F,p). Knight shifts of H differ between foils; H at octahedral sites vs at grain
  boundaries, dislocations and voids. Scaled linearly by Z1Z2 to d+d (estimate): ~0.95-2.1 keV hard, ~0.36 keV soft;
  with the ~Z^2 projectile scaling JSI report, ~0.2-0.3 keV.
- Lipoglavšek et al. 2010 (EPJA 44:71): 1H(7Li,α) in Pd shows large Ue only under tensile stress (1.9 ± 0.6 keV;
  Pd77Ag23 2.8 ± 0.7 keV; <0.3 keV without stress); no correlation with H content or Hall coefficient.
- Kasagi et al. 2002 and Lyakhov et al. 2008 (High Energy Chem. 42:319, cited by the paper as Lipson 2008): high Ue
  goes with low D density and high D mobility ("fluidity"); Lyakhov give a semi-empirical relation between Ue and the
  number of mobile D at the surface (formula not available to me).
- Lipoglavšek & Cvetinović 2026: thermal d+d in F-implanted Pd, Ue = 475 ± 10 eV inferred from
  sigma(0.025 eV) = (3.0 ± 1.0)e-16 b, assuming D trapped at F.

Material inputs common to these: vacancy and dislocation density, D trapping sites and energies, impurity
complexes, local DOS and spin density at the D nucleus, local stress. None gives a closed-form Ue.

## 8. Not screening, but read as Ue: D depth profile and yield floors

Thick-target yield with stopping ~ sqrt(E) (computed, `density_vs_ue.py`): a near-surface D excess that raises the
density probed at 2.5 keV relative to 6.5 keV by a factor F mimics a Ue shift in a 2.5-6.5 keV fit of
+0.16/0.25 keV (F = 1.5), +0.28/0.44 keV (F = 2), +0.61/0.96 keV (F = 4), for base Ue = 0.3/1.1 keV. The paper's
Pd -> PdD shift (+0.9 keV; +0.8 keV in the thick-target refit) is therefore equivalent to a factor ~2.2-4 change
in depth contrast. Brief fact 6 gives PdD/Pd absolute ratios of ~1.05 at 3 keV and ~0.45 at 6.5 keV, a contrast of
~2.3. A flat floor R at E_cm = 0.25 keV maps to per-energy Ue = 1.7 keV (R = 1e-4), 2.3 keV (1e-3),
3.2 keV (1e-2), 4.4 keV (5e-2). Huke 2008 argue the Bochum and Tohoku high values are density dynamics of this kind.
Chen et al. 2025 (Nature 644:640) report +15(2)% rate from electrochemical loading at E_cm ≈ 15 keV without
extracting Ue; read as screening at 15 keV that is ΔUe ≈ 0.5 keV (estimate), read as density it is +15% D.

## 9. DFT-based estimates of d-d screening in a metal

Not found as a direct DFT calculation of Ue in the sources I read. Closest: Huke, Czerski et al. 2008 EPJA 35:243
(ab initio Hartree-Fock of D in Li and Ta lattices; electrons migrate to D; Ue stays clearly below experiment);
Salzmann & Hass 2008 EPJA 36:359 (3D Thomas-Fermi of D in metal vs insulator; explains part of the metallic increase);
Czerski 2004 dielectric theory (above). The Czerski 2004 screening function shows Friedel oscillations and an
attractive tail, consistent with DFT dihydrogen calculations they cite. Kowalska 2025 explicitly call for DFT of
vacancy-D-impurity complexes. A DFT route for us is given in report.md.

## Which models fit which data

| Data | Fits | Does not fit |
|---|---|---|
| Gas 25 eV, insulators ≤30-60 eV, C ≈0 | adiabatic, TF, dielectric | Debye with n_eff > 0 |
| Ti/Zr/lanthanides at RT, high D/M: ≤30-125 eV | TF, dielectric (+resonance for Ti) | Debye unless n_eff(1-y) → 0 is imposed |
| HV accelerator metals 190-800 eV (Bochum, Berlin, Tohoku) | Debye by construction (n_eff chosen; Pd remeasured); density dynamics (Huke) | TF/dielectric (factor 2-6 low) |
| UHV clean Zr 105 eV; with defects/O 340-600 eV | dielectric + resonance (clean); defect/m* (damaged) | Debye (no T or n_eff change) |
| Glow discharges Ti 610, Pd 1000 eV (E_cm < 2 keV) | none quantitatively; ion energy spread and density dynamics are alternatives | TF, dielectric, Debye |
| JSI Pd hard foil, d+d-equivalent ~0.3-2 keV | lattice-position / defect picture (qualitative) | all static electron-density models |
| Paper: Pd 0.8, PdD 1.7, Ti 0.3, TiD 0.9 keV (2.5-6.5 keV fit) | density-profile mimic (factor 2-4) plus defect picture; Bochum Pd value 0.8 keV agrees | TF, dielectric (wrong sign for loading), Debye (wrong sign) |
| Paper plateau, per-energy Ue 2.2-3 keV at 0.25 keV | thermal-spike floor (Czerski) or any floor (section 8) | every constant-Ue model |
