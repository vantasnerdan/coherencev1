# transport: near-surface deuterium with beam damage

`transport.py` is a 1D finite-volume model over the full 250 um foil. Diffusion follows the gradient of the D chemical potential on Pd-D (alpha/plateau/beta) and Ti-D (alpha/plateau/delta) isotherms, so hydride fronts move. It carries vacancy traps (Pd: 0.23 eV, 6 D per vacancy, 1% saturation), vacancy diffusion and supply, surface desorption against 0.037 Pa chamber D2, sputter erosion, a chemisorbed D monolayer, an optional overlayer (sputtered, regrown, with its own permeance and D capacity), a back-side fugacity limited by electrochemical supply, and beam heating. Knobs are arguments to `Model(case)`: `overlayer`, `d_ol`, `j_rel`, `dT`, `strain`, `alloy`, `s_front`, `vsat`, `no_gas`, `eta_abs`, `v_init`.

`bca.py` is a binary-collision Monte Carlo (ZBL potential, Lindhard-Scharff stopping scaled to ICRU 49, NRT vacancy counting). Its tables rebuild on first run. It agrees with RustBCA within 4% in range (scripts in `lit/bca_scripts/`). `run_all.py` and `post.py` produce `summary.csv`, `ratios.csv`, `knobs.csv` and `ue_bias.csv`; `reach.py` and `reach_ti.py` give how far back-side loading gets (`reach_pd.csv`, `reach_ti.csv`); `checks.py` gives grid, mass-balance and ion-count checks (`checks.txt`).

Results:
- Clean Pd holds D/Pd ~0.01 in the top 70 nm, almost all in beam vacancies. Dissolved D is pinned at 3-4e-5 by the chamber gas.
- Ti and TiD both reach TiD ~1.85 under the beam, so the profiles give TiD/Ti = 1.00.
- Back-side loading reaches the beam side of Pd within minutes but stays alpha unless surface sticking is below ~1e-6. In Ti it never arrives within the experiment.
- 1-5 nm overlayers erode in 0.2-4 h at this flux.

Drops: 3D microstructure, blisters, ballistic mixing, vacancy clustering, run history, electrochemical kinetics, surface chemistry beyond one sticking coefficient, angular scattering in the yield integral.

Three biggest omissions: beam-side surface chemistry (carbon, sulfur, real adsorbed D), microstructure of the 200 nm sputtered film, run history and lateral beam profile.

Estimated inputs, to be replaced by computed ones: vacancy saturation (1%), H-vacancy binding (0.23 eV), Pd-D plateau pressure (5.2 kPa), oxide permeances and D capacities, alloy isotherms, Ti vacancy migration energy (0.5 eV), electrochemical absorption fraction (5%).
