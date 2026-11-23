# Beam — a genuinely energetic tail is possible only conditionally

**Mechanism (two sentences).** If the ordinary40eV barrier remains, an above-nominal incident tail could supply keV reacting pairs at the0.5keV-lab setting, whereas molecular ions divide the setting's energy among constituents. Neutralization preserves pre-conversion energy, but the documented low-setting optics do not establish a distinct keV beam subsequently decelerated to0.5keV.

**Species and monitor.** For Dn+ electrical-current fractions fn, yield/C=`Σ n fn Yn(L/n)/e` and atoms/C=`Σ n fn/e`:85/15 means1.15 D/charge, not85/15 atom fractions. The [Methods](https://doi.org/10.1038/s41467-026-74421-1) report approximate85/15 D+/D2+, not a measured spectrum at every setting; baseline D3+=0 is assumed. At **[.25,.5,1,2.5,6.5]CM keV**, pure D+ gives physical proton factors≈[1.1765,1.1765,1.1765,1.1762–1.1763,1.172–1.174]; raising D2+ current fraction to85% gives [.17647,.17647,.17648–.17649,.1771–.1777,.1885–.1968]. Pure D2+/D3+ suppress further. Profiles, sheets and stopping stay fixed; only incident paths change.

**Predicted versus required tail.** Adding **3e−4 neutral D per true charged unit** at13keV lab,5%FWHM, predicts these **physical proton Y/C factors**, using the same monitor response:

| Foil | Fp at the five CM energies |
|---|---|
| Pd | [1.33e15,6.80e8,1.01e4,1.331,1.000352] |
| PdD | [1.34e15,7.40e8,1.15e4,1.373,1.000352] |
| Ti | [6.59e16,8.57e9,4.59e4,1.628,1.000352] |
| TiD | [6.59e16,8.57e9,4.60e4,1.630,1.000352] |

Own-reference R divides by the changed neutron6.5 factor and is1 there; physical factors above are not forced to1. The.25 observations require **minimum local neutral D/charge** [0.01175,.05109,.000944,.002230], alternative-current [.001875,.008157,.000151,.000356], optimistically leaving the reference tail absent. All five-energy lower/paired requirements are in `beam_neutral_requirements.csv`; the physical baseline and required normalized excess remain separate.

Thus13keV can supply scale **if present**, but common13keV one-parameter fits give χ²96–122/16 nonreference points; the tested5/13keV constant-ratio/constant-rate families give χ²61–122, not an arbitrary-spectrum exclusion. At2keV lab the minimum.25 requirements are1.29e5–8.76e5 D/charge; at5keV,1.70–33.44. Ratios>1 are allowed for extra unmetered particles, but need a source/power budget beyond dilute chamber conversion. For δ=1, solved f is **at I_INT**, with local `f(E)=f·I_INT/Itrue(E)`, not a local fraction.

**What the source actually supplies.** Reported0.5–5kV gun delivery uses three grids/4cm flight; increasingly negative target bias begins only above5kV. Charge exchange or scattering of the low-setting beam cannot create>nominal energy; pre-bias neutrals at high settings are under-energetic. [Czerskiv3](https://arxiv.org/pdf/2409.02112v3) documents a different pre-deceleration neutral route (~3e−4), not a Davis measurement. The charged-only5% width was inferred over2.5–4.2kV, not from low-setting neutral spectroscopy.

All2.5–3e−4Torr chamber gas treated as D2 at300K over4cm gives column3.22–3.86e17m⁻². A **declared** σ=1e−20–1e−19m² bracket gives.32–3.79% conversion, not a measured fast fraction or whole-source bound. Torr is used because the printed4e−3mbar parenthesis is inconsistent. [Janev–Reiter–Samm](https://www.eirene.de/Documentation/report_4105.pdf), §5.2.1/Table31, gives a same-velocity **H+/H2 proxy**, not isotope-resolved D+/D2 data. Optical depth does not establish converting energy: no actual grid-voltage/pressure/neutral-acceptance map supplies the required13keV pathway at0.5keV.

**Coincident impacts/shared wake.** Conditional1Å separation at impact gives D2+ arrival gaps≤[.646,.457,.323,.204,.127]fs and transverse offsets≤1Å; orientation/breakup are unknown. Tiny independent-ion overlap does **not** bound this within-molecule correlation. Nevertheless, D2+ supplies only2.1e−9–3.9e−9 of baseline.25 proton yield; enhancing only the second constituent needs σ multipliers2.7e25–4.6e26 even with an unchanged reference. Granting it a rigid+.22eV throughout its path would give only total Fp(.25)=1+(5.6–10)e−11. That is an optimistic scale sensitivity: g6's actual one-ion2keV head-on result at0.1Å is not time-step checked, not a positive rigid shift and **not a molecular-wake theorem**.

**Absolute floor and non-yield discriminator.** Model-required.25 reactions/true C are [4898,21604,52355,123659], alternative [782,3449,8359,19744], conditional on the6.5 neutron anchor; all five energies are in `readout_floats.csv`. SI1 instead gives detected counts/cup C [1.075,1.535,.768,2.865], requiring `f_col/ηp` for emitted reactions/true C; no efficiency is invented. Neither these central values nor the four-energy table establish one common constant physical floor. A **foil-plane neutral/species energy spectrum versus voltage**, joined to actual grid potentials, would identify the missing incident-energy component. Fixed-composition molecular rate scales∝I; constant upstream neutral rate gives Y/C∝1/I, unlike a target state injected and probed by independent ions.

**Exactly three largest omissions:** (1) grid potentials and foil-plane neutral/species spectra/acceptance; (2) molecular breakup and the joint two-projectile electronic response; (3) same-run detector/charge calibration connecting relative curves to emitted rates.

**Finite proof:** atom-current identity1.11e−15, neutral changed-consumer agreement1.11e−16 relative; conditional timing calculated directly, exit0 (`beam_proof.json`). **Composition alone fails; an energetic tail is quantitatively possible but not established by the documented low-setting energy path.**

