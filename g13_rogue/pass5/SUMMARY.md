# Pass 5: a track-shaped launch on the same foil

With the declared **2-nm coherent-bath footprint and row-local thermalization**, the track construction launches **0.000260208 W**, or **0.000377112 of incident power**. The quarter-range hypothesis gives 0.000255055 W. The unchanged Gaussian comparisons give 9.46273e-06 and 1.01908e-06 W; the historical assumed-share arm gives 0.00862516 W on the same incident denominator. No interpolation between those values was imposed. The boost is mainly thermal: 0.000252729 W thermoelastic, 6.46845e-06 W interference, 3.71897e-07 W coherent mechanical, and 6.38311e-07 W recoil raw variance (`source_terms.csv`); it does not isolate a directional-recoil gain.

| Band | Track launch W | Share of track launch | Even/Cu track core parents | N_eff | Mean elastic eV/D |
|---|---:|---:|---:|---:|---:|
| 2.5e+05–2.5e+06 Hz | 6.228e-17 | 0.000% | 1.102e+09 | 1.253e+08 | 1.992e-23 |
| 2.5e+06–2.5e+07 Hz | 6.069e-15 | 0.000% | 3.06e+08 | 7.045e+05 | 3.447e-22 |
| 2.5e+07–2.5e+08 Hz | 6.078e-12 | 0.000% | 4.187e+07 | 371 | 2.446e-19 |
| 2.5e+08–2.5e+09 Hz | 6.007e-09 | 0.002% | 5.682e+05 | 0.2012 | 1.663e-16 |
| 2.5e+09–2.5e+10 Hz | 3.459e-06 | 1.329% | 1389 | 0.001837 | 7.509e-14 |
| 2.5e+10–2.5e+11 Hz | 0.000172 | 66.105% | 9.218 | 0.001777 | 2.621e-12 |
| 2.5e+11–2.5e+12 Hz | 8.473e-05 | 32.563% | 0.02972 | 0.0003329 | 5.811e-14 |

The largest unconditioned sampled packet is 0.0179016 J/m³ at 0.063 µm, 57.084 times its own same-point mean over 8.0954e-10 s. It has 6 core electrical ancestries and N_eff=1.000217; its inventory quotient is 8.9863e-11 eV/D. Paired shuffled energy at that same place/time is 0.0179034 J/m³.

The table separates finite first-five-image/eight-sigma front ancestry, the twelve-length energetic core, and energy-weighted participation. Every source, all three beams, both loss rules and all four depths are in `readout.csv`; `same_record_peaks.csv` and `tail_comparison.csv` bind maxima and controls to their actual records.

The two-bin low-edge transfer is:

| Band / loss, even steady | 30 µm J/m³ | 50 µm J/m³ | 200 µm J/m³ |
|---|---:|---:|---:|
| 4 / Cu_extrapolated | 2.3144e-08 | 5.2262e-09 | 2.3416e-12 |
| 4 / Al_transfer | 1.0232e-10 | 2.9293e-12 | 6.5386e-23 |
| 5 / Cu_extrapolated | 1.9006e-11 | 1.8973e-14 | 1.0309e-35 |
| 5 / Al_transfer | 9.9517e-23 | 9.3983e-33 | 2.7874e-106 |

The selected 30-µm low-edge coincidence peaks at 3.26378e-10 J/m³: **1.83×** the paired shuffle’s whole-window/nine-point maximum (1.78641e-10 J/m³). It has 6 supported parents, instantaneous N_eff=2.2072, and 2.799e-18 eV/local-D. At the originally selected point/time, the comparison is 11.99×; the shuffle’s same-point window maximum is 1.28765e-10 J/m³. At the old Gaussian peak’s identical location/time, the track source gives 2.98209e-10 versus 2.10712e-09 J/m³ on the same 1.4896e-09-s record, alongside the lower matched mean above. One selected event and one shuffle do not estimate an excess or rarity.

Depth reach uses a **declared reporting threshold**, mean elastic energy density ≥10⁻²⁰ J/m³ at a sampled depth. It is not a detection floor, physical zero, or bound on occasional coherent peaks. `depth_reach.csv` and `upper_band_depth_bound.csv` apply only to the retained paths and named loss rules.

On identical physical parent histories, the β=.35 conditional mark law changes realized source work by factors 1–1.079611; its conditional expected-work factors are 1–1.048079. The unamplified physical relative-speed correction reaches only |v_z/v_D|=1.813e-09. Phase shuffle, background-off and zero-prior controls are paired. Updated launches enter actual retarded receiver packets (`feedback_receiver.csv`).

The largest conditional change is band 5/Cu_extrapolated/hot_steady: 468 local source parents in 20 ns. Work changes from 9.58777e-17 to 1.03511e-16 J. Phase shuffle retains a factor 1.079583, whereas the physical relative-speed arm changes work by only 4.6e-12 fractionally. This is chiefly the stipulated mark-energy sensitivity. Source C counts 0/468. At its matched 63-nm receiver the quiet→conditional peak is 0.000733098→0.000733098 J/m³ (3.68001e-12→3.68001e-12 eV/D); B counts 1→1 initiations, 1→1 boundary-censored, in 0.8095 ns. These source and receiver denominators remain separate.

The upper energetic-band feedback has nonempty quiet packets in 0/24 bounded retarded receiver windows. Empty windows are exposure-limited; they are not assigned a physical zero or a rarity forecast.

The source and signed-background equations are in `README.md`; consumed inputs and the bounded Pd-loss search are in `inputs_ledger.md`. Those physical inputs, rather than further range refinements, set the interpretation.

Exact archive hashes, unpublished draft IDs, and remote member readback are supplied in the delivery receipt.
