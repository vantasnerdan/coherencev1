# Live foil, structured beam, peaks and rogue holes

## Built result and the first comparison

This is the corrected 4097 / 4103-as-corrected-by-4104 / 4109 construction. **No persistent damage or dose population is included.** The 48 live cases propagate a continuous, state-dependent *distribution* of signed launch marks through the cubic-grain vector foil, rather than putting one fixed small packet into a quiet medium. The conservative material-nonlinearity witness remains separately labeled in [fixed_source_readout.md](fixed_source_readout.md); it was not combined with this live-source sweep.

**Read the first A/C values with their actual resolution:** the local field is averaged over a 2.5-mm-square, 0–50-μm source neighborhood and frozen immediately **before each 50-ns bin**. Parents within that bin do not qualify themselves or one another. The retained **eigenmode band** is 5–250 kHz (53 modes, 28 classified flexural for attenuation), not a hard temporal Fourier filter of impulsive accelerations/flux. Neither a 54-ns bulk thickness front nor atomic splash structure is resolved by this low-band model. Two independent 100-μs records follow five maximum modal energy lifetimes of burn-in. The physical parent population is not compressed in time.

Primary structured example: **PdD-profile source, 2.5 keV, intrinsic energy lifetime 50 μs, L_E=.2 m, live feedback beta=.35, peaked/hot-spot spatial beam plus correlated bursts/ripple**. A is **93.429% of impacts**, C_any **1.0642%**, and resolved off-impact directional B occupancy **40.275%**. Uniform/steady comparison at identical recorded charge gives **94.932%, 1.4243%, 30.350%**. The construction therefore does not assume either that these classes are absent or that all forms of clumping change them monotonically.

- A: an impact into a pre-existing represented-band disturbance above .25 of the fixed quiet-reference rms.
- C: impact into either sign of a pre-existing epsilon_xx or u_z excursion beyond **the same ±3 causal-prior rms**. Positive peaks and negative holes are separate columns.
- B: inward Poynting flux from distinct parent-origin sectors and ≥90-degree inward-face directions at an entirely off-spot receiver, plus constructive quadratic interference. Multiple echoes of one sector cannot satisfy the distinct-origin requirement.

A/C are impact fractions/rates. B is receiving-window occupancy plus a connected-event rate per **receiver second**. They are not summed into a supposed whole-foil event rate. The source-cell and receiver banks are the explicitly represented sampling set; no overlapping-region tally is extrapolated over the disk.

## Spatial × temporal clumping at matched charge

Spatial levels: 0 uniform, 1 normalized Gaussian peak, 2 peak plus two local hot spots. Peak cell current density / uniform is **1 / 1.983 / 4.685**. Temporal levels: 0 steady, 1 .7-amplitude 100-kHz ripple, 2 ripple plus a random Markov burst process (mean ON/OFF 2/18 μs, added gain 9). These are hypotheses, not measured beam maps. Temporal rate CV is **0 / .495 / 1.361** for the canonical records; the burst arm's peak/mean is **6.117** and its 1-μs rate-modulation correlation **.526**. The final arm is genuinely random-intensity correlated timing, not deterministic ripple relabeled as bunching.

Each row below shares exactly **172,265,586,279 parent impacts** across the two measurement records. That is a common Poisson-drawn charge total allocated under different normalized spatial/time intensity paths, not a higher-current comparison. It is not 172 billion independent wave realizations.

| Space | Time | A impacts % | C_any impacts % | B receiver occupancy % | B starts / receiver s |
|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 94.932 | 1.4243 | 30.350 | 160,000 |
| 0 | 1 | 96.192 | .7050 | 16.475 | 60,000 |
| 0 | 2 | 97.720 | 1.1249 | 41.875 | 80,000 |
| 1 | 0 | 97.146 | 2.5677 | **56.000** | 145,000 |
| 1 | 1 | 97.969 | 1.6424 | 29.475 | 70,000 |
| 1 | 2 | 91.136 | 0 sampled | 26.175 | 95,000 |
| 2 | 0 | 96.328 | .3686 | 26.425 | 75,000 |
| 2 | 1 | 95.825 | 1.2072 | 38.450 | 100,000 |
| 2 | 2 | 93.429 | 1.0642 | 40.275 | 135,000 |

At space/time=2/2 the separate C fractions are strain-positive **.1146%**, strain-negative **.1446%**, displacement-positive **0 sampled**, displacement-negative **.8936%**. Their union is 1.0642%, not their arithmetic sum. A∩C equals 1.0642% here and A∪C equals 93.429%. Global temporal coincidence with remote B windows is separately labeled, not mistaken for co-location.

The impact-class launch heights are actual marked **self-energy per parent**, not energy per D or a collective regional peak: A/C_any mean **18.88 / 28.13 eV** in the 2/2 example, with sampled support maximum **88.55 eV per parent**. Uniform/steady gives **18.13 / 28.94 eV**. The quiet-law mean would be about 6.31 eV per charge parent at this source. All low/middle/high category fractions are retained in `launch_distribution.csv`; the propagator consumes their realized signed counts, not those means. Coherent work and the existing-wave energy loan are accounted once, with new work capped by the incoming host budget.

B's largest sampled occupancy is 56%, with center quadratic density peak **.01724 J/m³** and peak/mean about **4.26** in the peaked/steady case. B height uses a **receiver-center one-point box cubature**, not a fully volume-integrated regional mean; its matched point thermal mean/fluctuation columns are explicitly separate. Fully integrated region observables are in the signed/energy tables.

## Lifetime, range and standalone holes

Intrinsic lifetime is scanned at **2,50,200 μs**, independently of energy attenuation length **2 mm,2 cm,20 cm,2 m**. The table reports *first tested* 50%-common points, not an interpolated phase boundary. A first crosses 50% at **50 μs** for the tested ≥2-cm ranges in both uniform and fully structured examples; at 2 mm it remains below 50% (maximum about 41% in those grids). B crosses 50% at the peaked/steady 50-μs, .2-m point. C_any's largest observed impact fraction is **3.474%**, at uniform/steady, 200 μs, 2 m, beta=.35; no tested C row reaches 50%. This is a result of the built source/field law, not an opening statistical exclusion or a universal rarity claim.

Holes are also scored **without requiring an impact**, with the same signed ± thresholds, region sizes and duration conventions as positive peaks. `signed_tails.csv` includes both epsilon_xx and u_z, sides **2.5,10,21.269 mm**, depths **0–50,25–75,175–225 μm**, thresholds **1,3,5,10 causal-prior rms**, complete-event rates, censoring and duration tails. Near-zero supported displacement-average projectors are identified rather than divided by roundoff.

One actual pair of signed strain examples in the uniform-spatial/random-burst arm (50 μs,.2 m): extension **+2.059e−6**, standardized height **4.172**, at 25–75 μm; compression **−1.960e−6**, standardized depth **4.090**, at 175–225 μm. At the three-rms level each has one complete event over .2 ms, durations **13.1 / 13.7 μs**. Different depth locations are stated; this is not evidence of intrinsic positive/negative material asymmetry. All signs and locations, including zero sampled counts, remain in the tables. Large physical parent counts do not replace the two independent acoustic records.

## Time frames: arrivals continue while waves are in flight

The primary disk geometry is 250 μm thickness, **24 mm diameter**, beam area 1 cm² (disk-equivalent width **11.2838 mm**). The actual flange/clamp radius is not measured. The equal-area square model side is 21.269 mm; its own edge times are tabulated separately. Hill references use rho=12007 kg/m³ and the grounded near-300-K pure-Pd tuple. Thin-plate D_bend=**.201317 N m**, alpha=**.258972 m²/s**. Flexural speed is **group** speed 2 alpha k, not phase speed alpha k.

| Wave reference | Group speed m/s | Thickness ns | Spot crossing μs | Center→edge μs | Edge and return μs | Whole-spot parents during return, 38 / 138 μA |
|---|---:|---:|---:|---:|---:|---:|
| Bulk Hill L | 4620 | 54.11 | 2.442 | 2.597 | 5.195 | 1.232e9 / 4.474e9 |
| Bulk Hill T | 1988 | 125.79 | 5.677 | 6.038 | 12.075 | 2.864e9 / 1.040e10 |
| Flexural λ=5 mm | 651 | guided | 17.337 | 18.437 | 36.874 | 8.746e9 / 3.176e10 |
| Flexural λ=10 mm | 325 | guided | 34.673 | 36.874 | 73.748 | 1.749e10 / 6.352e10 |
| Flexural λ=15 mm | 217 | guided | 52.010 | 55.311 | 110.622 | 2.624e10 / 9.528e10 |

A guided flexural mode has **no separate thickness flight**. The optional λ=20-mm reference is below the selected 5-kHz cutoff and is labeled outside the band. The whole-spot totals above are raw charge-parent arrivals, not simultaneously coherent packets or source-age durations.

Within **one wavelength's radius** of the spot center, the neighborhood is clipped/normalized to the spot. Bulk wavelengths at the 100-kHz reference and flexural λ≥10 mm include the whole spot. For λ=5 mm, the uniform clipped area is .7854 cm²: return counts are **6.869e9 / 2.494e10** at 38/138 μA; the largest structured source-neighborhood count is about **7.537e9 / 2.737e10**. These structured values use the explicit coarse source-cell profile, not a measured microscope map. The full table contains both currents and every thickness/spot/span/edge/return path.

## Surviving lateral returns are actually modeled

Energy survival is `exp(−path/L_E−time/tau_E)` times the assumed edge **energy** reflection coefficient; amplitude survival is its square root. At intrinsic lifetime 200 μs and lossless edge, actual-disk-proxy round-trip energy fractions are:

| L_E | Bulk L | Flexural λ=10 mm |
|---:|---:|---:|
| 250 μm, earlier attenuation context | 1.98e−42 | 1.40e−42 |
| 2 mm | 5.99e−6 | 4.25e−6 |
| 2 cm | .2935 | .2083 |
| 20 cm | .8642 | .6134 |
| 2 m | .9627 | .6834 |

Thus the cm/m cases carry actual surviving lateral returns in the free variational field. The edge table labels ≥10% energy at edge and return; it does not just assume every mode survives. A separate R_edge=.1 modal-leakage arm preserves the reflected phase construction while losing 90% energy per nominal round trip. At 200 μs,2 m, structured beam, B occupancy changes **40.875→30.225%**, B peak center density **.05270→.004744 J/m³**, and connected starts **60,000→135,000/s**. Occupancy, height and fragmentation rate are different observables. Flange reflectivity, dispersion mapping and the square boundary are model parameters, not Pd/mount measurements. The D-mode 4.1-ps lifetime is **not** used for plate ringing.

## Backwards requirement readout and next choice

The retained **0.6423 meV/D** is full-loading context, not a material threshold at actual x<1. `output/live/region_peaks.csv` and its finite-event frontiers are rebuilt from the live fields. Their conversion uses unit transfer, 4.1-ps holding and a **median-modal-lifetime local-loss proxy**; it is not the true spatial dissipator or a new strain-to-D calibration.

The largest sampled conditional conversion proxy in this sweep is **4.662e−16 eV/D**, leaving factor **1.378e12**, at 50 μs,2 cm, structured beam, in the full-area 0–50-μm slice. Sustaining that peak indefinitely would require holding **5.65 s**. A different best finite-event envelope, at 50 μs,.2 m, hot-spot/steady beam, has an observed 7.65-μs event. Even putting its receiver peak throughout that duration and taking infinite holding gives **2.777e−10 eV/D**, residual concentration **2.313e6**. No holding-only solution exists for that **observed event envelope**; this says nothing universal about continuing-background storage.

`finite_event_joint_frontier.csv` reports the actual finite-duration concentration/holding tradeoff, rather than pretending that holding beyond an event grows its energy without bound. The separate fixed-source direct-family attenuation/nonlinearity/frontier calculation remains in the baseline readout and is not a result for the new classes. No fitted cheapest combination or optical Kerr coefficient is transplanted into Pd.

The decision-changing next input is the **joint space/time beam structure and the field-dependent launch response**, followed by frequency-resolved loss/mount response; the source-law and profile controls already move class populations differently. Connecting the acoustic reservoir to retained D storage remains the next physical consumer, not something more wave samples or an overlap tally supplies. No new MD, electronic, force, pair or damage campaign is started here.

## Evidence and scope

The live CLI completed naturally: **48 cases,576 impact-class rows,6912 signed-tail rows**, plus matched quadratic-energy tails and finite-event maps. Equal-charge factorial check: one unique total. Source/propagation energy closure is **1.64e−15 of incoming host energy**, maximum numerical source-cap excess **3.31e−24 J**, harmonic-reference covariance residual **1.06e−15**, and A/C union identity error **8.67e−18**. These are finite-model behavioral/self-consumer checks, not an independent continuum/material certificate.

The three classes use the actual pre-bin field. Standalone signed holes, thermal comparisons, birth rates, occupancy, censoring, source launch spread, same-parent molecular marks, and lateral survival each retain their units/denominators. Unknown beam, launch, loaded tensor, grain texture, mount and transfer inputs remain explicit. The delivery archive is unpublished and scope-only; its external byte/hash receipt is `delivery.json` rather than a self-referential hash inside the archive.

## Physical pre-impact class heights and matched consumer

For the primary structured 2/2 example, the additional impact-state aggregation uses **the identical parent counts and source trajectory**: all 48 cases match their saved source work with difference **0 J**. Existing signed/energy tails and convergence labels were reused, not reinterpreted.

| Tagged impacts | Mean absolute prior epsilon_xx | Prior epsilon_xx range | Mean absolute prior u_z | Prior u_z range |
|---|---:|---:|---:|---:|
| A | .9391e−6 | −3.503e−6 to +3.336e−6 | 91.84 nm | −355.8 to +315.4 nm |
| C_any | 2.1500e−6 | −3.083e−6 to +3.040e−6 | 260.6 nm | −355.8 to +205.5 nm |

The positive-strain C subgroup has mean absolute strain **2.350e−6**, negative-strain subgroup **2.422e−6**, and the negative-displacement subgroup mean absolute displacement **279.5 nm**. These are the **pre-existing** heights that the new impact encountered, separate from its marked launch self-energy. No positive-displacement C impact was sampled in this example, so its conditional mean/height fields are **empty, not zero-valued measurements**.

C_any fractions across the two independent records are **0 and 2.1284%**, averaging to the reported parent-weighted 1.0642%. A fractions are 90.05–96.81%. The large number of physical impacts does not hide the small number of independent wave histories or manufacture a material asymmetry.

At the peaked/steady B example, matched **point-observer** harmonic thermal mean is **4.515e−13 J/m³**, fluctuation sd **3.961e−13 J/m³**. The B peak is **3.82e10** thermal means or **4.35e10** thermal standard deviations in the same modal space. This is not an all-band lattice temperature, a regional integrated thermal energy, or retained-D optical excitation.

**Band/exterior boundary:** below-5-kHz global plate fundamentals and the true mounted-disk coherent net-traction response are not part of the live modal forcing. The separately retained whole-beam impulse/shot-noise plate budget is conditional on support/Q; the table does not claim a measured clamp. A downward band extension under a specified mount and launch spectrum could materially change displacement/hole memory. Burst correlation time relative to the flexural return/period is likewise a deciding source refinement; neither is replaced by a finer short-wave solver or more overlap counts in this draft.
