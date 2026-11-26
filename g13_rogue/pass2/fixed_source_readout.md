# Whole-foil coarse result

## What was actually propagated

The finite vector Ritz model represents the mass/area of a 24-mm Pd disk by a 21.269-mm free square, thickness 250 μm, with a **central** 1-cm² beam. It is not the measured mounted circular boundary. It retains 53–56 modes in the requested 5–250 kHz band, 324 vector polynomial coordinates before band selection, and synthetic 2.659/5.317-mm columnar grains. Cubic anisotropy is in the rotated tensor of the actual variational operator, not a phase-velocity emission histogram or grain gain. Short-wave caustics and the real texture remain unrepresented. The shear speed scale implies bulk wavelengths of roughly 6–10 mm at the upper band edge; grain/wavelength ratio is order .3–.9 there, not a geometric-ray limit. Collective flexure has its own dispersion.

Six independent **300-μs physical-time records** follow 250 μs burn-in; dt=50 ns. Each receiver/arm/case has 1.8 ms sampled exposure, not six times the foil area. Parents arrive at the chart I/e rate; exact signed Poisson count differences aggregate the very large arrival population without accelerating time. No sampling-window conditioning is used. Baseline band launch fraction is .01 and energy lifetime 50 μs, both unmeasured sensitivities. The latter is a 0.2-m loss-length proxy at 4 km/s, not a measured Pd loss. The nominal Pd and loaded-profile PdD runs use the same declared pure-Pd elastic/mass operator and the same count realizations, with their own accepted stopping and D populations.

## Linear tails and size/depth dependence

Representative **PdD, 2.5 keV**, 2.5-mm-square receiver, depth 25–75 μm. Quadratic strain energy is averaged over the receiver, not strain averaged before squaring. Counts are complete excursions above four times each arm's own beam mean.

| Arm | Mean / peak J m⁻³ | Peak / mean | Time fraction ≥4 mean | Complete excursions | Longest / median duration μs |
|---|---:|---:|---:|---:|---:|
| Hill isotropic | .05457 / .29788 | 5.458 | .02481 | 4 | 18.05 / 11.50 |
| Single cubic | .05338 / .27043 | 5.066 | .00939 | 3 | 14.40 / 1.60 |
| Cubic grains 2.659 mm | .04581 / .28257 | 6.169 | .02500 | 4 | 18.80 / 10.225 |
| Cubic grains 5.317 mm | .05861 / .27073 | 4.619 | .01289 | 4 | 8.20 / 5.525 |

At this receiver no arm sampled an eight-mean event. Under a Poisson excursion-rate approximation only, zero counts over 1.8 ms give a 95% rate-model upper bound of **1664 s⁻¹**. This is not a distribution-free bound for correlated crossings. At four mean the small-grain case has 4/2/0 complete events lasting ≥1/10/50 μs. Record-to-record variation is large: its four-mean time fraction spans 0–.0633 across the six histories (descriptive standard error .0126). These observations do not establish a tail-law ordering, much less a universal Rayleigh or exponential intensity distribution. Finite real multimode quadratic strain is not the complex diffuse single-mode limit.

For that same small-grain case, increasing the receiver side **2.5 → 5 → 10 → 21.269 mm**, in the 25–75-μm slice, gives peak/mean **6.169 → 4.874 → 4.075 → 2.215** and four-mean time fractions **.0250 → .00747 → .00114 → 0**. Complete counts are 4/3/1/0. Depth slices 5–55, 25–75, 175–225 μm give peak/mean 6.211/6.169/6.070 and longest four-mean events 18.75/18.80/17.20 μs. The low-band coarse field is not resolving the nm beam layer through these depth windows.

A 16-receiver union is computed directly as `max_R(u_R/mean_R)` at each time, not by adding overlapping spatial event counts. Its own mean is 1.4169 for this case; at four times **that bank mean** it gives one complete crossing over 1.8 ms. This bank statistic has no J/m³ unit or artificial thermal comparator and is not a survey of every possible hotspot center on the disk.

## Matched-band thermal and nominal material comparison

The harmonic thermal quadratic covariance uses the **same retained modes, strain components, tensor and region integration**. At the representative small-grain receiver, thermal mean is **3.9023e−13 J/m³**, fluctuation sd **3.0639e−13 J/m³**, and the correlated quadratic-rank proxy is **3.244**, not an assumed independent cell count. The peak is 7.24e11 thermal means or 9.22e11 thermal standard deviations in this exceptionally low band. Neither comparison is an optical-D energy or all-band lattice-temperature claim. The thermal mean does not vanish as the region grows: in the 25–75-μm slice its mean varies .390/.396/.489/1.035 ×10⁻¹² J/m³ across the four sides because modal participation is nonuniform, whereas its sd varies .306/.263/.212/.223 ×10⁻¹². Do not impose independent-volume scaling on these correlated low modes. The whole-volume receiver has rank 53 and the expected equipartition strain-energy trace.

Matched small-grain, 2.5-mm, 25–75-μm quadratic peaks:

| Nominal target / CM keV | Peak J/m³ | Region n_D m⁻³ | Acoustic peak eV/local D | Conditional conversion eV/D |
|---|---:|---:|---:|---:|
| Pd / .25 | .0099896 | 2.1878e24 | 2.8500e−8 | 2.3370e−15 |
| PdD / .25 | .0099890 | 2.1802e28 | 2.8596e−12 | 2.3449e−19 |
| Pd / 2.5 | .282888 | 2.1009e24 | 8.4043e−7 | 6.8915e−14 |
| PdD / 2.5 | .282572 | 1.9958e28 | 8.8368e−11 | 7.2462e−18 |

The conversion column is separately conditional on **unit** local transfer, quadratic acoustic dissipation `u/tau_energy`, and D energy lifetime 4.1 ps. It is not a calibration or a prediction of optical population. Different per-D values mostly reflect different denominators, not a different measured elastic operator. Saved 12-h/2-h loading profiles are not the 300-μs acoustic exposures. The .25/2.5-keV supported endpoints are not extrapolated to 6.5 keV, where this stopping source is absent.

## Backwards requirements, without a material verdict

The retained **0.6423 meV/D** is full-loading context, not a known threshold at actual x<1.

1. **Concentration alone, direct attenuation family.** At PdD 2.5 keV, depth 50 μm, L=100 μm, unit launch and unit transfer give 6.0549e−15 eV/D. Required factor is **1.0608e11**. The analogous loaded .25-keV value is 2.1949e−16 eV/D, factor **2.9263e12**. Those direct-family fractions differ deliberately from the modal .01 band launch.
2. **Attenuation length alone.** Conversion contains K(L,z)/L, not K alone. It tends to zero at both L→0 and L→∞. The finite scan maximum for 2.5 keV, 50 μm is at **L≈115.6 μm**, yielding **6.0896e−15 eV/D** and residual **1.0547e11**. Optima at depth 30/200 μm are about 69.3/460.4 μm. This is a scan optimum, not a certified global interval. No attenuation-only length is claimed to reach the reference.
3. **Holding alone, steady versus event.** Sustaining the direct-family L=100-μm conversion would require D energy lifetime **.4349 s**. The modal representative peak, at .01 band launch, would instead require **363.4 s** if sustained indefinitely. It is not sustained: the observed four-mean excursion lasts at most **18.8 μs**. Even setting its global receiver peak constant throughout that event and sending all dissipation to D, infinite holding stores at most **3.3226e−11 eV/D**, leaving factor **1.9331e7**. Thus no finite holding-only solution exists for **that observed event envelope**. Continuing-background accumulation is a separate mechanism, not ruled out by a finite excursion.
4. **Positive nonlinearity alone, fixed energy.** For the declared positive quartic model and total budget H0=P_band*tau=4.3561e−8 J, any nonnegative strength obeys `u_R≤lambda_max(A_R) H0`. At the representative receiver the bound is **3.3185 J/m³**, or conditional 4.1-ps converted **8.5100e−17 eV/D**, still factor **7.5476e12** from the reference. H0 is a reference budget, not a global bound on the unbounded Poisson ensemble. A larger coefficient cannot remove this fixed-energy inequality; it does not prove a bound for other energy/source models.

`joint_frontier.csv` gives the direct-family steady Pareto curve after its best L and sampled nonlinear envelope. With a symmetric squared-log concentration/holding cost proxy, equal log sharing is about **3.21e5 concentration × 3.21e5 holding**, holding near **1.32 μs**; this is an explicit arbitrary cost choice, not the universally cheapest physical combination. `finite_event_joint_frontier.csv` separately uses `C*tau_hold*(1-exp(-T/tau_hold))=reference/r`, saturating at the finite-event concentration asymptote. Neither curve silently combines the direct hemisphere's normalization with the foil's .01-band source.

## Whole-beam net traction

At 2.5 keV the absorbed mean momentum flux gives **2.1177 μN**; perfect reversal bounds this simple incident impulse by **4.2354 μN**. One-sided parent shot-noise force PSD is **1.0615e−26 N²/Hz**, retaining molecular-parent correlation. For illustration only, a unit-overlap 10-kHz supported mode with full-foil effective mass 1.358 g and Q=100 has static response 3.95e−13 m and shot-noise mode energy 3.11e−27 J. Support, mode overlap, effective mass and Q are not measured, so this modal value is **not a universal plate-energy upper bound**. The free finite operator excludes rigid motion; this separate directed-traction budget is not added a second time to its energy-normalized signed packets.

## Artifacts and next decision

`output/tails.csv`, `region_peaks.csv`, and `receiver_bank_union_tails.csv` carry the joint size/depth/height/duration and thermal readouts. `requirements.csv`, `attenuation_scan.csv`, `joint_frontier.csv`, `finite_event_requirements.csv`, and `finite_event_joint_frontier.csv` keep the conditional mechanisms separate. `whole_foil_histories.npz` stores all whole-volume strain-energy histories; the other region observables are in tables. The exact ledger and input hashes are in `inputs.json`, per-arm parameters in `models.csv`, and focused physical-metric, signed-cancellation and nonlinear-conservation evidence in their named files.

The decision-changing refinement is **the finite-duration D storage/holding kernel**, not more digits in an acoustic peak. To move the backwards requirement, the next mechanism must change energy supply/concentration or retained-D storage/transfer; a finer solver at the same fixed Hamiltonian energy cannot remove the positive quadratic bound. For the forward acoustic level, the first missing empirical input is band-resolved launch stress/spectrum and its acoustic loss, followed by mounted-boundary/grain texture. No new calibration or larger campaign is started here.

## Actual nonlinear redistribution and focused behavioral evidence

The final localized Hamiltonian uses **front-side bending epsilon_xx at depth 30 μm**, not weakly driven epsilon_zz at the foil's neutral plane. The latter preliminary projection barely changed the driven field; it was replaced rather than interpreted as a material-nonlinearity result. At the same operator, launch, counts and receivers, the final canonical strength labels .3/1 correspond to dimensionless quartic strain coefficients gamma≈4.38e11/1.46e12. These enormous unmeasured coefficients are normalized sensitivities, not inferred Pd third/fourth-order elastic constants. Gamma, unlike dimensionless chi scaled by source energy, is shared across the matched source cases.

At the representative PdD 2.5-keV receiver, chi=0/.3/1 gives peak **.28257/.27733/.26878 J/m³**, peak/mean **6.169/6.084/5.953**, four-mean time fraction **.02500/.02461/.02372**, and longest complete duration **18.8/18.7/18.5 μs**. The source-driven nonlinear projection variance is .01339/.01308/.01250. Other tested regions show increases as well as decreases: the **largest matched peak ratio anywhere in the finite scan is 1.02204**. It is an observed redistribution envelope, not a fitted concentration gain or an extrapolation to the required strength. It is far smaller than the conditional requirement factors and not evidence of a universal tail-shape change.

The unforced, unattenuated 80-μs witness changes modal quadratic-energy shares by L1 **.00703/.02183** for chi=.3/1; chi=0 changes only 1.05e−13. Its largest endpoint relative total-Hamiltonian drift is **2.84e−7**. Thus the nonlinear arm changes the trajectory and redistributes conserved energy; it does not simply add a quartic observable. A single actual effective D+ packet of **9.7342e−19 J** has peak 1.0675e−11 J/m³ at the named 2.5-mm, 5–55-μm receiver, changing only about 2e−12 fractionally at chi=1. The across-arrival result is therefore separate from single-parent packet change.

The documented full CLI exited naturally and wrote 24 model/source/strength rows, 1920 physical-receiver tail rows, 120 dimensionless union-tail rows, 240 finite-event requirement rows and 2160 finite-event joint rows. The actual downstream duration consumer also ran without replaying the producer while that consumer was introduced. The final default CLI includes it. Physical-metric error is **2.33e−15**, eigen residual **6.11e−12**, six rigid motions are removed, opposite signed propagated fields cancel with recorded residual zero, and the whole-volume receiver recovers modal elastic energy within the declared check. These are finite-model self-checks, not an independent continuum or material certificate. No parent numerical replay is needed.
