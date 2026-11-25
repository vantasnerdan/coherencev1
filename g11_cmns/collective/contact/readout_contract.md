# Conditional contact and G12 readout: PdD, 0.25 keV CM

## Preparation and meaning of contact

The ordinary concentrated preparation is `cold_localized_1.0` in [Part 1 results](../part1/results.json#L114-L129). [The producer](../part1/motion.py#L265-L278) displaces only the G10 source site (original index **D25**) along **+x**, starting from the model minimum. It starts at rest and evolves for 2 ps without drive or damping. The site displacement is **0.11988147186050971 Å**, matched to the **whole lattice** excess energy **0.004844094171855048 eV** of the uniform Γ control at 0.037968 Å. This is not the saved pump or the separate `cold_localized_2.0` sensitivity (Γ amplitude 0.075936 Å, site displacement 0.2223309455834531 Å, whole-lattice excess energy 0.021693574626444 eV).

The native zero-image vector from D25 to D63 is along **+[110]**, with length \(R_0=2.88782409436586\) Å. Projecting the D25 displacement onto that direction gives an inward relative displacement **\(\delta=0.08476900169119068\) Å** (\(s/\sqrt{2}\)), or equal/opposite path coordinate \(u=\delta/2=0.04238450084559534\) Å. The corresponding projected radial point is \(R_0-\delta=2.8030550926746693\) Å. The actual one-sided +x displacement also has a transverse component, and its exact distance to D63 is 2.804336576842573 Å, not the projected radial point. The archived ordinary global minimum uses the tied D25–D47 pair at \(t=0\); it does not select D25–D63 uniquely. For the doubled-Γ test the inward projection is 0.15721171928967695 Å, not the ordinary input.

A concentrated two-level radial preparation can match **only** the path work \(U_{\rm lat}(R_0-\delta)-U_{\rm lat}(R_0)\) via the real ground/first-excited mixture. When the work is between zero and the level gap, the excited-state weight is work/gap; choose the relative phase for an inward initial mean displacement. A real wavefunction gives zero initial current, but energy matching does not also fix its initial mean radius to \(R_0-\delta\). Such a state, and any 77/300-K radial thermal mixture, are one-path conditional model states, not the original 32-site trajectories or a reconstructed thermal ensemble. Part 1 thermal histories used classical canonical sampling and physical Maxwell momenta.

`rate_readout` takes the base-10 log of the three-dimensional contact density in Å⁻³ for **one normalized tagged pair** in the assumed pure spatial \(l=0\) sector. For a normalized reduced radial wavefunction \(u_R\) (\(\int_0^\infty |u_R(R)|^2\,dR=1\)), the 3D density is \(\rho=|u_R'(0)|^2/(4\pi)\); no second angular factor is applied here. This does not determine an \(l=0\) projection weight for an oriented lattice state. An ensemble or time-window density must be formed upstream, on the density scale, not by averaging logarithms.

## Fixed conversion and reference

For the allowed-spin statistical \(s\)-wave DD mixture (\(S=0:2\) in the ratio 1:5), \(\Gamma_p=A_p\rho\). [Part 1](../part1/motion.py#L346-L351) uses \(A_p=5.81220367098898\times10^7\) Å³/s, from \((3/4)S_p(0)/(\pi\alpha\mu c_{\rm light})\) with SI-to-Å³ conversion, \(\mu=m_D/2\), and Bosch–Hale \(S_p(0)=55.576\) keV barn. [Pairlaw `coefficient('DD')`](../../../g6_pairlaw/script/pairlaw.py#L41-L59) sums that intercept **and** the 53.701-keV-barn neutron intercept; its total coefficient is not a proton coefficient. No random-nine-spin, threshold-resonance, or attempt-frequency multiplier is silently added.

These are the archived [.25-keV PdD conversion rows](../part1/rates.csv#L14-L15), using the [fixed G12 inputs ledger](../../../g12_residual/inputs_ledger.md) and row `karahadian2026:fig2:PdD:0.25` in `g12_residual/regression/davis_residuals.csv`. The G12 own-foil modeled n(6.5 keV) effective readout is \(N_{65}=434741.4188244709\) reactions/C. The primary p/n join is 1. The retained focusing factor at .25 keV is \(f_{\rm cup}=I_{\rm cup}/I_{\rm true}=1\); the chart model gives \(I_{\rm cup}=I_{\rm true}=38\) µA (\(3.8\times10^{-5}\) C/s), versus \(I_{\rm cup}(6.5)=238\) µA. The chart-model current is not an actual run-current measurement.

| Archived scalar | Coefficient per unit \(\rho\) | Coefficient per unit \(\Xi\) |
|---|---:|---:|
| Emitted p/s | 58122036.7098898 | 1.5853014894375147e22 |
| Emitted p/true C | 1529527281839.2053 | 4.1718460248355654e26 |
| Added G12 ordinate \(\Delta R\), k0 | 3518246.0552643132 | 9.596154965211561e20 |
| Added G12 ordinate \(\Delta R\), k1 | 22035330.556655437 | 6.010223372948294e21 |

For a continuously active pair, \(Y_p=\Gamma_p/I_{\rm true}\) is conditional p/true C. For either hypothesis the addition is \(\Delta R_{k0}=Y_p/(f_{\rm cup}N_{65})\) and \(\Delta R_{k1}=(238/38)\,\Delta R_{k0}\), with focusing once, no second division by current and no change to the observed ordinate. These are model additions, not the adjacent residuals: the archived model ordinates are \(2.277535921907689\times10^{-19}\) (k0), \(1.4264567089842895\times10^{-18}\) (k1), versus \(R_{\rm obs}=0.05052\); both adjacent \(A_k\) values round to 0.05052.

`event_readout` is a **different hypothesis**. Supply the log of the positive excess \(\Xi=\int_0^{2\,{\rm ps}}[\rho_{\rm prepared}(t)-\rho_{\rm matched\ radial\ ground}(t)]\,dt\) in Å⁻³ s, after **one** subtraction. Then \(A_p\Xi\) is excess emitted p/event. For one such event per incident D, with the same excess in both beam components, \(Y_p=(0.85+2\cdot0.15)A_p\Xi/e\) p/true C (two D per D₂⁺ charge); p/s uses \(I_{\rm true}Y_p\). This does not infer a capture or preparation probability. If a modeled excess is nonpositive, do not take its positive-event log or present it as an added emission.

## API and limits

`rate_readout(log10_contact_A3)` and `event_readout(log10_excess_Xi_A3_s)` both return `log10_emitted_p_per_s`, `log10_emitted_p_per_true_C`, `log10_delta_R_k0`, and `log10_delta_R_k1`. They take finite base-10 log magnitudes of positive scalars in the stated units, not linear inputs. Convert from natural logs by division by \(\ln 10\). All operations are additive in log space; a tiny nonzero contact need not be exponentiated or rounded to zero. An exact zero has no finite log and is outside these positive-input conversions.

Neither hypothesis specifies a real active-pair/domain population, angular projection, preparation duty cycle or capture fraction. Neither estimates a foil rate. The independent [SI1 detected-count proxy](../../../g12_residual/regression/SI1_normalizations.csv#L8) is 5.833333333333333e−5 detected/s (1.5350877192982457 detected/model true C); detector acceptance and the absolute observed 6.5-keV reference are unknown. Do not compare by treating this proxy as emitted protons, subtracting it from the modeled emission, or using it for calibration. No inference is made for other foils.
