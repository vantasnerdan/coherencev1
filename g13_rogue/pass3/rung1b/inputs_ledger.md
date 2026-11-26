# Signed-source input ledger

All operators, parent semantics, stopping profiles, loaded-D populations and the random 1% law are inherited from accepted `../rung1/inputs_ledger.md`, `../rung1/rungs.py`, and `../../pass2/replication/{source_ledger.json,portable.c,replicate.py,inputs.npz}` at base `fe514a174214b01558c228fc9f8982ef3791fb06`. No measured conversion efficiency is introduced.

| Consumed input | Value / units | Publication / authority; assumption |
|---|---|---|
| Linear thermal expansion | 11.77e-6 K^-1 at 293.15 K | Arblaster, *Crystallographic Properties of Palladium*, Table I, [transcription](https://doczz.net/doc/7859361/crystallographic-properties-of-palladium), DOI 10.1595/147106712X646113. Pure Pd, not measured loaded PdD. |
| Heat capacity | 25.89 J mol^-1 K^-1 / 0.1064 kg mol^-1 = 243.327 J kg^-1 K^-1 at 298.15 K | Furukawa et al., NBS/J. Phys. Chem. Ref. Data, Table 6, [primary PDF](https://srd.nist.gov/jpcrdreprint/1.3253137.pdf). Constant near-room transfer; no temperature fit. |
| Conductivity | 0.718 W cm^-1 K^-1 = 71.8 W m^-1 K^-1 at 298.2 K | Ho, Powell and Liley, p.377 (PDF p.100), [primary PDF](https://srd.nist.gov/jpcrdreprint/1.3253100.pdf); well-annealed high-purity Pd, stated roughly ±5% near room temperature. |
| Density, stiffness | 12007 kg m^-3; C11/C12/C44=227.1/176.04/71.73 GPa | Accepted tensor/mass and seed-173 rotated grains; Rayne and Arblaster citations in inherited ledger. |
| Electrical parents | 138e-6 A, 5000 V; D+/D2+ probabilities .85/.15 | Accepted true-charge convention. Incident power .690 W. D2+ is one electrical parent, two correlated constituents at 2500 eV. |
| Incident momentum | p_D=sqrt(2 m_D 5000 eV); p_D2=sqrt(2) p_D | CODATA deuteron mass (SciPy constants); full-stop inward-transfer hypothesis, no reflection enhancement or random sign. Mean incident force about 2.1174 µN; transfer is not measured. |
| Deposited components | accepted species host, electron and direct-D cumulative tracks; terminal residual | Use actual cumulative increments and depth units, rather than assigning the host centroid to electrons. Primary assumption: all electronic, direct-D and residual energy retained; escape=0. Unresolved 100 eV/constituent is deposited at terminal depth as an explicit assumption. Alternative electron-escape effect is reported as a partition sensitivity, not fabricated data. |
| Heat / mechanical partition | host + electron + direct-D + residual minus actual random, pressure and thermoelastic work | Actual signed source work (including cross terms) is removed once from the heat reservoir; negative work returns energy. The random marked self-energy is a diagnostic, not a second heat subtraction. No bath/random acoustic energy is counted as new incident energy. |
| Work-debit location | J per bin, globally energy-consistent assumed apportionment | Random plus pressure work is distributed over mean beam-cell weights and the parent-probability-averaged normalized host-depth shape. Thermoelastic work is removed from the uniform thermal mode over those weights. Individual-parent/local work attribution is unresolved; this is not a measured local heat-partition law. |
| Thermal boundary | insulated front/back; thickness diffusion during finite exposure; negligible-lateral-diffusion approximation for beam-scale transient | Cosine spectral diffusion with exact cumulative-track source projection, including an undamped uniform thermal mode. No fictitious insulated DC solution. Separate lateral isothermal flange solution supplies a bounded mean-temperature/static-bias conditional. |
| Loss | retained gamma_E=1/(50 µs)+v_proxy/.2 m OR Cu gamma_E=4460 ln(10)/10 ×32 f/MHz | Separate alternatives, never added. Exact viscous modal propagation uses the same envelope energy-loss rates; Cu is explicitly an extrapolated prior attenuation proxy, not a material study or a measured Pd guided loss. |
| Static mechanical convention | force/torque-balanced elastic shape about CM; free rigid acceleration separate | The reported full-space static solve includes soft modes below 5 kHz. Only the bias subtracted from finite-band dynamic readouts is projected onto retained modes; it is not promoted to full DC compliance. An ideal isothermal thermal flange is not a mechanically calibrated clamp: a real clamped/sprung mount changes compliance, modes and ring-up. |
| D storage | unit transfer, holding time 4.1 ps | Accepted positive nonlocal sqrt(gamma)-weighted elastic/kinetic cycle-envelope proxy and fixed n_D(z); not instantaneous viscous heating, measured absorption or measured energy conversion. Actual work/loss accounting is separate. |

Computed diffusivity and thermal skin depth are dimensional diagnostics, not measurements. Temperatures, stopping retention, source size/duration and losses remain conditional. No PdD thermal refit, Ti, MD, electronic or force campaign is opened.

## Primary power partition

| Beam-averaged component | Nominal W at 138 µA | Credited disposition |
|---|---:|---|
| Full incident kinetic beam | 0.690000 | Ceiling, not an extra heat term |
| Electronic stopping | 0.586665 | Retained heat in the primary conditional model |
| Host recoil | 0.087123 | Shared reservoir for actual random/pressure work and remaining heat |
| Direct-D stopping | 0.000342 | Retained heat |
| Terminal cutoff remainder | 0.015870 | Unresolved; **assumed** locally retained, not measured heat |
| Electron/sputter escape | 0 in primary model | Explicit assumption, not an observed absence |
| Random marked self-energy | Per-record ledger | Diagnostic only; actual signed random work, not this diagnostic again, is removed from heat |
| Momentum / thermoelastic work | Per-record ledger | Actual signed work with interference; removed once from retained heat |
| Damping | Per-record ledger | Leaves the tracked mechanical-plus-heat reservoir; not re-credited as new beam power |

The electronic, host, direct-D and cutoff rows sum to the incident ceiling (rounding shown). The source normalizes central species depth shapes to the accepted integrated beam partition, then enforces exactly 5000 eV per electrical parent; it does not invent species-resolved measured heat-retention data.

If all electronic stopping escaped instead, the non-electronic available power would be about 0.103335 W before mechanical work, ~0.150 of the full-retention value. Thermal-only fields scale linearly with that retained heat under a **fixed random trajectory**; thermal-only quadratic energy scales as its square. This is a declared partition sensitivity, not a re-simulation of nonlinear beta .35 feedback.

The exact harmonic semi-infinite surface-**peak** temperature for a unit-peak sinusoidal 0.690-W modulation over 1 cm² is `q/(k sqrt(omega/D_th))`, about 0.601 mK at 100 kHz. `q delta/k≈0.850 mK` is the dimensional skin-length product, not that peak amplitude (it differs by sqrt(2)). The 64→128 cosine-mode check reports the **mechanically projected response** error separately from finite surface-temperature truncation; neither temperature is measured here.
