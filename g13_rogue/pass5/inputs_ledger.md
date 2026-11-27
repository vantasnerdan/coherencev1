# Pass-5 consumed inputs

| Input | Meaning and consumer | Authority / status |
|---|---|---|
| Base | `65dad6ce3e7920a4564b1456f8df2ab3e1c8fe5f`; only additive `g13_rogue/pass5/` | Axiom, Fleet Coherencev1 `g13_rogue` `entry_4184` |
| Plan | `coherencev1-g13-pass5-track-recoil-loss-20261127`; template `research-evaluation` at `14c5716f18cd5b1b8ee77624079ed31a0879c64388a625f4db88162dda02faf4` | Initial plan revision `b6ede3215fb447624d028c493415a5b4f400285e365089ea5fb4e3e60fe77e9e`; final record in `governing_plan.json` |
| Beam | 138 µA true electrical current, .690 W, 5 keV per electrical parent, D⁺/D₂⁺ .85/.15; 1.15 nuclei/parent | Frozen accepted beam composition and physical current |
| Tracks | D⁺ 63 rows, 128.944 nm, 408.153 fs; D₂⁺ constituent 58 rows, 81.871 nm, 327.218 fs | `g11_cmns/collective/part1/stopping.json`, PdD / 2.5-keV-CM central CSDA tracks |
| Row units | `depth_A` in Å and `time_fs` in fs; cumulative channels are finite-differenced including zero launch | The scalar `terminal_depth_A` is historically mislabeled and is already in metres |
| Channel energy | Electronic 4251.194688, host 631.324775, direct-D 2.480537, terminal 115 J/true C | Actual species energy vectors in exact Rung2 `signed_r2.npz`; summed to 5 keV/parent after accepted rescaling |
| Pd elasticity | ρ=12007 kg/m³, C11/C12/C44=227.1/176.04/71.73 GPa | Accepted near-300-K pure-Pd transfer, [Rayne](https://doi.org/10.1103/PhysRev.118.1545) |
| Expansion | αL=11.77×10⁻⁶ K⁻¹; eigenstress coefficient αL(C11+2C12) | [Arblaster 2012](https://doi.org/10.1595/147106712X646113), solid Pd near293 K |
| Heat capacity | 25.89 J/mol/K / .1064 kg/mol =243.327 J/kg/K | [Furukawa et al. 1974](https://doi.org/10.1063/1.3253137), Table6 |
| Conductivity | 71.8 W/m/K; Dth=κ/(ρCp) | [Ho, Powell, Liley1972](https://doi.org/10.1063/1.3253100), annealed pure Pd near300 K |
| Recoil kinematics | mD=3.3435837768×10⁻²⁷ kg; MPd=.1064 kg/mol divided by Avogadro's number; 1 eV=1.602176634×10⁻¹⁹ J | Two-body elastic identities, uniform azimuth; expectation normalization to actual row host energy |
| Screened marks / range | Existing `g4_beam/transport/bca.py` ZBL angle/stopping functions; Pd-on-Pd nuclear-only stopping, terminal .1 eV; range and quarter-range | Conditional direction/range model; no measured Pd recoil trajectories, azimuths, order or lateral track were supplied |
| Bath / thermalization | 2-nm Gaussian coherent-bath averaging radius; row-local electronic/direct-D heating, host work along finite recoil path | Declared mesoscopic source and instantaneous row-to-bath thermalization closure; CSDA transit is not the whole source/receiver duration |
| Foil and texture | 24-mm disk's equal-area square, 250 µm thick; illuminated 1 cm²; fixed8×8 seeded rotations | Accepted whole-foil geometry and source weights reused unchanged |
| Dispersion | Branch caps3.4/4.7/6.8 THz, sine interpolation preserving directional low-q speeds | Accepted conditional interpolation to pure-Pd branch evidence, [Miiller/Brockhouse](https://doi.org/10.1139/p71-087); 2.5 THz is the requested ladder end |
| Cu loss | αE=ln(10)/10 ×32(f/MHz) m⁻¹, applied to **actual path length** | [Ono2020](https://doi.org/10.3390/app10072230), Table3 N27, OFHC Cu longitudinal low-MHz data; unchanged extrapolated alternative |
| Al loss | τE=34 ps×240 GHz/f, applied to **actual transit time** | [Grossmann2017](https://doi.org/10.1088/1367-2630/aa6d05): Al longitudinal amplitude lifetime68 ps→energy34 ps at240 GHz, ~50–700 GHz thin-film fitted range; unchanged transfer alternative |
| Local D denominator | Fixed accepted loaded-D number density at63 nm,30/50/200 µm | Elastic energy density divided by local number density is an inventory quotient, not acquired microscopic D energy |
| Slow background | Exact accepted low-frequency histories and two loss alternatives | Original bytes copied/reused; high-frequency state does not refit the approved slow producer |
| Old source comparisons |20 nm/.5 ps and60 nm/2 ps Gaussian launch arrays; historical .1 host share's .99 high-frequency reservoir plus its separate .01 low-band permission | `old_launch.npz/csv` are exact unchanged Pass-4 bytes; fractions use the same incident-power denominator, not a per-band host fraction |
| Feedback | Exact Pass2 β=.35 pre-existing-state mark driver, now declared source velocity/displacement projection | `pass2/live_method.md:9–21`, `pass2/live.py:391–459`; velocity-aware signed work and physical relative-speed arm are separate |
| Thermal comparison | Accepted same-band300-K elastic Bose energy | `thermal_reference.csv` unchanged; not total thermal energy or a tail probability |

No transferable numerical Pd/PdHₓ GHz loss was found in the bounded inspected sources, so neither alternative is relabeled as measured Pd: [Tanigaki2011](https://www.jstage.jst.go.jp/article/jsmeatem/2011.10/0/2011.10__OS02F016-/_pdf/-char/en) gives arbitrary-unit echoes in95.1-nm Pd/Si without a fitted/deconvolved loss, [Youngblood1979](https://doi.org/10.1103/PhysRevB.19.6016) numerical figures were inaccessible (the [OSTI abstract](https://www.osti.gov/biblio/6007299) is qualitative), and [Leisure1983](https://doi.org/10.1088/0305-4608/13/10/013) concerns loaded PdH₀.₆₄ at10–270 MHz with no numerical attenuation in its accessible abstract.

For a future compatible Lorentzian **energy FWHM** ΔE, γE=ΔE/ℏ, τE=ℏ/ΔE and τA=2τE; HWHM must first be doubled. Only the actual measured branch, frequency window, loading, temperature and instrumental deconvolution would transfer. No linewidth or flat/zero continuation is inserted here.
