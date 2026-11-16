"""Write experiments.csv (one row per experiment). Values are transcribed from the sources named in each row.
'est.' marks my own estimates (see compare.py / compare_out.txt). 'verification' says how much of the source was read.
Run: python3 -I make_experiments_csv.py <output csv>
"""
import csv
import sys

COLS = ["id", "reference", "group", "reaction", "target", "method", "E_cm_lowest_keV", "E_cm_range_keV",
        "beam_species_purity", "energy_definition", "vacuum", "charge_measurement", "detection",
        "Ue_eV_reported", "low_energy_behaviour", "attributed_cause", "digitizable_data", "verification"]

R = []
R.append(dict(
    id="karahadian2026", reference="Karahadian, Colborne, Persaud, Schenkel, Munday, Nat. Commun. 17:8845 (2026); preprint arXiv:2512.06212v1 (5 Dec 2025)",
    group="UC Davis / LBNL", reaction="d(d,p)t; d(d,n)3He", target="250 um Pd or Ti foil, 200 nm sputtered Pd/Ti overlayer; beam-only vs electrochemically loaded (D2SO4/D2O)",
    method="microwave ECR gun (Tectra IonEtch), 3-grid single-aperture extraction, no mass/charge analysis, 3 mm aperture ~4 cm from target",
    E_cm_lowest_keV="0.25", E_cm_range_keV="0.25-6.5",
    beam_species_purity="85% D+ / 15% D2+ (stated from source literature; not measured per energy)",
    energy_definition="E_lab = gun voltage (<=5 kV, constant extraction voltage) or gun + negative target bias (>=5.5 kV); spread ~5% from retarding-field scan at 2.5-4.2 kV only",
    vacuum="base <=1e-6 Torr; D2 2.5-3e-4 Torr during runs",
    charge_measurement="charge computed from a separately measured Faraday-cup IV curve (not integrated on the foil during runs)",
    detection="Si diode behind 10 um Al (protons, E_cm 0.25-4.5); EJ-309 (neutrons, E_cm >=3-4.5); detector angle not stated",
    Ue_eV_reported="Pd 800+-200; PdD 1700+-200; Ti 300+-200; TiD 900+-200 (fit 2.5-6.5 keV); per-energy Ue rising to 3.0 keV (PdD), 2.2 keV (TiD) at 0.25 keV",
    low_energy_behaviour="floor below ~2 keV (v1: below 2.5 keV with 'slight increase' at lowest energies); est. floor/Y(6.5): Pd 4.1e-3, PdD 1.4e-2, Ti 1.2e-3, TiD 1.9e-3",
    attributed_cause="heterogeneous (surface vs bulk) defect-enhanced screening; low-energy resonances mentioned as alternative",
    digitizable_data="Fig 2 (relative yields); Fig 3d (hydride/metal ratios); SI Table 1 (S:B, currents); SI Table 2; Zenodo 10.5281/zenodo.20435774 embargoed to 2026-12-01",
    verification="full text + SI read"))
R.append(dict(
    id="colborne2025_poster", reference="Colborne, Karahadian, Unzueta, Johnston, Persaud, Munday, Schenkel, 'Quantifying Nuclear Reactions in Metal Hydrides at Low Energies', APS March Meeting 2025 poster; OSTI 2569797; eScholarship 6xk1h758 (2025-03-16)",
    group="LBNL / UC Davis (same apparatus)", reaction="d(d,n)3He; d(d,p)t", target="0.25 mm Pd and Ti foils; beam-loaded vs electrochemically loaded",
    method="same ion gun; beam <=15 kV; electrochemistry only up to 5 kV beam at that time",
    E_cm_lowest_keV="2.5", E_cm_range_keV="2.5-7.5 (neutrons)", beam_species_purity="not stated", energy_definition="not stated",
    vacuum="not stated", charge_measurement="plotted as neutron rate per hour (not per charge)",
    detection="EJ-309 (bg ~100 n/h; signal 10-20 /h at low energies); Si diode (bg 0.5 counts/h; proton data only below 9 kV because of x-ray noise)",
    Ue_eV_reported="none",
    low_energy_behaviour="neutron rates of all four materials converge to ~40-100 /h at E_cm 2.5-3.0 keV (my reading; probably background-dominated); text: Ti rate >> Pd; 'no increase' for electrochemically loaded Pd",
    attributed_cause="n/a",
    digitizable_data="OSTI PDF p.3 'Initial Results': neutron rate vs E_cm (Ti, TiDx, Pd, PdDx), proton spectra 5-9 kV (Ti); est. TiDx/Ti = 2.1 and PdDx/Pd = 0.43 at 6.5 keV; TiDx/Ti relative to 6.5 keV reproduces paper Fig 3d neutron points",
    verification="poster read and plot inspected"))
R.append(dict(
    id="czerski_zr_2409", reference="Czerski et al., arXiv:2409.02112 v1 (21 Aug 2024), v2 (11 Dec 2025), v3 (28 Sep 2026)",
    group="Univ. Szczecin (eLBRUS)", reaction="d(d,p)t", target="1 mm Zr implanted to ZrD2",
    method="ECR source, 90 deg analysing magnet, 2-lens deceleration (source at 12-20 kV), target at ground",
    E_cm_lowest_keV="0.625", E_cm_range_keV="0.625-12.5", beam_species_purity="D1+, D2+, D3+ magnetically separated",
    energy_definition="source potential; 10 eV definition, 5 eV stability",
    vacuum="target chamber 1e-10 mbar; beamline 5-8e-9 mbar; source 1e-6 mbar",
    charge_measurement="on target plus two Faraday cups; secondary electrons +10-15% checked with +60 V target bias",
    detection="PIPS Si at 135 deg, 5 cm, 4 um Al, 20 keV resolution",
    Ue_eV_reported="340+-30",
    low_energy_behaviour="plateau below E_d ~3 keV (E_cm 1.5) at 2.6e-4 counts/(uA ks); est. 3e-6 of Y(6.5), 5e-4 of Y(3.0)",
    attributed_cause="thermal D-D fusion in ion tracks (thermal spike) + screening + 0+ threshold resonance; with D1+/D2+ through the decelerator a neutral-beam component (fraction ~3e-4; 1e-4 in v1) made its own plateau, identified by lower proton energy at 135 deg and suppressed with D3+ at 12 kV",
    digitizable_data="v3 Fig 2a,b (yield and proton peak vs energy, D1+/D2+/D3+); Fig 7a,b (D1+ with neutral component); Table 2; Fig 1 spectra",
    verification="v1 and v3 full text read"))
R.append(dict(
    id="czerski_ti", reference="Czerski et al., arXiv:2605.27438 (22 May 2026); arXiv:2409.02112v3",
    group="Univ. Szczecin", reaction="d(d,p)t", target="1 mm Ti0.993Pd0.007 implanted to D/M ~2",
    method="as czerski_zr (deceleration used)", E_cm_lowest_keV="~0.75 (E_d 1.5 keV quoted)", E_cm_range_keV="~0.75-10",
    beam_species_purity="D1+, D2+, D3+ analysed; lowest points with D3+", energy_definition="as czerski_zr", vacuum="as czerski_zr",
    charge_measurement="as czerski_zr", detection="as czerski_zr",
    Ue_eV_reported="25+-10 (threshold-resonance excess visible below E_d 6 keV)",
    low_energy_behaviour="plateau below E_d 2-3.5 keV at 1.4e-4 counts/(uA ks); est. 2.4e-6 of Y(6.5), 7.7e-4 of Y(3.0)",
    attributed_cause="thermal spike; hot tracks assumed metallic with Ue ~370 eV; low cold Ue blamed on oxygen-rich, embrittled surface",
    digitizable_data="2409 v3 Fig 3a,b; 2605.27438 Fig 2; Table 2",
    verification="full text read"))
R.append(dict(
    id="czerski_pd", reference="Czerski et al., arXiv:2409.02112v3; arXiv:2605.27438",
    group="Univ. Szczecin", reaction="d(d,p)t (+ e+e- search)", target="1 mm Pd plate, D/Pd ~0.05 (beam-loaded)",
    method="ECR + magnet, no deceleration (so no above-nominal neutrals)", E_cm_lowest_keV="~1 (read from figure axis; not stated)", E_cm_range_keV="~1-10",
    beam_species_purity="D1+, D2+, D3+ analysed", energy_definition="source potential, 10 eV",
    vacuum="1e-10 mbar target chamber", charge_measurement="on target; beam power held at ~0.6 W (target 30-32 C; D density varied ~15%)",
    detection="1000 um Si at 135 deg, 6 cm, 1 um Al; two proton components resolved",
    Ue_eV_reported="400+-40",
    low_energy_behaviour="plateau below E_d ~6 keV (E_cm 3) at 6.9e-3 counts/(uA ks); est. 2.8e-3 of Y(6.5), 0.44 of Y(3.0) (cf. Karahadian Pd 4.1e-3 and 0.22)",
    attributed_cause="thermal fusion: a proton component at fixed (E_cm = 0) energy with energy-independent rate",
    digitizable_data="2409 v3 Fig 4a,b and Fig 1b; 2605.27438 Fig 1, Fig 3; Table 2",
    verification="full text read"))
R.append(dict(
    id="czerski2016_kaczmarski2017", reference="Czerski et al., EPL 113, 22001 (2016); Kaczmarski et al., Acta Phys. Pol. B 48, 489 (2017)",
    group="Univ. Szczecin / Berlin", reaction="d(d,p)t", target="deuteron-implanted Zr, UHV 1e-10 mbar, Ar-sputter cleaned",
    method="accelerator (eLBRUS predecessor)", E_cm_lowest_keV="not verified (E_d < 25 keV stated)", E_cm_range_keV="E_d < 25 keV",
    beam_species_purity="analysed", energy_definition="accelerator voltage", vacuum="1e-10 mbar",
    charge_measurement="not checked", detection="Si",
    Ue_eV_reported="Zr 490 (UHV, earlier) vs 300 (HV); ~105+-15 when resonance included",
    low_energy_behaviour="enhancement above screening fit below E_d ~10 keV (no plateau)",
    attributed_cause="0+ threshold resonance in 4He interfering with broad resonances; surface impurities raise Ue",
    digitizable_data="Kaczmarski 2017 Fig 1 (enhancement factors normalised at E_d = 15 keV)",
    verification="Kaczmarski full text; Czerski 2016 abstract only"))
R.append(dict(
    id="czerski2001_huke2008", reference="Czerski et al., EPL 54, 449 (2001); Huke et al., PRC 78, 015803 (2008), arXiv:0805.4538",
    group="TU Berlin", reaction="d(d,p)t; d(d,n)3He", target="implanted Al, Zr, Ta (2001); also Pd, Sr, Li, Na, C (2008)",
    method="RF source, 60 kV, magnetic analysis, magnetic steerer before target (removes neutrals), SE suppression cylinder",
    E_cm_lowest_keV="2.5", E_cm_range_keV="2.5-30", beam_species_purity="analysed D+", energy_definition="60 kV stabilised supply",
    vacuum="HV", charge_measurement="target with suppression; yield recorded vs implanted charge Y(q) to track D density",
    detection="Si detectors at several angles",
    Ue_eV_reported="Ta 322+-15; Zr 297+-8; Al 190+-15; Pd 313+-2; Sr 350-800; Li, Na <=150; C ~0",
    low_energy_behaviour="no floor in range",
    attributed_cause="screening larger than theory; Huke lists oxide layers, D-density dynamics and beam heating as sources of spurious screening in others' data",
    digitizable_data="Huke 2008 Table I; Fig 11 (Ue vs D density, all groups)",
    verification="Huke full text; Czerski 2001 abstract"))
R.append(dict(
    id="kasagi2002", reference="Kasagi, Yuki, Baba, Noda, Ohtsuki, Lipson, J. Phys. Soc. Jpn. 71, 2881 (2002); values tabulated in Kasagi, Prog. Theor. Phys. Suppl. 154, 365 (2004)",
    group="Tohoku", reaction="d(d,p)t", target="Ti, Fe, Pd, PdO, Au (beam-loaded, LN2-cooled holder)",
    method="low-energy high-current ion beam", E_cm_lowest_keV="1.25", E_cm_range_keV="1.25-5",
    beam_species_purity="D+ (analysis not checked)", energy_definition="accelerating voltage", vacuum="not checked",
    charge_measurement="beam dose; yields normalised to E_d = 10 keV with repeated monitor runs", detection="one Si detector at 90 deg",
    Ue_eV_reported="PdO 600+-30; Pd 310+-30; Fe 200+-20; Au 70+-20; Ti 65+-30",
    low_energy_behaviour="steep decrease to E_d 2.5 keV; PdO yield ~50x bare at E_d 2.5 keV; no floor",
    attributed_cause="screening plus an extra solid-state mechanism; Huke 2008 attributes part to D-density dynamics, oxides, heating and 10 keV normalisation",
    digitizable_data="excitation-function figure in JPSJ 71, 2881 (figure number not checked)",
    verification="abstract + secondary (Huke 2008, review table)"))
R.append(dict(
    id="yuki1997", reference="Yuki et al., J. Phys. G 23, 1459 (1997); J. Phys. Soc. Jpn. 66, 73 (1997); EXFOR O1084",
    group="Tohoku", reaction="d(d,p)t", target="thick Ti (pre-loaded with D2 gas); thick Yb (beam-loaded)",
    method="high-current ion beam generator", E_cm_lowest_keV="1.28 (Ti)", E_cm_range_keV="Ti 1.28-3.12; Yb 1.61-3.58",
    beam_species_purity="D+ (not checked)", energy_definition="accelerating voltage", vacuum="LN2 trap",
    charge_measurement="yield per dose normalised to E_d = 6.45 keV", detection="dE-E Si telescope (50 + 200 um)",
    Ue_eV_reported="(see Kasagi 2002 for Ti 65+-30)",
    low_energy_behaviour="no floor: Ti falls to 4.1e-5 of the E_d 6.45 keV yield at E_cm 1.28; est. Ti Y(1.53)/Y(3.12) = 5.0e-4, equal to bare thick-target model (Karahadian Ti: 0.12)",
    attributed_cause="n/a (authors note target D number changes during bombardment)",
    digitizable_data="EXFOR O1084 subentries 002 (Yb) and 003 (Ti): tabulated, digitised from Figs 1-2",
    verification="EXFOR entry read"))
R.append(dict(
    id="yuki1998", reference="Yuki et al., JETP Lett. 68, 823 (1998)", group="Tohoku", reaction="d(d,p)t",
    target="Pd; Au/Pd/PdO heterostructure (60 um)", method="as kasagi2002", E_cm_lowest_keV="~1.25 (not checked)", E_cm_range_keV="not checked",
    beam_species_purity="D+", energy_definition="accelerating voltage", vacuum="not checked", charge_measurement="not checked", detection="Si",
    Ue_eV_reported="Au/Pd/PdO 602+-23 (as quoted by Huke 2008)", low_energy_behaviour="no floor reported",
    attributed_cause="screening; Huke attributes high value to heating of thin target and density shifts",
    digitizable_data="not checked", verification="secondary only"))
R.append(dict(
    id="raiola2002_2004", reference="Raiola et al., Phys. Lett. B 547, 193 (2002); Eur. Phys. J. A 19, 283 (2004)",
    group="LUNA collaboration, Bochum 100 kV", reaction="d(d,p)t", target="58 samples: metals, insulators, semiconductors, beam-loaded",
    method="analysed beam", E_cm_lowest_keV="2.5", E_cm_range_keV="2.5-15 (E_d 5-30 keV per Huke 2008)",
    beam_species_purity="D+", energy_definition="accelerator voltage", vacuum="HV",
    charge_measurement="beam integration; D content from normalisation to gas cross section at E_d = 30 keV", detection="4 Si detectors at 130 deg",
    Ue_eV_reported="Pd 800+-90 (D/Pd ~0.03); large in metals, gas-like in insulators and lanthanides",
    low_energy_behaviour="no floor in range", attributed_cause="Debye screening by quasi-free electrons",
    digitizable_data="EPJA 19, 283 Table 1", verification="abstract + secondary"))
R.append(dict(
    id="raiola2005", reference="Raiola et al., J. Phys. G 31, 1141 (2005); EXFOR O1379",
    group="LUNA / Bochum 100 kV", reaction="d(d,p)t", target="Pt, Co, Ti, Sc, Y, Zr, Hf, lanthanides at -10 to 340 C",
    method="analysed beam", E_cm_lowest_keV="3.85", E_cm_range_keV="3.85-22", beam_species_purity="D+", energy_definition="0.7% cross-section uncertainty at E_d 4 keV from energy",
    vacuum="HV", charge_measurement="beam integration", detection="Si detectors cooled to 0 C",
    Ue_eV_reported="Ti: <30 (-10 C, D/Ti 2.1); <50 (50 C, 1.1); 250+-40 (100 C, 0.26); 295+-40 (150 C); 290+-65 (200 C, 0.20). Pt 675 (20 C) to 465 (300 C)",
    low_energy_behaviour="no floor", attributed_cause="Ue appears only when D solubility is low; Ue falls as T^-1/2",
    digitizable_data="EXFOR O1379: tabulated S-factors vs T (Pt, Hf) and summary table of Ue, solubility",
    verification="EXFOR entry read"))
R.append(dict(
    id="greife1995", reference="Greife, Gorris, Junker, Rolfs, Zahnow, Z. Phys. A 351, 107 (1995); EXFOR O1171",
    group="Bochum / LUNA 50 kV", reaction="d(d,p)t; d(d,n)3He", target="windowless D2 gas (also self-implanted Ti layer for n/p ratio)",
    method="analysed beams; molecular beams at lowest energies", E_cm_lowest_keV="1.62", E_cm_range_keV="1.62-130",
    beam_species_purity="analysed D+, D2+ (Coulomb-explosion correction applied)", energy_definition="effective E_cm corrected for loss, straggling, Doppler",
    vacuum="gas target", charge_measurement="beam calorimeter", detection="Si detectors (0.85 um Al)",
    Ue_eV_reported="~25 (gas; as quoted by Czerski/Huke)",
    low_energy_behaviour="no floor: S-factor 58-62 keV b from 2 to 15 keV, 74+-11 at 1.62 keV",
    attributed_cause="atomic screening",
    digitizable_data="EXFOR O1171 subentry 007 = paper Table 4 (S(E), 1.62-14.95 keV, tabulated); subentry 006 n/p ratio in self-implanted Ti, 2.4-13.1 keV",
    verification="EXFOR entry read"))
R.append(dict(
    id="bystritsky2012", reference="Bystritsky et al., Nucl. Phys. A 889, 93 (2012); Phys. At. Nucl. 75, 913 (2012)",
    group="JINR Dubna / Tomsk", reaction="d(d,n)3He", target="ZrD2, TiD2, TaD0.5 films",
    method="pulsed Hall plasma accelerator (broad ion energy distribution)", E_cm_lowest_keV="~2 (catalogue)", E_cm_range_keV="2-7",
    beam_species_purity="plasma beam, not analysed", energy_definition="time-of-flight / accelerator setting (not checked)", vacuum="not checked",
    charge_measurement="not checked", detection="plastic scintillators (2.5 MeV n)",
    Ue_eV_reported="ZrD2 205+-35; TiD2 125+-34; TaD0.5 313+-58; enhancement independent of T 20-200 C",
    low_energy_behaviour="no floor reported", attributed_cause="screening",
    digitizable_data="not checked", verification="catalogue abstracts only"))
R.append(dict(
    id="bystritsky_pd_d3he", reference="Bystritsky et al., arXiv:1306.6783, NIM A (doi 10.1016/j.nima.2013.11.072): d(p,g)3He in ZrD2; d(3He,p)4He in ZrD (Tomsk)",
    group="JINR / Tomsk", reaction="d(p,g)3He; d(3He,p)4He", target="ZrD2; ZrD", method="pulsed Hall accelerator",
    E_cm_lowest_keV="6.67 (p+d); 7.2 (d+3He)", E_cm_range_keV="6.67-12.67; 7.2-12.0",
    beam_species_purity="plasma beam", energy_definition="not checked", vacuum="not checked", charge_measurement="not checked", detection="not checked",
    Ue_eV_reported="d+3He in ZrD: 618+-155 (search snippet); p+d consistent with LUNA gas",
    low_energy_behaviour="no data below 5 keV", attributed_cause="n/a",
    digitizable_data="not checked", verification="abstracts only"))
R.append(dict(
    id="lipson2005", reference="Lipson, Rusetskii, Karabut, Miley, JETP 100, 1175 (2005)",
    group="IPC RAS Moscow / Illinois", reaction="d(d,p)t", target="Ti cathode",
    method="high-current pulsed glow discharge in D2 (300-600 mA/cm2)", E_cm_lowest_keV="0.4", E_cm_range_keV="0.4-1.23",
    beam_species_purity="not analysed (plasma ions; species mix unknown)", energy_definition="discharge voltage 0.8-2.45 kV",
    vacuum="D2 discharge", charge_measurement="discharge current", detection="3 MeV protons",
    Ue_eV_reported="Ti 610+-150",
    low_energy_behaviour="yield falls with voltage (fit by Ue); est. implied drop 6e-4 from E_cm 1.23 to 0.4 keV; no floor reported; ~1e9 enhancement at E_d 1 keV quoted in Lipson 2007 slides",
    attributed_cause="screening; enhancement depends on D diffusion in thin subsurface layer",
    digitizable_data="yield vs voltage figure (number not checked)", verification="abstract + secondary"))
R.append(dict(
    id="schenkel2019", reference="Schenkel et al., J. Appl. Phys. 126, 203302 (2019); arXiv:1905.03400",
    group="LBNL (Google-funded)", reaction="d(d,n)3He", target="Pd (and Ti) wire cathodes",
    method="pulsed glow discharge, 0.5-2 Torr D2", E_cm_lowest_keV="~1.2 nominal (lower after sheath losses)", E_cm_range_keV="~1.2-6.3",
    beam_species_purity="not measured; PIC: D+ with nearly equal D3+ flux", energy_definition="cathode voltage 2.4-12.6 kV minus sheath losses (~1.2 keV mean, PIC)",
    vacuum="discharge", charge_measurement="discharge current, not corrected for secondary electrons (factor-2 uncertainty)",
    detection="liquid scintillator with PSD (4 analyses differ x6 at lowest voltage); H2 control runs",
    Ue_eV_reported="Pd 1000+-250",
    low_energy_behaviour="yields 160-1000x bare at lowest voltage; no plateau claimed",
    attributed_cause="screening dependent on loading and defects",
    digitizable_data="Fig 9 (relative neutron yield vs cathode voltage)", verification="full arXiv text read"))
R.append(dict(
    id="chen2025", reference="Chen et al., Nature 644, 640 (2025); SI via Research Square rs-5046357",
    group="UBC (Berlinguette) with LBNL", reaction="d(d,n)3He", target="Pd foil, electrochemically loaded from back side (D/Pd ~0.7 electrolyte side)",
    method="microwave plasma + plasma-immersion, -30 kV sheath", E_cm_lowest_keV="15", E_cm_range_keV="15 (single energy)",
    beam_species_purity="D+ stated", energy_definition="sheath voltage -30 kV", vacuum="plasma chamber",
    charge_measurement="plasma current 0.5 mA", detection="EJ-309 at 12 cm, PSD, bg 0.2 n/s; beam on/off cycling",
    Ue_eV_reported="none", low_energy_behaviour="n/a; +15+-2% neutron rate with electrochemical loading",
    attributed_cause="higher D density", digitizable_data="main-text neutron rate vs time figures", verification="SI preprint read; main text not accessed"))
R.append(dict(
    id="cvetinovic2023", reference="Cvetinovic et al., Phys. Lett. B 838, 137684 (2023)",
    group="JSI Ljubljana et al.", reaction="1H(7Li,a)4He; 1H(19F,ag)16O; 2H(19F,p)20F (inverse kinematics)",
    target="two H/D-containing Pd foils", method="accelerator, heavy-ion beams", E_cm_lowest_keV="not checked", E_cm_range_keV="not checked",
    beam_species_purity="analysed", energy_definition="accelerator", vacuum="not checked", charge_measurement="not checked", detection="not checked",
    Ue_eV_reported="one foil: no large enhancement; other: up to ~10x adiabatic limit (talk slides: 2H(19F,p) 18.2+-3.3 keV 'hard' vs 3.2+-1.9 keV 'soft' Pd; unverified)",
    low_energy_behaviour="not d+d", attributed_cause="lattice position of hydrogen (non-equilibrium sites, defects)",
    digitizable_data="not checked", verification="abstract only"))
R.append(dict(
    id="dubey2025", reference="Dubey et al., Phys. Rev. X 15, 041004 (2025); arXiv:2408.07567",
    group="Univ. Szczecin", reaction="d+d -> p, n, e+e- (4He 0+ threshold resonance)", target="ZrD2",
    method="eLBRUS accelerator", E_cm_lowest_keV="2.5 (E_d 5 keV)", E_cm_range_keV=">=2.5",
    beam_species_purity="analysed", energy_definition="source potential", vacuum="UHV", charge_measurement="target",
    detection="Si of several thicknesses, NaI(Tl), HPGe",
    Ue_eV_reported="n/a", low_energy_behaviour="e+e- branching rises as energy falls; Gamma_ee >= 10 Gamma_p",
    attributed_cause="0+ threshold resonance", digitizable_data="branching-ratio figures (not checked)", verification="abstract only"))
R.append(dict(
    id="kitamura1996", reference="Kitamura, Saitoh, Itoh, Fusion Sci. Technol. 29, 372 (1996)",
    group="Kobe Univ.", reaction="d(d,p)t", target="Pd, Ti", method="keV D implantation with in-situ ERDA",
    E_cm_lowest_keV="not checked", E_cm_range_keV="keV", beam_species_purity="not checked", energy_definition="not checked", vacuum="not checked",
    charge_measurement="not checked", detection="charged-particle detectors",
    Ue_eV_reported="n/a", low_energy_behaviour="yield depends on beam current (flux) through temperature-dependent D concentration",
    attributed_cause="flux-dependent target D density", digitizable_data="not checked", verification="abstract only"))
R.append(dict(
    id="theory_rows", reference="Czerski, PRC 106, L011601 (2022); Kilic et al., EPJ Web Conf. 154, 01024 (2017)",
    group="Univ. Szczecin", reaction="d(d,p)t analysis", target="re-analysis of metal and gas data", method="theory (T-matrix, threshold resonance)",
    E_cm_lowest_keV="n/a", E_cm_range_keV="n/a", beam_species_purity="n/a", energy_definition="n/a", vacuum="n/a", charge_measurement="n/a", detection="n/a",
    Ue_eV_reported="n/a", low_energy_behaviour="no new data", attributed_cause="0+ threshold resonance; material dependence via local electron density",
    digitizable_data="none", verification="abstracts only"))
R.append(dict(
    id="gas_ref_recalled", reference="Casella et al., Nucl. Phys. A 706, 203 (2002) [d(p,g)3He gas]; Aliotta et al., Nucl. Phys. A 690, 790 (2001) [d+3He gas]",
    group="LUNA", reaction="d(p,g)3He; d(3He,p)4He", target="gas", method="LUNA underground", E_cm_lowest_keV="2.5; ~5 (recalled)",
    E_cm_range_keV="2.5-22; 5-60 (recalled)", beam_species_purity="analysed", energy_definition="accelerator", vacuum="gas target",
    charge_measurement="calorimeter", detection="BGO; Si", Ue_eV_reported="d+3He gas 219+-7 (as quoted in search results)",
    low_energy_behaviour="no floor", attributed_cause="atomic screening", digitizable_data="paper tables",
    verification="recalled, not re-checked this session"))

with open(sys.argv[1], "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    for r in R:
        w.writerow(r)
print(len(R), "rows")
