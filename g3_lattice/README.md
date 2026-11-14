# g3_lattice: lattice simulation with loading, arrangement, strain and pump as sliders

A lattice of hydrogen sites (the octahedral interstitials of the fcc metal). Each occupied site carries a three-dimensional displacement in the computed anharmonic well and is coupled to its neighbours by the computed force constants (`g2_input`, `g2_refit`). Zero-point motion enters once, as a stiffening of the well; thermal motion comes from Langevin noise. NumPy only.

    python3 lattice.py lattice_inputs.json --batch cases.json --out results
    python3 lattice.py lattice_inputs.json --batch pump_cases.json --out pump_phase
    python3 map_consumer.py lattice_inputs.json --source map_inputs --out map_results

`map_inputs/` is rebuilt from `native_map_inputs.json` (see the first lines of `delivery.json`'s companion command in the Fleet entry, or load each case's fields into `map_inputs/<case>/fcc_octahedral_grid.npz`). The main batch takes about a minute; a million sites for 32 steps takes 12 seconds.

Changing loading, arrangement, strain, temperature or pump costs nothing beyond a rerun. The wells and couplings are still the full-loading ones: occupancy removes sites and bonds but does not yet change the well of a site with empty neighbours. That law is being computed in `g3_loading`.

## NiH: loading against strain (13,824 sites, random occupancy, 0 K)

| Loading | No strain | 8.5% expansion | 10% expansion |
|---|---|---|---|
| 0.5 | disordered | disordered | static order |
| 0.75 | disordered | disordered | static order |
| 1.0 | disordered | patterned static order | patterned static order |

With the pump on, the unstrained cases at loading 0.5 and 0.75 become driven oscillations; the strained cases keep their static state. The strains here are far beyond the 1 percent range the stiffness was computed over, and the largest displacements leave the range the well was fitted on, so this is a map of the model, not of the material.

## Arrangement matters

At loading 0.75 and 8.5% expansion:

| Arrangement of the vacancies | Result | RMS static displacement |
|---|---|---:|
| Random | no ordered region | 0.000004 A |
| Clustered | 62 active sites, in singles, pairs and one triple | 0.011 A |
| Ordered superstructure | 12 periodic sheets of 576 sites each | 0.134 A |

The ordered case has a preferred wavelength of 0.53 nm in its displacement pattern. A strain gradient of 0.005 per nm also gives a patterned state.

## Pump and spectrum

A cold pumped lattice oscillates with RMS amplitude 0.06 A. Its spectrum peaks at 22.8 THz, half the stiffness-modulation frequency, with a weak line at three halves of it (1.7e-5 of the main line). At 300 K the same pump is not enough unless its frequency is raised by 10 percent. The line does not sharpen between 13,824 and 32,768 sites within the frequency resolution of these runs, so no size-dependent coherence is shown. No stable family of sidebands was resolved.

## On the Ni-on-Pd pattern (`g3_skew`)

The per-site stiffness map from the mismatch pattern was applied as a change relative to bulk. No static order forms. Under pump, the response sits almost entirely in the Pd layer at the interface (over 96 percent of the mode's weight) and appears as one hot site per repeat of the pattern: 25 sites for bare Pd, 81 for Pd at the PdD spacing. The response inherits the pattern's period. It stays as isolated sites; the largest connected cluster is 13 sites. This is conditional on the first-estimate map, which extrapolates the strain response and contains no hydrogen in its geometry.

## What this does and does not show

It shows that the model forms spatial patterns, that vacancy arrangement changes the outcome at fixed loading, that a pump produces a half-frequency oscillation, and that an interface pattern localises the pumped response. It does not show a large phase-coherent domain, and the PdD inputs remain provisional.

Figures: `results/*.png` (occupied sites and displacement in one layer), `map_results/*.png` (stiffness contrast and pumped response at the interface and in depth). `delivery.json` holds all 33 cases with spectra.

## The decay rate is not measured (`loss_table.py`, `loss_table.json`)

The pump threshold is proportional to how fast the vibration decays. No measured decay rate for a single hydrogen vibration was found for NiH, PdH or PdD: the published neutron data give the width of a whole band of vibrations, or the instrument resolution, and neither is a lifetime. The threshold is therefore unknown. Two labelled scenarios show how far it could range, without implying the real value lies between them: the exploratory loss used in these runs gives a stiffness-modulation threshold of 1.7 percent; a published calculated value for PdH at 295 K (Paulatto et al., arXiv:1411.5628; calculated, not measured) gives about 23 percent.

## Does the pumped state stay in step over distance? (`coherence/`)

    cd coherence && python3 replot_phase.py     (redraws the figures from saved data; no simulation)

Cold (0 K), full loading, pump at 8 percent stiffness modulation. The oscillation is not uniform: neighbouring sheets of hydrogen move in opposite directions, in stripes with a period of 0.53 nm. The measure below is how far that striped pattern keeps the same phase.

| Sites | Phase length at 41 ps |
|---:|---:|
| 13,824 | 1.44 nm |
| 110,592 | 2.01 nm |
| 1,000,000 | 2.16 nm |

- **It does not grow with the size of the system.** A million sites give the same length as a hundred thousand. The lattice breaks into patches a few nanometres across, each striped along a different diagonal, with walls between them (`coherence/cold_x1_h008_L100_t40_phase.png`).
- **It is still growing in time.** Between the windows 21-41 ps and 41-83 ps the fraction of the lattice in the leading pattern rises from 0.16 to 0.29: the patches are coarsening. The runs do not resolve how fast, or where it stops.
- **At 300 K there is no order beyond nearest neighbours.** A stronger or retuned pump raises the amplitude of the motion but not its range, which stays comparable to an unpumped lattice.

The length quoted is an integral of the phase correlation, not a fitted decay length. All of this is for NiH with the assumed decay rate, which is not measured.

## Rerun on the corrected wells (`rerun/`)

The wells this directory first used were recomputed in `g3_loading/` (PdD on-site stiffness 0.59 instead of 1.09 eV/A^2, coupling 0.39 instead of 0.064; NiH coupling 0.68 instead of 0.23). `rerun/` holds the results on the corrected inputs. Where they differ from the sections above, these replace them.

### Part 1: Pd oscillator and the vacancy force law (`rerun/part1/`)

Reproduce: `cd rerun/part1/pd_x1_quantum && python3 pd_x1_refit.py`, and `fit_vacancy_bias.py` as given in `rerun/part1/updated/part1.md`. Both were rerun here and reproduce the delivered numbers.

| | PdH | PdD |
|---|---:|---:|
| First excitation, one atom in its well (meV) | 74.0 | 48.9 |
| Collective mode, all atoms in step (meV) | 62.3 | 39.8 |
| Measured neutron peak (meV) | 55.8 | 37 |
| Static ordering strain at 0 K / 300 K (percent) | 7.5 / 8.8 | 5.3 / 7.1 |
| Collective mode frequency (THz) | 15.1 | 9.6 |
| Stiffness-modulation (pump) frequency (THz) | 30.1 | 19.3 |
| Modulation threshold at resonance (percent) | 1.0 | 1.6 |

H to D ratio: 1.51 for the single-atom excitation, 1.56 for the collective mode (measured peaks: 1.51). The measured peaks are averages over the vibration band, not this one mode. The ordering strains extrapolate from calculations at plus and minus 1 percent strain. The pump thresholds assume a loss rate of 0.48 per picosecond, which is not measured.

Vacancy force law: each vacant nearest-neighbour site pulls a D toward itself with 0.0707 eV/A along the line joining them; the third-shell term is negligible (-0.0001); the second-shell term cannot be determined from cells of this size and is set to zero. Fitted to 80 force vectors from the three partial-loading cells, RMS error 0.0041 eV/A; predicting each cell from the other two, 0.0030, 0.0086 and 0.0043 eV/A. Using the law in place of the computed forces changes the static displacement pattern in those cells by 2.4 to 3.6 percent.

### Part 2: lattice, loading, interface and NiH (`rerun/part2/`)

`updated_batch.py` ran 30 PdD cases, `ni_epoch.py` the NiH case, `updated_figures.py` redraws the three figures from the saved data (`updated/figures/`). Full write-up: `rerun/part2/updated/part2.md`. Temperatures of every case: `updated/temperatures.json`.

**NiH, pumped, 0 K, 110,592 sites.** With the corrected coupling the striped patches reported above do not appear. The hydrogens move in step across the whole simulated box: 99.8 percent of the pumped motion is in the uniform (all-in-phase) pattern. The in-step length, 3.85 nm by the same measure as before, is set by the size of the box, so this run puts no upper limit on the coherent region. It was not repeated at a larger size or at 300 K. NiH oscillator on the corrected inputs: single-atom first excitation 100.6 meV, collective mode 85.6 meV (measured peak 89); static ordering strain 8.1 percent at 0 K, 8.2 at 300 K; collective mode 20.7 THz, pump at 41.4 THz, threshold about 1.9 percent stiffness modulation for an assumed loss of 1.26 per picosecond. This run used the earlier NiH quartics (18.3, 1.13); with the recomputed 100 quartic (16.2) the single-atom excitation moves by -1.3 percent.

**PdD, pumped (8 percent stiffness modulation), 0 K.** Static displacement with no pump; pumped oscillation amplitude; in-step length; share of the pumped motion in the single strongest pattern:

| Loading, arrangement | Static (A) | Pumped (A) | In-step length (nm) | Strongest-pattern share |
|---|---:|---:|---:|---:|
| 1.0 | 0.000 | 0.038 | 1.38 | 0.12 |
| 0.875 random | 0.072 | 0.025 | 1.38 | 0.12 |
| 0.75 random | 0.093 | 0.028 | 1.29 | 0.10 |
| 0.6 random | 0.105 | 0.032 | 0.88 | 0.03 |
| 0.875 ordered | 0.000 | 0.043 | 0.82 | 0.03 |
| 0.875 clustered | 0.067 | 0.033 | 1.75 | 0.24 |

PdD responds to the pump at every loading but stays patchy, with an in-step length near 1 to 2 nm. Vacancies displace their neighbours statically (vacancy-driven displacement, not ordering); an ordered vacancy superlattice cancels that displacement by symmetry. Clustered vacancies give the longest in-step length, ordered vacancies the largest amplitude but the shortest. At 300 K the pump adds nothing measurable (0.00915 A without, 0.00919 with). Stretching by 5 percent or raising the pump frequency by 10 percent increases the amplitude and shortens the in-step length. Loadings below 0.75 use the loading law outside the range it was fitted on.

**Second-shell vacancy force.** Setting it to plus or minus half the first-shell value changes which pattern leads (the uniform share goes from 0.005 to 0.001 or 0.155), so the unknown coefficient affects the texture, not only the size.

**Ni on Pd.** With the strain law and loading law applied to the mismatch map, the pumped amplitude at the PdD-spacing interface rises about tenfold (0.0013 to 0.0146 A) and concentrates in the soft regions (2.9 times enrichment), but the pumped pattern does not copy the period of the mismatch pattern (correlation 0.25). At the bare-Pd spacing the effect is negligible. Every site in this map is outside the strain range that was computed.
