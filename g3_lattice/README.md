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
