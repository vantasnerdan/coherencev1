# g15_track pass 4: travelling fields on the whole foil (lead's summary)

Passes 2 and 3 (`../pass2/LEAD.md`, `../pass3/LEAD.md`) had narrowed to one patch around one impact, with every carrier except sound held as a stock in one place. This pass (Marvin) replaces that with one model of the whole foil in space and time, in which each carrier travels and each loss is the source of the next field where and when it happens. It takes the hot-electron stiffening and pressure of `../../g16_hot_electrons/LEAD.md` as its electron-to-deuterium and electron-to-sound coupling, carried two ways (all in the site stiffness, or all in the neighbour coupling).

`README.md` is Marvin's technical page, `inputs_ledger.md` the inputs, `run.py` the script. Pictures: `maps_*.png` (temperature, strain and deuterium stiffness over the foil at 0.01, 0.5 and 2 s), `first_travel.png` (one impact's sound arriving 5 mm away and coming back from the edge), `carriers_*.png`, `acoustic_*.png`, `medium_*.png`, `yield_comparison.png`.

## What is built

- **The foil.** 21.3 mm square, 250 µm thick, under two seconds of the 1 cm² beam as a random rain of impacts, in three beam shapes (even, hot spots, hot spots with bursts), at six beam energies, with three sets of rates and the two ways of carrying the coupling: 108 cases.
- **The grid.** 17 by 17 cells across (1.25 mm each) and three layers in depth (the top 0.25 µm, where the beam stops, and two below). What happens inside a cell near a track is carried as travelling segments with their own positions, times and phase.
- **The carriers.** Hot and thermal electrons (free flight, then diffusion), the deuterium vibration as a wave with its own frequency, wavelength and group speed, upper palladium vibrations, struck deuterons and recoils, sound in seven bands from 250 kHz to 2.5 THz, the foil's 303 lowest standing modes below that, and heat by conduction to the faces, edges and mount.
- **The couplings, both ways.** Temperature and strain change the deuterium stiffness; hot electrons change it and press on the palladium; the deuterium vibration presses back; the work done one way is taken from the other.
- **Energy.** Followed over all fields to 5 parts in 1e13.

## What came out

1. **Sound crosses the foil.** One impact's low-frequency sound reaches a point 5 mm away after 3.2 µs and returns from the edge, with its sign reversed, at 13.5 µs (group speed 1,574 m/s for the plate wave followed). Waves from different impacts do meet in the standing modes.
2. **How far the other carriers travel.** The deuterium vibration travels 0.3 nm to 0.2 µm from where it is excited before its energy has passed on (0.4 µm at most), depending on the set of rates. Hot electrons and upper palladium vibrations reach about the same distance. All are far below one grid cell.
3. **The largest change to the medium is slow heat.** After 2 s the beam spot is 8.5 K warmer than at the start and still warming. That raises the deuterium vibration frequency by up to 0.56 percent, unevenly across the foil, following the beam's hot spots. Strain from the warm spot shifts it by about 0.01 percent.
4. **The travelling sound barely changes it.** The low-frequency strain is at most 1.8e-8 (burst beam; 1e-12 for the steady beams) and shifts the deuterium frequency by 2.5e-8. The supported sound spectrum ends below twice the deuterium frequency, so sound does not pump it at the lattice layer's resonance.
5. **Hot electrons, near a track.** The electron energy reaches 0.03 eV per site close to a track, where `g16_hot_electrons/` gives a strong stiffening. Averaged over the deuterium wave an impact actually excites, the electronic part is at most 0.065 percent of the stiffness, and changes the wave's group speed by 0.15 m/s.
6. **Meetings of deuterium waves from different impacts.** In the 3.2 µs of rain sampled at full detail (2.5 keV), none occurred; 0.03 were expected. Over the whole foil in two seconds there are 30 to 175 million such meetings in the model. Their phases are followed for single impacts (spread 0.03 radian typical, 0.08 largest); the meetings themselves are not yet built.
7. **Yield.** The built fields add 2.3e-7 to 2.8e-7 of the unexplained signal of `g12_residual/` at 2.5 keV, falling steeply at lower energy (1e-11 at 1 keV, 1e-24 at 0.25 keV). The signal itself is nearly independent of beam energy.

## What the model cannot yet express

1. **Signed sound above 250 kHz.** Only the 303 standing modes carry wave shape, so interference is resolved only at wavelengths of millimetres. From 250 kHz to 2.5 THz, where nearly all of an impact's sound is launched, the bands carry energy and direction but not crests and troughs, so they cannot add up, focus or form tall events.
2. **Between 0.4 µm and 1.25 mm.** Nothing is resolved there. That is the range of the beam's fine structure and of sound wavelengths from megahertz to gigahertz.
3. **The deuterium wave away from the collective motion.** Its speed depends on the neighbour coupling as a function of wavelength, which is still the older `g3_lattice/` shape joined to the new collective value; and its lifetime in loaded PdD is unmeasured (three values carried).

Marvin's own three: the wavelength dependence and site/neighbour split of the new coupling; measured transport and loss of phase in loaded PdD and the phase of the launch; and the step from the closest computed pair distance to contact.

## Lead's note

Core archive SHA-256 9903fa2d65ef711c56c785c250768484905ca0e30aad2da54b3ea2576d067b59, all 64 file hashes in `manifest_core.json` checked, no machine paths. `energy_history.csv` (16 MB) is left out of the repository for size. The evidence archive (378 files) has SHA-256 83732e31ac468610609562f5a0368d63a7324e15dc31448b9c5d9142e16f90ad and is not in the repository.

The lead reran `run.py --consume` with NumPy, SciPy and Matplotlib on a clean copy of both archives alone. This rebuilds the yield tables from the saved fields; they agree with the delivery to 5e-15. The 108 foil runs themselves were not rerun.
