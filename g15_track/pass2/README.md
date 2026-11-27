# G15 Pass 2: does an impact make the deuterons move in step?

We kept Pass 1’s carrier budgets, fixed profiles and three beam shapes. The ion supplies a signed screened kick. The electron pulse supplies a signed force only where the energy density varies or the cage breaks local symmetry. Its strength and sign are inputs to the sensitivity, not measured PdD coefficients. Coherent work is taken from the existing moving-D and optical channels. The largest kicks and later chart exits remain in the moving-D bookkeeping instead of being kept inside a bound oscillator.

## What the kick did

At 2.5 keV, the unit electronic-gradient branch put at most 4.81e−8 eV into the full-footprint retarded axial-clock mode, 2.50e−8 eV into the uniform mode and 2.55e−15 eV into an alternating-sheet projection, across the three closures and molecular geometries. The sheet result uses a continuum average over track registration; a specific alignment of the track with FCC sites was not determined. This calculation estimates the linear electronic part over the full footprint. The local nonlinear motifs include both kicks, the actual finite pair energies and signed work limits. They do not assume that a large common displacement brings two neighbours closer.

In the nonlinear motifs, the largest resolved relative twice-frequency stiffness modulation was 0.00214 of the local threshold, even in the artificially strong, work-limited branches. This comparison uses a trailing fit to actual curvature, with a drift term; windows where the harmonic could not be separated were left unscored. For example, the 2.5 keV D+ nominal bulk case with no electron force and with the unit positive gradient both reach about 0.023 eV of excess coordinate energy, then fall to 0.017% of their sampled peaks by 4.41 ps. Some deliberately very strong cage cases retain displaced states while the applied force persists; that is not free growth of a collective mode.

| Local environment | Resting pair distance (Å) | Closest, electron force absent (Å) | Closest, all settings (Å) | Largest relative kinetic energy (eV) |
|---|---:|---:|---:|---:|
| Bulk PdD | 2.888 | 2.579 | 2.573 | 0.0530 |
| Six-D Pd vacancy | 1.966 | 1.835 | 1.734 | 0.0362 |

The distance and kinetic extrema need not occur together. `pair_consumers.csv` records synchronized transits at each native separation, with its Born–Oppenheimer cut energy and modeled relative kinetic energy identified separately. The finite shell does not determine the probability of reaching nuclear contact.

## Whole-foil readout

The sampled peak-level rates refer to an occupied target drawn from each modeled footprint and are conditional on the location and environment shown. They are not any-site patch probabilities or measured rates. Bulk and vacancy are separate endmembers; the physical fraction of available cages is not supplied. Input profiles were not changed.

At 2.5 keV in the nominal bulk, unit-gradient case, all 228, 257 and 244 sampled centre marks passed 0.01 eV per sampled D for the even, hot-spot and burst beams, respectively. None passed 0.1 eV in those windows. The conservative nearby-parent overlap limits for that support were 5.41e−8, 5.23e−8 and 2.95e−7, respectively; the last uses the maximum imposed burst intensity. The threshold includes the optical population and a local harmonic motion reference, not just the phase-locked component. A zero observed count does not settle the tail frequency.

Counting each modeled **patch** once, rather than drawing one site, changes that comparison. The nominal-bulk unit-gradient maxima are 0.0401 eV (D+) and 0.1058 eV (D₂+). In this fixed-registration model, both species reach 0.01 eV somewhere in the represented patch, while only the molecular parent reaches 0.1 eV. Before overlap, the two rates are therefore 8.61e14 and 1.29e14 patches/s at the 2.5 keV setting; neither reaches 1 eV. `patch_level_rates.csv` keeps both species maxima and the finite-horizon overlap bounds for each beam, separate from the single-site trials.

On the accepted proton reference, the combined D+/D₂+ supported beam-on-moving-target envelope reaches 2.54e−8 of the positive residual; the largest supported struck-D estimate is 4.31e−8 of that residual. Both maxima are at 2.5 keV. The beam number is a magnitude envelope, not the signed central estimate. The signed estimate, matched random comparison, unresolved molecular/struck target exposure and all six beam settings are in `yield_consumers.csv` and `yield_mixture.csv`. For coherent pair trajectories, the innermost *supplied* shells crossed were 2.648 Å in bulk and 1.846 Å in the cage. More extreme between-shell minima in the local table do not determine the contact transmission.

Across the 360 realized force cases the largest final signed-energy residual was 1.50e−8 eV per parent; the peak integrator defect was 1.63e−8 eV. Phase routing leaves the original whole-foil heat and boundary model in place. The deciding missing material input is a transient signed electron force and curvature response in PdD; the vacancy pair additionally needs its cage couplings and contact transmission.

## Replay

From repository base `c037c8c479bc08f6ae3ed01916cf8582cafdcf38`, extract both archives into the same root directory. Pass 1, G3, G10 and G12 inputs remain at their accepted paths, with file hashes in `source_index.json`. Python, NumPy, SciPy and Matplotlib are required.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 g15_track/pass2/run.py --consume
```

`--consume` checks the finished 72-case source fingerprint and uses the saved trajectories without repeating the source run. `--source-only` rebuilds those cases if needed. The main capsule contains the compact result and source code. The evidence archive contains the site-level arrays, trajectory and peak samples, shared electrical-parent marks and detailed mode diagnostics. `inputs_ledger.md` gives the force law, loss convention, pair support and cutoffs.
