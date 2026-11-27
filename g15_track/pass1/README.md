# G15 Pass 1: following one impact into the whole foil

This is one coarse conditional model, not a new atomic calculation. It keeps the accepted deuterium profile fixed and follows incoming kinetic energy through named moving particles, nonthermal electrons, deuterium vibration, the full palladium vibration band, transported thermal energy and three exits. It reuses accepted G11/G13 bytes without running or editing their producers.


## What came out

At the nominal0.69W beam, the three conditional closures send0.156,0.304 and0.553W through D optical vibration. Their simultaneous standing optical inventories are only1.25e-14,1.52e-13 and1.14e-12J. Transfer power is not stored energy or phase locking.

| Closure | Peak single-parent optical energy per target D | omega / energy-loss rate | h / threshold if1% is phase locked |
|---|---:|---:|---:|
| fast_acoustic | 0.0168–0.0186eV | 5.80 | 0.0121–0.0133 |
| nominal_mixed | 0.0102–0.0111eV | 36.16 | 0.0459–0.0495 |
| optical_moment_long | 0.00219–0.00231eV | 149.55 | 0.0405–0.0428 |

These are different radius/lifetime/branch scenarios, not a single-parameter lifetime sweep. The calculated fully locked ratios are1.21–4.95; actual phase locking is unknown, and the declared random-phase population is not that drive. In the finite receiver windows6,925 actual nearby parents were drawn;221,179 of221,184 uniform observation times were empty. Sampled conditional peaks coexist with nonzero geometric overlap probabilities,2.6e-10 to8.1e-6, rather than proving the absence of rare events.

Nominal-contact steady temperature rises are11.49K for the even beam and12.23K for hot spots, summing both species. Strong contacts give6.76/7.57K; weak contacts give169.9/170.8K in the explicitly linearized300K boundary law. The latter is a large-excursion contact sensitivity, not an accurate nonlinear radiative equilibrium. Boundary contact knowledge matters to the slow structural/thermal state separately from the optical lifetime.

The joined run conserved energy to9.44e-11 relative, and steady exits to4.42e-9W. Every single parent drained below the stated1e-8 remaining fraction. Kinematic tables consume the time-weighted causal optical spectrum, not the optical transfer-power numbers in this paragraph.

The refined optical velocity consumer adds at most2.0e-10 of the G12 proton residual in these scenarios, including the first-constituent field seen by the second D2+ constituent. Actual arrival-conditioned mean optical stock is4.3e-12 to7.1e-9eV per target D, much smaller than an event-conditioned packet peak. The struck-D source reaches at most5.8e-8 of that residual in the stated transfer/range model. Pair calculations instead expose the needed abundance times contact efficiency: ordinary inward crossings at supplied angstrom-scale cuts can be common within a packet, but none of those crossings is nuclear contact. Their strong depth/phase/energy-sharing dependence is kept in the table, without inventing the missing transmission or vacancy population.

`consumer_spectrum_join.csv` compares36 actual causal optical time integrals with the exact carrier residence; maximum relative defect is1.56e-15. The beam velocity consumer averages oscillator phase and direction with mean kinetic energy E/2 and sums both background source species against the actual traversed bulk profile. Pair relative-mode tails use finite packet energy, joint phase and every supplied separation/approach-energy cut. These are declared conditional state/geometry closures, not measured hot-cascade correlations.

For a scoped change confined to yield consumers, `run.py --yields-only` reuses the unchanged accepted carrier/foil result rather than replaying it. The complete `--all` path calls the same consumer functions.
## Run and open

From the repository root, run:

```
python3 g15_track/pass1/run.py --all
```

The capsule is also self-contained: the same command on its extracted `g15_track/pass1/run.py` needs Python3, NumPy, SciPy and Matplotlib. `bca.py` is the exact accepted stopping helper, with its original docstring. `inputs.json` contains only frozen selected scientific inputs and their source hashes. Nothing is fitted to the observed yields. The source writes only beside itself.

Open `single_impact_channels.png` for the one-impact picture. Its upper panel is incoming/stored energy; its lower panel is cumulative source and cumulative exits, not more stored energy. Logarithmic curves show small coherent contributions without hiding the other carriers. `results.json` and `conservation.csv` state the actual closure defects and stopping tolerance. `inputs_ledger.md` contains the kept/dropped input, all assumptions and citations.

## What the tables mean

- `single_channels.csv`: both electrical-parent species and all three rate scenarios, causal source rows sampled at most0.25fs, then exact logarithmic propagation until remaining energy is below1e-8 of the parent. All energies are eV per electrical parent; D2+ has the same5keV total as D+, not twice that energy.
- `moving_particles.csv`: continuing cutoff-ion trajectory/energy and source-normalized secondary kick/range scenarios. A stopped structural state holds explicit excess energy before its later release.
- `local_space_time.csv`: local radial packet readouts, with their own stated diffusion/radius/depth condition. These nanometre fields are not the millimetre FV grid. `second_constituent_prior.csv` gives only the first constituent's field before the second contributes; axial and transverse launch geometries are distinct.
- `whole_foil.csv` and `foil_*.npz`: three unchanged beam shapes, both species, three rate closures and three boundary scenarios. Temperatures are excess above300K. Species contributions add; free expansion is positive alpha*deltaT, fully constrained hydrostatic stress is negative/compressive. Gradient and clamping limits are not an elastic texture solve.
- `steady_impulse_join.csv`: exact residence-time identity joining each single impulse to the standing optical inventory, separately from its arrival rate and energy loss. `joined_beam_history.npz` propagates an actual100microsecond charge-normalized Cox/no-change history through the same carrier and foil operator, at50ns, starting from the standing state. Six comparable temperature snapshots span that history.
- `carrier_bands.csv`: incoming source energy, visited channel flow, residence and standing stock. A visited flow may sum to more than incident energy when energy passes sequentially through carriers; it is not an extra power source. Incoherent Pd high modes extend to6.8THz. D optical support30–65meV is not a measured homogeneous linewidth.
- `fluctuations.csv`, `packet_marks.csv`, `finite_windows.npz`:64 actual-flux receiver windows per case, selected to contain256 expected geometrically nearby parents. Actual durations, counts, empty records and empty uniform observation draws are retained. Event-conditioned amplitudes and quantiles answer a different question from uniform time averages. A zero finite overlap count does not imply no rare tail: the nonzero geometric probability and observed exposure are separate columns. Signed displacement holes/troughs are not negative site populations.
- `project_ratios.csv`: E/kBT, frequency/energy loss, distinct amplitude loss, neighbour and anharmonic ratios, and actual fractional stiffness modulation divided by the resonant Mathieu threshold. The stiffness harmonic is derived from the same-Kii oscillator's q-squared anharmonicity, with explicit locked fraction0,.01,1; the entire broadband optical stock is not silently called a pump. Static thermal strain changes stiffness, but is not a resonant pump.
- `optical_spectra.csv`, `yield_consumers.csv`, `pair_consumers.csv`, `low_energy_scenarios.csv`: all six available nominal PdD profile settings from0.25 to2.5keV, both branches, both species, three rate closures and all three beam histories. Beam-plus-moving-target adds only the changed velocity distribution relative to the quiet baseline. Struck-D reactions use a distinct source-normalized moving spectrum. Neighbour-pair output is conditional per occupied candidate and shell-to-contact efficiency; it does not invent vacancy abundance or extend the computed potential below its shortest separation. G12 observed/model/residual ordinates and the raw own-reference denominator remain explicit. Per-incident-D predictions are not per charge; the native mixture has1.15 incident nuclei per electrical parent.

## Reconstruct the arrays

All NPZ members are ordinary numeric arrays, readable with `numpy.load(...,allow_pickle=False)`. Foil arrays have axes(x,y,z) with coordinates `x_m,x_m,z_m`; their names include beam, rate, species and snapshot time. Species suffix D1 means D+ electrical parents; D2 means D2+ electrical parents.

`single_response.npz` stores thermal modal histories and the two small separable bases, avoiding a copied dense basis. Reconstruct it as:

```
Vxy = data['thermal_Vxy']; Vz = data['thermal_Vz']
U = np.einsum('ia,jb,kc->ijkabc',Vxy,Vxy,Vz).reshape(2601,2601)
delta_T = (U @ theta) / np.sqrt(data['thermal_volume_m3'])
```

A row `theta` comes from the matching scenario/species `thermal_modes` member; its time is in `thermal_time_s`. Reshape the result to(17,17,9). These fields are Kelvin, not energy. Total thermal energy is CV*sum(deltaT*volume), with CV from the script. A single source occupies one lateral cell; whole-beam fields are spatial superpositions with the same positive operator. Face/mount/radiation integrated exits are separate in the table. Macroscopic beam snapshots retain the same integrated parent and upstream carrier lags but do not resolve the sub-.5ps row-arrival width.

The decisive next physical measurement is the branch-resolved optical occupation decay and phase after an ion/electron impulse, together with a loaded-foil contact/temperature map; the model separates their effects rather than fitting them to the residual. Three largest omitted inputs are stated once in the ledger. This is a private draft, not a commit, publication or permission to restart another line.
