# g11_cmns: deuterons moving together, and the lattice-to-nucleus coupling (lead's summary)

`g12_residual` showed that the gap between model and measurement is an added signal that barely depends on beam energy. This node builds the one thing no earlier node built: many deuterons moving together, away from equilibrium. `survey/` is the literature table. `collective/part1/` is the motion of the deuteron sublattice. `coupling/` is the strength of the coupling between lattice motion and the nucleus. `collective/part2/` is a small many-body model with loss.

## Collective part 1: what collective motion does to the distance between neighbours

A cell of 32 deuterons in the PdD wells and neighbour couplings of `g10_quantum_site`, followed for 2 picoseconds. The same total energy (4.8 meV, set by the displacement in the pumped lattice run of `g3_lattice`) is put in with three different patterns.

| Pattern, same energy | Closest two neighbours get (A) | Peak energy of their relative motion | Energy held on one site |
|---|---:|---:|---|
| All deuterons in step | 2.888 (unchanged) | none | no |
| A wave across the cell | 2.860 | 0.6 meV | no |
| Concentrated on one site | 2.804 | 2.7 meV | a quarter or more for 0.41 ps |
| Concentrated, at twice the displacement | 2.735 | 14 meV | a quarter or more for 1.07 ps |

- **Motion in step does nothing to a pair.** Equal displacements cancel in the distance between two deuterons. What brings a pair closer is the pattern, not the amplitude. This is why a larger common amplitude was never going to show up in the reaction.
- **Concentration holds.** Energy placed on one site stays there for about a picosecond in the anharmonic lattice. That is long compared with the 0.15 to 0.5 fs arrival transient in `g9_edge`.
- **The closest approach and the largest relative speed do not coincide.** In every cold case the closest approach is at a turning point, where the relative speed is zero. They cannot be combined into one event.
- **The pumped state itself is not represented.** Its pattern has a wavelength of 7.5 nm, and the cell is 0.8 nm across.

**Heat, a beam ion, and a stopped deuteron.**
- At 300 K (classical motion) two neighbours come as close as 1.98 A. A stopped beam deuteron placed on an empty site comes within 1.93 A of a neighbour. Both are outside the range where the wells were computed (0.24 A from the site), so these two numbers rest on an extrapolated well. Only about half of the motion at 300 K stays inside the computed range.
- Where a 0.25 keV (centre of mass) ion's energy goes in PdD, with the fixed profile: 59 percent to electrons, 18 percent to recoil of palladium atoms, 0.014 percent directly to recoil of lattice deuterons, and 23 percent not followed below the transport cutoff. How much of the electron and palladium share ends up in deuteron motion is not known and was carried as a stated sensitivity.

**The rate for a pair in the lattice, now computed** (`collective/contact/`). With the pair energy from `g10_quantum_site/pdd-nearest-pair/` joined at short range to the earlier fast-ion pair calculation, the relative motion of two neighbouring deuterons is solved as a quantum problem, and the chance of finding them at the same point is turned into a rate in the units of the `g12_residual` signal (PdD at 0.25 keV, where the signal is 0.05).

| State of one pair | Added signal from that one pair, as a power of ten |
|---|---:|
| Cold (77 K) | -105 |
| 300 K | -100 |
| Concentrated motion of part 1, averaged over 2 ps | -105 |
| The pair held at 1 eV of relative energy | -84 |
| The pair held at 8 eV of relative energy | -58 |

- **Per pair it is about a hundred powers of ten below the signal.** No count of pairs in a foil closes that. The concentrated motion of part 1 raises the contact chance 2.7 to 2.9 times over a cold pair: it carries about 1 meV of relative energy on the pair, and the rate needs far more than that.
- **The rate climbs steeply with the energy on the pair,** about six powers of ten per additional 0.4 eV near 1 eV and about one power per 0.6 eV near 8 eV. On this curve a lattice source would have to put hundreds of eV or more onto single pairs. That is the same requirement `g9_edge` reached from the stored-energy side (about 1.7 million pumped deuterons' motion arriving on one pair for 1 keV).
- **How the short range is joined matters little when all five computed pair points are used** (7 percent at 300 K between joining at 1.29 A and at 1.6 A). If the two flagged short points are left out and the curve is continued by its last slope, the result moves by up to a factor of about a thousand. Against a hundred powers of ten, neither changes the reading.

This is one pair on one line, with the neighbours held still, relative motion with no angular momentum, and a smooth assumed form below 0.1 A. The thermal rows are a quantum treatment of that one line, not the classical lattice runs above.

## Focusing: how much shared energy can arrive on one pair

`collective/focusing/`. The same lattice model with the 0.8 A wells and the computed pair term, in cells of 32 to 131,072 deuterons. Every deuteron starts with the same average energy and the chosen pair starts at rest. The starting pattern is then searched for the one that puts the most energy of relative motion on that pair at one instant. "Gain" is that peak energy divided by the average energy per deuteron.

| Deuterons | Average energy each | Peak on the pair, converging waves (eV) | Peak on the pair, best pattern found (eV) | Gain, best pattern | Share of all the energy that reached the pair |
|---:|---:|---:|---:|---:|---:|
| 32 | 0.6 meV | 0.004 | 0.009 | 15 | 46% |
| 256 | 0.6 meV | 0.010 | 0.021 | 36 | 14% |
| 2048 | 0.6 meV | 0.008 | 0.132 | 221 | 11% |
| 16,384 | 0.6 meV | 0.007 | 0.125 | 208 | 1.3% |
| 131,072 | 0.6 meV | 0.006 | 0.351 | 585 | 0.45% |
| 2048 | 4.8 meV | 0.034 | 0.478 | 100 | 4.9% |
| 2048 | 25 meV | 0.157 | 0.83 (pair moving apart at the peak) | 33 | 1.6% |

- **Focusing is real in the model and the peak keeps growing with the number of deuterons.** At the pumped level (0.6 meV each) the best pattern found puts 0.13 eV on a pair from 2048 deuterons, 0.12 eV from 16,384 and 0.35 eV from 131,072, a gain of 585. The step from 2048 to 16,384 is flat and the next step rises again, so no growth law can be read off. Simple converging waves do far less (gain 6 to 16) and do not improve with size.
- **The share of the energy that arrives keeps falling:** 11 percent, 1.3 percent and 0.45 percent at the three largest sizes. The peak grows because there is more energy in the box, not because more of it is gathered.
- **The search is small, so these are lower figures for the model.** At each size 11 starting patterns are tried, all built from the small-motion (harmonic) picture and then adjusted once. The run is 2.4 ps long, in which the fastest lattice wave covers about 93 A, less than the distance to the corners of the largest box (113 A). In that box some of the deuterons cannot contribute within the run.
- **It is held for about 20 fs** (18 to 30 fs at half the peak), about a fifth of one vibration period.
- **The gain falls as the energy per deuteron rises:** 221, 100 and 33 at 0.6, 4.8 and 25 meV with 2048 deuterons.
- **The highest case, redone with the finer step throughout:** all 11 patterns rerun at 0.25 fs give 0.83 eV as the best, and at that instant the pair is moving apart. The 0.98 eV found earlier with the coarser step, with the pair closing, is not reproduced. The closest the pair comes anywhere in that search is 1.55 A.
- **The edge of the computed wells does not set these peaks.** The searches at 4.8 and 25 meV were repeated with the stop at 0.8 A removed and the wells continued beyond it in two different ways, one soft and one stiff. Both give the same best patterns and the same peaks (0.478 and 0.826 eV), and those runs never leave the computed range. Three of the eleven patterns at 25 meV do leave it, reaching 0.84 A on a deuteron that is not in the chosen pair, and none of them beats the best one. What holds the peak near 1 eV at the higher energies is therefore inside the model: the motion is no longer close to harmonic there, and patterns designed from the harmonic picture lose their timing. A search that works on the full motion has not been made.
- **With loss switched on** (one value, amplitude decay time 4.1 ps) the coarse-step peak of 0.98 eV fell to 0.32 eV. The larger boxes and the fine-step case were not run with loss.
- **As a rate:** a pair at 0.8 to 1.0 eV gives 1e-82 to 1e-88 of the signal units per pair, against 1e-105 for a cold pair. Focusing buys about twenty powers of ten and leaves about eighty. The new peaks were not converted to rates.
- **Energy needed, before any inefficiency:** 100 eV on one pair is all the energy of 167,000 deuterons at the pumped level, and 1 keV that of 1.7 million. The largest box run holds 79 eV in all. At the 0.45 percent share found there, 100 eV on a pair would take the energy of 37 million deuterons, and the share is still falling with size.

### Patterns adjusted against the full motion

The starting pattern is now improved step by step against the actual, non-harmonic motion, with the same total energy, the chosen pair at rest and the same limit on any one deuteron's starting speed. Three things are read separately: the peak energy of relative motion on the pair, the same counted only while the pair is closing, and the closest the pair comes.

| Deuterons, energy each | Peak on the pair before (eV) | Peak now (eV) | Peak while closing, before (eV) | Peak while closing, now (eV) | Closest approach before (A) | Closest approach now (A) |
|---|---:|---:|---:|---:|---:|---:|
| 2048, 4.8 meV | 0.478 | 0.569 | 0.478 | 0.569 | 1.79 | 1.70 |
| 2048, 25 meV | 0.826 | 1.027 | 0.809 | 0.969 | 1.55 | 1.46 |
| 16,384, 0.6 meV | 0.125 | 0.295 | 0.116 | 0.295 | 2.22 | 1.92 |

- **The gain from a better search is modest at 2048 deuterons (a fifth) and large at 16,384 (2.4 times),** and the 16,384 searches were still improving when they were stopped. The neighbour distance at rest is 2.89 A.
- **About 1 eV while closing takes a pair to 1.46 A,** which is where the computed pair energy says it should turn round. The pair closes, stops and rebounds within about 30 fs.
- **The 25 meV runs now reach the edge of the computed wells** (0.80 A) and leave it after the peaks counted here. Going higher than about 1 eV on a pair in this model needs wells and pair energies beyond what is computed.

**With loss designed in** (amplitude decay time 4.1 ps, the pattern optimised with the loss acting):

| Deuterons, energy each | Old pattern, with loss (eV) | Pattern above, with loss (eV) | Pattern optimised with loss (eV) | Without loss (eV) |
|---|---:|---:|---:|---:|
| 2048, 4.8 meV | 0.412 | 0.480 | 0.482 | 0.569 |
| 2048, 25 meV | 0.331 | 0.561 | 0.777 | 1.027 |
| 16,384, 0.6 meV | 0.082 | 0.217 | 0.268 | 0.295 |

- **Loss at this value costs 10 to 25 percent of the peak when the pattern is designed for it,** not the factor of three seen earlier with a pattern designed without it. While closing, the lossy peaks are 0.48, 0.75 and 0.24 eV.

### The closing pulses as a reaction chance

`collective/focusing/pulse/`. Each closing-and-rebound of the chosen pair (about 30 fs) is turned into a chance of one proton-branch reaction, by driving the quantum model of the pair's separation (the one used for the contact rate above) along the separation the lattice run gives. Four versions of the pair energy inside the computed range and two choices of how many levels to keep are carried, and the spread across them is given, as powers of ten.

| Deuterons, energy each | Energy while closing (eV) | Closest approach (A) | Chance of a proton per pair per pulse | In signal units, if every beam deuteron met one such pulse |
|---|---:|---:|---:|---:|
| 2048, 25 meV | 0.97 | 1.46 | 1e-95 to 1e-91 | 1e-80 to 1e-76 |
| 2048, 25 meV, loss designed in | 0.75 | 1.58 | 1e-101 to 1e-97 | 1e-87 to 1e-82 |
| 2048, 4.8 meV | 0.57 | 1.71 | 1e-104 to 1e-100 | 1e-89 to 1e-85 |
| 2048, 4.8 meV, loss designed in | 0.48 | 1.77 | 1e-106 to 1e-102 | 1e-91 to 1e-87 |
| 16,384, 0.6 meV | 0.29 | 1.92 | 1e-108 to 1e-106 | 1e-94 to 1e-91 |
| 16,384, 0.6 meV, loss designed in | 0.24 | 2.06 | 1e-110 to 1e-106 | 1e-95 to 1e-91 |

- **The best pulse built is about 75 powers of ten short,** under the generous assumption of one such pulse for every beam deuteron. The signal is about 0.05 in these units. Nothing in the model says how often such a pulse would form.
- **These are examples, not settled numbers.** The quantum model's average separation does not follow the lattice path (it is off by 0.13 to 0.39 A), and keeping more levels moves the result by up to 2.6 powers of ten. That uncertainty is small against the 75.
- They agree with the earlier steady-level estimate (1e-82 to 1e-88 for a pair held at 0.8 to 1.0 eV) to within what the two different assumptions allow.

### Larger regions, searched for closing energy

| Deuterons at 0.6 meV each | Energy while closing (eV) | Closest approach (A) | With loss (eV) |
|---|---:|---:|---:|
| 16,384 | 0.355 (was 0.295) | 1.86 | 0.270 (pattern optimised with loss) |
| 131,072 | 0.365 (was 0.351) | 1.88 | 0.248 (same pattern, loss switched on) |

- **Eight times more deuterons gave almost the same pulse.** The 16,384 case had 25 improvement steps and was stopped when further steps were adding little; the 131,072 case had one. So the two are not searched to the same depth and no size law follows, but there is no sign here of the peak growing in proportion to the energy available (9.8 eV against 79 eV in the box).
- **Every best pattern sits on the limit set for the start.** The searches require the starting energy to be spread over at least a quarter of the deuterons. That limit is a choice made to keep the start "shared", not a property of the material, and the peaks depend on it.

Left out: the forces off the computed lines, moving palladium and electrons apart from the one loss value, and how such a starting pattern would be produced.

## Coupling: how strongly lattice motion couples to the nucleus

The coupling is Hagelstein's boost term (doi 10.1088/1361-6455/acf3be, Eq. 27).

- **From the head-on (s-wave) entrance it is exactly zero at first order** for any 0+ state of helium-4. It needs the next entrance state up (p-wave), or two steps.
- **Inside the nucleus it is eV-sized.** Between internal states of the four-nucleon system, one deuteron moving in the 37 meV PdD mode gives a mixing of 6 eV; the pumped-lattice displacement gives 1.6 eV. The value falls 20 percent when the basis is enlarged from the smaller to the larger size tried, so it is not converged.
- **For two deuterons arriving at beam energy it is tiny, because they still have to tunnel.** With the repulsion included, the same coupling into the p-wave entrance is about 1e-18 eV at 250 eV and 1e-8 eV at 2.5 keV. It carries the same steep fall with energy as the ordinary reaction.
- **The two-step route** is about 4e-7 eV in size, but it connects states 34 MeV apart, and one mode can supply at most two quanta (0.074 eV) at a time. The page also shows that taking the first-order term twice is not the complete second-order operator; the complete one is not known.

## Collective part 2: one shared vibration and a nuclear excitation, with loss

A small quantum model: one 37 meV vibration shared by one, two or four cells, coupled with the 6 eV example value above to a nuclear excitation 34 MeV away, with loss to the surroundings. Run with QuTiP and quimb, which agree.

- Starting from four vibration quanta, the energy drains to the surroundings and the nuclear excitation is not populated.
- Starting from the nuclear excitation, it stays put unless nuclear decay is switched on, and then it leaves by the decay channels put in, not into the vibration.
- Going from one to four cells at fixed occupation changes nothing; going from zero to four quanta raises the (very small) mixing in proportion to 2n+1.

This model is far from the regime the collective-transfer proposals describe. It holds at most four quanta and four cells, and the gap is about a billion quanta wide. Part 2b below is the scaled family.

## Collective part 2b: exchanging many quanta at once, and what loss does to it

The same kind of model, now with the gap set to 3, 5, 9, 15 and 21 vibration quanta, 2 to 8 cells, and tens of thousands of quanta in the vibration. The coupling is a chosen model value, not the nuclear one. "Transfer" means the whole gap is exchanged in one step, in either direction.

- **Without loss, many quanta do move at once.** With a gap of 15 quanta, two cells and 30,000 quanta in the vibration, the chance of finding the gap exchanged reaches 0.17 within 210 ps (taking the vibration quantum as 41.6 meV). With a gap of 21 quanta it reaches 0.015.
- **Realistic loss removes it.** With loss on the intermediate states at 100 times the vibration frequency, placed on both sides, nothing arrives and everything is lost. Placed only on the low-energy side, which is the condition the published lossy models rely on, 2e-8 arrives and 95 percent is lost. Only in the limit of infinitely fast loss does transfer come back: 0.003 in the same 210 ps, and 0.32 if the model is followed for 7.5 ns.
- **Up and down are the same.** The chance of arriving is identical in the two directions. What is lost on the way differs.
- **More cells sharing the vibration makes it slower, not faster,** when the coupling of each cell is divided as the model divides it. Fitted from the runs without loss:

| Gap, in quanta | Rate goes as (number of cells) to the power | as (quanta in the vibration) to the power | as (coupling) to the power |
|---|---:|---:|---:|
| 3 | -1.1 | 1.1 | 2.6 |
| 5 | -2.0 | 2.2 | 4.2 |
| 9 | -3.9 | 3.6 | 7.3 |

Each power grows in step with the gap. For gaps of 15 and above with four or more cells the runs did not hold a clean branch, and no law is given there. A separate exact calculation at weak coupling shows that one-sided loss in the infinite limit removes the penalty for adding cells.

**Carried to the real gap.** The real gap (23.85 MeV, two deuterons to helium-4) is 573 million quanta of this vibration. The pumped lattice run holds 0.6 meV per deuteron, so the 32,768-deuteron box holds 472 quanta and a million deuterons hold 14,000. Neither holds enough energy to supply one gap. About 1e11 deuterons moving as one would. At that size the weak-coupling law gives an amplitude with an exponent of minus several billion before the unknown coupling is put in. That is a formula carried far outside where it was fitted, and the page says so; it is not a predicted number.

## What this says so far

The collective part shows that the pattern of motion, not its size, is what changes the distance between deuterons, and that concentrated motion lasts. With the pair energy now computed, that motion at the energies built first (millielectronvolts per pair) gives a rate about a hundred powers of ten below the added signal. The rate rises steeply with the energy on a pair, and the focusing part shows that shared motion can be gathered onto one pair: at the pumped level the peak on a pair reaches 0.35 eV from 131,072 deuterons and is still growing with size, while the share of the energy that arrives falls below one percent. At higher energy per deuteron the peaks stay near 1 eV, and that is not caused by the edge of the computed wells. A search on the full motion adds a fifth at 2048 deuterons and more than doubles the peak at 16,384, and designing the pattern with the loss in it keeps 75 to 90 percent of the peak. About 1 eV on a pair, which brings it from 2.89 to 1.46 A, gives, pulse by pulse, a reaction chance about 75 powers of ten below the signal even if every beam deuteron met one such pulse. Past 1 eV the model runs out of computed wells. At the pumped level, going from 16,384 to 131,072 deuterons left the pulse at about 0.36 eV. Shared motion of the size the lattice run holds, gathered as well as these searches can gather it, does not come near the signal through the ordinary barrier. What is not built is any way for the collective state to act other than by pushing two deuterons together: a change in the barrier or in the nuclear step itself. How such a starting pattern would arise is also not built. The coupling, as built, acts through the same barrier as the ordinary reaction, so by itself it follows the beam energy and does not give a flat signal.

The many-quanta route works in the model only without loss or with infinitely fast one-sided loss, gets steeply harder with the size of the gap, and needs more stored energy in one shared motion than the pumped lattice run holds by a factor of about a million. Neither route, as built, reaches the signal. Focusing is the part that moved the number, and it is not finished. Sites other than the bulk one (the surface, vacancies, titanium) are not built.

## Lead's note

Collective part 1 archive SHA-256 54e98dc8425a7c28b36927b8c071c23a6ef75596974259c677d0816bccab00b9. The lead reran `python g11_cmns/collective/part1/motion.py` (NumPy and SciPy only): all 19 cases reproduce, the largest differences being in the energy-balance checks at 1e-9 eV or less; the committed files are the delivered ones. Coupling archive SHA-256 e0e5b2f47389681ff94d6ccac32a582492605822dfce0eb36d01bee1c1784e35; not rerun (the nuclear code is not installed on the lead's machine). The percentages of the ion's energy and the comparison with the deuterium molecule are the lead's arithmetic. Collective part 2 archive SHA-256 4c2a255493188736f16f178ca1d80fb6174a4931e9e567bea34540b5d90448ac; not rerun (QuTiP and quimb are not installed on the lead's machine). Titanium has no wells yet, so nothing here applies to the titanium foils. Collective part 2b archive SHA-256 19e5928a67c9372d0e04c8d5252f1819acbb94283dd762dc5535dac61b54e0d9; not rerun (QuTiP and quimb are not installed on the lead's machine). `results.json` (64 MB) is committed compressed as `results.json.gz`; unpack it before running `figures.py`. The lead refitted the nine exponents in the table from `transfer.csv` (no loss, supported and uncensored rows, the detuned control rows left out): all nine reproduce to the three decimals given on the page. Part 2b uses a 41.6 meV quantum and the 23.85 MeV gap; part 2 used 37 meV and 34 MeV. Contact archive SHA-256 f1583ca5612247a7b988075995514e596dd20b8c3de17665b6e582823c8e90ff; all file sums check. The lead reran `python contact.py --refine` in a separate copy: every number agrees with the delivered `results.json` to 5 parts in 1e7. The rows at 1 eV and 8 eV are the lead's reading of the per-level contact values in `results.json` (five-point curve, join at 1.29 A) with the delivered conversion. The contact calculation was made before the checks of the pair energy were finished; they move the pair term by 1 percent, which does not change the reading. Focusing archive SHA-256 3fb12ab835d463697f53afc7dfbbe643072ea8cb549a4b807c850b098b8df1e4; all file sums check; not rerun (Numba is not installed on the lead's machine). The focusing table is read from `comparison.csv`; the statement that the selected runs were not cut off by the edge of the computed wells is Marvin's, on the lead's question. Focusing continuation archive SHA-256 e647f78288dab7eab4c51b4f5bfedc56d288297b58e81311394bad849be6af83; all 37 manifest hashes check; not rerun (Numba). The lead recomputed every case's best peak, its time, the largest displacement and the closest pair distance from the delivered traces (`continuation_traces.npz`); all agree with `continuation_results.json` exactly. The archive also carried copies of earlier inputs, which are identical to the committed ones and were not added again. The 37 million figure and the reading that harmonic design loses its timing are the lead's. Full-motion search archive SHA-256 5e07c466a4d4c59750f4c6575677021c6626eb94da0496b5a923da1b995d7add; all 34 manifest hashes check; not rerun (Numba). The lead recomputed the peaks, the peaks while closing, the closest approaches and the lossy peaks from the delivered traces; they agree with the delivered tables. Marvin marks the 16,384 results and some of the 2048 ones as best found, not as settled optima. Pulse and size archive SHA-256 32f184fc936038bc5bd4060256a4d08df8af4950daff16964c619ebcfa60caa1; all 55 manifest hashes check. The lead reran `pulse/replay.py` (NumPy and SciPy only) in a separate copy: its built-in comparison of recomputed levels with the stored ones stops on the lead's machine (energy difference 1e-9 eV against a limit of 5e-10), and with that limit loosened to 5e-9 every result agrees with the delivered `pulse/results.json` to 4 parts in 1e5. The closing energies and closest approaches of the 16,384 and 131,072 runs were recomputed from the delivered traces and agree. The pulse table gives the lowest and highest of the four versions and two level choices, rounded to whole powers of ten. The larger-size searches were not rerun (Numba). Site-pair patch archive SHA-256 9fc9fb3ea608912e818214078215fb0a623a8c1ca0794abcdff42b285f525d1f: `contact.py` and `pulse/replay.py` now take a description of any pair of sites (`contact/bulk_pair.json` is the bulk one). The lead reran both on the bulk pair in a clean copy: the contact results agree with the committed ones to 1 part in 1e6 and the pulse results to 6 parts in 1e5. The example output paths in the delivered pages were replaced by placeholders.
