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

`collective/focusing/`. The same lattice model with the 0.8 A wells and the computed pair term, in cells of 32, 256 and 2048 deuterons. Every deuteron starts with the same average energy and the chosen pair starts at rest. The starting pattern is then searched for the one that puts the most energy of relative motion on that pair at one instant. "Gain" is that peak energy divided by the average energy per deuteron.

| Deuterons | Average energy each | Peak on the pair, converging waves (eV) | Peak on the pair, best pattern found (eV) | Gain, best pattern | Share of all the energy that reached the pair |
|---:|---:|---:|---:|---:|---:|
| 32 | 0.6 meV | 0.004 | 0.009 | 15 | 46% |
| 256 | 0.6 meV | 0.010 | 0.021 | 36 | 14% |
| 2048 | 0.6 meV | 0.008 | 0.132 | 221 | 11% |
| 2048 | 4.8 meV | 0.034 | 0.478 | 100 | 4.9% |
| 2048 | 25 meV | 0.157 | 0.83 to 0.98 | 33 to 39 | 1.6 to 1.9% |

- **Focusing is real in the model and grows with the number of deuterons.** At the pumped level (0.6 meV each) the best pattern found takes a pair from 0.6 meV to 0.13 eV, a gain of 221 with 2048 deuterons, and the gain was still rising between 256 and 2048. Simple converging waves do far less (gain 6 to 16) and do not improve with size.
- **It is held for about 20 fs** (18 to 30 fs at half the peak), about a fifth of one vibration period.
- **The gain falls as the energy per deuteron rises:** 221, 100 and 33 to 39 at 0.6, 4.8 and 25 meV with 2048 deuterons. The peaks found top out near 1 eV. In the highest case the most-displaced deuteron reaches 0.79 to 0.80 A, against wells computed to 0.8 A, but the selected runs stayed inside the computed range for their whole length, so they were not cut off by it. What holds the peak down there is not separated yet: the wall toward palladium is up to 3.5 times stiffer than a spring at those displacements and spreads the motion, the search and the run length are finite, and patterns that leave the computed range are not counted. The 1 eV is neither a limit set by the edge nor a limit of the lattice.
- **The highest case is not settled numerically.** The same start gives 0.98 eV with a 0.5 fs step and 0.83 eV with a 0.25 fs step, at different times.
- **With loss switched on** (one value, amplitude decay time 4.1 ps) the highest peak falls from 0.98 to 0.32 eV.
- **As a rate:** a pair at 0.8 to 1.0 eV gives 1e-82 to 1e-88 of the signal units per pair, against 1e-105 for a cold pair. Focusing buys about twenty powers of ten and leaves about eighty.
- **Energy needed, before any inefficiency:** 100 eV on one pair is all the energy of 167,000 deuterons at the pumped level, and 1 keV that of 1.7 million. At the 11 percent share found for 2048 deuterons, if it held at larger sizes (not shown), those become 1.5 million and 15 million.

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

The collective part shows that the pattern of motion, not its size, is what changes the distance between deuterons, and that concentrated motion lasts. With the pair energy now computed, that motion at the energies built first (millielectronvolts per pair) gives a rate about a hundred powers of ten below the added signal. The rate rises steeply with the energy on a pair, and the focusing part shows that shared motion can be gathered onto one pair: a gain of 221 from 2048 deuterons at the pumped level, still rising with size, reaching about 1 eV on a pair in the largest case, with the deuterons then close to the edge of the computed wells. That is still about eighty powers of ten short in rate. What is not known is whether the gain keeps growing in larger regions and past the computed range. The coupling, as built, acts through the same barrier as the ordinary reaction, so by itself it follows the beam energy and does not give a flat signal.

The many-quanta route works in the model only without loss or with infinitely fast one-sided loss, gets steeply harder with the size of the gap, and needs more stored energy in one shared motion than the pumped lattice run holds by a factor of about a million. Neither route, as built, reaches the signal. Focusing is the part that moved the number, and it is not finished. Sites other than the bulk one (the surface, vacancies, titanium) are not built.

## Lead's note

Collective part 1 archive SHA-256 54e98dc8425a7c28b36927b8c071c23a6ef75596974259c677d0816bccab00b9. The lead reran `python g11_cmns/collective/part1/motion.py` (NumPy and SciPy only): all 19 cases reproduce, the largest differences being in the energy-balance checks at 1e-9 eV or less; the committed files are the delivered ones. Coupling archive SHA-256 e0e5b2f47389681ff94d6ccac32a582492605822dfce0eb36d01bee1c1784e35; not rerun (the nuclear code is not installed on the lead's machine). The percentages of the ion's energy and the comparison with the deuterium molecule are the lead's arithmetic. Collective part 2 archive SHA-256 4c2a255493188736f16f178ca1d80fb6174a4931e9e567bea34540b5d90448ac; not rerun (QuTiP and quimb are not installed on the lead's machine). Titanium has no wells yet, so nothing here applies to the titanium foils. Collective part 2b archive SHA-256 19e5928a67c9372d0e04c8d5252f1819acbb94283dd762dc5535dac61b54e0d9; not rerun (QuTiP and quimb are not installed on the lead's machine). `results.json` (64 MB) is committed compressed as `results.json.gz`; unpack it before running `figures.py`. The lead refitted the nine exponents in the table from `transfer.csv` (no loss, supported and uncensored rows, the detuned control rows left out): all nine reproduce to the three decimals given on the page. Part 2b uses a 41.6 meV quantum and the 23.85 MeV gap; part 2 used 37 meV and 34 MeV. Contact archive SHA-256 f1583ca5612247a7b988075995514e596dd20b8c3de17665b6e582823c8e90ff; all file sums check. The lead reran `python contact.py --refine` in a separate copy: every number agrees with the delivered `results.json` to 5 parts in 1e7. The rows at 1 eV and 8 eV are the lead's reading of the per-level contact values in `results.json` (five-point curve, join at 1.29 A) with the delivered conversion. The contact calculation was made before the finer-mesh and all-electron checks of the pair energy were finished. Focusing archive SHA-256 3fb12ab835d463697f53afc7dfbbe643072ea8cb549a4b807c850b098b8df1e4; all file sums check; not rerun (Numba is not installed on the lead's machine). The focusing table is read from `comparison.csv`; the statement that the selected runs were not cut off by the edge of the computed wells is Marvin's, on the lead's question.
