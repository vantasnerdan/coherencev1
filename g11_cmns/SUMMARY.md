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

**No reaction rate yet.** Nothing in the repository gives the energy of two lattice deuterons between 2.4 A apart (the edge of the computed wells) and 1.6 A (the outermost point of the earlier pair calculations), and those earlier points were computed for a fast ion on one ray, not for two bound deuterons. The page gives the conversion from "how much the two deuterons overlap at contact" to a rate in the units of the `g12_residual` signal, and leaves the overlap as the unknown. For scale: the closest approaches found so far are about 2 A, and the two nuclei in a deuterium molecule sit 0.74 A apart. The pair cut now being computed by giuliano reaches about 1.3 to 1.6 A and joins the two ranges along one path.

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

The collective part shows that the pattern of motion, not its size, is what changes the distance between deuterons, and that concentrated motion lasts. Whether that gives a reaction rate of the size of the added signal cannot be said until the pair energy between 2.4 and 1.3 A is computed. The coupling, as built, acts through the same barrier as the ordinary reaction, so by itself it follows the beam energy and does not give a flat signal.

The many-quanta route works in the model only without loss or with infinitely fast one-sided loss, gets steeply harder with the size of the gap, and needs more stored energy in one shared motion than the pumped lattice run holds by a factor of about a million. The open route for the flat signal is therefore the first one: concentrated motion bringing pairs closer, which waits on the pair energy.

## Lead's note

Collective part 1 archive SHA-256 54e98dc8425a7c28b36927b8c071c23a6ef75596974259c677d0816bccab00b9. The lead reran `python g11_cmns/collective/part1/motion.py` (NumPy and SciPy only): all 19 cases reproduce, the largest differences being in the energy-balance checks at 1e-9 eV or less; the committed files are the delivered ones. Coupling archive SHA-256 e0e5b2f47389681ff94d6ccac32a582492605822dfce0eb36d01bee1c1784e35; not rerun (the nuclear code is not installed on the lead's machine). The percentages of the ion's energy and the comparison with the deuterium molecule are the lead's arithmetic. Collective part 2 archive SHA-256 4c2a255493188736f16f178ca1d80fb6174a4931e9e567bea34540b5d90448ac; not rerun (QuTiP and quimb are not installed on the lead's machine). Titanium has no wells yet, so nothing here applies to the titanium foils. Collective part 2b archive SHA-256 19e5928a67c9372d0e04c8d5252f1819acbb94283dd762dc5535dac61b54e0d9; not rerun (QuTiP and quimb are not installed on the lead's machine). `results.json` (64 MB) is committed compressed as `results.json.gz`; unpack it before running `figures.py`. A plain refit by the lead of the three exponents from `transfer.csv` gives the same signs and the same growth with the gap; the cell-number powers come out less steep at gaps of 5 and 9 (-1.7 and -3.0), with a row selection that is not the delivered one. Part 2b uses a 41.6 meV quantum and the 23.85 MeV gap; part 2 used 37 meV and 34 MeV.
