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

This model is far from the regime the collective-transfer proposals describe. It holds at most four quanta and four cells, and the gap is about a billion quanta wide. It shows the small-number behaviour and gives no scaling law for large numbers. A scaled family (a gap of a few to a few tens of quanta, then size, occupation and loss varied) has been asked for, so that a scaling law comes from runs.

## What this says so far

The collective part shows that the pattern of motion, not its size, is what changes the distance between deuterons, and that concentrated motion lasts. Whether that gives a reaction rate of the size of the added signal cannot be said until the pair energy between 2.4 and 1.3 A is computed. The coupling, as built, acts through the same barrier as the ordinary reaction, so by itself it follows the beam energy and does not give a flat signal.

## Lead's note

Collective part 1 archive SHA-256 54e98dc8425a7c28b36927b8c071c23a6ef75596974259c677d0816bccab00b9. The lead reran `python g11_cmns/collective/part1/motion.py` (NumPy and SciPy only): all 19 cases reproduce, the largest differences being in the energy-balance checks at 1e-9 eV or less; the committed files are the delivered ones. Coupling archive SHA-256 e0e5b2f47389681ff94d6ccac32a582492605822dfce0eb36d01bee1c1784e35; not rerun (the nuclear code is not installed on the lead's machine). The percentages of the ion's energy and the comparison with the deuterium molecule are the lead's arithmetic. Collective part 2 archive SHA-256 4c2a255493188736f16f178ca1d80fb6174a4931e9e567bea34540b5d90448ac; not rerun (QuTiP and quimb are not installed on the lead's machine). Titanium has no wells yet, so nothing here applies to the titanium foils.
