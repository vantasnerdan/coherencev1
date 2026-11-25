# Collective motion and lattice–nuclear transfer

[Part 1](part1/README.md) models how thermal, track-driven and coherent motion changes the separation and momentum of neighbours. At the same modeled energy, the localized preparation brings deuterons closer than uniform motion. The nearest approach occurs before the kinetic peak. Turning those trajectories into reaction rates needs a matched interaction down to nuclear contact.

[Part 2](part2/README.md) follows energy in one shared vibration and a set of nuclear states. In the small model, the optical bath drains the initial phonons. A high-energy excitation decays through the nuclear channels when those channels are enabled.

[Part 2b](part2b/README.md) is the scaled family: gaps of 3 to 21 quanta, 2 to 8 cells, loss on both sides or one side. Many quanta move at once without loss; finite loss removes the transfer, and the rate falls with each added cell.

[Contact](contact/README.md) turns the computed pair energy into a contact chance and a rate for one pair on the nearest-neighbour line, cold, at 300 K and for the concentrated motion of part 1. All are about a hundred powers of ten below the added signal per pair.

[Focusing](focusing/README.md) searches for starting patterns that gather shared energy onto one pair. At the pumped level the best pattern found puts 0.13 eV on a pair from 2048 deuterons and 0.35 eV from 131,072 (gain 585), while the share of the energy that arrives falls to 0.45 percent. At higher energy per deuteron the peaks stay near 1 eV, and removing the edge of the computed wells does not change them. Adjusting the pattern against the full motion raises the peaks to 0.57 and 1.03 eV at 2048 deuterons and from 0.12 to 0.29 eV at 16,384; a pair closing with 0.97 eV reaches 1.46 A; with loss designed in, 75 to 90 percent of the peak is kept. Turned into a reaction chance pulse by pulse (`focusing/pulse/`), the best pulse gives 1e-95 to 1e-91 per pair, about 75 powers of ten below the signal even with one pulse per beam deuteron. At the pumped level 16,384 and 131,072 deuterons give almost the same closing energy (0.36 eV).

The separate [coupling report](../coupling/README.md) gives the first allowed matrix elements and their limits. Neither calculation changes the foil profile, loading or G12 baseline.
