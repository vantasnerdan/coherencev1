# p2_input: two real inputs for the chain

Pass 2 replaces placeholders in the pass 0 chain with real values, chosen from the sensitivity ranking in `p1_rank/`. Two were replaced: the heat capacity by lookup, and the bare gap by one electronic-structure calculation.

## Heat capacity, by lookup

    python3 heat_capacity.py

**Cv = 4.1e6 J/(m^3 K)** for Pd/D at loading 0.8 near 300 K, working range 3.5e6 to 4.2e6.

Dulong-Petit for the Pd sublattice (2.43e6) plus an Einstein term for the D optic mode at 37 meV, Einstein temperature 429 K (1.65e6). The range is a spread between models: the lower end uses a Debye Pd host at 276 K with the 65 meV shoulder, the upper end uses the tabulated Pd molar Cp scaled to the loaded cell. It is an estimate. No measured PdD heat capacity near 300 K was found. The placeholder it replaces spanned 1e5 to 1e7.

Sources: Antonov et al., Phys. Rev. B 80, 134302 (2009), pp. 1 and 3 (cell, 37 meV peak, 65 meV shoulder); de Podesta, *Understanding the Properties of Matter*, Table 7.10 (Pd Cp 25.98 J/(mol K)); ISSP Yamamuro group 2012 activity report (Debye temperatures 246 K for PdH, 276 K for Pd).

## Bare gap, by one GPAW calculation

    python3 frozen_phonon.py        (needs GPAW and ASE)

**Harmonic optic frequency at Gamma: 27.1 meV.** The measured PdD optical peak is 37 meV, so the calculation is 9.9 meV, or 27 percent, low.

Setup: two-atom stoichiometric rocksalt PdD at a = 4.084 Angstrom (the loading 1.0 cell, not 0.8), PBE, PAW, plane waves at 500 eV, 6x6x6 k-points, Fermi smearing 0.1 eV, GPAW 26.7.0. The D atom is displaced by plus and minus 0.02 Angstrom along a cubic axis; the force constant is 0.3465 eV/Angstrom^2 from the central difference of the force, and the frequency uses the Pd-D reduced mass (1.977 u). Three SCF runs, about 1.6 core-minutes. No convergence study.

The mismatch with experiment is not resolved into causes by this run. Anharmonicity (PdH and PdD are known to be strongly anharmonic), the difference between a Gamma harmonic frequency and a measured spectral peak, the functional, and the unconverged basis and k-mesh can all contribute. The result has not been re-run on another machine.

Quantum ESPRESSO was tried first, twice, and aborted at startup both times with a buffer overflow before any calculation began. That is a fault in that machine's QE 6.7 installation and says nothing about the material.
