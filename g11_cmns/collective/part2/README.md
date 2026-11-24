# Part 2 — what one shared vibration transfers

Does one shared vibration exchange energy with a nuclear excitation, or does each lose its energy separately to the surroundings? We ran the same model in both directions.

### What we put in

The vibration quantum is 37 meV. The nuclear gap is 33.98 MeV, taken from an illustrative odd-parity calculation in the [coupling report](../../coupling/README.md). The coupling in one cell is 6.06 eV for this example. It is not a measured D–D coupling. The leading even-wave D–D term is zero. The physical odd and full second-order terms remain unknown.

We used a single vibration shared by one, two or four cells. Its momentum per cell falls as $1/\sqrt N$ when the total vibration is held fixed. The calculation keeps the sign of each cell’s contribution. Two deuterons moving in opposite directions can cancel in the pair’s total momentum.

The starting states already include the weak mixing caused by the interaction. This avoids injecting extra energy by switching the interaction on suddenly. The optical bath removes low-energy quanta. A separate comparison adds nuclear decay and assigns 25% to photons and 75% to particles. Those fractions are model inputs, not foil measurements.

### What came out

| Starting state or model change | Result |
|---|---|
| Four phonons, no nuclear excitation | The optical bath received 3.97 quanta. The nuclear high-energy band stayed unpopulated. |
| One nuclear excitation, no nuclear decay | The nuclear excitation was retained. |
| One nuclear excitation, nuclear decay enabled | Over 10.34 ps, the baths received 0.248 photon and 0.745 particle quanta in expectation. |
| Increase from one to four cells at fixed mode occupation | The small virtual nuclear admixture stayed the same. |
| Increase from zero to four phonons | The virtual admixture rose from $3.18\times10^{-14}$ to $2.86\times10^{-13}$, following $2n+1$ here. This is mixing, not population of the high-energy band. |

QuTiP and quimb agreed on state populations to within $6.21\times10^{-9}$. The runs used at most four initial phonons and four cells. They do not determine scaling at larger occupation or size.

Four phonons at 37 meV hold 0.148 eV. One keV corresponds to 27,028 of those quanta. The archived pump model holds about 0.6 meV per deuteron, so one keV is the energy inventory of about 1.7 million deuterons at that level. An energy inventory does not specify how many deuterons transfer energy to one pair.

The separately calculated even-wave two-phonon term is $4.37\times10^{-7}$ eV in its finite basis. Its energy mismatch is 34.3 MeV. We did not use it as an on-shell coupling.

**Three omissions:** the physical continuum and higher-order nuclear matrix elements; transfer at much larger phonon occupation and cell count; measured foil mode populations and loss spectra.

Run `python g11_cmns/collective/part2/quantum.py` in the installed scientific environment. `results.json` contains the time histories and energy spectra. `scaling.csv` gives the finite-size comparisons. `coupling_proof.json` binds the input and the ordered two-phonon calculation.

The beam comparisons to make next are the plateau proton-to-neutron ratio, the proton line shape, and changes with loading, temperature or applied vibration.
