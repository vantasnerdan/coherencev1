# g16_hot_electrons: what hot electrons do to the deuterium's stiffness in PdD (lead's summary)

`g15_track` follows an impact's energy through the foil. About 85 percent of it goes first to the electrons, and passes 2 and 3 there (`../g15_track/pass3/LEAD.md`) carried the effect of hot electrons on the deuterium sites as a guessed setting. This directory (Giuliano) computes it for PdD from first principles: the two-atom cell of `g2_input/`, atoms and volume held fixed, electrons at 300 to 30,000 K.

`README.md` is the technical page, `response.csv` the table, `native/` the raw calculations.

## Result

The stiffness here is the collective one: all deuterons displaced together against the palladium. It is on-site stiffness minus the summed neighbour coupling, and this cell cannot separate the two.

| Electrons at | Extra electron energy per site, eV | Collective stiffness, eV/Å² (16³ mesh) | Change |
|---:|---:|---:|---:|
| 300 K | 0 | 0.1767 | |
| 1,000 K | 0.004 | 0.1732 | −2 % |
| 3,000 K | 0.074 | 0.2531 | +43 % |
| 10,000 K | 0.83 | 0.5793 | +228 % |
| 30,000 K | 7.4 | 1.8846 | +967 % |

- **Hot electrons stiffen the deuterium's collective motion strongly.** The rise is about 1 eV/Å² per eV of electron energy up to 0.07 eV per site, 0.4 from there to 0.8 eV per site, and 0.2 from there to 7 eV per site. Against the collective stiffness that is several hundred percent per eV; against the on-site stiffness of about 1.5 eV/Å² used in `g15_track` it is about 70 percent per eV near 0.07 eV per site. Pass 3 of `g15_track` carried 10 percent per eV at most.
- **Mildly heated electrons soften it slightly.** At 0.004 eV per site the stiffness is down by 1 to 2 percent on both meshes. Somewhere between 0.004 and 0.07 eV per site the change turns from softening to stiffening.
- **Below about 0.001 eV per site the change is too small for these meshes.** At the electron energies of the coarse cells of `g15_track` pass 1 (2e-5 to 2e-4 eV per site) the two meshes give opposite signs. Near a track the electron energy per site is far higher than in those cells.
- **No directed force.** A deuteron at the centre of its site feels no net force from heating the electrons, by symmetry. This says nothing about a vacancy cage or uneven heating.
- **Pressure.** At fixed volume the pressure rises by 0.8, 8.8 and 64.7 GPa at 3,000, 10,000 and 30,000 K.
- **Shape of the well.** Along the cube edge the quartic term rises from 11.2 to 14.3 eV/Å⁴ between 300 and 30,000 K; along the cube diagonal it falls from 3.4 to about zero (12³ mesh only).

## Follow-up: is the stiffening in the site or in the neighbour coupling?

`followup/` (Giuliano; `followup/README.md` is his page) repeats the calculation in a four-atom cell, which also holds the motion where the two deuterium sublattices move against each other along the cube edge. The difference between that and the collective motion gives the coupling between the two sublattices. 16³-equivalent mesh, motion along the wave direction:

| Electrons at | Collective stiffness | Opposite-phase stiffness | Diagonal part | Coupling between sublattices |
|---:|---:|---:|---:|---:|
| 300 K | 0.177 | 0.674 | 0.426 | 0.249 |
| 3,000 K | 0.253 | 0.841 | 0.547 (+29 %) | 0.294 (+18 %) |
| 10,000 K | 0.579 | 1.600 | 1.090 (+156 %) | 0.510 (+105 %) |

All in eV/Å².

- **Both rise.** Hot electrons stiffen the site and strengthen the coupling between neighbours. The vibration runs faster at every wavelength, and the spread of frequencies between long and short waves widens, so a deuterium wave also travels faster where the electrons are hot. Pass 4 of `g15_track` carried the two extremes (all in the site, all in the coupling); the computed answer lies between, about three parts diagonal to one part coupling at 3,000 K.
- **Where the softening turns to stiffening.** Two more electron temperatures place the turn between 1,000 and 1,500 K (0.004 eV per site at the lower end). At 2,000 K both meshes show stiffening (+9 to +13 percent).
- **Not yet a wave speed.** Two wavelengths do not fix how the coupling is shared among the neighbour shells, and the "diagonal part" still contains couplings to neighbours that move in the same phase in this cell. The speed of long deuterium waves is therefore not pinned by this.
- **Sideways motion.** For motion across the wave direction the room-temperature values differ between the two meshes (opposite-phase stiffness 0.065 and 0.129), so only the rise with heating is used: about the same size as along the wave.

Archive SHA-256 b79eb13a390b4d0b3e8b9407c722a71cb33af7e4f4ad6bb5e440a224b5d0c2f1, all 105 file hashes checked, no machine paths. The lead checked the table against `followup/folded-response.json` by hand; the replay scripts need ASE and were not run here.

## What it rests on

1. Electrons with a thermal distribution stand in for the non-thermal electrons of the first femtoseconds.
2. Two wavelengths only (all in step, and neighbouring sublattices opposed); the sharing among neighbour shells and the long-wave speed are open.
3. Fixed volume and fixed palladium; fully loaded PdD.

## Lead's note

Archive SHA-256 6c9aa40a153b54738a2aa7a1b14ae4caaa711996d4a8bd5fa0f49b2f526d3c2f, all 165 file hashes in `SHA256SUMS` checked, no machine paths. The lead recomputed the collective stiffness at 300, 1,000, 3,000, 10,000 and 30,000 K on both meshes directly from the raw forces in `native/`; it agrees with `response.csv` to the digits shown. The first-principles runs themselves were not repeated. `inputs_ledger.md` here holds Giuliano's proposed ledger rows; `g12_residual/inputs_ledger.md` is unchanged.
