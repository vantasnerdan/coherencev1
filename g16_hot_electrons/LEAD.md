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

- **Hot electrons stiffen the deuterium's collective motion strongly.** Between 0.07 and 7 eV per site the stiffness rises by roughly 1 eV/Å² for the first eV and 0.23 eV/Å² per eV after that. Against the collective stiffness that is several hundred percent per eV; against the on-site stiffness of about 1.5 eV/Å² used in `g15_track` it is about 70 percent per eV near 0.07 eV per site. Pass 3 of `g15_track` carried 10 percent per eV at most.
- **Mildly heated electrons soften it slightly.** At 0.004 eV per site the stiffness is down by 1 to 2 percent on both meshes. Somewhere between 0.004 and 0.07 eV per site the change turns from softening to stiffening.
- **Below about 0.001 eV per site the change is too small for these meshes.** At the electron energies of the coarse cells of `g15_track` pass 1 (2e-5 to 2e-4 eV per site) the two meshes give opposite signs. Near a track the electron energy per site is far higher than in those cells.
- **No directed force.** A deuteron at the centre of its site feels no net force from heating the electrons, by symmetry. This says nothing about a vacancy cage or uneven heating.
- **Pressure.** At fixed volume the pressure rises by 0.8, 8.8 and 64.7 GPa at 3,000, 10,000 and 30,000 K.
- **Shape of the well.** Along the cube edge the quartic term rises from 11.2 to 14.3 eV/Å⁴ between 300 and 30,000 K; along the cube diagonal it falls from 3.4 to about zero (12³ mesh only).

## What it rests on

1. Electrons with a thermal distribution stand in for the non-thermal electrons of the first femtoseconds.
2. The collective stiffness, not the on-site stiffness and neighbour coupling separately.
3. Fixed volume and fixed palladium; fully loaded PdD.

## Lead's note

Archive SHA-256 6c9aa40a153b54738a2aa7a1b14ae4caaa711996d4a8bd5fa0f49b2f526d3c2f, all 165 file hashes in `SHA256SUMS` checked, no machine paths. The lead recomputed the collective stiffness at 300, 1,000, 3,000, 10,000 and 30,000 K on both meshes directly from the raw forces in `native/`; it agrees with `response.csv` to the digits shown. The first-principles runs themselves were not repeated. `inputs_ledger.md` here holds Giuliano's proposed ledger rows; `g12_residual/inputs_ledger.md` is unchanged.
