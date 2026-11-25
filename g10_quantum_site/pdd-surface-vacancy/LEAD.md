# Wells at the covered PdD(111) face and in a palladium vacancy (lead's note)

Giuliano's capsule. One deuteron at a time is moved by 0.12 and 0.36 A along three fixed directions, with every other atom held, at three sites. The energies are changes in the electronic free energy from the same cell with nothing moved.

## The slab

A 2x2 PdD(111) slab, four palladium layers with the bottom one held, one adsorbed deuteron per four surface palladium atoms, relaxed until the largest force on a free atom is 0.048 eV/A.

- The adsorbed deuteron sits 0.81 A above its three nearest palladium atoms.
- The first deuteron below the surface sits 0.55 A further in than the ideal position used for the arrival run in `g9_edge/first_principles`.
- The adsorbed deuteron and the one below it are 2.44 A apart. The next neighbours of the subsurface deuteron are at 2.60 and 2.62 A. In the bulk every neighbour is at 2.89 A.

| Site | Direction | Cost of 0.12 A (eV), toward minus / plus | Cost of 0.36 A (eV), toward minus / plus |
|---|---|---:|---:|
| Adsorbed | along the surface normal (plus is outward) | 0.027 / 0.020 | 0.193 / 0.198 |
| Adsorbed | in the surface, first axis | 0.019 / 0.019 | 0.207 / 0.179 |
| Adsorbed | in the surface, second axis | 0.015 / 0.020 | 0.102 / 0.271 |
| Subsurface | along the surface normal | 0.009 / 0.010 | 0.083 / 0.076 |
| Subsurface | in the surface, first axis | 0.010 / 0.017 | 0.136 / 0.191 |
| Subsurface | in the surface, second axis | 0.011 / 0.009 | 0.194 / 0.065 |

The subsurface deuteron is held about half as firmly along the normal as the adsorbed one, and its softest directions cost 0.065 to 0.083 eV at 0.36 A.

## The vacancy

The cell is the one from `g7_defect`: seven palladium atoms, one missing, one deuteron. The deuteron's starting position is not a resting position. Moving it 0.36 A further from the centre of the empty site lowers the energy by 0.034 eV and the force still points outward. Sideways it is held (0.06 to 0.40 eV at 0.36 A). So this is a cut from a point that is still sliding, not a well, and Giuliano says so. There is only one deuteron in this cell, so it says nothing about pairs of deuterons in a vacancy.

## What this does and does not give

- It gives one-deuteron wells at two surface sites and the distances between neighbours there. The closest pair found is 2.44 A, against 2.89 A in the bulk.
- It does not give the energy of two deuterons approaching each other at the surface (no pair cut), and the wells stop at 0.36 A.
- No site here has two deuterons starting much closer than in the bulk. Two or more deuterons in one vacancy have not been computed.

## Checks

- Finer real-space grid, one displacement at each of two sites: the energy change moves by 0.0001 eV (adsorbed, 0.12 A along the normal) and 0.0007 eV (vacancy, 0.12 A sideways).
- The k-point mesh was not varied in this capsule.

## Verification

Archive SHA-256 5fb6d577d96413c2bc12ee057f4381c123b7a43e36a0530d6bb5c41fb3b9cd1f, asset 792605823; the delivered `SHA256SUMS` checks (kept as `SHA256SUMS.delivered`). The lead reran `replay/readout_site.py` on the delivered records: `site-wells.json` is reproduced byte for byte. The table and the neighbour distances are the lead's reading of `site-wells.json`.

The delivered records carried the directory names of the machine they were computed on. In the committed copies those are replaced by `<run>/`, so the committed files differ from the delivered ones in those strings only (every number is identical; `SHA256SUMS` is regenerated for the committed files). The replay's grid check compares those strings and stops on the committed copies; it runs on the delivered archive.
