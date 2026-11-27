# g13_rogue pass 4: the same whole foil carried up to the fast, short waves (lead's summary)

Continued in `../pass5/LEAD.md`, which rebuilds the launch from the recorded impact track (260 µW against the 1 to 9.5 µW here).

Pass 3 finished the foil's slow vibrations (5 to 250 kHz) and showed that every energy level in them rested on an assumed launch. An impact's energy starts in much faster, shorter waves. This pass (Marvin) keeps the whole foil, the same beam and the same three beam shapes, and carries the waves upward in seven steps of ten, from 0.25 MHz to 2.5 THz. In each step it counts how many impacts' waves are present at one point at one time, and what energy they hold at four depths (the stopping depth of 63 nm, then 30, 50 and 200 µm).

`README.md` is the technical page, `READOUT.md` the full readout, `LEAD.md` Marvin's own short readout, `inputs_ledger.md` the inputs and their sources, `transport.py`, `receiver.py`, `joint.py` and `readout.py` the scripts.

## What an impact launches

The sources are built from the beam itself: the momentum of the deuteron and the sudden heating along its track, with published palladium constants. Per impact they put 1.0e-20 J into sound waves, about one part in 80,000 of the 5 keV. Over the whole beam that is 9.5 µW of the 0.69 W if the heated spot is 20 nm wide and forms in 0.5 ps, and 1.0 µW if it is 60 nm and 2 ps. The spot size and time are assumed; nothing measured fixes them. Nearly all of it lands between 2.5 and 250 GHz, waves roughly 10 nm to 1 µm long.

A second scenario gives sound waves a tenth of the energy that goes into recoiling palladium atoms (about 1.2 percent of the impact). It is an assumption carried for comparison, as the 1 percent launch was in passes 2 and 3.

## Crowded where there is no energy, alone where there is

Even steady beam, stopping depth, built sources, 20 nm and 0.5 ps. Two loss rules are carried because none is measured for palladium at these frequencies: one extrapolated from copper (Ono 2020), one transferred from an aluminium film (Grossmann 2017).

| Waves | Launched, W | Impacts present at once (copper rule / aluminium rule) | Tallest over mean | Mean energy per deuteron, eV | Against heat in the same waves |
|---|---:|---:|---:|---:|---:|
| 0.25–2.5 MHz | 6e-17 | 1.1e9 / 1.1e9 | 1.0 | 2e-23 | 4e-5 |
| 2.5–25 MHz | 6e-15 | 3e8 / 3e8 | 1.0 | 3e-22 | 7e-6 |
| 25–250 MHz | 6e-12 | 4e7 / 2e7 | 1.1 | 2e-19 | 5e-6 |
| 0.25–2.5 GHz | 6e-9 | 6e5 / 9e3 | 9 | 1.5e-16 | 3e-6 |
| 2.5–25 GHz | 2.6e-6 | 290 / 11 | 47 | 1.8e-14 | 4e-7 |
| 25–250 GHz | 6.8e-6 | 3.4 / 0.1 | 200 | 2.6e-14 | 6e-10 |
| 0.25–2.5 THz | 1.2e-7 | 0.03 / 0.001 | (no event sampled) | 4e-16 | 9e-15 |

The picture is hail on a pool seen at two magnifications. In the slow, long waves a billion ripples overlap at every point, and their sum is a smooth, steady hum: the tallest moment equals the mean. In the fast, short waves, where the launched energy is, each impact is a splash that dies away before another lands near it. A point sees a few at a time or none, and the tall moments there are a single splash passing close by.

- **The tallest event sampled** is in the 25–250 GHz step at the stopping depth: 4.4e-5 J/m³, 57 times the mean at that point over the record. Six impacts were within reach at that moment, but one of them carried all of it (weighted count 1.00006), and scrambling the timing of the others leaves 4.43e-5. It is one impact's wave, not several adding up. It corresponds to 2e-13 eV per deuteron.
- **Uneven beam:** hot spots change the levels by about ten percent. Bursts raise the tallest-over-mean of the slow waves from 1.0 to 3.6, which is the burst itself arriving, and leave the fast waves as they were.
- **Depth:** under the copper rule the 2.5–25 GHz waves reach 30 µm with a twentieth of their energy and the 25–250 GHz waves reach it only at their lowest edge (6.6e-11 J/m³, three to four impacts at a time). Under the aluminium rule the fast waves stay within the first micrometres.
- **Against heat:** in every step the launched waves are a small addition to the thermal motion already in the same waves at 300 K, at most 4e-5 of it with the built sources and at most 4e-3 with the assumed share.
- **Assumed share:** the same pattern at about 3,000 to 4,000 times the energy; mean 1e-10 eV per deuteron in the 25–250 GHz step.

Before any of this was run, the fast-wave method was checked against the resolved slow-wave model of pass 3 on a common band. The guided-wave form reproduces it to 6e-15. A plain bulk-ray form gives only 5 to 9 percent of it, so the guided form is used wherever the foil's thickness matters.

## What the result rests on

1. **Loss of fast sound in loaded palladium.** Neither rule is a palladium measurement; they differ by up to a factor of thirty in how many impacts overlap and decide how deep the waves reach.
2. **Size and duration of the launch.** The assumed pair (20 nm and 0.5 ps against 60 nm and 2 ps) moves the launched power ninefold, and the built sources and the assumed share differ by three thousand. The built sources treat the recoiling palladium atoms purely as heat.
3. **Passage of energy from the palladium waves to the deuterons.** Energy per deuteron here is the wave energy divided by the deuterons present, not energy a deuteron has taken up.

The waves are also carried one way and linearly: a wave does not change the launch of the next impact, so the state-dependent launch of pass 2 has no counterpart above 250 kHz yet.

## Lead's note

Main archive SHA-256 f371e37e4ab21d9015523adfaadeac5cef4c61c14ea34ee110c4c911ad060722, all 89 file hashes in `manifest.json` checked, no machine paths. Evidence archive (168 sampled wave records and two distribution files) SHA-256 4833f8451d1f8a5e6e7990a9d333ec3da9ab7c51c010be2517ef81354c088941; it is not in the repository and `evidence_manifest.json` lists its contents. One array of the main archive, `band_depth_full_support.npz` (13 MB), is also left out; its table `band_depth_full_support.csv` is here.

The lead reran every step listed in `READOUT.md` with NumPy, SciPy and a C compiler on a clean copy of the repository at commit 2afd47f plus this directory and the pass 3 step 2 evidence it reads. The launch, the common-band check, the sampled wave records and the tallest event agree to 1e-9 or better (the records to 1e-15). The guided-wave shapes agree to 3 parts in 10,000 and the counts of impacts to 3 parts in a million.
