# g13_rogue pass 5: the launch rebuilt from the recorded impact track (lead's summary)

Pass 4 carried the whole foil up to the fast, short waves and found its levels decided by two unmeasured inputs: the size and duration of what an impact launches, and the loss of fast sound in palladium. This pass (Marvin) keeps the same foil, beam, bands and loss rules and replaces the assumed round heated spot with the impact's own recorded track: energy laid down row by row along the 80 to 130 nm path in 0.3 to 0.4 ps, with each recoiling palladium atom treated as a push with its own direction and range. It also lets an impact's launch depend on the motion already present where it lands, in the two bands that hold the energy.

`SUMMARY.md` is Marvin's readout, `README.md` the construction, `inputs_ledger.md` the inputs and their sources, `source.py`, `refine.py`, `feedback.py`, `readout.py` and `run.py` the scripts.

## What the track launches

| Launch | Sound launched by the whole beam | Share of the 0.69 W |
|---|---:|---:|
| Round spot, 60 nm and 2 ps (pass 4) | 1.0 µW | 1.5e-6 |
| Round spot, 20 nm and 0.5 ps (pass 4) | 9.5 µW | 1.4e-5 |
| Recorded track, heat spread over 2 nm | 260 µW | 3.8e-4 |
| Assumed share (a tenth of the recoil energy) | 8,600 µW | 1.25e-2 |

The track launches 27 times more than the tighter round spot. Almost all of it, 253 of the 260 µW, is the sudden heating: the track is a thin line, and a thin hot line launches sound more strongly than a broad warm spot. The recoiling palladium atoms as pushes add about 1 µW, and quartering their range changes the total by 2 percent. The 2 nm width over which the heat is taken to spread is assumed, and it is now the one number that sets the scale of the launch.

The sound moves upward in frequency: 66 percent between 25 and 250 GHz, 33 percent between 0.25 and 2.5 THz.

## The same pattern at a higher level

Even steady beam, stopping depth, copper loss rule.

| Waves | Impacts present at once | Tallest over mean | Mean energy per deuteron, eV |
|---|---:|---:|---:|
| 0.25–2.5 MHz | 1.1e9 | 1.0 | 2e-23 |
| 2.5–25 MHz | 3e8 | 1.0 | 3e-22 |
| 25–250 MHz | 4e7 | 1.1 | 2e-19 |
| 0.25–2.5 GHz | 6e5 | 7 | 1.7e-16 |
| 2.5–25 GHz | 1,400 | 22 | 7.5e-14 |
| 25–250 GHz | 9 | 190 | 2.6e-12 |
| 0.25–2.5 THz | 0.03 | (no event sampled) | 6e-14 |

- **Tallest event sampled:** 0.018 J/m³ at the stopping depth in the 25–250 GHz waves, 57 times the mean at its point and 9e-11 eV per deuteron. As in pass 4 it is one impact's wave: six impacts within reach, one carrying all of it, and scrambling the timing of the others leaves it unchanged.
- **Against heat:** the mean in the 25–250 GHz waves is 6e-8 of the thermal motion already in those waves at 300 K.
- **Depth:** the extra sound is in shorter waves that die sooner. At 30 µm the 25–250 GHz waves hold 1.9e-11 J/m³, less than the 6.6e-11 of pass 4 despite the larger launch. One event there is two impacts adding up (1.8 times the tallest the same record gives with scrambled timing), at 3e-18 eV per deuteron.

## An impact landing on moving ground

In the two bands that hold the energy, each impact was launched into the motion left by the impacts before it, over 20 ns of real beam history in a 30 µm patch.

- **The physical effect is tiny.** The ground under an impact moves at most 2e-9 of the deuteron's speed. The work an impact exchanges with the existing motion is about one part in 100,000 of what it launches, and it is positive for half the impacts and negative for the other half (258 against 260).
- **The pass 2 rule** (motion already present raises the chance of a tenfold launch, strength 0.35) changes the launched work by up to 8 percent. Scrambling the timing of the earlier waves gives the same 8 percent, so the change comes from the rule's own odds and not from where the waves happen to be. The tallest event at the matched point is unchanged.
- **The 0.25–2.5 THz waves** die so fast that none of the 24 sampled windows contained a wave to land on, so this band is untested rather than quiet.

## Loss

A search of published work found no usable figure for the loss of GHz sound in palladium or its hydride (Tanigaki 2011 shows echoes without a fitted loss; Leisure 1983 covers 10 to 270 MHz without a number we could reach). The copper and aluminium rules both stay.

## What the result rests on

1. **The width over which the track's heat appears as lattice heat**, taken as 2 nm and immediate. It sets the launch.
2. **The loss of fast sound in loaded palladium.** It sets how far the sound goes and how many impacts overlap.
3. **Passage of energy from the palladium waves to the deuterons.** Energy per deuteron is wave energy divided by the deuterons present.

## Lead's note

Main archive SHA-256 88718b16e625e32031ced576896dfe2b52069d5a2021d14eee8d088019a68605, all 67 file hashes in `manifest.json` checked, no machine paths. Evidence archive (186 MB) SHA-256 8ad82eaa76f7cc8d328aa7bf702ac37dab99c7a9a13a9ef80598d26dac6cbe31 as reported by Marvin; the lead did not download it and it is not in the repository. One table of the main archive, `feedback_events.csv` (10 MB), is also left out.

The lead reran the launch, the low-edge depth transfer, the state-dependent launch and the readout on a clean copy of the repository at commit 65dad6c with NumPy, SciPy and a C compiler. Launch powers, depth transfer and state-dependent work agree to 1e-11 or better. The step that samples wave records at the four depths was not rerun (it reads the 1 GB history files of pass 3), so the tallest-event figures are as delivered; that step uses the pass 4 receiver code, which the lead did rerun for pass 4.
