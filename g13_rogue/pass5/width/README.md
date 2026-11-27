# Pass 5 addendum: heating width and handover

Launch only, on accepted Pass 5 (`13ac21d20bba353cc7f8df7edd9923fc56f79f1c`). Width is the Gaussian standard deviation of the **heat source**; the mechanical bath form, recoil range and track deposits stay fixed. The electron-to-lattice handover is a causal exponential with the listed mean time; host, direct-D and cutoff heat retain the original timing. Vary one input at a time from 2 nm, instantaneous. Sub-nanometre results use the same continuum thermal model and are a sensitivity, not a calibrated prediction.

Launched power in µW. Each spectral column spans one decade.

| Width (nm) | Handover (ps) | 7-band total | 0.25–2.5 MHz | 2.5–25 MHz | 25–250 MHz | 0.25–2.5 GHz | 2.5–25 GHz | 25–250 GHz | 0.25–2.5 THz |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 0 | 260.2 | 6.228e−11 | 6.069e−9 | 6.078e−6 | 0.006007 | 3.459 | 172.0 | 84.73 |
| 0.5 | 0 | 654.4 | 6.221e−11 | 6.062e−9 | 6.071e−6 | 0.006001 | 3.470 | 203.0 | 448.0 |
| 1 | 0 | 443.5 | 6.225e−11 | 6.066e−9 | 6.075e−6 | 0.006004 | 3.470 | 196.1 | 243.9 |
| 5 | 0 | 92.76 | 6.231e−11 | 6.072e−9 | 6.081e−6 | 0.006007 | 3.378 | 86.23 | 3.150 |
| 10 | 0 | 31.42 | 6.232e−11 | 6.073e−9 | 6.082e−6 | 0.005998 | 3.118 | 27.51 | 0.7889 |
| 2 | 0.1 | 255.2 | 6.228e−11 | 6.069e−9 | 6.079e−6 | 0.006007 | 3.459 | 170.7 | 81.04 |
| 2 | 1 | 126.8 | 6.231e−11 | 6.071e−9 | 6.081e−6 | 0.006009 | 3.415 | 104.5 | 18.88 |
| 2 | 10 | 14.66 | 6.233e−11 | 6.073e−9 | 6.082e−6 | 0.005926 | 1.653 | 9.179 | 3.826 |

With instantaneous heating, widening from 0.5 to 10 nm cuts the seven-band power 20.8-fold. At 2 nm, a 10 ps handover cuts it 17.7-fold against instantaneous heating; 0.1 ps cuts it only 1.9%. Most of the change is above 25 GHz. The 0.5 nm case also launches 4.287 µW above 2.5 THz.

`launch.csv` gives watts in each of the seven bands, their sum, and power above 2.5 THz. Above-ladder power is included when solving the acoustic and residual heat budget, but not in the seven-band sum.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 g13_rogue/pass5/width/launch.py \
  --signed-reference /path/to/verified/rung2/signed_r2.npz
```

Requires NumPy and SciPy. The default signed input is the root-relative Rung 2 path in the accepted Pass 5 dependency ledger. The driver reads only its species energy vector and checks the file hash. The 2 nm, instantaneous case is compared against the accepted `g13_rogue/pass5/source.npz` by species, band and source term before the table is written.
