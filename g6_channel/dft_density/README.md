# Pristine PdD ⟨110⟩ valence density

PBE, Pd10/Hae, k2 (8 actual k points), charge0, four Pd plus four D, relaxed cubic a=4.134141684668574Å. The density is GPAW `get_all_electron_density(gridrefinement=2, skip_core=True)` **as returned**, electrons/Å³; frozen Pd core excluded. No new SCF, projectile, vacancy, rescale, atom-core correction or thermal broadening.

## Results

| Quantity | h=0.20Å | h=0.16Å |
|---|---:|---:|
| Native valence grid |40³|48³|
| Integrated electrons (nominal44) |44.18036228667709|44.14122948147666|
| Declared open-hole centre density, e/Å³ |0.23191552|0.23274610|
| Exterior best bPd, Å |0.28511204|0.28571058|
| Exterior best bD, Å |0.25792810|0.26374670|
| bPd / current0.45Å |0.633582|0.634912|
| bD / current0.2645885Å |0.974827|0.996818|
| Exterior relative RMS |2.06295%|2.53245%|
| Full-plane relative RMS of exterior fit |25.2765%|24.5622%|
| Centre overshoot of exterior fit |7.48459%|7.98915%|

The **Native grid** centre h change is at most0.3568594% across all eight preregistered geometric open-hole centres (<5%). This is NOT a5%accuracy claim for the two-width approximation: it overshoots the centre by about8%. Use the real projected grid when that difference matters. The model is a restricted effective fit, not an atomic charge decomposition or a new physical optimum.

## Files and replay

- `density/h020.npz`, `density/h016.npz`: original Native grid bytes, cell, positions, symbols and zero origin. Native uses H for the isotope's electronic setup; its nuclear mass does not change this static electronic density.
- `results/*-transverse.npz`: `n_perp_e_A3`, Fourier coefficients, supported-mode and exterior-fit masks, marked centres and densities, transverse lengths. Array axes are [1−10],[001]; beam is[110]. 128² is a display/interpolation grid, not finer DFT sampling.
- `results/geometry.json`: exact Pd/D string coordinates, line number densities and all centres. `transverse-density.png` marks Pd circles, D squares, centres and the exterior-mask boundary.
- `results/readout.json`: both integrated counts, mean preservation, centre h comparison and Nyquist sensitivity.
- `results/width-fit.json`: both primary fits, nine starts each, Jacobian conditioning, full-plane-fit diagnostic and residuals. Two-width static kernel exp(−r/b)/(8πb³), fixed10/1valence weights, no fitted amplitudes. Fourier bandwidth matches the respective Native projection. Full-plane fitting is not the primary estimator.
- `results/projection-qualification.json`: direct continuous-line quadrature on both actual grids, analytic longitudinal/Nyquist controls and wrong-grid/nonfinite refusals. Kernel radial Fourier quadrature is in the fit report.
- `native/`: immutable source metadata, setup-radius receipt and privacy-redacted GPAW logs. `manifest.json` separates ORIGINAL from EXPORTED byte hashes. GPWs remain private; their original hashes are carried by `inputs.json`.

Run with NumPy only:

```sh
python reader.py CAPSULE_DIRECTORY
```

It checks every exported byte, reprojects both original grids and forward-replays saved fit widths/residuals. It does not rerun GPAW, optimization or channel trajectories. Optional refitting needs SciPy; plotting needs Matplotlib. A separately attached actual-archive portability proof binds the final tar SHA256 and demonstrates blocked Native/SciPy imports plus altered-grid refusal.

## Representation and geometry qualifications

The original h0.20 intake failed a strict44±1e−6 count gate. GPAW explicitly warns skip_core does not generally integrate to nominal valence. Axiom accepted the approximately0.4%discrepancy for this use and directed no correction. That original `passed=false` record is preserved, not rewritten. h0.16 count is also reported, not normalized.

Projection integrates the canonical real trigonometric interpolant continuously along[110], splitting even Nyquist coefficients into Hermitian partners before selecting zero longitudinal modes. Sampled Nyquist aliases cannot identify an unseen continuous field. Omitting Nyquist planes changes centre density by0.002533/0.004030e/Å³ for coarse/fine grids; this convention is explicit. The orthogonal consumer cell is(a√2,a,a√2), twice the original cubic volume, so its projected-cell electron count is twice the original44.18/44.14. Density units remain electrons/Å³, not electrons/Å².

The actual Pd PAW matching radius is max(rcut_j)=2.5Bohr=1.3229430264Å. The5Bohr partial-wave spline truncation is larger and must NOT define the fit mask. Original h0.16 metadata called the spline cutoff `augmentation_support_radii_A`; the corrected distinction is bound in `inputs.json` and the Native setup receipt. Hae has no PAW augmentation. Primary fitting uses the1860/16384 transverse samples outside all Pd matching cylinders.

Declared centres are all rectangular holes in the full Pd+D string lattice, not atom positions. Their nearest Pd string is1.2658172065Å: slightly INSIDE the maximum matching cylinder. Accordingly centre model values are extrapolations from the exterior fit. The source consumer's old lattice is4.084Å; do not confuse it with this relaxed4.1341416847Å cell. Apply any downstream300K Debye-Waller factor separately; no channel yield or refinement was computed here.
