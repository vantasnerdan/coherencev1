## Result

**The updated cold NiH case is Γ-dominated, not striped.** On 110,592 sites, the completed 41.4508 ps trajectory has x-channel amplitude-weighted Γ fraction **0.998434**; the matched equal-time 16-cycle estimator gives **0.982759**. The unit-phase Γ fraction is 0.999952. Its phase-integral equivalent length **3.8525 nm** is essentially the whole periodic box: **box-limited phase locking, not a selected patch radius or infinite-system result**. The old frozen-kernel 2.005 nm/finite-q result remains an old-input result. This is an all-input-profile and retargeted-clock comparison, not an isolated J intervention.

**Temperatures:** the NiH pumped trajectory has Langevin bath **T=0 K**. Of the 30 PdD cases, only `thermal_x0.875_h0` and `thermal_x0.875_h0.08` have **T=300 K**; **all other 28 PdD cases have T=0 K**. `temperatures.json` enumerates every case. Old-kernel NiH 300 K evidence must not be presented as an updated-NiH calculation.

Ni uses all adopted893 harmonic/shell inputs, but its trajectory/static/clock calculations retain **401 quartics in both directions**, frozen before924 arrived. The frozen local first excitation is **100.601694 meV**, coupled pole **85.587939 meV**, versus an **89 meV band peak, not identified Γ**. The additional949-requested **oscillator-only** row uses new100=16.16569694 with old111=1.12836572: the local first excitation is **99.298889 meV**, a **−1.302806 meV (−1.2950%)** shift. Parallel 1→2 lines move **111.095/120.555 →109.285/117.584 meV**, perpendicular **95.634→94.880 meV**. No trajectory, static threshold or clock was recomputed with this late row; no matched-profile new111 exists.

For the frozen Ni epoch the RPA frequency is 20.6951 THz (stiffness clock41.3901 THz), while the actual classical-Hartree phase run uses **21.2416/42.4831 THz**. Preserved physical passive damping is κ=1.259414 ps⁻¹, amplitude decay1.58804 ps; bare dimensionless γ=0.00983398. The simple on-resonance pole thresholds are effective stiffness **h_eff=0.0193710** and bare-code **h=0.0184508**, not incident-field thresholds or measured losses. Frozen-variance static onsets are **8.0652/8.1987% at0/300 K**, explicitly extrapolated Ni-specific contrast closure, not SSCHA boundaries. Finite-record frequency resolution is not a physical linewidth.

### PdD loading/arrangement result

All rows below have T=0 K and pump h=.08. DC is the unpumped vector RMS; AΩ is the pumped vector RMS, each normalized as sqrt(mean over occupied sites of the three-component squared magnitude). ξ_int and the phase peak refer separately to the selected x-channel unit-phase field.

| loading/arrangement | DC, Å | AΩ, Å | ξ_int, nm | phase peak fraction |
|---|---:|---:|---:|---:|
|1.0|0.00000173|0.037968|1.3773|0.11766|
|.875 random|0.071562|0.025432|1.3754|0.11913|
|.75 random|0.092662|0.028151|1.2901|0.09996|
|.6 random|0.105205|0.032259|0.8799|0.03237|
|.875 ordered|0.00000176|0.043201|0.8192|0.02749|
|.875 clustered|0.067420|0.033234|1.7456|0.24366|

Vacancies create **vacancy-driven displacement** when their central forces do not cancel; a centrosymmetric ordered vacancy superlattice cancels that bias rather than automatically making a stronger phase domain. Here ordered occupancy raises amplitude but weakens phase concentration; the generic clustered mask helps finite-time phase concentration. This clustered mask is not the commissioned DFT bunched4/8 motif. Repeated b/c/d measured-force templates test their imposed periodicity, not a measured full Hessian or relaxed structure. Their exact-versus-fit comparisons are in Part1.

At 300 K, PdD AΩ is **0.00914749 without pump /0.00918874 with pump**: no comparable gain. Increasing the cold clock by10% gives AΩ=0.0630669 Å but phase peak0.00701; +5% stretch gives AΩ=0.0486602 Å but ξ_int=0.1721 nm. Larger amplitude is not established phase coherence.

### Bias sensitivity and interface map

On the single prescribed fixed32³ random x=.875 mask, C2=0,+C1/2,−C1/2 gives **amplitude-weighted x-channel Γ fractions 0.00482/0.00148/0.15487**. The minus scenario becomes Γ-led in this channel; the generic three-vector peak can still remain at7.50 nm. Thus the scenario changes selected polarization/carrier texture, not merely size, while remaining the vacancy-displaced/pump-amplified regime. These are scenarios, not calibrated bounds or a thermodynamic phase switch.

The Ni-on-Pd maps use **R[A(E_abs)−A(E_bulk(a(x)))]Rᵀ**, never absolute collective KΓ as an onsite well, with additional elastic ΔΣJ frozen. The expanded pd4084 map raises vector AΩ from0.00133890 to0.0145683 Å; the pd389 map remains tiny (0.00006446→0.00009683 Å) with a stronger forced2Ω response. In the expanded interface layer, Ω energy prefers soft regions (quartile enrichment2.86), but translation mismatch1.35 and soft-map correlation0.248 do **not** show faithful motif-period inheritance. These are per-layer/in-plane finite-map results; an open-depth 3DFFT qz is not a bulk Bloch wavevector. Every EMT site extrapolates PartA, and the hydrogen-free pure-metal geometry is not a validated PdD interface.

The exact two-term864 law and all shell coefficients were retained, including the poor sparse outer-shell held-out error1.535 eV/Å². At random x=.875 about17.96% of sites are outside observable hulls; at x=.75 the outermost-shell fraction is80.44%; x=.6 is entirely extrapolated. Below.75 uses the authorized continuation of the same reference-spacing formula, not an extra measurement. These limitations were not clipped or tuned away.

## Assumptions
Zero-point-renormalized classical FCC Langevin dynamics with frozen x1 local quartics/variance, scalar image-folded reciprocal bonds and vacancy bias, explicit reference-subtracted elastic contrasts/extrapolation, and phenomenological loss; the late100 sensitivity changes only the local oscillator row.

## Script

`updated_batch.py` executed all30 PdD cases; `ni_epoch.py` executed the one matched110592-site NiH case; `ni_oscillator924.py` adds only the requested cheap oscillator row. The packet retains all scalar summaries/configurations, source/input epochs and exact execution evidence; unplotted long curves and bulk per-site trajectories remain on the producing host. All plotted numeric curves and display-layer data are included, with only declared display-slice dtype conversion.

Re-render all three figures **without simulation**:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python updated_figures.py \
  --input updated/figure_input.npz --ni-input updated/ni_figure_input.npz \
  --out replot
```

The fixed-third-primitive-index slices are (001), with raster axes [110]/[-110].
