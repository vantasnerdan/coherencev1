---
description: Corrected finite-cell PdD local force column and directional quartics
producer: dense_pd.py
consumer: dense_pd_readout.py
visualization: plot_dense_wells.py
data_directory: dense_Pd/k3
---

# PdD local x=1 input

One D moves while the other31 D and32 Pd nuclei stay at their registered sites in a periodic2x2x2 conventional rocksalt cell, a=4.084A. This is a LOCAL periodic-site response, not the primitive collective Gamma oscillator. The source diagonal retains source-image contributions; it is not an isolated-defect limit.

The adopted electronic profile is PW500/PBE/Pd10, neutral nonmagnetic, Fermi0.1eV, nbands-20, Davidson2, Mixer(.05,5,50), max80, gamma-centred k3, density1e-4 and energy1e-7, with maximum all-atom force change1e-4eV/A. Spatial point-group reduction is disabled; all ordinary displaced states restart the same actual zero-density checkpoint with full SCF. The mesh is campaign-fixed k3, not a new mesh-convergence claim. Nuclear zero-point motion is absent from these Born coefficients.

The looser force-change.003 experiment is not the adopted input. Its force/energy curvature disagreement and large raw acoustic residual were diagnosed before quartic production. A bounded stricter-force continuation of the actual zero/±.03 states preserved geometry, electronic physics, software and PAW bytes; it updated the numerical stopping policy explicitly, without density transplantation. The force-derived on-site curvature and independent energy curvature agree within the declared5percent gate. The full-ion acoustic residual is reported unmodified, not patched tozero. These checks do not supply an error bar for every weak neighbour or prove mesh convergence.

## Actual column and normalization

The linear response uses independent -0.03/+0.03A source-x SCFs and all64 forces: Phi_(i alpha,jx)=-(F_i alpha(+h)-F_i alpha(-h))/(2h). No negative sign is reconstructed. Offsite pulling Jxx=-Phi_(ix,jx); source Kii=Phi_(jx,jx). Kii-sum(all actual offsite Jxx) is the H-only q=0 column sum with metals fixed, not an all-ion acoustic zero or a replacement for local Kii.

Every H row retains its original index, Cartesian response vector, nearest-image offset and image translation. Full distinct-site shell multiplicities are12,3,12,3,1, not unfolded infinite-FCC coordinations. Each occupied periodic site counts once. Do not multiply farther folded shells by an infinite-lattice coordination or add source-image bonds again.

For a scalar lattice kernel, use the mean per existing bond: actual occupied-shell sum divided by occupied multiplicity. Keep the actual sum for the same-column H-only Gamma calculation. At reduced loading, apply occupancy once; an already diluted occupied sum is not a full undiluted kernel. Directional/orientation averaging and extension from this measured source-x column to a full bond tensor are model closures.

## Local directional potential

The100 and111 curves each contain actual positive0.12/0.24A energies and projected forces, with the same zero reference. The final readout uses

V(u)-V(0)=K2*u^2/2+K4*u^4/4,

fixing K2 from the accepted actual central pair and fitting K4 jointly to energy and force, with force residual rows scaled by0.24A. K4 is four times the coefficient multiplying u^4. Energy-only/force-only comparisons and point residuals are retained. Cubic reference symmetry supplies an isotropic quadratic coefficient; even extension of each local potential is an assumption, not extra negative large-displacement measurements. The positive sampled range ends at0.24A; no larger-amplitude validity is demonstrated.

The reported harmonic D quantum is derived from local Kii and the stated deuteron mass. It is not a computed anharmonic excitation gap. Retaining these x1 quartics at lower loading is an explicit downstream approximation; there are no reduced-loading quartic measurements.

## Frozen and Native use

A fresh Native calculation requires the normal GPAW/setups environment and admitted systemd/Herdr resources. Stage order: zero, probe, well100, well111 through run-dense-pd.sh. The optional precision seed flag is restricted to the bounded zero/probe tightening control, validates relocated checkpoint bytes and disallows physical/software/PAW changes or weaker stopping. Complete point JSON/GPW publication is immutable and atomic; orphan checkpoints require explicit disposition. All failures and prior work remain charged.

The frozen readout needs source/config/point JSON files and NumPy/ASE/SciPy, but no private GPW/PAW bytes. Actual isolated consumption without those files and refusal of the wrong3.89A reference are recorded. GPW hashes in the numerical records bind private Native restart material; GPWs are deliberately not in the public small archive. All final quantities and raw measured arrays live in the shared observation archive, not private messages or rounded prose.

The4.084A cell is the campaign's77K beta-PdD author-extrapolated source reference, not PBE's own centred zero-pressure spacing. The separately measured PBE reference and film-negative Gamma curvatures belong to their own geometries/operators; none is substituted as local Kii here.
