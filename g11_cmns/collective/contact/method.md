# Conditional radial model and preparation

This added, unpublished calculation does not alter original Part 1. It asks how a **declared one-path, pure spatial l=0 bound-pair model** responds to moving its G10/G6 join. It does not infer angular projection, pair occupancy, active domains, capture, or physical foil emission. All lengths are Å, energies eV; reduced mass is μ=mD/2. The normalized reduced amplitude obeys ∫|uR|²dR=1 and −[ℏ²/(2μ)]uR″+UuR=EuR, with uR(0)=0. Contact is |uR′(0)|²/(4π), in Å⁻³.

## Geometry, force sign and one-line join

Original indices D25 and D63 have zero-image separation R0=2.88782409436586 along +[110], n=(1,1,0)/√2. The symmetric path is x25=+un, x63=−un, R=R0−2u. The lattice model is

**U(R)=Ulat(R) for R≥Rj; U(R)=VG6(R)+Ulat(Rj)−VG6(Rj) for R<Rj**, with **Ulat=2V.80,110((R0−R)/2)+W((R0−R)/2)** on approach, and Rj=1.28782409436586 or 1.600.

W=EAB−EA−EB+E0 is interpolated by connected cubic Hermite segments through u=0, W=W′=0 and every retained native E/force row. Its derivative is −n·(FAB,25−FA,25)+n·(FAB,63−FB,63). The exported all-atom mixed-force response is not this derivative. There is **no additional old +JAB u²**: Ulat replaces that path term, rather than adding it twice. The adopted one-body is exactly the accepted .80-Å [110] `wide_quadratic` model, including signed-native averaging at .12/.24/.48 and its old polynomial below .03. It is identical in both flag treatments; the .36 one-body value is not a new refit anchor.

All-five uses the five native k3 connected rows (.12/.24/.36/.48/.80); .48 and .80 remain explicitly PAW-overlap flagged, not deleted as invalid data. Three-row excludes only those connected rows and uses **W(u)=W(.36)+W′(.36)(u−.36)** beyond .36: one declared no-fit tangent continuation, not a bound. Its last supported R=2.167824094365859 exceeds both joins, so **neither excluded-row splice has measured distance overlap**. Nominal overlap starts near R=2.146856; flags are partial-wave geometry indicators, not electronic-error estimates. The received finer-mesh zero and projector-free zero/A lack the closest AB calculation; no second-treatment sensitivity is inferred or awaited.

For R>R0 the explicitly chosen cut is **W=0** and Ulat=2V.80,110(|u|). This is C1 at R0 because W(0)=W′(0)=0, but its outward connected force is unmeasured. Beyond |u|=.80, the one-body continuation is E.80+s.80 t+K t²/2, t=|u|−.80, with positive K; it confines the radial exterior. Thus the outer measured-distance limit is R0+1.60, not infinity. These signs/cuts are assumptions, not claims about other lattice directions or a full geometry surface.

G6 is the old **occupied/vacant total interaction**, from a different frozen single-projectile [110] trajectory. It is not the same energy/force chart as G10. The splice is an explicit conditional model, not a validated matching calculation. Both variants use PCHIP on the **electronic remainder** VG6−C/R at the same seven .10–1.60 samples, preserving C=αℏc in eV Å. The single additive shift above sets the same energy-zero rule. A finite derivative mismatch is allowed and recorded in results; no blend is fitted. Below a=.10, the electronic remainder is ε(a)+ε′(a)(R²−a²)/(2a), smooth at the origin and C1 at a. **This is not measured contact data.** A common smooth-core choice is held fixed rather than turning this into another uncertainty campaign. The sampled inner potential is not monotone: PCHIP produces a 1.976-meV shelf between 1.2160 and 1.2796 Å. Its .513-eV localization kinetic scale greatly exceeds that depth, no retained level lies in its energy interval, and the inner floor remains at least .658 eV above the lattice minimum. Thus it does not introduce a dominant inner bound-state population; rare-state contact weights are nevertheless retained and refined.

## Thermal and concentrated states

77/300 K are canonical **radial quantum** eigenstate ensembles, ρT=Σn exp[−(En−E0)/(kBT)]ρn/Z. They are not the original Part 1 classical thermal trajectories or the full 3D orientational ensemble. Logs are formed by log-sum-exp of the density-weighted Boltzmann terms; neither contacts nor rates are averaged as logarithms. Rare-state support is checked at the contact numerator, not merely Z.

Concentrated is a **two-level, projected-turning-work energy-matched analogue**, not a replay of the original 32-site state. The ordinary Part 1 localized D25 +x turning displacement s=.11988147186050971 projects onto D25→D63 as δ=s/√2=.08476900169119068. Its one-sided off-axis actual separation is 2.804336576842573, while the projected radial point is R0−δ=2.8030550926746693. The full-lattice .004844094171855048-eV energy is **not** assigned to this pair.

Define ΔEturn=Ulat(R0−δ)−Ulat(R0), p=ΔEturn/(E1−E0), and Ψ(0)=√(1−p)u0+σ√p u1. This normalized state has exactly that excess energy and real coefficients, hence zero initial current. σ is chosen so 2σ√[p(1−p)]⟨0|R|1⟩ is inward. Its mean displacement is reported, **not asserted equal to δ**. Restriction to levels 0/1 defines the actual preparation; it is not an after-the-fact Fourier cutoff. No translated far tail or poorly resolved high-energy launch coefficients contribute contact.

The time-window average is exact in this two-level model: with both origin slopes gauged positive,

ρavg=(1−p)ρ0+pρ1+2σ√[p(1−p)ρ0ρ1] sin(ωT)/(ωT), ω=(E1−E0)/ℏ, T=2 ps.

The matched reference is the **same joined Hamiltonian's cold radial ground**, not thermal or original-lattice data. Event excess is Ξ=T(ρavg−ρ0), subtracted once. A continuously active tagged pair uses total ρavg in the separate stationary comparison; it is not the event rate. The readout contract gives the conditional one-event-per-incident-D conversion, chart-model current and proton-only coefficient. No capture/all-pair multiplier is inferred.

## Numerical decision and replay

Run `python contact.py --refine` from this folder (NumPy/SciPy; one thread each). `inputs.json` is a compact, frozen exact-binary input projection with source bindings in `provenance.json`. No electronic jobs or old-well edits occur. The finite-difference half-line eigensystem uses uR=vector/√h, so its normalization is the radial integral, not a 3D volume sum. Tiny inner eigenvector tails are **not** trusted. The regular cusp solution uR(R)=R[1+CR/(2t)+O(R²)], t=ℏ²/(2μ), is propagated outward in logarithmic radius using y=(log φ)′, log φ and uR=√Rφ, then matched to the resolved normalized inner-barrier eigen-amplitude. This directly gives log |uR′(0)|²/(4π), without underflow.

The `--refine` comparison changes h=.001→.0005, Rmax=9→10, Ecut=8→10 eV, forbidden match height .15→.25 eV and ODE tolerance 2e−9→3e−10. It reports changes in the **final log contact/event excess**, which are also changes in final log G12 readouts, plus contact-weighted energy support. This is a finite numerical convergence comparison, not a rigorous enclosure or electronic-model uncertainty band. `consumer_proof.json` records the settled source/normalization/readout regression. `results.json` retains all levels/contact logs and all four model variants; the short README table is the decision readout.
