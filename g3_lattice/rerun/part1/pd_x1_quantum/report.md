# Adopted762 Pd local-well refit — staged, not a partial Fleet publication

## Result

The exact adopted **fixed-k3** force curvature and new same-profile quartics lower the Pd spectrum without fitting any output. No DFT or lattice/map dynamics was run.

|isotope|local0→1 before→after,meV|parallel1→2 singlet/doublet,meV|perpendicular1→2,meV|coupledΓ before→after,meV|new per-axis zero-point variance,Å²|new zero-point Hartree shift,eV/Å²|
|---|---|---|---|---|---|---|
|PdH|82.261244→74.034667|88.237035/93.594211|71.438427|80.643686→62.264359|0.0278968602|0.834485576|
|PdD|55.733571→48.883453|57.386734/60.492638|47.336258|54.537036→39.826501|0.0211533583|0.632765561|

Local H/D=1.514513866; coupledΓ H/D=1.563390116. The55.8/37 INS band-descriptor ratio is1.508108108; it is not an identifiedΓ-mode pair. Local-ratio agreement must not hide the differentΓ ratio or remaining absolute mismatch.

**Static branch:** use r(T,ε)=Kii−J0+[K_A(εI)−K_A(0)]+(3A+2B)σ_axis²(T), with the actual PartA polynomial contrast ΔK=−20.4999254225ε+88.5697793051ε² eV/Å². The local variance and quartics are frozen versus strain; the extra elastic change ofΣJ is assumed zero.

|T,K|PdH onset,linear lattice strain|PdD onset|H−D,percentage points|
|---:|---:|---:|---:|
|0|7.456467%|5.270446%|2.186021|
|300|8.769870%|7.050573%|1.719297|

There is **no uniform Hartree onset below5%** in this fitted PartA closure. All listed onsets extrapolate the computed hydrostatic±1% family; the mixed-strain coordinate box is not validation of those hydrostatic strains. The stored JSON supplies soft-mode approaches and distinguishes that sampled DFT family from the coordinate box; neither is an experimental phase boundary.

Do not attribute the old3–4%→new5–9% change to the new quartics alone. Published478 used the primitiveΓ intercept and linear464 law. The saved matched-old401/PartA comparator has no positive-strain polynomial crossing. With NEW matched762 data but the OLD linear464 contrast, H thresholds are5.050547/5.442886% and D4.067278/4.899171% at0/300K: law/reference choice materially changes a5% ordering statement. Historical published values and both clean comparisons remain in result.json.

|isotope|mode/half-pump frequency,THz|stiffness modulation,THz|conditional h_eff threshold|bare lattice h projected through RPA pole|
|---|---:|---:|---:|---:|
|PdH|15.055454941|30.110909882|1.022355%|0.927438%|
|PdD|9.630005194|19.260010387|1.598340%|1.470889%|

These are weak-drive, on-resonance, single-pole thresholds for held physical friction, not incident-field thresholds or measured/nonlinear-lattice onsets. The E² carrier can be at the half/modulated-mode frequency; a different Raman-beat route has its own carriers. No dynamic transduction, transmission or active Z2-switching linewidth is invented.

**Bath/units:** physical κ=0.483554674984ps⁻¹ and amplitude decay4.13603694363ps are held, not a literal gamma=.01. H/D bare harmonic units are49.5097083055/35.0173402220meV; corresponding time units0.013294603815/0.018796743348ps. Atomic velocity-friction coefficients become0.006428667827/0.009089253121. Literal atomic gamma=.01 would instead give amplitude decays2.658921/3.759349ps: a changed physical bath. Reduced-pole bare-force friction additionally carries the spectral mass1.006150745/1.005024402; it is not a replacement for the atomic mass in the full nonlinear lattice.

**Assumptions:** adopted fixed-k3 periodic local cubic well and same-column scalar-J RPA at4.084Å, H as isotope-mass substitution; frozen local variance/quartics with PartA contrast and fixed extraelasticΣJ; phenomenological physical decay held.

## Script

`pd_x1_refit.py` imports the unchanged local3D/RPA implementation in `engine/`. `input.json` is the exact762 force/quartic epoch, not an energy/force average. `result.json` contains full raw strengths, split transitions, pole residues/masses, variance0/300K, static comparisons and bath conventions. The sole changed-Hamiltonian probe compares parity basis8 versus10: maximum transition changes0.000699meV(H)/0.000265meV(D), below the declared0.01meV reporting resolution; this is a finite comparison, not a material uncertainty certificate. The used execution exited0 in4.37s, one thread. The initial scalar-serialization failure was repaired locally, not hidden.

Run with OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=1: `python pd_x1_refit.py > result.json`.

The ONE updated lattice/map batch is still held for actual loading-law delivery plus the producer’s conditional Ni-control budget decision/input784. No partial quantum Fleet observation is posted; this staged result will enter the one combined g3_lattice observation.
