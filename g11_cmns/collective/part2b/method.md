# Model and observables

## Hamiltonian and preparation

Set ℏ=ω=1. The model is

H=δ j + k +(g/√N)(J₊+J₋)(a+a†).

Here j=0,…,N counts nuclear excitations in the exact permutation-symmetric J=N/2 sector. The oscillator number is n₀+k. Each matrix element retains the actual square root of this number. The spin matrix element from j to j+1 is √[(N−j)(j+1)]. The common 1/√N mode normalization is applied once.

Up starts in the bare product state |j=0,n₀>. Its stored energy is n₀ quanta. Down starts in |j=1,n₀−q>. Its stored energy is δ+n₀−q. The excitation is the symmetric bright state. No half-inverted nuclear population is assumed. No displaced-state gauge is used. The large number n₀ is the physical initial Fock occupation in this model.

The target is a joint projector. Up measures |j=1,n₀−q>. Down measures |j=0,n₀>. A smaller off-resonant change in a mean energy is not counted as q-quantum transfer. The outputs also retain the other q-grid endpoints. The largest reported exchange requires at least 1% probability at such an endpoint. All endpoint probabilities remain available below that reporting threshold.

Combined spin/oscillator parity is exact. An even-q destination lies in the other parity block. The q=4 and g=0 controls require no parameter sweep.

The oscillator window extends to |k|≤Nq+max(12,q). Negative physical oscillator numbers are excluded. The retained source parity block is propagated. Boundary probability is recorded for every case.

## Finding the intended resonance

The requested q is the dressed exchange order. It is not the bare gap δ. A large-n guide solves

mean over phase of √[δ²+16g²n₀ cos²(phase)/N]=q.

This supplies the intended branch. A finite Brillouin–Wigner calculation then matches the real diagonal energies of the two endpoints. Every q-grid endpoint remains in the retained block. The other states are eliminated only for this resonance calculation. The time evolution uses the full finite Hamiltonian.

The source benchmark in [Hagelstein and Chaudhary, §6.4](https://arxiv.org/html/0709.3557) uses a bare gap 11 and a dressed order 15. It gives a critical occupation near 27000 for its chosen parameters. This work does not claim that exact benchmark as a replication.

A preliminary continuation selected a bare gap near 14.325 at the nominal 15 case. Its dressed order was 17. That branch was rejected. A nonzero bare 15-photon bin alone was not enough to identify the exchange order. The accepted branch begins near bare 11.918. `branch_proof.json` retains the spectrum and the continued endpoint cluster.

## What loss means here

The source is [Energy Exchange in the Lossy Spin-Boson Model, §3 and Appendix A](https://jcmns.org/api/v1/articles/72149-energy-exchange-in-the-lossy-spin-boson-model.pdf). Its loss acts on intermediate product states. It is not ordinary damping at the optical frequency.

For each q/N/n₀/g case, the product-state energy reference is fixed at its lossless tracked resonance. The lower-only mask contains intermediate states below that reference. The symmetric mask gives the same width to intermediates on both sides. Every q-grid endpoint is stable in both directions. A state within 10⁻¹⁰ quantum of the boundary is stable. The smallest boundary separation is reported.

The mask is not recomputed as its width changes. This keeps the bath fixed. The resonance is retuned along the same branch. Adaptive width steps follow that branch. A lost branch remains unsupported. Its failure is stored instead of selecting a different root.

Finite loss adds −iΓ_j/2 to the Hamiltonian. Lost amplitude enters absorbing sinks. It does not return to the active system. The singular Γ→∞ case removes the lossy intermediate states. That projection is reported separately from every finite width.

## Population, sink and rate

The active-sector wavefunction is not renormalized. The destination probability is |ψ_target|². Survival is ||ψ||². Sink probability is 1−||ψ||². The conditional ratio |ψ_target|²/||ψ||² is retained only to expose the difference.

Let P_max be the largest unconditioned destination probability in the stated window. Let t₅₀ be the first time reaching P_max/2. The reported event metric is r=P_max/(2t₅₀). It has units of ω. Multiplying by ω gives s⁻¹ for that prepared-state protocol. It is not a stationary event count, deposited-power rate or fusion rate.

The default window is 2π/|V_BW|, capped at 10⁷ in ωt. These are splitting-scaled windows. The decisive 15-quantum comparison also uses one common time window. A capped window is marked censored. A matrix-element splitting is reported separately from r.

This Hamiltonian is complex symmetric. Therefore its propagator is symmetric. Forward and reverse destination probabilities are equal for the paired endpoints. Survival and sink can differ. The two directions are retained explicitly.

QuTiP evolves a trace-preserving model with explicit sink jumps. Quimb evolves the matching active-sector wavefunction. The changed-consumer comparison also checks that survival conditioning would give a different answer.

## Fits and extrapolation

The full family varies one coordinate at a time. The orders are 3, 5, 9, 15, 21. N is 2, 4, 8. Occupation is q times 500, 1000, 2000, 4000. The reference coupling is .085√(q/150). Its multipliers are .5, .75, 1, 1.25. Finite widths are 1, 10, 100, 1000. The lower-only projected limit is additional.

`transfer.csv` retains the bare-bin records. Its order-labelled rate is blank when the tracked cluster is unsupported. The raw metric remains separate. The q≥15, N≥4 region is not used for an order-resolved law. Local fits use only retained clusters and uncensored windows. They also exclude boundary probability above 10⁻⁵ and destination peaks below 10⁻¹². These are reporting resolutions, not universal bounds. Large fit residuals are reported.

A separate 100-digit Decimal recurrence sums the leading q-step paths. The common oscillator factor is √[n₀!/(n₀−q)!]. The common coupling factor is g^q. This calculation protects the cancellation control at large N. It is not substituted for the unconditioned time evolution.

The q=5 projected/no-loss amplitude ratio reproduces (N−1)(15N−18)/18. The no-loss coefficient scales as N^{−(q−1)/2} in this weak corner. The projected coefficient approaches an N-independent limit. Finite-width sink depletion is absent from those amplitude statements.

Local gap fits use ln r=a+bq+cq ln q. One-at-a-time exponents are fitted separately. These rate summaries are not carried to the real gap. The only foil-scale extension is the named weak-regime amplitude law:

ln C_q(N)=a+bq+cq ln q+(d+eq)ln(N/1024).

It multiplies g^q√[n₀!/(n₀−q)!]. Its coefficients are fitted to the precise leading-path grid at q=3–21 and N=16–1024. `fits.json` gives the residuals and ranges. This formal amplitude extension gives no destination probability after finite loss.

The foil table maps the archived .6-meV-per-D energy inventory to a conditional Fock occupation. The actual pumped state is not known to be a Fock state. If n<q, that particular Fock source cannot supply the transition. This does not bound an unknown coherent-state tail or external work. Physical odd coupling and active participation remain unknown. The conditional amplitude extension and the absence of observed real-gap population support are separate columns.

**Three omissions:** a physical nuclear coupling and bath spectrum; a measured mode/preparation map to the foil; support for the enormous joint extrapolation beyond the computed coordinates.

## Reproduction

The numerical versions used here are Python 3.12.3, NumPy 2.5.5, quimb 1.16.0 and QuTiP 5.3.1. The figures use Matplotlib 3.11.2.

From the repository root, use an interpreter with those packages. The producer creates its local `runs.jsonl` cache. Analysis selects the declared family and classifies its branches. The release includes those selected results, not the cache's rejected precursor cases.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
python g11_cmns/collective/part2b/scaled.py --suite full
python g11_cmns/collective/part2b/analyze.py
python g11_cmns/collective/part2b/figures.py
```
