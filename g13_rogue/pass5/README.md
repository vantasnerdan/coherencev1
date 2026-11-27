# Pass 5: tracks, recoils, and launch on a moving source

Read `SUMMARY.md` for the result; `readout.csv` has every source/beam/loss/band/depth, `same_record_peaks.csv` binds peaks to their own exposure and phase shuffle, and `frequency_refinement.csv` resolves the 2.5–250-GHz low-edge depth transfer. The track source puts its energy mainly in 25–250 GHz and 0.25–2.5 THz, so those are the state-dependent bands in `feedback.csv`, `feedback_receiver.csv`, and `feedback_debits.csv`. Units and source ownership are in `inputs_ledger.md`.

## Construction

The accepted central CSDA timelines give mean cumulative deposition, not recoil events. Each row is differenced from the zero launch datum, with timeline Å and fs converted to m and s. D₂⁺ is one electrical parent with two constituent tracks and the accepted total channel energies. The straight track is retained because transverse projectile deflection was not recorded.

The force Fourier convention is `exp(+iωt−ik·r)`. There is one external axial stopping impulse per electrical parent. A ballistic Pd recoil is outside the coherent bath until stopping, and contributes

`q [exp(iω(t+τ)−ik·(r+ℓ)) − exp(iωt−ik·r)]`.
The receiver's conjugate carrier is `exp(+iω(t−t_parent)−ikR)`, so it consumes the **complex conjugate** source phasor. A row at `(z_j,t_j)` therefore shifts its carrier as `ω(t−t_parent−t_j)−k(R−n_z z_j)`. Its frequency/angle-integrated work is retained under the common finite packet filter.

Thus recoil pairs vanish at `(ω,k)=0`; the external momentum is neither added twice nor replaced by a sum of positive recoil magnitudes. Elastic two-body kinematics give `q²=2M_Pd T`, `q_z=(m_D+M_Pd)T/p_D`, and `cosθ=√(T/Tmax)`. ZBL impact-area quadrature and uniform azimuth specify the unresolved marks. Their Poisson intensities reproduce each row's host-energy expectation, including the independent-kick raw second moment `Σλq²=2M_Pd ΔE_host`. The spectrum contains coherent mean power **plus** that recoil variance. Nuclear-only ZBL Pd-on-Pd stopping to 0.1 eV supplies a finite upper-range hypothesis, compared once with one-quarter of that range; constant-deceleration flight gives `τ=2ℓ/v`.
The residual bath force is the declared axial CSDA momentum profile minus the Pd initiation kicks, plus their stopping returns; it can have either sign. This is the explicitly assigned unresolved projectile-deflection/electronic/bath pathway, not an independently recorded electronic force.

A Gaussian 2-nm coherent-bath averaging radius defines this mesoscopic source; it is a model resolution, not a measured lateral track. Electron/direct-D heat stays at its row; host work is deposited along the recoil range and time. The same heat capacity, expansion, conductivity and elastic constants as Pass 4 set diffusion and eigenstress. The remaining heat scale is solved from `hQ+E_wave(h)=Q`, including mechanical, thermal and interference work. Above-ladder acoustic work is included in the debit. Low frequencies use the accepted free-face thickness basis evaluated at actual row depths and times; sub-10-MHz recoil-pair corrections are omitted at their small `ωτ+kℓ` scale.
`kinetic_heat_debits.csv` assigns all mechanical wave work to host kinetic energy first, proportionally over its rows; the remaining host acoustic debit is thermal/interference work. Host residual heat is exactly `h ΔE_host`, and the other deposited channels retain `h` times their energy. This is a nonunique bookkeeping allocation with one total debit, not a second acoustic source.

Bulk outgoing work is `k² dω dΩ / [2ρ(2π)³v_g]` times the projected source quadratic. Its angular energy, coherent phase and recoil variance feed the receiver packets, rather than relabeling a normal-impulse scalar spectrum. A shared real two-point centered mark per electrical parent/branch preserves the computed raw second moment; its finite event powers and source signs are retained. This is a moment closure of unresolved collision marks, not an atomistic trajectory. The receiver retains Pass 4's full foil, seeded grains, first-five-image paths, conservative scalar interface splits, 20-nm/subwavelength readout smoothing and eight-sigma filter. Weak-front phase, late cavity paths and grain mode conversion are not reconstructed by this finite, mode-diagonal envelope.

`readout.csv` uses ensemble effective occupancy `N_eff=mean(packet)²/Var(packet)`, which can be below one. Peak rows use instantaneous positive-parent participation `(Σw_parent)²/Σw_parent²`, with each parent's retained echoes combined first. These differ from geometric ancestry; a paired phase control measures interference in a sampled record.

## State-dependent source

The exact predecessor is `pass2/live_method.md` lines 9–21 and `pass2/live.py` lines 391–459: pre-existing local strain/displacement changes a **probability of launch-energy marks**, not a measured relative-speed coefficient. Here the same bounded-state idea uses signed velocity and displacement projected over actual host-deposition rows before the new parent enters. Quiet same-band built-wave covariance sets fixed references. Each electrical parent has low/mid/high self-energy factors `(.1,1,10)/.865`, probabilities `.25`, `1−.25−p_high`, `p_high=.01+β signal`, and β=0 or .35. Both bands share the actual physical local parent history, with no compressed current or superparents.

A 30-µm source-history patch and 20-ns burn precede the 20-ns observation. This is a source sampling region, unrelated to sub-nm recoil range. Earlier parents are propagated chronologically to each new track row with the shared local Green limit. Only the background **outside** that explicit patch is represented by signed plane-mode covariance; the no-background arm tests this closure. The phase-shuffled and zero-prior arms are paired controls. The separately computed physical correction is `δp/p=−v_z/v_D`, with no fitted amplification. Launch work is `2Re(a_pre*δa)+|δa|²` in the stated normalized source-footprint mode; negative cross work is retained and the two-band incident/heat debit is joined by physical parent ID. Updated launch coefficients actually propagate through the receiver/split kernel. Receiver B uses a fixed same-record quiet threshold and its own exposure; source C uses signed pre-impact velocity and a physical-parent denominator. They are separate populations.
Incoming polarizations and directional group speeds enter the signed projection: a horizontal longitudinal wave is not counted as vertical target motion. The physical source-footprint projection orients eigenvectors consistently before Cartesian projection. Outside-patch covariance is the stationary spatial-beam mean (including for the burst sensitivity); its unknown signed realization is compared with the background-off arm, while the explicit parent history retains the recorded temporal intensity.
Each marked parent's outgoing angular quadratic is normalized explicitly. Its realized branch self-energy, including recoil-variance interference, pays the source debit; the identical normalization divides the packet Green amplitude, so it cannot become a second gain factor.

## Replay

Start at exact accepted repository commit `65dad6ce3e7920a4564b1456f8df2ab3e1c8fe5f`; extract both Pass-5 archives there. NumPy, SciPy and a C99 compiler are needed. Obtain the exact accepted Rung2 numeric/history dependency described in `dependencies.json` (Pass-3 draft 460035582, evidence asset 798321285, archive SHA-256 `e61ab2af46343f68cf682cab4ace136e6fdc8fbddbf5f37148605d87e5845638`). Its two numeric identities are:

- `model.npz`: `59697212d397a6233f7aee9a4d02bd07ac773c8d3d380d50eaba116dffbcef8d`
- `signed_r2.npz`: `55d3d2b85f66d12d77bee585e1dcb385e7e4f307fb2372646f85426ec7849548`

The exact Pass-4 draft is 460102194: capsule asset 798660474, evidence asset 798660627. Its unchanged comparison bytes and finite-front census are reused; no accepted low-band producer or normalization control is rebuilt.

From repository root, with the verified Rung2 extraction at `REF`:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 g13_rogue/pass5/run.py launch --reference-dir REF
python3 g13_rogue/pass5/run.py receivers --reference-dir REF
python3 g13_rogue/pass5/run.py refine --reference-dir REF --previous-evidence-dir PREV
python3 g13_rogue/pass5/run.py feedback --reference-dir REF
python3 g13_rogue/pass5/run.py check --reference-dir REF
python3 g13_rogue/pass5/run.py readout --reference-dir REF
```
`PREV` is the verified Pass-4 evidence extraction. The delivered archives already preserve its exact low-edge catalog bytes, so `refine` needs that argument only on an initial producer reconstruction. Both old Gaussian shapes and the updated source use the same frozen filters and parent catalogs; `tail_comparison.csv` includes energy at the old source's identical peak location/time, not merely separately selected maxima. `replay` regenerates signed fields from existing catalogs without drawing a new source history. To inspect the delivered frozen result, run `check` directly; no producer replay is needed.

`check` reads settled outputs, exercises actual source/packet and loss-unit paths, and uses stored causal kernels to test current-parent exclusion. Numerical spectral arrays are tolerance-portable across supported NumPy/SciPy/BLAS versions; original packed NPZ bytes are preserved in the archives. Packaging never reruns scientific producers. The separate delivery receipt names the unpublished draft and exact remotely read-back archive hashes. No commit, push, merge or publication is performed.
