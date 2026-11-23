# Equations, source identity and replay

At base `59a56b0cd55e601f2e95988106550f3e98bcba5d`, supply the original read-only matrix directory as an argument. The capsule does not contain the private NPZ files:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python g9_edge/edge.py --matrix-dir /path/to/original/data/Pd --output g9_edge/replay
```

Use the existing NumPy/SciPy environment. Output directory must not exist. Single thread; the requested bound is 2 GB RSS and 15 min total numerical CPU. `result.json` records actual CPU, wall time and RSS. Source archive SHA256: `1be3de744451f7938d44bd3dced26854849643e80e8c84f1f35e77991a384c10`. No DFT or Native is used.

## Reduction, support and clocks

Read **saved full-G Coulomb-symmetrized inverse**, then `(epsilon^-1)[G0,G0]`; do not replace it with `1/epsilon[G0,G0]`. Use the original actual energy grid and one matrix version: k12, G200 eV, 96 bands, eta=.1 eV, `primary-N3` and `lambda-G` specimens. Only after that inversion take the bulk isotropic equal-|q| mean. `reduction.csv` includes the wrong-order scalar witness, raw loss and the contracted full-G inverse identity. `dielectric-curve.csv` is static and is not used as a time spectrum. Static imaginary roundoff is recorded separately and L(0)=0 is the explicit spectral choice; positive-frequency losses are not clipped, including near-zero ones.

Between actual energy samples interpolate L=-Im(epsilon^-1) piecewise linearly. Where q is outside the sampled radial support, or energy exceeds 96 eV, use causal zero-temperature Lindhard/Mermin response. Its density is fixed by the smallest-q **finite-support first loss moment**, not by foil fitting; moment clocks are estimates from that realization, not measured plasmons. Mermin collision width .5 eV is a model parameter. The RPA sparse-grid static inverse and its exact piecewise-linear Kramers–Kronig moment need not agree; `static_moments.csv` retains both and the realized hybrid value, with **no static normalization correction**. The moment-derived free-electron continuation is itself a physical choice, not repaired ab initio information.

The finite realization is

`h(q,t) = -theta(t) 2/(pi hbar) exp(-eta_r t/hbar) integral_0^2048 L(q,E) sin(Et/hbar) dE`,

with explicit extra eta_r=.1 eV. Hence the same-model frozen coefficient is

`g0(q) = -2/pi integral E L(q,E)/(E^2+eta_r^2) dE`.

This regulator is not a measured lifetime. The bulk induced interaction for an arbitrary fixed observation point x is

`V_ind(x,t) = (2C/pi) integral_0^16 dq integral_0^t du h(q,u) sinc(q |x-v(t-u) ez|)`.

Thus t<0 is exactly zero, and both geometries have the same source gate. The infinite-prehistory moving bulk wake is the same spectral dispersion at Doppler energy hbar q.mu.v, with that same regulator. The actual G200 short-range convergence problem is inherited, not converted into a numerical error bar. q<=16 Å^-1 and E<=2048 eV are explicit support. qcap=12 sensitivity and a leading Lindhard q^-4 tail estimate are saved separately; neither is a certified tail bound.

## Surface reduction and its limits

The candidate specular-reflection nonlocal half-space kernel, z>0 inside, is

`E_s(k,z)=k/pi integral dq_z exp(i q_z z)/[(k^2+q_z^2) epsilon(Q,omega)]`,

`W_in,in/(2pi C/k)=E_s(z-z')+E_s(z+z')-2 E_s(z)E_s(z')/[1+E_s(0)]`,

`W_in,vac/(2pi C/k)=2 E_s(z) exp(k z')/[1+E_s(0)]`.

Reference: [Pitarke et al., Eqs. 7.30 and 7.48–7.49](https://arxiv.org/abs/cond-mat/0611257), with their inside-z convention reversed. In the local limit E_s=exp(-k|z|)/epsilon; epsilon=1 gives the bare uniform-space interaction, and the two kernels agree at the plane.

This replay takes the cheaper **bulk nonlocal plus local dielectric image** reduction, not the full Q-dependent half-space formula. Let `S=omega(omega+i gamma)/[omega(omega+i gamma)-Omega^2]`. The local image coefficient is `B=S(1-S)/(1+S)=-K_Omega+K_(Omega/sqrt(2))`; its exact damped-oscillator time kernel is used. This retains separate bulk and surface electronic clocks. Integrate `C integral_0^kcap dk J0(k rho) exp(-k z) integral_0^t du h_B(u) exp[-k v(t-u)]`, kcap=1 Å^-1, and show kcap=2 peak/yield sensitivity. Frozen B(0)=0 for this local conductor. At k Å-scale this reduction misses spatial dispersion of the image and microscopic boundary/local-field reconstruction; cutoff sensitivity is not an exact nonlocal-surface error bar. No claim of exact surface ab initio response is made.

## Work, populations and reaction normalization

The force is `-gradient_x V(x,t)` with x held fixed during differentiation. For stationary on-axis sites, interchange source time and lag to integrate impulse analytically in source position; neighbor diagnostics differentiate the actual rho,z geometry. `fixed_sites.csv` is for the finite edge, not a co-moving receiver. Kick energy is p_delta^2/(2m), while pair kinetic correction includes the signed term linear in incident v: `delta E_rel=m[(v+delta v_P-delta v_T)^2-v^2]/4`. These are different quantities.

For a collision at Z, accumulate the impulse until t_match=max[(Z-1 Å)/v,0]. Separately use `delta U= -[delta V(R_turn,t_collision)-delta V(1 Å,t_collision)]` for the inner-turn approximation. This defines a matched, first-order perturbative consumer, not a self-consistent nuclear trajectory. CSDA already includes settled electronic drag: subtract the same spectral model's **formed moving bulk self-force** from the causal projectile self-force before accumulating its kinetic correction. The target perturbation and inner-barrier correction remain causal-minus-frozen. Thus settled stopping is not applied twice. The molecular fragments have their actual per-D incident energies and native five-point spread; profile keys stay at their nominal run energy.

The edge window is from switch-on through the end of the first opposite-sign pair-barrier lobe, reported with both crossings; a missing crossing is explicitly flagged, not silently called a measured duration. Depth quadrature is split at the window boundary. The sheet is exactly t=0; `sheet_formed_bulk_wake` is a deliberately different prehistory sensitivity only for that sheet, not a solved vacuum approach. Outside the selected finite window the common U40 subtotal is retained, not an alleged proof that the later wake equals the adiabatic value. Raw p/n yields and raw ratios are primary; no own-reference normalization is substituted.

The stored-energy table does not turn equilibrium cloud/dipole energy, work function, chemical-potential shift, cohesion or ZPE into an excited release. Unknown release/coupling remains unknown in the physical interpretation; it is not a fitted additive yield contribution. The quantitative rate column is the charged-source electronic/impulse conditional model alone.
