"""One-body momentum from *closed* equilibrium paths of distinguishable nuclei.

Write k=p/hbar. For a primitive P-bead sample of H=T+V, the estimator is
Tr[exp(i r.k) B**P]/Tr[B**P], B=exp(-dt T/2) exp(-dt V) exp(-dt T/2),
where dt=beta/P. The translation goes into the kinetic link at the cut;
shifting bead j by y_j*r, y_j=(j+1/2)/P-1/2, makes the extra spring action
m*r**2/(2 beta hbar**2). The remaining factor is exp(-dt sum_j delta V_j).
This is the midpoint displaced-path estimator (Lin et al., PRL 105, 110602,
2010, https://arxiv.org/abs/1004.0868; i-PI 3.3.1 `displacedpath`). No bead
velocities, thermostatted momenta, or position-to-momentum conversion enter.

The finite-P operator is positive but not the continuum Gibbs operator. Its
outer free half steps contribute an exp(-dt E_k) factor to the momentum
diagonal before normalization. A small displayed dt E_k is not evidence of
P convergence. This estimator is for primitive canonical paths, not for a
higher-order, contracted, exchanged, or approximate configurational sampler.
"""

from dataclasses import dataclass
from typing import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.constants import atomic_mass, electron_volt, hbar, k as boltzmann
from scipy.integrate import lebedev_rule

Array = NDArray[np.float64]
SiteCoefficients = tuple[ArrayLike, ArrayLike, float, float, ArrayLike]
Energy = Callable[[Array], ArrayLike]
HBAR2_OVER_2AMU_EV_A2 = hbar**2 * 1e20 / (2 * atomic_mass * electron_volt)
KB_EV_PER_K = boltzmann / electron_volt


@dataclass(frozen=True)
class Characteristic:
    displacement_A: Array
    value: Array
    se: Array
    blocks: Array
    frame_weight_ess: float
    block_weight_ess: float
    numerator_block_ess: Array
    frames_used: int
    direction: Array | None
    angular_order: int | None
    kinetic_factor_eV_A2: float
    beta_over_beads_eV_inv: float
    # [<k^2>, <k^4>] are component moments, or magnitude moments if spherical.
    k2_k4: Array | None
    k2_k4_se: Array | None


@dataclass(frozen=True)
class MomentumInversion:
    k_A_inv: Array
    value: Array
    se: Array
    cutoff_delta: Array
    mesh_delta: Array
    resolution_A_inv: float
    comparison_limit_A_inv: float
    free_end_action: Array

    @property
    def detection_floor(self) -> Array:
        """Empirical absolute floor, not an enclosure or a continuum error."""
        return 3 * self.se + self.cutoff_delta + self.mesh_delta

    @property
    def positive_signal(self) -> NDArray[np.bool_]:
        """Pointwise measured signal against the mesh/window/SE floor only."""
        return (np.abs(self.k_A_inv) < self.comparison_limit_A_inv) & (
            self.value > self.detection_floor
        )


def _directions(direction: ArrayLike | None, order: int) -> tuple[Array, Array]:
    if direction is not None:
        axis = np.asarray(direction, dtype=float)
        if axis.shape != (3,) or not np.all(np.isfinite(axis)) or not np.any(axis):
            raise ValueError("direction must be a finite, nonzero 3-vector")
        axis = axis / np.max(np.abs(axis))
        return (axis / np.linalg.norm(axis))[None, :], np.ones(1)
    xyz, weights = lebedev_rule(order)
    if order < 5 or np.any(weights <= 0):
        raise ValueError("use a positive-weight Lebedev rule of order at least 5")
    x, y, z = xyz
    keep = (x > 0) | ((x == 0) & ((y > 0) | ((y == 0) & (z > 0))))
    half = weights[keep]
    return xyz[:, keep].T, half / half.sum()


def _bridge(p: int, beta: float, kinetic: float, reference: float) -> tuple[Array, float]:
    y = (np.arange(p) + 0.5) / p - 0.5
    prefactor = 1 / (4 * beta * kinetic)
    if reference == 0:
        return y, prefactor
    laplacian = (2 * p / (beta * kinetic)) * np.sin(np.pi * np.fft.rfftfreq(p))**2
    h = np.fft.irfft(np.fft.rfft(y) * (laplacian / (laplacian + beta / p * reference)), n=p)
    # z=h-y is a periodic change of variables. The reference's linear terms
    # vanish by L z + dt*reference*h = 0. Its added action is dt*ref*y.h/2.
    return h, prefactor + 0.5 * beta / p * reference * float(y @ h)


def _potential_terms(x: Array, stiffness: Array, force: Array, a: float, b: float,
                     field: Array, h: Array, u: Array) -> tuple[Array, Array, Array, Array]:
    """Coefficients of sum_j [V(q+h_j*r*u)-V(q)] in powers 1..4 of r."""
    x2 = x * x
    grad = x @ stiffness - force - field + (a - b) * x * x2
    if b:
        grad += b * x * x2.sum(axis=-1, keepdims=True)
    c1 = np.einsum('fbi,b->fi', grad, h) @ u.T
    h2 = h * h
    u2 = u * u
    x2_weighted = np.einsum('fbi,b->fi', x2, h2)
    c2 = 1.5 * (a - b) * (x2_weighted @ u2.T)
    c2 += 0.5 * h2.sum() * np.einsum('di,ij,dj->d', u, stiffness, u)
    if b:
        xx = np.einsum('fbi,fbj,b->fij', x, x, h2)
        c2 += b * np.einsum('di,fij,dj->fd', u, xx, u)
        c2 += 0.5 * b * x2_weighted.sum(axis=-1, keepdims=True)
    c3 = np.einsum('fbi,b->fi', x, h2 * h) @ ((a - b) * u * u2 + b * u).T
    c4 = 0.25 * ((a - b) * (u2 * u2).sum(axis=-1) + b) * (h2 @ h2)
    return c1, c2, c3, c4


def _polynomial_estimate(x: Array, params: tuple[Array, Array, float, float, Array],
                         h: Array, u: Array, direction_weights: Array, r: Array,
                         dt: float, prefactor: float, lw: Array) -> tuple[Array, Array]:
    stiffness, force, a, b, field = params
    c1, c2, c3, c4 = _potential_terms(x, stiffness, force, a, b, field, h, u)
    l1, l2, l3, l4 = -dt * c1, -dt * c2 - prefactor, -dt * c3, -dt * c4
    r2 = r * r
    even = (l2[..., None] + l4[:, None] * r2) * r2
    odd = (l1[..., None] + l3[..., None] * r2) * r
    sample = np.exp(np.logaddexp(even + odd, even - odd) - np.log(2) + lw[:, None, None])
    # Exact Taylor derivatives of the finite-P characteristic, not a PDF fit.
    l1_sq = l1 * l1
    m2 = -2 * l2 - l1_sq
    m4 = 24 * (l4 + l1 * l3) + 12 * l2 * (l2 + l1_sq) + l1_sq * l1_sq
    moments = np.column_stack((m2 @ direction_weights, m4 @ direction_weights))
    return np.einsum('fdr,d->fr', sample, direction_weights), moments * np.exp(lw[:, None])


def _energy_estimate(q: Array, particle: int, energy: Energy, base: Array, h: Array,
                     dirs: Array, weights: Array, r: Array, dt: float, prefactor: float,
                     reference: float, lw: Array) -> Array:
    sample = np.zeros((len(q), len(r)))
    shifted = q.copy()
    site = q[:, :, particle]
    for j, radius in enumerate(r):
        if radius == 0:
            sample[:, j] = np.exp(lw)
            continue
        for u, w in zip(dirs, weights):
            for sign in (-1, 1):
                delta = sign * radius * (h[:, None] * u)
                shifted[:, :, particle] = site + delta
                changed = np.asarray(energy(shifted), dtype=float)
                if changed.shape != q.shape[:2]:
                    raise ValueError("energy_eV must return one total energy per bead")
                dv = (changed - base).sum(axis=1)
                if reference:
                    dv -= reference * (np.einsum('fbi,bi->f', site, delta) + 0.5 * np.sum(delta**2))
                sample[:, j] += 0.5 * w * np.exp(lw - dt * dv - prefactor * radius**2)
    return sample


def _ratio_blocks(numerator: Array, denominator: Array) -> tuple[Array, Array, Array]:
    value = numerator.sum(axis=0) / denominator.sum()
    # Ratio delta-method influence, with the same block denominators used in
    # the central estimate. Linear transforms of these rows preserve covariance.
    blocks = value + (numerator - denominator[:, None] * value) / denominator.mean()
    return value, blocks.std(axis=0, ddof=1) / np.sqrt(len(blocks)), blocks


def _ess(weights: Array) -> Array:
    """Participation of nonnegative weights; not an autocorrelation estimate."""
    peak = weights.max(axis=0)
    scaled = weights / np.where(peak > 0, peak, 1)
    square_sum = np.sum(scaled**2, axis=0)
    return scaled.sum(axis=0)**2 / np.where(square_sum > 0, square_sum, 1)


def displaced_path(
    beads_A: ArrayLike, mass_amu: float, temperature_K: float, particle: int,
    displacement_A: ArrayLike, *, direction: ArrayLike | None = None,
    energy_eV: Energy | None = None, site_coefficients: SiteCoefficients | None = None,
    block_size: int, frame_batch: int = 128, cuts: int = 4, angular_order: int = 7,
    reference_stiffness_eV_per_A2: float = 0.0, log_frame_weights: ArrayLike | None = None,
) -> Characteristic:
    """Estimate the even 1D characteristic along `direction`, or its sphere mean.

    `beads_A`: (ordered frames, beads, occupied sites, 3), using the unchanged
    fixed graph and unwrapped well coordinates. Displace the selected site's
    path only. Supply either a vectorized `(F,P,N,3) -> (F,P)` total bead energy
    or `(K, f, a, b, field)` with eV/Å units and field=sum_j J_ij*x_j, shape
    (F,P,3). The affected energy is x.K.x/2 - (f+field).x + a*sum(x_i^4)/4
    + b*sum_{i<j} x_i^2*x_j^2/2; the neighbor shift is counted once. The latter
    route also yields direct second/fourth moments. In the spherical case the
    angular factors are 3 and 5; multiply by kinetic_factor and its square to
    obtain <E_k> and <E_k^2>. Each positive true moment yields a Markov tail
    bound for its own Gibbs model; the reported statistical SE is not a bound.

    A positive `reference_stiffness_eV_per_A2` changes integration variables
    to a harmonic bridge, exact for any choice (zero recovers i-PI's linear
    shift); it is not a Gaussian continuation for the momentum distribution.
    Antithetic signs and bead cuts are averaged *within* frames. `block_size`
    is in consecutive frames and must be long enough for trajectory and weight
    correlations; a trailing partial block is discarded. `se` is block SE,
    not a discretization, PES, angular, or free-energy-overlap uncertainty.
    `log_frame_weights`, when present, are whole-frame log importance ratios
    to the same V used in the shift. ESS diagnostics are weight participation,
    at frame or block level; even high ESS cannot detect absent overlap.
    """
    q = np.asarray(beads_A, dtype=float)
    r = np.asarray(displacement_A, dtype=float)
    if q.ndim != 4 or q.shape[-1] != 3 or not q.shape[1] or not 0 <= particle < q.shape[2]:
        raise ValueError("expected (frames, beads, sites, 3) Å and a valid occupied particle index")
    if r.ndim != 1 or not r.size or not np.all(np.isfinite(r)) or np.any(r < 0):
        raise ValueError("displacement_A must contain nonnegative finite separations")
    reference = reference_stiffness_eV_per_A2
    if not np.all(np.isfinite([mass_amu, temperature_K, reference])) or min(mass_amu, temperature_K) <= 0 or reference < 0:
        raise ValueError("mass and temperature must be positive; reference stiffness must be nonnegative")
    if block_size < 1 or frame_batch < 1 or cuts < 1 or len(q) < 2 * block_size:
        raise ValueError("need at least two complete blocks and positive frame_batch and cuts")
    if (energy_eV is None) == (site_coefficients is None):
        raise ValueError("supply exactly one of energy_eV or site_coefficients")
    kinetic = HBAR2_OVER_2AMU_EV_A2 / mass_amu
    beta = 1 / (KB_EV_PER_K * temperature_K)
    p = q.shape[1]
    dt = beta / p
    dirs, angular_weights = _directions(direction, angular_order)
    h, prefactor = _bridge(p, beta, kinetic, reference)
    cut_count = min(cuts, p)
    bridges = [np.roll(h, j * p // cut_count) for j in range(cut_count)]
    nblocks = len(q) // block_size
    used = block_size * nblocks
    lw = np.zeros(used) if log_frame_weights is None else np.asarray(log_frame_weights, dtype=float)
    expected = (used,) if log_frame_weights is None else (len(q),)
    if lw.shape != expected or np.any(np.isnan(lw)) or np.any(lw == np.inf) or not np.any(np.isfinite(lw[:used])):
        raise ValueError("expected usable log weights, one per input frame (-inf is allowed)")
    lw = lw[:used] - lw[:used].max()
    weight = np.exp(lw)
    denominator = weight.reshape(nblocks, block_size).sum(axis=1)
    sums = np.zeros((nblocks, len(r)))
    if site_coefficients is not None:
        stiffness, force, a, b, field = site_coefficients
        stiffness = np.asarray(stiffness, dtype=float)
        force = np.asarray(force, dtype=float)
        field = np.asarray(field, dtype=float)
        if stiffness.shape != (3, 3) or not np.allclose(stiffness, stiffness.T, rtol=0, atol=1e-12):
            raise ValueError("site stiffness must be a symmetric (3,3) matrix")
        if force.shape != (3,) or field.shape != (len(q), p, 3):
            raise ValueError("site force and neighbor field must be shaped (3,) and (frames,beads,3)")
        residual_stiffness = stiffness - reference * np.eye(3)
        moment_sums = np.zeros((nblocks, 2))
    for block in range(nblocks):
        end = (block + 1) * block_size
        for start in range(end - block_size, end, frame_batch):
            span = slice(start, min(start + frame_batch, end))
            if not np.any(weight[span]):
                continue
            if energy_eV is not None:
                base = np.asarray(energy_eV(q[span]), dtype=float)
                if base.shape != q[span].shape[:2]:
                    raise ValueError("energy_eV must return one total energy per bead")
                for hcut in bridges:
                    est = _energy_estimate(q[span], particle, energy_eV, base, hcut, dirs,
                                           angular_weights, r, dt, prefactor, reference, lw[span])
                    sums[block] += est.sum(axis=0)
            else:
                local = (residual_stiffness, force, a, b, field[span])
                for hcut in bridges:
                    est, mom = _polynomial_estimate(q[span, :, particle], local, hcut, dirs,
                                                     angular_weights, r, dt, prefactor, lw[span])
                    sums[block] += est.sum(axis=0)
                    moment_sums[block] += mom.sum(axis=0)
    sums /= cut_count
    sums[:, r == 0] = denominator[:, None]
    if not np.all(np.isfinite(sums)):
        raise FloatingPointError("non-finite displaced-path weights: check PES or reweight overlap")
    value, se, blocks = _ratio_blocks(sums, denominator)
    k2_k4 = k2_k4_se = None
    if site_coefficients is not None:
        moment_sums /= cut_count
        if direction is None:
            moment_sums *= [3, 5]
        if not np.all(np.isfinite(moment_sums)):
            raise FloatingPointError("non-finite displaced-path moment weights")
        k2_k4, k2_k4_se, _ = _ratio_blocks(moment_sums, denominator)
    return Characteristic(r, value, se, blocks, float(_ess(weight)), float(_ess(denominator)),
                          _ess(sums), used, dirs[0] if direction is not None else None,
                          None if direction is not None else angular_order, kinetic, dt,
                          k2_k4, k2_k4_se)


def _fourier(result: Characteristic, k_A_inv: ArrayLike, mode: str) -> MomentumInversion:
    r = result.displacement_A
    k = np.asarray(k_A_inv, dtype=float)
    if k.ndim != 1 or not np.all(np.isfinite(k)):
        raise ValueError("k_A_inv must contain finite wavenumbers")
    if (mode == 'component') != (result.direction is not None):
        raise ValueError("component inversion needs a direction; radial inversion needs a sphere mean")
    if mode != 'component' and np.any(k < 0):
        raise ValueError("radial wavenumbers must be nonnegative")
    if len(r) < 9 or r[0] != 0 or r[1] <= 0 or not np.allclose(np.diff(r), r[1], rtol=1e-10, atol=0):
        raise ValueError("inversion requires a uniform displacement grid from zero with at least 9 nodes")
    # Compare with the same last node on the coarser mesh, also for even len(r).
    coarse = np.r_[np.arange(0, len(r) - 1, 2), len(r) - 1]
    shorter = slice(0, 1 + 3 * (len(r) - 1) // 4)

    def kernel(s: Array) -> Array:
        ks = k[:, None] * s
        if mode == 'component':
            f = np.cos(ks) / np.pi
        elif mode == 'radial':
            f = (2 / np.pi) * ks * np.sin(ks)
        else:
            f = (-2 / np.pi) * k[:, None] * (np.sinc(ks / np.pi) - np.cos(ks))
        ds = np.diff(s) / 2
        f[:, :-1] *= np.r_[ds[0], ds[:-1] + ds[1:]]
        f[:, -1] *= ds[-1]
        return f

    transform = kernel(r)
    offset = 1 if mode == 'survival' else 0
    value = offset + transform @ result.value
    block_pdf = result.blocks @ transform.T
    se = block_pdf.std(axis=0, ddof=1) / np.sqrt(len(block_pdf))
    cutoff_delta = np.abs(offset + kernel(r[shorter]) @ result.value[shorter] - value)
    mesh_delta = np.abs(offset + kernel(r[coarse]) @ result.value[coarse] - value)
    return MomentumInversion(k, value, se, cutoff_delta, mesh_delta, np.pi / r[-1],
                             np.pi / (2 * r[1]), result.beta_over_beads_eV_inv * result.kinetic_factor_eV_A2 * k**2)


def radial_momentum(result: Characteristic, k_A_inv: ArrayLike) -> MomentumInversion:
    """Unwindowed finite-transform dP/d|k| [Å], not dP/d³k or dP/dE.

    The sphere-averaged characteristic yields the magnitude distribution even
    for an anisotropic site. Angular-rule sensitivity must be checked separately.
    Oscillatory cancellation can yield negative bins: never clip or normalize
    them. Mesh/window deltas and block SE are observed sensitivities only, not
    bounds on the unsampled characteristic. A floor at or above the value means
    the tail is unresolved there. Momentum has no pointwise displacement cutoff
    corresponding to the end of the fitted well: compare PES continuations.
    """
    return _fourier(result, k_A_inv, 'radial')


def momentum_component(result: Characteristic, k_A_inv: ArrayLike) -> MomentumInversion:
    """Unwindowed finite-transform dP/d(k.direction) [Å], on a signed axis."""
    return _fourier(result, k_A_inv, 'component')


def radial_survival(result: Characteristic, k_A_inv: ArrayLike) -> MomentumInversion:
    """Unclipped finite-transform estimate of P(|p|/hbar > k); cancellation-limited."""
    return _fourier(result, k_A_inv, 'survival')
