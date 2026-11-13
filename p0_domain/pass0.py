#!/usr/bin/env python3
"""Coherencev1 pass 0: hypothetical non-RWA Dicke mean field, NOT LENR evidence.

Units: tau=omega_c*t, Ec=hbar*omega_c (free), r=Egap/Ec, theta=kBT/Egap.
Egap is a bare local matter level splitting, never a condensation energy.
q+ip=<a>/sqrt(N), s=<Pauli>, |s|<=1. Per-oscillator internal energy / Ec:
  u=p*p+(1+4*d)*q*q+r*sz/2+2*g*q*sx.
  q'=p; p'=-(1+4*d)*q-g*sx-2*kappa*p;
  s'=cross((4*g*q,0,r),s)-gamma*(s-s_eq(q)),
  s_eq=-tanh(B/(2*r*theta))*(4*g*q,0,r)/B, B=|(4*g*q,0,r)|.
The bath relaxes toward the INSTANTANEOUS local field, not the bare spin.
It is a phenomenological passive thermal bath, not a microscopic metal bath.
Field viscosity removes 4*kappa*p*p from u per tau; spin bath exchanges
-gamma*dot((4*g*q,0,r),s-s_eq)/2. Without either bath, u is conserved.
No thermal noise/finite-N switching is included; equilibria are static.

rho=omega_p/omega_c, omega_p^2=n*e^2/(eps0*m), g=rho*sqrt(f*r)/2.
f is the selected ground-state directional oscillator-strength fraction.
A2 on: d=rho^2/4=g^2/(f*r), the full charged-species TRK coefficient;
A2 off: d=0 is a diagnostic, not consistent electric-dipole minimal coupling.
H/D changes the nuclear mass; default bare gap follows common force constant,
Egap_D/Egap_H=sqrt(mH/mD), unless gap_eV is explicitly supplied.
No host screening, cavity geometry, nuclear reaction, or pumping is inferred.

Primary normalization: Emary/Brandes, doi:10.1103/PhysRevE.67.066203.
TRK: Nataf/Ciuti, https://arxiv.org/abs/1006.1801, Eqs 8,11,13.
Run with NumPy/SciPy and OMP/BLAS threads=1. JSON null means not defined,
not zero. All examples and losses are exploratory unless labeled sourced.
"""
from dataclasses import dataclass, replace, asdict
import argparse
import json
import time
import numpy as np
from scipy.constants import c, e, epsilon_0, hbar, k, m_p, physical_constants
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

MASS = {"H": m_p, "D": physical_constants["deuteron mass"][0]}


@dataclass(frozen=True)
class Model:
    density_m3: float
    species: str = "H"
    Ec_eV: float = 0.1
    gap_H_eV: float = 0.1
    gap_eV: float | None = None
    T_K: float = 0.025 * e / k
    kappa: float = 0.15
    gamma: float = 0.1
    f: float = 1.0
    a2: bool = False

    def __post_init__(self):
        if self.species not in MASS or self.density_m3 < 0 or self.T_K < 0:
            raise ValueError("Require H/D, nonnegative density and temperature")
        if min(self.Ec_eV, self.gap_H_eV) <= 0 or not 0 < self.f <= 1:
            raise ValueError("Require positive energies and 0<f<=1")
        if self.gap_eV is not None and self.gap_eV <= 0:
            raise ValueError("Require positive bare gap")
        if min(self.kappa, self.gamma) < 0:
            raise ValueError("Losses cannot be negative")

    @property
    def omega_c(self):
        return self.Ec_eV * e / hbar

    @property
    def gap(self):
        return (self.gap_eV if self.gap_eV is not None else
                self.gap_H_eV * np.sqrt(MASS["H"] / MASS[self.species]))

    @property
    def r(self):
        return self.gap / self.Ec_eV

    @property
    def temperature(self):
        return k * self.T_K / (e * self.Ec_eV)

    @property
    def rho(self):
        return np.sqrt(self.density_m3 * e**2 / (epsilon_0 * MASS[self.species])) / self.omega_c

    @property
    def g(self):
        return self.rho * np.sqrt(self.f * self.r) / 2

    @property
    def d(self):
        return self.rho**2 / 4 if self.a2 else 0.0

    @property
    def stiffness(self):
        return 1 + 4 * self.d


def thermal_weight(B, m):
    return 1.0 if m.temperature == 0 else np.tanh(B / (2 * m.temperature))


def spin_target(q, m):
    field = np.array([4 * m.g * q, 0.0, m.r])
    B = np.linalg.norm(field)
    return -thermal_weight(B, m) * field / B


def rhs(tau, y, m):
    q, p = y[:2]
    spin = y[2:]
    return np.r_[p, -m.stiffness*q - m.g*spin[0] - 2*m.kappa*p,
                 np.cross([4*m.g*q, 0, m.r], spin)
                 - m.gamma*(spin-spin_target(q, m))]


def energy(y, m):
    q, p, sx, sy, sz = y
    return p*p + m.stiffness*q*q + m.r*sz/2 + 2*m.g*q*sx


def normal_matrix(m):
    w = thermal_weight(m.r, m)
    return np.array([[0, 1, 0, 0],
                     [-m.stiffness, -2*m.kappa, -m.g, 0],
                     [-m.gamma*4*m.g*w/m.r, 0, -m.gamma, -m.r],
                     [4*m.g*w, 0, m.r, -m.gamma]], dtype=float)


def equilibrium(m):
    w = thermal_weight(m.r, m)
    gain = m.f*m.rho**2*w/m.stiffness
    q = 0.0
    if gain > 1:
        upper = 4*m.g**2/m.stiffness
        B = upper if m.temperature == 0 else brentq(
            lambda B: B-upper*thermal_weight(B, m), m.r, upper)
        q = np.sqrt((B-m.r)*(B+m.r))/(4*m.g)
    y = np.r_[q, 0.0, spin_target(q, m)]
    normal = np.array([0, 0, 0, 0, -w])
    delta_u = energy(y, m)-energy(normal, m)
    B = np.hypot(m.r, 4*m.g*q)
    if m.temperature == 0:
        delta_f = m.stiffness*q*q-(B-m.r)/2
    else:
        logz = lambda B: np.logaddexp(B/(2*m.temperature), -B/(2*m.temperature))
        delta_f = m.stiffness*q*q-m.temperature*(logz(B)-logz(m.r))
    return y, float(delta_u), float(delta_f), float(gain)


def density_for_rho(rho, Ec_eV=0.1, species="H"):
    return rho*rho * (Ec_eV*e/hbar)**2 * epsilon_0*MASS[species]/e**2


def domain(m):
    y, du, df, gain = equilibrium(m)
    poles = np.linalg.eigvals(normal_matrix(m))
    growth = float(poles.real.max())
    w = thermal_weight(m.r, m)
    # Frozen thermal polarization, conservative small-excitation references.
    disc = np.sqrt((m.stiffness-m.r**2)**2+16*m.g**2*m.r*w)
    squared = [(m.stiffness+m.r**2-disc)/2,
               (m.stiffness+m.r**2+disc)/2]
    lossless = [float(np.sqrt(x)) if x >= 0 else None for x in squared]
    # For a2=True, (f*w-1)*rho^2>1 is impossible for 0<f,w<=1.
    threshold = None if m.a2 else float(1/np.sqrt(m.f*w)) if w > 0 else None
    ordered = gain > 1
    return {"inputs": asdict(m), "rho": float(m.rho), "r": float(m.r),
            "theta": float(m.temperature/m.r), "g": float(m.g), "d": float(m.d),
            "rho_threshold": threshold,
            "threshold_status": "none at finite density under TRK" if m.a2 else "f*rho^2*tanh(1/(2*theta))>1",
            "normal_seed_exponent_per_tau": growth,
            "ordered": ordered, "equilibrium_q": float(y[0]),
            "locked_omega_over_omega_c": 0.0 if ordered else None,
            "frequency_status": "static Z2 order" if ordered else "no ordered state",
            "domain_size_m": None, "size_status": "no selected geometry/propagation law",
            "bare_vacuum_wavelength_m_NOT_domain_size": float(2*np.pi*c/m.omega_c),
            "delta_U_eV_per_oscillator": float(du*m.Ec_eV),
            "delta_F_eV_per_oscillator": float(df*m.Ec_eV),
            "domain_lifetime_s": None,
            "lifetime_status": "no finite deterministic decay of bath-stabilized equilibrium" if ordered else "no domain to decay",
            "normal_seed_ringdown_s_NOT_domain_lifetime": float(-1/(growth*m.omega_c)) if growth < 0 else None,
            "normal_lossless_mode_omega_over_omega_c": lossless,
            "normal_damped_poles_per_tau": [[float(x.real), float(x.imag)] for x in poles],
            "normal_damped_mode_omega_over_omega_c": sorted(float(x.imag) for x in poles if x.imag > 1e-10)}


def integrate(m, initial, end, times=None):
    sol = solve_ivp(rhs, (0, end), initial, args=(m,), method="DOP853",
                    t_eval=times, rtol=1e-8, atol=1e-13)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def seed_observation(m):
    vals, vecs = np.linalg.eig(normal_matrix(m))
    j = int(np.argmax(vals.real))
    # Biorthogonal modal projection works for a damped oscillatory seed too.
    inverse = np.linalg.inv(vecs)
    initial = 1e-7*vecs[:, j].real/np.max(np.abs(vecs[:, j].real))
    y0 = np.r_[initial, -thermal_weight(m.r, m)]
    sol = integrate(m, y0, 8.0)
    a0 = inverse[j] @ initial
    a1 = inverse[j] @ sol.y[:4, -1]
    return {"observed_seed_exponent_per_tau": float(np.log(abs(a1/a0))/8),
            "linear_seed_exponent_per_tau": float(vals[j].real),
            "modal_amplitude_ratio_at_tau8": float(abs(a1/a0)), "rhs_evaluations": sol.nfev}


def domain_example():
    # Known case answers only whether the actual integrator reproduces its
    # uncoupled damped-field equation; not a numerical convergence study.
    m0 = Model(0)
    times = np.linspace(0, 8, 33)
    known = integrate(m0, [1, 0, 0, 0, -thermal_weight(m0.r, m0)], 8, times)
    omega = np.sqrt(1-m0.kappa**2)
    exact = np.exp(-m0.kappa*times)*(np.cos(omega*times)+m0.kappa/omega*np.sin(omega*times))
    rows = []
    for rho in (0.8, 0.98, 1.05, 1.3):
        m = Model(density_for_rho(rho))
        rows.append(domain(m) | seed_observation(m))
    m = Model(density_for_rho(1.3))
    eq, _, _, _ = equilibrium(m)
    # This is an equilibrium-relaxation observation, NOT a domain decay time.
    initial = eq.copy()
    initial[0] *= 1.01
    settled = integrate(m, initial, 200)
    return {"known_case_uncoupled_q_max_abs_error": float(np.max(abs(known.y[0]-exact))),
            "table": rows,
            "ordered_relaxation": {"rho": 1.3, "q_equilibrium": float(eq[0]),
                                   "q_tau200": float(settled.y[0, -1]),
                                   "rhs_evaluations": settled.nfev},
            "isotope_same_density": [domain(replace(m, species=s)) for s in ("H", "D")]}


def a2_example():
    rows = []
    for rho in (0.8, 0.98, 1.05, 1.3):
        m = Model(density_for_rho(rho))
        rows.append({"rho": rho, "off": domain(m),
                     "on": domain(replace(m, a2=True)) | seed_observation(replace(m, a2=True))})
    return {"comparison": rows,
            "criterion": "off: f*rho^2*w>1; on: (f*w-1)*rho^2>1, impossible for f,w<=1",
            "density_relation": "n=rho^2*epsilon0*m*omega_c^2/e^2"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["domain", "a2"])
    args = parser.parse_args()
    start = time.perf_counter()
    result = domain_example() if args.mode == "domain" else a2_example()
    result["elapsed_s"] = time.perf_counter()-start
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
